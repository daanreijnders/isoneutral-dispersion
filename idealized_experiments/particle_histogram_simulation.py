import math
import random
import numpy as np
import xarray as xr

from datetime import timedelta
import sys
import warnings


import parcels
from parcels import (
    AdvectionDiffusionEM,
    AdvectionRK4_3D,
    AdvectionDiffusionM1,
    DiffusionUniformKh,
    FieldSet,
    JITParticle,
    ScipyParticle,
    ParticleSet,
    Variable,
    ParcelsRandom,
)
import importlib
import pickle
import pandas as pd
from scipy.integrate import dblquad

sys.path.append("/Users/4302001/surfdrive/diffusion-hydro-mod/kernels")
sys.path.append("/nethome/4302001/diffusion-hydro-mod/kernels")
import idealized_isopycnal_field
import elements_analytical
import K_Le_Sommer
import K_Redi
import K_Redi_smallslope
import EM_3D_BC
import M1_3D_BC
import idealized_isopycnal_field
import Markov1_3D_BC_taper
import Markov1_3D_BC_taper_init
importlib.reload(idealized_isopycnal_field)

def create_analytical_fieldset(
    rhoField,
    Ki=1000,
    Tl=1 * 24 * 60 * 60,
    nu_squared=None,  # RMS eddy velocity, following Koszalka
    epsilon=None,
    eta=None,
    expansion_terms=None,
):
    """
    Prepares a fieldset of an idealized isopycnal field.

    The fieldset can be used for numerically computing trajectories using Markov-1 or Markov-0 (diffusive) parameterizations. 
    Note that in case of Markov-1, `nu_squared` is computed using the assymptotic behavior: Ki = nu_squared * Tl,
    so only `Ki` and `Tl` need to be specified. 
    """
    
    data = {
        "U": np.zeros(1),
        "V": np.zeros(1),
        "W": np.zeros(1),
        "boundaryMask": np.ones(1),
    }
    
    dims = {"lon": 1, "lat": 1, "depth": 1}

    fieldset = parcels.FieldSet.from_data(
        data,
        dims,
        mesh="flat",
        allow_time_extrapolation=True,
    )

    # For analytical solution
    fieldset.add_constant("gval", rhoField.g)
    fieldset.add_constant("Nsquaredval", rhoField.Nsquared)
    fieldset.add_constant("alphaxval", rhoField.alphax)
    fieldset.add_constant("alphayval", rhoField.alphay)
    fieldset.add_constant("kappaxval", rhoField.kappax)
    fieldset.add_constant("kappayval", rhoField.kappay)
    fieldset.add_constant("rho0val", rhoField.rho0)

    # For Milstein
    if expansion_terms:
        fieldset.add_constant("expansion_terms", expansion_terms)

    # For Markov-0
    fieldset.add_constant("kappaval", Ki)

    # For Markov-1
    fieldset.add_constant("TL", Tl)  # Lagrangian timescale
    if Ki:
        fieldset.add_constant("nusquared", Ki / Tl)
    else:
        fieldset.add_constant("nusquared", nu_squared)

    if epsilon:
        fieldset.add_constant("epsilon", epsilon)
    else:
        fieldset.add_constant("epsilon", 0)
        
    if eta:
        fieldset.add_constant("etaval", eta)

    return fieldset

class Markov1Particle(parcels.JITParticle):
    u_prime = parcels.Variable('u_prime', initial=0)
    v_prime = parcels.Variable('v_prime', initial=0)
    w_prime = parcels.Variable('w_prime', initial=0)
    
def get_test_particles(
    rhoField, fieldset, rho=1027.5, nx=25, ny=50, xbounds=(250000, 750000), ybounds=(250000, 1750000), pclass=Markov1Particle
):
    XX, YY = np.meshgrid(
        np.linspace(xbounds[0], xbounds[1], nx), np.linspace(ybounds[0], ybounds[1], ny), indexing='ij'
    )
    lons = XX.flatten()
    lats = YY.flatten()
    depth = rhoField.isopycnal_array(lons, lats, rho_iso=rho)
    pset = parcels.ParticleSet.from_list(
        fieldset,
        pclass=pclass,
        lon=lons,
        lat=lats,
        depth=depth,
        time=np.zeros(nx * ny),
        lonlatdepth_dtype=np.float64,
    )

    return pset

def get_test_particles_dense(
    rhoField, fieldset, rho=1027.5, nx=80, ny=160, xbounds=(0, 1_000_000), ybounds=(0, 2_000_000), pclass=Markov1Particle
):
    """Get a set of evenly spaced testparticles, specified within xbounds and ybounds"""
    
    x_shift = (xbounds[1]-xbounds[0])/nx/2
    y_shift = (ybounds[1]-ybounds[0])/ny/2
    
    if x_shift != y_shift:
        warnings.warn('dx and dy not the same')
    
    XX, YY = np.meshgrid(
        np.linspace(xbounds[0] + x_shift, xbounds[1] - x_shift, nx), np.linspace(ybounds[0] + y_shift, ybounds[1] - y_shift, ny), indexing='ij'
    )
    lons = XX.flatten()
    lats = YY.flatten()
    depth = rhoField.isopycnal_array(lons, lats, rho_iso=rho)
    pset = parcels.ParticleSet.from_list(
        fieldset,
        pclass=pclass,
        lon=lons,
        lat=lats,
        depth=depth,
        time=np.zeros(nx * ny),
        lonlatdepth_dtype=np.float64,
    )

    return pset

def rms_z_error(x, y, z, rho=1027.5):
    """Determines the error (average distance from the isopycnal)"""
    
    true = myField.isopycnal_array(x, y, rho)
    rms_error = np.sqrt(np.mean((z - true)**2))
    return rms_error

def diapycnal_flux(x, y, z, t, rho=1027.5, tmode="seconds"):
    """Determines the diapycnal flux over time"""
    
    rms_error = rms_z_error(x, y, z, rho)
    
    if tmode == "hours":
        multip = 60 * 60
    elif tmode == "seconds":
        multip = 1
    return 0.5 * rms_error**2 / (multip * t)

def plotTestParticles(testParticles):
    fig = go.Figure(data=[go.Surface(x=XX, y=YY, z=ZZ_iso1027, opacity=0.4),
                          go.Scatter3d(x=testParticles.lon, y=testParticles.lat, z=testParticles.depth,
                                       mode='markers', marker=dict(size=1, color='blue'), opacity=0.8)])
    fig.update_layout(
        title="Isopycnal",
        autosize=False,
        width=1200,
        height=600,
        margin=dict(l=65, r=50, b=65, t=90),
        scene=dict(
            aspectratio=dict(x=1, y=2, z=0.5),
            xaxis=dict(title="x [km]"),
            yaxis=dict(title="y [km]"),
            zaxis=dict(title="z [m]", nticks=3),
        ),
    )

    fig.show()
    
def deleteParticle(particle, fieldset, time):
    particle.delete()
    
def periodic_BC(particle, fieldset, time):
    if particle.lon > 1_000_000:
        particle.lon -= 1_000_000
    elif particle.lon < 0:
        particle.lon += 1_000_000
    if particle.lat > 2_000_000:
        particle.lat -= 2_000_000
    elif particle.lat < 0:
        particle.lat += 2_000_000
        
        
myField = idealized_isopycnal_field.densityField()
myField.create_interpolated_density_grid(nx=201, ny=401, nz=901, H=6000)

ZZ_iso1027, XX, YY = myField.isopycnal_grid(1027.5, myField.X, myField.Y)

fieldset = create_analytical_fieldset(myField, Ki=1000)

def add_dummy_settings(fieldset):
    fieldset.add_constant("upperBound", 1e9)
    fieldset.add_constant("lowerBound", -1e9)
    fieldset.add_constant("northBound", 1e9)
    fieldset.add_constant("southBound", -1e9)
    fieldset.add_constant("Sc", 1)
    fieldset.add_constant("Sd", 0.001)
    
# ----- SETTINGS ------
outpath = "/scratch/daanr/markov1_90d_dense_dt40_hist.nc"
runtimedays = 90
outputdays = 90
dt_min = 40
multip = 4
# ---------------------

fieldset = create_analytical_fieldset(myField, Ki=1000, Tl = 20 * 24 * 60 * 60, epsilon=dt_min * 60 / (20 * 24 * 60 * 60), eta=1e-8)
add_dummy_settings(fieldset)


testParticles = get_test_particles_dense(myField, fieldset, nx= multip*80, ny=160*multip, rho=1027.5)
ParcelsRandom.seed(1636)
testParticles.execute(
            testParticles.Kernel(elements_analytical.elements_analytical) 
            + testParticles.Kernel(Markov1_3D_BC_taper_init.Markov1_3D_BC_taper_init),
            dt=0,
            verbose_progress=True,
)
output_file = testParticles.ParticleFile(name = outpath, outputdt=timedelta(days=outputdays))

testParticles.execute(
    testParticles.Kernel(elements_analytical.elements_analytical) 
    + testParticles.Kernel(Markov1_3D_BC_taper.Markov1_3D_BC_taper),
    runtime=timedelta(days=runtimedays),
    dt=timedelta(minutes=40),
    verbose_progress=True,
    output_file=output_file
)

output_file.close()

print(rms_z_error(testParticles.lon, testParticles.lat, testParticles.depth))
diapycnal_flux(testParticles.lon, testParticles.lat, testParticles.depth, 90*24*60*60)