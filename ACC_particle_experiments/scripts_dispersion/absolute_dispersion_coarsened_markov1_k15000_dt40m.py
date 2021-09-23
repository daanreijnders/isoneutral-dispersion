import importlib
import sys
from datetime import timedelta as delta

import numpy as np
import parcels
import xarray as xr

sys.path.append("/nethome/4302001/diffusion-hydro-mod/kernels")
import Markov1_3D_BC_taper
import Markov1_3D_BC_taper_init
import density_elements_sampled

dt = 40*60 # seconds
Tl = 20*24*60*60
kappa = 15000
nusquared = kappa / Tl
eps = dt/Tl #
eta = 1e-5/(eps*nusquared*Tl)

# Random seed for reproducibility
parcels.rng.seed(1636)


inpath = "/data/oceanparcels/input_data/MITgcm/ACC_channel/"
dataFile = "ACC_ridge_fine_1y_locb_coarsened_tave_fixtemp_with_derivative.nc"
outFile = f"absolute_dispersion_ACC_ridge_coarsened_180d_Markov1_K{kappa}_dt{int(dt/60)}m.nc"
outPath = "/scratch/daanr/trajectories/dispersion/"

ds = xr.open_dataset(inpath + dataFile)

def initParticleGrid():
    # First, we create the particles at the path centers:
    center_lats = np.arange(90_000, 2000_000 - 50_000 + 1, 20_000)
    center_lons = np.arange(50_000, 1000_000 - 50_000 + 1, 20_000)
    depths = np.arange(-200, -1600 - 1, -100)

    # create a mesh of these centers
    grid_lats, grid_lons, grid_depths = np.meshgrid(center_lats, center_lons, depths)

    # flatten these
    plats = grid_lats.flatten()
    plons = grid_lons.flatten()
    pdepths = grid_depths.flatten()
    
    return plats, plons, pdepths


def create_fieldset(inpath, dataFile):
    filenames = {
        "U": inpath + dataFile,
        "V": inpath + dataFile,
        "W": inpath + dataFile,
        "THETA": inpath + dataFile,
        "boundaryMask": inpath + dataFile,
        "dRhodX" : inpath + dataFile,
        "dRhodY" : inpath + dataFile,
        "dRhodZ" : inpath + dataFile,
        "dRhodXX" : inpath + dataFile,
        "dRhodXY" : inpath + dataFile,
        "dRhodXZ" : inpath + dataFile,
        "dRhodYY" : inpath + dataFile,
        "dRhodYZ" : inpath + dataFile,
        "dRhodZZ" : inpath + dataFile,
    }
    variables = {"U": "UVEL", 
                 "V": "VVEL", 
                 "W": "WVEL", 
                 "THETA": "THETA",
                 "boundaryMask" : "boundaryMask",
                 "dRhodX" : "dTHETAdX",
                 "dRhodY" : "dTHETAdY",
                 "dRhodZ" : "dTHETAdZ",
                 "dRhodXX" : "dTHETAdXX",
                 "dRhodXY" : "dTHETAdXY",
                 "dRhodXZ" : "dTHETAdXZ",
                 "dRhodYY" : "dTHETAdYY",
                 "dRhodYZ" : "dTHETAdYZ",
                 "dRhodZZ" : "dTHETAdZZ",
                 }
    dimensions = {
        "U": {"lon": "XG", "lat": "YG", "depth": "Zl"},
        "V": {"lon": "XG", "lat": "YG", "depth": "Zl"},
        "W": {"lon": "XG", "lat": "YG", "depth": "Zl"},
        "THETA": {"lon": "XC", "lat": "YC", "depth": "Z"},
        "boundaryMask": {"lon": "XC", "lat": "YC", "depth": "Z"},
    }
    for derivdim in ['X', 'Y', 'Z', 'XX', 'XY', 'XZ', 'YY', 'YZ', 'ZZ']:
        dimensions[f"dRhod{derivdim}"] = dict(zip(["depth", "lat", "lon"], ds[f"dTHETAd{derivdim}"].dims))
    
    fieldset = parcels.FieldSet.from_c_grid_dataset(filenames, variables, dimensions, mesh="flat", tracer_interp_method='linear', gridindexingtype='mitgcm', deferred_load=False)
    
    
    fieldset.add_periodic_halo(zonal=True, meridional=True)
    fieldset.add_constant("domain_width", 1000000)
    fieldset.add_constant("northBound", 2000000)
    fieldset.add_constant("southBound", 50000)
    fieldset.add_constant("upperBound", -10)
    fieldset.add_constant("lowerBound", -3907.58)
    fieldset.add_constant("Sc", 0.008)
    fieldset.add_constant("Sd", 0.0005)
    fieldset.add_constant("TL", Tl)
    fieldset.add_constant("nusquared", nusquared)
    fieldset.add_constant("epsilon", eps)
    fieldset.add_constant("etaval", eta)
    
    return fieldset
    
    
class samplingParticle(parcels.JITParticle):
    pot_temp = parcels.Variable('pot_temp', initial=0)
    lon_adjustment = parcels.Variable('lon_adjustment', initial=0)
    inst_U = parcels.Variable('inst_U', initial=0)
    inst_V = parcels.Variable('inst_V', initial=0)
    inst_W = parcels.Variable('inst_W', initial=0)
    u_prime = parcels.Variable('u_prime', initial=0)
    v_prime = parcels.Variable('v_prime', initial=0)
    w_prime = parcels.Variable('w_prime', initial=0)
    
    
def sample(particle, fieldset, time):
    particle.pot_temp = fieldset.THETA[time, particle.depth, particle.lat, particle.lon]
    particle.inst_U = fieldset.U[time, particle.depth, particle.lat, particle.lon]
    particle.inst_V = fieldset.V[time, particle.depth, particle.lat, particle.lon]
    particle.inst_W = fieldset.W[time, particle.depth, particle.lat, particle.lon]
    
    
def periodicBC(particle, fieldset, time):
    if particle.lon < 0:
        particle.lon += fieldset.domain_width
        particle.lon_adjustment -= fieldset.domain_width
    elif particle.lon > fieldset.domain_width:
        particle.lon -= fieldset.domain_width
        particle.lon_adjustment += fieldset.domain_width
        
        
def deleteParticle(particle, fieldset, time):
    particle.delete()
    
    
if __name__ == '__main__': 
    fieldset = create_fieldset(inpath, dataFile)
    
    plat, plon, pdepth = initParticleGrid()
    
    pset = parcels.ParticleSet.from_list(
        fieldset=fieldset, 
        pclass=samplingParticle, 
        lon=plon, 
        lat=plat, 
        depth=pdepth,
        lonlatdepth_dtype=np.float64,
    )
    
    sampleKernel = pset.Kernel(sample)
    pset.execute(sampleKernel + \
                 density_elements_sampled.density_elements_sampled + \
                 Markov1_3D_BC_taper_init.Markov1_3D_BC_taper_init,
                 dt=0)
    
    output_file = pset.ParticleFile(name = outPath+outFile, outputdt=delta(hours=6))

    pset.execute(
        pset.Kernel(density_elements_sampled.density_elements_sampled)
        + pset.Kernel(Markov1_3D_BC_taper.Markov1_3D_BC_taper)
        + pset.Kernel(periodicBC)
        + sampleKernel,
        runtime=delta(days=180),
        dt=delta(seconds=dt),
        output_file=output_file,
        recovery={ErrorCode.ErrorInterpolation: deleteParticle,
                  ErrorCode.ErrorOutOfBounds: deleteParticle}
    )

    output_file.close()
