import importlib
import sys
from datetime import timedelta as delta

import numpy as np
import parcels
import xarray as xr

sys.path.append("/nethome/4302001/diffusion-hydro-mod/kernels")
import EM_3D_BC
import M1_3D_BC
import Markov1_3D_BC_taper
import density_elements_sampled

dt = 20*60 # seconds

inpath = "/data/oceanparcels/input_data/MITgcm/ACC_channel/"
dataFile = "ACC_ridge_fine_2y_loca.nc"
meanFile = "ACC_ridge_fine_1y_loca_tave.nc"
outFile = f"absolute_dispersion_ACC_ridge_fine_advection_180d_dt{int(dt/60)}m.nc"
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


def create_fieldset(inpath, dataFile, meanFile):
    filenames = {
    "U": inpath + dataFile,
    "V": inpath + dataFile,
    "W": inpath + dataFile,
    "THETA": inpath + dataFile,
    }
    variables = {"U": "UVEL", "V": "VVEL", "W": "WVEL", "THETA": "THETA"}
    dimensions = {
        "U": {"lon": "XG", "lat": "YG", "depth": "Zl", "time": "time"},
        "V": {"lon": "XG", "lat": "YG", "depth": "Zl", "time": "time"},
        "W": {"lon": "XG", "lat": "YG", "depth": "Zl", "time": "time"},
        "THETA": {"lon": "XC", "lat": "YC", "depth": "Z", "time": "time"},
    }
    fieldset = parcels.FieldSet.from_c_grid_dataset(filenames, variables, dimensions, mesh="flat", tracer_interp_method='linear', gridindexingtype='mitgcm')
    
    U_mean = parcels.Field.from_netcdf(inpath+meanFile, ("U_mean", "UVEL"), {"lon": "XG", "lat": "YG", "depth": "Zl"}, mesh='flat', timestamps=None, allow_time_extrapolation=True, interp_method = 'cgrid_velocity')
    V_mean = parcels.Field.from_netcdf(inpath+meanFile, ("V_mean", "VVEL"), {"lon": "XG", "lat": "YG", "depth": "Zl"}, mesh='flat', timestamps=None, allow_time_extrapolation=True, interp_method = 'cgrid_velocity')
    W_mean = parcels.Field.from_netcdf(inpath+meanFile, ("W_mean", "WVEL"), {"lon": "XG", "lat": "YG", "depth": "Zl"}, mesh='flat', timestamps=None, allow_time_extrapolation=True, interp_method = 'cgrid_velocity')
    
    for field in [U_mean, V_mean, W_mean]:
        fieldset.add_field(field)
    
    fieldset.add_periodic_halo(zonal=True)
    fieldset.add_constant("domain_width", 1_000_000)
    
    return fieldset
    
    
class samplingParticle(parcels.JITParticle):
    pot_temp = parcels.Variable('pot_temp', initial=0)
    lon_adjustment = parcels.Variable('lon_adjustment', initial=0)
    inst_U = parcels.Variable('inst_U', initial=0)
    inst_V = parcels.Variable('inst_V', initial=0)
    inst_W = parcels.Variable('inst_W', initial=0)
    mean_U = parcels.Variable('mean_U', initial=0)
    mean_V = parcels.Variable('mean_V', initial=0)
    mean_W = parcels.Variable('mean_W', initial=0)
    
    
def sample(particle, fieldset, time):
    particle.pot_temp = fieldset.THETA[time, particle.depth, particle.lat, particle.lon]
    particle.inst_U = fieldset.U[time, particle.depth, particle.lat, particle.lon]
    particle.inst_V = fieldset.V[time, particle.depth, particle.lat, particle.lon]
    particle.inst_W = fieldset.W[time, particle.depth, particle.lat, particle.lon]
    particle.mean_U = fieldset.U_mean[time, particle.depth, particle.lat, particle.lon]
    particle.mean_V = fieldset.V_mean[time, particle.depth, particle.lat, particle.lon]
    particle.mean_W = fieldset.W_mean[time, particle.depth, particle.lat, particle.lon]
    
    
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
    fieldset = create_fieldset(inpath, dataFile, meanFile)
    
    plat, plon, pdepth = initParticleGrid()
    
    pset = parcels.ParticleSet.from_list(fieldset=fieldset, 
                                         pclass=samplingParticle, 
                                         lon=plon, 
                                         lat=plat, 
                                         depth=pdepth)
    
    sampleKernel = pset.Kernel(sample)
    pset.execute(sampleKernel, dt=0)
    
    output_file = pset.ParticleFile(name = outPath+outFile, outputdt=delta(hours=6))

    pset.execute(
        pset.Kernel(parcels.AdvectionRK4_3D)
        + pset.Kernel(periodicBC)
        + sampleKernel,
        runtime=delta(days=180),
        dt=delta(seconds=dt),
        output_file=output_file,
        recovery={ErrorCode.ErrorInterpolation: deleteParticle,
                  ErrorCode.ErrorOutOfBounds: deleteParticle}
    )

    output_file.close()