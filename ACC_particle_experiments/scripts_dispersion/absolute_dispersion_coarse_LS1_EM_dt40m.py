import importlib
import sys
from datetime import timedelta as delta

import numpy as np
import parcels
import xarray as xr

sys.path.append("/nethome/4302001/diffusion-hydro-mod/kernels")
import density_elements_sampled
import Le_Sommer_elements
import EM_3D_BC
import M1_3D_BC
import K_Le_Sommer
import K_Redi_smallslope

dt = 40*60
kappa = 1500
C = 1
eps = 1e-5/kappa
K_kernel = "Le_Sommer"
AD_kernel = "EM"

# Random seed for reproducibility
parcels.rng.seed(1636)


inpath = "/data/oceanparcels/input_data/MITgcm/ACC_channel/"
dataFile = "ACC_ridge_coarse_1y_locb_fixtemp_with_derivative.nc"
outFile = f"absolute_dispersion_ACC_ridge_coarse_180d_{AD_kernel}_{K_kernel}_C{C}_dt{int(dt/60)}m.nc"
outPath = "/scratch/daanr/trajectories/dispersion/"


K_kernels = {"Redi_smallslope": K_Redi_smallslope.K_Redi_smallslope,
             "Le_Sommer" : K_Le_Sommer.K_Le_Sommer,
             "None" : None
            }
AD_kernels = {"EM": EM_3D_BC.EM_3D_BC,
              "M1": M1_3D_BC.M1_3D_BC,
              "advection" : parcels.kernels.AdvectionRK4_3D}

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


def create_fieldset(inpath, dataFile, derivatives= False, Le_Sommer=False):
    filenames = {
        "U": inpath + dataFile,
        "V": inpath + dataFile,
        "W": inpath + dataFile,
        "THETA": inpath + dataFile,
        "boundaryMask": inpath + dataFile,
    }
    variables = {"U": "UVEL", 
                 "V": "VVEL", 
                 "W": "WVEL", 
                 "THETA": "THETA",
                 "boundaryMask" : "boundaryMask",
                 }
    dimensions = {
        "U": {"lon": "XG", "lat": "YG", "depth": "Zl", "time": "time"},
        "V": {"lon": "XG", "lat": "YG", "depth": "Zl", "time": "time"},
        "W": {"lon": "XG", "lat": "YG", "depth": "Zl", "time": "time"},
        "THETA": {"lon": "XC", "lat": "YC", "depth": "Z", "time": "time"},
        "boundaryMask": {"lon": "XC", "lat": "YC", "depth": "Z"},
    }
    if derivatives == True:
        for derivdim in ['X', 'Y', 'Z', 'XX', 'XY', 'XZ', 'YY', 'YZ', 'ZZ']:
            filenames[f"dRhod{derivdim}"] = inpath + dataFile
            variables[f"dRhod{derivdim}"] = f"dTHETAd{derivdim}"
            dimensions[f"dRhod{derivdim}"] = dict(zip(["time", "depth", "lat", "lon"], ds[f"dTHETAd{derivdim}"].dims))
    
    if Le_Sommer == True:
        for var in ['Delta', 'P', 'Q', 'R']:
            filenames[var] = inpath + dataFile
            variables[var] = var.lower()
            dimensions[var] = dict(zip(["time", "depth", "lat", "lon"], ds[var.lower()].dims))
            for derivdim in ['X', 'Y', 'Z']:
                filenames[f"d{var.lower()}d{derivdim}"] = inpath + dataFile
                variables[f"d{var.lower()}d{derivdim}"] = f"d{var.lower()}d{derivdim}"
                dimensions[f"d{var.lower()}d{derivdim}"] = dict(zip(["time", "depth", "lat", "lon"], ds[f"d{var.lower()}d{derivdim}"].dims))
    
    fieldset = parcels.FieldSet.from_c_grid_dataset(filenames, variables, dimensions, mesh="flat", tracer_interp_method='linear', gridindexingtype='mitgcm', deferred_load=False, allow_time_extrapolation=True)

    fieldset.add_periodic_halo(zonal=True, meridional=True)
    fieldset.add_constant("expansion_terms", 10)
    fieldset.add_constant("domain_width", 1000000)
    fieldset.add_constant("northBound", 2000000)
    fieldset.add_constant("southBound", 50000)
    fieldset.add_constant("upperBound", -10)
    fieldset.add_constant("lowerBound", -3907.58)
    fieldset.add_constant("Sc", 0.008)
    fieldset.add_constant("Sd", 0.0005)
    if Le_Sommer == True:
        fieldset.add_constant("hsquared", C * 50_000**2)
    fieldset.add_constant("kappaval", kappa)
    fieldset.add_constant("epsilon", eps)
    return fieldset
    
    
class samplingParticle(parcels.JITParticle):
    pot_temp = parcels.Variable('pot_temp', initial=0)
    lon_adjustment = parcels.Variable('lon_adjustment', initial=0)
    
    
def sample(particle, fieldset, time):
    particle.pot_temp = fieldset.THETA[time, particle.depth, particle.lat, particle.lon]

    
    
def periodicBC(particle, fieldset, time):
    if particle.lon < 0:
        particle.lon += fieldset.domain_width
        particle.lon_adjustment -= fieldset.domain_width
    elif particle.lon > fieldset.domain_width:
        particle.lon -= fieldset.domain_width
        particle.lon_adjustment += fieldset.domain_width

        
def deleteParticle(particle, fieldset, time):
    particle.delete()

    

if AD_kernel == 'advection':
    fieldset = create_fieldset(inpath, dataFile)
else:
    if K_kernel == "Le_Sommer":
        fieldset = create_fieldset(inpath, dataFile, derivatives=True, Le_Sommer=True)
    else:
        fieldset = create_fieldset(inpath, dataFile, derivatives=True, Le_Sommer=False)

ds.close()

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
pset.execute(sampleKernel, dt=0)

if AD_kernel == "advection":
    kernels = parcels.kernels.AdvectionRK4_3D + pset.Kernel(periodicBC) + SampleTHETAKernel
else:
    if K_kernel == "Le_Sommer":
        kernels = pset.Kernel(density_elements_sampled.density_elements_sampled) \
        + pset.Kernel(Le_Sommer_elements.Le_Sommer_elements) \
        + pset.Kernel(K_kernels[K_kernel]) \
        + pset.Kernel(AD_kernels[AD_kernel]) \
        + pset.Kernel(periodicBC) \
        + sampleKernel
    else:
        kernels = pset.Kernel(density_elements_sampled.density_elements_sampled) \
        + pset.Kernel(K_kernels[K_kernel]) \
        + pset.Kernel(AD_kernels[AD_kernel]) \
        + pset.Kernel(periodicBC) \
        + sampleKernel


output_file = pset.ParticleFile(name = outPath+outFile, outputdt=delta(hours=24))

pset.execute(
    kernels,
    runtime=delta(days=180),
    dt=delta(seconds=dt),
    output_file=output_file,
    recovery={ErrorCode.ErrorInterpolation: deleteParticle,
              ErrorCode.ErrorOutOfBounds: deleteParticle}
)

output_file.close()


