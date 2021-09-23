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
import Markov1_3D_BC_taper
import Markov1_3D_BC_taper_init


nparts = 100_000
dt = 40*60 # seconds
Tl = 20*24*60*60
kappa = 1500
nusquared = kappa / Tl
eps = dt/Tl #
eta = 1e-5

# Random seed for reproducibility
parcels.rng.seed(1636)

inpath = "/data/oceanparcels/input_data/MITgcm/ACC_channel/"
refTracer = "ACC_ridge_fine_2y_loca.nc"
dataFile = "ACC_ridge_coarse_1y_locb_fixtemp_with_derivative.nc"
outFile = f"trelease_coarse_locb_Markov1_K{kappa}_Tl{Tl}_p{nparts}_dt{int(dt//60)}m.nc"
outPath = "/scratch/daanr/trajectories/ptracer/"

ds = xr.open_dataset(inpath + dataFile)
ds_ref = xr.open_dataset(inpath + refTracer)


def initializeParticles(reference_ds, level, nparts=50_000):
    """
    Initializes particles in a similar distribution as the initial tracer distribution.
    `nparts` is initial approximation. Actual number of particles used depends on initial
    tracer distribution
    
    Still a work in progress. To change: 
     - ability to specify (X, Y) center of release
     - custom domains, independent from reference dataset
    
    Returns
    -------
    tuple of 3 arrays
        Arrays describing initial X, Y, Z coordinates
    """
    nx = 200  # number of gridpoints in x-direction - 1000km/5km
    ny = 400  # number of gridpoints in y-direction - 2000km/5km
    nz = 30  # number of vertical levels
    
    SALT_3D = np.zeros((nz, ny, nx))
    
    X_release = int(reference_ds.XC[50].data)
    Y_release = int(reference_ds.YC[200].data)
    XX, YY = np.meshgrid(reference_ds.XC.data, reference_ds.YC.data)
    dist_to_release = np.sqrt(np.square(X_release - XX) + np.square(Y_release - YY))
    
    SALT_3D[level, :, :] = np.where(
        dist_to_release / 50000 <= 1,
        (1 / np.power(2 * np.exp(1), dist_to_release / 50000) - 1 / (2 * np.exp(1)))
        / (1 - 1 / (2 * np.exp(1))),
        0,
    )
    
    totalTracer = np.sum(SALT_3D[level, :, :])
    particlesPerUnitTracer = nparts / totalTracer
    particlesPerBin = np.ceil(particlesPerUnitTracer * SALT_3D[level, :, :]).astype("int")
    actualNParticles = np.sum(particlesPerBin)
    
    lon = np.array([])
    lat = np.array([])
    depth = np.ones(actualNParticles) * reference_ds.Z.data[level]
    coord_pairs = np.vstack(
        (XX[dist_to_release < 50_000], YY[dist_to_release < 50_000])
    )
    particlesPerBin_masked = particlesPerBin[dist_to_release < 50_000]
    # Random initialization of particles in each bin
    for idx in range(particlesPerBin_masked.shape[0]):
        lon = np.append(
            lon,
            np.random.uniform(
                coord_pairs[0, idx] - 2500,
                coord_pairs[0, idx] + 2500,
                particlesPerBin_masked[idx],
            ),
        )
        lat = np.append(
            lat,
            np.random.uniform(
                coord_pairs[1, idx] - 2500,
                coord_pairs[1, idx] + 2500,
                particlesPerBin_masked[idx],
            ),
        )

    return np.array(lon), np.array(lat), depth
#


plon, plat, pdepth = initializeParticles(ds_ref, 24, nparts=nparts)
ds_ref.close()


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
                filenames[f"d{var}d{derivdim}"] = inpath + dataFile
                variables[f"d{var}d{derivdim}"] = f"d{var.lower()}d{derivdim}"
                dimensions[f"d{var}d{derivdim}"] = dict(zip(["time", "depth", "lat", "lon"], ds[f"d{var.lower()}d{derivdim}"].dims))
    
    fieldset = parcels.FieldSet.from_c_grid_dataset(filenames, variables, dimensions, mesh="flat", tracer_interp_method='linear', gridindexingtype='mitgcm', deferred_load=False)

    fieldset.add_periodic_halo(zonal=True, meridional=True)
    fieldset.add_constant("expansion_terms", 10)
    fieldset.add_constant("domain_width", 1_000_000)
    fieldset.add_constant("northBound", 2_000_000)
    fieldset.add_constant("southBound", 50_000)
    fieldset.add_constant("upperBound", -10)
    fieldset.add_constant("lowerBound", -3907.58)
    fieldset.add_constant("Sc", 0.008)
    fieldset.add_constant("Sd", 0.0005)
    fieldset.add_constant("TL", Tl)
    fieldset.add_constant("nusquared", nusquared)
    fieldset.add_constant("epsilon", eps)
    fieldset.add_constant("etaval", eta)
    return fieldset


fieldset = create_fieldset(inpath, dataFile, derivatives=True, Le_Sommer=False)
    
    
class MyParticle(parcels.JITParticle):
    pot_temp = parcels.Variable('pot_temp', initial=0)
    lon_adjustment = parcels.Variable('lon_adjustment', initial=0)
    u_prime = parcels.Variable('u_prime', initial=0)
    v_prime = parcels.Variable('v_prime', initial=0)
    w_prime = parcels.Variable('w_prime', initial=0)

pset = parcels.ParticleSet.from_list(fieldset=fieldset, pclass=MyParticle, lon=plon, lat=plat, depth=pdepth)

def SampleTHETA(particle, fieldset, time):
         particle.pot_temp = fieldset.THETA[time, particle.depth, particle.lat, particle.lon]
        
SampleTHETAKernel = pset.Kernel(SampleTHETA)
pset.execute(SampleTHETAKernel + \
             density_elements_sampled.density_elements_sampled + \
             Markov1_3D_BC_taper_init.Markov1_3D_BC_taper_init,
             dt=0)

def periodicBC(particle, fieldset, time):
    if particle.lon < 0:
        particle.lon += fieldset.domain_width
        particle.lon_adjustment -= fieldset.domain_width
    elif particle.lon > fieldset.domain_width:
        particle.lon -= fieldset.domain_width
        particle.lon_adjustment += fieldset.domain_width

def deleteParticle(particle, fieldset, time):
    particle.delete()
    
output_file = pset.ParticleFile(
    name=outPath+outFile, outputdt=delta(hours=24)
)

pset.execute(
    pset.Kernel(density_elements_sampled.density_elements_sampled) \
    + pset.Kernel(Markov1_3D_BC_taper.Markov1_3D_BC_taper)
    + pset.Kernel(periodicBC) \
    + SampleTHETAKernel,
    runtime=delta(days=180),
    dt=delta(seconds=dt),
    output_file=output_file,
    recovery={ErrorCode.ErrorInterpolation: deleteParticle,
              ErrorCode.ErrorOutOfBounds: deleteParticle}
)

output_file.close()