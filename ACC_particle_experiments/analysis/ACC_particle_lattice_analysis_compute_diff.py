import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
import sys
sys.path.append('/nethome/4302001/diffusion-hydro-mod/tools/')
import ACC_tools
import pickle

trajdir = '/data/oceanparcels/output_data/data_Daan/project_diffusion/dispersion/'
tdays = 180

experiment_files = {
    "fine_advection" : "absolute_dispersion_ACC_ridge_fine_advection_180d_dt40m.nc",
    "coarsened_advection" : "absolute_dispersion_ACC_ridge_coarsened_advection_180d_dt40m.nc",
    "coarsened_EM_Le_Sommer_C0.3" : "absolute_dispersion_ACC_ridge_coarsened_180d_EM_Le_Sommer_C0.3_dt40m.nc",
    "coarsened_EM_Le_Sommer_C1" : "absolute_dispersion_ACC_ridge_coarsened_180d_EM_Le_Sommer_C1_dt40m.nc",
    "coarsened_EM_Le_Sommer_C3" : "absolute_dispersion_ACC_ridge_coarsened_180d_EM_Le_Sommer_C3_dt40m.nc",
    "coarsened_EM_Redi_smallslope_K1500" : "absolute_dispersion_ACC_ridge_coarsened_180d_EM_Redi_smallslope_K1500_dt40m.nc",
    "coarsened_EM_Redi_smallslope_K5000" : "absolute_dispersion_ACC_ridge_coarsened_180d_EM_Redi_smallslope_K5000_dt40m.nc",
    "coarsened_EM_Redi_smallslope_K15000" : "absolute_dispersion_ACC_ridge_coarsened_180d_EM_Redi_smallslope_K15000_dt40m.nc",
    "coarsened_M1_Redi_smallslope_K15000" : "absolute_dispersion_ACC_ridge_coarsened_180d_M1_Redi_smallslope_K15000_dt40m.nc",
    "coarsened_Markov1_K1500" : "absolute_dispersion_ACC_ridge_coarsened_180d_Markov1_K1500_dt40m.nc",
    "coarsened_Markov1_K5000" : "absolute_dispersion_ACC_ridge_coarsened_180d_Markov1_K5000_dt40m.nc",
    "coarsened_Markov1_K15000" : "absolute_dispersion_ACC_ridge_coarsened_180d_Markov1_K15000_dt40m.nc",
    "coarse_advection" : "absolute_dispersion_ACC_ridge_coarse_advection_180d_dt40m.nc",
    "coarse_EM_Le_Sommer_C0.3" : "absolute_dispersion_ACC_ridge_coarse_180d_EM_Le_Sommer_C0.3_dt40m.nc",
    "coarse_EM_Le_Sommer_C1" : "absolute_dispersion_ACC_ridge_coarse_180d_EM_Le_Sommer_C1_dt40m.nc",
    "coarse_EM_Le_Sommer_C3" : "absolute_dispersion_ACC_ridge_coarse_180d_EM_Le_Sommer_C3_dt40m.nc",
    "coarse_EM_Redi_smallslope_K1500" : "absolute_dispersion_ACC_ridge_coarse_180d_EM_Redi_smallslope_K1500_dt40m.nc",
    "coarse_EM_Redi_smallslope_K5000" : "absolute_dispersion_ACC_ridge_coarse_180d_EM_Redi_smallslope_K5000_dt40m.nc",
    "coarse_EM_Redi_smallslope_K15000" : "absolute_dispersion_ACC_ridge_coarse_180d_EM_Redi_smallslope_K15000_dt40m.nc",
    "coarse_Markov1_K1500" : "absolute_dispersion_ACC_ridge_coarse_180d_Markov1_K1500_dt40m.nc",
    "coarse_Markov1_K5000" : "absolute_dispersion_ACC_ridge_coarse_180d_Markov1_K5000_dt40m.nc",
    "coarse_Markov1_K15000" : "absolute_dispersion_ACC_ridge_coarse_180d_Markov1_K15000_dt40m.nc",
}

experiments = dict()
for experiment, filename in experiment_files.items():
    experiments[experiment] = xr.open_dataset(trajdir + str(filename), decode_times=False)
    
ds_field_fine = xr.open_dataset('/data/oceanparcels/input_data/MITgcm/ACC_channel/ACC_ridge_fine_1y_locb.nc')
ds_field_coarsened = xr.open_dataset("/data/oceanparcels/input_data/MITgcm/ACC_channel/ACC_ridge_fine_1y_locb_coarsened_fixtemp_with_derivative.nc")
ds_field_coarse = xr.open_dataset('/data/oceanparcels/input_data/MITgcm/ACC_channel/ACC_ridge_coarse_1y_locb_tave_fixtemp_with_derivative.nc')
    
diffusivity_dict = dict()

diffFile =  open("diffusivities_lattice.txt", "a")
diffFile.write(f"Computing diffusivity after {tdays} days of advection \n")
diffFile.close()
    
for exp_name, traj_ds in experiments.items():
    print(f"working on {exp_name}")
    
    if exp_name in ["coarse_Markov1_K15000", "coarsened_Markov1_K15000", "fine_advection", "coarsened_advection", "coarse_advection"]:
        obs = 720
    else:
        obs = tdays
    
    
    if "coarse_" in exp_name:
        theta = ds_field_coarse.THETA
    elif "coarsened_" in exp_name:
        theta = ds_field_coarsened.THETA
    else:
        theta = ds_field_fine.THETA
        
    diffusivity_dict[exp_name] = {}
    
    diffusivity_dict[exp_name]["full"] = dict(zip(["diffusivityAvg", "diffusivityArray", "diapycnalDist"], ACC_tools.compute_vertical_diffusivity(theta, traj_ds.where((traj_ds.isel(obs=0).lat > 60_000) \
                                                                                                     & (traj_ds.isel(obs=0).lat < 2_000_000 - 10_000) \
                                                                                                     & (traj_ds.z < - 50)).dropna('traj'), tdays, return_full=True, initTempWithSpline=True)))
    print(f"Diffusivity full: {diffusivity_dict[exp_name]['full']['diffusivityAvg']}")
    
    diffusivity_dict[exp_name]["top"] = dict(zip(["diffusivityAvg", "diffusivityArray", "diapycnalDist"], ACC_tools.compute_vertical_diffusivity(theta, traj_ds.where((traj_ds.isel(obs=0).lat > 60_000) \
                                                                                                    & (traj_ds.isel(obs=0).lat < 2_000_000 - 10_000) \
                                                                                                    & (traj_ds.isel(obs=0).z <= -200) \
                                                                                                    & (traj_ds.isel(obs=0).z > -600) \
                                                                                                    & (traj_ds.z < - 50)).dropna('traj'), tdays, return_full=True, initTempWithSpline=True)))
    print(f"Diffusivity top: {diffusivity_dict[exp_name]['top']['diffusivityAvg']}")
    
    diffusivity_dict[exp_name]["center"] = dict(zip(["diffusivityAvg", "diffusivityArray", "diapycnalDist"], ACC_tools.compute_vertical_diffusivity(theta, traj_ds.where((traj_ds.isel(obs=0).lat > 60_000) \
                                                                                                    & (traj_ds.isel(obs=0).lat < 2_000_000 - 10_000) \
                                                                                                    & (traj_ds.isel(obs=0).z <= -600) \
                                                                                                    & (traj_ds.isel(obs=0).z > -1200) \
                                                                                                    & (traj_ds.z < - 50)).dropna('traj'), tdays, return_full=True, initTempWithSpline=True)))
    print(f"Diffusivity center: {diffusivity_dict[exp_name]['center']['diffusivityAvg']}")
    
    diffusivity_dict[exp_name]["deep"] = dict(zip(["diffusivityAvg", "diffusivityArray", "diapycnalDist"], ACC_tools.compute_vertical_diffusivity(theta, traj_ds.where((traj_ds.isel(obs=0).lat > 60_000) \
                                                                                                    & (traj_ds.isel(obs=0).lat < 2_000_000 - 10_000) \
                                                                                                    & (traj_ds.isel(obs=0).z <= -1200) \
                                                                                                    & (traj_ds.isel(obs=0).z >= -1600) \
                                                                                                    & (traj_ds.z < - 50)).dropna('traj'), tdays, return_full=True, initTempWithSpline=True)))
    print(f"Diffusivity deep: {diffusivity_dict[exp_name]['deep']['diffusivityAvg']}")
    
    
    diffFile =  open("diffusivities_lattice.txt", "a")
    diffFile.write(f"{exp_name}, {diffusivity_dict[exp_name]['full']['diffusivityAvg']}, {diffusivity_dict[exp_name]['center']['diffusivityAvg']}, {diffusivity_dict[exp_name]['top']['diffusivityAvg']}, {diffusivity_dict[exp_name]['deep']['diffusivityAvg']} \n")
    diffFile.close()
    
    
    
with open("diffusivityDict_lattice.pickle", "wb") as pickFile:
    pickle.dump(diffusivity_dict, pickFile)