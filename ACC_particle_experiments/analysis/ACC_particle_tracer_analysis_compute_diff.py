import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
import sys
sys.path.append('/nethome/4302001/diffusion-hydro-mod/tools/')
import ACC_tools

trajdir = '/data/oceanparcels/output_data/data_Daan/project_diffusion/trelease/'
tdays = 180

experiment_files = {
    "fine_advection" : "trelease_fine_locb_advection_None_K15000_p100000_dt40m.nc",
    "coarsened_advection" : "trelease_coarsened_locb_advection_None_K15000_p100000_dt40m.nc",
    "coarsened_EM_Le_Sommer_C0.3" : "trelease_coarsened_locb_EM_Le_Sommer_C0.3_p100000_dt40m.nc",
    "coarsened_EM_Le_Sommer_C1" : "trelease_coarsened_locb_EM_Le_Sommer_C1_p100000_dt40m.nc",
    "coarsened_EM_Redi_smallslope_K1500" : "trelease_coarsened_locb_EM_Redi_smallslope_K1500_p100000_dt40m.nc",
    "coarsened_EM_Redi_smallslope_K5000" : "trelease_coarsened_locb_EM_Redi_smallslope_K5000_p100000_dt40m.nc",
    "coarsened_EM_Redi_smallslope_K15000" : "trelease_coarsened_locb_EM_Redi_smallslope_K15000_p100000_dt40m.nc",
    "coarse_advection" : "trelease_coarse_locb_advection_None_K15000_p100000_dt40m.nc",
    "coarse_EM_Le_Sommer_C0.3": "trelease_coarse_locb_EM_Le_Sommer_C0.3_p100000_dt40m.nc",
    "coarse_EM_Le_Sommer_C1":   "trelease_coarse_locb_EM_Le_Sommer_C1_p100000_dt40m.nc",
    "coarse_EM_Le_Sommer_C3":   "trelease_coarse_locb_EM_Le_Sommer_C3_p100000_dt40m.nc",
    "coarse_EM_Redi_smallslope_K1500": "trelease_coarse_locb_EM_Redi_smallslope_K1500_p100000_dt40m.nc",
    "coarse_EM_Redi_smallslope_K5000": "trelease_coarse_locb_EM_Redi_smallslope_K5000_p100000_dt40m.nc",
    "coarse_EM_Redi_smallslope_K15000": "trelease_coarse_locb_EM_Redi_smallslope_K15000_p100000_dt40m.nc",
    "coarse_M1_Redi_smallslope_K1500": "trelease_coarse_locb_M1_Redi_smallslope_K1500_p100000_dt40m.nc",
    "coarse_M1_Redi_smallslope_K5000": "trelease_coarse_locb_M1_Redi_smallslope_K5000_p100000_dt40m.nc",
    "coarse_M1_Redi_smallslope_K15000": "trelease_coarse_locb_M1_Redi_smallslope_K15000_p100000_dt40m.nc",
}

experiments = dict()
for experiment, filename in experiment_files.items():
    experiments[experiment] = xr.open_dataset(trajdir + str(filename), decode_times=False)
    
ds_field_fine = xr.open_dataset('/data/oceanparcels/input_data/MITgcm/ACC_channel/ACC_ridge_fine_1y_locb.nc')
ds_field_coarsened = xr.open_dataset("/data/oceanparcels/input_data/MITgcm/ACC_channel/ACC_ridge_fine_1y_locb_coarsened.nc")
ds_field_coarse = xr.open_dataset('/data/oceanparcels/input_data/MITgcm/ACC_channel/ACC_ridge_coarse_1y_locb_with_derivative.nc')
    
diffusivity_dict = dict()

diffFile =  open("diffusivities.txt", "a")
diffFile.write(f"Computing diffusivity after {tdays} days of advection")
diffFile.close()
    
for exp_name, exp_ds in experiments.items():
    print(f"working on {exp_name}")
    if "coarse_" in exp_name:
        diffusivity_dict[exp_name] = ACC_tools.compute_vertical_diffusivity(ds_field_coarse.THETA, exp_ds, tdays)
    elif "coarsened" in exp_name:
         diffusivity_dict[exp_name] = ACC_tools.compute_vertical_diffusivity(ds_field_coarsened.THETA, exp_ds, tdays)
    elif "fine" in exp_name:
         diffusivity_dict[exp_name] = ACC_tools.compute_vertical_diffusivity(ds_field_fine.THETA, exp_ds, tdays)
    print(f"Diffusivity: {diffusivity_dict[exp_name]}")
    diffFile =  open("diffusivities_trelease.txt", "a")
    diffFile.write(f"{exp_name}: {diffusivity_dict[exp_name]} \n")
    diffFile.close()