import numpy as np
import xarray as xr
import xgcm
import xmitgcm as xm
import MITgcmutils as mutil
from scipy.interpolate import InterpolatedUnivariateSpline
import sys
import cftime

# Writing utility
def write_binary_field(fname, data, dtype=np.dtype("float64")):
    dcopy = data.copy()
    print("Write to file: " + fname)
    os.path.exists(path)
    if dcopy.dtype != dtype:
        dcopy = dcopy.astype(dtype)
        print("changed dtype to " + str(dcopy))
    if sys.byteorder == "little":
        dcopy.byteswap(True)
    fid = open(fname, "wb")
    dcopy.tofile(fid)
    fid.close()
    
def read_binary_field(fname, shape, dtype=np.dtype("float64")):
    data = np.fromfile(fname, dtype)
    if sys.byteorder == "little":
        data.byteswap(True)
    return data.reshape(shape)

def mitgcm_stat_xr(statFile, depthFile, dt=None, t0=None):
    """Load MITgcm diagnostic statistics into xarray.
    
    Args:
        statFile (string): relative path to the statistics file
        depthFile (path) : relative path to binary `Depth` file which contains level depths
        dt (int)         : model timestep, for proper timestamping
        t0 (string)      : reference date for calendar. Format `yyyy-mm-dd hh:mm:ss`

    Returns:
        xarray.DataSet
    """
    
    statTypes = ['avg', 'std', 'min', 'max']
    
    ds = xr.Dataset()
    
    locs, totals, iters = mutil.diagnostics.readstats(statFile)
    RC = mutil.mds.rdmds(depthFile)[:, 0, 0]
    
    varNames = list(locs.keys())
    
    for varIdx, varName in enumerate(varNames):
        if locs[varName].shape[1] == 1:
            layered = False
        elif locs[varName].shape[1] == RC.shape[0]:
            layered = True
        else:
            raise RuntimeError("Vertical dimension must match that of layers or be unity.")
        
        for statIdx, statType in enumerate(statTypes):
            if layered:
                ds[f"{varName}_{statType}"] = xr.DataArray(locs['THETA'][:, :, statIdx], dims=('iterNum', 'RC'), coords=[iters['THETA'], RC])
            else: 
                ds[f"{varName}_{statType}"] = xr.DataArray(locs['THETA'][:, 0, statIdx], dims=('iterNum'), coords=[iters['THETA']])
                
    ds['volume'] = xr.DataArray(locs[varNames[0]][0, :, 4], dims=['RC'], coords=[RC])
    ds = ds.transpose()
    if dt and t0:
        ds['time'] = xr.DataArray(cftime.num2date(ds['iterNum'].values*dt, f"seconds since {t0}", calendar="360_day"), dims=['iterNum'], coords=[ds.iterNum])
        ds = ds.swap_dims({'iterNum': 'time'})
    return ds


def compute_vertical_diffusivity(EulerianField, LagrangianTrajectories, time, obs=None, interpMethod='linear', return_full=False, timestepDays=1, initTempWithSpline=False):
    """Compute a vertical diffusivity by localy interpolating z-coordinate of
    original isopycnal (here just the theta-level).

    Args:
        EulerianField (xarray.DataArray): Eulerian density/temperature field
        LagrangianTrajectories (xarray.DataSet): Particleset output from parcels
        time (int)           : time index to use in the Eulerian density/temperature field
        obs (int, optional)  : time/obs index to use in Particleset
        interpMethod (string): horizontal interpolation method for computing temperature/density profile 
        return_full (bool)   : if True, returns `diffusivityAvg, diffusivityArray, diapycnalDist`
        timestepDays (int)   : specifies how many days one unit of `time` is, for computing diffusivity
        initTempWithSpline (bool): if True, use spline interpolation to determine the initial temperatur

    Returns:
        int: average diffusivity (in m^2/s)

    """
    
    if not obs:
        obs = time
    
    if 'time' in EulerianField.coords:
        EulerianFieldStart = EulerianField.isel(time=0)
        EulerianFieldEnd = EulerianField.isel(time=time)
    else:
        EulerianFieldStart = EulerianFieldEnd = EulerianField
    
    if initTempWithSpline == True:
        theta_per_part_start = EulerianFieldStart.interp(XC=LagrangianTrajectories.isel(obs=0).lon, 
                                                         YC=LagrangianTrajectories.isel(obs=0).lat,
                                                         method=interpMethod)
    
    theta_per_part_end = EulerianFieldEnd.interp(XC=LagrangianTrajectories.isel(obs=obs).lon, 
                                                 YC=LagrangianTrajectories.isel(obs=obs).lat,
                                                 method=interpMethod)
    
    diapycnalDist = np.zeros(theta_per_part_end.traj.shape[0])
    z_coords = -theta_per_part_end.Z.data # Z normally runs negative but scipy requires monotonic increase

    for particleIdx in range(theta_per_part_end.traj.shape[0]):
        if initTempWithSpline == True:
            theta_orig = InterpolatedUnivariateSpline(z_coords, theta_per_part_start.isel(traj=particleIdx))(-LagrangianTrajectories.isel(obs=0, traj=particleIdx).z.data)
        else:
            theta_orig = LagrangianTrajectories.isel(traj=particleIdx, obs=0).pot_temp
        
        isopycnalZ = InterpolatedUnivariateSpline(z_coords, theta_per_part_end.isel(traj=particleIdx) - theta_orig).roots()
        currZ = LagrangianTrajectories.isel(traj=particleIdx, obs=obs).z.data
        
        if isopycnalZ.shape[0] == 1:
            diapycnalDist[particleIdx] = isopycnalZ[0] + currZ
        elif isopycnalZ.shape[0] > 0:
            diapycnalDist[particleIdx] = (isopycnalZ + currZ).min()
        else:
            diapycnalDist[particleIdx] = np.nan
    
    diffusivityArray = 0.5 * diapycnalDist**2 / (time * timestepDays * 24 * 60 * 60)
    diffusivityAvg = np.nanmean(diffusivityArray)
    
    if return_full:
        return diffusivityAvg, diffusivityArray, diapycnalDist
    
    else:
        return diffusivityAvg


