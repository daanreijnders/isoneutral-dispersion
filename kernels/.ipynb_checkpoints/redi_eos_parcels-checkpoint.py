#!/usr/bin/env python3

import math
from datetime import timedelta as delta
from glob import glob
from os import path

import numpy as np

import parcels

pclass = {'scipy': parcels.ScipyParticle,
          'jit': parcels.JITParticle}


def create_fieldset():
    workdir = path.realpath(
        "/projects/0/topios/hydrodynamic_data/MITgcm/ACC/ACC_30d/")

    variables = {"U": "UVEL", "V": "VVEL", "W": "WVEL", "T": "THETA", "S": "SALT"}
    dimensions = {
        "U": {"lon": "XG", "lat": "YG", "depth": "Zl", "time": "time"},
        "V": {"lon": "XG", "lat": "YG", "depth": "Zl", "time": "time"},
        "W": {"lon": "XG", "lat": "YG", "depth": "Zl", "time": "time"},
        "S": {"lon": "XC", "lat": "YC", "depth": "Z", "time": "time"},
        "T": {"lon": "XC", "lat": "YC", "depth": "Z", "time": "time"},
    }

    fieldset = parcels.FieldSet.from_mitgcm(
        [workdir + "/ACC_30d_diagnostics_with_salinity.nc"],
        variables,
        dimensions,
        mesh="flat",
    )
    fieldset.add_constant("domain_width", 1000000)
    fieldset.add_constant("dres_horiz", 50)
    fieldset.add_constant("dres_vert", 0.1)
    fieldset.S.interp_method = "linear"
    fieldset.T.interp_method = "linear"

    fieldset.add_periodic_halo(zonal=True)
    return fieldset


def AdvectionDiffusionEM_EOSredi_3D(particle, fieldset, time):
#     Do: - minimum horizontal diffusion
#         - mixed layer parameterization
#         - boundary conditions
    sqrt2 = 1.41421356237
#     Wiener increment with zero mean and std of sqrt(dt)
    dWx = random.uniform(-1.0, 1.0) * math.sqrt(math.fabs(particle.dt) * 3)
    dWy = random.uniform(-1.0, 1.0) * math.sqrt(math.fabs(particle.dt) * 3)
    dWz = random.uniform(-1.0, 1.0) * math.sqrt(math.fabs(particle.dt) * 3)
    
#     Compute isopycnal / neutral surface slopes
##    equation of state
##    In linear EOS, this is simple. In other case, use TEOS
#     rho0 = 1028
#     tRef = 10
#     sRef = 35.0000
#     alpha=1.7e-4
#     beta=0 # (In our case. Normal representative value: 7.6e-4)
    
    isopycK = 200.
    
    T = fieldset.T[time, particle.depth, particle.lat, particle.lon]
#     Gradient calculations. Note that normally, we should compute density at the +1/-1 points
#     through T and S (using an EOS). Then compute the gradient in rho through this.
    T_xm1 = fieldset.T[time, particle.depth, particle.lat, particle.lon - fieldset.dres_horiz]
    T_xp1 = fieldset.T[time, particle.depth, particle.lat, particle.lon + fieldset.dres_horiz]
#     S_xm1 = fieldset.S[time, particle.depth, particle.lat, particle.lon - fieldset.dres_horiz]
#     S_xp1 = fieldset.S[time, particle.depth, particle.lat, particle.lon + fieldset.dres_horiz]
    dRhodx = (T_xp1 - T_xm1) / (2 * fieldset.dres_horiz)
    
    T_ym1 = fieldset.T[time, particle.depth, particle.lat - fieldset.dres_horiz, particle.lon]
    T_yp1 = fieldset.T[time, particle.depth, particle.lat + fieldset.dres_horiz, particle.lon]
#     S_ym1 = fieldset.S[time, particle.depth, particle.lat - fieldset.dres_horiz, particle.lon]
#     S_yp1 = fieldset.S[time, particle.depth, particle.lat + fieldset.dres_horiz, particle.lon]
    dRhody = (T_yp1 - T_ym1) / (2 * fieldset.dres_horiz)

    T_zm1 = fieldset.T[time, particle.depth - fieldset.dres_vert, particle.lat, particle.lon]
    T_zp1 = fieldset.T[time, particle.depth + fieldset.dres_vert, particle.lat, particle.lon]
#     S_ym1 = fieldset.S[time, particle.depth, particle.lat - fieldset.dres_horiz, particle.lon]
#     S_yp1 = fieldset.S[time, particle.depth, particle.lat + fieldset.dres_horiz, particle.lon]
    dRhodz = (T_zp1 - T_zm1) / (2 * fieldset.dres_vert)
    
    Sx = -dRhodx/dRhodz
    Sy = -dRhody/dRhodz
    # Rudimentary slope clipping
    # look into more sophisticated options later
    maxSlope=1e-02
    if math.fabs(Sx) > maxSlope:
        Sx = math.copysign(maxSlope, Sx)
    if math.fabs(Sy) > maxSlope:
        Sy = math.copysign(maxSlope, Sy)
    Sabs2 = math.sqrt(Sx**2 + Sy**2)
    
#     Note: dKxx and dKyy are zero (except in case of variable diffusivity ~ Visbeck).
#           dKxz, dKyz, and dKzz are often zero 
#           (when [-2*dres, 2*dres] falls between two density gridpoints)
#           Computing these implies computing slopes for all adjacent points (6)
#           for which we need to know the density at even more neighboring points 
#           (6 * 6 = 36, but without overlap 25 total interpolations vs 6)
#           
    
#     see comment about variable diffusivity
    dk11dx = 0
    dk22dy = 0
#     Gradients of diffusivity tensor elements. This is currently very expensive; 18 more interpolations are needed (for just density)
    T_xp2 = fieldset.T[time, particle.depth, particle.lat, particle.lon + 2 * fieldset.dres_horiz]
    dRhodx_xp1 = (T_xp2 - T_xm1) / (2 * fieldset.dres_horiz)
    
    T_xp1_ym1 = fieldset.T[time, particle.depth, particle.lat - fieldset.dres_horiz, particle.lon + fieldset.dres_horiz]
    T_xp1_yp1 = fieldset.T[time, particle.depth, particle.lat + fieldset.dres_horiz, particle.lon + fieldset.dres_horiz]
    dRhody_xp1 = (T_xp1_yp1 - T_xp1_ym1) / (2 * fieldset.dres_horiz)
    
    T_xp1_zm1 = fieldset.T[time, particle.depth - fieldset.dres_vert, particle.lat, particle.lon + fieldset.dres_horiz]
    T_xp1_zp1 = fieldset.T[time, particle.depth + fieldset.dres_vert, particle.lat, particle.lon + fieldset.dres_horiz]
    dRhodz_xp1 = (T_xp1_zp1 - T_xp1_zm1) / (2 * fieldset.dres_vert)
    
    T_xm2 = fieldset.T[time, particle.depth, particle.lat, particle.lon - 2 * fieldset.dres_horiz]
    dRhodx_xm1 = (T - T_xm2) / (2 * fieldset.dres_horiz)
    
    T_xm1_ym1 = fieldset.T[time, particle.depth, particle.lat - fieldset.dres_horiz, particle.lon - fieldset.dres_horiz]
    T_xm1_yp1 = fieldset.T[time, particle.depth, particle.lat + fieldset.dres_horiz, particle.lon - fieldset.dres_horiz]
    dRhody_xm1 = (T_xm1_yp1 - T_xm1_ym1) / (2 * fieldset.dres_horiz)
    
    T_xm1_zm1 = fieldset.T[time, particle.depth - fieldset.dres_vert, particle.lat, particle.lon - fieldset.dres_horiz]
    T_xm1_zp1 = fieldset.T[time, particle.depth + fieldset.dres_vert, particle.lat, particle.lon - fieldset.dres_horiz]
    dRhodz_xm1 = (T_xm1_zp1 - T_xm1_zm1) / (2 * fieldset.dres_vert)
    
    dRhodx_yp1 = (T_xp1_yp1 - T_xm1_yp1) / (2 * fieldset.dres_horiz)
    
    T_yp2 = fieldset.T[time, particle.depth, particle.lat + 2 * fieldset.dres_horiz, particle.lon]
    dRhody_yp1 = (T_yp2 - T) / (2 * fieldset.dres_horiz)
    
    T_yp1_zm1 = fieldset.T[time, particle.depth - fieldset.dres_vert, particle.lat + fieldset.dres_horiz, particle.lon]
    T_yp1_zp1 = fieldset.T[time, particle.depth + fieldset.dres_vert, particle.lat + fieldset.dres_horiz, particle.lon]
    dRhodz_yp1 = (T_yp1_zp1 - T_yp1_zm1) / (2 * fieldset.dres_vert)
    
    dRhodx_ym1 = (T_xp1_ym1 - T_xm1_ym1) / (2 * fieldset.dres_horiz)
    
    T_ym2 = fieldset.T[time, particle.depth, particle.lat - 2 * fieldset.dres_horiz, particle.lon]
    dRhody_ym1 = (T - T_ym2) / (2 * fieldset.dres_horiz)
    
    T_ym1_zm1 = fieldset.T[time, particle.depth - fieldset.dres_vert, particle.lat - fieldset.dres_horiz, particle.lon]
    T_ym1_zp1 = fieldset.T[time, particle.depth + fieldset.dres_vert, particle.lat - fieldset.dres_horiz, particle.lon]
    dRhodz_ym1 = (T_ym1_zp1 - T_ym1_zm1) / (2 * fieldset.dres_vert)
    
    dRhodx_zm1 = (T_xp1_zm1 - T_xm1_zm1) / (2 * fieldset.dres_vert)
    dRhodx_zp1 = (T_xp1_zp1 - T_xm1_zp1) / (2 * fieldset.dres_vert)

    dRhody_zm1 = (T_yp1_zm1 - T_ym1_zm1) / (2 * fieldset.dres_vert)
    dRhody_zp1 = (T_yp1_zp1 - T_ym1_zp1) / (2 * fieldset.dres_vert)
    
    T_zp2 = fieldset.T[time, particle.depth + 2 * fieldset.dres_vert, particle.lat, particle.lon]
    dRhodz_zp1 = (T_zp2 - T) / (2* fieldset.dres_vert)
    
    T_zm2 = fieldset.T[time, particle.depth - 2 * fieldset.dres_vert, particle.lat, particle.lon]
    dRhodz_zm1 = (T - T_zm2) / (2* fieldset.dres_vert)

    Sx_xp1 = -dRhodx_xp1/dRhodz_xp1
    Sy_xp1 = -dRhody_xp1/dRhodz_xp1
    if math.fabs(Sx_xp1) > maxSlope:
        Sx_xp1 = math.copysign(maxSlope, Sx_xp1)
    if math.fabs(Sy_xp1) > maxSlope:
        Sy_xp1 = math.copysign(maxSlope, Sy_xp1)
#     Sabs2_xp1 = math.sqrt(Sx_xp1**2 + Sy_xp1**2)
    
    Sx_xm1 = -dRhodx_xm1/dRhodz_xm1
    Sy_xm1 = -dRhody_xm1/dRhodz_xm1
    if math.fabs(Sx_xm1) > maxSlope:
        Sx_xm1 = math.copysign(maxSlope, Sx_xm1)
    if math.fabs(Sy_xm1) > maxSlope:
        Sy_xm1 = math.copysign(maxSlope, Sy_xm1)
#     Sabs2_xm1 = math.sqrt(Sx_xm1**2 + Sy_xm1**2)
    
    Sx_yp1 = -dRhodx_yp1/dRhodz_yp1
    Sy_yp1 = -dRhody_yp1/dRhodz_yp1
    if math.fabs(Sx_yp1) > maxSlope:
        Sx_yp1 = math.copysign(maxSlope, Sx_yp1)
    if math.fabs(Sy_yp1) > maxSlope:
        Sy_yp1 = math.copysign(maxSlope, Sy_yp1)
#     Sabs2_yp1 = math.sqrt(Sx_yp1**2 + Sy_yp1**2)
    
    Sx_ym1 = -dRhodx_ym1/dRhodz_ym1
    Sy_ym1 = -dRhody_ym1/dRhodz_ym1
    if math.fabs(Sx_ym1) > maxSlope:
        Sx_ym1 = math.copysign(maxSlope, Sx_ym1)
    if math.fabs(Sy_ym1) > maxSlope:
        Sy_ym1 = math.copysign(maxSlope, Sy_ym1)
#     Sabs2_ym1 = math.sqrt(Sx_ym1**2 + Sy_ym1**2)
    
    Sx_zp1 = -dRhodx_zp1/dRhodz_zp1
    Sy_zp1 = -dRhody_zp1/dRhodz_zp1
    if math.fabs(Sx_zp1) > maxSlope:
        Sx_zp1 = math.copysign(maxSlope, Sx_zp1)
    if math.fabs(Sy_zp1) > maxSlope:
        Sy_zp1 = math.copysign(maxSlope, Sy_zp1)
    Sabs2_zp1 = math.sqrt(Sx_zp1**2 + Sy_zp1**2)
    
    Sx_zm1 = -dRhodx_zm1/dRhodz_zm1
    Sy_zm1 = -dRhody_zm1/dRhodz_zm1
    if math.fabs(Sx_zm1) > maxSlope:
        Sx_zm1 = math.copysign(maxSlope, Sx_zm1)
    if math.fabs(Sy_zm1) > maxSlope:
        Sy_zm1 = math.copysign(maxSlope, Sy_zm1)
    Sabs2_zm1 = math.sqrt(Sx_zm1**2 + Sy_zm1**2)
    
    dk31dx = isopycK * (Sx_xp1 - Sx_xm1) / (2 * fieldset.dres_horiz)
    dk32dy = isopycK * (Sy_yp1 - Sy_ym1) / (2 * fieldset.dres_horiz)
    dk33dz = isopycK * (Sabs2_zp1 - Sabs2_zm1) / (2 * fieldset.dres_vert)
    
#     Filling out the Redi tensor
    k11 = isopycK
#     k12 = 0
#     k13 = isopycK * Sx
#     k21 = 0
    k22 = isopycK
#     k23 = isopycK * Sy
    k31 = isopycK * Sx
    k32 = isopycK * Sy
    k33 = isopycK * Sabs2
    
#     Cholesky decomposed tensor elements
    sigma11 = math.sqrt(k11)
#     sigma12 = 0
#     sigma13 = sigma31
#     sigma21 = 0
    sigma22 = math.sqrt(k22)
#     sigma23 = sigma32
    sigma31 = k31 / sigma11
    sigma32 = k32 / sigma22
    sigma33 = math.sqrt(k33 - sigma31 ** 2 - sigma32 ** 2)

    u, v, w = fieldset.UVW[time, particle.depth, particle.lat, particle.lon]

    ax = u #+ dk11dx
    ay = v #+ dk22dy
    az = w #+ dk31dx + dk32dy + dk33dz

#     Particle positions are updated only after evaluating all terms.
#     Note: sqrt2 comes from 2 * K = sigma * sigma.T
    particle.lon += ax * particle.dt + sqrt2 * sigma11 * dWx
    particle.lat += ay * particle.dt + sqrt2 * sigma22 * dWy
    particle.depth += (
        az * particle.dt
        + sqrt2 * sigma31 * dWx
        + sqrt2 * sigma32 * dWy
        + sqrt2 * sigma33 * dWz
    )
    
#     for debugging
#     if time < 1800:
#         print('-------')
#         print("u * particle.dt", "{0:.6g}".format(u * particle.dt))
#         print("v * particle.dt", "{0:.6g}".format(v * particle.dt))
#         print("w * particle.dt", "{0:.6g}".format(w * particle.dt))
#         print("(dk31dx + dk32dy + dk33dz) * particle.dt", "{0:.6g}".format((dk31dx + dk32dy + dk33dz) * particle.dt))
#         print("sqrt2 * sigma11 * dWx", "{0:.6g}".format(sqrt2 * sigma11 * dWx))
#         print("sqrt2 * sigma22 * dWy", "{0:.6g}".format(sqrt2 * sigma22 * dWy))
#         print("sqrt2 * sigma31 * dWx", "{0:.6g}".format(sqrt2 * sigma31 * dWx))
#         print("sqrt2 * sigma32 * dWy", "{0:.6g}".format(sqrt2 * sigma32 * dWy))
#         print("sqrt2 * sigma33 * dWz", "{0:.6g}".format(sqrt2 * sigma33 * dWz))


def initialize(fieldset, nparts=1000, z_init=-1000, partType='jit', run="Redi", length=""):
    lons = np.ones(nparts) * fieldset.U.grid.lon[100]
    lats = np.ones(nparts) * fieldset.V.grid.lat[200]
    depth = z_init * np.ones(nparts)

    if run == "control":
        lons = lons[0]
        lats = lats[0]
        depth = depth[0]
    pset = parcels.ParticleSet.from_list(
        fieldset=fieldset, pclass=pclass[partType], lon=lons, lat=lats, depth=depth
    )

    outfile = pset.ParticleFile(
        name=f"MIT_test_{run}_{partType}_n{nparts}_d{length}.nc", outputdt=delta(hours=1))
    return pset, outfile


def periodicBC(particle, fieldset, time):
    if particle.lon < 0:
        particle.lon += fieldset.domain_width
    elif particle.lon > fieldset.domain_width:
        particle.lon -= fieldset.domain_width


if __name__ == "__main__":
    fieldset = create_fieldset()
    for run in ["RediEOS", "control"]:
        lenght = 10
        pset, outfile = initialize(fieldset, run=run, length=length)
        if run == "control":
            AD_kernel = parcels.kernels.AdvectionRK4_3D
        elif run == 'RediEOS':
            AD_kernel = AdvectionDiffusionEM_RediEOS_3D
        else:
            raise("Expected `run` to be either `Redi` or `control`.")
        try:
            pset.execute(
                pset.Kernel(AD_kernel) + pset.Kernel(periodicBC),
                # kernel_veloprint_3D_data,
                runtime=delta(days=length),
                #    runtime=delta(minutes=5),
                dt=delta(minutes=5),
                output_file=outfile,
            )
            pset.close()
        except (NameError, TypeError) as e:
            pass
