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

    filenames = {
        "U": workdir + "/ACC_30d_diagnostics.nc",
        "V": workdir + "/ACC_30d_diagnostics.nc",
        "W": workdir + "/ACC_30d_diagnostics.nc",
    }
    variables = {
        "U": "UVEL",
        "V": "VVEL",
        "W": "WVEL",
        "K11": "GM_Kux",
        "K22": "GM_Kvy",
        #     "K13": "GM_Kuz",
        #     "K23": "GM_Kvz",
        "K31": "GM_Kwx",
        "K32": "GM_Kwy",
        "K33": "GM_Kwz",
    }
    dimensions = {
        "U": {"lon": "XG", "lat": "YG", "depth": "Zl", "time": "time"},
        "V": {"lon": "XG", "lat": "YG", "depth": "Zl", "time": "time"},
        "W": {"lon": "XG", "lat": "YG", "depth": "Zl", "time": "time"},
        "K11": {"lon": "XG", "lat": "YC", "depth": "Z", "time": "time"},
        "K22": {"lon": "XC", "lat": "YG", "depth": "Z", "time": "time"},
        #     "K13": {"lon": "XG", "lat": "YC", "depth": "Z", "time": "time"},
        #     "K23": {"lon": "XC", "lat": "YG", "depth": "Z", "time": "time"},
        "K31": {"lon": "XC", "lat": "YC", "depth": "Zl", "time": "time"},
        "K32": {"lon": "XC", "lat": "YC", "depth": "Zl", "time": "time"},
        "K33": {"lon": "XC", "lat": "YC", "depth": "Zl", "time": "time"},
    }

    fieldset = parcels.FieldSet.from_mitgcm(
        [workdir + "/ACC_30d_diagnostics.nc"], variables, dimensions, mesh="flat"
    )
    fieldset.add_constant("domain_width", 1000000)
    fieldset.add_constant("dres_horiz", 500)
    fieldset.add_constant("dres_vert", 0.1)

    rediElemns = [
        "K11",
        "K22",
        #     "K13",
        #     "K23",
        "K31",
        "K32",
        "K33",
    ]
    for elemn in rediElemns:
        getattr(fieldset, elemn).interp_method = "linear"

    fieldset.add_periodic_halo(zonal=True)
    return fieldset


def AdvectionDiffusionEM_RediInterp_3D(particle, fieldset, time):
    # Wiener increment with zero mean and std of sqrt(dt)
    sqrt2 = 1.41421356237

    dWx = random.uniform(-1., 1.) * math.sqrt(math.fabs(particle.dt) * 3)
    dWy = random.uniform(-1., 1.) * math.sqrt(math.fabs(particle.dt) * 3)
    dWz = random.uniform(-1., 1.) * math.sqrt(math.fabs(particle.dt) * 3)

    k11 = fieldset.K11[time, particle.depth, particle.lat, particle.lon]
    k11xp1 = fieldset.K11[time, particle.depth,
                          particle.lat, particle.lon + fieldset.dres_horiz]
    k11xm1 = fieldset.K11[time, particle.depth,
                          particle.lat, particle.lon - fieldset.dres_horiz]
    dk11dx = (k11xp1 - k11xm1) / (2 * fieldset.dres_horiz)
    # ^Zero bc isopycnal diffusion (except in case of Visbeck)

    # K12 = 0
    # K21 = 0

    k22 = fieldset.K22[time, particle.depth, particle.lat, particle.lon]
    k22yp1 = fieldset.K22[time, particle.depth,
                          particle.lat + fieldset.dres_horiz, particle.lon]
    k22ym1 = fieldset.K22[time, particle.depth,
                          particle.lat - fieldset.dres_horiz, particle.lon]
    dk22dy = (k22yp1 - k22ym1) / (2 * fieldset.dres_horiz)

    # K13 & K23 are (in a continuous field) equal to K31 and K32

    k31 = fieldset.K31[time, particle.depth, particle.lat, particle.lon]
    k31xp1 = fieldset.K31[time, particle.depth,
                          particle.lat, particle.lon + fieldset.dres_horiz]
    k31xm1 = fieldset.K31[time, particle.depth,
                          particle.lat, particle.lon - fieldset.dres_horiz]
    dk31dx = (k31xp1 - k31xm1) / (2 * fieldset.dres_horiz)

    k32 = fieldset.K32[time, particle.depth, particle.lat, particle.lon]
    k32yp1 = fieldset.K32[time, particle.depth,
                          particle.lat + fieldset.dres_horiz, particle.lon]
    k32ym1 = fieldset.K32[time, particle.depth,
                          particle.lat - fieldset.dres_horiz, particle.lon]
    dk32dy = (k32yp1 - k32ym1) / (2 * fieldset.dres_horiz)

    k33 = fieldset.K33[time, particle.depth, particle.lat, particle.lon]
    k33zp1 = fieldset.K33[time, particle.depth
                          + fieldset.dres_vert, particle.lat, particle.lon]
    k33zm1 = fieldset.K33[time, particle.depth
                          - fieldset.dres_vert, particle.lat, particle.lon]
    dk33dz = (k33zp1 - k33zm1) / (2 * fieldset.dres_horiz)

    sigma11 = math.sqrt(k11)
    # sigma21 = 0
    sigma22 = math.sqrt(k22)
    sigma31 = k31 / sigma11
    sigma32 = k32 / sigma22
    sigma33 = math.sqrt(k33 - sigma31**2 - sigma32**2)

    u, v, w = fieldset.UVW[time, particle.depth, particle.lat, particle.lon]

    ax = u + dk11dx
    ay = v + dk22dy
    az = w + dk31dx + dk32dy + dk33dz

    # Particle positions are updated only after evaluating all terms.
    particle.lon += ax * particle.dt + sqrt2 * sigma11 * dWx
    particle.lat += ay * particle.dt + sqrt2 * sigma22 * dWy
    particle.depth += az * particle.dt + sqrt2 * sigma31 * \
        dWx + sqrt2 * sigma32 * dWy + sqrt2 * sigma33 * dWz

#     if time <= 1800:
#         print("-------------")
#         print("Time: ", time)
#         print("U * particle.dt: ", U * particle.dt)
#         print("V * particle.dt: ", V * particle.dt)
#         print("W * particle.dt: ", W * particle.dt)
#         print("ax * particle.dt: ", ax * particle.dt)
#         print("ay * particle.dt: ", ay * particle.dt)
#         print("az * particle.dt: ", az * particle.dt)
#         print("bx [sqrt2 * sigma11 * dWx]", sqrt2 * sigma11 * dWx)
#         print("by [sqrt2 * sigma22 * dWy]", sqrt2 * sigma22 * dWy)
#         print("bz [sqrt2 * sigma31 * dWx + sqrt2 * sigma32 * dWy + sqrt2 * sigma33 * dWz]", sqrt2 * sigma31 * dWx + sqrt2 * sigma32 * dWy + sqrt2 * sigma33 * dWz)


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
    for run in ["RediInterp", "control"]:
        lenght = 10
        pset, outfile = initialize(fieldset, run=run, length=length)
        if run == "control":
            AD_kernel = parcels.kernels.AdvectionRK4_3D
        elif run == 'RediInterp':
            AD_kernel = AdvectionDiffusionEM_RediInterp_3D
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
