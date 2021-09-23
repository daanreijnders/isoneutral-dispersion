import math
import parcels.rng as ParcelsRandom


def K_elemns_EOSslope(particle, fieldset, time):
# fmt: off
    #     Do: - minimum horizontal diffusion
    #         - mixed layer parameterization
    #         - boundary conditions

    #     Compute isopycnal / neutral surface slopes
    ##    equation of state
    ##    In linear EOS, this is simple. In other case, use TEOS
    #     rho0 = 1028
    #     tRef = 10
    #     sRef = 35.0000
    #     alpha=1.7e-4
    #     beta=0 # (In our case. Normal representative value: 7.6e-4)

    rho = fieldset.Rho[time, particle.depth, particle.lat, particle.lon]
    #     Gradient calcuparticle.lations. Note that normally, we should compute density at the +1/-1 x,y,z
    #     points through T and S (using an EOS). Then compute the gradient in rho through this.
    rho_xm1 = fieldset.Rho[time, particle.depth, particle.lat, particle.lon - fieldset.dsamp_xy]
    rho_xp1 = fieldset.Rho[time, particle.depth, particle.lat, particle.lon + fieldset.dsamp_xy]
    dRhodx = (rho_xp1 - rho_xm1) / (2 * fieldset.dsamp_xy)

    rho_ym1 = fieldset.Rho[time, particle.depth, particle.lat - fieldset.dsamp_xy, particle.lon]
    rho_yp1 = fieldset.Rho[time, particle.depth, particle.lat + fieldset.dsamp_xy, particle.lon]
    dRhody = (rho_yp1 - rho_ym1) / (2 * fieldset.dsamp_xy)

    rho_zm1 = fieldset.Rho[time, particle.depth - fieldset.dsamp_z, particle.lat, particle.lon]
    rho_zp1 = fieldset.Rho[time, particle.depth + fieldset.dsamp_z, particle.lat, particle.lon]
    dRhodz = (rho_zp1 - rho_zm1) / (2 * fieldset.dsamp_z)

    Sx = -dRhodx / dRhodz
    Sy = -dRhody / dRhodz
    
    # Rudimentary slope clipping
    # look into more sophisticated options particle.later
#     maxSlope = 1e-02
#     if math.fabs(Sx) > maxSlope:
#         Sx = math.copysign(maxSlope, Sx)
#     if math.fabs(Sy) > maxSlope:
#         Sy = math.copysign(maxSlope, Sy)
    Sabs2 = Sx ** 2 + Sy ** 2

    # Computing gradients in K requires knowing the slope at adjacent points
    #     Note: dKxx and dKyy are zero (except in case of variable diffusivity ~ Visbeck).
    #           dKzx, dKzy, and dKzz are often zero
    #           (when [-2*dres, 2*dres] falls between two density gridpoints)
    #           Computing these implies computing slopes for all adjacent points (6)
    #           for which we need to know the density at even more neighboring points
    #           (6 * 6 = 36, but without overlap 25 total interpoparticle.lations vs 6)
    rho_xp2 = fieldset.Rho[time, particle.depth, particle.lat, particle.lon + 2 * fieldset.dsamp_xy]
    dRhodx_xp1 = (rho_xp2 - rho) / (2 * fieldset.dsamp_xy)

    rho_xp1_ym1 = fieldset.Rho[time, particle.depth, particle.lat - fieldset.dsamp_xy, particle.lon + fieldset.dsamp_xy,]
    rho_xp1_yp1 = fieldset.Rho[time, particle.depth, particle.lat + fieldset.dsamp_xy, particle.lon + fieldset.dsamp_xy,]
    dRhody_xp1 = (rho_xp1_yp1 - rho_xp1_ym1) / (2 * fieldset.dsamp_xy)

    rho_xp1_zm1 = fieldset.Rho[time, particle.depth - fieldset.dsamp_z, particle.lat, particle.lon + fieldset.dsamp_xy,]
    rho_xp1_zp1 = fieldset.Rho[time, particle.depth + fieldset.dsamp_z, particle.lat, particle.lon + fieldset.dsamp_xy,]
    dRhodz_xp1 = (rho_xp1_zp1 - rho_xp1_zm1) / (2 * fieldset.dsamp_z)

    rho_xm2 = fieldset.Rho[time, particle.depth, particle.lat, particle.lon - 2 * fieldset.dsamp_xy]
    dRhodx_xm1 = (rho - rho_xm2) / (2 * fieldset.dsamp_xy)

    rho_xm1_ym1 = fieldset.Rho[time, particle.depth, particle.lat - fieldset.dsamp_xy, particle.lon - fieldset.dsamp_xy,]
    rho_xm1_yp1 = fieldset.Rho[time, particle.depth, particle.lat + fieldset.dsamp_xy, particle.lon - fieldset.dsamp_xy,]
    dRhody_xm1 = (rho_xm1_yp1 - rho_xm1_ym1) / (2 * fieldset.dsamp_xy)

    rho_xm1_zm1 = fieldset.Rho[time, particle.depth - fieldset.dsamp_z, particle.lat, particle.lon - fieldset.dsamp_xy,]
    rho_xm1_zp1 = fieldset.Rho[time, particle.depth + fieldset.dsamp_z, particle.lat, particle.lon - fieldset.dsamp_xy,]
    dRhodz_xm1 = (rho_xm1_zp1 - rho_xm1_zm1) / (2 * fieldset.dsamp_z)
    dRhodx_yp1 = (rho_xp1_yp1 - rho_xm1_yp1) / (2 * fieldset.dsamp_xy)

    rho_yp2 = fieldset.Rho[time, particle.depth, particle.lat + 2 * fieldset.dsamp_xy, particle.lon]
    dRhody_yp1 = (rho_yp2 - rho) / (2 * fieldset.dsamp_xy)

    rho_yp1_zm1 = fieldset.Rho[time, particle.depth - fieldset.dsamp_z, particle.lat + fieldset.dsamp_xy, particle.lon,]
    rho_yp1_zp1 = fieldset.Rho[time, particle.depth + fieldset.dsamp_z, particle.lat + fieldset.dsamp_xy, particle.lon,]
    dRhodz_yp1 = (rho_yp1_zp1 - rho_yp1_zm1) / (2 * fieldset.dsamp_z)
    dRhodx_ym1 = (rho_xp1_ym1 - rho_xm1_ym1) / (2 * fieldset.dsamp_xy)

    rho_ym2 = fieldset.Rho[time, particle.depth, particle.lat - 2 * fieldset.dsamp_xy, particle.lon]
    dRhody_ym1 = (rho - rho_ym2) / (2 * fieldset.dsamp_xy)

    rho_ym1_zm1 = fieldset.Rho[time, particle.depth - fieldset.dsamp_z, particle.lat - fieldset.dsamp_xy, particle.lon,]
    rho_ym1_zp1 = fieldset.Rho[time, particle.depth + fieldset.dsamp_z, particle.lat - fieldset.dsamp_xy, particle.lon,]
    dRhodz_ym1 = (rho_ym1_zp1 - rho_ym1_zm1) / (2 * fieldset.dsamp_z)

    dRhodx_zm1 = (rho_xp1_zm1 - rho_xm1_zm1) / (2 * fieldset.dsamp_xy)
    dRhodx_zp1 = (rho_xp1_zp1 - rho_xm1_zp1) / (2 * fieldset.dsamp_xy)

    dRhody_zm1 = (rho_yp1_zm1 - rho_ym1_zm1) / (2 * fieldset.dsamp_xy)
    dRhody_zp1 = (rho_yp1_zp1 - rho_ym1_zp1) / (2 * fieldset.dsamp_xy)

    rho_zp2 = fieldset.Rho[time, particle.depth + 2 * fieldset.dsamp_z, particle.lat, particle.lon]
    dRhodz_zp1 = (rho_zp2 - rho) / (2 * fieldset.dsamp_z)

    rho_zm2 = fieldset.Rho[time, particle.depth - 2 * fieldset.dsamp_z, particle.lat, particle.lon]
    dRhodz_zm1 = (rho - rho_zm2) / (2 * fieldset.dsamp_z)

    Sx_xp1 = -dRhodx_xp1 / dRhodz_xp1
    Sy_xp1 = -dRhody_xp1 / dRhodz_xp1
#     if math.fabs(Sx_xp1) > maxSlope:
#         Sx_xp1 = math.copysign(maxSlope, Sx_xp1)
#     if math.fabs(Sy_xp1) > maxSlope:
#         Sy_xp1 = math.copysign(maxSlope, Sy_xp1)
#     Sabs2_xp1 = Sx_xp1**2 + Sy_xp1**2

    Sx_xm1 = -dRhodx_xm1 / dRhodz_xm1
    Sy_xm1 = -dRhody_xm1 / dRhodz_xm1
#     if math.fabs(Sx_xm1) > maxSlope:
#         Sx_xm1 = math.copysign(maxSlope, Sx_xm1)
#     if math.fabs(Sy_xm1) > maxSlope:
#         Sy_xm1 = math.copysign(maxSlope, Sy_xm1)
#     Sabs2_xm1 = Sx_xm1**2 + Sy_xm1**2

    Sx_yp1 = -dRhodx_yp1 / dRhodz_yp1
    Sy_yp1 = -dRhody_yp1 / dRhodz_yp1
#     if math.fabs(Sx_yp1) > maxSlope:
#         Sx_yp1 = math.copysign(maxSlope, Sx_yp1)
#     if math.fabs(Sy_yp1) > maxSlope:
#         Sy_yp1 = math.copysign(maxSlope, Sy_yp1)
    Sabs2_yp1 = Sx_yp1**2 + Sy_yp1**2

    Sx_ym1 = -dRhodx_ym1 / dRhodz_ym1
    Sy_ym1 = -dRhody_ym1 / dRhodz_ym1
#     if math.fabs(Sx_ym1) > maxSlope:
#         Sx_ym1 = math.copysign(maxSlope, Sx_ym1)
#     if math.fabs(Sy_ym1) > maxSlope:
#         Sy_ym1 = math.copysign(maxSlope, Sy_ym1)
#     Sabs2_ym1 = Sx_ym1**2 + Sy_ym1**2

    Sx_zp1 = -dRhodx_zp1 / dRhodz_zp1
    Sy_zp1 = -dRhody_zp1 / dRhodz_zp1
#     if math.fabs(Sx_zp1) > maxSlope:
#         Sx_zp1 = math.copysign(maxSlope, Sx_zp1)
#     if math.fabs(Sy_zp1) > maxSlope:
#         Sy_zp1 = math.copysign(maxSlope, Sy_zp1)
    Sabs2_zp1 = Sx_zp1 ** 2 + Sy_zp1 ** 2

    Sx_zm1 = -dRhodx_zm1 / dRhodz_zm1
    Sy_zm1 = -dRhody_zm1 / dRhodz_zm1
#     if math.fabs(Sx_zm1) > maxSlope:
#         Sx_zm1 = math.copysign(maxSlope, Sx_zm1)
#     if math.fabs(Sy_zm1) > maxSlope:
#         Sy_zm1 = math.copysign(maxSlope, Sy_zm1)
    Sabs2_zm1 = Sx_zm1 ** 2 + Sy_zm1 ** 2

    # Filling out the Redi tensor
    k11 = fieldset.Ki
    k21 = 0
    k22 = fieldset.Ki
    k31 = fieldset.Ki * Sx
    k32 = fieldset.Ki * Sy
    k33 = fieldset.Ki * Sabs2
    
    dk11dx = 0
    dk12dy = 0
    dk13dz = fieldset.Ki * (Sx_zp1 - Sx_zm1) / (2 * fieldset.dsamp_z)
    dk21dx = 0
    dk22dy = 0
    dk23dz = fieldset.Ki * (Sy_zp1 - Sy_zm1) / (2 * fieldset.dsamp_z)
    dk31dx = fieldset.Ki * (Sx_xp1 - Sx_xm1) / (2 * fieldset.dsamp_xy)
    dk31dy = fieldset.Ki * (Sx_yp1 - Sx_ym1) / (2 * fieldset.dsamp_xy)
    dk31dz = fieldset.Ki * (Sx_zp1 - Sx_zm1) / (2 * fieldset.dsamp_z)
    dk32dx = fieldset.Ki * (Sy_xp1 - Sy_xm1) / (2 * fieldset.dsamp_xy)
    dk32dy = fieldset.Ki * (Sy_yp1 - Sy_ym1) / (2 * fieldset.dsamp_xy)
    dk32dz = fieldset.Ki * (Sy_zp1 - Sy_zm1) / (2 * fieldset.dsamp_z)
    dk33dz = fieldset.Ki * (Sabs2_zp1 - Sabs2_zm1) / (2 * fieldset.dsamp_z)

    # fmt: on
        