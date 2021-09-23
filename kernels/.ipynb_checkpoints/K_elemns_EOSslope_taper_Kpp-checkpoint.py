import math
import parcels.rng as ParcelsRandom


def K_elemns_EOSslope_taper_Kpp(particle, fieldset, time):
    
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
    rho_xm1 = fieldset.Rho[time, particle.depth, particle.lat, particle.lon - dx]
    rho_xp1 = fieldset.Rho[time, particle.depth, particle.lat, particle.lon + dx]
    dRhodx = (rho_xp1 - rho_xm1) / (2 * dx)

    rho_ym1 = fieldset.Rho[time, particle.depth, particle.lat - dy, particle.lon]
    rho_yp1 = fieldset.Rho[time, particle.depth, particle.lat + dy, particle.lon]
    dRhody = (rho_yp1 - rho_ym1) / (2 * dy)

    rho_zm1 = fieldset.Rho[time, particle.depth - dz, particle.lat, particle.lon]
    rho_zp1 = fieldset.Rho[time, particle.depth + dz, particle.lat, particle.lon]
    dRhodz = (rho_zp1 - rho_zm1) / (2 * dz)

    Sx = -dRhodx / dRhodz
    Sy = -dRhody / dRhodz

    Sabs2 = Sx ** 2 + Sy ** 2

    # Computing gradients in K requires knowing the slope at adjacent points
    #     Note: dKxx and dKyy are zero (except in case of variable diffusivity).
    #           dKzx, dKzy, and dKzz are often zero
    #           (when [-2*dres, 2*dres] falls between two density gridpoints)
    #           Computing these implies computing slopes for all adjacent points (6)
    #           for which we need to know the density at even more neighboring points
    #           (6 * 6 = 36, but without overlap 25 total interpoparticle.lations vs 6)
    rho_xp2 = fieldset.Rho[time, particle.depth, particle.lat, particle.lon + 2 * dx]
    dRhodx_xp1 = (rho_xp2 - rho) / (2 * dx)

    rho_xp1_ym1 = fieldset.Rho[time, particle.depth, particle.lat - dy, particle.lon + dx]
    rho_xp1_yp1 = fieldset.Rho[time, particle.depth, particle.lat + dy, particle.lon + dx]
    dRhody_xp1 = (rho_xp1_yp1 - rho_xp1_ym1) / (2 * dy)

    rho_xp1_zm1 = fieldset.Rho[time, particle.depth - dz, particle.lat, particle.lon + dx]
    rho_xp1_zp1 = fieldset.Rho[time, particle.depth + dz, particle.lat, particle.lon + dx]
    dRhodz_xp1 = (rho_xp1_zp1 - rho_xp1_zm1) / (2 * dz)

    rho_xm2 = fieldset.Rho[time, particle.depth, particle.lat, particle.lon - 2 * dx]
    dRhodx_xm1 = (rho - rho_xm2) / (2 * dx)

    rho_xm1_ym1 = fieldset.Rho[time, particle.depth, particle.lat - dy, particle.lon - dx]
    rho_xm1_yp1 = fieldset.Rho[time, particle.depth, particle.lat + dy, particle.lon - dx]
    dRhody_xm1 = (rho_xm1_yp1 - rho_xm1_ym1) / (2 * dy)

    rho_xm1_zm1 = fieldset.Rho[time, particle.depth - dz, particle.lat, particle.lon - dx]
    rho_xm1_zp1 = fieldset.Rho[time, particle.depth + dz, particle.lat, particle.lon - dx]
    dRhodz_xm1 = (rho_xm1_zp1 - rho_xm1_zm1) / (2 * dz)
    dRhodx_yp1 = (rho_xp1_yp1 - rho_xm1_yp1) / (2 * dx)

    rho_yp2 = fieldset.Rho[time, particle.depth, particle.lat + 2 * dy, particle.lon]
    dRhody_yp1 = (rho_yp2 - rho) / (2 * dy)

    rho_yp1_zm1 = fieldset.Rho[time, particle.depth - dz, particle.lat + dy, particle.lon]
    rho_yp1_zp1 = fieldset.Rho[time, particle.depth + dz, particle.lat + dy, particle.lon]
    dRhodz_yp1 = (rho_yp1_zp1 - rho_yp1_zm1) / (2 * dz)
    dRhodx_ym1 = (rho_xp1_ym1 - rho_xm1_ym1) / (2 * dx)

    rho_ym2 = fieldset.Rho[time, particle.depth, particle.lat - 2 * dy, particle.lon]
    dRhody_ym1 = (rho - rho_ym2) / (2 * dy)

    rho_ym1_zm1 = fieldset.Rho[time, particle.depth - dz, particle.lat - dy, particle.lon]
    rho_ym1_zp1 = fieldset.Rho[time, particle.depth + dz, particle.lat - dy, particle.lon]
    dRhodz_ym1 = (rho_ym1_zp1 - rho_ym1_zm1) / (2 * dz)

    dRhodx_zm1 = (rho_xp1_zm1 - rho_xm1_zm1) / (2 * dx)
    dRhodx_zp1 = (rho_xp1_zp1 - rho_xm1_zp1) / (2 * dx)

    dRhody_zm1 = (rho_yp1_zm1 - rho_ym1_zm1) / (2 * dy)
    dRhody_zp1 = (rho_yp1_zp1 - rho_ym1_zp1) / (2 * dy)

    rho_zp2 = fieldset.Rho[time, particle.depth + 2 * dz, particle.lat, particle.lon]
    dRhodz_zp1 = (rho_zp2 - rho) / (2 * dz)

    rho_zm2 = fieldset.Rho[time, particle.depth - 2 * dz, particle.lat, particle.lon]
    dRhodz_zm1 = (rho - rho_zm2) / (2 * dz)

    Sx_xp1 = -dRhodx_xp1 / dRhodz_xp1
    Sy_xp1 = -dRhody_xp1 / dRhodz_xp1
    Sabs2_xp1 = Sx_xp1**2 + Sy_xp1**2

    Sx_xm1 = -dRhodx_xm1 / dRhodz_xm1
    Sy_xm1 = -dRhody_xm1 / dRhodz_xm1
    Sabs2_xm1 = Sx_xm1**2 + Sy_xm1**2

    Sx_yp1 = -dRhodx_yp1 / dRhodz_yp1
    Sy_yp1 = -dRhody_yp1 / dRhodz_yp1
    Sabs2_yp1 = Sx_yp1**2 + Sy_yp1**2

    Sx_ym1 = -dRhodx_ym1 / dRhodz_ym1
    Sy_ym1 = -dRhody_ym1 / dRhodz_ym1
    Sabs2_ym1 = Sx_ym1**2 + Sy_ym1**2

    Sx_zp1 = -dRhodx_zp1 / dRhodz_zp1
    Sy_zp1 = -dRhody_zp1 / dRhodz_zp1
    Sabs2_zp1 = Sx_zp1 ** 2 + Sy_zp1 ** 2

    Sx_zm1 = -dRhodx_zm1 / dRhodz_zm1
    Sy_zm1 = -dRhody_zm1 / dRhodz_zm1
    Sabs2_zm1 = Sx_zm1 ** 2 + Sy_zm1 ** 2

    # Tapering
    if Sabs2 > 0:
        taper = 0.5 * (1. + math.tanh((fieldset.Sc - math.sqrt(Sabs2)) / fieldset.Sd))
    else:
        taper = 1
    
    if Sabs2_xp1 > 0:
        taper_xp1 = 0.5 * (1. + math.tanh((fieldset.Sc - math.sqrt(Sabs2_xp1)) / fieldset.Sd))
    else:
        taper_xp1 = 1
    
    if Sabs2_xm1 > 0:
        taper_xm1 = 0.5 * (1. + math.tanh((fieldset.Sc - math.sqrt(Sabs2_xm1)) / fieldset.Sd))
    else:
        taper_xm1 = 1
        
    if Sabs2_yp1 > 0:
        taper_yp1 = 0.5 * (1. + math.tanh((fieldset.Sc - math.sqrt(Sabs2_yp1)) / fieldset.Sd))
    else:
        taper_yp1 = 1
        
    if Sabs2_ym1 > 0:
        taper_ym1 = 0.5 * (1. + math.tanh((fieldset.Sc - math.sqrt(Sabs2_ym1)) / fieldset.Sd))
    else:
        taper_ym1 = 1
        
    if Sabs2_zp1 > 0:
        taper_zp1 = 0.5 * (1. + math.tanh((fieldset.Sc - math.sqrt(Sabs2_zp1)) / fieldset.Sd))
    else:
        taper_zp1 = 1
        
    if Sabs2_zm1 > 0:
        taper_zm1 = 0.5 * (1. + math.tanh((fieldset.Sc - math.sqrt(Sabs2_zm1)) / fieldset.Sd))
    else:
        taper_zm1 = 1
    
    # Filling out the Redi tensor
    k11 = taper * fieldset.Ki
    k21 = 0
    k22 = taper * fieldset.Ki
    k31 = taper * fieldset.Ki * Sx
    k32 = taper * fieldset.Ki * Sy
    k33 = taper * fieldset.Ki * Sabs2
    
    kv = fieldset.Kv[time, particle.depth, particle.lat, particle.lon]
    kv_zp1 = fieldset.Kv[time, particle.depth + dz, particle.lat, particle.lon] 
    kv_zm1 = fieldset.Kv[time, particle.depth - dz, particle.lat, particle.lon] 
    dkvdz = (kv_zp1 - kv_zm1) / (2 * dz)
    kv_xp1 = fieldset.Kv[time, particle.depth, particle.lat, particle.lon + dx] 
    kv_xm1 = fieldset.Kv[time, particle.depth, particle.lat, particle.lon - dx] 
    dkvdx = (kv_xp1 - kv_xm1) / (2 * dx)
    kv_yp1 = fieldset.Kv[time, particle.depth, particle.lat + dy, particle.lon] 
    kv_ym1 = fieldset.Kv[time, particle.depth, particle.lat - dy, particle.lon] 
    dkvdy = (kv_yp1 - kv_ym1) / (2 * dy)
    
    dk11dx = 0
    dk12dy = 0
    dk13dz = fieldset.Ki * ((taper_zp1 - taper_zm1) / (2 * dz) * Sx + taper * (Sx_zp1 - Sx_zm1) / (2 * dz))
    dk21dx = 0
    dk22dy = 0
    dk23dz = fieldset.Ki * ((taper_zp1 - taper_zm1) / (2 * dz) * Sy + taper * (Sy_zp1 - Sy_zm1) / (2 * dz))
    dk31dx = fieldset.Ki * ((taper_xp1 - taper_xm1) / (2 * dx) * Sx + taper * (Sx_xp1 - Sx_xm1) / (2 * dx))
    dk31dy = fieldset.Ki * ((taper_yp1 - taper_ym1) / (2 * dy) * Sx + taper * (Sx_yp1 - Sx_ym1) / (2 * dy))
    dk31dz = fieldset.Ki * ((taper_zp1 - taper_zm1) / (2 * dz) * Sx + taper * (Sx_zp1 - Sx_zm1) / (2 * dz))
    dk32dx = fieldset.Ki * ((taper_xp1 - taper_xm1) / (2 * dx) * Sy + taper * (Sy_xp1 - Sy_xm1) / (2 * dx))
    dk32dy = fieldset.Ki * ((taper_yp1 - taper_ym1) / (2 * dy) * Sy + taper * (Sy_yp1 - Sy_ym1) / (2 * dy))
    dk32dz = fieldset.Ki * ((taper_zp1 - taper_zm1) / (2 * dz) * Sy + taper * (Sy_zp1 - Sy_zm1) / (2 * dz))
    dk33dz = fieldset.Ki * ((taper_zp1 - taper_zm1) / (2 * dz) * Sabs2 + taper * (Sabs2_zp1 - Sabs2_zm1) / (2 * dz))

    # fmt: on
        