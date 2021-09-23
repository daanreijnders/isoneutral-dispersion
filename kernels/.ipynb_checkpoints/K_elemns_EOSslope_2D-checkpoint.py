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

    particle.depth = particle.particle.depth
    particle.lat = particle.particle.lat
    particle.lon = particle.particle.lon
    fieldset.dsamp_xy = fieldset.fieldset.dsamp_xy
    fieldset.dsamp_z = fieldset.fieldset.dsamp_z

    rho = fieldset.Rho[time, particle.depth, particle.lat, particle.lon]
    rho_xm1 = fieldset.Rho[time, particle.depth, particle.lat, particle.lon - fieldset.dsamp_xy]
    rho_xp1 = fieldset.Rho[time, particle.depth, particle.lat, particle.lon + fieldset.dsamp_xy]
    rho_ym1 = fieldset.Rho[time, particle.depth, particle.lat - fieldset.dsamp_z, particle.lon]
    rho_yp1 = fieldset.Rho[time, particle.depth, particle.lat + fieldset.dsamp_z, particle.lon]
    rho_xm2 = fieldset.Rho[time, particle.depth, particle.lat, particle.lon - 2 * fieldset.dsamp_xy]
    rho_xp2 = fieldset.Rho[time, particle.depth, particle.lat, particle.lon + 2 * fieldset.dsamp_xy]
    rho_ym2 = fieldset.Rho[time, particle.depth, particle.lat - 2 * fieldset.dsamp_z, particle.lon]
    rho_yp2 = fieldset.Rho[time, particle.depth, particle.lat + 2 * fieldset.dsamp_z, particle.lon]
    rho_xp1_ym1 = fieldset.Rho[time, particle.depth, particle.lat - fieldset.dsamp_z, particle.lon + fieldset.dsamp_xy,]
    rho_xp1_yp1 = fieldset.Rho[time, particle.depth, particle.lat + fieldset.dsamp_z, particle.lon + fieldset.dsamp_xy,]
    rho_xm1_ym1 = fieldset.Rho[time, particle.depth, particle.lat - fieldset.dsamp_z, particle.lon - fieldset.dsamp_xy,]
    rho_xm1_yp1 = fieldset.Rho[time, particle.depth, particle.lat + fieldset.dsamp_z, particle.lon - fieldset.dsamp_xy,]
    
    dRhodx = (rho_xp1 - rho_xm1) / (2 * fieldset.dsamp_xy)
    dRhody = (rho_yp1 - rho_ym1) / (2 * fieldset.dsamp_z)
    dRhodx_xm1 = (rho - rho_xm2) / (2 * fieldset.dsamp_xy)
    dRhodx_xp1 = (rho_xp2 - rho) / (2 * fieldset.dsamp_xy)
    dRhody_xp1 = (rho_xp1_yp1 - rho_xp1_ym1) / (2 * fieldset.dsamp_z)
    dRhody_xm1 = (rho_xm1_yp1 - rho_xm1_ym1) / (2 * fieldset.dsamp_z)
    dRhodx_yp1 = (rho_xp1_yp1 - rho_xm1_yp1) / (2 * fieldset.dsamp_xy)
    dRhody_yp1 = (rho_yp2 - rho) / (2 * fieldset.dsamp_z)
    dRhodx_ym1 = (rho_xp1_ym1 - rho_xm1_ym1) / (2 * fieldset.dsamp_xy)
    dRhody_ym1 = (rho - rho_ym2) / (2 * fieldset.dsamp_z)

    S = -dRhodx / dRhody
    S_xp1 = -dRhodx_xp1 / dRhody_xp1
    S_xm1 = -dRhodx_xm1 / dRhody_xm1
    S_yp1 = -dRhodx_yp1 / dRhody_yp1
    S_ym1 = -dRhodx_ym1 / dRhody_ym1

    # Filling out the Redi tensor
    k11 = fieldset.Ki
    k21 = fieldset.Ki * S
    k22 = fieldset.Ki * S**2
    
    dk11dx = 0
    dk12dy = fieldset.Ki * (S_yp1 - S_ym1) / (2 * fieldset.dsamp_z)
    dk21dx = fieldset.Ki * (S_xp1 - S_xm1) / (2 * fieldset.dsamp_xy)
    dk22dy = fieldset.Ki * (S_yp1**2 - S_ym1**2) / (2 * fieldset.dsamp_z)
    # fmt: on
        
        
        