import math
import parcels.rng as ParcelsRandom

def K_elemns_sampled_S_deriv(particle, fieldset, time):
# fmt: off
    sx = fieldset.Sx[time, particle.depth, particle.lat, particle.lon]
    sy = fieldset.Sy[time, particle.depth, particle.lat, particle.lon]

    k11 = fieldset.Ki
    k21 = 0
    k22 = fieldset.Ki
    k31 = fieldset.Ki * sx
    k32 = fieldset.Ki * sy
    k33 = fieldset.Ki * (sx**2 + sy**2)
    
    dsxdx = fieldset.dSxdx[time, particle.depth, particle.lat, particle.lon]
    dsxdy = fieldset.dSxdy[time, particle.depth, particle.lat, particle.lon]
    dsxdz = fieldset.dSxdz[time, particle.depth, particle.lat, particle.lon]
    dsydx = fieldset.dSxdx[time, particle.depth, particle.lat, particle.lon]
    dsydy = fieldset.dSydy[time, particle.depth, particle.lat, particle.lon]
    dsydz = fieldset.dSydz[time, particle.depth, particle.lat, particle.lon]

    dk11dx = 0
    dk12dy = 0
    dk13dz = fieldset.Ki * dsxdz
    dk21dx = 0
    dk22dy = 0
    dk23dz = fieldset.Ki * dsydz
    dk31dx = fieldset.Ki * dsxdx
    dk31dy = fieldset.Ki * dsxdy
    dk31dz = fieldset.Ki * dsxdz
    dk32dx = fieldset.Ki * dsydx
    dk32dy = fieldset.Ki * dsydy
    dk32dz = fieldset.Ki * dsydz
    dk33dz = Ki * 2 * (sx * dsxdz + sy * dsydz)

    # fmt: on