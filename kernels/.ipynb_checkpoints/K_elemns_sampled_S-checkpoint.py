import math
import parcels.rng as ParcelsRandom

def K_elemns_sampled_S(particle, fieldset, time):
# fmt: off
    sx = fieldset.Sx[time, particle.depth, particle.lat, particle.lon]
    sy = fieldset.Sy[time, particle.depth, particle.lat, particle.lon]

    k11 = fieldset.Ki
    k21 = 0
    k22 = fieldset.Ki
    k31 = fieldset.Ki * sx
    k32 = fieldset.Ki * sy
    k33 = fieldset.Ki * (sx**2 + sy**2)
    
    dsxdx = (fieldset.Sx[time, particle.depth, particle.lat, particle.lon + fieldset.dsamp_xy] - fieldset.Sx[time, particle.depth, particle.lat, particle.lon - fieldset.dsamp_xy]) / (2 * fieldset.dsamp_xy)
    dsxdy = (fieldset.Sx[time, particle.depth, particle.lat + fieldset.dsamp_xy, particle.lon] - fieldset.Sx[time, particle.depth, particle.lat - fieldset.dsamp_xy, particle.lon]) / (2 * fieldset.dsamp_xy)
    dsxdz = (fieldset.Sx[time, particle.depth + fieldset.dsamp_z, particle.lat, particle.lon] - fieldset.Sx[time, particle.depth - fieldset.dsamp_z, particle.lat, particle.lon]) / (2 * fieldset.dsamp_z)
    dsydx = (fieldset.Sy[time, particle.depth, particle.lat, particle.lon + fieldset.dsamp_xy] - fieldset.Sy[time, particle.depth, particle.lat, particle.lon - fieldset.dsamp_xy]) / (2 * fieldset.dsamp_xy)
    dsydy = (fieldset.Sy[time, particle.depth, particle.lat + fieldset.dsamp_xy, particle.lon] - fieldset.Sy[time, particle.depth, particle.lat - fieldset.dsamp_xy, particle.lon]) / (2 * fieldset.dsamp_xy)
    dsydz = (fieldset.Sy[time, particle.depth + fieldset.dsamp_z, particle.lat, particle.lon] - fieldset.Sy[time, particle.depth - fieldset.dsamp_z, particle.lat, particle.lon]) / (2 * fieldset.dsamp_z)
    
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