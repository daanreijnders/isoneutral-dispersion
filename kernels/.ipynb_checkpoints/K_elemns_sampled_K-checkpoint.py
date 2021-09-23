import math
import parcels.rng as ParcelsRandom

def K_elemns_sampled_K(particle, fieldset, time):
# fmt: off
    k11 = fieldset.K11[time, particle.depth, particle.lat, particle.lon]
    k21 = fieldset.K21[time, particle.depth, particle.lat, particle.lon]
    k22 = fieldset.K22[time, particle.depth, particle.lat, particle.lon]
    k31 = fieldset.K31[time, particle.depth, particle.lat, particle.lon]
    k32 = fieldset.K32[time, particle.depth, particle.lat, particle.lon]
    k33 = fieldset.K33[time, particle.depth, particle.lat, particle.lon]
    
    dk11dx = 0
    dk12dy = 0
    dk13dz = (fieldset.K31[time, particle.depth + fieldset.dsamp_z, particle.lat, particle.lon] - fieldset.K31[time, particle.depth - fieldset.dsamp_z, particle.lat, particle.lon]) / (2 * fieldset.dsamp_z)
    dk21dx = 0
    dk22dy = 0
    dk23dz = (fieldset.K32[time, particle.depth + fieldset.dsamp_z, particle.lat, particle.lon] - fieldset.K32[time, particle.depth - fieldset.dsamp_z, particle.lat, particle.lon]) / (2 * fieldset.dsamp_z)
    dk31dx = (fieldset.K31[time, particle.depth, particle.lat, particle.lon + fieldset.dsamp_xy] - fieldset.K31[time, particle.depth, particle.lat, particle.lon - fieldset.dsamp_xy]) / (2 * fieldset.dsamp_xy)
    dk31dy = (fieldset.K31[time, particle.depth, particle.lat + fieldset.dsamp_xy, particle.lon] - fieldset.K31[time, particle.depth, particle.lat - fieldset.dsamp_xy, particle.lon]) / (2 * fieldset.dsamp_xy)
    dk31dz = (fieldset.K31[time, particle.depth + fieldset.dsamp_z, particle.lat, particle.lon] - fieldset.K31[time, particle.depth - fieldset.dsamp_z, particle.lat, particle.lon]) / (2 * fieldset.dsamp_z)
    dk32dx = (fieldset.K32[time, particle.depth, particle.lat, particle.lon + fieldset.dsamp_xy] - fieldset.K32[time, particle.depth, particle.lat, particle.lon - fieldset.dsamp_xy]) / (2 * fieldset.dsamp_xy)
    dk32dy = (fieldset.K32[time, particle.depth, particle.lat + fieldset.dsamp_xy, particle.lon] - fieldset.K32[time, particle.depth, particle.lat - fieldset.dsamp_xy, particle.lon]) / (2 * fieldset.dsamp_xy)
    dk32dz = (fieldset.K32[time, particle.depth + fieldset.dsamp_z, particle.lat, particle.lon] - fieldset.K32[time, particle.depth - fieldset.dsamp_z, particle.lat, particle.lon]) / (2 * fieldset.dsamp_z)
    dk33dz = (fieldset.K33[time, particle.depth + fieldset.dsamp_z, particle.lat, particle.lon] - fieldset.K33[time, particle.depth - fieldset.dsamp_z, particle.lat, particle.lon]) / (2 * fieldset.dsamp_z)

    # fmt: on