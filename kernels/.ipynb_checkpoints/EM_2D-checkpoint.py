import math
import parcels.rng as ParcelsRandom


def EM_2D(particle, fieldset, time):
    # Wiener increment with zero mean and std of sqrt(dt)
    dWx = ParcelsRandom.normalvariate(0, math.sqrt(math.fabs(particle.dt)))
    # dWy = ParcelsRandom.normalvariate(0, math.sqrt(math.fabs(particle.dt)))
    
    # Cholesky-Banachiewicz decomposed diffusivity tensor elements
    sigma11 = math.sqrt(2 * k11)
    sigma21 = 2 * k21 / sigma11
    sigma22 = 0 # math.sqrt(2 * particle.k22 - sigma21**2) = 0 independently of rho
    
    u, v = fieldset.UV[time, particle.depth, particle.lat, particle.lon]
    ax = u + dk11dx + dk12dy
    ay = v + dk21dx + dk22dy
#     dVxx = comp_Vxx(particle.lon)

    # Particle positions are updated only after evaluating all terms.  
    particle.lon += ax * particle.dt + sigma11 * dWx
    particle.lat += ay * particle.dt + sigma21 * dWx #+ sigma22 * dWy
    