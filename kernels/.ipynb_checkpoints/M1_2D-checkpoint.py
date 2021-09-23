import math
import parcels.rng as ParcelsRandom

def M1_2D(particle, fieldset, time):
    # Wiener increment with zero mean and std of sqrt(dt)
    dWx = ParcelsRandom.normalvariate(0, math.sqrt(math.fabs(particle.dt)))
    # dWy = random.normalvariate(0, math.sqrt(math.fabs(particle.dt)))
    
    # Cholesky-Banachiewicz decomposed diffusivity tensor elements
    sigma11 = math.sqrt(2 * k11)
    sigma21 = 2 * k21 / sigma11
    sigma22 = 0 # math.sqrt(2 * k22 - sigma21**2) = 0 independently of rho
    
    dsigma11dx = dk11dx / sigma11
    dsigma11dy = dk11dy / sigma11
    dsigma21dx = (sigma11 * dk21dx - k21 * dsigma11dx) / k11
    dsigma21dy = (sigma11 * dk12dy - k21 * dsigma11dy) / k11
    
    u, v = fieldset.UV[time, particle.depth, particle.lat, particle.lon]
    ax = u + particle.dk11dx + particle.dk12dy
    ay = v + particle.dk21dx + particle.dk22dy
#     dVxx = comp_Vxx(particle.lon)

    # Particle positions are updated only after evaluating all terms.  
    particle.lon += ax * particle.dt + sigma11 * dWx + sigma11 * dsigma11dx * 0.5 * (dWx ** 2 - particle.dt) \
                                                     + sigma21 * dsigma11dy * 0.5 * (dWx ** 2 - particle.dt)
    particle.lat += ay * particle.dt + sigma21 * dWx + sigma11 * dsigma21dx * 0.5 * (dWx ** 2 - particle.dt) \
                                                     + sigma21 * dsigma21dy * 0.5 * (dWx ** 2 - particle.dt)
    