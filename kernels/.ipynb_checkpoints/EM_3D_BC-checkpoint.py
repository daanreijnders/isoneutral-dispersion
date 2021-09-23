import math
import parcels.rng as ParcelsRandom

def EM_3D_BC(particle, fieldset, time):
    # fmt: off  
    # Independent Wiener increments with zero mean and std of sqrt(dt)
    dWx = ParcelsRandom.normalvariate(0, math.sqrt(math.fabs(particle.dt)))
    dWy = ParcelsRandom.normalvariate(0, math.sqrt(math.fabs(particle.dt)))
    # dWz = random.normalvariate(0, math.sqrt(math.fabs(particle.dt)))
    
    # Cholesky-Banachiewicz decomposed diffusivity tensor elements
    sigma11 = math.sqrt(2 * k11)
    sigma21 = 0 #2 * k21 / sigma11
    sigma22 = math.sqrt(2 * k22) # math.sqrt(2 * k22 - sigma21**2)
    sigma31 = 2 * k31 / sigma11
    sigma32 = 2 * k32 / sigma22 #(2 * k32 - sigma31 * sigma21) / sigma22
    # sigma33 should cancel out to zero in small slope approximation
    # sigma33 = math.sqrt(2 * k33 - sigma31 ** 2 - sigma32 ** 2)

    u, v, w = fieldset.UVW[time, particle.depth, particle.lat, particle.lon]
    ax = u + dk11dx + dk12dy + dk13dz
    ay = v + dk21dx + dk22dy + dk23dz
    az = w + dk31dx + dk32dy + dk33dz

    dLon = ax * particle.dt + sigma11 * dWx
    dLat = ay * particle.dt + sigma21 * dWx + sigma22 * dWy
    dDepth = az * particle.dt + sigma31 * dWx + sigma32 * dWy
  
    if (particle.lat + dLat < fieldset.northBound) and (particle.lat + dLat > fieldset.southBound) and (particle.depth + dDepth < fieldset.upperBound):
        particle.lon += dLon
        particle.lat += dLat
        particle.depth += dDepth 
        if particle.dt < fieldset.max_dt:
            particle.update_next_dt(particle.dt * 2)
    else:
        particle.dt /= 2
        if particle.dt < fieldset.min_dt:
            if particle.lat + dLat > fieldset.northBound:
                particle.lon += dLon
                particle.lat -= dLat
                particle.depth += dDepth 
            elif particle.depth + dDepth > fieldset.upperBound:
                particle.lon += dLon
                particle.lat += dLat
                particle.depth -= dDepth 
        else:
            return OperationCode.Repeat
    # fmt: on
