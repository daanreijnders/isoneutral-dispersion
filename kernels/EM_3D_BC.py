import math
import parcels.rng as ParcelsRandom

def EM_3D_BC(particle, fieldset, time):
    if k11 > 0:
        # Independent Wiener increments with zero mean and std of sqrt(dt)
        dWx = ParcelsRandom.normalvariate(0, math.sqrt(math.fabs(particle.dt)))
        dWy = ParcelsRandom.normalvariate(0, math.sqrt(math.fabs(particle.dt)))
        dWz = ParcelsRandom.normalvariate(0, math.sqrt(math.fabs(particle.dt)))

        # Cholesky-Banachiewicz decomposed diffusivity tensor elements
        # `math.fabs` to assure computation when the result is -epsilon
        b11 = math.sqrt(2. * k11)
        b21 = 2. * k21 / b11
        b22 = math.sqrt(math.fabs(2. * k22 - b21**2))
        b31 = 2. * k31 / b11
        b32 = (2. * k32 - b31 * b21) / b22
        b33 = math.sqrt(math.fabs(2. * k33 - b31**2 - b32**2))

        u, v, w = fieldset.UVW[time, particle.depth, particle.lat, particle.lon]
        ax = u + dk11dx + dk21dy + dk31dz
        ay = v + dk21dx + dk22dy + dk32dz
        az = w + dk31dx + dk32dy + dk33dz

        particle.lon += ax * particle.dt + b11 * dWx
        particle.lat += ay * particle.dt + b21 * dWx + b22 * dWy
        particle.depth += az * particle.dt + b31 * dWx + b32 * dWy + b33 * dWz
        
        # Boundary condition
        if particle.lat > fieldset.northBound:
            particle.lat = 2. * fieldset.northBound - particle.lat
        elif particle.lat < fieldset.southBound:
            particle.lat = 2. * fieldset.southBound - particle.lat
        if particle.depth > fieldset.upperBound:
            particle.depth = 2. * fieldset.upperBound - particle.depth
        elif particle.depth < fieldset.lowerBound:
            particle.depth = 2. * fieldset.lowerBound - particle.depth 
            
    else:
        u, v, w = fieldset.UVW[time, particle.depth, particle.lat, particle.lon]
        particle.lon += u * particle.dt
        particle.lat += v * particle.dt
        particle.depth += w * particle.dt
   
