import math
import parcels.rng as ParcelsRandom

def EX1_3D(particle, fieldset, time):
    # fmt: off  
    # Independent Wiener increments with zero mean and std of sqrt(dt)
    dWx = ParcelsRandom.normalvariate(0, math.sqrt(math.fabs(particle.dt)))
    dWy = ParcelsRandom.normalvariate(0, math.sqrt(math.fabs(particle.dt)))
    # dWz = ParcelsRandom.normalvariate(0, math.sqrt(math.fabs(particle.dt)))
    sqrt2 = math.sqrt(2)
    
    # Cholesky-Banachiewicz decomposed diffusivity tensor elements
    sigma11 = math.sqrt(2 * k11)
    sigma21 = 0 #2 * k21 / sigma11
    sigma22 = math.sqrt(2 * k22) # math.sqrt(2 * k22 - sigma21**2)
    sigma31 = 2 * k31 / sigma11
    sigma32 = 2 * k32 / sigma22 #(2 * k32 - sigma31 * sigma21) / sigma22
    # sigma33 should cancel out to zero in small slope approximation
    # sigma33 = math.sqrt(2 * k33 - sigma31 ** 2 - sigma32 ** 2)
    
    dsigma31dx = 2 * dk31dx / sigma11
    dsigma31dy = 2 * dk31dy / sigma11
    dsigma31dz = 2 * dk31dz / sigma11
    dsigma32dx = 2 * dk32dx / sigma11
    dsigma32dy = 2 * dk32dy / sigma11
    dsigma32dz = 2 * dk32dz / sigma11

    u, v, w = fieldset.UVW[time, particle.depth, particle.lat, particle.lon]
    ax = u + dk11dx + dk12dy + dk13dz
    ay = v + dk21dx + dk22dy + dk23dz
    az = w + dk31dx + dk32dy + dk33dz
    
    I11 = 0.5 * (dWx**2 - particle.dt)
    I22 = 0.5 * (dWy**2 - particle.dt)
    
    xi1 = 1/math.sqrt(particle.dt) * dWx
    xi2 = 1/math.sqrt(particle.dt) * dWy
    #p = 10
    r = 1
    rqp = 0
    sum12 = 0
    sum21 = 0
    mu1 = ParcelsRandom.normalvariate(0, 1)
    mu2 = ParcelsRandom.normalvariate(0, 1)
    while r <= fieldset.expansion_terms:
        zeta1 = ParcelsRandom.normalvariate(0, 1)
        zeta2 = ParcelsRandom.normalvariate(0, 1)
        eta1 = ParcelsRandom.normalvariate(0, 1)
        eta2 = ParcelsRandom.normalvariate(0, 1)
        rqp += 1 / (r**2)
        sum12 += 1/r * (zeta1 * (sqrt2 * xi2 + eta2) - zeta2 * (sqrt2 * xi1 + eta1))
        sum21 += 1/r * (zeta2 * (sqrt2 * xi1 + eta1) - zeta1 * (sqrt2 * xi2 + eta2))
        r += 1
    sqrt_rho_p = math.sqrt(1./12 - 1 / (2 * math.pi**2) * rqp)
    I12 = particle.dt * (0.5 * xi1 * xi2 + sqrt_rho_p * (mu1 * xi2 - mu2 * xi1)) \
        + particle.dt / (2 * math.pi) * sum12
    I21 = particle.dt * (0.5 * xi2 * xi1 + sqrt_rho_p * (mu2 * xi1 - mu1 * xi2)) \
        + particle.dt / (2 * math.pi) * sum21
    #     Particle positions are updated only after evaluating all terms.
    particle.lon += ax * particle.dt + sigma11 * dWx
    particle.lat += ay * particle.dt + sigma21 * dWx + sigma22 * dWy
    particle.depth += az * particle.dt + sigma31 * dWx + sigma32 * dWy + (
                      sigma11 * dsigma31dx * I11 + sigma11 * dsigma32dx * I12 \
                    + sigma22 * dsigma31dy * I21 + sigma22 * dsigma32dy * I22 \
                    + sigma31 * dsigma31dz * I11 + sigma31 * dsigma32dz * I12 \
                    + sigma32 * dsigma31dz * I21 + sigma32 * dsigma32dz * I22
    )
    
    # fmt: on