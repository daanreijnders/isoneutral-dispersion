import math
import parcels.rng as ParcelsRandom

def M1_3D_BC_KPP(particle, fieldset, time):
    # fmt: off  
    # Independent Wiener increments with zero mean and std of sqrt(dt)
    dWx = ParcelsRandom.normalvariate(0, math.sqrt(math.fabs(particle.dt)))
    dWy = ParcelsRandom.normalvariate(0, math.sqrt(math.fabs(particle.dt)))
    dWz = ParcelsRandom.normalvariate(0, math.sqrt(math.fabs(particle.dt)))
    sqrt2 = math.sqrt(2)
    
    
    # Cholesky-Banachiewicz decomposed diffusivity tensor elements
    sigma11 = math.sqrt(2 * k11)
    sigma21 = 0 #2 * k21 / sigma11
    sigma22 = math.sqrt(2 * k22) # math.sqrt(2 * k22 - sigma21**2)
    sigma31 = 2 * k31 / sigma11
    sigma32 = 2 * k32 / sigma22 #(2 * k32 - sigma31 * sigma21) / sigma22
    # sigma33's isopycnal elements cancel out, but its vertical K contribution does not:
    sigma33 = math.sqrt(2 * kv)
    
    dsigma31dx = 2 * dk31dx / sigma11
    dsigma31dy = 2 * dk31dy / sigma11
    dsigma31dz = 2 * dk31dz / sigma11
    dsigma32dx = 2 * dk32dx / sigma11
    dsigma32dy = 2 * dk32dy / sigma11
    dsigma32dz = 2 * dk32dz / sigma11
    dsigma33dx = dkvdx / sigma33
    dsigma33dy = dkvdy / sigma33
    dsigma33dz = dkvdz / sigma33

    u, v, w = fieldset.UVW[time, particle.depth, particle.lat, particle.lon]
    ax = u + dk11dx + dk12dy + dk13dz
    ay = v + dk21dx + dk22dy + dk23dz
    az = w + dk31dx + dk32dy + dk33dz + dkvdz
    
    I11 = 0.5 * (dWx**2 - particle.dt)
    I22 = 0.5 * (dWy**2 - particle.dt)
    I33 = 0.5 * (dWz**2 - particle.dt)
    
    xi1 = 1/math.sqrt(particle.dt) * dWx
    xi2 = 1/math.sqrt(particle.dt) * dWy
    xi3 = 1/math.sqrt(particle.dt) * dWz
    r = 1
    rqp = 0
    sum12 = 0
    sum13 = 0
    sum21 = 0
    sum23 = 0
    sum31 = 0
    sum32 = 0
    mu1 = ParcelsRandom.normalvariate(0, 1)
    mu2 = ParcelsRandom.normalvariate(0, 1)
    mu3 = ParcelsRandom.normalvariate(0, 1)
    while r <= fieldset.expansion_terms:
        zeta1 = ParcelsRandom.normalvariate(0, 1)
        zeta2 = ParcelsRandom.normalvariate(0, 1)
        zeta3 = ParcelsRandom.normalvariate(0, 1)
        eta1 = ParcelsRandom.normalvariate(0, 1)
        eta2 = ParcelsRandom.normalvariate(0, 1)
        eta3 = ParcelsRandom.normalvariate(0, 1)
        rqp += 1 / (r**2)
        sum12 += 1/r * (zeta1 * (sqrt2 * xi2 + eta2) - zeta2 * (sqrt2 * xi1 + eta1))
        sum13 += 1/r * (zeta1 * (sqrt2 * xi3 + eta3) - zeta3 * (sqrt2 * xi1 + eta1))
        sum21 += 1/r * (zeta2 * (sqrt2 * xi1 + eta1) - zeta1 * (sqrt2 * xi2 + eta2))
        sum23 += 1/r * (zeta2 * (sqrt2 * xi3 + eta3) - zeta3 * (sqrt2 * xi2 + eta2))
        sum31 += 1/r * (zeta3 * (sqrt2 * xi1 + eta1) - zeta1 * (sqrt2 * xi3 + eta3))
        sum32 += 1/r * (zeta3 * (sqrt2 * xi2 + eta2) - zeta2 * (sqrt2 * xi3 + eta3))
        r += 1
    sqrt_rho_p = math.sqrt(1./12 - 1 / (2 * math.pi**2) * rqp)
    I12 = particle.dt * (0.5 * xi1 * xi2 + sqrt_rho_p * (mu1 * xi2 - mu2 * xi1)) \
        + particle.dt / (2 * math.pi) * sum12
    I13 = particle.dt * (0.5 * xi1 * xi3 + sqrt_rho_p * (mu1 * xi3 - mu3 * xi1)) \
        + particle.dt / (2 * math.pi) * sum13
    I21 = particle.dt * (0.5 * xi2 * xi1 + sqrt_rho_p * (mu2 * xi1 - mu1 * xi2)) \
        + particle.dt / (2 * math.pi) * sum21
    I23 = particle.dt * (0.5 * xi2 * xi3 + sqrt_rho_p * (mu2 * xi3 - mu3 * xi2)) \
        + particle.dt / (2 * math.pi) * sum23
    I31 = particle.dt * (0.5 * xi3 * xi1 + sqrt_rho_p * (mu3 * xi1 - mu1 * xi3)) \
        + particle.dt / (2 * math.pi) * sum31
    I32 = particle.dt * (0.5 * xi3 * xi2 + sqrt_rho_p * (mu3 * xi2 - mu2 * xi3)) \
        + particle.dt / (2 * math.pi) * sum32
    #     Particle positions are updated only after evaluating all terms.
    
    dLon += ax * particle.dt + sigma11 * dWx
    dLat += ay * particle.dt + sigma21 * dWx + sigma22 * dWy
    dDepth += az * particle.dt + sigma31 * dWx + sigma32 * dWy + sigma33 * dWz + (
        (sigma11 * dsigma31dx + sigma31 * dsigma31dz) * I11 + \
        (sigma11 * dsigma32dx + sigma31 * dsigma32dz) * I12 + \
        (sigma11 * dsigma33dx + sigma31 * dsigma33dz) * I13 + \
        (sigma22 * dsigma31dy + sigma32 * dsigma31dz) * I21 + \
        (sigma22 * dsigma32dy + sigma32 * dsigma32dz) * I22 + \
        (sigma22 * dsigma33dy + sigma32 * dsigma33dz) * I23 + \
        sigma33 * dsigma31dz * I31 + \
        sigma33 * dsigma32dz * I32 + \
        sigma33 * dsigma33dz * I33
    )

    
    if (particle.lat + dLat < fieldset.northBound) and (particle.lat + dLat > fieldset.southBound) and (particle.depth + dDepth < fieldset.upperBound) and (particle.depth + dDepth > fieldset.lowerBound):
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
            if particle.lat + dLat < fieldset.southBound:
                particle.lon += dLon
                particle.lat -= dLat
                particle.depth += dDepth 
            elif particle.depth + dDepth > fieldset.upperBound:
                particle.lon += dLon
                particle.lat += dLat
                particle.depth -= dDepth 
            elif particle.depth + dDepth < fieldset.lowerBound:
                particle.lon += dLon
                particle.lat += dLat
                particle.depth -= dDepth 
        else:
            return OperationCode.Repeat
    # fmt: on