import math
import parcels.rng as ParcelsRandom

def M1_3D_BC(particle, fieldset, time):
    if k11 > 0:
        # Independent Wiener increments with zero mean and std of sqrt(dt)
        dWx = ParcelsRandom.normalvariate(0, math.sqrt(math.fabs(particle.dt)))
        dWy = ParcelsRandom.normalvariate(0, math.sqrt(math.fabs(particle.dt)))
        dWz = ParcelsRandom.normalvariate(0, math.sqrt(math.fabs(particle.dt)))
        sqrt2 = math.sqrt(2.)
    
        # Cholesky-Banachiewicz decomposed diffusivity tensor elements
        b11 = math.sqrt(2. * k11)
        b21 = 2. * k21 / b11
        b22 = math.sqrt(math.fabs(2. * k22 - b21**2))
        b31 = 2. * k31 / b11
        b32 = (2. * k32 - b31 * b21) / b22
        b33 = math.sqrt(math.fabs(2. * k33 - b31**2 - b32**2))

        db11dx = dk11dx / b11
        db11dy = dk11dy / b11
        db11dz = dk11dz / b11 

        db21dx = 2. * (b11 * dk21dx - k21 * db11dx) / (b11**2)
        db21dy = 2. * (b11 * dk21dy - k21 * db11dy) / (b11**2)
        db21dz = 2. * (b11 * dk21dz - k21 * db11dz) / (b11**2)

        db22dx = (dk22dx - b21 * db21dx) / b22
        db22dy = (dk22dy - b21 * db21dy) / b22
        db22dz = (dk22dz - b21 * db21dz) / b22

        db31dx = 2. * (b11 * dk31dx - k31 * db11dx) / (b11**2)
        db31dy = 2. * (b11 * dk31dy - k31 * db11dy) / (b11**2)
        db31dz = 2. * (b11 * dk31dz - k31 * db11dz) / (b11**2)

        db32dx = (b22 * (2. * dk32dx - b31 * db21dx - db31dx * b21) - db22dx * (2 * k32 - b31 * b21)) / (b22**2)
        db32dy = (b22 * (2. * dk32dy - b31 * db21dy - db31dy * b21) - db22dy * (2 * k32 - b31 * b21)) / (b22**2)
        db32dz = (b22 * (2. * dk32dz - b31 * db21dz - db31dz * b21) - db22dz * (2 * k32 - b31 * b21)) / (b22**2)

        db33dx = (dk33dx - b31 * db31dx - b32 * db32dx) / b33
        db33dy = (dk33dy - b31 * db31dy - b32 * db32dy) / b33
        db33dz = (dk33dz - b31 * db31dz - b32 * db32dz) / b33

        u, v, w = fieldset.UVW[time, particle.depth, particle.lat, particle.lon]
        ax = u + dk11dx + dk21dy + dk31dz
        ay = v + dk21dx + dk22dy + dk32dz
        az = w + dk31dx + dk32dy + dk33dz

        I11 = 0.5 * (dWx**2 - particle.dt)
        I22 = 0.5 * (dWy**2 - particle.dt)
        I33 = 0.5 * (dWz**2 - particle.dt)

        # Compute double stochastic integrals
        xi1 = 1./math.sqrt(particle.dt) * dWx
        xi2 = 1./math.sqrt(particle.dt) * dWy
        xi3 = 1./math.sqrt(particle.dt) * dWz
        w = 1. # normally p, but already in use in Le Sommer
        rqp = 0.
        sum12 = 0.
        sum13 = 0.
        sum21 = 0.
        sum23 = 0.
        sum31 = 0.
        sum32 = 0.
        mu1 = ParcelsRandom.normalvariate(0, 1)
        mu2 = ParcelsRandom.normalvariate(0, 1)
        mu3 = ParcelsRandom.normalvariate(0, 1)
        while w <= fieldset.expansion_terms:
            zeta1 = ParcelsRandom.normalvariate(0, 1)
            zeta2 = ParcelsRandom.normalvariate(0, 1)
            zeta3 = ParcelsRandom.normalvariate(0, 1)
            eta1 = ParcelsRandom.normalvariate(0, 1)
            eta2 = ParcelsRandom.normalvariate(0, 1)
            eta3 = ParcelsRandom.normalvariate(0, 1)
            rqp += 1. / (w**2)
            sum12 += 1./w * (zeta1 * (sqrt2 * xi2 + eta2) - zeta2 * (sqrt2 * xi1 + eta1))
            sum13 += 1./w * (zeta1 * (sqrt2 * xi3 + eta3) - zeta3 * (sqrt2 * xi1 + eta1))
            sum21 += 1./w * (zeta2 * (sqrt2 * xi1 + eta1) - zeta1 * (sqrt2 * xi2 + eta2))
            sum23 += 1./w * (zeta2 * (sqrt2 * xi3 + eta3) - zeta3 * (sqrt2 * xi2 + eta2))
            sum31 += 1./w * (zeta3 * (sqrt2 * xi1 + eta1) - zeta1 * (sqrt2 * xi3 + eta3))
            sum32 += 1./w * (zeta3 * (sqrt2 * xi2 + eta2) - zeta2 * (sqrt2 * xi3 + eta3))
            w += 1
        sqrt_rho_p = math.sqrt(1./12 - 1. / (2. * math.pi**2) * rqp)
        I12 = particle.dt * (0.5 * xi1 * xi2 + sqrt_rho_p * (mu1 * xi2 - mu2 * xi1)) \
            + particle.dt / (2. * math.pi) * sum12
        I13 = particle.dt * (0.5 * xi1 * xi3 + sqrt_rho_p * (mu1 * xi3 - mu3 * xi1)) \
            + particle.dt / (2. * math.pi) * sum13
        I21 = particle.dt * (0.5 * xi2 * xi1 + sqrt_rho_p * (mu2 * xi1 - mu1 * xi2)) \
            + particle.dt / (2. * math.pi) * sum21
        I23 = particle.dt * (0.5 * xi2 * xi3 + sqrt_rho_p * (mu2 * xi3 - mu3 * xi2)) \
            + particle.dt / (2. * math.pi) * sum23
        I31 = particle.dt * (0.5 * xi3 * xi1 + sqrt_rho_p * (mu3 * xi1 - mu1 * xi3)) \
            + particle.dt / (2. * math.pi) * sum31
        I32 = particle.dt * (0.5 * xi3 * xi2 + sqrt_rho_p * (mu3 * xi2 - mu2 * xi3)) \
            + particle.dt / (2. * math.pi) * sum32
        #     Particle positions are updated only after evaluating all terms.
        particle.lon += ax * particle.dt + b11 * dWx + (
            (b11 * db11dx + b21 * db11dy + b31 * db11dz) * I11 + \
            (b22 * db11dy + b32 * db11dz) * I21 + \
            (b33 * db11dz) * I31
        )
        particle.lat += ay * particle.dt + b21 * dWx + b22 * dWy + (
            (b11 * db21dx + b21 * db21dy + b31 * db21dz) * I11 + \
            (b11 * db22dx + b21 * db22dy + b31 * db22dz) * I12 + \
            (b22 * db21dy + b32 * db21dz) * I21 + \
            (b22 * db22dy + b32 * db22dz) * I22 + \
            (b33 * db21dz) * I31 + \
            (b33 * db22dz) * I32
        )
        particle.depth += az * particle.dt + b31 * dWx + b32 * dWy + b33 * dWz + (
            (b11 * db31dx + b21 * db31dy + b31 * db31dz) * I11 + \
            (b11 * db32dx + b21 * db32dy + b31 * db32dz) * I12 + \
            (b11 * db33dx + b21 * db33dy + b31 * db33dz) * I13 + \
            (b22 * db31dy + b32 * db31dz) * I21 + \
            (b22 * db32dy + b32 * db32dz) * I22 + \
            (b22 * db33dy + b32 * db33dz) * I23 + \
            (b33 * db31dz) * I31 + \
            (b33 * db32dz) * I32 + \
            (b33 * db33dz) * I33
        )

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