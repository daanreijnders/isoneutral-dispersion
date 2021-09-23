import math
import parcels.rng as ParcelsRandom

def Markov1_3D_BC_taper(particle, fieldset, time):
    # Independent Wiener increments with zero mean and std of sqrt(dt)
    dWx = ParcelsRandom.normalvariate(0, math.sqrt(math.fabs(particle.dt)))
    dWy = ParcelsRandom.normalvariate(0, math.sqrt(math.fabs(particle.dt)))
    dWz = ParcelsRandom.normalvariate(0, math.sqrt(math.fabs(particle.dt)))

    # Model parameters
    Tl = fieldset.TL
    nusq = fieldset.nusquared
    eps = fieldset.epsilon
    eta = fieldset.etaval

    # Velocities
    (u1, v1, w1) = fieldset.UVW[time, particle.depth, particle.lat, particle.lon]

    # (double) derivatives in density field are assumed available 
    # Computation of taper function
    Sx = - drhodx / drhodz
    Sy = - drhody / drhodz
    Sabs = math.sqrt(Sx**2 + Sy**2)
    # Tapering based on Danabasoglu and McWilliams, J. Clim. 1995
    # Tapering is discrete outside of a given range, to prevent
    # it from having its own decay timescale
    if Sabs < fieldset.Sc - 3. * fieldset.Sd:
        taper1 = 1.
    elif Sabs > fieldset.Sc + 3. * fieldset.Sd:
        taper1 = 0.
    else:
        taper1 = 0.5 * (1. + math.tanh((fieldset.Sc - Sabs)/fieldset.Sd))
    # Near boundaries, drop to 0
    if fieldset.boundaryMask[particle] < 0.95:
        taper2 = 0.
    else:
        taper2 = 1.
    
    # Computation of inverse Lagrangian timescale tensor, (inverse) velocity variance tensor (and derivative), and their product.
    A = (drhodx**2 + drhody**2 + drhodz**2)
    Arcp = 1. / (drhodx**2 + drhody**2 + drhodz**2)

    thetainv11 = Arcp / Tl * (1./eps * drhodx**2 + (drhody**2 + drhodz**2))
    thetainv12 = Arcp / Tl * (1/eps - 1) * drhodx * drhody
    thetainv13 = Arcp / Tl * (1./eps - 1) * drhodx * drhodz
    thetainv21 = Arcp / Tl * (1./eps - 1) * drhodx * drhody
    thetainv22 = Arcp / Tl * (1./eps * drhody**2 + (drhodx**2 + drhodz**2))
    thetainv23 = Arcp / Tl * (1./eps - 1) * drhody * drhodz
    thetainv31 = Arcp / Tl * (1./eps - 1) * drhodx * drhodz
    thetainv32 = Arcp / Tl * (1./eps - 1) * drhody * drhodz
    thetainv33 = Arcp / Tl * (1./eps * drhodz**2 + (drhodx**2 + drhody**2))

    sigma11 = nusq * Arcp * (eta * drhodx**2 + (drhody**2 + drhodz**2))
    sigma12 = nusq * Arcp * (eta - 1) * drhodx * drhody
    sigma13 = nusq * Arcp * (eta - 1) * drhodx * drhodz 
    sigma21 = nusq * Arcp * (eta - 1) * drhodx * drhody
    sigma22 = nusq * Arcp * (eta * drhody**2 + (drhodx**2 + drhodz**2))
    sigma23 = nusq * Arcp * (eta - 1) * drhody * drhodz
    sigma31 = nusq * Arcp * (eta - 1) * drhodx * drhodz 
    sigma32 = nusq * Arcp * (eta - 1) * drhody * drhodz
    sigma33 = nusq * Arcp * (eta * drhodz**2 + (drhodx**2 + drhody**2))
    
    dsigma11dx = nusq * Arcp**2 * 2 * (A * (eta * drhodx*drhodxx + (drhody*drhodyx + drhodz*drhodzx)) - (eta * drhodx**2 + (drhody**2 + drhodz**2)) * (drhodx*drhodxx + drhody*drhodyx + drhodz*drhodzx))
    dsigma21dy = nusq * Arcp**2 * (A * (eta - 1) * (drhodxy * drhody + drhodx * drhodyy) - 2 * (eta - 1) * drhodx * drhody * (drhodx*drhodxy + drhody*drhodyy + drhodz*drhodzy))
    dsigma31dz = nusq * Arcp**2 * (A * (eta - 1) * (drhodxz * drhodz + drhodx * drhodzz) - 2 * (eta - 1) * drhodx * drhodz * (drhodx*drhodxz + drhody*drhodyz + drhodz*drhodzz))
    dsigma21dx = nusq * Arcp**2 * (A * (eta - 1) * (drhodxx * drhody + drhodx * drhodyx) - 2 * (eta - 1) * drhodx * drhody * (drhodx*drhodxx + drhody*drhodyx + drhodz*drhodzx))
    dsigma22dy = nusq * Arcp**2 * 2 * (A * (eta * drhody*drhodyy + (drhodx*drhodxy + drhodz*drhodzy)) - (eta * drhody**2 + (drhodx**2 + drhodz**2)) * (drhodx*drhodxy + drhody*drhodyy + drhodz*drhodzy))
    dsigma32dz = nusq * Arcp**2 * (A * (eta - 1) * (drhodyz * drhodz + drhody * drhodzz) - 2 * (eta - 1) * drhody * drhodz * (drhodx*drhodxz + drhody*drhodyz + drhodz*drhodzz))
    dsigma31dx = nusq * Arcp**2 * (A * (eta - 1) * (drhodxx * drhodz + drhodx * drhodzx) - 2 * (eta - 1) * drhodx * drhodz * (drhodx*drhodxx + drhody*drhodyx + drhodz*drhodzx))
    dsigma32dy = nusq * Arcp**2 * (A * (eta - 1) * (drhodyy * drhodz + drhody * drhodzy) - 2 * (eta - 1) * drhody * drhodz * (drhodx*drhodxy + drhody*drhodyy + drhodz*drhodzy))
    dsigma33dz = nusq * Arcp**2 * 2 * (A * (eta * drhodz*drhodzz + (drhodx*drhodxz + drhody*drhodyz)) - (eta * drhodz**2 + (drhodx**2 + drhody**2)) * (drhodx*drhodxz + drhody*drhodyz + drhodz*drhodzz))

    dsigmainv11dx = Arcp**2 / nusq * 2 * (A * (1/eta * drhodx*drhodxx + (drhody*drhodyx + drhodz*drhodzx)) - (1/eta * drhodx**2 + (drhody**2 + drhodz**2)) * (drhodx*drhodxx + drhody*drhodyx + drhodz*drhodzx))
    dsigmainv12dx = Arcp**2 / nusq * (A * (1/eta - 1) * (drhodxx * drhody + drhodx * drhodyx) - 2 * (1/eta - 1) * drhodx * drhody * (drhodx*drhodxx + drhody*drhodyx + drhodz*drhodzx))
    dsigmainv13dx = Arcp**2 / nusq * (A * (1/eta - 1) * (drhodxx * drhodz + drhodx * drhodzx) - 2 * (1/eta - 1) * drhodx * drhodz * (drhodx*drhodxx + drhody*drhodyx + drhodz*drhodzx))
    dsigmainv21dx = Arcp**2 / nusq * (A * (1/eta - 1) * (drhodxx * drhody + drhodx * drhodyx) - 2 * (1/eta - 1) * drhodx * drhody * (drhodx*drhodxx + drhody*drhodyx + drhodz*drhodzx))
    dsigmainv22dx = Arcp**2 / nusq * 2 * (A * (1/eta * drhody*drhodyx + (drhodx*drhodxx + drhodz*drhodzx)) - (1/eta * drhody**2 + (drhodx**2 + drhodz**2)) * (drhodx*drhodxx + drhody*drhodyx + drhodz*drhodzx))
    dsigmainv23dx = Arcp**2 / nusq * (A * (1/eta - 1) * (drhodyx * drhodz + drhody * drhodzx) - 2 * (1/eta - 1) * drhody * drhodz * (drhodx*drhodxx + drhody*drhodyx + drhodz*drhodzx))
    dsigmainv31dx = Arcp**2 / nusq * (A * (1/eta - 1) * (drhodxx * drhodz + drhodx * drhodzx) - 2 * (1/eta - 1) * drhodx * drhodz * (drhodx*drhodxx + drhody*drhodyx + drhodz*drhodzx))
    dsigmainv32dx = Arcp**2 / nusq * (A * (1/eta - 1) * (drhodyx * drhodz + drhody * drhodzx) - 2 * (1/eta - 1) * drhody * drhodz * (drhodx*drhodxx + drhody*drhodyx + drhodz*drhodzx))
    dsigmainv33dx = Arcp**2 / nusq * 2 * (A * (1/eta * drhodz*drhodzx + (drhodx*drhodxx + drhody*drhodyx)) - (1/eta * drhodz**2 + (drhodx**2 + drhody**2)) * (drhodx*drhodxx + drhody*drhodyx + drhodz*drhodzx))
    
    dsigmainv11dy = Arcp**2 / nusq * 2 * (A * (1/eta * drhodx*drhodxy + (drhody*drhodyy + drhodz*drhodzy)) - (1/eta * drhodx**2 + (drhody**2 + drhodz**2)) * (drhodx*drhodxy + drhody*drhodyy + drhodz*drhodzy))
    dsigmainv12dy = Arcp**2 / nusq * (A * (1/eta - 1) * (drhodxy * drhody + drhodx * drhodyy) - 2 * (1/eta - 1) * drhodx * drhody * (drhodx*drhodxy + drhody*drhodyy + drhodz*drhodzy))
    dsigmainv13dy = Arcp**2 / nusq * (A * (1/eta - 1) * (drhodxy * drhodz + drhodx * drhodzy) - 2 * (1/eta - 1) * drhodx * drhodz * (drhodx*drhodxy + drhody*drhodyy + drhodz*drhodzy))
    dsigmainv21dy = Arcp**2 / nusq * (A * (1/eta - 1) * (drhodxy * drhody + drhodx * drhodyy) - 2 * (1/eta - 1) * drhodx * drhody * (drhodx*drhodxy + drhody*drhodyy + drhodz*drhodzy))
    dsigmainv22dy = Arcp**2 / nusq * 2 * (A * (1/eta * drhody*drhodyy + (drhodx*drhodxy + drhodz*drhodzy)) - (1/eta * drhody**2 + (drhodx**2 + drhodz**2)) * (drhodx*drhodxy + drhody*drhodyy + drhodz*drhodzy))
    dsigmainv23dy = Arcp**2 / nusq * (A * (1/eta - 1) * (drhodyy * drhodz + drhody * drhodzy) - 2 * (1/eta - 1) * drhody * drhodz * (drhodx*drhodxy + drhody*drhodyy + drhodz*drhodzy))
    dsigmainv31dy = Arcp**2 / nusq * (A * (1/eta - 1) * (drhodxy * drhodz + drhodx * drhodzy) - 2 * (1/eta - 1) * drhodx * drhodz * (drhodx*drhodxy + drhody*drhodyy + drhodz*drhodzy))
    dsigmainv32dy = Arcp**2 / nusq * (A * (1/eta - 1) * (drhodyy * drhodz + drhody * drhodzy) - 2 * (1/eta - 1) * drhody * drhodz * (drhodx*drhodxy + drhody*drhodyy + drhodz*drhodzy))
    dsigmainv33dy = Arcp**2 / nusq * 2 * (A * (1/eta * drhodz*drhodzy + (drhodx*drhodxy + drhody*drhodyy)) - (1/eta * drhodz**2 + (drhodx**2 + drhody**2)) * (drhodx*drhodxy + drhody*drhodyy + drhodz*drhodzy))
    
    dsigmainv11dz = Arcp**2 / nusq * 2 * (A * (1/eta * drhodx*drhodxz + (drhody*drhodyz + drhodz*drhodzz)) - (1/eta * drhodx**2 + (drhody**2 + drhodz**2)) * (drhodx*drhodxz + drhody*drhodyz + drhodz*drhodzz))
    dsigmainv12dz = Arcp**2 / nusq * (A * (1/eta - 1) * (drhodxz * drhody + drhodx * drhodyz) - 2 * (1/eta - 1) * drhodx * drhody * (drhodx*drhodxz + drhody*drhodyz + drhodz*drhodzz))
    dsigmainv13dz = Arcp**2 / nusq * (A * (1/eta - 1) * (drhodxz * drhodz + drhodx * drhodzz) - 2 * (1/eta - 1) * drhodx * drhodz * (drhodx*drhodxz + drhody*drhodyz + drhodz*drhodzz))
    dsigmainv21dz = Arcp**2 / nusq * (A * (1/eta - 1) * (drhodxz * drhody + drhodx * drhodyz) - 2 * (1/eta - 1) * drhodx * drhody * (drhodx*drhodxz + drhody*drhodyz + drhodz*drhodzz))
    dsigmainv22dz = Arcp**2 / nusq * 2 * (A * (1/eta * drhody*drhodyz + (drhodx*drhodxz + drhodz*drhodzz)) - (1/eta * drhody**2 + (drhodx**2 + drhodz**2)) * (drhodx*drhodxz + drhody*drhodyz + drhodz*drhodzz))
    dsigmainv23dz = Arcp**2 / nusq * (A * (1/eta - 1) * (drhodyz * drhodz + drhody * drhodzz) - 2 * (1/eta - 1) * drhody * drhodz * (drhodx*drhodxz + drhody*drhodyz + drhodz*drhodzz))
    dsigmainv31dz = Arcp**2 / nusq * (A * (1/eta - 1) * (drhodxz * drhodz + drhodx * drhodzz) - 2 * (1/eta - 1) * drhodx * drhodz * (drhodx*drhodxz + drhody*drhodyz + drhodz*drhodzz))
    dsigmainv32dz = Arcp**2 / nusq * (A * (1/eta - 1) * (drhodyz * drhodz + drhody * drhodzz) - 2 * (1/eta - 1) * drhody * drhodz * (drhodx*drhodxz + drhody*drhodyz + drhodz*drhodzz))
    dsigmainv33dz = Arcp**2 / nusq * 2 * (A * (1/eta * drhodz*drhodzz + (drhodx*drhodxz + drhody*drhodyz)) - (1/eta * drhodz**2 + (drhodx**2 + drhody**2)) * (drhodx*drhodxz + drhody*drhodyz + drhodz*drhodzz))
                                      
    SigmaTheta11 = 2 * nusq / Tl * Arcp * (eta/eps * drhodx**2 + (drhody**2 + drhodz**2))
    SigmaTheta21 = 2 * nusq / Tl * Arcp * (eta/eps - 1) * drhodx * drhody
    SigmaTheta22 = 2 * nusq / Tl * Arcp * (eta/eps * drhody**2 + (drhodx**2 + drhodz**2))
    SigmaTheta31 = 2 * nusq / Tl * Arcp * (eta/eps - 1) * drhodx * drhodz
    SigmaTheta32 = 2 * nusq / Tl * Arcp * (eta/eps - 1) * drhody * drhodz
    SigmaTheta33 = 2 * nusq / Tl * Arcp * (eta/eps * drhodz**2 + (drhodx**2 + drhody**2))
    
    # Cholesky-Banachiewicz decomposed diffusivity tensor elements
    b11 = math.sqrt(SigmaTheta11)
    b21 = SigmaTheta21 / b11
    b22 = math.sqrt(math.fabs(SigmaTheta22 - b21**2))
    b31 = SigmaTheta31 / b11
    b32 = (SigmaTheta32 - b31*b21)/b22
    b33 = math.sqrt(math.fabs(SigmaTheta33 - b31**2 - b32**2)) 
    
    # Drift correction terms
    a1x = 0.5 * (dsigma11dx + dsigma21dy + dsigma31dz)
    a1y = 0.5 * (dsigma21dx + dsigma22dy + dsigma32dz)
    a1z = 0.5 * (dsigma31dx + dsigma32dy + dsigma33dz)
    
    a2x = sigma11/2. * (u1 + particle.u_prime) * dsigmainv11dx * particle.u_prime + \
           sigma12/2. * (u1 + particle.u_prime) * dsigmainv12dx * particle.u_prime + \
           sigma13/2. * (u1 + particle.u_prime) * dsigmainv13dx * particle.u_prime + \
           sigma11/2. * (v1 + particle.v_prime) * dsigmainv11dy * particle.u_prime + \
           sigma12/2. * (v1 + particle.v_prime) * dsigmainv12dy * particle.u_prime + \
           sigma13/2. * (v1 + particle.v_prime) * dsigmainv13dy * particle.u_prime + \
           sigma11/2. * (w1 + particle.w_prime) * dsigmainv11dz * particle.u_prime + \
           sigma12/2. * (w1 + particle.w_prime) * dsigmainv12dz * particle.u_prime + \
           sigma13/2. * (w1 + particle.w_prime) * dsigmainv13dz * particle.u_prime + \
           sigma11/2. * (u1 + particle.u_prime) * dsigmainv21dx * particle.v_prime + \
           sigma12/2. * (u1 + particle.u_prime) * dsigmainv22dx * particle.v_prime + \
           sigma13/2. * (u1 + particle.u_prime) * dsigmainv23dx * particle.v_prime + \
           sigma11/2. * (v1 + particle.v_prime) * dsigmainv21dy * particle.v_prime + \
           sigma12/2. * (v1 + particle.v_prime) * dsigmainv22dy * particle.v_prime + \
           sigma13/2. * (v1 + particle.v_prime) * dsigmainv23dy * particle.v_prime + \
           sigma11/2. * (w1 + particle.w_prime) * dsigmainv21dz * particle.v_prime + \
           sigma12/2. * (w1 + particle.w_prime) * dsigmainv22dz * particle.v_prime + \
           sigma13/2. * (w1 + particle.w_prime) * dsigmainv23dz * particle.v_prime + \
           sigma11/2. * (u1 + particle.u_prime) * dsigmainv31dx * particle.w_prime + \
           sigma12/2. * (u1 + particle.u_prime) * dsigmainv32dx * particle.w_prime + \
           sigma13/2. * (u1 + particle.u_prime) * dsigmainv33dx * particle.w_prime + \
           sigma11/2. * (v1 + particle.v_prime) * dsigmainv31dy * particle.w_prime + \
           sigma12/2. * (v1 + particle.v_prime) * dsigmainv32dy * particle.w_prime + \
           sigma13/2. * (v1 + particle.v_prime) * dsigmainv33dy * particle.w_prime + \
           sigma11/2. * (w1 + particle.w_prime) * dsigmainv31dz * particle.w_prime + \
           sigma12/2. * (w1 + particle.w_prime) * dsigmainv32dz * particle.w_prime + \
           sigma13/2. * (w1 + particle.w_prime) * dsigmainv33dz * particle.w_prime
    
    a2y = sigma21/2. * (u1 + particle.u_prime) * dsigmainv11dx * particle.u_prime + \
           sigma22/2. * (u1 + particle.u_prime) * dsigmainv12dx * particle.u_prime + \
           sigma23/2. * (u1 + particle.u_prime) * dsigmainv13dx * particle.u_prime + \
           sigma21/2. * (v1 + particle.v_prime) * dsigmainv11dy * particle.u_prime + \
           sigma22/2. * (v1 + particle.v_prime) * dsigmainv12dy * particle.u_prime + \
           sigma23/2. * (v1 + particle.v_prime) * dsigmainv13dy * particle.u_prime + \
           sigma21/2. * (w1 + particle.w_prime) * dsigmainv11dz * particle.u_prime + \
           sigma22/2. * (w1 + particle.w_prime) * dsigmainv12dz * particle.u_prime + \
           sigma23/2. * (w1 + particle.w_prime) * dsigmainv13dz * particle.u_prime + \
           sigma21/2. * (u1 + particle.u_prime) * dsigmainv21dx * particle.v_prime + \
           sigma22/2. * (u1 + particle.u_prime) * dsigmainv22dx * particle.v_prime + \
           sigma23/2. * (u1 + particle.u_prime) * dsigmainv23dx * particle.v_prime + \
           sigma21/2. * (v1 + particle.v_prime) * dsigmainv21dy * particle.v_prime + \
           sigma22/2. * (v1 + particle.v_prime) * dsigmainv22dy * particle.v_prime + \
           sigma23/2. * (v1 + particle.v_prime) * dsigmainv23dy * particle.v_prime + \
           sigma21/2. * (w1 + particle.w_prime) * dsigmainv21dz * particle.v_prime + \
           sigma22/2. * (w1 + particle.w_prime) * dsigmainv22dz * particle.v_prime + \
           sigma23/2. * (w1 + particle.w_prime) * dsigmainv23dz * particle.v_prime + \
           sigma21/2. * (u1 + particle.u_prime) * dsigmainv31dx * particle.w_prime + \
           sigma22/2. * (u1 + particle.u_prime) * dsigmainv32dx * particle.w_prime + \
           sigma23/2. * (u1 + particle.u_prime) * dsigmainv33dx * particle.w_prime + \
           sigma21/2. * (v1 + particle.v_prime) * dsigmainv31dy * particle.w_prime + \
           sigma22/2. * (v1 + particle.v_prime) * dsigmainv32dy * particle.w_prime + \
           sigma23/2. * (v1 + particle.v_prime) * dsigmainv33dy * particle.w_prime + \
           sigma21/2. * (w1 + particle.w_prime) * dsigmainv31dz * particle.w_prime + \
           sigma22/2. * (w1 + particle.w_prime) * dsigmainv32dz * particle.w_prime + \
           sigma23/2. * (w1 + particle.w_prime) * dsigmainv33dz * particle.w_prime
    
    a2z = sigma31/2. * (u1 + particle.u_prime) * dsigmainv11dx * particle.u_prime + \
           sigma32/2. * (u1 + particle.u_prime) * dsigmainv12dx * particle.u_prime + \
           sigma33/2. * (u1 + particle.u_prime) * dsigmainv13dx * particle.u_prime + \
           sigma31/2. * (v1 + particle.v_prime) * dsigmainv11dy * particle.u_prime + \
           sigma32/2. * (v1 + particle.v_prime) * dsigmainv12dy * particle.u_prime + \
           sigma33/2. * (v1 + particle.v_prime) * dsigmainv13dy * particle.u_prime + \
           sigma31/2. * (w1 + particle.w_prime) * dsigmainv11dz * particle.u_prime + \
           sigma32/2. * (w1 + particle.w_prime) * dsigmainv12dz * particle.u_prime + \
           sigma33/2. * (w1 + particle.w_prime) * dsigmainv13dz * particle.u_prime + \
           sigma31/2. * (u1 + particle.u_prime) * dsigmainv21dx * particle.v_prime + \
           sigma32/2. * (u1 + particle.u_prime) * dsigmainv22dx * particle.v_prime + \
           sigma33/2. * (u1 + particle.u_prime) * dsigmainv23dx * particle.v_prime + \
           sigma31/2. * (v1 + particle.v_prime) * dsigmainv21dy * particle.v_prime + \
           sigma32/2. * (v1 + particle.v_prime) * dsigmainv22dy * particle.v_prime + \
           sigma33/2. * (v1 + particle.v_prime) * dsigmainv23dy * particle.v_prime + \
           sigma31/2. * (w1 + particle.w_prime) * dsigmainv21dz * particle.v_prime + \
           sigma32/2. * (w1 + particle.w_prime) * dsigmainv22dz * particle.v_prime + \
           sigma33/2. * (w1 + particle.w_prime) * dsigmainv23dz * particle.v_prime + \
           sigma31/2. * (u1 + particle.u_prime) * dsigmainv31dx * particle.w_prime + \
           sigma32/2. * (u1 + particle.u_prime) * dsigmainv32dx * particle.w_prime + \
           sigma33/2. * (u1 + particle.u_prime) * dsigmainv33dx * particle.w_prime + \
           sigma31/2. * (v1 + particle.v_prime) * dsigmainv31dy * particle.w_prime + \
           sigma32/2. * (v1 + particle.v_prime) * dsigmainv32dy * particle.w_prime + \
           sigma33/2. * (v1 + particle.v_prime) * dsigmainv33dy * particle.w_prime + \
           sigma31/2. * (w1 + particle.w_prime) * dsigmainv31dz * particle.w_prime + \
           sigma32/2. * (w1 + particle.w_prime) * dsigmainv32dz * particle.w_prime + \
           sigma33/2. * (w1 + particle.w_prime) * dsigmainv33dz * particle.w_prime

    # compute change in velocity perturbation 
    du_prime = (-(thetainv11 * particle.u_prime + thetainv12 * particle.v_prime + thetainv13 * particle.w_prime) + a1x - a2x) * particle.dt + b11 * dWx
    dv_prime = (-(thetainv21 * particle.u_prime + thetainv22 * particle.v_prime + thetainv23 * particle.w_prime) + a1y - a2y) * particle.dt + b21 * dWx + b22 * dWy
    dw_prime = (-(thetainv31 * particle.u_prime + thetainv32 * particle.v_prime + thetainv33 * particle.w_prime) + a1z - a2z) * particle.dt + b31 * dWx + b32 * dWy + b33 * dWz
    
    # Update velocity perturbation and apply taper
    particle.u_prime = taper1 * taper2 * (particle.u_prime + du_prime)
    particle.v_prime = taper1 * taper2 * (particle.v_prime + dv_prime)
    particle.w_prime = taper1 * taper2 * (particle.w_prime + dw_prime)
    
    # Prevent excursions in the velocity perturbations larger than 4*nu (standard deviation in velocity variance)
    if math.sqrt(particle.u_prime**2 + particle.v_prime**2 + particle.w_prime**2) < 4 * math.sqrt(nusq):
        taper3 = 1.
    else:
        taper3 = 0.
    particle.u_prime = taper3 * particle.u_prime
    particle.v_prime = taper3 * particle.v_prime
    particle.w_prime = taper3 * particle.w_prime
    
    # RK4 advection of the mean flow; Euler for the perturbation.
    lon1 = particle.lon + u1*.5*particle.dt
    lat1 = particle.lat + v1*.5*particle.dt
    dep1 = particle.depth + w1*.5*particle.dt
    (u2, v2, w2) = fieldset.UVW[time + .5 * particle.dt, dep1, lat1, lon1, particle]
    lon2 = particle.lon + u2*.5*particle.dt
    lat2 = particle.lat + v2*.5*particle.dt
    dep2 = particle.depth + w2*.5*particle.dt
    (u3, v3, w3) = fieldset.UVW[time + .5 * particle.dt, dep2, lat2, lon2, particle]
    lon3 = particle.lon + u3*particle.dt
    lat3 = particle.lat + v3*particle.dt
    dep3 = particle.depth + w3*particle.dt
    (u4, v4, w4) = fieldset.UVW[time + particle.dt, dep3, lat3, lon3, particle]
    
    # Compute displacement
    particle.lon += ((u1 + 2*u2 + 2*u3 + u4) / 6. + particle.u_prime) * particle.dt
    particle.lat += ((v1 + 2*v2 + 2*v3 + v4) / 6. + particle.v_prime) * particle.dt
    particle.depth += ((w1 + 2*w2 + 2*w3 + w4) / 6. + particle.w_prime) * particle.dt
    
    # Boundary condition
    if particle.lat > fieldset.northBound:
        particle.lat = 2. * fieldset.northBound - particle.lat
        particle.v_prime = - particle.v_prime
    elif particle.lat < fieldset.southBound:
        particle.lat = 2. * fieldset.southBound - particle.lat
        particle.v_prime = - particle.v_prime
    if particle.depth > fieldset.upperBound:
        particle.depth = 2. * fieldset.upperBound - particle.depth
        particle.w_prime = - particle.w_prime
    elif particle.depth < fieldset.lowerBound:
        particle.depth = 2. * fieldset.lowerBound - particle.depth
        particle.w_prime = - particle.w_prime

