import math
import parcels.rng as ParcelsRandom

def K_Redi(particle, fieldset, time):
    kappa = fieldset.kappaval
    eps = fieldset.epsilon
    
    denom = 1 / (drhodx**2 + drhody**2 + drhodz**2)
    denom_recip = drhodx**2 + drhody**2 + drhodz**2  
    
    k11 = kappa * denom * (eps * drhodx**2 + (drhody**2 + drhodz**2))
    k21 = kappa * denom * (eps - 1) * drhodx * drhody
    k22 = kappa * denom * (eps * drhody**2 + (drhodx**2 + drhodz**2))
    k31 = kappa * denom * (eps - 1) * drhodx * drhodz 
    k32 = kappa * denom * (eps - 1) * drhody * drhodz
    k33 = kappa * denom * (eps * drhodz**2 + (drhodx**2 + drhody**2))
    
    Sx = - drhodx / drhodz
    Sy = - drhody / drhodz
    Sabs = math.sqrt(Sx**2 + Sy**2)
    Sabs2 = Sx**2 + Sy**2
    
    if Sabs < fieldset.Sc - 3 * fieldset.Sd:
        taper1 = 1.
    elif Sabs > fieldset.Sc + 3 * fieldset.Sd:
        taper1 = 0.
    else:
        taper1 = 0.5 * (1. + math.tanh((fieldset.Sc - Sabs)/fieldset.Sd))
    if fieldset.boundaryMask[particle] < 0.95:
        taper2 = 0
    else:
        taper2 = 1.
    
    dk11dx = taper1 * taper2 * kappa * denom**2 * 2 * (denom_recip * (eps * drhodx*drhodxx + (drhody*drhodyx + drhodz*drhodzx)) - (eps * drhodx**2 + (drhody**2 + drhodz**2)) * (drhodx*drhodxx + drhody*drhodyx + drhodz*drhodzx))
    dk11dy = taper1 * taper2 * kappa * denom**2 * 2 * (denom_recip * (eps * drhodx*drhodxy + (drhody*drhodyy + drhodz*drhodzy)) - (eps * drhodx**2 + (drhody**2 + drhodz**2)) * (drhodx*drhodxy + drhody*drhodyy + drhodz*drhodzy))
    dk11dz = taper1 * taper2 * kappa * denom**2 * 2 * (denom_recip * (eps * drhodx*drhodxz + (drhody*drhodyz + drhodz*drhodzz)) - (eps * drhodx**2 + (drhody**2 + drhodz**2)) * (drhodx*drhodxz + drhody*drhodyz + drhodz*drhodzz))
    
    dk21dx = taper1 * taper2 * kappa * denom**2 * (denom_recip * (eps - 1) * (drhodxx * drhody + drhodx * drhodyx) - 2 * (eps - 1) * drhodx * drhody * (drhodx*drhodxx + drhody*drhodyx + drhodz*drhodzx))
    dk21dy = taper1 * taper2 * kappa * denom**2 * (denom_recip * (eps - 1) * (drhodxy * drhody + drhodx * drhodyy) - 2 * (eps - 1) * drhodx * drhody * (drhodx*drhodxy + drhody*drhodyy + drhodz*drhodzy))
    dk21dz = taper1 * taper2 * kappa * denom**2 * (denom_recip * (eps - 1) * (drhodxz * drhody + drhodx * drhodyz) - 2 * (eps - 1) * drhodx * drhody * (drhodx*drhodxz + drhody*drhodyz + drhodz*drhodzz))
    
    dk22dx = taper1 * taper2 * kappa * denom**2 * 2 * (denom_recip * (eps * drhody*drhodyx + (drhodx*drhodxx + drhodz*drhodzx)) - (eps * drhody**2 + (drhodx**2 + drhodz**2)) * (drhodx*drhodxx + drhody*drhodyx + drhodz*drhodzx))
    dk22dy = taper1 * taper2 * kappa * denom**2 * 2 * (denom_recip * (eps * drhody*drhodyy + (drhodx*drhodxy + drhodz*drhodzy)) - (eps * drhody**2 + (drhodx**2 + drhodz**2)) * (drhodx*drhodxy + drhody*drhodyy + drhodz*drhodzy))
    dk22dz = taper1 * taper2 * kappa * denom**2 * 2 * (denom_recip * (eps * drhody*drhodyz + (drhodx*drhodxz + drhodz*drhodzz)) - (eps * drhody**2 + (drhodx**2 + drhodz**2)) * (drhodx*drhodxz + drhody*drhodyz + drhodz*drhodzz))

    dk31dx = taper1 * taper2 * kappa * denom**2 * (denom_recip * (eps - 1) * (drhodxx * drhodz + drhodx * drhodzx) - 2 * (eps - 1) * drhodx * drhodz * (drhodx*drhodxx + drhody*drhodyx + drhodz*drhodzx))
    dk31dy = taper1 * taper2 * kappa * denom**2 * (denom_recip * (eps - 1) * (drhodxy * drhodz + drhodx * drhodzy) - 2 * (eps - 1) * drhodx * drhodz * (drhodx*drhodxy + drhody*drhodyy + drhodz*drhodzy))
    dk31dz = taper1 * taper2 * kappa * denom**2 * (denom_recip * (eps - 1) * (drhodxz * drhodz + drhodx * drhodzz) - 2 * (eps - 1) * drhodx * drhodz * (drhodx*drhodxz + drhody*drhodyz + drhodz*drhodzz))
    
    dk32dx = taper1 * taper2 * kappa * denom**2 * (denom_recip * (eps - 1) * (drhodyx * drhodz + drhody * drhodzx) - 2 * (eps - 1) * drhody * drhodz * (drhodx*drhodxx + drhody*drhodyx + drhodz*drhodzx))
    dk32dy = taper1 * taper2 * kappa * denom**2 * (denom_recip * (eps - 1) * (drhodyy * drhodz + drhody * drhodzy) - 2 * (eps - 1) * drhody * drhodz * (drhodx*drhodxy + drhody*drhodyy + drhodz*drhodzy))
    dk32dz = taper1 * taper2 * kappa * denom**2 * (denom_recip * (eps - 1) * (drhodyz * drhodz + drhody * drhodzz) - 2 * (eps - 1) * drhody * drhodz * (drhodx*drhodxz + drhody*drhodyz + drhodz*drhodzz))
    
    dk33dx = taper1 * taper2 * kappa * denom**2 * 2 * (denom_recip * (eps * drhodz*drhodzx + (drhodx*drhodxx + drhody*drhodyx)) - (eps * drhodz**2 + (drhodx**2 + drhody**2)) * (drhodx*drhodxx + drhody*drhodyx + drhodz*drhodzx))
    dk33dy = taper1 * taper2 * kappa * denom**2 * 2 * (denom_recip * (eps * drhodz*drhodzy + (drhodx*drhodxy + drhody*drhodyy)) - (eps * drhodz**2 + (drhodx**2 + drhody**2)) * (drhodx*drhodxy + drhody*drhodyy + drhodz*drhodzy))
    dk33dz = taper1 * taper2 * kappa * denom**2 * 2 * (denom_recip * (eps * drhodz*drhodzz + (drhodx*drhodxz + drhody*drhodyz)) - (eps * drhodz**2 + (drhodx**2 + drhody**2)) * (drhodx*drhodxz + drhody*drhodyz + drhodz*drhodzz))