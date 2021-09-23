import math
import parcels.rng as ParcelsRandom

def K_Le_Sommer(particle, fieldset, time):
    hsq = fieldset.hsquared
    eps = fieldset.epsilon
    kappa = fieldset.kappaval
    
    Sx = - drhodx / drhodz
    Sy = - drhody / drhodz
    Sabs = math.sqrt(Sx**2 + Sy**2)
    
    dSxdx = (drhodz * - drhodxx + drhodx * drhodzx) / (drhodz**2)
    dSxdy = (drhodz * - drhodxy + drhodx * drhodzy) / (drhodz**2)
    dSxdz = (drhodz * - drhodxz + drhodx * drhodzz) / (drhodz**2)
    
    dSydx = (drhodz * - drhodyx + drhody * drhodzx) / (drhodz**2)
    dSydy = (drhodz * - drhodyy + drhody * drhodzy) / (drhodz**2)
    dSydz = (drhodz * - drhodyz + drhody * drhodzz) / (drhodz**2)
    
    if Sabs < fieldset.Sc - 3. * fieldset.Sd:
        taper1 = 1.
    elif Sabs > fieldset.Sc + 3. * fieldset.Sd:
        taper1 = 0.
    else:
        taper1 = 0.5 * (1 + math.tanh((fieldset.Sc - Sabs)/fieldset.Sd))
    if fieldset.boundaryMask[particle] < 0.95:
        taper2 = 0.
    else:
        taper2 = 1.
    
    k11 = taper1 * taper2 * hsq/2 * (1 + delta) * p
    k21 = taper1 * taper2 * hsq/2 * (1 + delta) * r
    k22 = taper1 * taper2 * hsq/2 * (1 + delta) * q
    k31 = taper1 * taper2 * hsq/2 * (1 + delta) * (p * Sx + r * Sy)
    k32 = taper1 * taper2 * hsq/2 * (1 + delta) * (r * Sx + q * Sy)
    k33 = taper1 * taper2 * (hsq/2 * (1 + delta) * (p * Sx**2 + q * Sy**2 + 2 * r * Sx * Sy) + eps * kappa)
    
    dk11dx = taper1 * taper2 * hsq/2 * ((1 + delta) * dpdx + ddeltadx * p)
    dk11dy = taper1 * taper2 * hsq/2 * ((1 + delta) * dpdy + ddeltady * p)
    dk11dz = taper1 * taper2 * hsq/2 * ((1 + delta) * dpdz + ddeltadz * p)
    
    dk21dx = taper1 * taper2 * hsq/2 * ((1 + delta) * drdx + ddeltadx * r)
    dk21dy = taper1 * taper2 * hsq/2 * ((1 + delta) * drdy + ddeltady * r)
    dk21dz = taper1 * taper2 * hsq/2 * ((1 + delta) * drdz + ddeltadz * r)
    
    dk22dx = taper1 * taper2 * hsq/2 * ((1 + delta) * dqdx + ddeltadx * q)
    dk22dy = taper1 * taper2 * hsq/2 * ((1 + delta) * dqdy + ddeltady * q)
    dk22dz = taper1 * taper2 * hsq/2 * ((1 + delta) * dqdz + ddeltadz * q)

    dk31dx = taper1 * taper2 * hsq/2 * ((ddeltadx) * (p * Sx + r * Sy) + (1 + delta) * (p * dSxdx + dpdx * Sx + r * dSydx + drdx * Sy))
    dk31dy = taper1 * taper2 * hsq/2 * ((ddeltady) * (p * Sx + r * Sy) + (1 + delta) * (p * dSxdy + dpdy * Sx + r * dSydy + drdy * Sy))
    dk31dz = taper1 * taper2 * hsq/2 * ((ddeltadz) * (p * Sx + r * Sy) + (1 + delta) * (p * dSxdz + dpdz * Sx + r * dSydz + drdz * Sy))
    
    dk32dx = taper1 * taper2 * hsq/2 * ((ddeltadx) * (r * Sx + q * Sy) + (1 + delta) * (r * dSxdx + drdx * Sx + q * dSydx + dqdx * Sy))
    dk32dy = taper1 * taper2 * hsq/2 * ((ddeltady) * (r * Sx + q * Sy) + (1 + delta) * (r * dSxdy + drdy * Sx + q * dSydy + dqdy * Sy))
    dk32dz = taper1 * taper2 * hsq/2 * ((ddeltadz) * (r * Sx + q * Sy) + (1 + delta) * (r * dSxdz + drdz * Sx + q * dSydz + dqdz * Sy))
    
    dk33dx = taper1 * taper2 * hsq/2 * ((ddeltadx) * (p * Sx**2 + q * Sy**2 + 2 * r * Sx * Sy) + (1 + delta) * (p * 2 * Sx * dSxdx + dpdx * Sx**2 + q * 2 * Sy * dSydx + dqdx * Sy**2 + 2 * r * Sx * dSydx + 2 * r * dSxdx * Sy + 2 * drdx * Sx * Sy))
    dk33dy = taper1 * taper2 * hsq/2 * ((ddeltady) * (p * Sx**2 + q * Sy**2 + 2 * r * Sx * Sy) + (1 + delta) * (p * 2 * Sx * dSxdy + dpdy * Sx**2 + q * 2 * Sy * dSydy + dqdy * Sy**2 + 2 * r * Sx * dSydy + 2 * r * dSxdy * Sy + 2 * drdy * Sx * Sy))
    dk33dz = taper1 * taper2 * hsq/2 * ((ddeltadz) * (p * Sx**2 + q * Sy**2 + 2 * r * Sx * Sy) + (1 + delta) * (p * 2 * Sx * dSxdz + dpdz * Sx**2 + q * 2 * Sy * dSydz + dqdz * Sy**2 + 2 * r * Sx * dSydz + 2 * r * dSxdz * Sy + 2 * drdz * Sx * Sy))
    