import math
import parcels.rng as ParcelsRandom

def K_Redi_smallslope(particle, fieldset, time):
    kappa = fieldset.kappaval
    eps = fieldset.epsilon
    
    Sx = - drhodx / drhodz
    Sy = - drhody / drhodz
    Sabs = math.sqrt(Sx**2 + Sy**2)
    Sabs2 = Sx**2 + Sy**2
    
    dSxdx = (drhodz * - drhodxx + drhodx * drhodzx) / (drhodz**2)
    dSxdy = (drhodz * - drhodxy + drhodx * drhodzy) / (drhodz**2)
    dSxdz = (drhodz * - drhodxz + drhodx * drhodzz) / (drhodz**2)
    
    dSydx = (drhodz * - drhodyx + drhody * drhodzx) / (drhodz**2)
    dSydy = (drhodz * - drhodyy + drhody * drhodzy) / (drhodz**2)
    dSydz = (drhodz * - drhodyz + drhody * drhodzz) / (drhodz**2)
    
    dSabsdx = (Sx*dSxdx + Sy*dSydx)/Sabs
    dSabsdy = (Sx*dSxdy + Sy*dSydy)/Sabs
    dSabsdz = (Sx*dSxdz + Sy*dSydz)/Sabs

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
    
    k11 = taper1 * taper2 * kappa
    k21 = taper1 * taper2 * kappa * - Sx * Sy
    k22 = taper1 * taper2 * kappa
    k31 = taper1 * taper2 * kappa * Sx
    k32 = taper1 * taper2 * kappa * Sy
    k33 = taper1 * taper2 * kappa * (Sabs**2 + eps)
    
    dk11dx = 0
    dk11dy = 0
    dk11dz = 0
    
    dk21dx = taper1 * taper2 * kappa * - (Sx * dSydx + dSxdx * Sy)
    dk21dy = taper1 * taper2 * kappa * - (Sx * dSydy + dSxdy * Sy)
    dk21dz = taper1 * taper2 * kappa * - (Sx * dSydz + dSxdz * Sy)
    
    dk22dx = 0
    dk22dy = 0
    dk22dz = 0

    dk31dx = taper1 * taper2 * kappa * dSxdx
    dk31dy = taper1 * taper2 * kappa * dSxdy
    dk31dz = taper1 * taper2 * kappa * dSxdz
    
    dk32dx = taper1 * taper2 * kappa * dSydx
    dk32dy = taper1 * taper2 * kappa * dSydy
    dk32dz = taper1 * taper2 * kappa * dSydz
    
    dk33dx = taper1 * taper2 * kappa * 2 * Sabs * dSabsdx
    dk33dy = taper1 * taper2 * kappa * 2 * Sabs * dSabsdy
    dk33dz = taper1 * taper2 * kappa * 2 * Sabs * dSabsdz