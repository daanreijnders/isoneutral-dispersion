import math
import parcels.rng as ParcelsRandom

def elements_analytical(particle, fieldset, time):
    # fmt: off
    g = fieldset.gval
    Nsquared = fieldset.Nsquaredval
    alphax = fieldset.alphaxval
    alphay = fieldset.alphayval
    kappax = fieldset.kappaxval
    kappay = fieldset.kappayval
    rho0 = fieldset.rho0val

    drhodx = rho0 * alphax * kappax * math.cos(kappax * particle.lon)
    drhody = rho0 * alphay * kappay * math.cos(kappay * particle.lat)
    drhodz = -rho0 * Nsquared / g
    
    drhodxx = rho0 * alphax * kappax**2 * -math.sin(kappax * particle.lon)
    drhodxy = drhodyx = 0
    drhodxz = drhodzx = 0
    drhodyy = rho0 * alphay * kappay**2 * -math.sin(kappay * particle.lat)
    drhodyz = drhodzy = 0
    drhodzz = drhodzz = 0
    
    # Le Sommer
    delta = 1
    p = 1
    q = 1
    r = 0.1
    
    ddeltadx = 0
    ddeltady = 0
    ddeltadz = 0
    
    dpdx = 0
    dpdy = 0
    dpdz = 0
    
    dqdx = 0
    dqdy = 0
    dqdz = 0
    
    drdx = 0
    drdy = 0
    drdz = 0
