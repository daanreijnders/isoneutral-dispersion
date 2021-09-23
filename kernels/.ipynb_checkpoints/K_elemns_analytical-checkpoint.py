import math
import parcels.rng as ParcelsRandom

def K_elemns_analytical(particle, fieldset, time):
    # fmt: off
    g = fieldset.gval
    Nsquared = fieldset.Nsquaredval
    alphax = fieldset.alphaxval
    alphay = fieldset.alphayval
    kappax = fieldset.kappaxval
    kappay = fieldset.kappayval
    rho0 = fieldset.rho0val

    x = particle.lon
    y = particle.lat
    
    drhodx = rho0 * alphax * kappax* math.cos(kappax * x)
    drhody = rho0 * alphay * kappay * math.cos(kappay * y)
    drhodz = -rho0 * Nsquared / g
    
    drhodxx = - rho0 * alphax * kappax**2 * math.sin(kappax * x)
    drhodyy = - rho0 * alphay * kappay**2 * math.sin(kappay * y)
    
    denom = (drhodx**2 + drhody**2 + drhodz**2)
    k_const = fieldset.Ki / denom
    
    # Lower triangular K elements only, since K is symmetric
    k11 = k_const * (drhody**2 + drhodz**2)
    k21 = k_const * (-drhodx * drhody)
    k22 = k_const * (drhodx**2 + drhodz**2)
    k31 = k_const * (-drhodz*drhodx)
    k32 = k_const * (-drhodz*drhody)
    k33 = k_const * (drhodx**2 + drhody**2)

    dk11dx = -fieldset.Ki * (drhody**2 + drhodz**2) * 2 * drhodx * drhodxx / denom**2
    dk12dy = fieldset.Ki * (-denom * drhodx * drhodyy + 2 * drhodx * drhody**2 * drhodyy) / denom**2
    dk13dz = 0
    dk21dx = fieldset.Ki * (-denom * drhody * drhodxx + 2 * drhody * drhodx**2 * drhodxx) / denom**2
    dk22dy = -fieldset.Ki * (drhodx**2 + drhodz**2) * 2 * drhody * drhodyy / denom**2
    dk23dz = 0
    dk31dx = fieldset.Ki * (-denom * drhodz * drhodxx + 2 * drhodz * drhodx**2 * drhodxx) / denom**2
    dk31dy = fieldset.Ki * drhodz * drhodx * 2 * drhody * drhodyy / denom**2
    dk31dz = 0
    dk32dx = fieldset.Ki * drhodz * drhody * 2 * drhodx * drhodxx / denom**2
    dk32dy = fieldset.Ki * (-denom * drhodz * drhodyy + 2 * drhodz * drhody**2 * drhodyy) / denom**2
    dk32dz = 0
    dk33dz = 0

    # fmt: on