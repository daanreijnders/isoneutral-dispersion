import math
import parcels.rng as ParcelsRandom


def K_elemns_analytical_smallslope(particle, fieldset, time):
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
    
    k11 = fieldset.Ki 
    k21 = 0
    k22 = fieldset.Ki
    k31 = fieldset.Ki * (-drhodx / drhodz)
    k32 = fieldset.Ki * (-drhody / drhodz)
    k33 = fieldset.Ki * ((drhodx / drhodz)**2 + (drhody / drhodz)**2)

    dk11dx = 0
    dk12dy = 0
    dk13dz = 0
    dk21dx = 0
    dk22dy = 0
    dk23dz = 0
    dk31dx = fieldset.Ki * (-drhodxx / drhodz)
    dk31dy = 0
    dk31dz = 0
    dk32dx = 0
    dk32dy = fieldset.Ki * (-drhodyy / drhodz)
    dk32dz = 0
    dk33dz = 0
    
    # fmt: on