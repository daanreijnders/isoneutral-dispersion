import math
import parcels.rng as ParcelsRandom


def K_elemns_analytical(particle, fieldset, time):
    # fmt: off
    g = fieldset.gval
    Nsquared = fieldset.Nsquaredval
    alpha = fieldset.alphaval
    kappa = fieldset.kappaval
    rho0 = fieldset.rho0val

    x = particle.lon
    y = particle.lat
    
    drhodx = rho0 * alpha * kappa * math.cos(kappa * x) * (y / 6000 + 0.25)
    drhody = -rho0 * Nsquared / g + alpha * math.sin(kappa * x) / 6000
    
    drhodxx = - rho0 * alpha * kappa**2 * math.sin(kappa * x) * (y / 6000 + 0.25)
    drhodxy = rho0 * alpha * kappa * math.cos(kappa * x) / 6000
    
    denom = (drhodx**2 + drhody**2)
    k_const = fieldset.Ki / denom
    
    # Lower triangular K elements only, since K is symmetric
    k11 = k_const * (drhody**2)
    k21 = k_const * (-drhodx * drhody)
    k22 = k_const * (drhodx**2)

    dk11dx = -fieldset.Ki * drhody**2 * 2 * drhodx * drhodxx / denom**2
    dk12dy = -fieldset.Ki * (denom * drhody * drhodxy - 2 * drhodx**2 * drhody * drhodxy) / denom**2
    dk21dx = fieldset.Ki * (-denom * drhody * drhodxx + 2 * drhody * drhodx**2 * drhodxx) / denom**2
    dk22dy = fieldset.Ki * (drhody**2 * 2 * drhodx * drhodxy) / denom**2
    # fmt: on