def Markov1_3D_BC_taper_init(particle, fieldset, time):
    # Initialize random velocity magnitudes in the isopycnal plane
    u_iso_prime_abs = ParcelsRandom.normalvariate(0, math.sqrt(fieldset.nusquared))
    v_iso_prime_abs = ParcelsRandom.normalvariate(0, math.sqrt(fieldset.nusquared))
    
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
        taper1 = 0.5 * (1 + math.tanh((fieldset.Sc - Sabs)/fieldset.Sd))
    # Near boundaries, drop to 0
    if fieldset.boundaryMask[particle] < 0.95:
        taper2 = 0.
    else:
        taper2 = 1.
    
    # Compute the gradient vector, which is perpendicular to the local isoneutral
    grad = [drhodx, drhody, drhodz]
    
    # Compute u in the isopycnal plane (arbitrary direction)
    # by ensuring dot(u_iso, grad) = 0
    u_iso = [-(drhody+drhodz)/drhodx, 1, 1]
    u_iso_norm = math.sqrt(u_iso[0]**2 + u_iso[1]**2 + u_iso[2]**2)
    u_iso_normed = [u_iso[0]/u_iso_norm, u_iso[1]/u_iso_norm, u_iso[2]/u_iso_norm]
    
    # Fix v_iso by computing cross(u_iso, grad)
    v_iso = [grad[1]*u_iso[2] - grad[2]*u_iso[1],
             grad[2]*u_iso[0] - grad[0]*u_iso[2],
             grad[0]*u_iso[1] - grad[1]*u_iso[0]]
    v_iso_norm = math.sqrt(v_iso[0]**2 + v_iso[1]**2 + v_iso[2]**2)
    v_iso_normed = [v_iso[0]/v_iso_norm, v_iso[1]/v_iso_norm, v_iso[2]/v_iso_norm]
    
    # Compute the initial isopycnal velocity vector, which is a linear combination
    # of u_iso and v_iso
    vel_init = [u_iso_prime_abs * u_iso_normed[0] + v_iso_prime_abs * v_iso_normed[0],
                u_iso_prime_abs * u_iso_normed[1] + v_iso_prime_abs * v_iso_normed[1],
                u_iso_prime_abs * u_iso_normed[2] + v_iso_prime_abs * v_iso_normed[2]]
    
    
    particle.u_prime = taper1 * taper2 * vel_init[0]
    particle.v_prime = taper1 * taper2 * vel_init[1]
    particle.w_prime = taper1 * taper2 * vel_init[2]