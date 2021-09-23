import math
import parcels.rng as ParcelsRandom

def K_elemns_sampled_S_taper_Kpp(particle, fieldset, time):
# fmt: off
    sx = fieldset.Sx[time, particle.depth, particle.lat, particle.lon]
    sy = fieldset.Sy[time, particle.depth, particle.lat, particle.lon]
    sabs2 = sx ** 2 + sy ** 2
    
    sx_xp1 = fieldset.Sx[time, particle.depth, particle.lat, particle.lon + dx]
    sx_xm1 = fieldset.Sx[time, particle.depth, particle.lat, particle.lon - dx]
    sx_yp1 = fieldset.Sx[time, particle.depth, particle.lat + dy, particle.lon]
    sx_ym1 = fieldset.Sx[time, particle.depth, particle.lat - dy, particle.lon]
    sx_zp1 = fieldset.Sx[time, particle.depth + dz, particle.lat, particle.lon]
    sx_zm1 = fieldset.Sx[time, particle.depth - dz, particle.lat, particle.lon]
    sy_xp1 = fieldset.Sy[time, particle.depth, particle.lat, particle.lon + dx]
    sy_xm1 = fieldset.Sy[time, particle.depth, particle.lat, particle.lon - dx]
    sy_yp1 = fieldset.Sy[time, particle.depth, particle.lat + dy, particle.lon]
    sy_ym1 = fieldset.Sy[time, particle.depth, particle.lat - dy, particle.lon]
    sy_zp1 = fieldset.Sy[time, particle.depth + dz, particle.lat, particle.lon]
    sy_zm1 = fieldset.Sy[time, particle.depth - dz, particle.lat, particle.lon]
    
    sabs2_xp1 = sx_xp1**2 + sy_xp1**2
    sabs2_xm1 = sx_xm1**2 + sy_xm1**2
    sabs2_yp1 = sx_yp1**2 + sy_yp1**2
    sabs2_ym1 = sx_ym1**2 + sy_ym1**2
    sabs2_zp1 = sx_zp1**2 + sy_zp1**2
    sabs2_zm1 = sx_zm1**2 + sy_zm1**2
    
    if sabs2 > 0:
        taper = 0.5 * (1. + math.tanh((fieldset.Sc - math.sqrt(sabs2)) / fieldset.Sd))
    else:
        taper = 1
    
    if sabs2_xp1 > 0:
        taper_xp1 = 0.5 * (1. + math.tanh((fieldset.Sc - math.sqrt(sabs2_xp1)) / fieldset.Sd))
    else:
        taper_xp1 = 1
    
    if sabs2_xm1 > 0:
        taper_xm1 = 0.5 * (1. + math.tanh((fieldset.Sc - math.sqrt(sabs2_xm1)) / fieldset.Sd))
    else:
        taper_xm1 = 1
        
    if sabs2_yp1 > 0:
        taper_yp1 = 0.5 * (1. + math.tanh((fieldset.Sc - math.sqrt(sabs2_yp1)) / fieldset.Sd))
    else:
        taper_yp1 = 1
        
    if sabs2_ym1 > 0:
        taper_ym1 = 0.5 * (1. + math.tanh((fieldset.Sc - math.sqrt(sabs2_ym1)) / fieldset.Sd))
    else:
        taper_ym1 = 1
        
    if sabs2_zp1 > 0:
        taper_zp1 = 0.5 * (1. + math.tanh((fieldset.Sc - math.sqrt(sabs2_zp1)) / fieldset.Sd))
    else:
        taper_zp1 = 1
        
    if sabs2_zm1 > 0:
        taper_zm1 = 0.5 * (1. + math.tanh((fieldset.Sc - math.sqrt(sabs2_zm1)) / fieldset.Sd))
    else:
        taper_zm1 = 1
    
    k11 = taper * fieldset.Ki
    k21 = 0
    k22 = taper * fieldset.Ki
    k31 = taper * fieldset.Ki * sx
    k32 = taper * fieldset.Ki * sy
    k33 = taper * fieldset.Ki * sabs2
    
    kv = fieldset.Kv[time, particle.depth, particle.lat, particle.lon]
    kv_zp1 = fieldset.Kv[time, particle.depth + dz, particle.lat, particle.lon] 
    kv_zm1 = fieldset.Kv[time, particle.depth - dz, particle.lat, particle.lon] 
    dkvdz = (kv_zp1 - kv_zm1) / (2 * dz)
    kv_xp1 = fieldset.Kv[time, particle.depth, particle.lat, particle.lon + dx] 
    kv_xm1 = fieldset.Kv[time, particle.depth, particle.lat, particle.lon - dx] 
    dkvdx = (kv_xp1 - kv_xm1) / (2 * dx)
    kv_yp1 = fieldset.Kv[time, particle.depth, particle.lat + dy, particle.lon] 
    kv_ym1 = fieldset.Kv[time, particle.depth, particle.lat - dy, particle.lon] 
    dkvdy = (kv_yp1 - kv_ym1) / (2 * dy)
    
    dk11dx = 0
    dk12dy = 0
    dk13dz = fieldset.Ki * ((taper_zp1 - taper_zm1) / (2 * dz) * sx + taper * (sx_zp1 - sx_zm1) / (2 * dz))
    dk21dx = 0
    dk22dy = 0
    dk23dz = fieldset.Ki * ((taper_zp1 - taper_zm1) / (2 * dz) * sy + taper * (sy_zp1 - sy_zm1) / (2 * dz))
    dk31dx = fieldset.Ki * ((taper_xp1 - taper_xm1) / (2 * dx) * sx + taper * (sx_xp1 - sx_xm1) / (2 * dx))
    dk31dy = fieldset.Ki * ((taper_yp1 - taper_ym1) / (2 * dy) * sx + taper * (sx_yp1 - sx_ym1) / (2 * dy))
    dk31dz = fieldset.Ki * ((taper_zp1 - taper_zm1) / (2 * dz) * sx + taper * (sx_zp1 - sx_zm1) / (2 * dz))
    dk32dx = fieldset.Ki * ((taper_xp1 - taper_xm1) / (2 * dx) * sy + taper * (sy_xp1 - sy_xm1) / (2 * dx))
    dk32dy = fieldset.Ki * ((taper_yp1 - taper_ym1) / (2 * dy) * sy + taper * (sy_yp1 - sy_ym1) / (2 * dy))
    dk32dz = fieldset.Ki * ((taper_zp1 - taper_zm1) / (2 * dz) * sy + taper * (sy_zp1 - sy_zm1) / (2 * dz))
    dk33dz = fieldset.Ki * ((taper_zp1 - taper_zm1) / (2 * dz) * sabs2 + taper * (sabs2_zp1 - sabs2_zm1) / (2 * dz))

    # fmt: on