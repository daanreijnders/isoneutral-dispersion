import math
import parcels.rng as ParcelsRandom
    
def density_elements_sampled(particle, fieldset, time):
    drhodx = fieldset.dRhodX[particle]
    drhody = fieldset.dRhodY[particle]
    drhodz = fieldset.dRhodZ[particle]
    
    drhodxx = fieldset.dRhodXX[particle]
    drhodxy = drhodyx = fieldset.dRhodXY[particle]
    drhodxz = drhodzx = fieldset.dRhodXZ[particle]
    drhodyy = fieldset.dRhodYY[particle]
    drhodyz = drhodzy = fieldset.dRhodYZ[particle]
    drhodzz = drhodzz = fieldset.dRhodZZ[particle]
    