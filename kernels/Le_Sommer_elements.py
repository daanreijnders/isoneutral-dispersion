import math
import parcels.rng as ParcelsRandom

def Le_Sommer_elements(particle, fieldset, time):    
    delta = fieldset.Delta[particle]
    p = fieldset.P[particle]
    q = fieldset.Q[particle]
    r = fieldset.R[particle]
    
    ddeltadx = fieldset.ddeltadX[particle]
    ddeltady = fieldset.ddeltadY[particle]
    ddeltadz = fieldset.ddeltadZ[particle]
    
    dpdx = fieldset.dpdX[particle]
    dpdy = fieldset.dpdY[particle]
    dpdz = fieldset.dpdZ[particle]
    
    dqdx = fieldset.dqdX[particle]
    dqdy = fieldset.dqdY[particle]
    dqdz = fieldset.dqdZ[particle]
    
    drdx = fieldset.drdX[particle]
    drdy = fieldset.drdY[particle]
    drdz = fieldset.drdZ[particle]