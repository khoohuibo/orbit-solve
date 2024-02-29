import orekit
import numpy as np
import csv
import math
import string
import matplotlib.pyplot as plt
from datetime import datetime
import pytz
from scipy.signal import find_peaks

vm = orekit.initVM()

from orekit.pyhelpers import setup_orekit_curdir, absolutedate_to_datetime,datetime_to_absolutedate
setup_orekit_curdir()

from org.orekit.orbits import KeplerianOrbit, OrbitType
from org.orekit.frames import FramesFactory, TopocentricFrame
from org.orekit.bodies import OneAxisEllipsoid, GeodeticPoint
from org.orekit.time import TimeScalesFactory, AbsoluteDate, DateComponents, TimeComponents
from org.orekit.utils import IERSConventions, Constants

from org.orekit.propagation.analytical.tle import TLE, TLEPropagator
from org.orekit.propagation.analytical.tle.generation import FixedPointTleGenerationAlgorithm

from lib.helper_functions import get_abs_vel, mean_motion, distance_between_two, get_solar_declination, get_ecliptic_longitude_of_sun
from lib.propagator import set_up_prop, get_ecef

def checksum(line):
    L = line.strip()
    cksum = 0
    for i in range(68):
        c = L[i]
        if c == ' ' or c == '.' or c == '+' or c in string.ascii_letters:
            continue
        elif c == '-':
            cksum = cksum + 1
        else:
            cksum = cksum + int(c)

    cksum %= 10
    
    return cksum

with open('tle/input_2.tle', 'r') as f:
    output = f.read()

output = output.splitlines()
s = list(output[1])
s[9:15] = ['2', '3', '1', '0', '9', 'G']
s[-1] = str(checksum(''.join(s)))
s = ''.join(s)
#print(s)
t = output[2]

mytle = TLE(s,t)

print (mytle)
print ('Epoch :',mytle.getDate())

ITRF = FramesFactory.getITRF(IERSConventions.IERS_2010, True)
earth = OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS, 
                         Constants.WGS84_EARTH_FLATTENING, 
                         ITRF)


longitude = math.radians(103.8198)
latitude  = math.radians(1.3521)
altitude  = 125.0
station = GeodeticPoint(latitude, longitude, altitude)
station_frame = TopocentricFrame(earth, station, "CRISP")

inertialFrame = FramesFactory.getEME2000()

propagator = TLEPropagator.selectExtrapolator(mytle)

extrapDate = mytle.getDate()
finalDate = extrapDate.shiftedBy(60.0*60*24) #seconds

state_list = []
tle_state_list = []

pv_list = []
pos = []
vel = []
vel_abs = []
acc = []

date = []
dist_list = []

el_store_1=[]
az_store_1 = []
t_store_1 = []

while (extrapDate.compareTo(finalDate) <= 0.0):  
    pv_state = propagator.propagate(extrapDate)
    state_list.append(pv_state)
    pv = pv_state.getPVCoordinates()
    pv_list.append(pv)
    #dist = distance_between_two(pv_state.getPosition(), tle_state.getPosition())
    #dist_list.append(dist)
    pos_tmp = pv.getPosition()

    vel_tmp = pv.getVelocity()

    acc_tmp = pv.getAcceleration()

    pos.append([pos_tmp.getX(),pos_tmp.getY(),pos_tmp.getZ()])
    vel.append([vel_tmp.getX(), vel_tmp.getY(), vel_tmp.getZ()])
    acc.append([acc_tmp.getX(), acc_tmp.getY(), acc_tmp.getZ()])
    vel_abs.append(get_abs_vel(vel_tmp))

    date.append(absolutedate_to_datetime(extrapDate))
    
    el_tmp = station_frame.getElevation(pv.getPosition(),
                    inertialFrame,
                    extrapDate)*180.0/math.pi
    print(el_tmp, extrapDate)
    el_store_1.append(el_tmp)

    az_tmp = station_frame.getAzimuth(pv.getPosition(),
                    inertialFrame,
                    extrapDate)#*180.0/math.pi
    az_store_1.append(az_tmp)
    #print extrapDate, pos_tmp, vel_tmp
    t_store_1.append(extrapDate)
    extrapDate = extrapDate.shiftedBy(1.0)

#print(date)
derivedOrbit_alt, seedOrbit, newEpoch = get_ecef(date, absolutedate_to_datetime(mytle.getDate()), vel_abs, pos, vel, pv_list, inertialFrame, ITRF, inclination=mytle.getI(), t_delta_0=60)

#print(derivedOrbit_alt)
#print(seedOrbit)
#print(OrbitType.KEPLERIAN.convertType(state_list[-1].getOrbit()))
#print(newEpoch)

extrapDate = mytle.getDate()
finalDate = extrapDate.shiftedBy(60.0*60*24)
derivedpropagator = set_up_prop(0,0,0,0,0,0,newEpoch,inertialFrame,ITRF,initialOrbit=derivedOrbit_alt)

derived_state_list = []
dist_list = []

while (extrapDate.compareTo(finalDate) <= 0.0):
    derived_pv_state = derivedpropagator.propagate(extrapDate)
    derived_state_list.append(derived_pv_state)

    dist = distance_between_two(pv_state.getPosition(), derived_pv_state.getPosition())/1000
    dist_list.append(dist)





newTLE = TLE.stateToTLE(state_list[-1], mytle, FixedPointTleGenerationAlgorithm())
print(newTLE)

with open("tle/output_guess.tle", "w") as f:
    f.write("NEW_GUESS\r\n")
    f.write(str(newTLE))
    f.close()


peaks, _ =find_peaks(el_store_1, height=0)

if len(peaks) == 0:
    raise Exception

width=300

"""
fig = plt.figure()

ax = plt.subplot(1,1,1, polar=True)

ax.set_theta_zero_location('N')
ax.set_theta_direction(-1)
ax.set_rlim(bottom=90, top=0)

ax.plot(az, el)
ax.grid(True)
"""

#print(el_store_1)
print(len(peaks))

fig = plt.figure()

for index in range(len(peaks)-4, len(peaks)):
    #print(az_store_1[peaks[index]-width:peaks[index]+width])
    #print(el_store_1[peaks[index]-width:peaks[index]+width])
    #print(t_store_1[peaks[index]-width:peaks[index]+width])
    if index < 0:
        continue
    ax = plt.subplot(3, 2, index-(len(peaks)-4)+3, polar=True)
    
    ax.set_theta_zero_location('N')

    ax.set_theta_direction(-1)

    ax.set_rlim(bottom=90, top=0)

    #ax.plot(az_store_2[peaks_2[index]-width:peaks_2[index]+width], el_store_2[peaks_2[index]-width:peaks_2[index]+width], color='b', zorder=1, label='OLD TLE')
    
    ax.plot(az_store_1[peaks[index]-width:peaks[index]+width], el_store_1[peaks[index]-width:peaks[index]+width], linestyle='--', color='r', linewidth=2, zorder=2, label='NEW TLE')
    
    ax.legend()
    
    ax.grid(True)
    
    ax.set_title("%s \r\n %s" % (t_store_1[peaks[index]], t_store_1[peaks[index]]))


plt.show()