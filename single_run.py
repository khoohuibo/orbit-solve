import orekit
import numpy as np
import csv
import sys, getopt



from orekit.pyhelpers import setup_orekit_curdir, absolutedate_to_datetime,datetime_to_absolutedate
from org.orekit.orbits import OrbitType, KeplerianOrbit
from org.orekit.time import AbsoluteDate, TimeScalesFactory
from org.orekit.utils import  IERSConventions, PVCoordinatesProvider
from org.orekit.frames import FramesFactory
from org.orekit.propagation.analytical.tle import TLE
from org.orekit.propagation.analytical.tle.generation import FixedPointTleGenerationAlgorithm
from org.orekit.bodies import CelestialBodyFactory
from org.hipparchus.geometry.euclidean.threed import Vector3D

from lib.helper_functions import get_abs_vel, mean_motion, distance_between_two, get_beta_angle_alternate, get_beta_angle
from lib.propagator import set_up_prop, get_ecef

from orekit import JArray_double
import math
import matplotlib.pyplot as plt

vm = orekit.initVM()
setup_orekit_curdir()
utc = TimeScalesFactory.getUTC()

ra = 421 *  1000         #  Apogee
rp = 414 * 1000         #  Perigee
i = math.radians(51.6)      # inclination
omega = math.radians(0.0)   # perigee argument
raan = math.radians(0.0)  # right ascension of ascending node
lv = math.radians(0.0)    # True anomaly

epochDate = AbsoluteDate(2025, 7, 1, 0, 0, 00.000, utc)
initialDate = epochDate

inertialFrame = FramesFactory.getEME2000()

ITRF = FramesFactory.getITRF(IERSConventions.IERS_2010, True)

propagator = set_up_prop(rp, ra, i, omega, raan, lv, epochDate, inertialFrame, ITRF)

print(propagator.getInitialState().getOrbit())

extrapDate = AbsoluteDate(2025, 7, 1, 0, 0, 00.000, utc)
finalDate = extrapDate.shiftedBy(60.0*5)

state_list = []
tle_state_list = []

pv_list = []
pos = []

vel = []

vel_abs = []
acc = []

date = []
dist_list = []


while (extrapDate.compareTo(finalDate) <= 0.0):  
    
    #pv = propagator.getPVCoordinates(extrapDate, inertialFrame)
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
    print(extrapDate, end="\r")
    extrapDate = extrapDate.shiftedBy(1.0)


tle_first_guess = TLE(99997, 
                        "0",
                        2025,
                        1,
                        "F",
                        0,
                        999,
                        finalDate,
                        mean_motion(16, 600),
                        0.0, 
                        0.0,
                        propagator.getInitialState().getE(),
                        i,
                        omega,
                        raan,
                        lv,
                        100,
                        0.0)
print("HERE HERE HERE!")
print(OrbitType.KEPLERIAN.convertType(state_list[1].getOrbit()))
myTLE = TLE.stateToTLE(state_list[138], tle_first_guess, FixedPointTleGenerationAlgorithm())
print(myTLE)
#print(OrbitType.KEPLERIAN.convertType(state_list[3].getOrbit()))
#print(date)
derivedOrbit_alt, seedOrbit, newEpoch, t_delta = get_ecef(date, absolutedate_to_datetime(epochDate), vel_abs, pos, vel, pv_list, inertialFrame, ITRF, inclination=i)

print(derivedOrbit_alt)
#print(seedOrbit)

print(t_delta)

sma = derivedOrbit_alt.getA()
ecc = derivedOrbit_alt.getE()
inc = derivedOrbit_alt.getI()
aop = derivedOrbit_alt.getPerigeeArgument()
raan_alt = derivedOrbit_alt.getRightAscensionOfAscendingNode()
lv_new = seedOrbit.getMeanAnomaly() + t_delta * mean_motion(16, 600)
#derivedOrbit_alt.getTrueAnomaly()#
print(lv_new)
print(t_delta * mean_motion(16, 600))
propagator = set_up_prop(rp, ra, i, omega, raan, lv, epochDate, inertialFrame, ITRF, max_drag=False)

raan_list = []
max_dist_list = []


raan_new = raan_alt 
derivedpropagator = set_up_prop(rp, ra, inc, aop, raan_new, lv_new, newEpoch, inertialFrame, ITRF, a=sma, e=ecc, max_drag=True)
#derivedpropagator = set_up_prop(rp, ra, i, omega, raan, lv, epochDate, inertialFrame, ITRF, max_drag=True)
print(propagator.getInitialState().getOrbit())
print(derivedpropagator.getInitialState().getOrbit())

date = []
dist_list = []

state_list = []
derived_state_list = []

beta_list_sun = []
beta_list = []
beta_list_alt = []
derived_beta_list_sun=[]

raan_list= []

extrapDate = AbsoluteDate(2025, 7, 1, 0, 0, 00.000, utc)
finalDate = extrapDate.shiftedBy(60.0*60*24*365)
first = False
sun = CelestialBodyFactory.getSun();
sun = PVCoordinatesProvider.cast_(sun)  # But we want the PVCoord interface

while (extrapDate.compareTo(finalDate) <= 0.0):  

    #pv = propagator.getPVCoordinates(extrapDate, inertialFrame)
    pv_state = propagator.propagate(extrapDate)
    state_list.append(pv_state)
    pv = pv_state.getPVCoordinates(inertialFrame)

    pos_sun = sun.getPVCoordinates(extrapDate, inertialFrame).getPosition()
    print(pos_sun)
    #pos_sun_new = sun.getPosition(extrapDate, inertialFrame)
    beta_angle_new = 0.5 * math.pi - Vector3D.angle(pos_sun, pv.getMomentum())
    print("vector_angle: %f" % math.degrees(Vector3D.angle(pos_sun, pv.getMomentum())))

    beta_list_sun.append(math.degrees(beta_angle_new))

    int_orbit = KeplerianOrbit(OrbitType.KEPLERIAN.convertType(pv_state.getOrbit()))
    print(int_orbit)
    raan_list.append(math.degrees(int_orbit.getRightAscensionOfAscendingNode()))
    #print(math.degrees(int_orbit.getRightAscensionOfAscendingNode()))
    beta_angle = get_beta_angle(absolutedate_to_datetime(extrapDate), int_orbit.getRightAscensionOfAscendingNode(), int_orbit.getI())
    #beta_angle_alt = get_beta_angle_alternate(absolutedate_to_datetime(extrapDate), int_orbit.getRightAscensionOfAscendingNode(), int_orbit.getI())
    #beta_list.append(math.degrees(beta_angle))
    #beta_list_alt.append(math.degrees(beta_angle_alt))

    derived_pv_state = derivedpropagator.propagate(extrapDate)
    derived_state_list.append(derived_pv_state)

    derived_pv = derived_pv_state.getPVCoordinates(inertialFrame)

    pos_sun = sun.getPVCoordinates(extrapDate, inertialFrame).getPosition()
    #pos_sun_new = sun.getPosition(extrapDate, inertialFrame)
    beta_angle_new = 0.5 * math.pi - Vector3D.angle(pos_sun, derived_pv.getMomentum())

    derived_beta_list_sun.append(math.degrees(beta_angle_new))

    dist = distance_between_two(pv_state.getPosition(), derived_pv_state.getPosition())/1000
    dist_list.append(dist)

    date.append(absolutedate_to_datetime(extrapDate))
    print(extrapDate, end="\r\n")
    extrapDate = extrapDate.shiftedBy(3600.0)

header = ['date', 'distance', 'beta_angle_sun', 'derived_beta_angle_sun']
with open('output_sso/beta_three_sso_july_2025_true_delayed.csv', 'w') as f:
    writer = csv.writer(f)

    writer.writerow(header)
    for j in range(len(date)):
        writer.writerow([date[j], dist_list[j], beta_list_sun[j], derived_beta_list_sun[j]])

    f.close()

print(OrbitType.KEPLERIAN.convertType(state_list[0].getOrbit()))

print(OrbitType.KEPLERIAN.convertType(derived_state_list[0].getOrbit()))
    
print(OrbitType.KEPLERIAN.convertType(state_list[-1].getOrbit()))

print(OrbitType.KEPLERIAN.convertType(derived_state_list[-1].getOrbit()))


val, idx = min((val, idx) for (idx, val) in enumerate(beta_list_sun))

print(date[idx])
print(OrbitType.KEPLERIAN.convertType(state_list[idx].getOrbit()))
print(beta_list_sun[idx])

val, idx = max((val, idx) for (idx, val) in enumerate(beta_list_sun))

print(date[idx])
print(OrbitType.KEPLERIAN.convertType(state_list[idx].getOrbit()))
print(beta_list_sun[idx])

for i in range(len(date)):
    if dist_list[i] <= 1000:
        print("date: %s, dist: %f" % (date[i], dist_list[i]))

#plt.plot(date, dist_list)
plt.plot(date, beta_list_sun)
#plt.plot(date, beta_list)
#plt.plot(date, beta_list_alt)
plt.plot(date, derived_beta_list_sun)
plt.legend(['1', '2', '3', '4']) 
#plt.plot(date, raan_list)
plt.show()