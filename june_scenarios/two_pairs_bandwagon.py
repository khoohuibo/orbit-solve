import orekit
import numpy as np
import csv
import sys, getopt



from orekit.pyhelpers import setup_orekit_curdir, absolutedate_to_datetime,datetime_to_absolutedate
from org.orekit.orbits import OrbitType, KeplerianOrbit
from org.orekit.time import AbsoluteDate, TimeScalesFactory
from org.orekit.utils import  IERSConventions, PVCoordinatesProvider, Constants
from org.orekit.frames import FramesFactory
from org.orekit.propagation.analytical.tle import TLE
from org.orekit.propagation.analytical.tle.generation import FixedPointTleGenerationAlgorithm
from org.orekit.bodies import CelestialBodyFactory, OneAxisEllipsoid
from org.hipparchus.geometry.euclidean.threed import Vector3D

from lib.helper_functions import get_abs_vel, mean_motion, distance_between_two, get_beta_angle_alternate, get_beta_angle, get_abs_acc
from lib.propagator import set_up_prop, get_ecef, add_dsst_force_models

from orekit import JArray_double
import math
import matplotlib.pyplot as plt

vm = orekit.initVM()
setup_orekit_curdir()
utc = TimeScalesFactory.getUTC()
satellite_mass = 16
start_date = AbsoluteDate(2025, 6, 1, 0, 0, 00.000, utc)
pair_2_start_date = AbsoluteDate(2026, 1, 1, 0, 0, 00.000, utc)

ra = 600.14 * 1000        #  Apogee
rp = 600.139999 * 1000         #  Perigee
i = math.radians(97.7532)      # inclination (SSO)
omega = math.radians(0.0)   # perigee argument
raan = math.radians(167.272)  # right ascension of ascending node (SSO)
lv = math.radians(0.0)    # True anomaly

bandwagon_i = math.radians(45.1) # inclination (Bandwagon)
bandwagon_raan = math.radians(167.272) # RAAN (Bandwagon)

epochDate = start_date
#epochDate = AbsoluteDate(2024, 4, 1, 0, 0, 00.000, utc)
initialDate = epochDate

inertialFrame = FramesFactory.getEME2000()

ITRF = FramesFactory.getITRF(IERSConventions.IERS_2010, True)

propagator, _ = set_up_prop(i, omega, raan, lv, epochDate, inertialFrame, ITRF, rp=rp, ra=ra)

print(propagator.getInitialState().getOrbit())

extrapDate = start_date
#extrapDate = AbsoluteDate(2024, 4, 1, 0, 0, 00.000, utc)
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
                        mean_motion(satellite_mass, 600),
                        0.0, 
                        0.0,
                        propagator.getInitialState().getE(),
                        i,
                        omega,
                        raan,
                        lv,
                        100,
                        0.0)

myTLE = TLE.stateToTLE(state_list[138], tle_first_guess, FixedPointTleGenerationAlgorithm())

derivedOrbit_alt, seedOrbit, newEpoch, t_delta = get_ecef(date, absolutedate_to_datetime(epochDate), vel_abs, pos, vel, pv_list, inertialFrame, ITRF, inclination=i)

sma = derivedOrbit_alt.getA()
ecc = derivedOrbit_alt.getE()
inc = derivedOrbit_alt.getI()
aop = derivedOrbit_alt.getPerigeeArgument()
raan_alt = derivedOrbit_alt.getRightAscensionOfAscendingNode()
lv_new = seedOrbit.getMeanAnomaly() + t_delta * mean_motion(satellite_mass, 600)
propagator, og_driver = set_up_prop(i, omega, raan, lv, epochDate, inertialFrame, ITRF,rp=rp, ra=ra, max_drag=False)

derivedpropagator, dervied_driver = set_up_prop(inc, aop, raan_alt, lv_new, newEpoch, inertialFrame, ITRF, rp=rp, ra=ra, a=sma, e=ecc, max_drag=False)

pair_2_propagator, _ = set_up_prop(bandwagon_i, omega, bandwagon_raan, lv, pair_2_start_date, inertialFrame, ITRF, rp=rp, ra=ra)

extrapDate = pair_2_start_date
#extrapDate = AbsoluteDate(2024, 4, 1, 0, 0, 00.000, utc)
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
    pv_state = pair_2_propagator.propagate(extrapDate)
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
                        2026,
                        1,
                        "F",
                        0,
                        999,
                        finalDate,
                        mean_motion(satellite_mass, 600),
                        0.0, 
                        0.0,
                        propagator.getInitialState().getE(),
                        i,
                        omega,
                        raan,
                        lv,
                        100,
                        0.0)

myTLE = TLE.stateToTLE(state_list[138], tle_first_guess, FixedPointTleGenerationAlgorithm())

derived2Orbit_alt, seed2Orbit, newEpoch, t_delta = get_ecef(date, absolutedate_to_datetime(epochDate), vel_abs, pos, vel, pv_list, inertialFrame, ITRF, inclination=bandwagon_i)

sma_2 = derived2Orbit_alt.getA()
ecc_2 = derived2Orbit_alt.getE()
raan_2 = derived2Orbit_alt.getRightAscensionOfAscendingNode()
lv_2 = seed2Orbit.getMeanAnomaly() + t_delta * mean_motion(satellite_mass, 600)

pair_2_1, _ = set_up_prop(bandwagon_i, omega, bandwagon_raan, lv, pair_2_start_date, inertialFrame, ITRF,rp=rp, ra=ra, max_drag=False)

pair_2_2, _ = set_up_prop(bandwagon_i, omega, raan_2, lv_2, pair_2_start_date, inertialFrame, ITRF,rp=rp, ra=ra, max_drag=False)

print("SEED ORBIT-------------------------------")
print(propagator.getInitialState().getOrbit())

print("DERIVED ORBIT-------------------------------")
print(derivedpropagator.getInitialState().getOrbit())

print("SEED 2 ORBIT-----------------------")
print(pair_2_1.getInitialState().getOrbit())

print("DERIVED 2 ORBIT-----------------------")
print(pair_2_2.getInitialState().getOrbit())

date = []
dist_list = []
A_C_list = []
A_D_list = []
C_D_list = []
pair_2_date = []

state_list = []
derived_state_list = []
seed_2_state_list = []
derived_2_state_list = []

beta_list_sun = []
beta_list = []
beta_list_alt = []
derived_beta_list_sun=[]

raan_list= []
i_list = []
sma_list = []
ecc_list = []
aop_list = []
lv_list = []

pair_2_sma_list = []
B_D_list = []
B_C_list = []

extrapDate = start_date
#extrapDate = AbsoluteDate(2024, 4, 1, 0, 0, 00.000, utc)
finalDate = extrapDate.shiftedBy(60.0*60*24*365*3)
first = False
sun = CelestialBodyFactory.getSun();
sun = PVCoordinatesProvider.cast_(sun)  # But we want the PVCoord interface
earth = OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS, 
                        Constants.WGS84_EARTH_FLATTENING, 
                        ITRF)
manu_1 = False
manu_2 = False
ahead = False
previous_lv_diff = 0

while (extrapDate.compareTo(finalDate) <= 0.0):  

    #pv = propagator.getPVCoordinates(extrapDate, inertialFrame)
    pv_state = propagator.propagate(extrapDate)
    state_list.append(pv_state)
    pv = pv_state.getPVCoordinates(inertialFrame)

    pos_sun = sun.getPVCoordinates(extrapDate, inertialFrame).getPosition()
    #print(pos_sun)
    #pos_sun_new = sun.getPosition(extrapDate, inertialFrame)
    beta_angle_new = 0.5 * math.pi - Vector3D.angle(pos_sun, pv.getMomentum())
    #print("vector_angle: %f" % math.degrees(Vector3D.angle(pos_sun, pv.getMomentum())))

    beta_list_sun.append(math.degrees(beta_angle_new))

    int_orbit = KeplerianOrbit(OrbitType.KEPLERIAN.convertType(pv_state.getOrbit()))
    #print(int_orbit)
    appending_rann = math.degrees(int_orbit.getRightAscensionOfAscendingNode())
    if appending_rann < 0:
        appending_rann = appending_rann + 360
    raan_list.append(appending_rann)
    i_list.append(math.degrees(int_orbit.getI()))
    sma_list.append(int_orbit.getA())
    ecc_list.append(int_orbit.getE())
    aop_list.append(math.degrees(int_orbit.getPerigeeArgument()))
    lv_list.append(math.degrees(int_orbit.getMeanAnomaly()))

    derived_pv_state = derivedpropagator.propagate(extrapDate)
    derived_state_list.append(derived_pv_state)

    derived_int_orbit = KeplerianOrbit(OrbitType.KEPLERIAN.convertType(derived_pv_state.getOrbit()))

    derived_pv = derived_pv_state.getPVCoordinates(inertialFrame)

    pos_sun = sun.getPVCoordinates(extrapDate, inertialFrame).getPosition()
    #pos_sun_new = sun.getPosition(extrapDate, inertialFrame)
    beta_angle_new = 0.5 * math.pi - Vector3D.angle(pos_sun, derived_pv.getMomentum())

    derived_beta_list_sun.append(math.degrees(beta_angle_new))
    dist = Vector3D.distance(pv_state.getPosition(), derived_pv_state.getPosition())/1000
    dist_list.append(dist)

    if extrapDate.compareTo(pair_2_start_date) >= 0.0:
        seed_2_state = pair_2_1.propagate(extrapDate)
        seed_2_state_list.append(seed_2_state)


        pair_2_orbit = KeplerianOrbit(OrbitType.KEPLERIAN.convertType(seed_2_state.getOrbit()))
        pair_2_sma_list.append(pair_2_orbit.getA())

        derived_2_state = pair_2_2.propagate(extrapDate)
        derived_2_state_list.append(derived_2_state)

        A_C_list.append(Vector3D.distance(seed_2_state.getPosition(), pv_state.getPosition())/1000)
        A_D_list.append(Vector3D.distance(derived_2_state.getPosition(), pv_state.getPosition())/1000) # D-A

        C_D_list.append(Vector3D.distance(derived_2_state.getPosition(), seed_2_state.getPosition())/1000) #D-C
        B_D_list.append(Vector3D.distance(derived_pv_state.getPosition(), derived_2_state.getPosition())/1000) #D-B
        B_C_list.append(Vector3D.distance(derived_pv_state.getPosition(), seed_2_state.getPosition())/1000) #C-B

    else:
        A_C_list.append(0) # C - A
        B_C_list.append(0) # C - B
        B_D_list.append(0) # D - B
        A_D_list.append(0) # D -A
        C_D_list.append(0) # C-D
        pair_2_sma_list.append(0)
    
    pair_2_date.append(absolutedate_to_datetime(extrapDate))
    date.append(absolutedate_to_datetime(extrapDate))

    print(extrapDate, end="\r")
    extrapDate = extrapDate.shiftedBy(24*3600.0)

header = ['date', 'distance(A-B)', 'distance(A-C)', 'distance(A-D)', 'distance(B-C)', 'distance(B-D)', 'distance(C-D)', 'beta_angle_sun', 'derived_beta_angle_sun', 'sma', 'inc', 'ecc', 'raan', 'aop', 'lv']
with open('output_sso/pair_2_bandwagon.csv', 'w') as f:
    writer = csv.writer(f)

    writer.writerow(header)
    for j in range(len(date)):
        writer.writerow([date[j], dist_list[j], A_C_list[j], A_D_list[j], B_C_list[j], B_D_list[j], C_D_list[j], beta_list_sun[j], derived_beta_list_sun[j], sma_list[j], i_list[j], ecc_list[j], raan_list[j], aop_list[j], lv_list[j]])

    f.close()
#print("--------------------------")
#print(OrbitType.KEPLERIAN.convertType(state_list[0].getOrbit()))

#print(OrbitType.KEPLERIAN.convertType(derived_state_list[0].getOrbit()))
#print("--------------------------")    
#print(OrbitType.KEPLERIAN.convertType(state_list[-1].getOrbit()))

#print(OrbitType.KEPLERIAN.convertType(derived_state_list[-1].getOrbit()))

print("SEED 2 ORBIT-----------------------")
print(OrbitType.KEPLERIAN.convertType(seed_2_state_list[1].getOrbit()))

print("DERIVED 2 ORBIT-----------------------")
print(OrbitType.KEPLERIAN.convertType(derived_2_state_list[1].getOrbit()))

val, idx = min((val, idx) for (idx, val) in enumerate(beta_list_sun))

#print(date[idx])
#print(OrbitType.KEPLERIAN.convertType(state_list[idx].getOrbit()))
#print(beta_list_sun[idx])

val, idx = max((val, idx) for (idx, val) in enumerate(beta_list_sun))

#print(date[idx])
#print(OrbitType.KEPLERIAN.convertType(state_list[idx].getOrbit()))
#print(beta_list_sun[idx])

#for i in range(len(date)):
#    if dist_list[i] <= 1000:
#        print("date: %s, dist: %f" % (date[i], dist_list[i]))

fig , axs = plt.subplots(4,3)
fig.suptitle(str(start_date) + "three years sso")
axs[0,0].plot(date, beta_list_sun)
axs[0,0].plot(date, derived_beta_list_sun)
axs[0,0].set_title("beta angle")
axs[0,1].plot(date, sma_list)
axs[0,1].set_title("sma_list")
axs[0,2].plot(date, i_list)
axs[0,2].set_title("i_list")
axs[1,0].plot(date, ecc_list)
axs[1,0].set_title("ecc_list")
axs[1,1].plot(date, raan_list)
axs[1,1].set_title("raan_list")
axs[1,2].plot(date, aop_list)
axs[1,2].set_title("aop_list")
axs[2,0].plot(date, lv_list)
axs[2,0].set_title("lv_list")
axs[2,1].plot(date, dist_list)
axs[2,1].set_title("dist_list")
axs[2,2].plot(pair_2_date, A_C_list)
axs[2,2].set_title("A_C_list")
axs[3,0].plot(pair_2_date, B_C_list)
axs[3,0].set_title("pair_2_sma_list")
axs[3,1].plot(pair_2_date, A_D_list)
axs[3,1].set_title("A_D_list")
axs[3,2].plot(pair_2_date, B_D_list)
axs[3,2].set_title("pair_2_dist_list_intra")
plt.show()