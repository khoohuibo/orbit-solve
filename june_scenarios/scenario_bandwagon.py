import orekit
import numpy as np
import csv
import sys, getopt
from poliastro.earth.util import raan_from_ltan
from astropy.time import Time
from astropy import units as u


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
satellite_mass = 16
start_date = AbsoluteDate(2026, 1, 1, 0, 0, 00.000, utc)
astro_star_date = Time('2026-01-01T00:00:00.000000000', format='isot', scale='utc')

ra = 600.14 * 1000        #  Apogee
rp = 600.139999 * 1000         #  Perigee
i = math.radians(45.1)      # inclination (SSO)
omega = math.radians(0.0)   # perigee argument
#accurate_raan_sso_1 = raan_from_ltan(astro_star_date, ltan=10.5 * u.hourangle).value - 360
raan = math.radians(0.0)  # right ascension of ascending node (SSO)
lv = math.radians(0.0)    # True anomaly

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
#print("HERE HERE HERE!")
#print(OrbitType.KEPLERIAN.convertType(state_list[1].getOrbit()))
myTLE = TLE.stateToTLE(state_list[138], tle_first_guess, FixedPointTleGenerationAlgorithm())
#print(myTLE)
#print(OrbitType.KEPLERIAN.convertType(state_list[3].getOrbit()))
#print(date)
derivedOrbit_alt, seedOrbit, newEpoch, t_delta = get_ecef(date, absolutedate_to_datetime(epochDate), vel_abs, pos, vel, pv_list, inertialFrame, ITRF, inclination=i)

#print(derivedOrbit_alt)
#print(seedOrbit)

#print(t_delta)

sma = derivedOrbit_alt.getA()
ecc = derivedOrbit_alt.getE()
inc = derivedOrbit_alt.getI()
aop = derivedOrbit_alt.getPerigeeArgument()
raan_alt = derivedOrbit_alt.getRightAscensionOfAscendingNode()
lv_new = seedOrbit.getMeanAnomaly() + t_delta * mean_motion(satellite_mass, 600)
#derivedOrbit_alt.getTrueAnomaly()#
propagator, _ = set_up_prop(i, omega, raan, lv, epochDate, inertialFrame, ITRF,rp=rp, ra=ra, max_drag=False)

max_dist_list = []


raan_new = raan_alt 
derivedpropagator, _ = set_up_prop(inc, aop, raan_new, lv_new, newEpoch, inertialFrame, ITRF, rp=rp, ra=ra, max_drag=False)
#derivedpropagator, _ = set_up_prop(i, omega, raan, lv, epochDate, inertialFrame, ITRF, rp=rp, ra=ra, max_drag=True)
print("SEED ORBIT-------------------------------")
print(propagator.getInitialState().getOrbit())

print("DERIVED ORBIT-------------------------------")
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
i_list = []
sma_list = []
ecc_list = []
aop_list = []
lv_list = []

extrapDate = start_date
#extrapDate = AbsoluteDate(2024, 4, 1, 0, 0, 00.000, utc)
finalDate = extrapDate.shiftedBy(60.0*60*24*365*3)
first = False
sun = CelestialBodyFactory.getSun();
sun = PVCoordinatesProvider.cast_(sun)  # But we want the PVCoord interface

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
    sma_list.append(math.degrees(int_orbit.getA()))
    ecc_list.append(math.degrees(int_orbit.getE()))
    aop_list.append(math.degrees(int_orbit.getPerigeeArgument()))
    lv_list.append(math.degrees(int_orbit.getMeanAnomaly()))
    #print(math.degrees(int_orbit.getRightAscensionOfAscendingNode()))
    #beta_angle = get_beta_angle(absolutedate_to_datetime(extrapDate), int_orbit.getRightAscensionOfAscendingNode(), int_orbit.getI())
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
    print(extrapDate, end="\r")
    extrapDate = extrapDate.shiftedBy(24*3600.0)

header = ['date', 'distance', 'beta_angle_sun', 'derived_beta_angle_sun', 'sma', 'inc', 'ecc', 'raan', 'aop', 'lv']
with open('output_sso/scenario_3.csv', 'w') as f:
    writer = csv.writer(f)

    writer.writerow(header)
    for j in range(len(date)):
        writer.writerow([date[j], dist_list[j], beta_list_sun[j], derived_beta_list_sun[j], sma_list[j], i_list[j], ecc_list[j], raan_list[j], aop_list[j], lv_list[j]])

    f.close()
print("--------------------------")
print(OrbitType.KEPLERIAN.convertType(state_list[0].getOrbit()))

print(OrbitType.KEPLERIAN.convertType(derived_state_list[0].getOrbit()))
print("--------------------------")    
print(OrbitType.KEPLERIAN.convertType(state_list[-1].getOrbit()))

print(OrbitType.KEPLERIAN.convertType(derived_state_list[-1].getOrbit()))


val, idx = min((val, idx) for (idx, val) in enumerate(derived_beta_list_sun))

print(date[idx])
print(OrbitType.KEPLERIAN.convertType(derived_state_list[idx].getOrbit()))
print(derived_beta_list_sun[idx])

val, idx = max((val, idx) for (idx, val) in enumerate(derived_beta_list_sun))

print(date[idx])
print(OrbitType.KEPLERIAN.convertType(derived_state_list[idx].getOrbit()))
print(derived_beta_list_sun[idx])

for i in range(len(date)):
    if dist_list[i] <= 1000:
        print("date: %s, dist: %f" % (date[i], dist_list[i]))

fig , axs = plt.subplots(3,3)
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

#plt.plot(date, dist_list)
#plt.plot(date, beta_list_sun)
#plt.plot(date, beta_list)
#plt.plot(date, beta_list_alt)
#plt.plot(date, derived_beta_list_sun)
#plt.legend(['1', '2', '3', '4']) 
#plt.plot(date, raan_list)
plt.show()