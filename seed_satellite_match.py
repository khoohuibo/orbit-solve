import orekit
import numpy as np
import csv

vm = orekit.initVM()

from orekit.pyhelpers import setup_orekit_curdir, absolutedate_to_datetime,datetime_to_absolutedate
setup_orekit_curdir()

from org.orekit.orbits import KeplerianOrbit, OrbitType
from org.orekit.time import AbsoluteDate, TimeScalesFactory
from org.orekit.utils import Constants, IERSConventions, PVCoordinates, TimeStampedPVCoordinates
from org.orekit.frames import FramesFactory
from org.hipparchus.geometry.euclidean.threed import Vector3D
from org.orekit.propagation.analytical.tle import TLE
from org.orekit.propagation.analytical.tle import TLEPropagator
from org.orekit.propagation.analytical.tle.generation import FixedPointTleGenerationAlgorithm

from lib.helper_functions import get_abs_vel, mean_motion, distance_between_two
from lib.propagator import set_up_prop, get_ecef

from orekit import JArray_double
import math
import matplotlib.pyplot as plt

utc = TimeScalesFactory.getUTC()

ra = 600.01 *  1000         #  Apogee
rp = 600 * 1000         #  Perigee
i = math.radians(45.4)      # inclination
omega = math.radians(0.0)   # perigee argument
raan = math.radians(0.0)  # right ascension of ascending node
lv = math.radians(0.0)    # True anomaly

epochDate = AbsoluteDate(2025, 6, 1, 0, 0, 00.000, utc)
initialDate = epochDate

inertialFrame = FramesFactory.getEME2000()

ITRF = FramesFactory.getITRF(IERSConventions.IERS_2010, True)

propagator = set_up_prop(rp, ra, i, omega, raan, lv, epochDate, inertialFrame, ITRF)

print(propagator.getInitialState().getOrbit())

extrapDate = AbsoluteDate(2025, 6, 1, 0, 0, 00.000, utc)
finalDate = extrapDate.shiftedBy(60.0*10)

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

#print(date)
derivedOrbit_alt, seedOrbit, newEpoch = get_ecef(date, absolutedate_to_datetime(epochDate), vel_abs, pos, vel, pv_list, inertialFrame, ITRF, inclination=i)

#print(derivedOrbit_alt)
print(seedOrbit)

sma = derivedOrbit_alt.getA()
ecc = derivedOrbit_alt.getE()
inc = derivedOrbit_alt.getI()
aop = derivedOrbit_alt.getPerigeeArgument()
raan_alt = derivedOrbit_alt.getRightAscensionOfAscendingNode()
lv_new = derivedOrbit_alt.getTrueAnomaly()

propagator = set_up_prop(rp, ra, i, omega, raan, lv, epochDate, inertialFrame, ITRF)
def raan_flex(k):
    for k in range(-20, 100, 1):
    raan_alt = raan_alt + math.radians(-10)
    derivedpropagator = set_up_prop(rp, ra, inc, aop, raan_alt, lv_new, newEpoch, inertialFrame, ITRF, a=sma, e=ecc)

    print(derivedpropagator.getInitialState().getOrbit())

    date = []
    dist_list = []

    state_list = []
    derived_state_list = []

    extrapDate = AbsoluteDate(2025, 6, 1, 0, 0, 00.000, utc)
    finalDate = extrapDate.shiftedBy(60.0*60*24*365)

    while (extrapDate.compareTo(finalDate) <= 0.0):  

        #pv = propagator.getPVCoordinates(extrapDate, inertialFrame)
        pv_state = propagator.propagate(extrapDate)
        state_list.append(pv_state)

        derived_pv_state = derivedpropagator.propagate(extrapDate)
        derived_state_list.append(derived_pv_state)
        dist = distance_between_two(pv_state.getPosition(), derived_pv_state.getPosition())/1000
        dist_list.append(dist)

        date.append(absolutedate_to_datetime(extrapDate))
        print(extrapDate, end="\r")
        extrapDate = extrapDate.shiftedBy(24*3600.0)

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

    myTLE = TLE.stateToTLE(state_list[-1], tle_first_guess, FixedPointTleGenerationAlgorithm())
    print(myTLE)
    print(OrbitType.KEPLERIAN.convertType(state_list[-1].getOrbit()))

    myTLE = TLE.stateToTLE(derived_state_list[-1], tle_first_guess, FixedPointTleGenerationAlgorithm())
    print(myTLE)
    print(OrbitType.KEPLERIAN.convertType(derived_state_list[-1].getOrbit()))

    
    header = ['date', 'Distance']
    with open('output_45/raan_offset_%d_45_june_2025.csv' % k, 'w') as f:
        writer = csv.writer(f)

        writer.writerow(header)
        for j in range(len(date)):
            writer.writerow([date[j], dist_list[j]])

        f.close()

    np.savez("output_45/raan_offset_%d_45" % k, date, dist_list)

fig = plt.figure()


#print(dist_list)

ax2 = fig.add_subplot(222)
ax2.plot(date, dist_list)
ax2.title.set_text('derived distance')

plt.legend()
plt.show()