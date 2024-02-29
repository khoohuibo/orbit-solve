import orekit
import numpy as np
import csv
import sys, getopt



from orekit.pyhelpers import setup_orekit_curdir, absolutedate_to_datetime,datetime_to_absolutedate
from org.orekit.orbits import OrbitType, KeplerianOrbit
from org.orekit.time import AbsoluteDate, TimeScalesFactory
from org.orekit.utils import  IERSConventions
from org.orekit.frames import FramesFactory
from org.orekit.propagation.analytical.tle import TLE
from org.orekit.propagation.analytical.tle.generation import FixedPointTleGenerationAlgorithm

from lib.helper_functions import get_abs_vel, mean_motion, distance_between_two, get_beta_angle, get_beta_angle_alternate
from lib.propagator import set_up_prop, get_ecef

from orekit import JArray_double
import math
import matplotlib.pyplot as plt

vm = orekit.initVM()
setup_orekit_curdir()
utc = TimeScalesFactory.getUTC()

ra = 600.01 *  1000         #  Apogee
rp = 600 * 1000         #  Perigee
i = math.radians(97.6)      # inclination
omega = math.radians(0.0)   # perigee argument
raan = math.radians(0.0)  # right ascension of ascending node
lv = math.radians(0.0)    # True anomaly

epochDate = AbsoluteDate(2025, 4, 1, 0, 0, 00.000, utc)
initialDate = epochDate

inertialFrame = FramesFactory.getEME2000()

ITRF = FramesFactory.getITRF(IERSConventions.IERS_2010, True)

propagator = set_up_prop(rp, ra, i, omega, raan, lv, epochDate, inertialFrame, ITRF)

print(propagator.getInitialState().getOrbit())

extrapDate = AbsoluteDate(2025, 4, 1, 0, 0, 00.000, utc)
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

raan_list = []
max_dist_list = []

def raan_flex(k, factor):
    raan_new = raan_alt + math.radians(k/factor)
    derivedpropagator = set_up_prop(rp, ra, inc, aop, raan_new, lv_new, newEpoch, inertialFrame, ITRF, a=sma, e=ecc)

    print(derivedpropagator.getInitialState().getOrbit())

    date = []
    dist_list = []

    state_list = []
    derived_state_list = []

    beta_list = []
    beta_list_alt = []

    extrapDate = AbsoluteDate(2025, 4, 1, 0, 0, 00.000, utc)
    finalDate = extrapDate.shiftedBy(60.0*60*24*365*3)

    while (extrapDate.compareTo(finalDate) <= 0.0):  

        #pv = propagator.getPVCoordinates(extrapDate, inertialFrame)
        pv_state = propagator.propagate(extrapDate)
        state_list.append(pv_state)

        int_orbit = KeplerianOrbit(OrbitType.KEPLERIAN.convertType(pv_state.getOrbit()))
        beta_angle = get_beta_angle(absolutedate_to_datetime(extrapDate), int_orbit.getRightAscensionOfAscendingNode(), int_orbit.getI())
        beta_angle_alt = get_beta_angle_alternate(absolutedate_to_datetime(extrapDate), int_orbit.getRightAscensionOfAscendingNode(), int_orbit.getI())
        beta_list.append(math.degrees(beta_angle))
        beta_list_alt.append(math.degrees(beta_angle_alt))
        print("beta_angle: %f, beta_angle_alt: %f" % (beta_angle, beta_angle_alt))

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

    myTLE = TLE.stateToTLE(state_list[3], tle_first_guess, FixedPointTleGenerationAlgorithm())
    print(myTLE)
    print(OrbitType.KEPLERIAN.convertType(state_list[3].getOrbit()))

    myTLE = TLE.stateToTLE(derived_state_list[3], tle_first_guess, FixedPointTleGenerationAlgorithm())
    print(myTLE)
    print(OrbitType.KEPLERIAN.convertType(derived_state_list[3].getOrbit()))

    header = ['date', 'Distance']
    with open('output_sso/raan_offset_%d_%d.csv' % (k, inc), 'w') as f:
        writer = csv.writer(f)

        writer.writerow(header)
        for j in range(len(date)):
            writer.writerow([date[j], dist_list[j]])

        f.close()

    np.savez("output_sso/raan_offset_%d_sso" % (k), date, dist_list)

    return raan_new, max(dist_list)


if __name__ == "__main__":
    try:
        opts, args = getopt.getopt(sys.argv[1:], "sef", ["start=", "end=", "factor="])
    except getopt.GetoptError:
        print("thread_dsst.py --start=<deg*10> --end=<deg*10> --factor=<0..100>")

    start = 0
    end = 1
    factor = 1

    for opt, arg in opts:
        if opt in ('-s', "--start"):
            start = int(arg)
            print("start = %d" % start)
        elif opt in ('-e', "--end"):
            end = int(arg)
            print("end = %d" % end)
        elif opt in ('-f', "--factor"):
            factor = int(arg)
            print("factor = %d" % factor)
    for k in range(start, end, 1):
        raan_output, maximum_dist = raan_flex(k, factor)
        raan_list.append(raan_output)
        max_dist_list.append(maximum_dist)

    header = ['RAAN', 'Maximum Distance']

    with open('output/raan_all_sso.csv', 'w', encoding='UTF8') as f:
        writer = csv.writer(f)

        # write the header
        writer.writerow(header)

        # write the data
        for l in range(len(raan_list)):
            writer.writerow([raan_list[l], max_dist_list[l]])

        f.close()

    np.savez("output/raan_all_sso", raan_list, max_dist_list)

    plt.scatter(raan_list, max_dist_list)
    
    plt.legend()
    plt.show()