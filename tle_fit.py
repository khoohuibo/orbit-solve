import orekit
import numpy as np
import csv
import math
import string
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import getopt
import sys

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

from lib.helper_functions import get_abs_vel, distance_between_two, checksum
from lib.propagator import set_up_prop, get_ecef, process_tle

def main(argv):

    verbose = False
    input_file = None
    compare_file = None
    diff = None
    # we need the location of the input TLE file, followed by 
    try:
        opts, args = getopt.getopt(argv,"hvi:c:d:",["help", "verbose", "input=", "compare=", "diff="])
    except getopt.GetoptError:
        print("tle_fit.py --input=<TLE/FILELOCATION> --compare=<TLE/FILELOCATION - default empty is propagate> --diff=<seconds-offset expected> --help --verbose")
        sys.exit(2)
    
    for opt, arg in opts:
        if opt in ('-h', "--help"):
            print("""
            TLE fitting helper
            1. Input TLE file location
            2. TLE to compare to. If empty, will propagate and create in default location
            3. diff - used to check any compared TLE. If propagated, then diff determines offset along track
            """)
            sys.exit()
        elif opt in ('-i', "--input"):
            input_file = arg
            print("input selected = %s" % arg)
        elif opt in ('-c', "--compare"):
            compare_file = arg
            print("compare selected = %s" % arg)
        elif opt in ('-d', "--diff"):
            diff = arg
            print("time difference = %s" % arg)
        elif opt in ('-v', "--verbose"):
            verbose = True
            print("verbosity increased!")
        else:
            print("tle_fit.py --input=<TLE/FILELOCATION> --compare=<TLE/FILELOCATION - default empty is propagate> --diff=<seconds-offset expected> --help")
            sys.exit()
    
    if input_file is None:
        print("defaulting to standard TLE file!")
        input_file = 'tle/input_2.tle'

    if verbose is False:
        print("Print Messages will not be visible if you do not enable verbosity!")
    
    if compare_file is not None:
        print("executing Comparison Mode! Input TLE will be compared with Compare TLE!")

    if diff is None:
        print("defaulting to standard time difference of 60 seconds")
        diff = 60 
    

    # open the old TLE
    propagator, station_frame, inertialFrame, tle_epoch, mytle = process_tle(input_file)

    extrapDate = tle_epoch
    time_span_sec = 24.0*60*60
    finalDate = extrapDate.shiftedBy(time_span_sec) #seconds

    state_list = []

    pv_list = []
    pos = []
    vel = []
    vel_abs = []
    acc = []

    date = []

    el_store_1=[]
    az_store_1 = []
    t_store_1 = []

    while (extrapDate.compareTo(finalDate) <= 0.0):  
        pv_state = propagator.propagate(extrapDate)
        state_list.append(pv_state)
        pv = pv_state.getPVCoordinates()
        pv_list.append(pv)

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
        #print(el_tmp, extrapDate)
        el_store_1.append(el_tmp)

        az_tmp = station_frame.getAzimuth(pv.getPosition(),
                        inertialFrame,
                        extrapDate)#*180.0/math.pi
        az_store_1.append(az_tmp)
        #print extrapDate, pos_tmp, vel_tmp
        t_store_1.append(absolutedate_to_datetime(extrapDate))
        extrapDate = extrapDate.shiftedBy(30.0)

    if compare_file is None:
        ITRF = FramesFactory.getITRF(IERSConventions.IERS_2010, True)
        derivedOrbit_alt, _, newEpoch, _ = get_ecef(date, absolutedate_to_datetime(mytle.getDate()), vel_abs, pos, vel, pv_list, inertialFrame, ITRF, inclination=mytle.getI(), t_delta_0=diff)

        derivedpropagator, _ = set_up_prop(0,0,0,0,newEpoch,inertialFrame,ITRF,rp=0, ra=0, initialOrbit=derivedOrbit_alt)
    
    else:
        derivedpropagator, _, _, _, _ = process_tle(compare_file)

    extrapDate = tle_epoch
    finalDate = extrapDate.shiftedBy(time_span_sec) # seconds
    derived_state_list = []
    el_store_2=[]
    az_store_2 = []
    t_store_2 = []

    while (extrapDate.compareTo(finalDate) <= 0.0):
        derived_pv_state = derivedpropagator.propagate(extrapDate)
        derived_state_list.append(derived_pv_state)
        derived_pv = derived_pv_state.getPVCoordinates()

        date.append(absolutedate_to_datetime(extrapDate))
        
        el_tmp = station_frame.getElevation(derived_pv.getPosition(),
                        inertialFrame,
                        extrapDate)*180.0/math.pi
        #print(el_tmp, extrapDate)
        el_store_2.append(el_tmp)

        az_tmp = station_frame.getAzimuth(derived_pv.getPosition(),
                        inertialFrame,
                        extrapDate)#*180.0/math.pi
        az_store_2.append(az_tmp)
        #print extrapDate, pos_tmp, vel_tmp
        t_store_2.append(absolutedate_to_datetime(extrapDate))

        print(extrapDate, end="\r")
        extrapDate = extrapDate.shiftedBy(30.0)

    peaks_2, _ = find_peaks(el_store_2, height=0)
    peaks, _ = find_peaks(el_store_1, height=0)

    width = 10
    error_index = []
    max_az_offset = math.radians(0.5)
    max_el_offset = 0.2
    max_t_az = 0
    max_t_el = 0


    # comparison of track with each other
    for index in range(len(peaks)):
        # the TCA point will be different for both of them. 
        # generally, 1 is ahead of 2 in track. e.g------sat 2 ----------------> sat 1
        for i in range(width):
            try:
                az_offset = math.pi - abs(abs(az_store_1[peaks[index] + i] - az_store_2[peaks_2[index] + i]) - math.pi)
                el_offset = el_store_1[peaks[index] + i] - el_store_2[peaks_2[index]+i]
                        #print(az_offset, el_offset, t_store_1[peaks[index]+i].value, t_store_1[peaks_2[index]+i].value)
                if az_offset > max_az_offset:
                    max_az_offset = az_offset
                    max_t_az = t_store_1[peaks[index]+i]
                if abs(el_offset) > max_el_offset:
                    max_el_offset = abs(el_offset) 
                    max_t_el = t_store_1[peaks[index]+i]
                if az_offset > 0.015 or el_offset > 0.2:
                    error_index.append(index+i)
                    print("---------Beyond Tolerance-----------")
                    print("AZ-OFFSET, AZ_OLD, AZ_NEW")
                    print(az_offset, az_store_1[peaks[index] +i], az_store_2[peaks_2[index]+i], i)
                    print("EL_OFFSET, EL_OLD, EL_NEW")
                    print(el_offset, el_store_1[peaks[index] +i], el_store_2[peaks_2[index]+i], i)
                    print("T_OLD, T_NEW")
                    print(t_store_1[peaks[index]+i], t_store_1[peaks_2[index]+i])
                    print("---------%s----------" % t_store_1[peaks[index]+i])
            except IndexError:
                continue

    if len(error_index) > 0:
        print("Divergence found! Check gpredict and compare! Total number of divergences: %d" % len(error_index))
        print("list of error points------")
        print("Fraction over simulation: %0.2f" % (len(error_index)/time_span_sec))
        print("Maximum Azimuth Offset: %f %s\r\n" % (max_az_offset, max_t_az))
        print("maximum elevation offset: %f  %s\r\n" % (max_el_offset, max_t_el))

        #for i in error_index:
        #    print("error_index:[%d] - %s" % (i, t_store_1[i].value))
    """
    Tolerance for automated TLE
        - 0.5 degrees azimuth
        - 0.2 degrees elevation
        - Number of divergences/Number of visible seconds = less than 5% 
    """
    newTLE = TLE.stateToTLE(derived_state_list[-1], mytle, FixedPointTleGenerationAlgorithm())
    print(newTLE)

    # this where plotting begins
    fig = plt.figure(figsize=(10,10))
    ax_2 = plt.subplot(3, 2, (1,2))

    ax_2.plot(t_store_2, el_store_1, '-gD', markevery=peaks)
    ax_2.set_title("All peaks/passes")

    for index in range(len(peaks)-4, len(peaks)):
        ax = plt.subplot(3, 2, index-(len(peaks)-4)+3, polar=True)
        
        ax.set_theta_zero_location('N')

        ax.set_theta_direction(-1)

        ax.set_rlim(bottom=90, top=0)

        ax.plot(az_store_2[peaks_2[index]-width:peaks_2[index]+width], el_store_2[peaks_2[index]-width:peaks_2[index]+width], color='b', zorder=1, label='OLD TLE')

        ax.plot(az_store_1[peaks[index]-width:peaks[index]+width], el_store_1[peaks[index]-width:peaks[index]+width], linestyle='--', color='r', linewidth=2, zorder=2, label='NEW TLE')
        
        ax.legend()
        
        ax.grid(True)
        
        ax.set_title("%s \r\n" % (t_store_1[peaks[index]]))


    plt.subplots_adjust(bottom=0.038, left=0.022 ,right=0.992, top=0.967, hspace=0.39, wspace=0.01)
    plt.show()


if __name__ == "__main__":
    main(sys.argv[1:])