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
from org.orekit.utils import  IERSConventions, PVCoordinatesProvider, Constants
from org.orekit.frames import FramesFactory, TopocentricFrame
from org.orekit.propagation.analytical.tle import TLE
from org.orekit.propagation.analytical.tle.generation import FixedPointTleGenerationAlgorithm
from org.orekit.bodies import CelestialBodyFactory, GeodeticPoint, OneAxisEllipsoid
from org.hipparchus.geometry.euclidean.threed import Vector3D
from org.orekit.propagation.events import ElevationDetector, EventsLogger
from org.orekit.propagation.events.handlers import ContinueOnEvent

from lib.helper_functions import get_abs_vel, mean_motion, distance_between_two, get_beta_angle_alternate, get_beta_angle
from lib.propagator import set_up_prop, get_ecef

from orekit import JArray_double
import math
import matplotlib.pyplot as plt

vm = orekit.initVM()
setup_orekit_curdir()
utc = TimeScalesFactory.getUTC()
satellite_mass = 16
start_date = AbsoluteDate(2025, 10, 1, 0, 0, 00.000, utc)
astro_star_date = Time('2025-10-01T00:00:00.000000000', format='isot', scale='utc')

ra = 600.14 * 1000        #  Apogee
rp = 600.139999 * 1000         #  Perigee
i = math.radians(97.7532)      # inclination (SSO)
omega = math.radians(0.0)   # perigee argument
accurate_raan_sso_1 = raan_from_ltan(astro_star_date, ltan=10.5 * u.hourangle).value - 360
raan = math.radians(accurate_raan_sso_1)  # right ascension of ascending node (SSO)
lv = math.radians(0.0)    # True anomaly

epochDate = start_date
#epochDate = AbsoluteDate(2024, 4, 1, 0, 0, 00.000, utc)
initialDate = epochDate

inertialFrame = FramesFactory.getEME2000()

ITRF = FramesFactory.getITRF(IERSConventions.IERS_2010, True)

propagator, _ = set_up_prop(i, omega, raan, lv, epochDate, inertialFrame, ITRF, rp=rp, ra=ra)

print(propagator.getInitialState().getOrbit())

ITRF = FramesFactory.getITRF(IERSConventions.IERS_2010, True)
earth = OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS, 
                        Constants.WGS84_EARTH_FLATTENING, 
                        ITRF)

longitude = math.radians(121.534608)
latitude  = math.radians(25.04366)
altitude  = 125.0
station = GeodeticPoint(latitude, longitude, altitude)
station_frame = TopocentricFrame(earth, station, "RPTK-8U")

elevation_detector = ElevationDetector(station_frame).withConstantElevation(10.0).withHandler(ContinueOnEvent())
logger = EventsLogger()
logged_detector = logger.monitorDetector(elevation_detector)

propagator.addEventDetector(logged_detector)

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

el_store_1 = []
az_store_1 = []
t_store_1 = []

extrapDate = start_date
#extrapDate = AbsoluteDate(2024, 4, 1, 0, 0, 00.000, utc)
finalDate = extrapDate.shiftedBy(60.0*60*24)
first = False
sun = CelestialBodyFactory.getSun()
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

    pos_sun = sun.getPVCoordinates(extrapDate, inertialFrame).getPosition()
    #pos_sun_new = sun.getPosition(extrapDate, inertialFrame)

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
    date.append(absolutedate_to_datetime(extrapDate))
    print(extrapDate, end="\r")
    extrapDate = extrapDate.shiftedBy(60.0)

events = logger.getLoggedEvents()
if (events.size() == 0):
    print("huh?")
print(events.size())
print(el_store_1)
raise Exception

header = ['date', 'beta_angle_sun', 'derived_beta_angle_sun', 'sma', 'inc', 'ecc', 'raan', 'aop', 'lv']
with open('output_sso/scenario_beta.csv', 'w') as f:
    writer = csv.writer(f)

    writer.writerow(header)
    for j in range(len(date)):
        writer.writerow([date[j], beta_list_sun[j], derived_beta_list_sun[j], sma_list[j], i_list[j], ecc_list[j], raan_list[j], aop_list[j], lv_list[j]])

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