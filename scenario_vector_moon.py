import orekit
import numpy as np
import csv
import sys, getopt
from poliastro.earth.util import raan_from_ltan
from astropy.time import Time
from astropy import units as u


from orekit.pyhelpers import setup_orekit_curdir, absolutedate_to_datetime,datetime_to_absolutedate
from org.orekit.orbits import OrbitType, KeplerianOrbit
from org.orekit.bodies import OneAxisEllipsoid, GeodeticPoint
from org.orekit.time import AbsoluteDate, TimeScalesFactory
from org.orekit.utils import  IERSConventions, PVCoordinatesProvider, Constants
from org.orekit.frames import FramesFactory, TopocentricFrame
from org.orekit.propagation.analytical.tle import TLE
from org.orekit.propagation.analytical.tle.generation import FixedPointTleGenerationAlgorithm
from org.orekit.bodies import CelestialBodyFactory
from org.hipparchus.geometry.euclidean.threed import Vector3D
from org.orekit.estimation.iod import IodGooding
from org.orekit.geometry.fov import CircularFieldOfView
from org.orekit.attitudes import NadirPointing

from lib.helper_functions import get_abs_vel, mean_motion, distance_between_two, get_beta_angle_alternate, get_beta_angle
from lib.propagator import set_up_prop, get_ecef

from orekit import JArray_double
import math
import matplotlib.pyplot as plt

vm = orekit.initVM()
setup_orekit_curdir()
utc = TimeScalesFactory.getUTC()
satellite_mass = 16
start_date = AbsoluteDate(2025, 6, 1, 0, 0, 00.000, utc)
astro_star_date = Time('2025-06-01T00:00:00.000000000', format='isot', scale='utc')

ra = 600.14 * 1000        #  Apogee
rp = 600.139999 * 1000         #  Perigee
i = math.radians(97.7532)      # inclination (SSO)
#i = math.radians(23.5)
omega = math.radians(0.0)   # perigee argument
accurate_raan_sso_1 = raan_from_ltan(astro_star_date, ltan=10.5 * u.hourangle).value - 360
raan = math.radians(accurate_raan_sso_1)  # right ascension of ascending node (SSO)
#raan = math.radians(0.0)
lv = math.radians(0.0)    # True anomaly


# create point on earth frame
ITRF = FramesFactory.getITRF(IERSConventions.IERS_2010, True)
earth = OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS, 
                        Constants.WGS84_EARTH_FLATTENING, 
                        ITRF)

longitude = math.radians(103.8198)
latitude  = math.radians(1.3521)
altitude  = 125.0
station = GeodeticPoint(latitude, longitude, altitude)
station_frame = TopocentricFrame(earth, station, "CRISP")

cameraDirection = Vector3D(0.0, -1.0, 0.6)
camerafov = CircularFieldOfView(cameraDirection, math.radians(30), 0.0)

epochDate = start_date
#epochDate = AbsoluteDate(2024, 4, 1, 0, 0, 00.000, utc)
initialDate = epochDate

inertialFrame = FramesFactory.getEME2000()

ITRF = FramesFactory.getITRF(IERSConventions.IERS_2010, True)

propagator, _ = set_up_prop(i, omega, raan, lv, epochDate, inertialFrame, ITRF, rp=rp, ra=ra)

print(propagator.getInitialState().getOrbit())

propagator.setAttitudeProvider(NadirPointing(inertialFrame, earth))

extrapDate = start_date

state_list = []
tle_state_list = []

pv_list = []
pos = []

vel = []

vel_abs = []
acc = []

date = []
dist_list = []

print("SEED ORBIT-------------------------------")
print(propagator.getInitialState().getOrbit())

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

fov_list = []
fov_list_moon = []
fov_list_terra = []


finalDate = extrapDate.shiftedBy(60.0*60*24*30)
first = False
sun = CelestialBodyFactory.getSun()
sun = PVCoordinatesProvider.cast_(sun)
moon = CelestialBodyFactory.getMoon()
moon = PVCoordinatesProvider.cast_(moon)
terra = CelestialBodyFactory.getEarth()
terra = PVCoordinatesProvider.cast_(terra)

while (extrapDate.compareTo(finalDate) <= 0.0):  

    #pv = propagator.getPVCoordinates(extrapDate, inertialFrame)
    pv_state = propagator.propagate(extrapDate)

    pv = pv_state.getPVCoordinates(inertialFrame)

    #pos_earth = terra.getPVCoordinates(extrapDate, inertialFrame).getPosition()

    pos_sun = sun.getPVCoordinates(extrapDate, inertialFrame).getPosition()

    pos_moon = moon.getPVCoordinates(extrapDate, inertialFrame).getPosition()

    beta_angle_new = 0.5 * math.pi - Vector3D.angle(pos_sun, pv.getMomentum())

    beta_list_sun.append(math.degrees(beta_angle_new))

    fov_center = pv_state.toStaticTransform().getInverse().transformVector(camerafov.getCenter())
    #print(pv_state.toStaticTransform().getInverse())
    #print(fov_center)
    FOV_angle = Vector3D.angle(pos_sun, fov_center)
    FOV_angle_moon = Vector3D.angle(pos_moon, fov_center)
    #FOV_angle_terra = Vector3D.angle(pos_earth, fov_center)
    fov_list.append(math.degrees(FOV_angle))
    fov_list_moon.append(math.degrees(FOV_angle_moon))
    #fov_list_terra.append(math.degrees(FOV_angle_terra))

    date.append(absolutedate_to_datetime(extrapDate))
    print(extrapDate, end="\r")
    extrapDate = extrapDate.shiftedBy(60.0)

header = ['date', 'beta_angle_sun', 'fov-to-sun', 'fov-to-moon', 'fov-to-terra']
with open('output_sso/scenario_vector_moon.csv', 'w') as f:
    writer = csv.writer(f)

    writer.writerow(header)
    for j in range(len(date)):
        0
        writer.writerow([date[j], beta_list_sun[j], fov_list[j], fov_list_moon[j]])

    f.close()

print("--------------------------")
#print(OrbitType.KEPLERIAN.convertType(state_list[0].getOrbit()))

#print(OrbitType.KEPLERIAN.convertType(derived_state_list[0].getOrbit()))
print("--------------------------")    
#print(OrbitType.KEPLERIAN.convertType(state_list[-1].getOrbit()))

#print(OrbitType.KEPLERIAN.convertType(derived_state_list[-1].getOrbit()))


#val, idx = min((val, idx) for (idx, val) in enumerate(derived_beta_list_sun))

#print(date[idx])
#print(OrbitType.KEPLERIAN.convertType(derived_state_list[idx].getOrbit()))
#print(derived_beta_list_sun[idx])

#val, idx = max((val, idx) for (idx, val) in enumerate(derived_beta_list_sun))

#print(date[idx])
#print(OrbitType.KEPLERIAN.convertType(derived_state_list[idx].getOrbit()))
#print(derived_beta_list_sun[idx])
blind_counter = 0
for i in range(len(date)):
    if fov_list_moon[i] < 30:
        blind_counter = blind_counter + 1

print(blind_counter)


fig , axs = plt.subplots(4,1)
fig.suptitle(str(start_date) + "three years sso")
axs[0].plot(date, beta_list_sun)
axs[0].set_title("beta angle")
axs[1].plot(date, fov_list)
axs[1].set_title("fov angle")
axs[2].plot(date, fov_list_moon)
axs[2].set_title("fov angle moon")
#axs[3].plot(date, fov_list_terra)
#axs[3].set_title("fov angle terra")


plt.show()