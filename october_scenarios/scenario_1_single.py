import orekit
import numpy as np
import csv
import sys, getopt
from poliastro.earth.util import raan_from_ltan
from astropy.time import Time
from astropy import units as u


from orekit.pyhelpers import setup_orekit_curdir, absolutedate_to_datetime
from org.orekit.orbits import OrbitType, KeplerianOrbit, PositionAngleType, OrbitType
from org.orekit.time import AbsoluteDate, TimeScalesFactory
from org.orekit.utils import  IERSConventions, PVCoordinatesProvider, Constants
from org.orekit.frames import FramesFactory
from org.orekit.bodies import CelestialBodyFactory, OneAxisEllipsoid
from org.hipparchus.geometry.euclidean.threed import Vector3D

from org.orekit.propagation.semianalytical.dsst import DSSTPropagator
from org.orekit.propagation.semianalytical.dsst.forces import DSSTAtmosphericDrag, DSSTZonal, DSSTSolarRadiationPressure, DSSTTesseral
from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
from org.orekit.propagation import SpacecraftState, PropagationType

from org.orekit.utils import Constants, IERSConventions
from org.orekit.frames import FramesFactory
from org.orekit.forces.gravity.potential import GravityFieldFactory
from org.orekit.forces.drag import DragForce,IsotropicDrag
from org.orekit.models.earth.atmosphere import NRLMSISE00
from org.orekit.models.earth.atmosphere.data import CssiSpaceWeatherData
from org.orekit.forces.radiation import IsotropicRadiationSingleCoefficient
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
raan = math.radians(344.920)
lv = math.radians(0.0)    # True anomaly

inertialFrame = FramesFactory.getEME2000()

ITRF = FramesFactory.getITRF(IERSConventions.IERS_2010, True)

a = (rp + ra + 2 * Constants.WGS84_EARTH_EQUATORIAL_RADIUS) / 2.0 
e = 1.0 - (rp + Constants.WGS84_EARTH_EQUATORIAL_RADIUS) / a

earth = OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS, 
                        Constants.WGS84_EARTH_FLATTENING, 
                        ITRF)
 
initialOrbit = KeplerianOrbit(a, e, i, omega, raan, lv,
                                    PositionAngleType.TRUE,
                                    inertialFrame, start_date, Constants.WGS84_EARTH_MU)


satellite_mass = 10.0  # The models need a spacecraft mass, unit kg.
initialState = SpacecraftState(initialOrbit, satellite_mass) 

minStep = 1.0#initialState.getKeplerianPeriod()
maxStep = 100. * initialState.getKeplerianPeriod()#100. * minStep
tolerances = DSSTPropagator.tolerances(1.0, initialOrbit)

integrator = DormandPrince853Integrator(minStep, maxStep, 
    JArray_double.cast_(tolerances[0]),  # Double array of doubles needs to be casted in Python
    JArray_double.cast_(tolerances[1]))

propagator = DSSTPropagator(integrator)
propagator.setInitialState(initialState, PropagationType.MEAN)

cswl = CssiSpaceWeatherData("SpaceWeather-All-v1.2.txt")

atmosphere = NRLMSISE00(cswl, CelestialBodyFactory.getSun(), earth)

cross_section = 0.32
cd = 2.2
isotropic_drag = IsotropicDrag(cross_section, cd)

drag_drivers = isotropic_drag.getDragParametersDrivers()

drag_force = DragForce(atmosphere, isotropic_drag)

gravityProvider = GravityFieldFactory.getUnnormalizedProvider(10, 10)
propagator.addForceModel(DSSTZonal(gravityProvider))

propagator.addForceModel(DSSTAtmosphericDrag(drag_force, gravityProvider.getMu()))

propagator.addForceModel(DSSTTesseral(earth.getBodyFrame(), Constants.WGS84_EARTH_ANGULAR_VELOCITY, gravityProvider))

rad = IsotropicRadiationSingleCoefficient(0.24, 1.0)
propagator.addForceModel(DSSTSolarRadiationPressure(CelestialBodyFactory.getSun(), earth, rad, gravityProvider.getMu()))

print(propagator.getInitialState().getOrbit())

date = []
dist_list = []

state_list = []


beta_list_sun = []
raan_list= []
i_list = []
sma_list = []

ecc_list = []
aop_list = []
lv_list = []

ap_list = []
daily_flux = []

extrapDate = start_date
finalDate = extrapDate.shiftedBy(60.0*60*24*365*3)
sun = CelestialBodyFactory.getSun()
sun = PVCoordinatesProvider.cast_(sun)  # But we want the PVCoord interface
cswl = CssiSpaceWeatherData("SpaceWeather-All-v1.2.txt")

while (extrapDate.compareTo(finalDate) <= 0.0):  


    pv_state = propagator.propagate(extrapDate)
    state_list.append(pv_state)
    pv = pv_state.getPVCoordinates(inertialFrame)

    pos_sun = sun.getPVCoordinates(extrapDate, inertialFrame).getPosition()
    beta_angle_new = 0.5 * math.pi - Vector3D.angle(pos_sun, pv.getMomentum())
    beta_list_sun.append(math.degrees(beta_angle_new))

    int_orbit = KeplerianOrbit(OrbitType.KEPLERIAN.convertType(pv_state.getOrbit()))

    appending_rann = math.degrees(int_orbit.getRightAscensionOfAscendingNode())
    if appending_rann < 0:
        appending_rann = appending_rann + 360
    raan_list.append(appending_rann)

    i_list.append(math.degrees(int_orbit.getI()))
    sma_list.append(int_orbit.getA()-6371000)
    ecc_list.append(math.degrees(int_orbit.getE()))
    aop_list.append(math.degrees(int_orbit.getPerigeeArgument()))
    lv_list.append(math.degrees(int_orbit.getMeanAnomaly()))

    ap_list.append(cswl.getAp(extrapDate))
    daily_flux.append(cswl.getDailyFlux(extrapDate))

    date.append(absolutedate_to_datetime(extrapDate))
    print(extrapDate, end="\r")
    extrapDate = extrapDate.shiftedBy(3600.0)

header = ['date', 'beta_angle_sun', 'sma', 'inc', 'ecc', 'raan', 'aop', 'lv', 'ap_list', 'daily_flux']
with open('output/scenario_single.csv', 'w') as f:
    writer = csv.writer(f)

    writer.writerow(header)
    for j in range(len(date)):
        writer.writerow([date[j], beta_list_sun[j], sma_list[j], i_list[j], ecc_list[j], raan_list[j], aop_list[j], lv_list[j], ap_list[j], daily_flux[j]])

    f.close()
print("--------------------------")
print(OrbitType.KEPLERIAN.convertType(state_list[0].getOrbit()))


print("--------------------------")    
print(OrbitType.KEPLERIAN.convertType(state_list[-1].getOrbit()))


print("-------------min seed-------------")

val, idx = min((val, idx) for (idx, val) in enumerate(beta_list_sun))

print(date[idx])
print(OrbitType.KEPLERIAN.convertType(state_list[idx].getOrbit()))
print(beta_list_sun[idx])
print(ap_list[idx])
print(daily_flux[idx])

print("-------------max seed-------------")

val, idx = max((val, idx) for (idx, val) in enumerate(beta_list_sun))

print(date[idx])
print(OrbitType.KEPLERIAN.convertType(state_list[idx].getOrbit()))
print(beta_list_sun[idx])
print(ap_list[idx])
print(daily_flux[idx])


fig , axs = plt.subplots(3,3)
fig.suptitle(str(start_date) + "three years sso")
axs[0,0].plot(date, beta_list_sun)
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

plt.show()