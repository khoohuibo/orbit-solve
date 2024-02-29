import numpy as np
import math

from org.orekit.orbits import KeplerianOrbit, PositionAngleType, OrbitType
from org.orekit.propagation.numerical import NumericalPropagator
from org.orekit.propagation.semianalytical.dsst import DSSTPropagator
from org.orekit.propagation.semianalytical.dsst.forces import DSSTAtmosphericDrag, DSSTZonal, DSSTSolarRadiationPressure, DSSTTesseral
from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
from org.orekit.propagation import SpacecraftState, PropagationType

from org.orekit.utils import Constants, IERSConventions
from org.orekit.frames import FramesFactory
from org.orekit.bodies import OneAxisEllipsoid, CelestialBodyFactory
from org.hipparchus.geometry.euclidean.threed import Vector3D
from org.orekit.forces.gravity.potential import GravityFieldFactory
from org.orekit.forces.gravity import HolmesFeatherstoneAttractionModel
from org.orekit.forces.drag import DragForce,IsotropicDrag
from org.orekit.models.earth.atmosphere import NRLMSISE00
from org.orekit.models.earth.atmosphere.data import CssiSpaceWeatherData
from org.orekit.utils import PVCoordinates
from org.orekit.forces.radiation import IsotropicRadiationSingleCoefficient, SolarRadiationPressure, KnockeRediffusedForceModel
from orekit.pyhelpers import datetime_to_absolutedate
from orekit import JArray_double

from lib.helper_functions import get_angle_rez, alternate_delta_t

def get_ecef(date, dat_0, vel_abs, pos, vel, pv_list, inertialFrame, ITRF, inclination, t_delta_0=False):
    We = 7.272 * 10 ** -5
    trip = False
    if t_delta_0 == False:
        t_delta_0 = alternate_delta_t(veloc=vel_abs[0], inclination=inclination)
    else:
        trip = True

    for k in range(len(date)):

        newEpoch = datetime_to_absolutedate(date[k])
        #print(newEpoch)

        dat = date[k]

        #print(dat, dat_0)

        #pog = float(dat.hour + dat.minute/60 + dat.second/3600)

        #lon_g0 = get_GMST(dat, pog, 0)

        # find time difference required.
        if trip == True:
            t_delta = t_delta_0
        else:
            t_delta = alternate_delta_t(veloc=vel_abs[k], inclination=inclination) 

        #comp_a = - lon_g0 - We * t_delta 
        comp_a = - We * t_delta 

        r_3 = np.array([[math.cos(comp_a), -math.sin(comp_a), 0],
                    [math.sin(comp_a), math.cos(comp_a), 0],
                    [0,0,1]])

        pos_eci = np.array(pos[k])

        vel_eci = np.array(vel[k])

        eciToEcf = inertialFrame.getTransformTo(ITRF, newEpoch)
        pvECF = eciToEcf.transformPVCoordinates(pv_list[k])
        ecfToEci = ITRF.getTransformTo(inertialFrame, newEpoch)

        pos_ecef_new = pvECF.getPosition()
        vel_ecef_new = pvECF.getVelocity()

        pos_ecef_new_tmp = np.array([pos_ecef_new.getX(), pos_ecef_new.getY(),pos_ecef_new.getZ()])
        #pos_ecef_new_tmp = r_3.dot(pos_ecef_new_tmp)

        vel_ecef_new_tmp = np.array([vel_ecef_new.getX(), vel_ecef_new.getY(),vel_ecef_new.getZ()])
        #vel_ecef_new_tmp = r_3.dot(vel_ecef_new_tmp)

        pos_ecef_alt = Vector3D(float(pos_ecef_new_tmp[0]), float(pos_ecef_new_tmp[1]), float(pos_ecef_new_tmp[2]))
        vel_ecef_alt = Vector3D(float(vel_ecef_new_tmp[0]), float(vel_ecef_new_tmp[1]), float(vel_ecef_new_tmp[2]))

        pv_ecef_alt = PVCoordinates(pos_ecef_alt, vel_ecef_alt)

        pv_eci_alt = ecfToEci.transformPVCoordinates(pv_ecef_alt)

        if (dat - dat_0).total_seconds() > t_delta_0 :
            print("-------%s-------" % date[k])
            print(pos_eci)
            print(vel_eci)

            print("********ECEF_ALT********")
            print(pos_ecef_new_tmp)
            print(vel_ecef_new_tmp)

            print("&&&&&&ECI_ALT&&&&&&")
            print(pv_eci_alt.getPosition())
            print(pv_eci_alt.getVelocity())

            seedOrbit = KeplerianOrbit(pv_list[k], inertialFrame, newEpoch, Constants.WGS84_EARTH_MU)
            derivedOrbit_alt = KeplerianOrbit(pv_eci_alt, inertialFrame, newEpoch, Constants.WGS84_EARTH_MU)
            #print(newEpoch, t_delta)
            break
    
    return derivedOrbit_alt, seedOrbit, newEpoch

def add_dsst_force_models(propagator, earth):

    cswl = CssiSpaceWeatherData("SpaceWeather-All-v1.2.txt")

    atmosphere = NRLMSISE00(cswl, CelestialBodyFactory.getSun(), earth)

    cross_section = 0.02
    cd = 2.2
    isotropic_drag = IsotropicDrag(cross_section, cd)

    drag_force = DragForce(atmosphere, isotropic_drag)

    gravityProvider = GravityFieldFactory.getUnnormalizedProvider(10, 10)
    propagator.addForceModel(DSSTZonal(gravityProvider))

    propagator.addForceModel(DSSTAtmosphericDrag(drag_force, gravityProvider.getMu()))

    propagator.addForceModel(DSSTTesseral(earth.getBodyFrame(), Constants.WGS84_EARTH_ANGULAR_VELOCITY, gravityProvider))

    rad = IsotropicRadiationSingleCoefficient(0.24, 1.0)
    propagator.addForceModel(DSSTSolarRadiationPressure(CelestialBodyFactory.getSun(), earth, rad, gravityProvider.getMu()))

    return propagator
    

def add_force_models(propagator, earth, ra, gravity=True, drag=True, solar=True, albedo=True):

    if gravity:
        gravityProvider = GravityFieldFactory.getNormalizedProvider(10, 10)
        propagator.addForceModel(HolmesFeatherstoneAttractionModel(FramesFactory.getITRF(IERSConventions.IERS_2010, True), gravityProvider))

    if drag:
        cswl = CssiSpaceWeatherData("SpaceWeather-All-v1.2.txt")

        atmosphere = NRLMSISE00(cswl, CelestialBodyFactory.getSun(), earth)

        cross_section = 0.02
        cd = 2.2
        isotropic_drag = IsotropicDrag(cross_section, cd)

        drag_force = DragForce(atmosphere, isotropic_drag)

        propagator.addForceModel(drag_force)
    
    if solar:

        rad = IsotropicRadiationSingleCoefficient(0.24, 1.0)
        sol_pressure_force = SolarRadiationPressure(CelestialBodyFactory.getSun(), earth , rad)
        sol_pressure_force.addOccultingBody(CelestialBodyFactory.getMoon(), Constants.MOON_EQUATORIAL_RADIUS)

        propagator.addForceModel(sol_pressure_force)

    if albedo:

        angularResolution = get_angle_rez(ra)
        print("Angular Resolution for ra: %d = %f rad-> %f" % (ra, angularResolution, math.degrees(angularResolution)))
        albedo_pressure_force = KnockeRediffusedForceModel(CelestialBodyFactory.getSun(), rad, earth.getEquatorialRadius(), angularResolution)
        
        propagator.addForceModel(albedo_pressure_force)
    
    return propagator

def set_up_prop(rp, ra, i, omega, raan, lv, epochDate, inertialFrame, ITRF, a=False, e=False, initialOrbit=False, DSST=True):
    if a == False:
        a = (rp + ra + 2 * Constants.WGS84_EARTH_EQUATORIAL_RADIUS) / 2.0    
    if e == False:
        e = 1.0 - (rp + Constants.WGS84_EARTH_EQUATORIAL_RADIUS) / a

    earth = OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS, 
                            Constants.WGS84_EARTH_FLATTENING, 
                            ITRF)

    ## Orbit construction as Keplerian
    if initialOrbit == False:
        initialOrbit = KeplerianOrbit(a, e, i, omega, raan, lv,
                                    PositionAngleType.TRUE,
                                    inertialFrame, epochDate, Constants.WGS84_EARTH_MU)
    satellite_mass = 16.0  # The models need a spacecraft mass, unit kg.
    initialState = SpacecraftState(initialOrbit, satellite_mass) 
    
    
    if DSST == False:
        minStep = 0.001
        maxstep = 1000.0
        initStep = 60.0
        positionTolerance = 1.0 #1.0 * 10 ** (-3)
        tolerances = NumericalPropagator.tolerances(positionTolerance, 
                                                    initialOrbit, 
                                                    initialOrbit.getType())

        integrator = DormandPrince853Integrator(minStep, maxstep, 
            JArray_double.cast_(tolerances[0]),  # Double array of doubles needs to be casted in Python
            JArray_double.cast_(tolerances[1]))
        integrator.setInitialStepSize(initStep)

        propagator = NumericalPropagator(integrator)
        propagator.setOrbitType(OrbitType.CARTESIAN)
        propagator.setInitialState(initialState)

        propagator = add_force_models(propagator, earth, a/2, albedo=False)
    else:
        minStep = initialState.getKeplerianPeriod()
        maxStep = 100. * minStep
        tolerances = DSSTPropagator.tolerances(1.0, initialOrbit)

        integrator = DormandPrince853Integrator(minStep, maxStep, 
            JArray_double.cast_(tolerances[0]),  # Double array of doubles needs to be casted in Python
            JArray_double.cast_(tolerances[1]))

        propagator = DSSTPropagator(integrator)
        propagator.setInitialState(initialState, PropagationType.MEAN)
        propagator = add_dsst_force_models(propagator, earth)


    return propagator