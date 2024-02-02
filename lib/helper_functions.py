import math
import numpy as np
import datetime

def get_julian_datetime(date):
    """
    Convert a datetime object into julian float.
    Args:
        date: datetime-object of date in question

    Returns: float - Julian calculated datetime.
    Raises: 
        TypeError : Incorrect parameter type
        ValueError: Date out of range of equation
    """

    # Ensure correct format
    if not isinstance(date, datetime.datetime):
        raise TypeError('Invalid type for parameter "date" - expecting datetime')
    elif date.year < 1801 or date.year > 2099:
        raise ValueError('Datetime must be between year 1801 and 2099')

    # Perform the calculation
    julian_datetime = 367 * date.year - int((7 * (date.year + int((date.month + 9) / 12.0))) / 4.0) + int(
        (275 * date.month) / 9.0) + date.day + 1721013.5 + (
                          date.hour + date.minute / 60.0 + date.second / math.pow(60,
                                                                                  2)) / 24.0 - 0.5 * math.copysign(
        1, 100 * date.year + date.month - 190002.5) + 0.5

    return julian_datetime

def get_GMST(dat, utc=0, long=0):
    """
    Returns the siderial time in decimal hours. Longitude (long) is in 
    decimal degrees. If long=0, return value is Greenwich Mean Siderial Time 
    (GMST).
    """
    jd = get_julian_datetime(dat)

    t = (jd - 2451545.0)/36525
    # Greenwich siderial time at 0h UTC (hours)
    st = (24110.54841 + 8640184.812866 * t +
          0.093104 * t**2 - 0.0000062 * t**3) / 3600

    # Greenwich siderial time at given UTC
    st = st + 1.00273790935*utc

    # Local siderial time at given UTC (longitude in degrees)
    st = st + long/15
    st = st % 24
    st = st / 24
    st = st * 2*math.pi

    return(st)

def mean_motion(mass, altitude):
    grav = 6.67384 * 10 ** -11
    earth_mass = 5.972* 10 ** 24
    a = (6378 + altitude) * 1000
    return math.sqrt(grav * (earth_mass + mass)/ (a ** 3))


def delta_t(veloc):
    Vx = veloc/1000
    We = 7.272 * 10 ** -5
    Re = 6378

    i = math.radians(45.4)

    A = Vx
    B = We*Re

    comp_under_1 = (2 * A **2)
    comp_under_2 = (4*A*B*math.cos(i))
    comp_under_3 = (B * math.cos(i)) ** 2
    comp_under_4 = (B**2)

    comp_under_combined = (comp_under_1 - comp_under_2 + comp_under_3 + comp_under_4)

    delta_t = math.sqrt(1000**2/(comp_under_combined))


    delta_y = math.sqrt((Vx **2 + (We*Re)**2) * delta_t ** 2 - (2*Vx*We*Re*math.cos(i))*delta_t ** 2)
    delta_x = -Vx * delta_t + We*Re*math.cos(i)*delta_t

    delta_d = math.sqrt(delta_y**2 + delta_x**2)

    delta_RAAN = We * delta_t 

    delta_mean_motion = mean_motion(16, 605) * delta_t

    return delta_t

def alternate_delta_t(veloc, inclination):
    Vx = veloc/1000
    We = 7.272 * 10 ** -5
    Re = 6378

    i = inclination

    A = Vx
    B = We*Re

    comp_under_combined = (A ** 2) + (B ** 2) - (2 * A * B * math.cos(i))

    delta_t = math.sqrt(1000**2/(comp_under_combined))

    return delta_t

def r3(angle):
    return(np.array([[math.cos(angle), -math.sin(angle), 0],
            [math.sin(angle), math.cos(angle), 0],
            [0,0,1]]))

def get_abs_vel(vel_tmp):
    return math.sqrt((vel_tmp.getX())**2 + (vel_tmp.getY())**2 + (vel_tmp.getZ())**2)

def distance_between_two(p, p2):
    return math.sqrt((p2.getX() - p.getX())**2 + (p2.getY() - p.getY())**2 + (p2.getZ() - p.getZ())**2)

def get_angle_rez(ra):
    return (math.pi/2) - math.asin(6378/ra)

