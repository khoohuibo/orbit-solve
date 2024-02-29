
import math
from datetime import datetime
import pytz
from lib.helper_functions import get_julian_datetime

newDate = datetime(1986, 3, 10, 0, 0, 0, tzinfo=pytz.utc)

jd = get_julian_datetime(newDate)

D = jd - 2451545.0

g = 357.529 + 0.98560028 * D
q = 280.459 + 0.98564736 * D
L = math.radians(q + 1.915 * math.sin(math.radians(g)) + 0.020 * math.sin(math.radians(2 * g)))

R = 1.00015 - 0.01671 * math.cos(g) - 0.00014 * math.cos(2 * g)
e = math.radians(23.439 - 0.00000036 * D)
RA = math.atan((math.cos(e) * math.sin(L))/math.cos(L))
d = math.asin(math.sin(e) * math.sin(L))
print(jd, D, g, q, L, math.degrees(L), R, e)
print(360+math.degrees(RA), 360+math.degrees(d))
print(RA, d)

g = -4614.824413
q = -4692.131
L = -4690.377707
e = 23.4408162
