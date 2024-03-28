from lib.helper_functions import beta_angle_test
import math
import matplotlib.pyplot as plt

beta_list = []
R_list = []
for x in range(360):
    R = math.radians(x)
    R_list.append(math.degrees(R))
    b = beta_angle_test(RA_sun=math.radians(x), raan=math.radians(228+x*10), inc=math.radians(97.7))
    beta_list.append(math.degrees(b))

print(len(R_list))
print(len(beta_list))
plt.plot(R_list, beta_list)
plt.show()