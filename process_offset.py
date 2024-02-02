import numpy as np
import matplotlib.pyplot as plt
import math
import csv

npzfile = np.load("output/raan_offset_0_sso.npz", allow_pickle=True)
raan_list = list(npzfile['arr_0'])
max_dist_list = list(npzfile['arr_1'])

fig = plt.figure()

ax2 = fig.add_subplot(222)
#ax2.plot(df['RAAN'], df['Max Dist'])
ax2.plot(raan_list, max_dist_list)
ax2.title.set_text('derived distance')


plt.legend()
plt.show()