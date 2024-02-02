import numpy as np
import matplotlib.pyplot as plt
import math
import csv


npzfile = np.load("output/raan_all_-20_20.npz", allow_pickle=True)
npzfile_new = np.load("output/raan_all_-10_10.npz", allow_pickle=True)
npzfile_1 = np.load("output/raan_all_20_50.npz", allow_pickle=True)
npzfile_2 = np.load("output/raan_all_50_100.npz", allow_pickle=True)

raan_list = [math.degrees(x) for x in list(npzfile['arr_0'])]
raan_list_new = [math.degrees(x) for x in list(npzfile_new['arr_0'])]
raan_list_1 = [math.degrees(x) for x in list(npzfile_1['arr_0'])]
raan_list_2 = [math.degrees(x) for x in list(npzfile_2['arr_0'])]
raan_list = raan_list + raan_list_new + raan_list_1 + raan_list_2
#raan_list = list(npzfile['arr_0'])
max_dist_list = list(npzfile['arr_1'])
max_dist_list_new = list(npzfile_new['arr_1'])
max_dist_list_1 = list(npzfile_1['arr_1'])
max_dist_list_2 = list(npzfile_2['arr_1'])
max_dist_list = max_dist_list + max_dist_list_new + max_dist_list_1 + max_dist_list_2


print(raan_list)
print(max_dist_list)

import pandas as pd


data = {'RAAN': raan_list, 'Max Dist': max_dist_list}
df = pd.DataFrame(data).sort_values(by=['RAAN'])
print(df)

df.to_csv('out.csv', index=True)

"""
npzfile = np.load("output/raan_offset_0.npz", allow_pickle=True)
raan_list = list(npzfile['arr_0'])
max_dist_list = list(npzfile['arr_1'])
"""
fig = plt.figure()

ax2 = fig.add_subplot(222)
ax2.plot(df['RAAN'], df['Max Dist'])
#ax2.plot(raan_list, max_dist_list)
ax2.title.set_text('derived distance')


plt.legend()
plt.show()