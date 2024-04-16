import time
import matplotlib.pyplot as plt

current_timestamp = time.time()
current_mission_start_timestamp = current_timestamp + 10
current_duration = 65535

stop_signal = 0

second_list = []
current_duration_list = []

while (time.time() < current_mission_start_timestamp):
    0

counter = 0
while (stop_signal == 0):
    current_timestamp = time.time()

    if (current_timestamp - current_mission_start_timestamp > current_duration):
        stop_signal = 1
    else:
        current_duration = current_duration - (time.time() - current_mission_start_timestamp)

    print(current_duration)
    print(current_duration_list.append(current_duration))
    counter += 1
    print(counter)
    second_list.append(counter)
    time.sleep(1)

plt.plot(second_list, current_duration_list)
plt.show()