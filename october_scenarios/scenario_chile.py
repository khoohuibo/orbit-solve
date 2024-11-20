

from datetime import datetime
from datetime import timedelta


lines = []
    
with open('chiledeets.txt', 'r') as f:
    lines = f.read()
    f.close()



print([lines])
lines = lines.split("Observer:")
print(lines[1])
raise Exception
counter = 0
start_times = []
maximum_time_diff = timedelta()
for i in range(len(lines)):
    if "May" in lines[i]:
        dates = lines[i].split("   ")
        start_time = datetime.strptime(dates[0], "%d %b %Y %H:%M:%S.%f")
        start_times.append(start_time)

print(start_times)

for i in range(len(start_times)):
    if (i+1) >= len(start_times):
        break
    else:
        time_diff = start_times[i+1] - start_times[i]
        if time_diff > maximum_time_diff:
            maximum_time_diff = time_diff

print(maximum_time_diff)



        
