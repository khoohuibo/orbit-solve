import numpy as np
import matplotlib.pyplot as plt
import math
import csv


date = []
dist = []
with open("output_sso/beta_three_sso_july_2025_false_delayed.csv", newline='') as f:
    reader = csv.reader(f, delimiter=' ', quotechar="|")
    for row in reader:
        if "date" in row[0]:
            continue
        divided = row[1].split(",")
        #print(row)
        try:
            if float(divided[1]) <= 1000:
                date.append(row[0])
                dist.append(divided[1])
        except ValueError as e:
            print(e)
for i in range(len(date)):
    print(date[i], dist[i])