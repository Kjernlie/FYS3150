import numpy as np
import matplotlib.pylab as plt
import csv

f_1 = []
f_2 = []

f = open("warmup.csv", "r+")
try:
    reader = csv.reader(f)
    for row in reader:
        f_1.append(float(row[0]))
        f_2.append(float(row[1]))

    f.close()

f_1 = np.array(f_1)
f_2 = np.array(f_2)
