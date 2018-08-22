# This file is used to bar charts of the input data.
# In the input file, rows are the iteration number, columns are the core number and the value is the rutime. 
# Users need to store it in the same path as the input data's.
# The x-axis is the core number.
# The y-axis is the runtime and the unit is second. The programme will extract the maximum, minimum and average value in rows.

import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('out1.txt')
data_5000 = data[5000:9999,:]
data_mean = data_5000.sum(axis=0)/5000
data_max = np.max(data_5000,0)
data_min = np.min(data_5000,0)

rank = np.arange(len(data_mean))

fig = plt.figure()

l3 = plt.bar(rank,data_max,color='b')
l2 = plt.bar(rank,data_mean,color='g')
l1 = plt.bar(rank,data_min,color='r')
plt.rcParams.update({'font.size':22})

plt.xlabel('Thread Number')
plt.ylabel('Time (s)')
plt.title('Computation time in each iteration on every thread (Cirrus and multiple version)')
plt.legend(handles = [l1, l2, l3], labels = ['min', 'average','max'], loc = 'best')

plt.show()



