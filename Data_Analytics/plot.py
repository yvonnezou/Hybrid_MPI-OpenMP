import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('time_comp.txt')

data_mean = data.sum(axis=0)/505000

#data_max = np.max(data,0)

#data_min = np.min(data,0)

rank = np.arange(len(data_mean))

plt.bar(rank,data_mean)

plt.show()



