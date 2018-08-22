# This file is used to plot 3D images of the input data.
# In the input file, rows are the iteration number, columns are the core number and the value is the rutime. 
# Users need to store it in the same path as the input data's
# Due to the runtime, it only processes the first 5000 rows in the data file.
# The x-axis is the core number.
# The y-axis is the iteration number.
# The z-axis is the runtime and the unit is second. 

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

data = np.loadtxt('time_iter.txt')
data_5000 = data[0:4999,:]
nx = len(data_5000[0])
ny = len(data_5000)
x = range(nx)
y = range(ny)
hf = plt.figure()
ha = hf.add_subplot(111, projection='3d')
 
X, Y = np.meshgrid(x, y)

ha.plot_surface(X, Y, data_5000,cmap=plt.cm.coolwarm)

plt.show()





