import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
from matplotlib.patches import Polygon

data    = np.loadtxt('FLOWfieldGE.dat')
airfoil = np.loadtxt('GNUplot_tg_norm.dat')

pts = np.array([airfoil[:,0],airfoil[:,1]])
pts = np.transpose(pts)
p = Polygon(pts, closed=False)

ax = plt.gca()
ax.add_patch(p)

plt.plot(airfoil[:,0],airfoil[:,1])
plt.quiver(data[:,0], data[:,1], data[:,2], data[:,3])
plt.show()

x = data[:,0]
y = data[:,1]
z = data[:,5] 

# Set up a regular grid of interpolation points
xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
xi, yi = np.meshgrid(xi, yi)

# Interpolate
rbf = scipy.interpolate.Rbf(x, y, z, function='linear')
zi = rbf(xi, yi)

plt.imshow(zi, vmin=z.min(), vmax=z.max(), origin='lower',
           extent=[x.min(), x.max(), y.min(), y.max()])
plt.scatter(x, y, c=z)
plt.colorbar()
plt.show()
