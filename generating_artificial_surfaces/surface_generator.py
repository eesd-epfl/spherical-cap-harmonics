"""
Created on Mon Mar  8 16:26:19 2021

@author: mahmoudshaqfa

Citing
Tamaas is the result of a science research project. To give proper credit to Tamaas and the researchers who have developed the numerical methods that it implements, please cite Tamaas as:

Frérot , L., Anciaux, G., Rey, V., Pham-Ba, S., & Molinari, J.-F. Tamaas: a library for elastic-plastic contact of periodic rough surfaces. Journal of Open Source Software, 5(51), 2121 (2020). doi:10.21105/joss.02121

If you use the elastic-plastic contact capabilities of Tamaas, please cite:

Frérot, L., Bonnet, M., Molinari, J.-F. & Anciaux, G. A Fourier-accelerated volume integral method for elastoplastic contact. Computer Methods in Applied Mechanics and Engineering 351, 951–976 (2019) doi:10.1016/j.cma.2019.04.006.

If you use the adhesive contact capabilities of Tamaas, please cite:

Rey, V., Anciaux, G. & Molinari, J.-F. Normal adhesive contact on rough surfaces: efficient algorithm for FFT-based BEM resolution. Comput Mech 1–13 (2017) doi:10.1007/s00466-017-1392-5.

"""

import numpy as np
import tamaas as tm
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import pyvista as pv

n = 128;
scale = 0.05; # scale the height of the aspirities

# Create spectrum object
spectrum = tm.Isopowerlaw2D()
# Set spectrum parameters
spectrum.q0 = 10
spectrum.q1 = 10
spectrum.q2 = 30
spectrum.hurst = 0.2
generator = tm.SurfaceGeneratorFilter2D()
generator.setSizes([n, n])
generator.setFilter(spectrum)

generator.random_seed = 0
surface = generator.buildSurface()

u = np.linspace(0, n-1, endpoint=True, num=n)
v = np.linspace(0, n-1, endpoint=True, num=n)

u, v = np.meshgrid(u, v)

surface *= scale;


fig = plt.figure()
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})


surf = ax.plot_surface(u, v, surface)
plt.show()

points = np.c_[u.reshape(-1), v.reshape(-1), surface.reshape(-1)]
cloud = pv.PolyData(points)
#cloud.plot(point_size=15)
surf = cloud.delaunay_2d()
surf.plot(cpos="xy", show_edges=True)
#surf.save("mesh.stl")
pv.save_meshio("rec_artificial_H_" + str(spectrum.hurst) + ".stl", surf)