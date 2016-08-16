import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Find the x and y coordinate for each element in w

def get_data(filename):
    with open(filename, 'r') as f:
        x = []
        y = []
        z = []
        lines = f.readlines()
        points = lines[0].split()
        num_points = int(points[0])
        for line in lines[1:num_points+1]:
            words = line.split()
            x.append(float(words[0]))
            y.append(float(words[1]))
            z.append(float(words[2]))
        return (num_points,x,y,z)

# Find the nodes in each triangle. This relies
#   on the fact that region_mesh is structured.

def get_tri(filename, num_points):
    with open(filename, 'r') as f:
        i = []
        j = []
        k = []
        lines = f.readlines()
        tris = lines[num_points+1].split()
        num_tris = int(tris[0])
        for line in lines[num_points+2:(num_points+2)+(num_tris)]:
            words = line.split()
            i.append(float(words[0]))
            j.append(float(words[1]))
            k.append(float(words[2]))
        return (i,j,k)

xy = get_data("data/w_vector_d_2.100_rho_0.000_gamma_1.000.dat")
points = [xy[1],xy[2],xy[3]]
points = np.asarray(points)
x = np.asarray(points[0])
y = np.asarray(points[1])
z = np.asarray(points[2])

tris = get_tri("data/w_vector_d_2.100_rho_0.000_gamma_1.000.dat", xy[0])
tris = np.asarray(tris).T

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_trisurf(x,y,z, triangles=tris)
plt.show()

xy = get_data("data/w_vector_d_10.000_rho_0.000_gamma_1.000.dat")
points = [xy[1],xy[2],xy[3]]
points = np.asarray(points)
x = np.asarray(points[0])
y = np.asarray(points[1])
z = np.asarray(points[2])

tris = get_tri("data/w_vector_d_10.000_rho_0.000_gamma_1.000.dat", xy[0])
tris = np.asarray(tris).T

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_trisurf(x,y,z, triangles=tris)
plt.show()

xy = get_data("data/w_vector_d_5.510_rho_0.000_gamma_1.000.dat")
points = [xy[1],xy[2],xy[3]]
points = np.asarray(points)
x = np.asarray(points[0])
y = np.asarray(points[1])
z = np.asarray(points[2])

tris = get_tri("data/w_vector_d_5.510_rho_0.000_gamma_1.000.dat", xy[0])
tris = np.asarray(tris).T

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_trisurf(x,y,z, triangles=tris)
plt.show()
