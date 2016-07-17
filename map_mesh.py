import numpy as np
import matplotlib.pyplot as plt

def get_xy(filename):
    with open(filename, 'r') as f:
        x = []
        y = []
        lines = f.readlines()
        points = lines[0].split()
        num_points = int(points[1])
        for line in lines[1:num_points+1]:
            words = line.split()
            x.append(float(words[0]))
            y.append(float(words[1]))
        return (num_points,x,y)

def get_tri(filename, num_points):
    with open(filename, 'r') as f:
        i = []
        j = []
        k = []
        lines = f.readlines()
        tris = lines[num_points+2].split()
        num_tris = int(tris[1])
        for line in lines[num_points+3:(num_points+3)+(num_tris)]:
            words = line.split()
            i.append(float(words[0]))
            j.append(float(words[1]))
            k.append(float(words[2]))
        return (i,j,k)

xy = get_xy("region_mesh.txt")
points = [xy[1],xy[2]]
points = np.asarray(points)
x = np.asarray(points[0])
y = np.asarray(points[1])
tris = get_tri("region_mesh.txt", xy[0])
tris = np.asarray(tris).T

plt.figure()
plt.gca().set_aspect('equal')
plt.triplot(x,y,tris)
plt.show()
plt.close()
