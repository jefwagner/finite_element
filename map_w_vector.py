import numpy as np
import matplotlib.pyplot as plt

# Find the x and y coordinate for each element in w

def get_data(filename):
    with open(filename, 'r') as f:
        x = []
        y = []
        lines = f.readlines()
        points = lines[0].split()
        num_points = int(points[0])
        for line in lines[1:num_points+1]:
            words = line.split()
            x.append(float(words[0]))
            y.append(float(words[1]))
        return (num_points,x,y)

xy = get_data("w_vector.txt")
points = [xy[1],xy[2]]
points = np.asarray(points)
x = np.asarray(points[0])
y = np.asarray(points[1])

plt.figure()
# plt.gca().set_aspect('equal')
plt.plot(x,y)
plt.show()
plt.close()
