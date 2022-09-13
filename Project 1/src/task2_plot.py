import numpy as np
import matplotlib.pyplot as plt

x = []
u = []
filename = 'task2.txt'
with open(filename, 'r') as file:
    file.readline()
    for i in file:
        I = i.split()
        x.append(float(I[0]))
        u.append(float(I[1]))


x_val = np.array(x)
u_val = np.array(u)
    
# print (x, u)



def analytical(x):
    return 1 - (1 - np.exp(-10)) * x - np.exp(-10 * x)

X = np.linspace(0, 1, 1000)
U = analytical (X)
plt.plot(X, U, 'k', linewidth=0.8)
plt.plot(x, u, '.', color='r')