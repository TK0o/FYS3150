import matplotlib.pyplot as plt
import numpy as np

x = []
u = []
filename = 'data/task2.txt'
with open(filename, 'r') as file:
    file.readline()
    for i in file:
        I = i.split()
        x.append(float(I[0]))
        u.append(float(I[1]))
        
x_val = np.array(x)
u_val = np.array(u)

plt.plot(x_val, u_val, color='r', label='$u (x)$')
plt.title('Analytical solution u(x)')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.savefig('data/task2_plot.pdf')