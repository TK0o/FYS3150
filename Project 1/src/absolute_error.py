import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(1, 1, figsize=(9, 7))
C = ['r', 'b', 'g', 'darkorange', 'aqua', 'purple']

with open("data/absolute_error.txt", 'r') as data:
    c_index = 0 
    s = 0
    for i in data:
        if i.find ('n = ') != -1:
            I = i.split()
            data.readline()
            n = float(I[2])
            name = i
            
            int_n = int (n)
            arr = np.zeros((int_n - 1, 2), float)
            limit = int_n - 1
            s = -1
        
        elif s < limit:
            values = i.split()
            arr[s, 0] = values[0]
            arr[s, 1] = values[1]
        elif i.find('--') != -1:
            ax.plot(arr[:, 0], arr[:, 1], 
                    linewidth=0.75, label=name, color=C[c_index])
            c_index += 1
        s += 1

    ax.legend(loc = "lower center", ncol = 3)
    ax.set_title('Absolute error')
    ax.set_ylabel("$log_{10}$")
    ax.set_xlabel('x')
    
    plt.savefig('data/absolute_error.pdf')