import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(2, 3, figsize=(15, 10))

plt.rc('legend', fontsize=12)
plt.rc('axes', labelsize=15)
plt.rc('axes', titlesize=14)
plt.rc('font', size=13)

with open("data/task7_alter_b.txt", 'r') as data:
    p_index_x = 0
    p_index_y = 0
    s = 0
    
    for i in data:
        if i.find ('n + 1 = ') != -1:
            I = i.split()
            data.readline()
            n = float(I[4])
            int_n = int (n)
            name = i
            
            arr = np.zeros((int_n, 3), float)
            dot_index = int((n - 1) / 10)
            limit = int_n 
            s = -1           
            
        elif s < limit:
            values = i.split()
            arr[s, 0] = values[0]
            arr[s, 1] = values[1]
            arr[s, 2] = values[2]

        elif i.find('--') != -1:
            ax[p_index_y, p_index_x].plot(
                arr[:, 0], arr[:, 1], color='blue', 
                linewidth=0.75, label='Analytical')
                                    
            ax[p_index_y, p_index_x].plot(
                arr[:, 0], arr[:, 2], 
                color='k', linestyle='dashed', 
                linewidth=0.75, label='Numeric')
            
            ax[p_index_y, p_index_x].plot(
                arr[::dot_index, 0], arr[::dot_index, 2], 
                '.', color='r', label='Numeric')
            
            ax[p_index_y, p_index_x].set_title(f"n = {int_n - 1}")
            ax[p_index_y, p_index_x].set_xlabel("x")
            ax[p_index_y, p_index_x].set_ylabel("y")
            ax[p_index_y, p_index_x].legend(loc='upper right')
            
            p_index_x += 1 
            if p_index_x > 2:
                p_index_y += 1
                p_index_x = 0
        s += 1

    plt.suptitle('Comparison between Numeric '
                 'and Analytic, Specific Algorithm', fontsize=20)
                 
    plt.subplots_adjust(bottom=0.1, top=0.9, hspace=0.25, wspace=0.3)
    plt.savefig('data/alter_algorithm_comparison.pdf')
