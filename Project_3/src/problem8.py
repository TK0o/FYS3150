import numpy as np
import matplotlib.pyplot as plt

loc_str = "data/problem_8/"
file = lambda filename: np.loadtxt(loc_str + filename)

def problem_8_single_particle(filename):
    arr = file(filename)
    x, y, z = arr[:, 0], arr[:, 1], arr[:, 2]
    ax = plt.axes(projection ='3d')
    ax.plot(x, y, z)
    ax.set_xlabel("x $[\mu m]$")
    ax.set_ylabel("y $[\mu m]$")
    ax.set_zlabel("z $[\mu m]$")
    ax.set_title("rk4 single particle")
    plt.savefig('data/images/single_particle_3d.png')
    plt.show()

def problem_8_double_particle(filename):
    
    C = ['r', 'b']
    name = ['$p_{1}$', '$p_{2}$']

    fig, ax = plt.subplots(1, 2, figsize=(8, 4))
    for j in range (2):
        if j == 0:
            arr = file(filename)
            ax[j].set_title("rk4 double particle no interaction")
        else:
            arr = file(filename.replace('0', '1'))
            ax[j].set_title("rk4 double particle with interaction")
        
        for i in range (2):
            x, y = arr[:, i * 3], arr[:, 1 + i * 3]
            ax[j].plot(x, y, color=C[i], label=name[i])
            ax[j].set_xlabel("x $[\mu m]$")
            ax[j].set_ylabel("y $[\mu m]$")
            ax[j].legend()
            
    fig.tight_layout()
    plt.savefig('data/images/double_particle_comparison.png')
    plt.show()

def problem_8_phase(filename):
    C = ['r', 'b']
    name = ['$p_{1}$', '$p_{2}$']

    fig, ax = plt.subplots(2, 2, figsize=(10, 10))
    for j in range (2):
        for i in range (2):
            if j == 0:
                arr = file(filename)
                title_str1 = "rk4 2 particle no interaction "
                ax[0, j].set_title(title_str1 + "$(x, \ v_x)$")
                ax[1, j].set_title(title_str1 + "$(z, \ v_z)$")
            else:
                arr = file(filename.replace('0', '1'))
                title_str2 = "rk4 2 particle with interaction "
                ax[0, j].set_title(title_str2 + "$(x, \ v_x)$")
                ax[1, j].set_title(title_str2 + "$(z, \ v_z)$")
           
            x1, y1 = arr[:, i * 2], arr[:, 1 + i * 2]
            x2, y2 = arr[:, i * 2 + 4], arr[:, 5 + i * 2]
            ax[i, j].plot(x1, y1, color=C[0], label=name[0])
            ax[i, j].plot(x2, y2, color=C[1], label=name[1])
            if i == 0:
                ax[i, j].set_xlabel("x $[\mu m]$")
                ax[i, j].set_ylabel("$v_x \ [\mu m /\mu s]$")
            elif i == 1:
                ax[i, j].set_xlabel("z $[\mu m]$")
                ax[i, j].set_ylabel("$v_z \ [\mu m /\mu s]$")
            ax[i, j].legend(bbox_to_anchor=(0.475, -0.25), 
                            loc='lower center', ncol=2)
           
    fig.tight_layout()
    plt.savefig('data/images/phase.png')
    plt.show()
    
def problem_8_double_interact(filename, interact=0):
    C = ['r', 'b']
    name = ['$p_{1}$', '$p_{2}$']

    arr = file(filename)
    ax = plt.axes(projection ='3d')
    if interact == 0:
        ax.set_title("rk4 double particle with no interaction")
    else:
        ax.set_title("rk4 double particle with interaction")
        
    for i in range (2):
        x, y, z = arr[:, i * 2], arr[:, i * 2 + 1], arr[:, i * 2 + 2]
        
        ax.plot(x, y, z, label=name[i], color=C[i])
        ax.set_xlabel("x $[\mu m]$")
        ax.set_ylabel("y $[\mu m]$")
        ax.set_zlabel("z $[\mu m]$")
            
    ax.legend(bbox_to_anchor=(0.85, -0.05) ,loc='upper right', ncol=2)
    if interact == 0:
        plt.savefig('data/images/double_interaction_0_3d.png')
    else:
        plt.savefig('data/images/double_interaction_1_3d.png')
    
    plt.show()


def problem_8_error():
    filename1 = "euler_p8_relative_error{}_.txt"
    filename2 = "rk4_p8_relative_error{}_.txt"
    
    C = ['r', 'b']
    name = ['euler', 'rk4']
    s = 0
    fig, ax = plt.subplots(2, 2, figsize=(10, 10))
    for i in range (2):
        for j in range (2):
            n = int (4000 * 2**s)
            t = np.arange(n) + 1
            t = t * (50 / n)
            arr1 = file(filename1.format(n))
            arr2 = file(filename2.format(n))
            
            arr_num_euler = arr1[:, 0:3]
            arr_num_rk4 = arr2[:, 0:3]
            arr_analytical = arr1[:, 3:]
            
            arr_err_euler = abs(arr_analytical - arr_num_euler)
            arr_err_rk4 = abs(arr_analytical - arr_num_rk4)
            
            arr_euler = np.zeros(n)
            arr_rk4 = np.zeros(n)

            for q in range (n):
                arr_euler[q] = np.linalg.norm(arr_err_euler[q])
                arr_rk4[q] = np.linalg.norm(arr_err_rk4[q])
            ax[i, j].set_title(f'relative error for n = {n}')
            ax[i, j].plot(t, arr_euler, ls='--', lw=0.75, 
                          label=name[0], color=C[0])
            ax[i, j].plot(t, arr_rk4, ls='--', lw=0.75, 
                          label=name[1], color=C[1])
            ax[i, j].legend(loc='upper left')
            ax[i, j].set_xlabel('t $[\mu s]$')
 
            s += 1
    plt.savefig('data/images/error.png')
    plt.show()

def convergence_rate():
    filename1 = "euler_p8_relative_error{}_.txt"
    filename2 = "rk4_p8_relative_error{}_.txt"
    
    h_k = np.zeros(4)
    max_rk4 = np.zeros(4)
    max_euler = np.zeros(4)

    for i in range (4):
        n = int (4000 * 2**i)

        arr1 = file(filename1.format(n))
        arr2 = file(filename2.format(n))
        
        arr_num_euler = arr1[:, 0:3]
        arr_num_rk4 = arr2[:, 0:3]
        arr_analytical = arr1[:, 3:]
        
        arr_err_euler = abs(arr_analytical - arr_num_euler)
        arr_err_rk4 = abs(arr_analytical - arr_num_rk4)
        
        h_k[i] = n
        max_euler[i] = np.max(arr_err_euler)
        max_rk4[i] = np.max(arr_err_rk4)

    def error_sum(arr, h, i):
        return np.log10(arr[i + 1] / arr[i]) / np.log10(h[i + 1] / h[i])

    S_rk4 = 0
    S_euler = 0
    for i in range(3):
        S_rk4 += error_sum (max_rk4, h_k, i)
        S_euler += error_sum (max_euler, h_k, i)
        
    r_err_rk4 = S_rk4 / 3  
    r_err_euler = S_euler / 3
    print (f'convergence rate r_err for rk4 is r_err_rk4 = {r_err_rk4}')
    print (f'convergence rate r_err for auler is r_err_euler = {r_err_euler}')

problem_8_single_particle("rk4_p8_1_particle.txt")
problem_8_double_particle("rk4_p8_2_particle_interact_0_.txt")
problem_8_double_interact("rk4_p8_2_particle_interact_0_.txt")
problem_8_double_interact("rk4_p8_2_particle_interact_1_.txt", 1)
problem_8_phase("rk4_p8_phase_interact_0_.txt")
problem_8_double_interact("rk4_p8_2_particle_interact_1_.txt")
problem_8_error()
# convergence_rate()



















