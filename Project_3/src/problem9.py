import numpy as np
import matplotlib.pyplot as plt

loc_str = "data/problem_9/"
file = lambda filename: np.loadtxt(loc_str + filename)

def p9_general():
    C = ['r', 'b', 'g']
    arr = file('p9_time_dependent.txt')
    for i in range(1, 4):
        plt.plot(arr[:, 0], arr[:, i], label=f'f = {0.1 + 0.3 * (i - 1)}', 
                 color=C[i - 1])
    
    plt.ylabel('n particles')
    plt.xlabel('$W_V \ [MHz]$')
    plt.title('n particles in the trap depending on $W_V$')
    plt.legend(bbox_to_anchor=(0.475, -0.275), loc='lower center', ncol=3)
    plt.savefig('data/images/fraction_of_particle.png')
    plt.show()

def p9_fine_search():
    arr = file("p9_fine_search_interaction_.txt")
    
    C = ['r', 'b', 'g']
    name = ['no_interaction', 'interaction']
    plt.title('specific search for n particles in the trap depending on $W_V$')
    for i in range(1, 3):
        plt.plot(arr[:, 0], arr[:, i], lw=0.9,
                 label=name[i - 1], color=C[i - 1])           
    plt.legend()
    plt.ylabel('n particles')
    plt.xlabel('$W_V \ [MHz]$')
    plt.savefig('data/images/fraction_of_particle_fine_search.png')
    plt.show()

def p9_fine_search_last_attempt():
    arr = file("p9_fine_search_ 2_interaction_.txt")
    
    C = ['r', 'b', 'g']
    name = ['no_interaction', 'interaction']
    plt.title('specific search for n particles in the trap depending on $W_V$')
    for i in range(1, 3):
        plt.plot(arr[:, 0], arr[:, i], lw=0.9,
                 label=name[i - 1], color=C[i - 1])           
    plt.legend()
    plt.ylabel('n particles')
    plt.xlabel('$W_V \ [MHz]$')
    plt.savefig('data/images/fraction_finest_search.png')
    plt.show()

p9_fine_search_last_attempt()
p9_general()
p9_fine_search()
