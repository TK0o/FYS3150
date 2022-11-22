import numpy as np
import matplotlib.pyplot as plt


def problem_7():
    parallel = np.loadtxt('data/time_p_1.txt')
    parallel_off = np.loadtxt('data/time_p_0.txt')
    
    fig, ax = plt.subplots(1, 1, figsize=(8, 4))
    ax.plot(parallel[:, 1], parallel[:, 0], c='r', label='4 threads')
    ax.plot(parallel[:, 1], parallel_off[:, 0], c='b', label='single thread')
    ax.set_xlabel('N MCMC cycles')
    ax.set_ylabel('time [s]')
    ax.legend()
    plt.title('Timing tests')
    plt.savefig('data/timing_parallel.pdf')
    plt.savefig('data/timing_parallel.png')

def methods():
    arr = np.loadtxt('data/test_methods.txt')
    fig, ax = plt.subplots(1, 1, figsize=(8, 4))
    N_cycle = np.zeros(len(arr))
    name = ['if test', 'modulo', 'vector index', 'int division']
    for i in range(len(arr)):
        N_cycle[i] = 10**(i + 1)

    for i in range(4):
        ax.plot(N_cycle, arr[:, i], label=name[i])
    
    ax.set_xlabel('N MCMC cycles')
    ax.set_ylabel('tot time [s]')
    plt.legend()
    plt.savefig('data/method_timing.pdf')
    plt.savefig('data/method_timing.png')
    


def problem_5():
    arr = np.loadtxt('data/problem_5.txt')
    fig, ax = plt.subplots(2, 2, figsize=(16, 8), sharex=True)
    
    cycle = arr[0]
    name = ['T = 1 J/k unordered', 'T = 2.4 J/k unordered', 
            'T = 1 J/k ordered', 'T = 2.4 J/k ordered']
    
    Color = ['r', 'g', 'b', 'orange']
    unit = ['[ J ]', '[<m>]']
    
    ax[0, 1].plot(cycle[:12], arr[4, :12], c=Color[0], label=name[0])
    ax[0, 0].plot(cycle[:12], abs(arr[5, :12]), c=Color[0], label=name[0])
    ax[1, 1].plot(cycle[:12], arr[14, :12], c=Color[2], label=name[2])
    ax[1, 0].plot(cycle[:12], abs(arr[15, :12]), c=Color[2], label=name[2])
    
    for i in range (2):
        for j in range (2):
            ax[i, j].legend()
    
    ax[0, 0].set_ylabel(unit[1])
    ax[1, 0].set_ylabel(unit[1])
    ax[0, 1].set_ylabel(unit[0])
    ax[1, 1].set_ylabel(unit[0])
    
    ax[0, 0].set_title('<|m|>')
    ax[0, 1].set_title('$<\epsilon>$')
    
    fig.subplots_adjust(hspace=0)
    plt.savefig('data/p5_t1.pdf')
    plt.savefig('data/p5_t1.png')
    plt.show()
    
    fig, ax = plt.subplots(2, 2, figsize=(16, 8), sharex=True)
    ax[0, 0].plot(cycle[:12], arr[4, :12], c=Color[1], label=name[1])
    ax[1, 0].plot(cycle[:12], abs(arr[5, :12]), c=Color[1], label=name[1])
    ax[0, 1].plot(cycle[:12], arr[14, :12], c=Color[3], label=name[3])
    ax[1, 1].plot(cycle[:12], abs(arr[15, :12]), c=Color[3], label=name[3])
    
    for i in range (2):
        for j in range (2):
            ax[i, j].legend()
        
    ax[0, 0].set_title('<|m|>')
    ax[0, 1].set_title('$<\epsilon>$')
    fig.subplots_adjust(hspace=0)
    plt.savefig('data/p5_t2.pdf')
    plt.savefig('data/p5_t2.png')
    plt.show()
    

def problem_8():
    for i in range(4):
        N = 40 + 20 * i
        arr = np.loadtxt(f'data/problem8_N{N}.txt')
        
        fig, ax = plt.subplots(2, 2, figsize=(18, 10))
        
        # ax[0, 0].scatter(arr[:, 0], arr[:, 1], s=25, marker='o', color='b')
        # ax[0, 0].set_xlabel("Temperature (T)", fontsize=20);
        # ax[0, 0].set_ylabel("Energy ", fontsize=20);
        # ax[0, 0].axis('tight');
        
        # ax[0, 1].scatter(arr[:, 0], abs(arr[:, 2]), s=25, marker='o', color='b')
        # ax[0, 1].set_xlabel("Temperature (T)", fontsize=20);
        # ax[0, 1].set_ylabel("Magnetization ", fontsize=20);
        # ax[0, 1].axis('tight');
        
        # ax[1, 0].scatter(arr[:, 0], arr[:, 3], s=25, marker='o', color='b')
        # ax[1, 0].set_xlabel("Temperature (T)", fontsize=20);  
        # ax[1, 0].set_ylabel("Specific Heat ", fontsize=20);   
        # ax[1, 0].axis('tight');
        
        # ax[1, 1].scatter(arr[:, 0], arr[:, 4], s=25, marker='o', color='b')
        # ax[1, 1].set_xlabel("Temperature (T)", fontsize=20); 
        # ax[1, 1].set_ylabel("Susceptibility", fontsize=20); 
        # ax[1, 1].axis('tight');
        # plt.suptitle(f'L = {N}', fontsize=20)
        # plt.savefig(f'data/phase_transition_N{N}.pdf')
        
        
        ax[0, 0].plot(arr[:, 0], arr[:, 1], color='b')
        ax[0, 0].set_xlabel("Temperature (T)", fontsize=20);
        ax[0, 0].set_ylabel("Energy ", fontsize=20);
        ax[0, 0].axis('tight');
        
        ax[0, 1].plot(arr[:, 0], abs(arr[:, 2]), color='b')
        ax[0, 1].set_xlabel("Temperature (T)", fontsize=20);
        ax[0, 1].set_ylabel("Magnetization ", fontsize=20);
        ax[0, 1].axis('tight');
        
        ax[1, 0].plot(arr[:, 0], arr[:, 3], color='b')
        ax[1, 0].set_xlabel("Temperature (T)", fontsize=20);  
        ax[1, 0].set_ylabel("Specific Heat ", fontsize=20);   
        ax[1, 0].axis('tight');
        
        ax[1, 1].plot(arr[:, 0], arr[:, 4], color='b')
        ax[1, 1].set_xlabel("Temperature (T)", fontsize=20); 
        ax[1, 1].set_ylabel("Susceptibility", fontsize=20); 
        ax[1, 1].axis('tight');
        
        plt.suptitle(f'phase transition L = {N}', fontsize=20)
        plt.savefig(f'data/phase_transition_N{N}.pdf')
        plt.savefig(f'data/phase_transition_N{N}.png')
        plt.show()
    
    
    
def problem_8_v2():    
     for i in range(4):
        N = 40 + 20 * i
        arr = np.loadtxt(f'data/problem8_v2_N{N}.txt')
        
        fig, ax = plt.subplots(2, 2, figsize=(18, 10))
        
        # ax[0, 0].scatter(arr[:, 0], arr[:, 1], s=25, marker='o', color='b')
        # ax[0, 0].set_xlabel("Temperature (T)", fontsize=20);
        # ax[0, 0].set_ylabel("Energy ", fontsize=20);
        # ax[0, 0].axis('tight');
        
        # ax[0, 1].scatter(arr[:, 0], abs(arr[:, 2]), s=25, marker='o', color='b')
        # ax[0, 1].set_xlabel("Temperature (T)", fontsize=20);
        # ax[0, 1].set_ylabel("Magnetization ", fontsize=20);
        # ax[0, 1].axis('tight');
        
        # ax[1, 0].scatter(arr[:, 0], arr[:, 3], s=25, marker='o', color='b')
        # ax[1, 0].set_xlabel("Temperature (T)", fontsize=20);  
        # ax[1, 0].set_ylabel("Specific Heat ", fontsize=20);   
        # ax[1, 0].axis('tight');
        
        # ax[1, 1].scatter(arr[:, 0], arr[:, 4], s=25, marker='o', color='b')
        # ax[1, 1].set_xlabel("Temperature (T)", fontsize=20); 
        # ax[1, 1].set_ylabel("Susceptibility", fontsize=20); 
        # ax[1, 1].axis('tight');
        # plt.suptitle(f'L = {N}', fontsize=20)
        # plt.savefig(f'data/phase_transition_N{N}.pdf')
        
        
        ax[0, 0].plot(arr[:, 0], arr[:, 1], color='b')
        ax[0, 0].set_xlabel("Temperature (T)", fontsize=20);
        ax[0, 0].set_ylabel("Energy ", fontsize=20);
        ax[0, 0].axis('tight');
        
        ax[0, 1].plot(arr[:, 0], abs(arr[:, 2]), color='b')
        ax[0, 1].set_xlabel("Temperature (T)", fontsize=20);
        ax[0, 1].set_ylabel("Magnetization ", fontsize=20);
        ax[0, 1].axis('tight');
        
        ax[1, 0].plot(arr[:, 0], arr[:, 3], color='b')
        ax[1, 0].set_xlabel("Temperature (T)", fontsize=20);  
        ax[1, 0].set_ylabel("Specific Heat ", fontsize=20);   
        ax[1, 0].axis('tight');
        
        ax[1, 1].plot(arr[:, 0], arr[:, 4], color='b')
        ax[1, 1].set_xlabel("Temperature (T)", fontsize=20); 
        ax[1, 1].set_ylabel("Susceptibility", fontsize=20); 
        ax[1, 1].axis('tight');
        
        plt.suptitle(f'phase transition L = {N}', fontsize=20)
        plt.savefig(f'data/phase_transition_v2_N{N}.pdf')
        plt.savefig(f'data/phase_transition_v2_N{N}.png')
        plt.show()
    
    
    
methods()
problem_7() 
problem_8_v2()
problem_8() 
problem_5()
    
    
    
    
    