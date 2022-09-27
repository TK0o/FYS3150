import matplotlib.pyplot as plt
import numpy as np

x_axis_line = lambda x: 0 * x;
filename = "data/vectors_" + str(10) + ".txt"

for i in range (1, 3):
    filename = "data/vectors_" + str(10**i) + ".txt"    
    arr = np.loadtxt(filename)
    
    x = arr[:, 0]
    arr = np.delete(arr, 0, 1)
    analytic = arr[:, :3]
    jacobian = arr[:, 3:]
    
    fig, ax = plt.subplots(3, 1, figsize=(9,8), sharex=True)
    for i in range (3):
        if i < 1:
            ax[i].plot(x, analytic[:, i], color='r', 
                       linestyle='-.', label='analytic')
            ax[i].plot(x, jacobian[:, i], color='b', 
                       linestyle='--', label='jacobian')
        else:
            ax[i].plot(x, analytic[:, i], color='r', linestyle='-.')
            ax[i].plot(x, jacobian[:, i], color='b', linestyle='--')
        ax[i].set_ylabel('$V_{}$'.format(i + 1))
    
    fig.legend(bbox_to_anchor=(0.65, 0.16), ncol = 2)
    plt.subplots_adjust(bottom=0.2, top=0.94, wspace=1.5, hspace=0.2)
    fig.suptitle('corresponding eigenvectors for the 3 ' +
                 'smallest eigenvalue with n = {}'.format(len(x) - 1))
    plt.savefig("data/eigenvectors_n_{}.pdf".format(len(x) - 1))
    plt.show()

# , loc='lower center'
# v_1, v_2, v_3 = arr[:, 1], arr[:, 2], arr[:, 3]
# ana_v_1, ana_v_2, ana_v_3 = arr[:, 4], arr[:, 5], arr[:, 6]

# n = len(x) - 1
# v_list = np.array([v_1, v_2, v_3])
# ana_v_list = np.array([ana_v_1, ana_v_2, ana_v_3])
# v_colors = np.array(['black', 'royalblue', 'red'])

# fig, ax = plt.subplots(3, 1, figsize=(9,6), sharex=True)
# for i in range(len(v_list)):
#     ax[i].plot(x, np.zeros(len(x)), linestyle='dashed', color='dimgrey', lw=.8)
#     ax[i].plot(x, ana_v_list[i], color=v_colors[i], linestyle='dashed', label='Analytical')
#     ax[i].plot(x, v_list[i], color=v_colors[i], label='Numerical')
#     ax[i].set_ylabel(f'$v_{int(i+1)}$')
#     ax[i].legend()





