import matplotlib.pyplot as plt
from numba import vectorize
import numpy as np

@vectorize
def replace_zero(arr):
    if arr == 0:
        return np.nan
    else:
        return arr

arr = np.loadtxt("data/task5.txt")
a = np.log10(replace_zero(arr[:, 1]))

fig, ax = plt.subplots(2, 1, figsize=(8, 8), sharex=True)

ax[0].plot(arr[:, 0], arr[:, 1], 'r', label='transformations')
ax[0].set_title('transformations given N')
ax[0].set_ylabel('amount of transformations')
ax[0].legend(loc='lower right')

ax[1].plot(arr[:, 0], a, 'r', label='transformations')
ax[1].set_title('log10 of transformations given N')
ax[1].set_ylabel('amount of transformations log10')
ax[1].set_xlabel('matrix size N')
ax[1].legend(loc='lower right')

plt.subplots_adjust(bottom=0.2, hspace=0.2)
plt.savefig("data/transformations.pdf")
plt.show()