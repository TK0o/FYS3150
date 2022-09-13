import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 1, figsize=(6, 4))
with open("data/relative_error.txt", 'r') as data:
    data.readline()
    n = []
    y = []
    for i in data:
        I = i.split()        
        n.append(float(I[0]))
        y.append(float(I[1]))
        
    ax.plot(n, y, 'r', label='relative error')
    ax.set_title('relative error')
    ax.legend(loc="lower right")
    ax.set_ylabel("$log_{10}$")
    ax.set_xlabel('n')
    
    plt.savefig('data/relative_error.pdf')