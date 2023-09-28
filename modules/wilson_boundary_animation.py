import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import seaborn as sns
import os

def load_data(file):
    f = open(file,'r')
    data = []
    beginRead = False
    beta = None
    for line in f:
        if line.startswith('beta'):
           beta = str(line.split()[1])
        if line.startswith('twist_coeff'):
           twist_c = str(line.split()[1])
        if not beginRead:
            if line.startswith('MEASURE start'):
                beginRead = True
        elif line.startswith('MEASURE end'):
            beginRead = False
        else:
            if line.startswith('plaquette:'):
                line = line.split(':')[1]
                try:
                    data.append([float(x) for x in line.split()])
                except:
                    data.append(line.split())

    data = pd.DataFrame(data[1:],columns=data[0])
    name = "Beta: " + beta + ', Twist coefficient: ' + twist_c
    return name,data

def create_figure(data,cols=4,window=50,title=""):
    fig, ax = plt.subplots()
    line, = ax.plot(list(data.columns.values), data[0:window].mean(), 'o--', label="Step: 0")
    #ax.set_ylim(data.to_numpy().min(),data.to_numpy().max())
    ax.set_ylim(data.to_numpy().min() , data.to_numpy().max())
    legend=ax.legend(loc=1)
    def animate(i):
        line.set_ydata(data[i:i+window].mean())  # update the data.
        legend.get_texts()[0].set_text(f"Step: {i}")
        return line,
    ani = animation.FuncAnimation(
        fig, animate,frames=len(data), interval=4, blit=True, save_count=50)
    fig.suptitle(title)
    plt.show()

if __name__ == "__main__":
    

    name,data = load_data(sys.argv[1])
    print(name)
    data = data.drop(columns=['sum'])
    columns = list(data.columns.values)
    columns = columns[-len(columns)//2:] + columns[:-len(columns)//2]
    print(columns)
    data = data.reindex(columns=columns)
    data = data[10:]
    print(data)
    create_figure(data,window=100,title=name)