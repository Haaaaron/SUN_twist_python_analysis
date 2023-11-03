import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import seaborn as sns
import os

def create_figure(data,title="",name=""):
    fig, ax = plt.subplots()
    #data = data.set_index("Unnamed: 0")
    betas = list(data.index)
    x = list(data.columns.values)
    data = data.to_numpy()
    line, = ax.plot(x, data[0], 'o--', label=f"Step: 0, Beta: {betas[0]}")
    #ax.set_ylim(data.to_numpy().min(),data.to_numpy().max())
    ax.set_ylim(data.min() , data.max())
    ax.set_xlabel("z-index")
    ax.set_ylabel(r"Wilson action average")
    legend=ax.legend(loc=1)
    def animate(i):
        line.set_ydata(data[i])  # update the data.
        legend.get_texts()[0].set_text(f"Step: {i}, Beta: {betas[i]}")
        return line,
    ani = animation.FuncAnimation(
        fig, animate,frames=len(data), interval=20, blit=True, save_count=1000)
    fig.suptitle(title)
    plt.show()
    writervideo = animation.FFMpegWriter(fps=6)
    ani.save(f'{name}.mp4', writer=writervideo)
    plt.close()

if __name__ == "__main__":
    
    data = pd.read_csv(sys.argv[1],index_col=0)
    create_figure(data,title=sys.argv[2],name=sys.argv[3])