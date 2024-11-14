import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import os
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation 
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
FIG_SIZE=(10,20)
plt.rcParams.update({'font.size': 11})

def create_figure_polar(datas,cols=4, mean = True,both=False, title=""):
    max_radial = 0
    #Looking for maximal polar radius to have even plot
    for k,(name,data) in enumerate(datas.items()):
        

        data = data['sum']
        max_radial = max(max_radial,np.max(abs(data)))

    total = len(datas)
    rows = total // cols
    if total % cols != 0:
        rows += 1
        
    position = range(1,total+1)
    fig = plt.figure(1,figsize=FIG_SIZE)
    for k,(name,data) in enumerate(datas.items()):
        print(k)
        ax = fig.add_subplot(rows,cols,position[k],projection="polar")

        #delta = data.to_numpy().max()-data.to_numpy().min()
        delta = 0

        if mean:
            data = data['sum']
            angular = np.angle(data.mean())
            radial = abs(data.mean())
            ax.plot(angular,radial,'o')
        if both:
            data = data['sum']
            angular_m = np.angle(data.mean())
            radial_m = abs(data.mean())
            angular = np.angle(data)
            radial = abs(data)
            if k == 1: ax.plot(angular[0],radial[0],'b+',alpha=0.5,markersize=5,label=r"$P_i$")
            for i,(ang,rad) in enumerate(zip(angular[1:2000],radial[1:2000])):
                ax.plot(ang,rad,'b+',alpha=0.5,markersize=5)
            if k==1:
                ax.plot(angular_m,radial_m,'ro',label=r"$\langle P_i \rangle$")
            else:
                ax.plot(angular_m,radial_m,'ro')
        #ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5),fontsize="7")
        else:
            data = data['sum']
            angular = np.angle(data)
            radial = abs(data)
            for i,(ang,rad) in enumerate(zip(angular,radial)):
                ax.plot(ang,rad,'o',label=f'{i}')
        #ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5),fontsize="7")
        ax.set_title(r"$\beta=${}".format(name.split()[0]))
        ax.set_rmax(max_radial+max_radial/10)
    fig.legend("lower right")

    fig.suptitle(title)
    print("something")
    #plt.xticks([float(x) for x in datas[0].head()[1:-1]])
    fig.tight_layout(pad=0.4, w_pad=1, h_pad=1.0)
    #plt.show()
    plt.savefig("./output_polyakov.svg")

def create_figure_polar_zindex(datas,cols=4, mean = True,both=False, title="",skip=0,stride=1):
    max_radial = 0
    #Looking for maximal polar radius to have even plot
    for k,(name,data) in enumerate(datas.items()):
        
        #print(data["sum"])
        data = data.drop(columns=['sum'])
        max_radial = max(max_radial,np.max(abs(data)))

        total = len(data.keys())
    rows = total // cols
    if total % cols != 0:
        rows += 1
         
    position = range(1,total+1)
    fig = plt.figure(1,figsize=FIG_SIZE)
    for k,(name,data) in enumerate(datas.items()):
        data = data.drop(columns=['sum'])

        #delta = data.to_numpy().max()-data.to_numpy().min()
        delta = 0
        colors = []
        n_points = len(data)
        cmap = mpl.cm.get_cmap("viridis",n_points)
        c = np.arange(1., n_points + 1)
        norm = mpl.colors.BoundaryNorm(np.arange(len(c)+1)+0.5,len(c))
        for d in data:
            print(d)
            ax = fig.add_subplot(rows,cols,position[int(d)//stride],projection="polar")
            num_z=len(data.keys())
            angular_m = np.angle(data[d].mean())
            radial_m = abs(data[d].mean())
            angular = np.angle(data[d])
            radial = abs(data[d])
            
            if k == 1: ax.plot(angular[0],radial[0],'+',alpha=0.3,markersize=5)
            for i,(ang,rad) in enumerate(zip(angular[1:],radial[1:])):
                # if i == 0:
                #     #ax.plot(ang,rad,'+',color=((i+0.1)/len(data),0.7,0.5),alpha=0.5,markersize=5)
                #     ax.plot(ang,rad,'+',color=cmap(i),alpha=0.3,markersize=5)
                if i % skip == 0:
                    ax.plot(ang,rad,'+',color=cmap(i),alpha=0.5,markersize=5)
                

            ax.plot(angular_m,radial_m,'o',color='r')

        #ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5),fontsize="7")

        #ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5),fontsize="7")
            ax.set_title(rf"$x_3=${d}")
        ax.set_rmax(max_radial+max_radial/10)
    #fig.legend(loc="lower right")
    fig.suptitle(title)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sub_ax = plt.axes([0.1, -.02, 0.8, 0.02])
    sample_len = len(list(datas.values())[0])
    fig.colorbar(sm, ticks=list(range(0, sample_len, int(sample_len/10))),location="bottom",cax=sub_ax)
    #plt.xticks([float(x) for x in datas[0].head()[1:-1]])
    plt.tight_layout(pad=0.4, w_pad=1, h_pad=1.0)

    #plt.show()
    plt.savefig("./output_polyakov.svg")


def create_figure_real_imag(datas,cols=4, mean = True, title=""):
    total = len(datas)
    rows = total // cols
    if total % cols != 0:
        rows += 1
        
    position = range(1,total+1)
    fig = plt.figure(1,figsize=FIG_SIZE)
    for k,(name,data) in enumerate(datas.items()):
        ax = fig.add_subplot(rows,cols,position[k])

        #delta = data.to_numpy().max()-data.to_numpy().min()
        delta = 0
        y = data["sum"].to_numpy()[1000:]
        x = data.index.to_numpy()[1000:]
        ax.plot(x,y.real,label="real")
        ax.plot(x,y.imag,label="imag")
        ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5),fontsize="7")
        ax.set_title(f"Beta={name}")
    fig.suptitle(title)

    #plt.xticks([float(x) for x in datas[0].head()[1:-1]])
    fig.tight_layout(pad=0.4, w_pad=1, h_pad=1.0)
    #plt.show()
    plt.savefig("./output_polyakov.svg")
    
def animate_polar(data):
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    max_magnitude = np.abs(data.values).max()
    cmap = plt.cm.jet
    index = np.array(list(map(int,data.columns)))
    norm = Normalize(vmin=index.min(), vmax=index.max())
        
    def update(frame):
        ax.clear()
        complex_numbers = data.iloc[frame, :]
        angles = np.angle(complex_numbers)
        magnitudes = np.abs(complex_numbers)
        ax.set_ylim(0, max_magnitude)  # Set a fixed y-limit for better comparison across frames
        for i, (theta, r) in enumerate(zip(angles, magnitudes)):
            ax.plot(theta, r, marker='o',color=cmap(norm(i)), label=f'z-index: {index[i]}')
        
    ani = FuncAnimation(fig, update, frames=len(data), interval=100, repeat=True)
    
    # Save the animation
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=10, metadata=dict(artist='Me'), bitrate=1800)
    ani.save('/home/haaaaron/SUN_twist_python_analysis/videos/polar_animation.mp4', writer=writer, dpi=100)
    
    # Close the figure to suppress plotting in Jupyter Notebook
    plt.close(fig)

if __name__ == "__main__":
    datas = {}
    for root,dirs,files in os.walk(sys.argv[1]):
        for name in files:
            print(os.path.join(root,name))
            name,data = load_data(os.path.join(root,name))
            print(name)
            data = data.drop(columns=['sum'])
            columns = list(data.columns.values)
            columns = columns[-len(columns)//2:] + columns[:-len(columns)//2]
            data = data.reindex(columns=columns)
            #data = data[50:]
            datas[name] = data
            
    # Sort dataframe based on keys
    myKeys = list(datas.keys())
    myKeys.sort()
    datas = {i: datas[i] for i in myKeys}
    #print(datas)
    create_figure(datas,cols=3,mean=True,title=sys.argv[2])
    # plt.imshow(data, cmap ="RdYlBu",aspect='auto')
    # plt.xticks([float(x) for x in data.head()[1:-1]])
    # plt.colorbar()
    # plt.show()