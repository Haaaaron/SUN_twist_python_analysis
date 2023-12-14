import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import os

def create_figure_polar(datas,cols=4, mean = True, title=""):
    max_radial = 0
    #Looking for maximal polar radius to have even plot
    for k,(name,data) in enumerate(datas.items()):
        data = data.drop(columns=['sum'])
        if mean:
            max_radial = max(max_radial,np.max(abs(data.mean().to_numpy())))
        else:
            max_radial = max(max_radial,np.max(abs(data.to_numpy())))
    total = len(datas)
    rows = total // cols
    if total % cols != 0:
        rows += 1
        
    position = range(1,total+1)
    fig = plt.figure(1,figsize=(40,40))
    for k,(name,data) in enumerate(datas.items()):
        ax = fig.add_subplot(rows,cols,position[k],projection="polar")

        #delta = data.to_numpy().max()-data.to_numpy().min()
        delta = 0
        data = data.drop(columns=['sum'])
        if mean:
            angular = np.angle(data.mean().to_numpy())
            radial = abs(data.mean().to_numpy())
            for i,(ang,rad) in enumerate(zip(angular,radial)):
                ax.plot(ang,rad,'o',label=f'{i}')
            
        else:
            for column in data:
                angular = np.angle(data[column].to_numpy())
                radial = abs(data[column].to_numpy())
                #print(radial)
                ax.plot(angular,radial, 'o', label=column)
        #ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5),fontsize="7")
        ax.set_title(f"Beta={name}")
        ax.set_rmax(max_radial)

    fig.suptitle(title)
    print("something")
    #plt.xticks([float(x) for x in datas[0].head()[1:-1]])
    fig.tight_layout(pad=0.4, w_pad=1, h_pad=1.0)
    plt.show()

def create_figure_real_imag(datas,cols=4, mean = True, title=""):
    total = len(datas)
    rows = total // cols
    if total % cols != 0:
        rows += 1
        
    position = range(1,total+1)
    fig = plt.figure(1,figsize=(20,10))
    for k,(name,data) in enumerate(datas.items()):
        ax = fig.add_subplot(rows,cols,position[k])

        #delta = data.to_numpy().max()-data.to_numpy().min()
        delta = 0
        y = data["sum"].to_numpy()
        x = data.index.to_numpy()
        ax.plot(x,y.real,label="real")
        ax.plot(x,y.imag,label="imag")
        ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5),fontsize="7")
        ax.set_title(f"Beta={name}")
    fig.suptitle(title)

    #plt.xticks([float(x) for x in datas[0].head()[1:-1]])
    fig.tight_layout(pad=0.4, w_pad=1, h_pad=1.0)
    plt.show()

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