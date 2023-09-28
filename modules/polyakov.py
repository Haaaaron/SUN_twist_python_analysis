import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
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
            if line.startswith('polyakov:'):
                line = line.split(':')[1]
                try:
                    numbers = [x for x in line.split(",")][:-1]
                    complex_converted = []
                    for x in numbers:
                        num = x.split()
                        complex_converted.append(complex(float(num[0]),float(num[1])))
                        #complex(num[0])
                    data.append(complex_converted)
                except:
                    data.append(line.split(",")[:-1])
    #print(data[1])
    data[0] = data[0] + ["sum"]
    data = pd.DataFrame(data[1:],columns=data[0])
    name = beta + " " + twist_c
    return name,data

def create_figure(datas,cols=4, mean = True, title=""):
    total = len(datas)
    rows = total // cols
    if total % cols != 0:
        rows += 1
        
    position = range(1,total+1)
    fig = plt.figure(1)
    for k,(name,data) in enumerate(datas.items()):
        ax = fig.add_subplot(rows,cols,position[k],projection="polar")

        #delta = data.to_numpy().max()-data.to_numpy().min()
        delta = 0
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
        ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5),fontsize="7")
        ax.set_title(f"Beta={name}")
    fig.suptitle(title)

    #plt.xticks([float(x) for x in datas[0].head()[1:-1]])
    fig.tight_layout(pad=5.0, h_pad=5)
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