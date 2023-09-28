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
            if line.startswith('plaquette:'):
                line = line.split(':')[1]
                try:
                    data.append([float(x) for x in line.split()])
                except:
                    data.append(line.split())

    data = pd.DataFrame(data[1:],columns=data[0])
    name = "Beta: " + beta + ', Twist coefficient: ' + twist_c
    return name,data

def create_figure(datas,cols=4,mean=False):
    total = len(datas)
    rows = total // cols
    if total % cols != 0:
        rows += 1
        
    position = range(1,total+1)
    fig = plt.figure(1)
    fig.tight_layout(pad=5.0)
    for k,(name,data) in enumerate(datas.items()):
        ax = fig.add_subplot(rows,cols,position[k])
        #ax.imshow(datas[name], cmap ="RdYlBu",aspect='auto')
        delta = data.to_numpy().max()-data.to_numpy().min()
        #print(data.mean().to_numpy())
        #print(list(data.columns.values))
        if mean:
            ax.plot(data.columns.values,data.mean())
        else:
            sns.heatmap(data,ax=ax)
        ax.set_title(f"{name}, Delta={delta}")

    #plt.xticks([float(x) for x in datas[0].head()[1:-1]])
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
            data = data[50:]
            #data = data[50:]
            datas[name] = data
            
    # Sort dataframe based on keys
    myKeys = list(datas.keys())
    myKeys.sort()
    datas = {i: datas[i] for i in myKeys}
    
    create_figure(datas,mean=False)
