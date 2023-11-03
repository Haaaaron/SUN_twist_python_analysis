import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import os

def load_data(file,empty):
    f = open(file,'r')
    data = []
    beginRead = False
    beta = None
    for line in f:
        if line.startswith('beta'):
           beta = float(line.split()[1])
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
                    
    try:
        data = pd.DataFrame(data[1:],columns=data[0])
    except:
        empty +=1
    #name = "Beta: " + str(beta) + ', Twist coefficient: ' + twist_c
    return beta,data,empty

def create_figure(data):
    for (datas,name) in data:
        x = list(datas.keys())
        y = list(datas.values())
        plt.plot(x,y,'o--',label=name)
        #plt.xticks([float(x) for x in datas[0].head()[1:-1]])
    plt.title(r'Wilson action total average with respect to $\beta$')
    plt.xlabel(r'$\beta$')
    plt.ylabel(r'$\langle$ S $\rangle$')
    plt.legend()
    plt.show()

if __name__ == "__main__":
    datas_1 = {}
    datas_2 = {}
    name_1 = sys.argv[3]
    name_2 = sys.argv[4]
    empty = 0
    print(f"Computing {name_1}")
    for root,dirs,files in os.walk(sys.argv[1]):
        for name in files:
            #print(os.path.join(root,name))
            name,data,empty_ret = load_data(os.path.join(root,name),empty)
            if empty == empty_ret:
                data = data['sum'].to_numpy()
                #columns = list(data.columns.values)
                #columns = columns[-len(columns)//2:] + columns[:-len(columns)//2]
                #data = data.reindex(columns=columns)
                data = data[50:].mean()
                #data = data[50:]
                datas_1[name] = data
            else:
                empty = empty_ret
    print(f"Computing {name_2}")
    for root,dirs,files in os.walk(sys.argv[2]):
        for name in files:
            #print(os.path.join(root,name))
            
            name,data,empty_ret = load_data(os.path.join(root,name),empty)
            if empty == empty_ret:
                data = data['sum'].to_numpy()
                #columns = list(data.columns.values)
                #columns = columns[-len(columns)//2:] + columns[:-len(columns)//2]
                #data = data.reindex(columns=columns)
                data = data[50:].mean()
                #data = data[50:]
                datas_2[name] = data
            else:
                empty = empty_ret
            
    # Sort dataframe based on keys
    data = [[datas_1,name_1],[datas_2,name_2]]
    for i,(datas,name) in enumerate(data):
        myKeys = list(datas.keys())
        myKeys.sort()
        datas = {i: datas[i] for i in myKeys}
        data[i][0] = datas
    print(data[0],data[1])
    print(empty)
    create_figure(data)
