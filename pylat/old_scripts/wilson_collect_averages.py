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


if __name__ == "__main__":
    datas = {}
    empty = 0
    for root,dirs,files in os.walk(sys.argv[1]):
        for name in files:
            #print(os.path.join(root,name))
            
            name,data,empty_ret = load_data(os.path.join(root,name),empty)
            if empty == empty_ret:
                #print(data)
                data = data.drop(columns=['sum'])
                columns = list(data.columns.values)
                columns = columns[-len(columns)//2:] + columns[:-len(columns)//2]
                data = data.reindex(columns=columns)
                data = data[50:].mean()
                #data = data[50:]
                datas[name] = data
            else:
                empty = empty_ret
            
    # Sort dataframe based on keys

    myKeys = list(datas.keys())
    myKeys.sort()
    datas = {i: datas[i] for i in myKeys}
    df = pd.DataFrame.from_dict(datas).T
    #print(df.index)
    df.to_csv(sys.argv[2])