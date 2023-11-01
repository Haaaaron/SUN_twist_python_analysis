import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def create_average_action_figure(notwist,name_1,twist,name_2,thermalization=50):
    data = [[notwist,name_1],[twist,name_2]]
    for i,(datas,name) in enumerate(data):
        average_datas = {}
        for simulation in datas:
            average_sum = datas[simulation]['sum'].to_numpy().astype(float)
            average_datas[simulation] = average_sum[thermalization:].mean()
        data[i][0] = average_datas
    print(data)
    for (datas,name) in data:
        x = [float(x.split()[0]) for x in list(datas.keys())]
        print(x)
        y = list(datas.values())
        plt.plot(x,y,'o--',label=name)
        #plt.xticks([float(x) for x in datas[0].head()[1:-1]])
    plt.title(r'Wilson action total average with respect to $\beta$')
    plt.xlabel(r'$\beta$')
    plt.ylabel(r'$\langle$ S $\rangle$')
    plt.legend()
    plt.show()
    
def create_z_index_heat_map(datas,cols=4,mean=False):
    total = len(datas)
    rows = total // cols
    if total % cols != 0:
        rows += 1
        
    position = range(1,total+1)
    fig = plt.figure(1,figsize=(15,7))
    for k,(name,data) in enumerate(datas.items()):
        no_sum = data.drop(["sum"],axis=1)
        ax = fig.add_subplot(rows,cols,position[k])
        #ax.imshow(datas[name], cmap ="RdYlBu",aspect='auto')
        delta = no_sum.to_numpy().max()-no_sum.to_numpy().min()
        #print(data.mean().to_numpy())
        #print(list(data.columns.values))
        if mean:
            ax.plot(no_sum[50:].columns.values,no_sum[50:].mean())
        else:
            sns.heatmap(no_sum,ax=ax)
        ax.set_title(f"{name}")
    plt.tight_layout(pad=0.4, w_pad=3, h_pad=1.0)
    #plt.xticks([float(x) for x in datas[0].head()[1:-1]])
    plt.show()