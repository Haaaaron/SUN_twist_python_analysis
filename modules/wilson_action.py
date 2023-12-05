import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import seaborn as sns
from scipy import integrate

def create_average_action_figure(notwist,name_1,twist,name_2,errors_1=None,errors_2=None,thermalization=50):
    data = [[notwist,name_1,errors_1],[twist,name_2,errors_2]]
    for i,(datas,name,errors) in enumerate(data):
        average_datas = {}
        for simulation in datas:
            average_sum = datas[simulation]['sum'].to_numpy().astype(float)
            average_datas[simulation] = average_sum[thermalization:].mean()
        data[i][0] = average_datas
    for (datas,name,errors) in data:
        x = [float(x.split()[0]) for x in list(datas.keys())]
        y = list(datas.values())
        if (errors == None):
            plt.plot(x,y,'o--',label=name)
        else:
            error = []
            for key in datas.keys():
                error.append(float(errors[key]["error"].iloc[-1]))
            plt.errorbar(x,y,yerr=error,fmt='o--',label=name)
        #plt.xticks([float(x) for x in datas[0].head()[1:-1]])
    plt.title(r'Wilson action total average with respect to $\beta$')
    plt.xlabel(r'$\beta$')
    plt.ylabel(r'$\langle$ S $\rangle$')
    plt.legend()
    plt.show()
    
def create_average_action_figure_jackknife(notwist,name_1,twist,name_2,errors_1=None,errors_2=None,FS=False,plot_fmt="-"):
    data = [[notwist,name_1,errors_1],[twist,name_2,errors_2]]
    plt.figure(figsize=(12,10))

    for (datas,name,errors) in data:
        if (errors == None):
            for y in datas[1]:
                plt.plot(datas[0],y,plot_fmt,label=name)
                #break
        else:
            error = []
            for key in datas.keys():
                error.append(float(errors[key]["error"].iloc[-1]))
            plt.errorbar(x,y,yerr=error,fmt='o--',label=name)
        #plt.xticks([float(x) for x in datas[0].head()[1:-1]])
    if FS:
        plt.title(r'Wilson action total average with respect to $\beta$ produced with jackknife and FS reweighting')
    else:
        plt.title(r'Wilson action total average with respect to $\beta$ produced with jackknife')
    plt.xlabel(r'$\beta$')
    plt.ylabel(r'$\langle$ S $\rangle$')
    #plt.legend()
    plt.show()
    
def create_twist_notwist_difference_figure_jackknife(notwist,twist,errors_1=None,errors_2=None,plot_fmt="-"):
    plt.figure(figsize=(12,10))
    for (y_notwist,y_twist) in zip(notwist[1],twist[1]):
        plt.plot(notwist[0],y_twist-y_notwist,plot_fmt)
        #break

        #plt.xticks([float(x) for x in datas[0].head()[1:-1]])
    plt.title(r'Wilson action difference between twist and notwist with respect to $\beta$ using jackknife and FS reweighting')
    plt.xlabel(r'$\beta$')
    plt.ylabel(r'$\langle S_{t} - S_{nt} \rangle$')
    #plt.legend()
    plt.show()
    
def create_integral_figure_jackknife(notwist,twist,mean=True,plot_fmt='-',N=10):
    integrated_ys = []
    plt.figure(figsize=(12,10))
    for (y_notwist,y_twist) in zip(notwist[1],twist[1]):
        x=notwist[0]
        y = y_twist-y_notwist
        y_int = integrate.cumtrapz(y,x,initial=0)
        if mean:
            integrated_ys.append(y_int)
        else:
            plt.plot(x,y_int,plot_fmt)
        #y_int = integrate.simpson(y,x)
    y_mean = np.mean(integrated_ys)
    y_err = np.std(integrated_ys,axis=0)
    print(y_err)
    if mean:
        plt.errorbar(x,y_int,yerr=np.sqrt(N-1)*y_err,fmt=plot_fmt,ecolor="r")
    plt.title(r'Integral of difference between twist and no-twist wilson action. Averaged over set of jackknife data.')
    plt.xlabel(r'$\beta$')
    plt.ylabel(r'$\int \langle S_{t} - S_{nt} \rangle$')
    #plt.legend()
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
    
