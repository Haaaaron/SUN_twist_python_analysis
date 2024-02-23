import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import seaborn as sns
from scipy import integrate

def create_average_action_figure(notwist,name_1,twist,name_2,errors_1=None,errors_2=None,thermalization=50):
    plt.figure(figsize=(13,10))
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
    plt.savefig('./output.svg') 
    #plt.show()
    
def create_average_action_figure_jackknife(notwist,name_1,twist,name_2,errors_1=None,errors_2=None,FS=False,plot_fmt="-",fig_text=None):
    data = [[notwist,name_1,errors_1],[twist,name_2,errors_2]]
    plt.figure(figsize=(13,10))
    #plt.figure()

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
    if fig_text: plt.figtext(0,0,fig_text,fontstyle="italic")
    #plt.legend()
    plt.savefig('./output.svg') 

    
def create_twist_notwist_difference_figure_jackknife(notwist,twist,errors_1=None,errors_2=None,FS=False,plot_fmt="-",fig_text=None):
    plt.figure(figsize=(13,10))
    for (y_notwist,y_twist) in zip(notwist[1],twist[1]):
        plt.plot(notwist[0],y_twist-y_notwist,plot_fmt)
        #break

        #plt.xticks([float(x) for x in datas[0].head()[1:-1]])
    if FS:    
        plt.title(r'Wilson action difference between twist and notwist with respect to $\beta$ using jackknife and FS reweighting')
    else:
        plt.title(r'Wilson action total average with respect to $\beta$ produced with jackknife')
    plt.xlabel(r'$\beta$')
    plt.ylabel(r'$\langle S_{t} - S_{nt} \rangle$')
    if fig_text: plt.figtext(0,0,fig_text,fontstyle="italic")
    #plt.legend()
    plt.savefig('./output.svg') 

    
def create_integral_figure_jackknife(notwist,twist,dim,mean=True,plot_fmt='-',fig_text=None,N=10,FS=False):
    integrated_ys = []
    plt.figure(figsize=(13,10))
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
        volume=dim[0]*dim[1]*dim[2]*dim[3]
        plt.errorbar(x,y_int*volume*6/(dim[0]*dim[1]),yerr=np.sqrt(N-1)*y_err*volume*6/(dim[0]*dim[1]),fmt=plot_fmt,ecolor="r")
    if FS: 
        plt.title(r'Integral of difference between twist and no-twist wilson action. Averaged over set of jackknife data paired with FS reweighting.')
    else: 
        plt.title(r'Integral of difference between twist and no-twist wilson action. Averaged over set of jackknife data.')

    plt.xlabel(r'$\beta$')
    plt.ylabel(r'$\int \langle S_{t} - S_{nt} \rangle$')
    if fig_text: plt.figtext(0,0,fig_text,fontstyle="italic")
    #plt.legend()
    plt.savefig('./output.svg') 


def create_z_index_heat_map(datas,cols=4,mean=False,series_range=(10000,-1)):
    total = len(datas)
    rows = total // cols
    if total % cols != 0:
        rows += 1
        
    position = range(1,total+1)
    fig = plt.figure(1,figsize=(6,6))
    for k,(name,data) in enumerate(datas.items()):
        no_sum = data.drop(["sum"],axis=1)[series_range[0]:series_range[1]]
        ax = fig.add_subplot(rows,cols,position[k])
        #ax.imshow(datas[name], cmap ="RdYlBu",aspect='auto')
        #delta = no_sum.to_numpy().max()-no_sum.to_numpy().min()
        #print(data.mean().to_numpy())
        #print(list(data.columns.values))
        #print(no_sum)
        if mean:
            ax.plot(no_sum.columns.values,no_sum.mean())
        else:
            sns.heatmap(no_sum,ax=ax)
        ax.set_title(f"{name}")
    plt.tight_layout(pad=0.4, w_pad=3, h_pad=1.0)
    #plt.xticks([float(x) for x in datas[0].head()[1:-1]])
    plt.savefig('./output.svg') 

    
def plot_action_series(datas,cols=4,series_range=(10000,-1)):
    total = len(datas)
    rows = total // cols
    if total % cols != 0:
        rows += 1
        
    position = range(1,total+1)
    fig = plt.figure(1,figsize=(9,6))
    for k,(name,data) in enumerate(datas.items()):
        only_sum = data["sum"][series_range[0]:series_range[1]]
        ax = fig.add_subplot(rows,cols,position[k])
        #ax.imshow(datas[name], cmap ="RdYlBu",aspect='auto')
       # delta = no_sum.to_numpy().max()-no_sum.to_numpy().min()
        #print(data.mean().to_numpy())
        #print(list(data.columns.values))
        ax.plot(list(range(len(only_sum))),only_sum)
        ax.set_title(f"{name}")
    plt.tight_layout(pad=0.4, w_pad=3, h_pad=1.0)
    #plt.xticks([float(x) for x in datas[0].head()[1:-1]])
    plt.show()
    
def create_amplitude_graph(notwist,twist,extra,plot=True):
    f_1 = np.mean(np.array(twist[1]),axis=0)
    f_2 = np.mean(np.array(notwist[1]),axis=0)
    f = np.mean(np.array(twist[1] - notwist[1]),axis=0)
    x = twist[0]
    df = np.gradient(f,x)
    #plt.plot(x,df)
    figure, axis = plt.subplots(3,figsize=(9,10)) 
    index_df_min = np.argmax(df)
    index_df_max = np.argmin(df)
    if plot:
        f_1_plot = f_1[index_df_max-extra:index_df_min+extra]
        f_2_plot = f_2[index_df_max-extra:index_df_min+extra]
        f_plot = f[index_df_max-extra:index_df_min+extra]
        x_plot = x[index_df_max-extra:index_df_min+extra]
        df_plot = df[index_df_max-extra:index_df_min+extra]
        axis[0].plot(x_plot,f_1_plot,"-r",label=r"$S_{t=1}$")
        axis[0].plot(x_plot,f_2_plot,"-g",label=r"$S_{t=2}$")
        axis[0].plot(x_plot[-extra],f_2_plot[-extra],"om",)
        axis[0].plot(x_plot[extra],f_1_plot[extra],"oy")
        axis[0].set(ylabel=r"$S$")

        axis[1].plot(x_plot,f_plot)
        axis[1].plot(x_plot[-extra],f_plot[-extra],"om")
        axis[1].plot(x_plot[extra],f_plot[extra],"oy")
        axis[1].plot(x_plot[np.argmax(f_plot)],np.max(f_plot),"og")
        axis[1].set(ylabel=r"$\overline{S}(\beta)=S_{t=2}-S_{t=1}$")
        
        axis[2].plot(x_plot[-extra],df_plot[-extra],"om",label=r"$\max(\overline{S}'(\beta))$")
        axis[2].plot(x_plot[extra],df_plot[extra],"oy",label=r"$\min(\overline{S}'(\beta))$")
        axis[2].plot(x_plot,df_plot)
        axis[2].set(ylabel=r"$\overline{S}'(\beta)$")

        #plt.setp(axis[1],ylabel="$\bar{S}(\beta)$")
    figure.gca().set_xlabel(r"$\beta$")
    figure.suptitle(r"Phase transition 32,32,48,6")
    figure.legend()
    plt.savefig('./output.svg') 

    return x[index_df_max]-x[index_df_min],np.max(df),x[index_df_min],x[index_df_max]
