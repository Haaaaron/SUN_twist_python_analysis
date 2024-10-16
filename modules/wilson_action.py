import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import seaborn as sns
from scipy import integrate
FIG_SIZE=(7,5)
#standard font 11
plt.rcParams.update({'font.size': 11})

def create_average_action_figure(notwist,name_1,twist,name_2,errors_1=None,errors_2=None,thermalization=50):
    plt.figure(figsize=FIG_SIZE)
    #data = [[notwist,name_1,errors_1],[twist,name_2,errors_2]]
    data = [[notwist,name_1,errors_1]]

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
            plt.plot(x,y,'.--',label=name)
        else:
            error = []
            for key in datas.keys():
                error.append(float(errors[key]["error"].iloc[-1]))
            plt.errorbar(x,y,yerr=error,fmt='.--',ecolor="r",label=name)
        #plt.xticks([float(x) for x in datas[0].head()[1:-1]])
    plt.xlabel(r'$\beta$')
    plt.ylabel(r'$\langle$ S $\rangle / (6\cdot V \cdot\beta)$')
    #plt.legend()
    plt.savefig('./output.svg') 
    #plt.show()
    
def create_average_action_figure_difference(notwist,name_1,twist,name_2,errors_1=None,errors_2=None,covariance=None,thermalization=50):
    plt.figure(figsize=FIG_SIZE)
    data = [[notwist,name_1,errors_1],[twist,name_2,errors_2]]
    #data = [[notwist,name_1,errors_1]]

    for i,(datas,name,errors) in enumerate(data):
        average_datas = {}
        for simulation in datas:
            average_sum = datas[simulation]['sum'].to_numpy().astype(float)
            average_datas[simulation] = average_sum[thermalization:].mean()
        data[i][0] = average_datas
    x = [float(x.split()[0]) for x in list(data[0][0].keys())]
    y = np.array(list(data[1][0].values())) - np.array(list(data[0][0].values()))
    error=[]
    for i,key in enumerate(zip(data[0][0].keys(),data[1][0].keys())):
        error.append(np.sqrt(float(data[0][2][key[0]]["error"].iloc[-1])**2+float(data[1][2][key[1]]["error"].iloc[-1])**2))
    if (data[0][2] == None):
        plt.plot(x,y,'o--',label=name)
    else:        
        plt.errorbar(x,y,yerr=error,fmt='.--',ecolor="r")
        #plt.xticks([float(x) for x in datas[0].head()[1:-1]])
    plt.xlabel(r'$\beta$')
    plt.ylabel(r'$ \frac{1}{6\cdot V \cdot\beta}\langle S_{t} - S_{nt} \rangle $')
    plt.savefig('./output.svg') 
    #plt.show()
    
def create_average_action_figure_jackknife(notwist,name_1,twist,name_2,error=True,FS=False,TITLE=False,plot_fmt="-",fig_text=None):
    data = [[notwist,name_1],[twist,name_2]]
    plt.figure(figsize=FIG_SIZE)
    #plt.figure()

    for (datas,name) in data:
        if (error == False):
            # for i,y in enumerate(datas[1]):
            #     plt.plot(datas[0],y,linestyle="dotted",label=f"{i}")
            #     #break
            plt.plot(datas[0],np.mean(datas[1],axis=0),label=name)
        else:
            mean = np.mean(datas[1],axis=0)
            M = len(datas[1])
            errors= np.zeros(len(datas[0]))
            for block in datas[1]:
                errors += (block - mean)**2
            errors = np.sqrt((M-1)/M*errors)
            plt.errorbar(datas[0],mean,yerr=errors,ecolor="r",fmt='-',label=name)
        #plt.xticks([float(x) for x in datas[0].head()[1:-1]])
    if FS:
        plt.title(r'Wilson action total average with respect to $\beta$ produced with jackknife and FS reweighting')
    if TITLE:
        plt.title(r'Wilson action average with respect to $\beta$')
    plt.xlabel(r'$\beta$')
    plt.ylabel(r'$\langle$ S $\rangle / (6\cdot V \cdot\beta)$')
    if fig_text: plt.figtext(0.1,-0.01,fig_text,fontstyle="italic")
    plt.legend()
    plt.savefig('./output.svg',dpi=600, bbox_inches = "tight") 


    
def create_twist_notwist_difference_figure_jackknife(notwist,twist_list,jk=False,error=True,FS=False,TITLE=False,plot_fmt="-",fig_text=None):
    plt.figure(figsize=FIG_SIZE)
    for twist,name in twist_list:
        if jk == True:
            for i,(y_notwist,y_twist) in enumerate(zip(notwist[1],twist[1])):
                plt.plot(notwist[0],y_twist-y_notwist,linestyle="dotted",label=f"{i}")
        if error == True:
            difference = twist[1] - notwist[1]
            M = len(notwist[0])
            mean = np.mean(twist[1],axis=0)-np.mean(notwist[1],axis=0)
            errors= np.zeros(len(notwist[0]))
            for block in difference:
                print(block[-1])
                errors += (block - mean)**2
            errors = np.sqrt((M-1)/M*errors)
            print(np.max(errors))
            if len(twist_list) != 1: plt.errorbar(notwist[0],mean,yerr=errors,fmt='-',label=name)
            else: plt.errorbar(notwist[0],mean,yerr=errors,ecolor='r',fmt='-',label=name)
            #plt.xticks([float(x) for x in datas[0].head()[1:-1]])
    if FS:    
        plt.title(r'Wilson action difference between twist and notwist with respect to $\beta$ using jackknife and FS reweighting')
    if TITLE:
        plt.title(r'Wilson action difference between twist and notwist with respect to $\beta$')
    plt.xlabel(r'$\beta$')
    plt.ylabel(r'$\frac{1}{6\cdot V\cdot \beta}\langle S_{t} - S_{nt} \rangle$')
    if len(twist_list) != 1: plt.legend()
    if fig_text: plt.figtext(0.1,-0.01,fig_text,fontstyle="italic")
    plt.tight_layout()
    plt.savefig('./output.svg',dpi=600, bbox_inches = "tight") 

    
def create_integral_figure_jackknife(notwist,twist_list,dim,mean=True,plot_fmt='-',fig_text=None,N=10,FS=False,TITLE=False):
    
    plt.figure(figsize=FIG_SIZE)
    #print(twist[1])
    for twist,name in twist_list:
        integrated_ys = []
        for (y_notwist,y_twist) in zip(notwist[1],twist[1]):
            x=notwist[0]
            y = y_twist-y_notwist
            #y_int = integrate.cumtrapz(y,x,initial=0)
            y_int = integrate.cumulative_trapezoid(y,x,initial=0)

            if mean:
                integrated_ys.append(y_int)
            else:
                plt.plot(x,y_int,plot_fmt)
            #y_int = integrate.simpson(y,x)
        y_mean = np.mean(integrated_ys,axis=0)
        errors= np.zeros(len(notwist[0]))
        for block in integrated_ys:
            errors += (block - y_mean)**2
        errors = np.sqrt((N-1)/N*errors)
        if mean:
            volume=dim[0]*dim[1]*dim[2]*dim[3]
            if len(twist_list) != 1: plt.errorbar(x,y_mean*volume*6*dim[3]**2/(dim[0]*dim[1]),yerr=errors*volume*6*dim[3]**2/(dim[0]*dim[1]),fmt=plot_fmt,capsize=1,label=name)
            else: plt.errorbar(x,y_mean*volume*6*dim[3]**2/(dim[0]*dim[1]),yerr=errors*volume*6*dim[3]**2/(dim[0]*dim[1]),ecolor='r',fmt=plot_fmt,capsize=1,label=name)
        print(name,y_mean[-1]*volume*6*dim[3]**2/(dim[0]*dim[1]))
    if FS: 
        plt.title(r'Integral of difference between twist and no-twist wilson action. Averaged over set of jackknife data paired with FS reweighting.')
    if TITLE: 
        plt.title(r'Integral of difference between twist and no-twist wilson action')

    plt.xlabel(r'$\beta$')
    plt.ylabel(r'$\alpha_{o-o}/T^3$')
    if fig_text: plt.figtext(0.1,-0.01,fig_text,fontstyle="italic")
    if len(twist_list) != 1: plt.legend()
    
    #plt.legend()
    plt.savefig('./output.svg',dpi=600, bbox_inches = "tight") 

def create_integral_figure_ratio(notwist,twist_list,dim,mean=True,plot_fmt='-',fig_text=None,N=10,FS=False,TITLE=False):
    
    plt.figure(figsize=FIG_SIZE)
    #print(twist[1])
    ys = {}
    for twist,name in twist_list:
        integrated_ys = []
        for (y_notwist,y_twist) in zip(notwist[1],twist[1]):
            x=notwist[0]
            y = y_twist-y_notwist
            #y_int = integrate.cumtrapz(y,x,initial=0)
            y_int = integrate.cumulative_trapezoid(y,x,initial=0)

            integrated_ys.append(y_int)

        ys[name] = integrated_ys
            #y_int = integrate.simpson(y,x)
    volume=dim[0]*dim[1]*dim[2]*dim[3]
    ys_ratio = []
    for y_1,y_2 in zip(ys[twist_list[0][1]],ys[twist_list[1][1]]):
        ys_ratio.append(y_2/y_1)
    
    y_mean =np.mean(ys_ratio,axis=0)
    errors= np.zeros(len(notwist[0]))
    for block in ys_ratio:
        errors += (block - y_mean)**2
    errors = np.sqrt((N-1)/N*errors)
    start = np.where(x==10.787)[0][0] 
    beta_critical = 10.787 
    critical = convert_beta_to_T([beta_critical])[0]
    T = convert_beta_to_T(x[start:])
    #print(critical/T)
    print(y_mean[-1],errors[-1])
    if mean:
        plt.errorbar(critical/T,y_mean[start:],yerr=errors[start:],ecolor='r',fmt=plot_fmt,capsize=1)
    if FS: 
        plt.title(r'Integral of difference between twist and no-twist wilson action. Averaged over set of jackknife data paired with FS reweighting.')
    if TITLE: 
        plt.title(r'Integral of difference between twist and no-twist wilson action')
    plt.plot(critical/T,np.ones(len(T))*4/3,label=r'$\sigma_{k=2}/\sigma_1=4/3$')
    plt.xlabel(r'$T/T_c$')
    plt.ylabel(r'$ \alpha_{z_2}/\alpha_{z_1} $')
    if fig_text: plt.figtext(0.1,-0.01,fig_text,fontstyle="italic")

    plt.legend()
    plt.savefig('./output.svg',dpi=600, bbox_inches = "tight") 


def create_z_index_heat_map(datas,cols=2,mean=False,series_range=(10000,-1)):
    total = len(datas)
    rows = total // cols
    if total % cols != 0:
        rows += 1
        
    position = range(1,total+1)
    fig = plt.figure(1,figsize=FIG_SIZE)
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
        ax.set_title(r"$\beta$={}".format(name.split()[0]))
        ax.set_xlabel(r"$x_3$-index")
    plt.tight_layout(pad=0.4, w_pad=3, h_pad=1.0)
    
    #plt.xticks([float(x) for x in datas[0].head()[1:-1]])
    plt.savefig('./output.svg') 

    
def plot_action_series(datas,cols=4,series_range=(10000,-1)):
    total = len(datas)
    rows = total // cols
    if total % cols != 0:
        rows += 1
        
    position = range(1,total+1)
    fig = plt.figure(1,figsize=FIG_SIZE)
    for k,(name,data) in enumerate(datas.items()):
        only_sum = data["sum"][series_range[0]:series_range[1]]
        ax = fig.add_subplot(rows,cols,position[k])
        #ax.imshow(datas[name], cmap ="RdYlBu",aspect='auto')
       # delta = no_sum.to_numpy().max()-no_sum.to_numpy().min()
        #print(data.mean().to_numpy())
        #print(list(data.columns.values))
        ax.plot(list(range(len(only_sum))),only_sum)
        ax.set_title(r"$\beta=${}".format(name.split()[0]))
    plt.tight_layout(pad=0.4, w_pad=3, h_pad=1.0)
    #plt.xticks([float(x) for x in datas[0].head()[1:-1]])
    plt.savefig('./output.svg') 
    
def create_amplitude_graph(notwist,twist,extra,plot=True):
    f_1 = np.mean(np.array(twist[1]),axis=0)
    f_2 = np.mean(np.array(notwist[1]),axis=0)
    f = np.mean(np.array(twist[1] - notwist[1]),axis=0)
    x = twist[0]
    df = np.gradient(f,x)
    #plt.plot(x,df)
    figure, axis = plt.subplots(3,figsize=FIG_SIZE) 
    index_df_min = np.argmax(df)
    index_df_max = np.argmin(df)
    if plot:
        f_1_plot = f_1[index_df_max-extra:index_df_min+extra]
        f_2_plot = f_2[index_df_max-extra:index_df_min+extra]
        f_plot = f[index_df_max-extra:index_df_min+extra]
        x_plot = x[index_df_max-extra:index_df_min+extra]
        df_plot = df[index_df_max-extra:index_df_min+extra]
        axis[0].plot(x_plot,f_1_plot,"-r")
        axis[0].plot(x_plot,f_2_plot,"-g")
        axis[0].plot(x_plot[-extra],f_2_plot[-extra],"om",)
        axis[0].plot(x_plot[extra],f_1_plot[extra],"oy")
        axis[0].set(ylabel=r"S = $\frac{1}{6\cdot V \cdot\beta}\langle$ S $\rangle$")

        axis[1].plot(x_plot,f_plot)
        axis[1].plot(x_plot[-extra],f_plot[-extra],"om")
        axis[1].plot(x_plot[extra],f_plot[extra],"oy")
        axis[1].plot(x_plot[np.argmax(f_plot)],np.max(f_plot),"og",label=r"$\max(\overline{S}(\beta))$")
        axis[1].set(ylabel=r"$\overline{S}(\beta)=S_{t}-S_{nt}$")
        
        axis[2].plot(x_plot[-extra],df_plot[-extra],"om",label=r"$\max(\overline{S}'(\beta))$")
        axis[2].plot(x_plot[extra],df_plot[extra],"oy",label=r"$\min(\overline{S}'(\beta))$")
        axis[2].plot(x_plot,df_plot)
        axis[2].set(ylabel=r"$\overline{S}'(\beta)$")

        #plt.setp(axis[1],ylabel="$\bar{S}(\beta)$")
    figure.gca().set_xlabel(r"$\beta$")
    figure.suptitle(r"")
    figure.legend()
    plt.savefig('./output.svg') 

    return x[index_df_max]-x[index_df_min],np.max(df),x[index_df_min],x[index_df_max],df_plot[-extra],df_plot[extra]

def convert_beta_to_T(beta,beta_0=10.789,c_0=0.2702,c_1=3.5393,c_2=18.790,c_3=24,N=4):
    beta= np.array(beta)
    x = N**2*(1/beta - 1/beta_0)
    a_sigma = c_0*np.exp(-12*np.pi**2/11*(beta-beta_0)/N**2+c_1*x+c_2*x**2)

    return a_sigma