import numpy as np
import pandas as pd
import subprocess
from IPython.display import Markdown, display
import os


def display_markdown_title(title_info):
    title_info = title_info.split("/")
    group = title_info[3].split('-')[0][2]
    volume = title_info[3].split('-')[1:]
    temp = title_info[4].split('-')[1]
    twist_coeff = title_info[4].split('-')[3]
    title_markdown = r"# $\text{{SU}}({})$, $V$ = {}, $\beta=$ {}, twist coeff $=$ {}".format(group, volume, temp, twist_coeff)
    display(Markdown(title_markdown))

def select_subset(data,x,y,stride=1):
    #Select subset of simulations in plaquette data
    keys = [key for key in data]
    keys = keys[x:y:stride]
    subset_data = {}
    try:
        for key in keys:
            subset_data[key] = data[key]
    except:
        subset_data[keys] = data[keys]
        
    return subset_data

def compute_with_aa(data,thermalization=50):
    error_dict = {}
    for i,(name,datas) in enumerate(data.items()):
        np.savetxt("./modules/error_temp/temp_file.txt",datas.to_numpy()[thermalization:,:])
        error_data = subprocess.run(["/home/haaaaron/bin/aa /home/haaaaron/SUN_twist_python_analysis/modules/error_temp/temp_file.txt"],text=True,shell=True,capture_output=True).stdout.split("\n")
        empty_array = []
        for line in error_data:
            if len(line )!= 0:
                empty_array.append(line.split())
        error_data = np.array(empty_array)
        error_data = np.delete(error_data, 0,1)
        error_data = pd.DataFrame(error_data,columns=["average","error","atocorrelation","autocorrelation", "ind.meas"])
        error_dict[name] = error_data
    return error_dict

def compute_with_aa_jackknife(data,column,bins,only_sum=True,thermalization=1000):
    error_dict = {}
    for i,(name,datas) in enumerate(data.items()):
        #print(datas.to_numpy()[thermalization:])
        if only_sum: 
            np.savetxt("./modules/error_temp/temp_file.txt",datas["sum"].to_numpy()[thermalization:])
            column=1
        else:
            np.savetxt("./modules/error_temp/temp_file.txt",datas.to_numpy()[thermalization:])
        error_data = subprocess.run([f"/home/haaaaron/bin/aa -d {column} -j {bins} /home/haaaaron/SUN_twist_python_analysis/modules/error_temp/temp_file.txt"],text=True,shell=True,capture_output=True).stdout.split("\n")
        jackknife_data = subprocess.run([f"/home/haaaaron/bin/aa -d {column} -J {bins} /home/haaaaron/SUN_twist_python_analysis/modules/error_temp/temp_file.txt"],text=True,shell=True,capture_output=True).stdout.split("\n")        
        empty_array = []
        for line in jackknife_data:
            if len(line )!= 0:
                empty_array.append(line.split()[1:2])
        error_dict[name] = (error_data[0].split()[1:3],empty_array)
    return error_dict

def compute_with_aa_jackknife_fourier(fourier_profile,bins,only_sum=True,thermalization=1000):
    error_dict = {}

    # for i,frame in enumerate(data):
    #print(datas.to_numpy()[thermalization:])
    fourier_profile_average = []
    errors = []
    np.savetxt("/home/haaaaron/SUN_twist_python_analysis/modules/error_temp/temp_file.txt",fourier_profile[thermalization:])
    for i in range(1,np.shape(fourier_profile)[1]+1):
        error_data = subprocess.run([f"/home/haaaaron/bin/aa -d {i} -j {bins} /home/haaaaron/SUN_twist_python_analysis/modules/error_temp/temp_file.txt"],text=True,shell=True,capture_output=True).stdout.split("\n")
        jackknife_data = subprocess.run([f"/home/haaaaron/bin/aa -d {i} -J {bins} /home/haaaaron/SUN_twist_python_analysis/modules/error_temp/temp_file.txt"],text=True,shell=True,capture_output=True).stdout.split("\n")        
        empty_array = []
        for line in jackknife_data:
            if len(line )!= 0:
                empty_array.append(float(line.split()[1]))
        errors.append(float(error_data[0].split()[2]))
        fourier_profile_average.append(np.array(empty_array).mean())
    #error_dict[name] = (error_data[0].split()[1:3],empty_array)
    return (fourier_profile_average,errors)

def compute_with_aa_autocorrelation(data,thermalization=1000):
    autocorr_dict = {}
    for i,(name,datas) in enumerate(data.items()):
        #print(datas.to_numpy()[thermalization:])
        np.savetxt("./modules/autocorrelation_temp/temp_file.txt",datas["sum"].to_numpy()[thermalization:])
        column=1

        aa_output = subprocess.run([f"/home/haaaaron/bin/aa -d {column} /home/haaaaron/SUN_twist_python_analysis/modules/autocorrelation_temp/temp_file.txt"],text=True,shell=True,capture_output=True).stdout.split("\n")
        try:
            autocorrelation = float(aa_output[0].split()[3][:-1])
            autocorrelation_error = float(aa_output[0].split()[4][:-1])
            autocorr_dict[name] = (autocorrelation,autocorrelation_error)

        except:
            autocorrelation,autocorrelation_error=aa_output[0].split()[3].split("(")
            autocorr_dict[name] = (float(autocorrelation),float(autocorrelation_error[:-1]))
    return autocorr_dict

def compute_with_fsh_jackknife(data,column,bins,system_size,min_b,max_b,path="/home/haaaaron/SUN_twist_python_analysis/modules/fsh_temp/",acc="0.001",thermalization_range=(1000,2000)):
    fsh_list_lines = []
    subprocess.run(f"rm {path}*",shell=True)
    for i,(name,datas) in enumerate(data.items()):
        file_posfix = name.split()[0]
        file_name = f"r_{file_posfix}"
        number_of_lines = len(datas["sum"][thermalization_range[0]:thermalization_range[1]])
        np.savetxt(f"./modules/fsh_temp/{file_name}",datas["sum"][thermalization_range[0]:thermalization_range[1]],fmt='%.18f')
        fsh_list_lines.append(f"{file_posfix} {number_of_lines} {path}{file_name}")
   # print(np.array(fsh_list_lines))
    #np.savetxt(f"./modules/fsh_temp/fsh_list",np.array(fsh_list_lines))

    with open("./modules/fsh_temp/fsh_list", "w") as txt_file:
        txt_file.write(f"ascii {system_size} \n\n")
        for line in fsh_list_lines:
            txt_file.write(line + "\n")
    print("Runnin fsh")
    out = subprocess.run(f"/home/haaaaron/bin/fsh -t100 -d1 -b{min_b}:{acc}:{max_b} -x -j{bins} -J {path}jackknife_output {path}fsh_list ",text=True,shell=True,capture_output=True)
    jackknife,beta,action = np.loadtxt(f"{path}jackknife_output",usecols=(0,2,3),unpack=True)
    jackknife_data = []
    print("Packing data")
    for ind in np.unique(jackknife):
        mask = jackknife == ind
        jackknife_data.append(action[mask])
    return (beta[mask],np.array(jackknife_data), name.split()[1])

def construct_jackknife_functions(data):
    all_jackknife_data = []
    x = []
    twist_coeff = None
    for key in data:
        all_jackknife_data.append([float(x[0]) for x in data[key][1]])
        x.append(float(key.split()[0]))
        twist_coeff = key.split()[1]
 
    return (x,np.array(all_jackknife_data).T,twist_coeff)

def sort_plaquette_dict(plaquette_dict):
    keys = list(plaquette_dict.keys())
    keys = sorted(keys,key = lambda x: float(x.split()[0]))
    return {i : plaquette_dict[i] for i in keys}

def compute_covariance(twist,notwist,thermalization=50):
    cov_dict = {}
    for x in zip(twist,notwist):
        t_sample = twist[x[0]]["sum"][thermalization:]
        nt_sample = notwist[x[1]]["sum"][thermalization:]
        mean_t = np.mean(t_sample)
        mean_nt = np.mean(nt_sample)
        num_sample = len(t_sample)
        cov = np.sum([(val[0]-mean_t)*(val[1]-mean_nt) for val in zip(t_sample,nt_sample)])/num_sample
        cov_dict[x[0]] = cov
    return cov_dict

def print_df_as_markdown_fourier_modes(df):
    df["linear"] = df["linear"].apply(lambda x: r"{:.4g} $\pm$ {:.4g}".format(x[0], x[1]))
    df.rename(columns={"linear": r"Linear fit ($\sigma / T^3$)"}, inplace=True)
    table = df.to_markdown(index=False)
    display(Markdown(table))
    
def list_all_folders(root_folder, has_attribute):
    folder_list = []
    for dirpath, dirnames, _ in os.walk(root_folder):
        #print(dirpath,has_attribute in dirpath)
        if has_attribute in dirpath:
            for dirname in dirnames:
                #print(dirpath,dirname)
                folder_list.append(os.path.join(dirpath, dirname))
    folder_list.sort()
    for i, folder in enumerate(folder_list):
        print(folder,f", index: {i}")
    #print(folder_list)
    return folder_list

def write_surface_tension_dict(data_dict):
    with open("./Results/Surface_tension_smeared/Surface_tension_smeared.txt", "w") as file:
        for key, value in data_dict.items():
            file.write(f"{key}: {value}\n")

def read_surface_tension_dict(file_path):
    data_dict = {}
    with open(file_path, "r") as file:
        for line in file:
            key, value = line.strip().split(": ")
            data_dict[key] = float(value)
    return data_dict
