import numpy as np
import pandas as pd
import subprocess

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

def compute_with_aa(data):
    error_dict = {}
    for i,(name,datas) in enumerate(data.items()):
        np.savetxt("./modules/error_temp/temp_file.txt",datas.to_numpy())
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

def compute_with_fsh_jackknife(data,column,bins,system_size,min_b,max_b,path="/home/haaaaron/SUN_twist_python_analysis/modules/fsh_temp/",acc="0.001"):
    fsh_list_lines = []
    subprocess.run(f"rm {path}*",shell=True)

    for i,(name,datas) in enumerate(data.items()):
        file_posfix = name.split()[0]
        file_name = f"r_{file_posfix}"
        number_of_lines = len(datas["sum"])
        np.savetxt(f"./modules/fsh_temp/{file_name}",datas["sum"],fmt='%.18f')
        fsh_list_lines.append(f"{file_posfix} {number_of_lines} {path}{file_name}")
   # print(np.array(fsh_list_lines))
    #np.savetxt(f"./modules/fsh_temp/fsh_list",np.array(fsh_list_lines))

    with open("./modules/fsh_temp/fsh_list", "w") as txt_file:
        txt_file.write(f"ascii {system_size} \n\n")
        for line in fsh_list_lines:
            txt_file.write(line + "\n")
    print("Runnin fsh")
    out = subprocess.run(f"/home/haaaaron/bin/fsh -t400 -d1 -b{min_b}:{acc}:{max_b} -x -j{bins} -J {path}jackknife_output {path}fsh_list ",text=True,shell=True,capture_output=True)
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

    