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