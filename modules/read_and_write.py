import sys
import os
import io
import pandas as pd
import numpy as np
import ast

def load_data_file_complex(file, data_line, delim =",", dim=1,column=None,skip_lines=2):
    with open(file, 'r') as f:
        data = []
        beginRead = False
        beta = None
        for line in f:
            ##if not beginRead:
            if line.startswith('beta'):
                beta = str(line.split()[1])
            elif line.startswith('twist_coeff') or line.startswith('twist coeff'):
                twist_c = str(line.split()[1])
            #     if line.startswith('MEASURE start'):
            #         beginRead = True
            # elif line.startswith('MEASURE end'):
            #     beginRead = False
            # else:
            elif line.startswith(data_line): 
                if skip_lines != 0 or ('sum' in line):
                    skip_lines -= 1
                    continue
                line = line.split(':')[1].replace('\n','')
                complex_converted = []
                #try:
                if dim != 1 and not column:   
                    #print(line) 
                    #print(line.replace(",",""))                    
                    #complex_converted = np.fromstring(line.replace(",",""), sep=" ", dtype=np.float64)
                    #complex_converted = complex_converted[::2] + 1j*complex_converted[1::2]
                    #print(complex_converted)
                    #break
                    numbers = line.split(delim)
                    complex_converted = [complex(*map(float, x.split())) for x in numbers]
                elif dim !=1 and column:
                    numbers = line.split(delim)
                    complex_converted = complex(*map(float, numbers[column].split()))
                else:
                    #print(line)
                    line = line.replace(",","")
                    #print(line)
                    complex_converted = complex(*map(float, line.split()))
                #except:
                #    numbers[-1] = numbers[-1].replace('\n', '')
                #    complex_converted = [complex(*map(float, x.split())) for x in numbers]
                data.append(complex_converted)
        data = np.array(data[:-2])
        print("Data array shape:", np.shape(data))
        print("First row of data:", data[0])
        name = beta + " " + twist_c
        return name, data

def load_data_file_real(file, data_line, delim, dim=1, column=None,skip_lines=2,data_type=float):
    f = open(file,'r')
    data = []
    beginRead = False
    beta = None
    line_data = None
    for line in f:
        if line.startswith('beta'):
            beta = str(line.split()[1])
        elif line.startswith('twist_coeff') or line.startswith('twist coeff'):
            twist_c = str(line.split()[1])
        elif line.startswith(data_line):
            if (skip_lines != 0) or ('sum' in line):
                skip_lines -= 1
                continue
            line = line.split(':')[1].replace('\n','')
            if (dim != 1) and (not column):
                line_data = line.split(delim)
                data.append(np.array(line_data, dtype=float))
            elif dim !=1 and column != None:
                if delim != None:
                    line_data = line.split()[column]
                    data.append(np.array(line_data, dtype=float))
                else:
                    line_data = line.split(delim)[column]
                    data.append(np.array(line_data, dtype=float))
            else:
                line_data = line
                data.append(np.array(line_data, dtype=float))

            
            
            # elif len(line_data) == 1:
            #     data.append(float(line_data[0]))
            #else:
            
            #break

    #data = pd.DataFrame(data[1:],columns=data[0],dtype=np.float32)
    data = np.array(data[:-2])
    print("Data array shape:", np.shape(data))
    print("First row of data:", data[0])
    name = beta + " " + twist_c
    return name,data

def load_from_folder(folder,data_line,real_or_complex,delim=None,dim=1,column=None):
    print("Loading data from folder: ",folder)
    data_per_file = {}
    if real_or_complex == "real":
        load_func = load_data_file_real
    else:
        delim = ","
        load_func = load_data_file_complex
        
    for root,dirs,files in os.walk(folder):
        for name in files:
            #print(os.path.join(root,name))
            if ("out" in name) or  ("suN_" in name):
                print(name)
                name,data = load_func(os.path.join(root,name), data_line, delim,dim=dim,column=column)
                #data = data['sum'].to_numpy()
                #columns = list(data.columns.values)
                #columns = columns[-len(columns)//2:] + columns[:-len(columns)//2]
                #data = data.reindex(columns=columns)
                #data = data[50:].mean()
                #data = data[50:]
                #print(data)
                data_per_file[name] = data
    #print(data_per_file['5.081616161616162 0'])
    myKeys = list(data_per_file.keys())
    myKeys.sort()
    data_per_file = {i: data_per_file[i] for i in myKeys}
    return data_per_file

def write_reweight_data(data, file_name):
    data_to_write = np.concatenate(([data[0]],data[1]))
    np.savetxt(file_name,data_to_write)

def read_reweight_data(file_name):
    data = np.loadtxt(file_name)
    #print(data)
    return (data[0],data[1:],file_name.split("_")[-2])


def read_surface_data(file_path, file_name="surface_smooth", samples=(0, None)):
    # Initialize an empty list to store all dataframes read from the file
    all_arrays = []
    
    # Open the file and read line by line
    with open(file_path + "/" + file_name, 'r') as file:
        lines = file.readlines()
        
    # Initialize variables to track the start of each dataframe section
    header_found = False
    data_start_index = None
    volume_count = 0
    start_sample, end_sample = samples
    print("Constructing surface array")
    # Iterate over each line in the file
    for idx, line in enumerate(lines):
        line = line.strip()
        
        if idx == 0:
            volume = np.array(line.split(" ")[1:]).astype(int)
        # Check if the line starts with 'volume:'
        if line.startswith('volume:'):
            volume_count += 1
            if volume_count < start_sample:
                continue
            if end_sample is not None and volume_count > end_sample:
                break
            # If header is found and data_start_index is set (indicating previous dataframe section),
            # extract the dataframe section from data_start_index to current index (exclusive)
            if header_found and data_start_index is not None:
                array_section = lines[data_start_index:idx]
                array = np.loadtxt(io.StringIO('\n'.join(array_section)), delimiter=' ', skiprows=3)
                all_arrays.append(array)
            # Reset header_found and set new data_start_index
            header_found = True
            data_start_index = idx
        
    # After loop ends, extract the last dataframe section if any
    if header_found and data_start_index is not None:
        array_section = lines[data_start_index:]
        array = np.loadtxt(io.StringIO('\n'.join(array_section)), delimiter=' ', skiprows=3)
        all_arrays.append(array)
    return volume, np.array(all_arrays)

def read_fourier_profile(file_path, file_name="fourier_profile_*"):
    # Initialize an empty list to store all dataframes read from the file
    all_arrays = []
    
    # Open the file and read line by line
    with open(file_path + "/" + file_name, 'r') as file:
        lines = file.readlines()
        
    # Initialize variables to track the start of each dataframe section
    header_found = False
    data_start_index = None
    
    # Iterate over each line in the file
    for idx, line in enumerate(lines):
        line = line.strip()
        
        if idx == 0:
            volume = np.array(line.split(" ")[1:]).astype(int)
        # Check if the line starts with 'volume:'
        if line.startswith('volume:'):
            # If header is found and data_start_index is set (indicating previous dataframe section),
            # extract the dataframe section from data_start_index to current index (exclusive)
            if header_found and data_start_index is not None:
                array_section = lines[data_start_index:idx]
                array = np.loadtxt(io.StringIO('\n'.join(array_section)), delimiter=' ', skiprows=3)
                modes = array[:,0]
                freq = array[:,1]
                all_arrays.append(freq)
            # Reset header_found and set new data_start_index
            header_found = True
            data_start_index = idx
    # After loop ends, extract the last dataframe section if any
    if header_found and data_start_index is not None:
        array_section = lines[data_start_index:]
        array = np.loadtxt(io.StringIO('\n'.join(array_section)), delimiter=' ', skiprows=3)
        modes = array[:,0]
        freq = array[:,1]
        all_arrays.append(freq)
    
    return volume,modes,np.array(all_arrays)

def write_surface_tension_dict(data_dict,file_path="./Results/Surface_tension_smeared/Surface_tension_smeared.txt"):
    with open(file_path, "w") as file:
        for key, value in data_dict.items():
            file.write(f"{key}: {value}\n")

def read_surface_tension_dict(file_path="./Results/Surface_tension_smeared/Surface_tension_smeared.txt"):
    data_dict = {}
    with open(file_path, "r") as file:
        for line in file:
            key, value = line.strip().split(": ", 1)
            data_dict[key] = ast.literal_eval(value)
    return data_dict

def update_surface_tension_dict(new_data_dict, file_path="./Results/Surface_tension_smeared/Surface_tension_smeared.txt"):
    # Load the existing data
    existing_data_dict = read_surface_tension_dict(file_path)
    
    # Update the existing data with new data
    existing_data_dict.update(new_data_dict)
    print(existing_data_dict)
    # Write the updated data back to the file
    write_surface_tension_dict(existing_data_dict)

if __name__ == "__main__":
    print(load_from_folder("./current_output/notwist","plaquette:","real"))
    #load_data_file("../current_output/")
    pass