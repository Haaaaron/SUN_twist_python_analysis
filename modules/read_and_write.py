import sys
import os
import io
import pandas as pd
import numpy as np

def load_data_file_complex(file, data_line):
    f = open(file,'r')
    data = []
    beginRead = False
    beta = None
    for line in f:
        if not beginRead:
            if line.startswith('beta'):
                beta = str(line.split()[1])
            if line.startswith('twist_coeff'):
                twist_c = str(line.split()[1])
            if line.startswith('MEASURE start'):
                beginRead = True
        elif line.startswith('MEASURE end'):
            beginRead = False
        else:
            if line.startswith(data_line): 
                line = line.split(':')[1]
                numbers = [x for x in line.split(",")[:]]
                complex_converted = []
                try:
                    for x in numbers:
                        num = x.split()
                        complex_converted.append(complex(float(num[0]),float(num[1])))
                        #complex(num[0])
                except:
                    numbers[-1] = numbers[-1].replace('\n','')
                    
                    complex_converted = [x.replace(' ','') for x in numbers]
                data.append(complex_converted)
    if (data[1][-1] == "sum"):
        data = pd.DataFrame(data[1:],columns=data[0]+["sum"])
    else:
        data = pd.DataFrame(data[1:],columns=data[0])
    name = beta + " " + twist_c
    return name,data

def load_data_file_real(file, data_line):
    f = open(file,'r')
    data = []
    beginRead = False
    beta = None
    for line in f:
        if not beginRead:
            if line.startswith('beta'):
                beta = str(line.split()[1])
            if line.startswith('twist_coeff'):
                twist_c = str(line.split()[1])
            if line.startswith('MEASURE start'):
                beginRead = True
        elif line.startswith('MEASURE end'):
            beginRead = False
        else:
            if line.startswith(data_line):
                
                line = line.split(':')[1]
                #print(line.split())
                data.append(line.split())

    data = pd.DataFrame(data[1:],columns=data[0],dtype=np.float32)
    name = beta + " " + twist_c
    return name,data

def load_from_folder(folder,data_line,real_or_complex):
    data_per_file = {}
    if real_or_complex == "real":
        load_func = load_data_file_real
    else:
        load_func = load_data_file_complex
        
    for root,dirs,files in os.walk(folder):
        for name in files:
            #print(os.path.join(root,name))
            name,data = load_func(os.path.join(root,name), data_line)
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


def read_surface_data(file_path, file_name="surface_smooth"):
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
                all_arrays.append(array)
            # Reset header_found and set new data_start_index
            header_found = True
            data_start_index = idx
        
    # After loop ends, extract the last dataframe section if any
    if header_found and data_start_index is not None:
        array_section = lines[data_start_index:]
        array = np.loadtxt(io.StringIO('\n'.join(array_section)), delimiter=' ', skiprows=3)
        all_arrays.append(array)
    return volume,np.array(all_arrays)

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

if __name__ == "__main__":
    print(load_from_folder("./current_output/notwist","plaquette:","real"))
    #load_data_file("../current_output/")
    pass