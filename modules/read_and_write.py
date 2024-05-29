import sys
import os
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
                numbers = [x for x in line.split(",")[:-1]]
                complex_converted = []
                try:
                    for x in numbers:
                        num = x.split()
                        complex_converted.append(complex(float(num[0]),float(num[1])))
                        #complex(num[0])
                    #print(complex_converted)
                except:
                    complex_converted = numbers
                data.append(complex_converted)
                
    data = pd.DataFrame(data[1:],columns=data[0]+["sum"])
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

if __name__ == "__main__":
    print(load_from_folder("./current_output/notwist","plaquette:","real"))
    #load_data_file("../current_output/")
    pass