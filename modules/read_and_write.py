import sys
import os
import pandas as pd

def load_data_file(file, data_line):
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
                #If data is complex then do a little error handling
                try:
                    numbers = [x for x in line.split(",")[:-1]]
                    #print(numbers)
                    complex_converted = []
                    for x in numbers:
                        num = x.split()
                        complex_converted.append(complex(float(num[0]),float(num[1])))
                        #complex(num[0])
                    #print(complex_converted)
                    data.append(complex_converted)
                except:
                    data.append(line.split())

    data = pd.DataFrame(data[1:],columns=data[0])
    name = beta + " " + twist_c
    return name,data

def load_from_folder(folder,data_line):
    data_per_file = {}
    
    for root,dirs,files in os.walk(folder):
        for name in files:
            #print(os.path.join(root,name))
            name,data = load_data_file(os.path.join(root,name), data_line)
            #data = data['sum'].to_numpy()
            #columns = list(data.columns.values)
            #columns = columns[-len(columns)//2:] + columns[:-len(columns)//2]
            #data = data.reindex(columns=columns)
            #data = data[50:].mean()
            #data = data[50:]
            #print(data)
            data_per_file[name] = data
    #print(data_per_file['5.081616161616162 0'])
    
    return data_per_file

if __name__ == "__main__":
    print(load_from_folder("./current_output/notwist","polyakov:"))
    #load_data_file("../current_output/")
    pass