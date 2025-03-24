import pandas as pd
import time
import os
import sys
from matplotlib.animation import FuncAnimation
import argparse
import subprocess
import matplotlib.pyplot as plt
import numpy as np
plt.style.use('dark_background')

def rsync_file(remote_path="lumi:/projappl/project_462000781/HILA/applications/suN_gauge_twist/24-24-36-6-output-twist-0/out.txt"):
    local_path = '/home/haaaaron/SUN_twist_python_analysis/data/tmp/out.txt'
    rsync_command = ['rsync', '-avz', remote_path, local_path]
    subprocess.run(rsync_command)


def read_data(file_path, label='Order parameter:', skip_lines=0):
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            if label in line:
                if skip_lines > 0:
                    skip_lines -= 1
                    continue
                val = line.split(":")[-1].split(",")[-1]
                value = float(val)
                data.append(value)
    # Drop the first 1000 points for thermalization
    return data[100:]

def read_data_polyakov(file_path, label='polyakov:', skip_lines=0):
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            if label in line:
                if skip_lines > 0:
                    skip_lines -= 1
                    continue
                parts = line.split(":")[-1].split()
                real_part = float(parts[0])
                imag_part = float(parts[1])
                data.append(complex(real_part, imag_part))
    # Drop the first 1000 points for thermalization
    return np.abs(np.array(data[100:]))

def plot_ascii(data, height, width):
    os.system('clear')
    max_value = max(data)
    min_value = min(data)
    scale = (max_value - min_value) / height
    x_scale = len(data) / width

    # Print x-axis labels
    x_labels = ' ' * 10  # Offset for y-axis labels
    for j in range(width):
        if j % (width // 10) == 0:  # Print 10 labels evenly spaced
            x_labels += f'{int(j * x_scale):4d}'
        else:
            x_labels += ' '

    for i in range(height, -1, -1):
        line = f'{min_value + i * scale:8.5f} '  # Add yscale values on the first column
        for j in range(width):
            index = int(j * len(data) / width)
            value = data[index]
            if value >= min_value + i * scale and value < min_value + (i + 1) * scale:
                if j > 0 and data[int((j - 1) * len(data) / width)] < value:
                    line += '/'
                elif j > 0 and data[int((j - 1) * len(data) / width)] > value:
                    line += '\\'
                else:
                    line += '*'
            else:
                line += ' '
        print(line)

def animate(i, label,file, skip_lines, remote_file):
    if remote_file:
        rsync_file(remote_path=remote_file)
    if label=='polyakov:' or label=='polyakov muca:':
        data = read_data_polyakov(file, label=label,skip_lines=skip_lines)
    else:
        data = read_data(file, label,skip_lines)
        
    plt.subplot(2, 1, 1)
    plt.plot(data, color='gray')
    plt.xlabel('Index')
    plt.ylabel('Order Parameter')
    plt.title('Order Parameter vs Index')
    plt.grid(True)

    plt.subplot(2, 1, 2)
    plt.hist(data, bins=200, color='gray')
    plt.xlabel('Order Parameter')
    plt.ylabel('Frequency')
    plt.title('Histogram of Order Parameter')
    plt.grid(True)

def main():
    parser = argparse.ArgumentParser(description='Choose the plotting method.')
    parser.add_argument('--method', choices=['ascii', 'matplotlib'], required=True, help='Plotting method: ascii or matplotlib')
    parser.add_argument('--label', type=str, default='Order parameter:', help='Label to search for in the data file')
    parser.add_argument('--weight_file', type=str, default='/home/haaaaron/SUN_twist_python_analysis/data/tmp/out.txt', help='Path to the weight file')
    parser.add_argument('--skip_lines', type=int, default=100, help='Number of lines to skip for thermalization')
    parser.add_argument('--remote_file', type=str, default=None, help='Full remote file path to rsync before processing')
    args = parser.parse_args()


    weight_file = args.weight_file

    height = 40   # Height for higher resolution
    width = 180  # Width for higher resolution

    if args.method == 'matplotlib':
        ani = FuncAnimation(plt.gcf(), animate, fargs=(args.label,weight_file,args.skip_lines,args.remote_file,), interval=100)
        plt.show()
    elif args.method == 'ascii':
        while True:
            data = read_data(weight_file, args.label)
            print(data)
            plot_ascii(data, height, width)
            time.sleep(1)

if __name__ == "__main__":
    main()
