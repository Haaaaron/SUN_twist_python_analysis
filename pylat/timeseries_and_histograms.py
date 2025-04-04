import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import matplotlib.cm as cm

from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import make_interp_spline
from matplotlib_inline.backend_inline import set_matplotlib_formats
from scipy.optimize import curve_fit

set_matplotlib_formats('svg')


def create_histogram_animation(data_dict, filename, temp_min, temp_max):
    fig, ax = plt.subplots()

    # Filter the data_dict based on the temperature range
    filtered_data_dict = {k: v for k, v in data_dict.items() if temp_min <= float(k.split(" ")[0]) <= temp_max}

    # Calculate global min and max values for the x-axis using numpy operations

    all_data = np.concatenate([data[1000:] for data in filtered_data_dict.values()])
    global_min = np.min(all_data)
    global_max = np.max(all_data)

    # Function to update the histogram for each frame
    def update_hist(num, data, bins, ax):
        ax.clear()
        ax.hist(data[num][1000:], bins=bins, color='blue', alpha=0.7)
        #ax.set_xlim(global_min, global_max)
        ax.set_title(f"Temperature: {list(filtered_data_dict.keys())[num]}")
        ax.set_xlabel('Average Plaquette Action')
        ax.set_ylabel('Frequency')

    # Extract the data and number of bins
    data = list(filtered_data_dict.values())
    bins = 100

    # Create the animation
    ani = animation.FuncAnimation(fig, update_hist, frames=len(data), fargs=(data, bins, ax), repeat=False)

    # Save the animation
    ani.save(filename, writer='ffmpeg')

    # Close the plot to avoid displaying it in the notebook
    plt.close(fig)
    
def plot_3d_histogram_curves(data_dict, temp_min, temp_max,term=100, column_index=-1):
    fig = plt.figure(figsize=(12, 8))  # Adjust the figure size here
    ax = fig.add_subplot(111, projection='3d')

    # Filter the data_dict based on the temperature range
    filtered_data_dict = {k: v for k, v in data_dict.items() if temp_min <= float(k.split(" ")[0]) <= temp_max}

    # Create a colormap
    cmap = cm.get_cmap('plasma')
    norm = plt.Normalize(temp_min, temp_max)

    # Create the 3D plot
    for temp, data in reversed(filtered_data_dict.items()):
        if data.ndim == 2:
            data = data[:, column_index]
        if np.iscomplexobj(data):
            data = np.abs(data[term:])  # Compute the absolute values if data is complex
        else:
            data = data[term:]  # Use the data as is if it's not complex
        hist, bin_edges = np.histogram(data, bins=200)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        temp_value = float(temp.split(" ")[0])
        color = cmap(norm(temp_value))
        zorder = temp_value  # Use temperature value as zorder to preserve depth perception
        
        # Create a smooth curve for the histogram
        bin_centers_smooth = np.linspace(bin_centers.min(), bin_centers.max(), 300)
        hist_smooth = make_interp_spline(bin_centers, hist)(bin_centers_smooth)
        ax.plot(bin_centers_smooth, [temp_value] * len(bin_centers_smooth), hist_smooth, label=temp, color=color, solid_capstyle='round', zorder=zorder)
        ax.set_xlabel('Value')
        ax.set_ylabel('Temperature')
        ax.set_zlabel('Frequency')
        temp_values = [float(k.split(" ")[0]) for k in filtered_data_dict.keys()]
        ax.set_title(f'3D Histogram Curves (Temp Range: {max(temp_values) - min(temp_values):.4f})')
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.4f}'))
    fig.legend()
    plt.show()
    
def plot_3d_timeseries_curves(data_dict, temp_min, temp_max, column_index=-1):
    fig = plt.figure(figsize=(12, 8))  # Adjust the figure size here
    ax = fig.add_subplot(111, projection='3d')

    # Filter the data_dict based on the temperature range
    filtered_data_dict = {k: v for k, v in data_dict.items() if temp_min <= float(k.split(" ")[0]) <= temp_max}

    # Create a colormap
    cmap = cm.get_cmap('viridis')
    norm = plt.Normalize(temp_min, temp_max)

    # Create the 3D plot
    for temp, data in reversed(filtered_data_dict.items()):  # Reverse the order of plotting
        if len(np.shape(data)) != 1:
            if np.iscomplexobj(data):
                data = np.abs(data[1000:,column_index])  # Compute the absolute values if data is complex
            else:
                data = data[1000:,column_index]  # Use the data as is if it's not complex
        else:
            if np.iscomplexobj(data):
                data = np.abs(data[1000:])  # Compute the absolute values if data is complex
            else:
                data = data[1000:]  # Use the data as is if it's not complex
        #data = np.squeeze(data)  # Remove any extra dimensions
        temp_value = float(temp.split(" ")[0])
        color = cmap(norm(temp_value))
        zorder = temp_value  # Use temperature value as zorder to preserve depth perception

        # Create the timeseries plot
        timeseries_data = data
        time_points = np.arange(len(timeseries_data))
        ax.plot(time_points, [temp_value] * len(timeseries_data), timeseries_data, color=color, solid_capstyle='round', zorder=zorder, linewidth=0.04)
        ax.set_xlabel('Time')
        ax.set_ylabel('Temperature')
        ax.set_zlabel('Value')
        temp_values = [float(k.split(" ")[0]) for k in filtered_data_dict.keys()]
        ax.set_title(f'3D Timeseries Curves (Temp Range: {max(temp_values) - min(temp_values):.4f})')
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.4f}'))

    plt.show()
    
def plot_time_series_histograms(dicts, temp_min, temp_max, hist=True, log_scale=False, column=-1, thermalization=1000):
    # Filter the keys based on the temperature range
    filtered_keys = [key for key in dicts[0][0].keys() if temp_min <= float(key.split()[0]) <= temp_max]

    # Create subplots
    if len(filtered_keys) == 1:
        fig, axs = plt.subplots(1, len(dicts), figsize=(12, 4))
        axs = [axs]  # Wrap in a list to keep the indexing consistent
    else:
        fig, axs = plt.subplots(len(filtered_keys), len(dicts), figsize=(12, 2 * len(filtered_keys)))

    # Define a list of colors for each subplot
    colors = ['blue', 'green', 'red', 'purple', 'orange', 'brown', 'pink', 'gray', 'olive', 'cyan']

    for i, key in enumerate(filtered_keys):
        for j, (data_dict,name) in enumerate(dicts):
            data = data_dict[key]
            if np.iscomplexobj(data):
                data = np.abs(data[thermalization:])
            if len(data.shape) != 1:
                data = data[thermalization:, column]
            else:
                data = data[thermalization:]
            color = colors[j % len(colors)]  # Cycle through colors if there are more subplots than colors
            
            if hist:
                axs[i][j].hist(data, bins=100, color=color, alpha=0.7)
                axs[i][j].set_title(f"{name} at {key}")
                axs[i][j].set_ylabel('Frequency')
                axs[i][j].set_xlabel('Value')
                if log_scale:
                    axs[i][j].set_yscale('symlog')
            else:
                time_points = np.arange(len(data))
                axs[i][j].plot(time_points, data, color=color, alpha=0.7, linewidth=0.5)
                axs[i][j].set_title(f"{name} at {key}")
                axs[i][j].set_ylabel('Value')
                axs[i][j].set_xlabel('Time')

            # Add more ticks in x range
            #axs[i][j].xaxis.set_major_locator(plt.MaxNLocator(nbins=10))
            # Add a legend to the figure

    plt.tight_layout()
    plt.show()
    
def plot_time_series_derivative(dicts, temp_min, temp_max, chosen_dict=0, log_scale=False, column=-1, thermalization=1000):
    # Filter the keys based on the temperature range
    filtered_keys = [key for key in dicts[0][0].keys() if temp_min <= float(key.split()[0]) <= temp_max]

    # Create subplots
    if len(filtered_keys) == 1:
        fig, axs = plt.subplots(1, len(dicts), figsize=(12, 4))
        axs = [axs]  # Wrap in a list to keep the indexing consistent
    else:
        fig, axs = plt.subplots(len(filtered_keys), len(dicts)+1, figsize=(12, 2 * len(filtered_keys)))

    # Define a list of colors for each subplot
    colors = ['blue', 'green', 'red', 'purple', 'orange', 'brown', 'pink', 'gray', 'olive', 'cyan']

    for i, key in enumerate(filtered_keys):
        data_dict, name = dicts[chosen_dict]
        data = data_dict[key]
        if np.iscomplexobj(data):
            data = np.abs(data[thermalization:])
        if len(data.shape) != 1:
            data = data[thermalization:, column]
        else:
            data = data[thermalization:]

        time_points = np.arange(len(data))
        # Define the sigmoid function
        def sigmoid(x, L, x0, k, b):
            return L / (1 + np.exp(-k * (x - x0))) + b

        # Fit the sigmoid function to the data
        popt, _ = curve_fit(sigmoid, time_points, data, maxfev=10000)
        sigmoid_data = sigmoid(time_points, *popt)
        derivative_data = np.gradient(sigmoid_data, time_points)
        second_derivative_data = np.gradient(derivative_data, time_points)
        axs[i][3].plot(time_points, second_derivative_data, color='purple', alpha=0.7, linewidth=0.5, label='Second Derivative')
        axs[i][0].plot(time_points, data, color='blue', alpha=0.3, linewidth=0.5, label='Data')
        axs[i][0].plot(time_points, sigmoid_data, color='red', alpha=0.7, linewidth=0.5, label='Sigmoid Fit')
        distance_from_sigmoid = np.abs(data - sigmoid_data)
        axs[i][1].plot(time_points, distance_from_sigmoid, color='blue', alpha=0.7, linewidth=0.5, label='Distance from Sigmoid')
        axs[i][2].plot(time_points, derivative_data, color='green', alpha=0.7, linewidth=0.5, label='Derivative')
        axs[i][0].set_title(f"Derivative of {name} at {key}")
        axs[i][0].set_ylabel('Value')
        axs[i][0].set_xlabel('Time')
        axs[i][0].legend()

            # Add more ticks in x range
            #axs[i][j].xaxis.set_major_locator(plt.MaxNLocator(nbins=10))
            # Add a legend to the figure

    plt.tight_layout()
    plt.show()

def plot_polyakov_heatmaps(polyakov_data, columns=1, polar=False):
    # Calculate the number of rows needed
    rows = int(np.ceil(len(polyakov_data.keys()) / columns))
    
    # Create a figure with subplots
    fig, axs = plt.subplots(rows, columns, figsize=(8 * columns, 4 * rows), subplot_kw={'projection': 'polar'} if polar else {})
    axs = axs.flatten()  # Flatten the array of axes for easy iteration
    
    # Plot each set of complex numbers in a subplot
    for i, (key, values) in enumerate(polyakov_data.items()):
        real_parts = np.real(values)
        imag_parts = np.imag(values)
        
        if polar:
            # Convert to polar coordinates
            angles = np.angle(values)
            radii = np.abs(values)
            
            # Create a 2D histogram (heatmap) in polar coordinates
            heatmap, angle_edges, radius_edges = np.histogram2d(angles, radii, bins=50)
            
            # Plot the heatmap
            angle_centers = (angle_edges[:-1] + angle_edges[1:]) / 2
            radius_centers = (radius_edges[:-1] + radius_edges[1:]) / 2
            angle_grid, radius_grid = np.meshgrid(angle_centers, radius_centers)
            im = axs[i].pcolormesh(angle_grid, radius_grid, heatmap.T, shading='auto')
            axs[i].set_title(f'Polar Heatmap for {key}')
            fig.colorbar(im, ax=axs[i])
        else:
            # Create a 2D histogram (heatmap)
            heatmap, xedges, yedges = np.histogram2d(real_parts, imag_parts, bins=50)
            
            # Plot the heatmap
            im = axs[i].imshow(heatmap.T, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], origin='lower', aspect='auto')
            axs[i].set_title(f'Heatmap for {key}')
            axs[i].set_xlabel('Real Part')
            axs[i].set_ylabel('Imaginary Part')
            fig.colorbar(im, ax=axs[i])
    
    # Hide any unused subplots
    for j in range(i + 1, len(axs)):
        fig.delaxes(axs[j])
    
    plt.tight_layout()
    plt.show()