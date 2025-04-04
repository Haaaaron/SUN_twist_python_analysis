import pandas as pd
import ast
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import matplotlib.animation as animation 
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation, FFMpegWriter
from scipy.optimize import curve_fit
from itertools import cycle

FIG_SIZE=(10,15)
plt.rcParams.update({'font.size': 11})
def polar_plot(complex_values,index,title=None):
    magnitude = np.abs(complex_values)
    phase = np.angle(complex_values)
    cmap = plt.cm.jet

    # Normalize magnitudes to map to colormap range
    index = np.array(list(map(int,index)))
    norm = Normalize(vmin=index.min(), vmax=index.max())
    # Plotting
    plt.figure()
    ax = plt.subplot(111, projection='polar')
    for i, (theta, r) in enumerate(zip(phase, magnitude)):
        ax.plot(theta, r, marker='o',color=cmap(norm(i)), label=f'z-index: {index[i]}')
    if title:
        ax.set_title(title)
    ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    #plt.colorbar(sm, ax=ax,label='Magnitude')
    plt.legend(loc='upper left', bbox_to_anchor=(1.1, 1.1))
    plt.show()
    
def surface_in_3d(surface, artistic=False,extra=None):
    # Create meshgrid for x and y values
    x = np.unique(surface[:, 0])
    y = np.unique(surface[:, 1])
    
    z = surface[:, 2].reshape(len(y), len(x))
    x, y = np.meshgrid(x, y)

    # Create a 3D plot
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    z_min, z_max = z.min(), z.max()
    ax.set_zlim(z_min - 20, z_max + 20)

    # Plot the surface with artistic style
    surf = ax.plot_surface(x, y, z, cmap='YlGn', edgecolor='none', alpha=0.8)

    # Plot contours
    #ax.contour(x, y, z, zdir='z', offset=z_min - 20, cmap='magma')
    if artistic:

        if extra:
            z_1 = extra[0][:, 2].reshape(len(y), len(x))
            surf = ax.plot_surface(x, y, z_1-10, cmap='bone', edgecolor='none', alpha=0.8)
            z_2 = extra[1][:, 2].reshape(len(y), len(x))
            surf = ax.plot_surface(x, y, z_2-20, cmap='YlGnBu', edgecolor='none', alpha=0.8) 
        else:
            surf = ax.plot_surface(x, y, z-10, cmap='bone', edgecolor='none', alpha=0.8)
            surf = ax.plot_surface(x, y, z-20, cmap='YlGnBu', edgecolor='none', alpha=0.8) 

        ax.set_axis_off()

    else:
        cbar = fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5)
        cbar.set_label('Z Value')

        # Set labels with artistic font
        ax.set_xlabel('X', fontsize=15, fontweight='bold', color='darkblue')
        ax.set_ylabel('Y', fontsize=15, fontweight='bold', color='darkgreen')
        ax.set_zlabel('Z', fontsize=15, fontweight='bold', color='darkred')

        # Set background color
        ax.set_facecolor('lightgrey')
        

    plt.show()
def animate_surface_in_3d(frames, volume, output_file="../videos/surface_animation.mp4", fps=10, bitrate=1800):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    def update_surface(frame):
        ax.clear()
        ax.set_zlim(-20,volume[2]+20)
        # Extract data for the current frame
        surface = frames[frame]
        x = np.unique(surface[:, 0])
        y = np.unique(surface[:, 1])
        
        z = surface[:, 2].reshape(len(y), len(x))
        x,y = np.meshgrid(x,y)

        # Plot the surface
        ax.plot_surface(x, y, z, cmap='viridis')

        # Set labels
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        ax.set_title(f'Frame {frame + 1}/{len(frames)}')

        return ax
    plt.ioff()
    ani = FuncAnimation(fig, update_surface, frames=len(frames), interval=100)
    # Create animation
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=fps, metadata=dict(artist='Me'), bitrate=bitrate)
    ani.save(output_file, writer=writer, dpi=100)
    
def exponential_func(x, a, b, c):
    return a * np.exp(-b * x) + c

def linear_func(x,a,b):
    return -x*a + b

global_fig = None
markers = cycle(['o', 's', '^', 'D', 'x', 'v', 'p', '*', 'h', '+', '1', '2', '3', '4', '|', '_'])

def compute_fourier_profile(n_2, f_n, volume, errors=None, beta=None, twist=None, fit_range=3, smearing=None, show_plot=True):
    global global_fig,ax1,ax2
    n_2 = np.array(n_2[1:])
    f_n = np.array(f_n[1:])
    y_data = (f_n * n_2 * np.pi**2 * 4) / 36
    
    # Fitting curve
    coeff, cov = curve_fit(linear_func, n_2[:fit_range], y_data[:fit_range], sigma=errors[1:fit_range+1])
    perr = np.sqrt(np.diag(cov))
    
    # Compute y graph and error bands
    n_2_fit = np.concatenate(([0], n_2[:fit_range]))
    y = linear_func(n_2_fit, *coeff)
    y_low = linear_func(n_2_fit, *(coeff - perr))
    y_high = linear_func(n_2_fit, *(coeff + perr))
    
    # Compute error and y_0
    y_0 = linear_func(0, *coeff)
    partials = np.array([0,1])
    sigma_y0 = np.sqrt(partials.T @ cov @ partials)
    if global_fig is None:
        global_fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))
    if show_plot:
        fmt = next(markers)
        color = ax1._get_lines.get_next_color()
        ax1.errorbar(n_2, y_data, yerr=errors[1:], fmt="-"+fmt, label=f"smear = {smearing}", color=color, markersize=3,capsize=5)
        ax1.set_xlabel(r"$n$")
        ax1.set_ylabel(r"$f_n^2 4 \pi^2 n^2/ N_t^2$")
        ax1.set_title(r"Fourier profile V={}, beta={}, twist={}".format(volume, beta, twist))
        legend_text = r"$N_t^2 / f_0^2 4 \pi^2  = \sigma / T^3=$ {:.4g} $\pm$ {:.4g}, smear={}".format(1/y_0, sigma_y0/y_0**2, smearing)
        ax2.plot(n_2_fit,y, "-", label=legend_text, color=color)
        ax2.fill_between(n_2_fit, y_low, y_high, color=color, alpha=0.5)
        ax2.errorbar(n_2[:fit_range], y_data[:fit_range], yerr=errors[1:fit_range+1], fmt=fmt, color=color, markersize=3, capsize=5)
        ax2.set_xlabel(r"$n$")
        ax2.set_ylabel(r"Fit and Error Bands")
        
        global_fig.subplots_adjust(bottom=0.2)
        ax2.legend()
        ax1.legend()
        plt.tight_layout()
    return 1 / y_0, sigma_y0 / y_0**2
    
def compute_fourier_profile_exponential_fit(n_2, f_n, volume, errors=None, beta=None, twist=None, smearing=None, show_plot=True):
    global global_fig, ax1, ax2
    n_2 = np.array(n_2[1:])
    f_n = np.array(f_n[1:])
    y_data = (f_n * n_2 * np.pi**2 * 4) / 36
    
    # Fitting curve
    coeff, cov = curve_fit(exponential_func, n_2, y_data)
    y = exponential_func(n_2, *coeff)
    perr = np.sqrt(np.diag(cov))
    
    # Compute y graph and error bands
    y = exponential_func(n_2, *coeff)
    y_low = exponential_func(n_2, *(coeff - perr))
    y_high = exponential_func(n_2, *(coeff + perr))
    
    # Compute error and y_0 at x_0 = 0
    y_0 = exponential_func(0, *coeff)
    a, b, c = coeff
    partials = np.array([1, 0, 1])
    sigma_y0 = np.sqrt(partials.T @ cov @ partials)
    
    if global_fig is None:
        global_fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))
    
    if show_plot:
        fmt = next(markers)
        color = ax1._get_lines.get_next_color()
        ax1.errorbar(n_2, y_data, yerr=errors[1:], fmt="-"+fmt, label=f"smear = {smearing}", color=color, markersize=3, capsize=5)
        ax1.set_xlabel(r"$n$")
        ax1.set_ylabel(r"$f_n^2 4 \pi^2 n^2/ N_t^2$")
        ax1.set_title(r"Fourier profile V={}, beta={}, twist={}".format(volume, beta, twist))
        legend_text = r"$N_t^2 / f_0^2 4 \pi^2  = \sigma / T^3=$ {:.4g} $\pm$ {:.4g}, smear={}".format(1/y_0, sigma_y0/y_0**2, smearing)
        ax2.plot(n_2, y, "-", label=legend_text, color=color)
        ax2.fill_between(n_2, y_low, y_high, color=color, alpha=0.5)
        ax2.errorbar(n_2, y_data, yerr=errors[1:], fmt=fmt, color=color, markersize=3, capsize=5)
        ax2.set_xlabel(r"$n$")
        ax2.set_ylabel(r"Fit and Error Bands")
        
        global_fig.subplots_adjust(bottom=0.2)
        ax2.legend()
        ax1.legend()
        plt.tight_layout()
    
    return 1 / y_0, sigma_y0 / y_0**2

def save_global_figure(volume, beta,twist, fit_type, z_smear):
    global global_fig
    if global_fig is not None:
        filename = f"/home/haaaaron/SUN_twist_python_analysis/plots/fourier_profile/volume_{volume[0]}-{volume[1]}-{volume[2]}-{volume[3]}_beta_{beta}_twist_{twist}_fit_{fit_type}_zsmear_{z_smear}.pdf"
        global_fig.savefig(filename)
        print(f"Figure saved as {filename}")
    else:
        print("No global figure to save.")
        
def surface_tension_ratios(surface_tensions, twist_first, twist_second):
    # Extract volume from one of the keys
    sample_key = next(iter(surface_tensions))
    volume = sample_key.split("-")[-4:]
    grouped_surface_tensions = {}
    for key, value in surface_tensions.items():
        beta_value = key.split("-")[1]
        twist = key.split("-")[3]
        if beta_value not in grouped_surface_tensions:
            grouped_surface_tensions[beta_value] = {}
        grouped_surface_tensions[beta_value][twist] = value

    num_plots = len(grouped_surface_tensions)
    print(num_plots)
    if num_plots == 1:
        fig, axes = plt.subplots(1, 1, figsize=(9,9))
        axes = [axes]  # Make it iterable
    elif num_plots == 2:
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        axes = np.array([axes])  # Make it 2D
    else:
        fig, axes = plt.subplots(2, (num_plots + 1) // 2, figsize=(14, 5 * ((num_plots + 1) // 2)))
    
    for idx, (beta_value, twists) in enumerate(grouped_surface_tensions.items()):
        print(beta_value, twists)
        
        data_1 = twists[twist_first]
        data_2 = twists[twist_second]
        
        # Create a matrix by dividing each element with each other
        len_twist_first = len(data_1['linear'])
        len_twist_second = len(data_2['linear'])
        matrix = np.zeros(shape=(len_twist_first, len_twist_second))
        data_2_linear = np.array(data_2['linear'])
        data_1_linear = np.array(data_1['linear'])
        matrix = data_2_linear[:, 0][:, np.newaxis] / data_1_linear[:, 0]
        
        # Compute error propagation
        error_2 = data_2_linear[:, 1][:, np.newaxis]
        error_1 = data_1_linear[:, 1]
        error_matrix = matrix * np.sqrt((error_2 / data_2_linear[:, 0][:, np.newaxis])**2 + (error_1 / data_1_linear[:, 0])**2)
        print(error_matrix)
        N = 4
        k = 2
        # Create a heatmap
        center_value = k * (N - k) / (N - 1)
        vmin = center_value - 0.3
        vmax = center_value + 0.3
        ax = axes[idx // 2, idx % 2] if num_plots > 1 else axes[0]
        cax = ax.matshow(matrix, cmap='coolwarm', vmin=vmin, vmax=vmax)
        for (i, j), val in np.ndenumerate(matrix):
            ax.text(j, i, f'{val:.2f}\n±{error_matrix[i, j]:.2f}', ha='center', va='center', color='black')
        #for (i, j), val in np.ndenumerate(matrix):
        #    ax.text(j, i, f'{val:.4f}', ha='center', va='center', color='black')
        fig.colorbar(cax, ax=ax)
        ax.set_xticks(np.arange(len(data_1['smearing'])))
        ax.set_xticklabels(data_1['smearing'])
        ax.set_yticks(np.arange(len(data_2['smearing'])))
        ax.set_yticklabels(data_2['smearing'])
        ax.set_title(r'$\beta=$ {}'.format(beta_value))
        ax.set_xlabel(f'Smearing Level, twist={twist_first}')
        ax.set_ylabel(f'Smearing Level, twist={twist_second}')
    fig.suptitle(f'Surface tension ratio, volume={volume}, sigma_2/sigma_1 = {k*(N-k)}/{N-1} = {k*(N-k)/(N-1):.2f}')
    plt.tight_layout(pad=2.0)
    plt.subplots_adjust(hspace=0.3, wspace=0)  # Adjust space between subplots
    plt.savefig(f"./Results/tension_ratio/volume_{volume}_twist_{twist_first}_{twist_second}.pdf")
    #plt.show()
