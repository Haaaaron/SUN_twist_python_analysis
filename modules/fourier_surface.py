import pandas as pd
import ast
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import matplotlib.animation as animation 
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation, FFMpegWriter
from scipy.optimize import curve_fit
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
    
def surface_in_3d(df):
    pivot_df = df.pivot(index='y', columns='x', values='z')

    # Create meshgrid for x and y values
    x = pivot_df.columns.values
    y = pivot_df.index.values
    x, y = np.meshgrid(x, y)

    # Extract z-values from the pivot table
    z = pivot_df.values

    # Create a 3D plot
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_zlim(0,84)

    # Plot the surface
    ax.plot_surface(x, y, z, cmap='viridis')

    # Set labels
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    
def animate_surface_in_3d(frames, volume, output_file="../videos/surface_animation.mp4", fps=10, bitrate=1800):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    def update_surface(frame):
        ax.clear()
        ax.set_zlim(0,volume[2])
        # Extract data for the current frame
        df = frames[frame]
        pivot_df = df.pivot(index='y', columns='x', values='z')
        x = pivot_df.columns.values
        y = pivot_df.index.values
        x, y = np.meshgrid(x, y)
        z = pivot_df.values

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

def compute_fourier_profile(n_2, f_n, volume, errors=None, beta=None, fit_range=3, smearing=None, show_plot=True):
    
    n_2 = np.array(n_2[1:])
    f_n = np.array(f_n[1:])
    y_data = (f_n * n_2 * np.pi**2 * 4) / 36
    
    # Fitting curve
    coeff, cov = curve_fit(linear_func, n_2[:fit_range], y_data[:fit_range], sigma=errors[1:fit_range+1])
    perr = np.sqrt(np.diag(cov))
    
    # Compute y graph and error bands
    y = linear_func(n_2[:fit_range], *coeff)
    y_low = linear_func(n_2[:fit_range], *(coeff - perr))
    y_high = linear_func(n_2[:fit_range], *(coeff + perr))
    
    # Compute error and y_0
    y_0 = linear_func(0, *coeff)
    partials = np.array([0,1])
    sigma_y0 = np.sqrt(partials.T @ cov @ partials)
    
    if show_plot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.errorbar(n_2, y_data, yerr=errors[1:], fmt='o')
        
        fig.subplots_adjust(bottom=0.2)
        ax.set_xlabel(r"$n$")
        ax.set_ylabel(r"$f_n^2 4 \pi^2 n^2/ N_t^2$")
        ax.set_title(r"Fourier profile V={}, beta={}, smear={}".format(volume, beta, smearing))
        ax.plot(n_2[:fit_range], y, "-")
        ax.fill_between(n_2[:fit_range], y_low, y_high, color='gray', alpha=0.5)
        legend_text = r"$N_t^2 / f_0^2 4 \pi^2  = \sigma / T^3=$ {:.4g} $\pm$ {:.4g}".format(1/y_0, sigma_y0/y_0**2)
        plt.gcf().text(0.15, 0.05, legend_text, fontsize=12, bbox=dict(facecolor='white', alpha=0.5))
        plt.show()
    return 1 / y_0, sigma_y0 / y_0**2
    
def compute_fourier_profile_exponential_fit(n_2, f_n, volume, errors=None, beta=None, smearing=None, show_plot=True):
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
    a,b,c = coeff
    partials = np.array([1,0,1])
    sigma_y0 = np.sqrt(partials.T @ cov @ partials)

    if show_plot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.errorbar(n_2, y_data, yerr=errors[1:], fmt='o')
        fig.subplots_adjust(bottom=0.2)
        ax.set_xlabel(r"$n$")
        ax.set_ylabel(r"$f_n^2 4 \pi^2 n^2/ N_t^2$")
        ax.set_title(r"Fourier profile V={}, beta={}, smear={}".format(volume, beta, smearing))
        ax.plot(n_2, y, "-")
        ax.fill_between(n_2, y_low, y_high, color='gray', alpha=0.5)
        legend_text = r"$N_t^2 / f_0^2 4 \pi^2  = \sigma / T^3=$ {:.4g} $\pm$ {:.4g}".format(1/y_0, sigma_y0/y_0**2)
        plt.gcf().text(0.15, 0.05, legend_text, fontsize=12, bbox=dict(facecolor='white', alpha=0.5))
        plt.show()
    
    return 1 / y_0, sigma_y0 / y_0**2
