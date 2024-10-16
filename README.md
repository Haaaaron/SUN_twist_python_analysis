# SUN Analysis

Collection of data analysis methods for analysing the phase transitions and subsequent surface in pure SUN gauge theory. The studied boundary is produced by a "twist" operation imposed on the SUN gauge field. A more simple analog can be viewed for the [Ising model](https://haarti.net/blog/).

An example video of the surface in pure $\text{SU}(N)$ gauge theory:

![Surface Animation](/videos/example/surface_animation.mp4)

This collection of methods and analysis tools is in no means a consice package and simply ment to document the process of the user.

## Dependencies
| Dependencies|
|-------------|
| jupyter     |
| numpy       |
| scipy       |
| matplotlib  |
| fsh (included) |
| aa (included) |

## Notebooks

Most of the data manipulation is done in python notebooks where the custom data format is loaded into a notebook from the /data folder. After this the data is post processed and the outputted by various plots. 

| Filename | Descritpion |
| --- | --- |
| SUN_compute_fourier_modes.ipynb | Measures and plots the surface tension based on fourier modes |
| SUN_plot_surface.ipynb | Animate polyakov field in 3D |
| SUN_polyakov.ipynb | Methods to plot polyakov lines |
| SUN_wilson_action.ipynb | Computes surface tension by integrating wilson action |

## modules

The modules folder containts useful python functions. For example utility and problem specific analysis methods.