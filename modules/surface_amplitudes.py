import numpy as np 
import matplotlib.pyplot as plt

def surface_amplitudes(smooth_surfaces,return_threshold=False,plot_histogram=True, thermalization=100):
    
    surface_amplitudes = {}
    for smearing_level, surface in smooth_surfaces.items():
        print(f"Smearing Level: {smearing_level}")
        surface_amplitudes[smearing_level] = np.array([np.abs(instance[:,2].max() - instance[:,2].min()) for instance in surface[thermalization:]])
    
    if plot_histogram:
        fig, axs = plt.subplots(len(surface_amplitudes), figsize=(10, 10))
        fig.tight_layout(pad=3.0)

        for ax, (level, amplitudes) in zip(axs.flatten(), surface_amplitudes.items()):
            counts, bins, patches = ax.hist(amplitudes, bins=200, edgecolor='black', alpha=0.7)
            max_bin = bins[np.argmax(counts)]
            ax.axvline(max_bin, color='r', linestyle='dashed', linewidth=1)
            ax.set_xlabel('Amplitude')
            ax.set_ylabel('Frequency')
            ax.set_title(f'Histogram of Surface Amplitudes for Smearing Level {level}')
            ax.legend([f'Max Frequency at {max_bin:.2f}'])

        # Hide any unused subplots
        for i in range(len(surface_amplitudes), len(axs.flatten())):
            fig.delaxes(axs.flatten()[i])

        plt.show()
        
    average_surface_amplitudes = {level: (np.mean(amplitudes),np.min(amplitudes), np.max(amplitudes)) for level, amplitudes in surface_amplitudes.items()}
    print(average_surface_amplitudes)
    
    if return_threshold:
        level_indices = {}
        for level, amplitudes in surface_amplitudes.items():
            indices = np.where(amplitudes > return_threshold)[0]
            level_indices[level] = set(indices)
            
        return level_indices