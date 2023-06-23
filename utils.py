import matplotlib.pyplot as plt
from scipy.fftpack import fft, fftfreq
import numpy as np


def plot_obsesrvable_and_drive(x_array, drive_array = None, t_array = None):
    
    fig, ax1 = plt.subplots(figsize=(15,2), dpi = 150)
    color = "steelblue"
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Displacement', color=color)
    ax1.plot(t_array, x_array)
    ax1.tick_params(axis='y', labelcolor=color)

    if drive_array is not None:
        # Instantiate a second y-axis that shares the same x-axis
        ax2 = ax1.twinx()

        color = "r"
        # We already handled the x-label with ax1y
        ax2.set_ylabel('Drive function', color=color)
        ax2.plot(t_array, drive_array, color=color, alpha = 0.2)
        ax2.tick_params(axis='y', labelcolor=color)

        fig.tight_layout()  # Otherwise the right y-label is slightly clipped


def plot_power_spectrum(t, x_array):
    """
    Compute and plot the power spectrum of a signal.
    
    Parameters:
    t: Array of time points.
    x_array: Array of signal values.
    """
    # Compute the Fast Fourier Transform (FFT) of the signal
    fft_vals = fft(x_array)

    # Compute the power spectral density (PSD)
    psd = np.abs(fft_vals) ** 2

    # Get the frequencies corresponding to the values of the PSD
    fft_freq = fftfreq(len(psd), d=(t[1]-t[0]))

    # Only keep the positive frequencies (the spectrum is symmetric)
    i = (fft_freq > 0) & (fft_freq < 100)

    # Plot the power spectrum
    plt.figure(figsize=(10, 2))
    plt.loglog(fft_freq[i], psd[i])
    plt.xlabel('Frequency')
    plt.ylabel('Power Spectral Density')
    # plt.title('Power Spectrum')
    plt.grid(True)
    plt.show()