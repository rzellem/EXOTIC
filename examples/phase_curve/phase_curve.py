import numpy as np
import pandas as pd
from pylightcurve.models.exoplanet_lc import eclipse_mid_time, transit
import matplotlib.pyplot as plt
from exotic.api.elca import phase_curve

if __name__ == "__main__":
    # Read the CSV file into a DataFrame
    df = pd.read_csv('wasp-18b_lightcurve.csv')

    # Extract the relevant columns from the DataFrame
    time_values = [time + 2450000 for time in df['Time (JD)'].tolist()[:5000]]
    relative_flux_values = df['Relative Flux'].tolist()[:5000]
    relative_flux_error_values = df['Relative Flux Error'].tolist()[:5000]

    # Normalize the relative flux values
    relative_flux_array = np.array(relative_flux_values)
    mean_value = np.mean(relative_flux_array)
    normalized_relative_flux_values = relative_flux_array / mean_value
    relative_flux_values = normalized_relative_flux_values.tolist()

    # Convert to numpy arrays
    time_values = np.array(time_values)
    relative_flux_values = np.array(relative_flux_values)
    relative_flux_error_values = np.array(relative_flux_error_values)

    # Define the dictionary of transit parameters
    transit_params = {
        'u0': 2.8779848420193685,
        'u1': -5.214163845369339,
        'u2': 5.10549846519914,
        'u3': -1.7693193857713148,
        'rprs': 0.0969,
        'per': 0.941,
        'ars': 3.48,
        'ecc': 0,
        'inc': 83.5,
        'omega': 265,
        'tmid': 2456740.74189,
        'A': 250e-6,
        'B': 10e-6,
        'fpfs': 0.1  # Assuming a planet-to-star flux ratio of 0.1 (example value)
    }

    # Compute the phase curve model with eclipse model and transit model
    phase_curve_model = phase_curve(time_values, transit_params)

    # Compute phase values
    period = transit_params['per']
    tmid = transit_params['tmid']
    phases = (tmid - time_values / period) % 1 - 0.182

    # Plot the phase curve model and observed data
    plt.figure(figsize=(12, 6))
    plt.plot(phases, phase_curve_model, 'r.', label='Phase Curve Model', linestyle='none')
    plt.scatter(phases, relative_flux_values, c='b', s=2, alpha=0.5, label='Observed Data')
    plt.xlabel('Phase', fontsize=14)
    plt.ylabel('Relative Flux', fontsize=14)
    plt.title('Phase Curve Model for WASP-18 b', fontsize=16)
    plt.grid(True, ls='--', alpha=0.5)
    plt.legend()
    plt.show()