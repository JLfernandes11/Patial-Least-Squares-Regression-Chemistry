# Partial Least Squares Regression (PLS)

Partial Least Squares Regression is a latent-variable method that simultaneously projects the X (spectral data/absorbance) and Y (concentrations) matrices into a low-dimensional space of latent variables or principal component (PC) in a way that maximize the covariance between them. PLS is especially effective dealing with highly collinear and noisy spectral data obtained with chemical spectroscopy (UV, vis, infrared, etc). 

To use the model, run the "PLS.py" file. Then select the standards absorbances, wavelength/wavenumber values, known standards concentrations and sample absorbances (all must be .csv files!). 

The program will output the raw spectra, the first and second derivative plot, the mean squared error plot for each latent variables or PC by cross-validation, the regression plot for each analyte and the model's performance parameters (R² from calibration, R² from Cross-validation, MSE from calibration and MSE from Cross-validation). 

You can test the program with the files: 
- standards_spectral.csv
- wavelength.csv
- standard_concentrations.csv
- sample_simulation.csv
