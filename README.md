# Partial Least Squares Regression (PLS)

Partial Least Squares Regression is a latent-variable method that simultaneously projects the X (spectral data/absorbance) and Y (concentrations) matrices into a low-dimensional space of latent variables in a way that maximize the covariance between them. PLS is especially effective dealing with highly collinear and noisy spectral data obtained with chemical spectroscopy (UV, vis, infrared, etc). 

To use the model, run the "PLS.py" file. Then select the standards absorbances, wavelength/wavenumber values, known standards concentrations and sample absorbances (all must be .csv files!). 



You can test the program with the files: 
- standards_spectral.csv
- wavelength.csv
- standard_concentrations.csv
- sample_simulation.csv
