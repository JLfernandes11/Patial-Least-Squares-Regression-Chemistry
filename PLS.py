
# PLS regression with graphical user interface

from sys import stdout
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from sklearn.cross_decomposition import PLSRegression
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.preprocessing import StandardScaler
import os
import tkinter as tk
from tkinter import filedialog, simpledialog
import customtkinter
from tkinter import messagebox


# Initialize tkinter window
window = tk.Tk()
window.title('PLS Regression Chemical analysis')
window.geometry("1200x600")

# Create a frame for the select file buttons
frame_files = tk.Frame(window)
#Display the frame on the window and put it on nowth-west
frame_files.pack(anchor="nw") 

# Row 1: Standards X CSV
row1 = tk.Frame(frame_files)
row1.pack(anchor="w", pady=2)
tk.Label(row1, text="Select the standards (training) spectral data (absorbance) .CSV:").pack(side="left")
entry_x = tk.Entry(row1, width=50)
entry_x.pack(side="left", padx=5)
# Row 2: Wavelengths CSV
row2 = tk.Frame(frame_files)
row2.pack(anchor="w", pady=2)
tk.Label(row2, text="Select the wavelength or wavenumber values the spectral data were colected .CSV:").pack(side="left")
entry_wl = tk.Entry(row2, width=50)
entry_wl.pack(side="left", padx=5)
# Row 3: Select concentrations CSV
row3 = tk.Frame(frame_files)
row3.pack(anchor="w", pady=2)
tk.Label(row3, text="Select the known concentrations of the standards (training) .CSV:").pack(side="left")
entry_y = tk.Entry(row3, width=50)
entry_y.pack(side="left", padx=5)
# Row 4: Select unknowns X CSV
row4 = tk.Frame(frame_files)
row4.pack(anchor="w", pady=2)
tk.Label(row4, text="Select the spectral data (absorbance) of the sample .CSV:").pack(side="left")
entry_xs = tk.Entry(row4, width=50)
entry_xs.pack(side="left", padx=5)

X_path = wl_path = y_path = xs_path = None 

# Ask for the four CSVs
def openfile1():
    global X_path
    X_path  = filedialog.askopenfilename(title="Select standards X CSV", filetypes=[("CSV","*.csv")])
    if X_path:
        entry_x.delete(0, tk.END)
        entry_x.insert(0, X_path)
def openfile2():
    global wl_path
    wl_path = filedialog.askopenfilename(title="Select wavelengths CSV", filetypes=[("CSV","*.csv")])
    if wl_path:
        entry_wl.delete(0, tk.END)
        entry_wl.insert(0, wl_path)
def openfile3():
    global y_path
    y_path  = filedialog.askopenfilename(title="Select concentrations CSV", filetypes=[("CSV","*.csv")])
    if y_path:
        entry_y.delete(0, tk.END)
        entry_y.insert(0, y_path)
def openfile4():
    global xs_path
    xs_path = filedialog.askopenfilename(title="Select unknowns X CSV", filetypes=[("CSV","*.csv")])
    if xs_path:
        entry_xs.delete(0, tk.END)
        entry_xs.insert(0, xs_path)


# Creating the buttons that use the function
button_openfile1 = tk.Button(row1, text="Browse", command=openfile1)
button_openfile2 = tk.Button(row2, text="Browse", command=openfile2)
button_openfile3 = tk.Button(row3, text="Browse", command=openfile3)
button_openfile4 = tk.Button(row4, text="Browse", command=openfile4)

# Display the buttons on the window
button_openfile1.pack(side="top")
button_openfile2.pack(side="top")
button_openfile3.pack(side="top")
button_openfile4.pack(side="top")

# Creating a button to submit the data for analysis (closes the window)
close_button = tk.Button(window, text="Submit", command=window.destroy)
close_button.pack(side="bottom", pady=5) 

# Creating a text explaining the format of the data

frame_text1 = customtkinter.CTkScrollableFrame(window, width=1150, height=500, fg_color="#1b2d34")
frame_text1.pack()


label_text1 = customtkinter.CTkLabel(frame_text1, text="1. The standards (training) spectral data (absorbance) needs to be inserted as a .CSV file. \n"
                        "The spectral data must have one sample per row (no index column), absorbances as columns and the first row as the header (e.g. abs1, abs2, …). Example: ",
                        bg_color="#1b2d34",
                        fg_color="#1b2d34",
                        text_color="white",
                        justify="left")
label_text1.pack(anchor="nw", padx=2, pady=3) 
table_str1 = ("Abs1\tAbs2\tAbs3\n"
            "0.48\t0.49\t0.51\n"
            "0.49\t0.50\t0.53\n"
            "0.50\t0.52\t0.54\n"
            "0.52\t0.53\t0.55")
label_table1 = customtkinter.CTkLabel(frame_text1,
                    text=table_str1,
                    font=("Courier", 10),    # monospace so columns line up
                    bg_color="#1b2d34",
                    fg_color="#1b2d34",
                    text_color="white",
                    justify="left",
                    anchor="w").pack(fill="x", padx=20)

label_text2 = customtkinter.CTkLabel(frame_text1, text="2. The wavelength or wavenumber the spectral data were colected needs to be inserted as a .CSV file. \n"
                        "Each wavelength or wavenumber must be a row (first row is the header. e.g. wavenumber or wavelength) and only one column. Example:",
                        bg_color="#1b2d34",
                        fg_color="#1b2d34",
                        text_color="white",
                        justify="left")
label_text2.pack(anchor="nw", padx=2, pady=3)
table_str2 = ("Wavelengh\n"
            "2498.53638\n"
            "2494.81472\n"
            "2491.09306\n"
            "2487.3714\n")
label_table2 = customtkinter.CTkLabel(frame_text1,
                    text=table_str2,
                    font=("Courier", 10),    # monospace so columns line up
                    bg_color="#1b2d34",
                    fg_color="#1b2d34",
                    text_color="white",
                    justify="left",
                    anchor="w").pack(fill="x", padx=20)

label_text3 = customtkinter.CTkLabel(frame_text1, text="3. The known concentrations of the standards (training) needs to be inserted as a .CSV file. \n"
                        "The concentrations data must have one sample per row (no index column), concentrations of each analyte as columns and the first row as the header (e.g. analyte1, analyte2, …). \n" 
                        "Example:",
                        bg_color="#1b2d34",
                        fg_color="#1b2d34",
                        text_color="white",
                        justify="left")
label_text3.pack(anchor="nw", padx=2, pady=3)
table_str3 = ("Analyte1\tAnalyte2\tAnalyte3\n"
            "0.48      \t0.49    \t0.51\n"
            "0.49      \t0.50    \t0.53\n"
            "0.50      \t0.52    \t0.54\n"
            "0.52      \t0.53    \t0.55")
label_table3 = customtkinter.CTkLabel(frame_text1,
                    text=table_str3,
                    font=("Courier", 10),    # monospace so columns line up
                    bg_color="#1b2d34",
                    fg_color="#1b2d34",
                    text_color="white",
                    justify="left",
                    anchor="w").pack(fill="x", padx=20)

label_text4 = customtkinter.CTkLabel(frame_text1, text="4. The spectral data (absorbance) of the sample needs to be inserted as a .CSV file. \n"
                        "The spectral data must have one sample per row (no index column), absorbances as columns and the first row as the header (e.g. abs1, abs2, …). Example:",
                        bg_color="#1b2d34",
                        fg_color="#1b2d34",
                        text_color="white",
                        justify="left")
label_text4.pack(anchor="nw", padx=2, pady=3)
table_str4 = ("Abs1\tAbs2\tAbs3\n"
            "0.48\t0.49\t0.51\n"
            "0.49\t0.50\t0.53\n"
            "0.50\t0.52\t0.54\n"
            "0.52\t0.53\t0.55")
label_table4 = customtkinter.CTkLabel(frame_text1,
                    text=table_str4,
                    font=("Courier", 10),    # monospace so columns line up
                    bg_color="#1b2d34",
                    fg_color="#1b2d34",
                    text_color="white",
                    justify="left",
                    anchor="w").pack(fill="x", padx=20)



window.mainloop()

data_X        = pd.read_csv(X_path)
data_wl       = pd.read_csv(wl_path)
data_y        = pd.read_csv(y_path)
data_X_sample = pd.read_csv(xs_path) 

X = data_X.to_numpy()
y = data_y.to_numpy()
wl = data_wl.to_numpy().flatten()
X_sample = data_X_sample.to_numpy()


### Ask pre-processing SG parameters
derivative = simpledialog.askinteger("Derivative",
    "Do you want to apply the first or second derivative to the data? (type 1 or 2):", minvalue=1, maxvalue=2)
window     = simpledialog.askinteger("Window length",
    "Enter the window length to be applied to the Savgol filter (odd integer ≥ 3):", minvalue=3)
polyorder  = simpledialog.askinteger("Polyorder",
    f"Enter the polyorder for the baseline correction (1–{window-1}): \n"
    "(If you chose the second derivative, do not choose polyorder=1 here)"
    , minvalue=1, maxvalue=window-1)



#... Preprocessing (1): first or second derivative
match derivative:
        case 1:
            d1X = savgol_filter(X, window, polyorder, deriv=1)
        case 2:
            d2X = savgol_filter(X, window, polyorder, deriv=2) 
        case default:
            print("Error: Derivative values can only be 1 or 2") 
# for deriv=2, polyorder has to be > 1.
# this is because the Savgol filter devides the data into subsets. 
# Each subset data will have a first-order fit applied and then the seconde derivative. 
# The second derivative of a first-order (Ex: ax + b) is zero.
# So each subset will have second derivative result as zero.

# polyroder=1 and deriv=1 results in a consant but it's fine 
# because each subset will have a constant value that it's different compared to the other subsets. 


# Aplying the filter first or second derivative to the spectral data of unknown samples
match derivative:
        case 1:
            d1X_sample = savgol_filter(X_sample, window, polyorder, deriv=1)
        case 2:
            d2X_sample = savgol_filter(X_sample, window, polyorder, deriv=2)


#... Preprocess (2) Center the standards data by removing the mean
match derivative:
        case 1:
            scaler_X = StandardScaler(with_std=False).fit(d1X) # store the means of the standards spectral data so we can center the sample spectral data with it, so both unknown sample and standards (training) are centered to the same basis
            X_centered = scaler_X.transform(d1X)
        case 2:
            scaler_X = StandardScaler(with_std=False).fit(d2X) # store the means of the standards spectral data so we can center the sample spectral data with it, so both unknown sample and standards (training) are centered to the same basis
            X_centered = scaler_X.transform(d2X)
# Centering the Y matrix
scaler_Y = StandardScaler(with_std=False).fit(y)  # store the means of the standard concentrations Y matrix so we can de-center the result unknown sample determined concentration Y matrix later 
Y_centered = scaler_Y.transform(y) 
# Centering the X_sample matrix
match derivative:
    case 1:
        X_sam_cen = scaler_X.transform(d1X_sample)
    case 2: 
        X_sam_cen = scaler_X.transform(d2X_sample)



#... First derivative to make the graph
d1X_plot = savgol_filter(X, 18, polyorder=1, deriv=1)

#... Second derivative to make the graph
d2X_plot = savgol_filter(X, 18, polyorder=2, deriv=2) 


#... Plot absorbance spectra
plt.figure(figsize=(7,5.5))
with plt.style.context(('default')):
    plt.plot(wl, X.T)
    plt.xlabel(data_wl.columns[0], fontsize=14)
    plt.ylabel('Absorbance (a.u.)', fontsize=14)
plt.show()

#... Plot first derivative
plt.figure(figsize=(7,5.5))
with plt.style.context(('default')):
    plt.plot(wl, d1X_plot.T, linewidth=1)
    plt.xlabel(data_wl.columns[0], fontsize=14)
    plt.ylabel('$1^{st}$ Derivative of Absorbance', fontsize=14)
    plt.show()

#... Plot second derivative
plt.figure(figsize=(7,5.5))
with plt.style.context(('default')):
    plt.plot(wl, d2X_plot.T, linewidth=1)
    plt.xlabel(data_wl.columns[0], fontsize=14)
    plt.ylabel('$2^{nd}$ Derivative of Absorbance', fontsize=14)
    plt.show()


def optimise_pls_cv(X_centered, Y_centered, pc):
 
    '''Run PLS including a variable number of components,
       and calculate MSE'''
 
    pls = PLSRegression(n_components=pc) 

    #... Fit 
    pls.fit(X_centered, Y_centered) 

    #... Calibration
    y_c = pls.predict(X_centered)
 
    #... Cross-validation
    y_cv = cross_val_predict(pls, X_centered, Y_centered, cv=10)

    #... Calculate scores for calibration and cross-validation
    score_c = r2_score(Y_centered, y_c)
    score_cv = r2_score(Y_centered, y_cv)

    #... Calculate mean square error for calibration and cross validation
    mse_c = mean_squared_error(Y_centered, y_c)
    mse_cv = mean_squared_error(Y_centered, y_cv)

    return(y_cv, score_c, score_cv, mse_c, mse_cv)

#... test with 8 components
ncomp_test = 8   # Increase this number if you want more components to be tested
r2r_test = []
r2cv_test = []
mser_test = []
mscv_test = []
xticks = np.arange(1, ncomp_test)
for n_comp in xticks:
    predicted_pc, r2r_pc, r2cv_pc, mser_pc, mscv_pc = optimise_pls_cv(X_centered, Y_centered, pc=n_comp)
    r2r_test.append(r2r_pc)
    r2cv_test.append(r2cv_pc)
    mser_test.append(mser_pc)
    mscv_test.append(mscv_pc)
R2_CAL_CV = np.column_stack((r2r_test, r2cv_test))
MSE_CAL_CV = np.column_stack((mser_test, mscv_test))
minimum_in_MSE_cv = np.argmin(mscv_test)
maximum_in_R2_cv = np.argmax(r2cv_test)




#... Plot the MSE and r²
def plot_metrics(vals, ylabel, loc_best_MSE_R2):
    with plt.style.context('default'):
        plt.plot(xticks, np.array(vals[:,1]), '-o', color='blue', mfc='blue', 
                 markersize=7, label = "Cross Validation")
        plt.plot(xticks[loc_best_MSE_R2], vals[loc_best_MSE_R2,1], 'P', ms=12, mfc='k')
        plt.xlabel('Number of Principal Components', fontsize=14)
        plt.xticks = xticks
        plt.ylabel(ylabel, fontsize=14)
        plt.legend(loc="best")
        plt.show()

"""plot_metrics(MSE_CAL_CV, 'MSE', 'min')
plot_metrics(R2_CAL_CV, '$R^{2}$', 'max')"""
plot_metrics(MSE_CAL_CV, 'MSE', minimum_in_MSE_cv)
plot_metrics(R2_CAL_CV, '$R^{2}$', maximum_in_R2_cv)

#... Calculate and print the position of minimum in MSE
print("\nSuggested number of components:", minimum_in_MSE_cv+1)

# Ask which number of PCs to use
suggested = minimum_in_MSE_cv + 1
pc = simpledialog.askinteger(
    "Select number of PCs",
    f"Suggested number of principal components = {suggested} \n"
    f"Enter desired number of principal components to be used (1–{xticks[-1]}):",
    minvalue=1, maxvalue=xticks[-1])

# Calculate PLS with the desired number of pcs 
predicted, r2r, r2cv, mser, mscv = optimise_pls_cv(X_centered, Y_centered, pc)
print('R² from calibration: %5.3f'  % r2r)
print('R² from Cross-validation: %5.3f'  % r2cv)
print('MSE from calibration: %5.3f' % mser)
print('MSE from Cross-validation: %5.3f' % mscv, "\n")


#... Regression plot

for i in range(Y_centered.shape[1]):
    # Extract current column from both actual and predicted matrices
    y_col = Y_centered[:, i]
    predicted_col = predicted[:, i]
    
    # Perform linear regression for this column
    z = np.polyfit(y_col, predicted_col, 1)
    regression_line = z[1] + z[0] * y_col

    # Calculate R²
    r2 = r2_score(y_col, predicted_col)
    
    # Create regression plot
    with plt.style.context(('default')):
        fig, ax = plt.subplots()
        
        # Scatter plot of actual vs predicted
        ax.scatter(y_col, predicted_col, s=150, c='red', 
                   edgecolors='red', marker="*")
        
        # Regression line
        ax.plot(y_col, regression_line, c='blue', linewidth=2, 
                label=f'y = {z[0]:.4f}x + {z[1]:.4f}')
        
        # Labels and title
        plt.xlabel(f'Measured {data_y.columns[i]} Concentration (mg/mL)', fontsize=13)
        plt.ylabel(f'Predicted {data_y.columns[i]} Concentration (mg/mL)', fontsize=13)
        plt.title(f'Regression for {data_y.columns[i]}', fontsize=15)
        
        # Add R² to plot
        ax.text(0.01, 1.04, f'$R^2 = {r2:.4f}$', 
                transform=ax.transAxes, fontsize=12,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8))
        
        plt.show()




# Calculating the concentrations for unknown samples

def pls_sample(X_sam_cen, Y_centered, pc, X_centered, scaler_Y):

    # PLS witht the selected number of pc
    pls = PLSRegression(n_components=pc) 

    # Fit
    pls.fit(X_centered, Y_centered)

    # Getting the concentrations of the the unknown sample
    y_s = pls.predict(X_sam_cen)

    # De-centering the concentration of the unknown sample
    y_sample = scaler_Y.inverse_transform(y_s)

    return(y_sample) 

# Getting the result (y_sample) from the function
result = pls_sample(X_sam_cen, Y_centered, pc, X_centered, scaler_Y) 

# Making it pandas dataframe and inserting the columns label
result_df = pd.DataFrame(result, columns=data_y.columns) 

# Get current directory (same as this .py)
dir = os.path.dirname(os.path.abspath(__file__))

# Create output path in the same directory
output_path = os.path.join(dir, 'predicted_concentrations_pls.csv') 

# Saving results as .csv file
result_df.to_csv(output_path, index=False)


# Showing the statistical info and where the precdiced concentrations were saved 

window2 = tk.Tk()
window2.title('All done!')
window2.geometry("700x500")
window2.config(bg="black") 

frame2 = tk.Frame(window2, bg="#1b2d34")
frame2.pack(fill="both", expand=True, pady=60)

msg = (
    f"Your PLS run is complete.\n\n"
    f"Predicted concentrations saved to:\n{output_path}\n\n\n"
    f"R² from calibration:        {r2r:.3f}\n"
    f"R² from Cross-validation:   {r2cv:.3f}\n"
    f"MSE from calibration:       {mser:.3f}\n"
    f"MSE from Cross-validation:  {mscv:.3f}"
)

tk.Label(frame2, text=msg, bg="#1b2d34", fg="white", font=("Arial", 12)).pack(pady=20)

bottom_frame = tk.Frame(frame2, bg="#1b2d34")
bottom_frame.pack(side="bottom", pady=20)

close_button2 = tk.Button(bottom_frame, text="Done", command=window2.destroy)
close_button2.pack(pady=4)


window2.mainloop()