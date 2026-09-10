# Import python packages
# Only numpy and netCDF4 are absolutely required for reading and manipulating the H1 dataset
# matplotlib required for visualizations
import gc
import os
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt

# Open netcdf file and create netcdf object
# Change data_path variable to data of interest
data_path = './ncdata/2000.nc'
data = Dataset(data_path)

print('Printing a netcdf object shows the description, dimensions, and variables')
print(data)

# Load the actual data for the month of march
# Note that in python the slicing array syntax is necessary
# to load the data to memory
march_data = data['march'][:,:,:]

# Now that data is loaded to memory, close netcdf file
data.close()

# Filter out gridcells with WS 7
# Wherever a gridcell is equal to 7, make gridcell 1, otherwise make 0
ws = 7
ws_march_data = np.where(march_data==ws,1,0)

# Plot for single 3-hour period
# Two plots, one (360,180) the other rotated counter-clockwise 90 degrees
print('Plotting WS 7 for a single 3 hour period')
print('Right plot is rotated counter-clockwise 90 degrees')
print('NOTE: must close plot window for program to continue')
fig, axes = plt.subplots(ncols=2)
axes[0].imshow(ws_march_data[0,:,:])
axes[1].imshow(np.rot90(ws_march_data[0,:,:]))
plt.show()

print('Now making same plot for entire year')
# Create empty list to hold each variable, or 3D data matrix
data_list = []

# Open the file and create a netcdf object called data
# Use the 'with' command to have python automatically handle
# file opening and closing
with Dataset(data_path) as data:
    # Loop over all variables, print the variables, and load their data
    for var in data.variables:
        print(var)
        # Skip over WS histograms
        if var == 'ws':
            continue
        data_list.append(data[var][:,:,:])
        
# Choose the WS to filter out and create an empty array
# to continually add the results of each month to
ws = 7
year_plot = np.zeros((360,180))
for data in data_list:
    # Shape of each array will differ due to number of days in each month
    # Every array should only have values of -127 or 1-11
    print(data.shape, np.unique(data))
    wsdata = np.where(data==ws,1,0)
    # Add the number of counts each gridcell was assigned to the WS
    # Do this by adding along the 'temporal' or first dimension to create 2D map
    year_plot += np.sum(wsdata,axis=0)

fig, axes = plt.subplots(ncols=2)
axes[0].imshow(year_plot)
axes[1].imshow(np.rot90(year_plot))
plt.show()

# Next, load Weather States (WS) and plot
with Dataset(data_path) as data:
    ws = data['ws'][:,:]

# Due to how raw data was processed, in order to properly
# plot WS centroids, they must be reshaped from vector of length (42)
# to a small matrix of shape (6,7) then transposed 
# (the .T after reshape below does the transpose)
fig, axes = plt.subplots(nrows=4, ncols=3)
for i, ax in enumerate(axes.flatten()):
    # If statement for removing unwanted axes/empty subplots from showing
    if i>=10:
        ax.set_axis_off()
        continue
    ax.imshow(ws[i,:].reshape(6,7).T)
