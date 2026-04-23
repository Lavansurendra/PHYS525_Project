# Importing neccesary libraries for running Fortran/C++ files set up by make command in moose and flare directories
import sys
import os

# importing numpy and matplotlib for calculations and visualization of parralel connection length
import numpy as np
import matplotlib.pyplot as plt

# specifying the location of the .so files for flare and moose code 
    # NOTE: the insert command allows us to specify what order the directories should be searched in when we run the flare import
    #       since we specified an index of 0 this command ensures that the specified directories are searched first
sys.path.insert(0, os.path.expanduser('~/PHYS525/PHYS525_Project/moose/python'))
sys.path.insert(0, os.path.expanduser('~/PHYS525/PHYS525_Project/flare/python'))

# importing functions from the flare and moose directories that we will need to use to run the C++ code of moose and flare
import moose
import flare
from flare import model
from flare.tasks import fieldline_connection
from moose.grids import R3grid


# Specifying the path to the VMEC output file for HSX.
# NOTE: When other computers run this code they will have to change this path
vmec_file_path = "/home/lavan/PHYS525/PHYS525_Project/wout_HSX_ar6.nc"

# taking the VMEC output files found online and generating a continuous 3d mathematical model of the 
# magnetic field
# the flare code run here can calculate the magnetic vector B to a fraction of a millimeter at any point
# within the device using 3D interpolation
try:
    print("Loading HSX mgrid model...")
    model.load("HSX/mgrid")
    print("Successfully loaded MGRID data.")
except Exception as e:
    print(f"Error loading MGRID data: {e}")
    sys.exit(1)

# creating a boundary torus as a stopping condition for the magnetic field lines (placeholder until we 
# can find the real HSX mesh)
# vessel = model.Torosurf()

# # generating 50 evenly spaced starting points along the major radius between 1.35 meters and 1.54 meters 
# # to analyze how connection length changes as we move outward from the core of the plasma to the edge 
R_array = np.linspace(1.35, 1.54, 50)
# # fixing 50 starting points at z=0 and phi=0 so that all starting points are found on the outboard midplane of the plasma  
Z_array = np.array([0.0])
phi_val = 0.0              # A single Phi angle# # organizing the starting points to be provided to the flare code
# start_points = np.column_stack((R_start, Z_start, Phi_start))




# Generate the grid object and save it exactly how Fortran expects
my_grid = R3grid.rzmesh(R_array, Z_array, phi_val)
my_grid.savetxt("grid.dat")
print("Saved starting coordinates to 'grid.dat' using MOOSE R3grid.")

# providing the flare simulation code with the magnetic field, boundary and starting points set up above

# using the fieldline_connection function from the flare code to trace the magnetic field lines starting from the specified points
# NOTE: since a field line that is perfectly confined will be traced forever this sets a maximum limit for which any field line will be traced
      # if a field line travels 10 kilometers without hitting the boundary then we will assume it is confined and stop tracing
fieldline_connection(grid='grid.dat', lcmax=10000.0, output='lc.dat')

# once the C++ engine stops running lc.dat contains trace data which we can plot against our starting coordinates
if os.path.exists('lc.dat'):
    print("Trace complete. Loading results from 'lc.dat'.")
    # NOTE: np.loadtxt loads the data file but skips all the commented lines explaining what the columns contain
    data = np.loadtxt('lc.dat')

    # extracting the connection lengths from the data file
    if data.ndim > 1:
        # from printing the dataset we found that the backward and forward connection lengths are found in the first and second column of the lc.dat file so we will add them together to get the total connection length for each starting point
        L_c = data[:,0] + data[:,1]
        print(L_c)
    else:
        # Failsafe just in case it only reads one row
        L_c = data[0] + data[1]
        print(L_c)

# plotting the distances of each of the field lines starting from our specified starting R coordinates to see where each one collides with the boundary
    plt.figure(figsize=(8, 5))
    plt.plot(R_array, L_c, marker='o', linestyle='-')
    plt.title("HSX Parallel Connection Length (First Wall)")
    plt.xlabel("Starting Radius R (m)")
    plt.ylabel("Connection Length L_c (m)")
    plt.yscale('log') 
    plt.grid(True)

    plt.savefig("HSX_Connection_Length.png")

else:
    print("Error: The C++ engine failed to generate the 'lc.dat' output file.")
