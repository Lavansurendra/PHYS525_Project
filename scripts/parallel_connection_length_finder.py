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
vessel = model.Torosurf()

# generating 50 evenly spaced starting points along the major radius between 1.2 meters and 1.5 meters 
# to analyze how connection length changes as we move outward from the core of the plasma to the edge 
R_start = np.linspace(1.2, 1.5, 50)
# fixing 50 starting points at z=0 and phi=0 so that all starting points are found on the outboard midplane of the plasma  
Z_start = np.zeros(50)               
Phi_start = np.zeros(50)             
# organizing the starting points to be provided to the flare code
start_points = np.column_stack((R_start, Z_start, Phi_start))

# providing the flare simulation code with the magnetic field, boundary and starting points set up above
task = fieldline_connection.Task(
    boundary=vessel,
    points=start_points,

    # since a field line that is perfectly confined will be traced forever this sets a maximum limit for which any field line will be traced
    # if a field line travels 10 kilometers without hitting the boundary then we will assume it is confined and stop tracing
    max_length=10000.0  
)

# runs the C++ code within the flare module to simulate the field lines from each of the starting points
results = task.execute()

# once the C++ engine stops running results.connection_length contains an array of 50 numbers which we can plot against our starting coordinates
L_c = results.connection_length

# plotting the distances of each of the field lines starting from our specified starting R coordinates to see where each one collides with the boundary
plt.plot(R_start, L_c, marker='o', linestyle='-')

plt.figure(figsize=(8, 5))
plt.plot(R_start, L_c, marker='o', linestyle='-')
plt.title("HSX Parallel Connection Length (QHS)")
plt.xlabel("Starting Radius R (m)")
plt.ylabel("Connection Length L_c (m)")
plt.yscale('log') # Usually plotted in log scale since confined lines hit max_length
plt.grid(True)
plt.savefig("HSX_Connection_Length.png")
