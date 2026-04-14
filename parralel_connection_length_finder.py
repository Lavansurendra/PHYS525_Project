import sys
import os
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.expanduser('~/PHYS525/PHYS525_Project/moose/python'))
sys.path.append(os.path.expanduser('~/PHYS525/PHYS525_Project/flare/python'))

import moose
import flare
from flare.model import bfield, boundary
from flare.tasks import fieldline_connection

B_field = bfield.load("wout_HSX_ar6.nc")

vessel = boundary.Torosurf()

R_start = np.linspace(1.2, 1.5, 50)  
Z_start = np.zeros(50)               
Phi_start = np.zeros(50)             

start_points = np.column_stack((R_start, Z_start, Phi_start))

task = fieldline_connection.Task(
    bfield=B_field,
    boundary=vessel,
    points=start_points,
    max_length=10000.0  
)
results = task.execute()

L_c = results.connection_length

plt.plot(R_start, L_c, marker='o', linestyle='-')

plt.figure(figsize=(8, 5))
plt.plot(R_start, L_c, marker='o', linestyle='-')
plt.title("HSX Parallel Connection Length (QHS)")
plt.xlabel("Starting Radius R (m)")
plt.ylabel("Connection Length L_c (m)")
plt.yscale('log') # Usually plotted in log scale since confined lines hit max_length
plt.grid(True)
plt.savefig("HSX_Connection_Length.png")
print("Analysis complete! Saved plot to HSX_Connection_Length.png")
