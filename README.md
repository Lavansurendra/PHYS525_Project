# PHYS525 Project Part 1: Stellarator Parallel Connection Length & 1D Parallel Diffusion

This repository contains the Python scripts required to simulate 3D magnetic field line parallel connection lengths ($L_c$) for the Helically Symmetric eXperiment (HSX). These $L_c$ values are used to set the length of the magnetic field line over which a diffusion equation describing the propagation of a disturbance in the density of the plasma is solved.

Because this project relies on the **FLARE** physics engine and the **MOOSE** framework (which must be compiled from source), you cannot simply `pip install` the requirements. 

## 1. System Requirements & System Packages
You must be running Linux (or Ubuntu via WSL) with Miniconda installed. 

First, install the underlying C/Fortran compilers and the MPI/Math libraries. Open your terminal and run:
```
sudo apt update
sudo apt install build-essential cmake gfortran pkgconf \
                 libnetcdf-dev libnetcdff-dev \
                 liblapack-dev openmpi-bin libopenmpi-dev
```

## 2. Conda Environment Setup
We use Conda for the base Python environment, but you must install `mpi4py` via `pip` to ensure it compiles against the system's OpenMPI and `gfortran` libraries. **Do not use Conda to install MPI.**

```
# Create and activate the environment
conda create --name PHYS525 python=3.10
conda activate PHYS525

# Install the pre-compiled math libraries via Conda
conda install numpy scipy sympy matplotlib netcdf4

# Install MPI via pip (CRITICAL to avoid Fortran version clashes)
pip install mpi4py
```

## 3. Clone the Sub-Repositories
Make sure you are in the main project folder. You need to pull down the source code for Dr. Frerichs' MOOSE and FLARE libraries. 

```
git clone https://gitlab.com/hfrerichs/moose.git
git clone https://gitlab.com/hfrerichs/flare.git
```

## 4. Compile the MOOSE Framework
MOOSE must be compiled first and "installed" into a local hidden directory so FLARE can find its configuration blueprints.

```
cd moose
mkdir build
cd build

# Configure CMake (Force it to use the Conda Python path and install locally)
cmake -DPython3_EXECUTABLE=$(which python) -DCMAKE_INSTALL_PREFIX=$HOME/.local ..

# Compile and install
make
make install

# Return to the main project directory
cd ../..
```

## 5. Compile the FLARE Physics Engine
Now you must build FLARE and explicitly point it to the local MOOSE installation you just created.

```
cd flare
mkdir build
cd build

# Configure CMake (Point to MOOSE targets in .local)
cmake -DCMAKE_PREFIX_PATH=$HOME/.local -DCMAKE_INSTALL_PREFIX=$HOME/.local -DPython3_EXECUTABLE=$(which python) ..

# Compile the physics engine
make

# Install the Python bindings locally
make install

# Return to the main project directory
cd ../..
```

## 6. The Global Database Structure (Critical for FLARE Tracing)
FLARE does not read raw magnetic grid files via absolute paths in Python. It relies on a strictly defined global `DATABASE` directory in your Linux home folder to load "Machine Models." You must build this exact folder structure to run the HSX simulations.

**Step 1: Create the directory**
Run this command from anywhere in your WSL/Linux terminal:
```
mkdir -p ~/DATABASE/flare/HSX/mgrid
```

**Step 2: Add the MGRID File**
Place the massive 3D magnetic grid file (`mgrid_res2p5cm_180pln 1.nc` or another magnetic grid file if you have one) directly into the `~/DATABASE/flare/HSX/mgrid/` folder. 

*NOTE: The `mgrid_res2p5cm_180pln 1.nc` file was provided to us by Dr. Dieter Boeyaert.*

**Step 3: Create the `.bfield` Configuration**
Inside the `mgrid` folder, create a hidden file named `.bfield`:
```
nano ~/DATABASE/flare/HSX/mgrid/.bfield
```
Paste the following configuration, which specifies the file name and the exact coil current amplitudes (in Amperes) required to generate the QHS state:
```
[DEFAULT]

[equi3d_mgrid]
filename: mgrid_res2p5cm_180pln 1.nc
amplitudes: [-1.0722E+04,-1.0722E+04,-1.0722E+04,-1.0722E+04,-1.0722E+04,-1.0722E+04,0.00,0.00,0.00,0.00,0.00,0.0000E+00]
dtype: 'magnetic_field'
```
*NOTE: the mgrid_res2p5cm_180pln 1.nc and vessel_hsx_flare.txt.txt files and configuration specifications used throughout this code were provided by Dr. Dieter Boeyaert*

*NOTE: ensure the file name matches the file name you have in the corresponding folder*

**Step 4: Add the 3D Vessel Mesh**
Place the vacuum vessel coordinate file (`vessel_hsx_flare.txt` or another coordinate file if you have one) directly into the `~/DATABASE/flare/HSX/mgrid/` folder. 

*NOTE: The `vessel_hsx_flare.txt` file was provided to us by Dr. Dieter Boeyaert.*

Create a `.boundary` configuration file to define the 3D `torosurf` shape:
```
nano ~/DATABASE/flare/HSX/mgrid/.boundary
```
Paste the following text exactly as written:
```
[DEFAULT]

[torosurf]
filename: vessel_hsx_flare.txt
```

## 7. Running the Code

**Calculate Parallel Connection Lengths:**
Because the FLARE C++ wrapper operates via file I/O, this script uses `moose.grids` to generate a `grid.dat` input file for the Fortran engine, executes the field-line trace against the 3D vessel mesh, and outputs data to `lc.dat`.
```
python parallel_connection_length_finder.py
```

## 8. Run the 1D Diffusion Equation Solver
Set the parallel connection length, the desired electron temperature of the plasma in eV, the total amount of time you want the diffusion equation to be solved for, the number of spatial grid points you want the diffusion equation solved for along the magnetic field line with length equal to the parallel connection length, the number of time steps you want the diffusion equation solved for, and the baseline plasma density. Note that if your settings cause the diffusion equation solving method (Forward Time Center Space) to become numerically unstable, an error will be thrown suggesting a new number of timesteps to calculate over to make the solver stable. You can also modify the shape (amplitude, initial position, width) of the density disturbance. The solver will save the calculated density values for each time and spatial position along the field line to a CSV and then animate the propagation of the disturbance. Set the recalculate_data flag to be true if you want to recalculate the density at every point along the field line and at every time and then animate the new data. If you leave the flag set to be false, it will load data from an already saved csv and animate that instead.

```
python Par_Diffusion_eq_solver.py
```

# PHYS525 Project Part 2: Magnetic Field Strength and & 1D Perpendicular Diffusion

This repository also contains the python scripts required to evaluate the magnetic field strength ($|B|$) at a number of positions for the Helically Symmetric eXperiment (HSX). These $|B|$ values are used to determine the Larmor radius which is used to determine the perpendicular thermal diffusivity which is a part of the 1D heat diffusion equation. Then this equation is used to describe the diffusion of a heat disturbance across the magnetic field lines of the stellerator and solved.

The setup for the field strength calculation is the same as it was for the parallel connection length calculation. Please follow the numbered steps 1-6 above.

## 8. Run the 1D Diffusion Equation Solver
Running this solver will also run the Field_Strength_finder function in the field_strength_finder.py file which will automatically evaluate the magnetic field strength at specific points. The points are specified by setting the variable N in the 1D diffusion equation solver. You should also set the desired electron temperature of the plasma in eV, the total amount of time you want the diffusion equation to be solved for, the number of time steps you want the diffusion equation solved for, and the baseline plasma density. Note that if your settings cause the diffusion equation solving method (Forward Time Center Space) to become numerically unstable, an error will be thrown suggesting a new number of timesteps to calculate over to make the solver stable. You can also modify the shape (amplitude, initial position, width) of the initial heat disturbance. The solver will save the calculated heat values for each time and spatial position along the field line to a CSV and then animate the diffusion of the disturbance. Set the recalculate_data flag to be true if you want to recalculate the heat at every point along the field line and at every time and then animate the new data. If you leave the flag set to be false, it will load data from an already saved csv and animate that instead.


```
python Perp_Diffusion_eq_solver.py
```
