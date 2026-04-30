# PHYS525 Project: Stellarator Parallel Connection Length & 1D Radial Diffusion

This repository contains the Python scripts required to simulate 3D magnetic field line parallel connection lengths ($L_c$) for the Helically Symmetric eXperiment (HSX) using the MGRID format. These $L_c$ values are used to set the length of the magnetic field line over which a diffusion equation describing the propagation of a disturbance in the density of the plasma is solved.

Because this project relies on the **FLARE** physics engine and the **MOOSE** framework (which must be compiled from source), you cannot simply `pip install` the requirements. 

Please follow these exact instructions to replicate the Ubuntu/WSL environment.

## 1. System Requirements & System Packages
You must be running Linux (or Ubuntu via WSL) with Miniconda installed. 

First, install the underlying C/Fortran compilers and the MPI/Math libraries. Open your terminal and run:
```bash
sudo apt update
sudo apt install build-essential cmake gfortran pkgconf \
                 libnetcdf-dev libnetcdff-dev \
                 liblapack-dev openmpi-bin libopenmpi-dev
```

## 2. Conda Environment Setup
We use Conda for the base Python environment, but you must install `mpi4py` via `pip` to ensure it compiles against the system's OpenMPI and `gfortran` libraries. **Do not use Conda to install MPI.**

```bash
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

```bash
git clone https://gitlab.com/hfrerichs/moose.git
git clone https://gitlab.com/hfrerichs/flare.git
```

## 4. Compile the MOOSE Framework
MOOSE must be compiled first and "installed" into a local hidden directory so FLARE can find its configuration blueprints.

```bash
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

```bash
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
```bash
mkdir -p ~/DATABASE/flare/HSX/mgrid
```

**Step 2: Add the MGRID File**
Place the massive 3D magnetic grid file (`mgrid_res2p5cm_180pln.nc` or another magnetic grid file if you have one) directly into the `~/DATABASE/flare/HSX/mgrid/` folder. 

**Step 3: Create the `.bfield` Configuration**
Inside the `mgrid` folder, create a hidden file named `.bfield`:
```bash
nano ~/DATABASE/flare/HSX/mgrid/.bfield
```
Paste the following configuration, which specifies the file name and the exact coil current amplitudes (in Amperes) required to generate the QHS state:
```ini
[DEFAULT]

[equi3d_mgrid]
filename: mgrid_res2p5cm_180pln 1.nc
amplitudes: [-1.0722E+04,-1.0722E+04,-1.0722E+04,-1.0722E+04,-1.0722E+04,-1.0722E+04,0.00,0.00,0.00,0.00,0.00,0.0000E+00]
dtype: 'magnetic_field'
```
*NOTE: the files and configuration specifications used throughout this code were provided by Dr. Dieter Boeyaert*
*NOTE: ensure the file name matches the file name you have in the corresponding folder*
**Step 4: Add the 3D Vessel Mesh**
Place the vacuum vessel coordinate file (`vessel_hsx_flare.txt` or another coordinate file if you have one) directly into the `~/DATABASE/flare/HSX/mgrid/` folder. Create a `.boundary` configuration file to define the 3D `torosurf` shape:
```bash
nano ~/DATABASE/flare/HSX/mgrid/.boundary
```
Paste the following text exactly as written:
```
[DEFAULT]

[torosurf]
filename: vessel_hsx_flare.txt
```

## 7. Running the Code
You do not need to install MOOSE or FLARE into your Python environment. Our Python scripts handle this dynamically. As long as the folders are named `moose` and `flare` and sit in the root of this project directory, scripts will automatically append the build paths using `sys.path.insert()`.

**1. Calculate Parallel Connection Lengths:**
Because the FLARE C++ wrapper operates via file I/O, this script uses `moose.grids` to generate a `grid.dat` input file for the Fortran engine, executes the field-line trace against the 3D vessel mesh, and outputs data to `lc.dat`.
```bash
python parallel_connection_length_finder.py
```

## 8. Run the 1D Diffusion Equation Solver
Set the parallel connection length, the desired electron temperature of the plasma in eV, the total amount of time you want the diffusion equation to be solved for, the number of spatial grid points you want the diffusion equation solved for along the magnetic field line with length equal to the parallel connection length, the number of time steps you want the diffusion equation solved for, and the baseline plasma density. Note that if your settings cause the diffusion equation solving method (Forward Time Center Space) to become numerically unstable, an error will be thrown suggesting a new number of timesteps to calculate over to make the solver stable. You can also modify the shape (amplitude, initial position, width) of the density disturbance. The solver will save the calculated density values for each time and spatial position along the field line to a CSV and then animate the propagation of the disturbance. Set the recalculate_data flag to be true if you want to recalculate the density at every point along the field line and at every time and then animate the new data. If you leave the flag set to be false, it will load data from an already saved csv and animate that instead.
