# PHYS525 Project: Stellarator 1D Diffusion Analysis & Parallel Connection Length

This repository contains the Python scripts and data required to extract 1D magnetic geometry from VMEC equilibria for a Crank-Nicolson diffusion solver, as well as simulating 3D magnetic field line parallel connection lengths ($L_c$) using the MGRID format.

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
We use Conda for the base Python environment, but we must install `mpi4py` via `pip` to ensure it compiles against the system's OpenMPI and `gfortran` libraries. **Do not use Conda to install MPI.**

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

*(Note: These folders should be in the `.gitignore` so we don't accidentally push compiled binaries to the repo).*

```bash
git clone [https://gitlab.com/hfrerichs/moose.git](https://gitlab.com/hfrerichs/moose.git)
git clone [https://gitlab.com/hfrerichs/flare.git](https://gitlab.com/hfrerichs/flare.git)
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
Now we build FLARE, explicitly pointing it to the local MOOSE installation we just created.

```bash
cd flare
mkdir build
cd build

# Configure CMake (Point to MOOSE targets in .local)
cmake -DCMAKE_PREFIX_PATH=$HOME/.local -DPython3_EXECUTABLE=$(which python) ..

# Compile the physics engine
make

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
Place the massive 3D magnetic grid file (`mgrid_res2p5cm_180pln.nc`) directly into the `~/DATABASE/flare/HSX/mgrid/` folder. *(Note: Do not use the VMEC `wout_*.nc` files for the connection length trace; they lack vacuum field data).*

**Step 3: Create the `.bfield` Configuration**
Inside the `mgrid` folder, create a hidden file named `.bfield`:
```bash
nano ~/DATABASE/flare/HSX/mgrid/.bfield
```
Paste the following configuration, which specifies the file name and the exact coil current amplitudes (in Amperes) required to generate the QHS state:
```ini
[DEFAULT]

[equi3d_mgrid]
filename: mgrid_res2p5cm_180pln.nc
amplitudes: [-1.0722E+04,-1.0722E+04,-1.0722E+04,-1.0722E+04,-1.0722E+04,-1.0722E+04,0.00,0.00,0.00,0.00,0.00,0.0000E+00]
dtype: 'magnetic_field'
```

**Step 4: Create the `.boundary` Bypass**
FLARE's internal parser will refuse to load the model unless a boundary configuration file exists in the same directory. Until the physical CAD mesh of the HSX vessel is acquired, create an **empty** boundary file to bypass the security check without triggering a build error:
```bash
touch ~/DATABASE/flare/HSX/mgrid/.boundary
```

## 7. Running the Code
You do not need to install MOOSE or FLARE into your Python environment. Our Python scripts handle this dynamically. As long as the folders are named `moose` and `flare` and sit in the root of this project directory, scripts will automatically append the build paths using `sys.path.insert()`.

**To test the 1D VMEC extraction:**
```bash
python extract_hsx.py
```

**To test the 3D MGRID Parallel Connection Length trace:**
Because the FLARE C++ wrapper operates via file I/O, this script uses `moose.grids` to generate a `grid.dat` input file for the Fortran engine, executes the trace, and reads the output from `lc.dat`.
```bash
python parallel_connection_length_finder.py
```

## 8. Known Limitations & Next Steps
* **Missing Vessel Mesh:** The current parallel connection length simulation is operating "wall-free." Because there is no physical boundary to stop the trace, particles eventually crash into the physical electromagnets, causing the toroidal field to drop to zero and throwing a Runge-Kutta integrator failure (`Error 7`). 
* **Next Step:** Acquire the 3D first-wall mesh of the HSX vacuum vessel, update the `.boundary` config file to include it, and re-run to get physically accurate $L_c$ measurements.