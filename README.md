# PHYS525 Project: Stellarator 1D Diffusion Analysis

This repository contains the Python scripts and data required to extract 1D magnetic geometry from VMEC equilibria and run a Crank-Nicolson diffusion solver. 

Because this project relies on the **FLARE** physics engine and the **MOOSE** framework (which must be compiled from source), you cannot simply `pip install` the requirements. 

Please follow these exact instructions to replicate the Ubuntu/WSL environment.

## 1. System Requirements & System Packages
You must be running Linux (or Ubuntu via WSL) with Miniconda installed. 

First, install the underlying C/Fortran compilers and the MPI/Math libraries. Open your terminal and run:
`bash
sudo apt update
sudo apt install build-essential cmake gfortran pkgconf \
                 libnetcdf-dev libnetcdff-dev \
                 liblapack-dev openmpi-bin libopenmpi-dev
`

## 2. Conda Environment Setup
We use Conda for the base Python environment, but we must install `mpi4py` via `pip` to ensure it compiles against the system's OpenMPI and `gfortran` libraries. **Do not use Conda to install MPI.**

`bash
# Create and activate the environment
conda create --name PHYS525 python=3.10
conda activate PHYS525

# Install the pre-compiled math libraries via Conda
conda install numpy scipy sympy matplotlib netcdf4

# Install MPI via pip (CRITICAL to avoid Fortran version clashes)
pip install mpi4py
`

## 3. Clone the Sub-Repositories
Make sure you are in the main project folder. You need to pull down the source code for Dr. Frerichs' MOOSE and FLARE libraries. 

*(Note: These folders should be in the `.gitignore` so we don't accidentally push compiled binaries to the repo).*

`bash
git clone https://gitlab.com/hfrerichs/moose.git
git clone https://gitlab.com/hfrerichs/flare.git
`

## 4. Compile the MOOSE Framework
MOOSE must be compiled first and "installed" into a local hidden directory so FLARE can find its configuration blueprints.

`bash
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
`

## 5. Compile the FLARE Physics Engine
Now we build FLARE, explicitly pointing it to the local MOOSE installation we just created.

`bash
cd flare
mkdir build
cd build

# Configure CMake (Point to MOOSE targets in .local)
cmake -DCMAKE_PREFIX_PATH=$HOME/.local -DPython3_EXECUTABLE=$(which python) ..

# Compile the physics engine
make

# Return to the main project directory
cd ../..
`

## 6. Running the Code
You do not need to install MOOSE or FLARE into your Python environment. Our Python scripts handle this dynamically. 

As long as the folders are named `moose` and `flare` and sit in the root of this project directory, scripts like `extract_hsx.py` will automatically append the build paths using `sys.path.append()`.

To verify your environment is working, run:
`bash
python extract_hsx.py
`
If you see the success message, your machine is ready to perform field line traces!