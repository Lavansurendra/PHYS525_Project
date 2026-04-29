
# importing neccesary libraries for running Fortran/C++ files set up by make command in moose and flare directories
import sys                  # Lets us modify how Python searches for files
import os                   # Lets us interact with the computer's file system (folders/paths)

# importing numpy and matplotlib for calculations and visualization of magnetic field strength
import numpy as np          # The standard library for heavy array math
import matplotlib.pyplot as plt  # The standard library for plotting graphs
import csv                  # Lets us read and write spreadsheet files

def Field_Strength_Finder(num_minor_radius_points, phi_angle):
    # specifying the location of the .so files for flare and moose code 
        # NOTE: the insert command allows us to specify what order the directories should be searched in when we run the flare import
        #       since we specified an index of 0 this command ensures that the specified directories are searched first
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), 'moose', 'src', 'python')))
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), 'flare', 'src', 'python')))


    # creating a temporary system variable called 'DATABASE' pointing to your ~/DATABASE folder.
    os.environ["DATABASE"] = os.path.expanduser("~/DATABASE")

    # importing functions from the flare and moose directories that we will need to use to run the C++ code of moose and flare
    from flare import model
    from flare.analysis import bfield

    # using the HSX mgrid model given by Dr. Dieter Boeyaert and generating a continuous 3d mathematical model of the 
    # magnetic field
    # the flare code run here can calculate the magnetic vector B to a fraction of a millimeter at any point
    # within the device using 3D interpolation

    try:
        print("Loading HSX mgrid model...")
        model.load("HSX/mgrid")
        print("Successfully loaded MGRID data.")
    except Exception as e:
        
        # NOTE: accurately loading the model requires that the HSX magnetic field file is correctly placed in the DATABASE folder and that the flare code is correctly set up to read it
        print(f"Error loading MGRID data: {e}")
        sys.exit(1)


    # defining the grid of points where we want to calculate the magnetic field strength |B|
    # we will be calculating |B| across a 2D slice of the outboard midplane of the plasma (z=0) for a range of minor radii and toroidal angles 
    # this will allow us to visualize how the magnetic field strength changes as we move outward from the core of the plasma to the edge and as we move around the stellerator

    # defining the major radius of the magnetic axis in meters (from HSX documentation)
    R_axis = 1.20   

    # defining the minor radius (user specified number of points)
        # NOTE: from the hsx website we see that the average plasma minor radius is 0.12 meters. However since HSX is very twisted the plasma is not a perfect cylinder
        # correspondingly at at some toroidal angles the cross section of the plasma can be wider than 0.12 meters
        # to ensure we capture the full plasma cross section we will calculate |B| up to a minor radius of 0.15 meters
    r_grid = np.linspace(0.0, 0.15, num_minor_radius_points)

    # creating an empty 1D array to store our calculated field strengths
    B_strength = np.zeros(len(r_grid))


    print(f"Calculating |B| at phi = {phi_angle} rad...")

    for i, r in enumerate(r_grid):
        
        # translating 1D minor radius into a 3D coordinate on the midplane
        R = R_axis + r
        Z = 0.0  
        
        # using flare code to get the 3D vector [Br, Bz, Bphi] of the magnetic field at the specified point
        B_vector = bfield.eval(R, Z, phi_angle)
        
        # calculating the scalar magnitude |B| and saving it
        B_strength[i] = np.linalg.norm(B_vector)

    print("Data extraction complete!")


    #  creating a data folder if it doesn't already exist
    output_folder = os.path.join("..", "data", "perp_diffusion_solver")
    os.makedirs(output_folder, exist_ok=True)

    # naming the file based on the user specified angle
    csv_filename = os.path.join(output_folder, f"hsx_bfield_1d_phi_{phi_angle}.csv")

    print(f"Saving 1D data to {csv_filename}...")

    with open(csv_filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        
        # writing a simple two-column header
        writer.writerow(["r_m", f"B_tesla_phi_{phi_angle}"])
        
        # writing the radius and corresponding field strength row by row
        for i in range(len(r_grid)):
            writer.writerow([r_grid[i], B_strength[i]])


    # initializing a canvas for plotting the 1D profile of |B| vs minor radius r at the user specified toroidal angle
    fig, ax = plt.subplots(figsize=(8, 5))
    fig.patch.set_facecolor('#0f0f1a')
    ax.set_facecolor('#0f0f1a')

    # plotting the 1D profile
    ax.plot(r_grid, B_strength, color='#00cfff', linewidth=2.0, label=f'|B| at $\phi$ = {phi_angle}')

    # labelling and formatting
    ax.set_xlabel("Minor Radius  r  (m)", color='white', fontsize=12)
    ax.set_ylabel("Magnetic Field Strength  |B|  (Tesla)", color='white', fontsize=12)
    ax.set_title(f"HSX Radial Magnetic Field Profile ($\phi$ = {phi_angle})", color='white', fontsize=13)
    ax.tick_params(colors='white')

    # customizing the spines and legend
    for spine in ax.spines.values():
        spine.set_edgecolor('#444466')

    # placing the legend in the upper right corner with a transparent background and white text
    ax.legend(loc='upper right', framealpha=0.2, labelcolor='white')

    # naming the file based on the user specified angle
    plot_filename = os.path.join(output_folder, f"hsx_bfield_1d_phi_{phi_angle}.png")
    print(f"Saving plot image to {plot_filename}...")
    fig.savefig(plot_filename, dpi=300, facecolor=fig.get_facecolor(), edgecolor='none')

    plt.tight_layout()
    plt.show()

    return B_strength

if __name__ == "__main__":

    # Example usage: calculate and plot |B| vs minor radius at phi = 0.0 radians with 100 points
    # Field_Strength_Finder(num_minor_radius_points=1001, phi_angle=0)
    # Field_Strength_Finder(num_minor_radius_points=1001, phi_angle=0.7854)
    # Field_Strength_Finder(num_minor_radius_points=1001, phi_angle=1.5708)
    # Field_Strength_Finder(num_minor_radius_points=1001, phi_angle=2.3562)
    # Field_Strength_Finder(num_minor_radius_points=1001, phi_angle=3.1416)
    # Field_Strength_Finder(num_minor_radius_points=1001, phi_angle=3.9270)
    # Field_Strength_Finder(num_minor_radius_points=1001, phi_angle=4.7124)
    # Field_Strength_Finder(num_minor_radius_points=1001, phi_angle=5.4978)
    Field_Strength_Finder(num_minor_radius_points=1001, phi_angle=6.2832)