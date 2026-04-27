import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import csv
import os

from field_strength_finder import Field_Strength_Finder

# ─────────────────────────────────────────
# PARAMETERS — edit these freely
# ─────────────────────────────────────────

# simulation parameters
N = 1001            # number of spatial grid points (NOTE: this must be an odd number for the finite difference method in cylindrical coordinates to work properly)
total_sim_time = 5 * 10**(-8)   # total time simulator runs (modify)
n_steps = 10000   # total number of time steps to simulate (chosen values: {3.33: [12000, 19000, 27000, 37000, 49000, 62000], })
animate_every = 1000  # only render every Nth frame (keeps animation smooth)
obs_phi = np.pi         # the toroidal angle at which you want to observe the diffusion of temperature in HSX in radians

# constants
me = 9.11 * 10**(-31)    # electron mass in kg
e = 1.602 * 10**(-19)   # elementary charge
epsilon = 8.85 * 10**(-12)   # permittivity of free space (check units)

# plasma parameters
minor_r = 0.12       # minor radius of HSX (m)
Te = 2500              # electron temperature in eV (chosen values: [5,6,7,8,9,10])
n_baseline = 2 * 10**17      # baseline plasma density in m^(-3)
B = Field_Strength_Finder(N, obs_phi)

# Gaussian initial condition parameters
bump_center = 0   # center of the Gaussian bump
bump_width = 0.05     # standard deviation (controls how wide the spike is) (modify)
bump_height = Te / 10     # peak amplitude (modify)

d_min_thres = (1/np.e) * bump_height + Te   # disturbance height at center of disturbance used to measure how fast decay dissipates

# if the boolean is set to True the code will run the simulation, save the data to a csv file, and animate the diffusion. if the boolean is set to False the code will skip the simulation and animate using the data in the csv file 
recalculate_data = True

# this variable should be set to the name of the csv you want to create or save to if the recalculate_data boolean is set to True, or the name of the csv you want to read from if the recalculate_data boolean is set to False
csv_filename = "dif_151_5.csv"

# this variable should be set to the name of the gif you want to create corresponding to the animation
gif_filename = "dif_151_5_animation.gif"

# only recalculate data if the boolean is set to True, otherwise just read the data from the csv file and skip the simulation

if recalculate_data:
    # ─────────────────────────────────────────
    # GRID SETUP
    # ─────────────────────────────────────────

    s = np.linspace(-minor_r, minor_r, N)       # spatial grid along field line
    center = len(s) // 2
    ds = s[1] - s[0]               # grid spacing
    dt = total_sim_time / n_steps         # time step in seconds (modify)


    # ─────────────────────────────────────────
    # INITIAL CONDITIONS
    # ─────────────────────────────────────────

    T_fluc = bump_height * np.exp(-0.5 * ((s - bump_center) / bump_width)**2)
    T = Te + T_fluc
    T_initial = T.copy()

    L_debye = np.sqrt((epsilon * T_initial[1:-1]) / (n_baseline * e**2))
    Lambda = 12 * np.pi * n_baseline * L_debye**3
    nu_ei = (n_baseline * (e**4) * np.log(Lambda)) / (3 * np.pi**(3/2) * epsilon**2 * me**(1/2) * (2*T_initial[1:-1]*e)**(3/2))
    vth = np.sqrt((T_initial[1:-1] * e) / me)
    r_larmor = (me**2 * vth**2) / (e**2 * B[1:-1])
    chi_perp = r_larmor**2 * nu_ei

    # Stability check: diffusion number must be <= 0.5 for explicit scheme
    r_init = np.max(chi_perp) * dt / ds**2
    if r_init > 0.5:
        
        new_N_steps = math.ceil(np.max(chi_perp) * (total_sim_time / (0.5 * ds **2)) + 1)

        raise ValueError(
            f"Unstable! Diffusion number r = {r_init:.3f} > 0.5. "
            f"Increase number of time steps to {new_N_steps}"
        )
    else:
        print(f"Stability check passed: r = {r_init:.4f} (must be ≤ 0.5)")


    # ─────────────────────────────────────────
    # TIME STEPPING FUNCTION (explicit FTCS)
    # ─────────────────────────────────────────
    # Solves: dn/dt = 1/r(d/dr(r*chi(dT/dr)))
    # Using forward-time, centered-space (FTCS) finite differences with insulated boundary conditions

    def step(T_fluc):
        
        """Advance n by one time step using cylindrical coordinate FTCS scheme with zero-flux boundaries."""
        
        T = Te + T_fluc

        L_debye = np.sqrt((epsilon * T[1:-1]) / (n_baseline * e**2))
        Lambda = 12 * np.pi * n_baseline * L_debye**3
        nu_ei = (n_baseline * (e**4) * np.log(Lambda)) / (3 * np.pi**(3/2) * epsilon**2 * me**(1/2) * (2*T[1:-1]*e)**(3/2))
        vth = np.sqrt((T[1:-1] * e) / me)
        r_larmor = (me**2 * vth**2) / (e**2 * B[1:-1])
        chi_perp = r_larmor**2 * nu_ei
        r = chi_perp * dt
        
        T_fluc_new = T_fluc.copy()
        
        # Interior points
        T_fluc_new[1:-1] = T_fluc[1:-1] + r * ((T_fluc[2:] - 2*T_fluc[1:-1] + T_fluc[:-2]) / ds **2 + (1/s[1:-1]) * (T_fluc[2:] - T_fluc[:-2]) / (2*ds))
        
        # Insulated Boundary conditions (we assume no heat transfer out of the plasma to the walls)
        T_fluc_new[0] = T_fluc_new[1]
        T_fluc_new[-1] = T_fluc_new[-2]
        T_fluc_new[center] = T_fluc[center] + 2 * r * (T_fluc[center+1] - 2*T_fluc[center] + T_fluc[center-1]) / ds **2

    # ─────────────────────────────────────────
    # PRE-COMPUTE FRAMES FOR ANIMATION
    # ─────────────────────────────────────────

    frames = []
    times  = []

    T_fluc_current = T_fluc.copy()
    T_current = T_initial.copy()

    decay_time = False

    for i in range(n_steps):
        
        # animation code
        if i % animate_every == 0:
            frames.append(T_current.copy())
            times.append(i * dt)
        
        # calculate the next set of density values
        T_fluc_current = step(T_fluc_current)
        T_current = Te + T_fluc_current

        # if the current disturbance height is below the threshold, record the time at which it crossed the threshold
        if np.max(T_current) <= d_min_thres and decay_time == False:
            
            decay_time = i * dt


    print(f"Decay time was: {decay_time} seconds")

    # ─────────────────────────────────────────
    # EXPORT DATA TO CSV
    # ─────────────────────────────────────────
    
    # defining the output folder we will save our data to, and creating the folder if it doesn't already exist
    output_folder = os.path.join("..", "data", "perp_diffusion_solver")
    os.makedirs(output_folder, exist_ok=True)

    # defining the full path to the csv file we will save our data to
    path_name = os.path.join(output_folder, csv_filename)

    print(f"Writing data to {path_name}...")

    # Open the file in write mode
    with open(path_name, mode='w', newline='') as file:
        writer = csv.writer(file)
        
        # 1. Create and write the header row
        # Column 1 is "Position (m)", followed by columns for each recorded time step
        header = ["Position_m"] + [f"Temperature_t={t:.8f}s" for t in times]
        writer.writerow(header)
        
        # 2. Write the data rows
        # We loop through every spatial index 'i' (from 0 to N)
        for i in range(N):
            # Start the row with the spatial coordinate s[i]
            row = [s[i]]
            
            # Add the density value at this exact position for every saved time frame
            for frame in frames:
                row.append(frame[i])
                
            # Write the completed row to the CSV
            writer.writerow(row)

    print(f"Data successfully saved to {path_name}")

    # ─────────────────────────────────────────
    # ANIMATION
    # ─────────────────────────────────────────

    fig, ax = plt.subplots(figsize=(10, 5))
    fig.patch.set_facecolor('#0f0f1a')
    ax.set_facecolor('#0f0f1a')

    # Plot initial condition as faint reference
    ax.plot(s, T_initial, color='white', alpha=0.15, linewidth=1.2, linestyle='--', label='Initial condition')

    # Main evolving line
    line, = ax.plot(s, frames[0], color='#00cfff', linewidth=2.0, label='n(s, t)')

    # Filled area under curve
    fill = ax.fill_between(s, frames[0], alpha=0.15, color='#00cfff')

    # Labels and formatting
    ax.set_xlabel('Distance along field line  s  (m)', color='white', fontsize=12)
    ax.set_ylabel('Perturbation Temperature  T(r, phi, t)', color='white', fontsize=12)
    ax.tick_params(colors='white')
    for spine in ax.spines.values():
        spine.set_edgecolor('#444466')
    ax.set_xlim(0, 2*minor_r)
    ax.set_ylim(0.95*Te, 1.5*Te)

    legend = ax.legend(loc='upper right', framealpha=0.2, labelcolor='white')
    ax.set_title('1D Perpendicular Diffusion Across Field Lines', color='white', fontsize=13, pad=12)

    def update(frame_idx):
        global fill
        y = frames[frame_idx]

        line.set_ydata(y)

        # Redraw fill
        for coll in ax.collections:
            coll.remove()
        ax.fill_between(s, y, alpha=0.15, color='#00cfff')

        return line

    ani = animation.FuncAnimation(
        fig, update,
        frames=len(frames),
        interval=1,       # milliseconds between frames
        blit=False
    )
    
   # defining the full path to the animation file we will save our data to
    path_name = os.path.join(output_folder, gif_filename)

    # Save the animation as a GIF file using Pillow writer
    print(f"Saving animation to {path_name}... (this may take a minute or two)")
    ani.save(path_name, writer='pillow', fps=30)
    print(f"GIF successfully saved to {gif_filename}!")
    plt.tight_layout()
    plt.show()

else:
    print(f"Reading data from {csv_filename}...")

    s_vals = []
    frames = []
    times = []

    # defining the output folder we will save our data to, and creating the folder if it doesn't already exist
    output_folder = os.path.join("..", "data", "perp_diffusion_solver")
    os.makedirs(output_folder, exist_ok=True)

    # defining the full path to the csv file we will save our data to
    path_name = os.path.join(output_folder, csv_filename)

    try:
        with open(path_name, 'r') as file:
            reader = csv.reader(file)
            header = next(reader)
            
            # Parse the time step from the column headers (e.g., "Density_t=0.00000050s")
            for col in header[1:]:
                time_str = col.split('=')[1].replace('s', '')
                times.append(float(time_str))
            
            # Pre-allocate an empty list for every time frame
            for _ in times:
                frames.append([])
                
            # Read the rows: column 0 is position, the rest are densities at different times
            for row in reader:
                s_vals.append(float(row[0]))
                for i, val in enumerate(row[1:]):
                    frames[i].append(float(val))
                    
    except FileNotFoundError:
        print(f"Error: Could not find '{path_name}'. Make sure it's in the correct folder!")
        exit()

    # Convert standard Python lists to Numpy arrays for Matplotlib
    s = np.array(s_vals)
    frames = [np.array(f) for f in frames]

    print(f"Successfully loaded {len(frames)} frames across {len(s)} spatial points.")

    # ─────────────────────────────────────────
    # 2. INFER PLOTTING PARAMETERS FROM DATA
    # ─────────────────────────────────────────
    # Rather than hardcoding your parameters, we can calculate them from the data:
    L = s[-1]                      # The maximum distance along the field line
    T_initial = frames[0]          # The very first time step is our initial condition
    n_baseline = np.min(T_initial) # The edges of your initial Gaussian are the baseline

    # ─────────────────────────────────────────
    # 3. RENDER THE ANIMATION (Original Style)
    # ─────────────────────────────────────────
    print("Rendering animation...")
    fig, ax = plt.subplots(figsize=(10, 5))
    fig.patch.set_facecolor('#0f0f1a')
    ax.set_facecolor('#0f0f1a')

    # Plot initial condition as faint reference
    ax.plot(s, T_initial, color='white', alpha=0.15, linewidth=1.2, linestyle='--', label='Initial condition')

    # Main evolving line
    line, = ax.plot(s, frames[0], color='#00cfff', linewidth=2.0, label='n(s, t)')

    # Filled area under curve
    fill = ax.fill_between(s, frames[0], alpha=0.15, color='#00cfff')

    # Labels and formatting
    ax.set_xlabel('Distance along field line  s  (m)', color='white', fontsize=12)
    ax.set_ylabel('Perturbation density  n(s, t)', color='white', fontsize=12)
    ax.tick_params(colors='white')
    for spine in ax.spines.values():
        spine.set_edgecolor('#444466')

    ax.set_xlim(0, L)
    ax.set_ylim(0.95 * n_baseline, 1.5 * n_baseline)

    # Using scientific notation (.2e) for time since the total sim time is microscopic
    time_text = ax.text(0.02, 0.93, '', transform=ax.transAxes, color='#ffdd88', fontsize=11, fontfamily='monospace')

    legend = ax.legend(loc='upper right', framealpha=0.2, labelcolor='white')
    ax.set_title('1D Parallel Diffusion Along a Field Line (CSV Data)', color='white', fontsize=13, pad=12)

    def update(frame_idx):
        global fill
        y = frames[frame_idx]

        line.set_ydata(y)

        # Redraw fill
        for coll in ax.collections:
            coll.remove()
        fill = ax.fill_between(s, y, alpha=0.15, color='#00cfff')

        time_text.set_text(f't = {times[frame_idx]:.2e} s')
        return line, time_text

    ani = animation.FuncAnimation(
        fig, update,
        frames=len(frames),
        interval=1,       # milliseconds between frames
        blit=False
    )
    
    
    # defining the full path to the animation file we will save our data to
    path_name = os.path.join(output_folder, gif_filename)

    # Save the animation as a GIF file using Pillow writer
    print(f"Saving animation to {path_name}... (this may take a minute or two)")
    ani.save(path_name, writer='pillow', fps=30)
    print(f"GIF successfully saved to {path_name}!")

    plt.tight_layout()
    plt.show()