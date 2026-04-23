import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import csv

# ─────────────────────────────────────────
# PARAMETERS — edit these freely
# ─────────────────────────────────────────

L = 1.02          # parallel connection length (chosen values: [1.02, 1.14, 1.33, 1.51, 3.33])
Te = 5              # electron temperature in eV (chosen values: [5,6,7,8,9,10])
me = 9.11 * 10**(-31)    # electron mass in kg
e = 1.602 * 10**(-19)   # elementary charge
epsilon = 8.85 * 10**(-12)   # permittivity of free space (check units)
N = 1000            # number of spatial grid points
total_sim_time = 5 * 10**(-8)   # total time simulator runs (modify)
n_steps = 150000   # total number of time steps to simulate (chosen values: {3.33: [], })
animate_every = 1000  # only render every Nth frame (keeps animation smooth)

n_baseline = 2 * 10**17      # baseline plasma density in m^(-3)
alpha = 1                    # factor multiplied onto parallel diffusion parameter to account for 

# Gaussian initial condition parameters
bump_center = L / 2   # center of the Gaussian bump
bump_width = 0.05     # standard deviation (controls how wide the spike is) (modify)
bump_height = n_baseline / 3     # peak amplitude (modify)

d_min_thres = (1/np.e) * bump_height + n_baseline   # disturbance height at center of disturbance used to measure how fast decay dissipates

# if the boolean is set to True the code will run the simulation, save the data to a csv file, and animate the diffusion. if the boolean is set to False the code will skip the simulation and animate using the data in the csv file 
recalculate_data = True

# this variable should be set to the name of the csv you want to create or save to if the recalculate_data boolean is set to True, or the name of the csv you want to read from if the recalculate_data boolean is set to False
csv_filename = "diffusion_results.csv"

# this variable should be set to the name of the gif you want to create corresponding to the animation
gif_filename = "diffusion_animation.gif"

# only recalculate data if the boolean is set to True, otherwise just read the data from the csv file and skip the simulation

if recalculate_data:
    # ─────────────────────────────────────────
    # GRID SETUP
    # ─────────────────────────────────────────

    s = np.linspace(0, L, N)       # spatial grid along field line
    ds = s[1] - s[0]               # grid spacing
    dt = total_sim_time / n_steps         # time step in seconds (modify)


    # ─────────────────────────────────────────
    # INITIAL CONDITIONS
    # ─────────────────────────────────────────

    n_fluc = bump_height * np.exp(-0.5 * ((s - bump_center) / bump_width)**2)
    n = n_baseline + n_fluc
    n_initial = n.copy()

    L_debye = 7430 * np.sqrt((Te) / n_initial[1:-1])
    Lambda = 12 * np.pi * n_initial[1:-1] * L_debye**3
    vth = np.sqrt((Te * e) / me)
    nu_ei = (n_initial[1:-1] * (e**4) * np.log(Lambda)) / (3 * np.pi**(3/2) * epsilon**2 * me**(1/2) * (2*Te*e)**(3/2))
    D_parallel = alpha * vth**2 / nu_ei

    # Stability check: diffusion number must be <= 0.5 for explicit scheme
    r_init = np.max(D_parallel) * dt / ds**2
    if r_init > 0.5:
        
        new_N_steps = math.ceil(np.max(D_parallel) * (total_sim_time / (0.5 * ds **2)) + 1)

        raise ValueError(
            f"Unstable! Diffusion number r = {r_init:.3f} > 0.5. "
            f"Increase number of time steps to {new_N_steps}"
        )
    else:
        print(f"Stability check passed: r = {r_init:.4f} (must be ≤ 0.5)")


    # ─────────────────────────────────────────
    # TIME STEPPING FUNCTION (explicit FTCS)
    # ─────────────────────────────────────────
    # Solves: dn/dt = D * d²n/ds²
    # Using forward-time, centered-space (FTCS) finite differences with insulated boundary conditions

    def step(n_baseline, n_fluc):
        """Advance n by one time step using FTCS scheme with zero-flux boundaries."""
        
        n = n_baseline + n_fluc

        L_debye = 7430 * np.sqrt((Te) / n[1:-1])    # Debye Length
        Lambda = 12 * np.pi * n[1:-1] * L_debye**3  # Coulomb Factor (for Coulomb logarithm)
        nu_ei = (n_initial[1:-1] * (e**4) * np.log(Lambda)) / (3 * np.pi**(3/2) * epsilon**2 * me**(1/2) * (2*Te*e)**(3/2))     # collision frequency between electrons and ions
        D_parallel = alpha * vth**2 / nu_ei   # parallel diffusion coefficient (m^2/s or arbitrary units)
        r = (D_parallel * dt) / ds**2
        
        n_fluc_new = n_fluc.copy()
        
        # Interior points
        n_fluc_new[1:-1] = n_fluc[1:-1] + r * (n_fluc[2:] - 2*n_fluc[1:-1] + n_fluc[:-2])
        # Boundary conditions (flux equal at both ends)
        n_fluc_new[0] = n_fluc_new[1]
        n_fluc_new[-1] = ((3) * n_fluc_new[0] - 4 * n_fluc_new[1] + n_fluc_new[2] - n_fluc_new[-3] + 4 * n_fluc_new[-2]) / 3
        return n_fluc_new

    # ─────────────────────────────────────────
    # PRE-COMPUTE FRAMES FOR ANIMATION
    # ─────────────────────────────────────────

    frames = []
    times  = []

    n_fluc_current = n_fluc.copy()
    n_current = n_initial.copy()

    decay_time = False

    for i in range(n_steps):
        
        # animation code
        if i % animate_every == 0:
            frames.append(n_current.copy())
            times.append(i * dt)
        
        # calculate the next set of density values
        n_fluc_current = step(n_baseline, n_fluc_current)
        n_current = n_baseline + n_fluc_current

        # if the current disturbance height is below the threshold, record the time at which it crossed the threshold
        if np.max(n_current) <= d_min_thres and decay_time == False:
            
            decay_time = i * dt


    print(f"Decay time was: {decay_time} seconds")

    # ─────────────────────────────────────────
    # EXPORT DATA TO CSV
    # ─────────────────────────────────────────

    print(f"Writing data to {csv_filename}...")

    # Open the file in write mode
    with open(csv_filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        
        # 1. Create and write the header row
        # Column 1 is "Position (m)", followed by columns for each recorded time step
        header = ["Position_m"] + [f"Density_t={t:.8f}s" for t in times]
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

    print(f"Data successfully saved to {csv_filename}")

    # ─────────────────────────────────────────
    # ANIMATION
    # ─────────────────────────────────────────

    fig, ax = plt.subplots(figsize=(10, 5))
    fig.patch.set_facecolor('#0f0f1a')
    ax.set_facecolor('#0f0f1a')

    # Plot initial condition as faint reference
    ax.plot(s, n_initial, color='white', alpha=0.15, linewidth=1.2, linestyle='--', label='Initial condition')

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
    ax.set_ylim(0.95*n_baseline, 1.5*n_baseline)

    # time_text = ax.text(0.02, 0.93, '', transform=ax.transAxes, color='#ffdd88', fontsize=11, fontfamily='monospace')

    legend = ax.legend(loc='upper right', framealpha=0.2, labelcolor='white')
    ax.set_title('1D Parallel Diffusion Along a Field Line', color='white', fontsize=13, pad=12)

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
    
    # Save the animation as a GIF file using Pillow writer
    print(f"Saving animation to {gif_filename}... (this may take a minute or two)")
    ani.save(gif_filename, writer='pillow', fps=30)
    print(f"GIF successfully saved to {gif_filename}!")
    plt.tight_layout()
    plt.show()

else:
    print(f"Reading data from {csv_filename}...")

    s_vals = []
    frames = []
    times = []

    try:
        with open(csv_filename, 'r') as file:
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
        print(f"Error: Could not find '{csv_filename}'. Make sure it's in the same folder!")
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
    n_initial = frames[0]          # The very first time step is our initial condition
    n_baseline = np.min(n_initial) # The edges of your initial Gaussian are the baseline

    # ─────────────────────────────────────────
    # 3. RENDER THE ANIMATION (Original Style)
    # ─────────────────────────────────────────
    print("Rendering animation...")
    fig, ax = plt.subplots(figsize=(10, 5))
    fig.patch.set_facecolor('#0f0f1a')
    ax.set_facecolor('#0f0f1a')

    # Plot initial condition as faint reference
    ax.plot(s, n_initial, color='white', alpha=0.15, linewidth=1.2, linestyle='--', label='Initial condition')

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

    # Save the animation as a GIF file using Pillow writer
    print(f"Saving animation to {gif_filename}... (this may take a minute or two)")
    ani.save(gif_filename, writer='pillow', fps=30)
    print(f"GIF successfully saved to {gif_filename}!")

    plt.tight_layout()
    plt.show()