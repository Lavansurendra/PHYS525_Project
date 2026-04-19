import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# ─────────────────────────────────────────
# PARAMETERS — edit these freely
# ─────────────────────────────────────────

L = 100.0          # parallel connection length
N = 10000            # number of spatial grid points
kTe = 5              # electron temperature in eV
me = 9.11 * 10**(-31)    # electron mass in kg
epsilon = 8.85 * 10**(-12)   # permittivity of free space (check units)
vth = np.sqrt((2 * kTe) / me)
max_length = 10             # the length of the field line over which density values will be calculated
dt = 0.001          # time step in seconds
n_steps = 50000     # total number of time steps to simulate
animate_every = 5  # only render every Nth frame (keeps animation smooth)

n_baseline = 2 * 10**17      # baseline plasma density in m^(-3)

# Gaussian initial condition parameters
bump_center = max_length / 2   # center of the Gaussian bump
bump_width  = 0.05     # standard deviation (controls how wide the spike is)
bump_height = n_baseline / 3     # peak amplitude

# ─────────────────────────────────────────
# GRID SETUP
# ─────────────────────────────────────────

s = np.linspace(0, max_length, N)       # spatial grid along field line
ds = s[1] - s[0]               # grid spacing


# ─────────────────────────────────────────
# INITIAL CONDITION
# ─────────────────────────────────────────

n_fluc = bump_height * np.exp(-0.5 * ((s - bump_center) / bump_width)**2)
n = n_baseline + n_fluc
n_initial = n.copy()

L_debye = 7430 * np.sqrt((kTe) / n_initial[1:-1])
Lambda = 12 * np.pi * n_initial[1:-1] * L_debye**3
nu_ei = (n_initial[1:-1] * (np.e**4) * np.log(Lambda)) / (3 * epsilon**2 * me**(1/2) * (2*kTe)**(3/2))
D_parallel = vth**2 / nu_ei   # parallel diffusion coefficient (m^2/s or arbitrary units)

# Stability check: diffusion number must be <= 0.5 for explicit scheme
r_init = np.max(D_parallel) * dt / ds**2
if r_init > 0.5:
    raise ValueError(
        f"Unstable! Diffusion number r = {r_init:.3f} > 0.5. "
        f"Reduce dt or increase N."
    )
else:
    print(f"Stability check passed: r = {r_init:.4f} (must be ≤ 0.5)")

# ─────────────────────────────────────────
# TIME STEPPING FUNCTION (explicit FTCS)
# ─────────────────────────────────────────
# Solves: dn/dt = D * d²n/ds²
# Using forward-time, centered-space (FTCS) finite differences

def step(n_baseline, n_fluc):
    """Advance n by one time step using FTCS scheme with zero-flux boundaries."""
    
    n = n_baseline + n_fluc

    L_debye = 7430 * np.sqrt((kTe) / n[1:-1])
    Lambda = 12 * np.pi * n[1:-1] * L_debye**3
    nu_ei = (n[1:-1] * (np.e**4) * np.log(Lambda)) / (3 * epsilon**2 * me**(1/2) * (2*kTe)**(3/2))
    D_parallel = vth**2 / nu_ei   # parallel diffusion coefficient (m^2/s or arbitrary units)
    tau = L**2 / D_parallel

    r = D_parallel * dt / ds**2
    
    n_fluc_new = n_fluc.copy()
    # Interior points
    n_fluc_new[1:-1] = n_fluc[1:-1] + r * (n_fluc[2:] - 2*n_fluc[1:-1] + n_fluc[:-2]) #- n_fluc[1:-1] / tau
    # Neumann boundary conditions (zero flux at both ends)
    n_fluc_new[0] = n_fluc_new[1]
    n_fluc_new[-1] = n_fluc_new[-2]
    return n_fluc_new

# ─────────────────────────────────────────
# PRE-COMPUTE FRAMES FOR ANIMATION
# ─────────────────────────────────────────

frames = []
times  = []

n_fluc_current = n_fluc.copy()
n_current = n_initial.copy()
for i in range(n_steps):
    if i % animate_every == 0:
        frames.append(n_current.copy())
        times.append(i * dt)
    n_fluc_current = step(n_baseline, n_fluc_current)
    n_current = n_baseline + n_fluc_current

print(f"Computed {len(frames)} frames over t = 0 to {n_steps * dt:.4f} s")

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
ax.set_xlim(0, max_length)
ax.set_ylim(0.95*n_baseline, 1.5*n_baseline)

time_text = ax.text(0.02, 0.93, '', transform=ax.transAxes, color='#ffdd88', fontsize=11, fontfamily='monospace')

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

    time_text.set_text(f't = {times[frame_idx]:.2f} s    (D∥ = {D_parallel})')
    return line, time_text

ani = animation.FuncAnimation(
    fig, update,
    frames=len(frames),
    interval=1,       # milliseconds between frames
    blit=False
)

plt.tight_layout()
plt.show()