import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Constants
G = 6.67430e-6
dt = 0.01
num_steps = 2000

# Masses: central mass + 3 satellites
masses = np.array([10000, 1, 1, 1])

# Initial positions (x, y) in 2D
positions = np.array([
    [0.0, 0.0],    # central mass
    [1.5, 0.0],    # satellite 1
    [0.0, 2.0],    # satellite 2
    [-2.0, 0.0]    # satellite 3
], dtype=float)

# Initialize velocities array
velocities = np.zeros_like(positions)

# Calculate initial velocities for circular orbits
for i in range(1, len(masses)):
    r_vec = positions[i] - positions[0]
    r = np.linalg.norm(r_vec)
    speed = np.sqrt(G * masses[0] / r)
    # Velocity perpendicular to radius vector (in 2D)
    vel_dir = np.array([-r_vec[1], r_vec[0]])
    vel_dir /= np.linalg.norm(vel_dir)
    velocities[i] = speed * vel_dir

num_bodies = len(masses)
pos_history = np.zeros((num_bodies, num_steps, 2))

def compute_accelerations(pos, mass):
    acc = np.zeros_like(pos)
    softening = 1e-2
    for i in range(num_bodies):
        for j in range(num_bodies):
            if i != j:
                diff = pos[j] - pos[i]
                dist_sq = np.sum(diff**2) + softening**2
                inv_dist3 = 1.0 / (dist_sq * np.sqrt(dist_sq))
                acc[i] += G * mass[j] * diff * inv_dist3
    return acc

# Verlet integration setup
pos_old = positions - velocities * dt

for step in range(num_steps):
    acc = compute_accelerations(positions, masses)
    pos_new = 2 * positions - pos_old + acc * dt**2
    velocities = (pos_new - pos_old) / (2 * dt)

    pos_history[:, step, :] = pos_new

    pos_old = positions
    positions = pos_new

# Plot & animate in 2D
fig, ax = plt.subplots()
ax.set_xlim(-3, 3)
ax.set_ylim(-3, 3)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_title('2D N-Body Orbital Simulation')
ax.grid(True)
ax.set_aspect('equal')

colors = ['yellow', 'blue', 'green', 'red']
sizes = [150, 50, 50, 50]

scatters = [ax.scatter([], [], s=size, color=color) for size, color in zip(sizes, colors)]

def update(frame):
    for i, scatter in enumerate(scatters):
        scatter.set_offsets(pos_history[i, frame])
    return scatters

ani = FuncAnimation(fig, update, frames=num_steps, interval=20, blit=True)
plt.show()
