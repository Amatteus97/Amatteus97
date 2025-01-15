import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

def runge_kutta_system(f, yv0, t0, tf, h):
    t_values = [t0]
    yv_values = [yv0]
    t = t0
    yv = yv0

    steps = int((tf - t0) / h)
    for _ in tqdm(range(steps), desc="Processing"):
        k1 = h * np.array(f(t, yv))
        k2 = h * np.array(f(t + h/2, yv + k1/2))
        k3 = h * np.array(f(t + h/2, yv + k2/2))
        k4 = h * np.array(f(t + h, yv + k3))
        
        yv = yv + (k1 + 2*k2 + 2*k3 + k4) / 6
        t = t + h
        
        t_values.append(t)
        yv_values.append(yv)

    return t_values, yv_values

# Define the ODEs for the three-body problem
def three_body_problem(t, yv):
    G = 6.67430e-11  # gravitational constant

    # Unpack positions and velocities
    r1, r2, r3 = yv[:2], yv[2:4], yv[4:6]
    v1, v2, v3 = yv[6:8], yv[8:10], yv[10:12]

    # Masses of the bodies (in kg)
    m1, m2, m3 = 1.989e30, 8.681e25, 1.024e26  # Sun, Uranus, Neptune

    # Compute distances
    r12 = np.linalg.norm(r2 - r1)
    r13 = np.linalg.norm(r3 - r1)
    r23 = np.linalg.norm(r3 - r2)

    # Compute accelerations
    a1 = G * m2 * (r2 - r1) / r12**3 + G * m3 * (r3 - r1) / r13**3
    a2 = G * m1 * (r1 - r2) / r12**3 + G * m3 * (r3 - r2) / r23**3
    a3 = G * m1 * (r1 - r3) / r13**3 + G * m2 * (r2 - r3) / r23**3

    return np.concatenate((v1, v2, v3, a1, a2, a3))

# Initial conditions
# Positions (in meters)
r1_0 = np.array([0.0, 0.0])  # Sun
r2_0 = np.array([2.871e12, 0.0])  # Uranus
r3_0 = np.array([4.495e12, 0.0])  # Neptune

# Velocities (in meters per second)
v1_0 = np.array([0.0, 0.0])  # Sun
v2_0 = np.array([0.0, 6.8e3])  # Uranus
v3_0 = np.array([0.0, 5.4e3])  # Neptune

yv0 = np.concatenate((r1_0, r2_0, r3_0, v1_0, v2_0, v3_0))

# Time range
t0 = 0
tf = 365 * 24 * 3600 * 500  # 100 years
h = 24 * 3600  # 1 day

t_values, yv_values = runge_kutta_system(three_body_problem, yv0, t0, tf, h)

# Extract positions
r1_values = np.array([yv[:2] for yv in yv_values])
r2_values = np.array([yv[2:4] for yv in yv_values])
r3_values = np.array([yv[4:6] for yv in yv_values])

# Plot the trajectories
plt.plot(r1_values[:, 0], r1_values[:, 1], label='Sun')
plt.plot(r2_values[:, 0], r2_values[:, 1], label='Uranus')
plt.plot(r3_values[:, 0], r3_values[:, 1], label='Neptune')
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.title('Three-Body Problem Trajectories: Sun, Uranus, and Neptune')
plt.legend()
plt.grid(True)
plt.show()
