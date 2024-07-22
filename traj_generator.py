import numpy as np
import plotly.graph_objs as go



'''
def linear_motion(t, length=2, speed=1, reverse=False, y_speed = 0, z_speed=0):
    period = 2 * length / speed
    t = t % period
    if reverse:
        t = period - t
    if t < period / 2:
        x = speed * t - length / 2
    else:
        x = length / 2 - speed * (t - period / 2)
    y = y_speed * t
    z = z_speed * t
    return x, y, z



def linear_motion(t, length=2, speed=1, reverse=False, z_speed=0, y_speed=0):
    period = 2 * length / speed
    t = t % period
    direction = 1 if t < period / 2 else -1
    if reverse:
        direction *= -1
    t = t if t < period / 2 else period - t
    x = direction * (speed * t - length / 2)
    y = direction * y_speed * t
    z = direction * z_speed * t
    return x, y, z
'''

def circular_motion(t, radius=1, speed=1, clockwise=True, z_speed=0, perturbation=0):
    angle = speed * t * (1 if clockwise else -1)
    x = radius * np.cos(angle) + np.random.uniform(-perturbation, perturbation)
    y = radius * np.sin(angle) + np.random.uniform(-perturbation, perturbation)
    z = z_speed * t + np.random.uniform(-perturbation, perturbation)
    return x, y, z



def linear_motion(t, length=2, y_length=2, z_length=2, speed=1, y_speed=1, z_speed=1, reverse=False, perturbation=0):
    x_length = length
    x_speed = speed


    x_period = 2 * x_length / x_speed
    y_period = 2 * y_length / y_speed
    z_period = 2 * z_length / z_speed
    
    # Calculate x position
    x_t = t % x_period
    if x_t < x_period / 2:
        x = x_speed * x_t - x_length / 2
    else:
        x = x_length / 2 - x_speed * (x_t - x_period / 2)
    if reverse:
        x = -x
    
    # Calculate y position
    y_t = t % y_period
    if y_t < y_period / 2:
        y = y_speed * y_t - y_length / 2
    else:
        y = y_length / 2 - y_speed * (y_t - y_period / 2)
    
    # Calculate z position
    z_t = t % z_period
    if z_t < z_period / 2:
        z = z_speed * z_t - z_length / 2
    else:
        z = z_length / 2 - z_speed * (z_t - z_period / 2)
    
    # perturbations
    x += np.random.uniform(-perturbation, perturbation)
    y += np.random.uniform(-perturbation, perturbation)
    z += np.random.uniform(-perturbation, perturbation)
    
    return x, y, z


def generate_trajectories(num_particles, num_timesteps, motion_func, group_label, **kwargs):
    trajectories = []
    for p in range(num_particles):
        particle_trajectory = []
        initial_x, initial_y, initial_z = np.random.uniform(-0.5, 0.5, 3)
        for t in range(num_timesteps):
            x, y, z = motion_func(t, **kwargs)
            x += initial_x
            y += initial_y
            z += initial_z + np.random.uniform(-0.05, 0.05)
            particle_trajectory.append([x, y, z, group_label])
        trajectories.append(particle_trajectory)
    return trajectories

def generate_all_trajectories(num_particles_per_group, num_timesteps):
    all_trajectories = []
    group_label = 1
    
    for num_particles, motion_func, kwargs in num_particles_per_group:
        group_trajectories = generate_trajectories(num_particles, num_timesteps, motion_func, group_label, **kwargs)
        all_trajectories.extend(group_trajectories)
        group_label += 1
    
    return all_trajectories