  # -*- coding: utf-8 -*-
"""
Created on Tue Sep  9 17:58:16 2025

@author: ziadm
"""

# scientific and numerical imports
from scipy.integrate import solve_ivp
import numpy as np

# matplotlib imports
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# main code
if __name__ == '__main__':
    # assign constants
    g = 9.81
    
    # error-handling functions
    def nonfloat_error_handling(variable, variable_name, units):
        while variable is not float:
            try:
                variable = float(variable)
                break
            except:
                print(f'INVALID ENTRY: {variable_name} must be a number.\n')
                variable = input(f'Re-enter the {variable_name} (in {units}): ')
        return variable
    
    def nonpositive_error_handling(variable, variable_name, units):
        while variable <= 0:
            print(f'INVALID ENTRY: {variable_name} must be greater than zero.\n')
            variable = input(f'Re-enter the {variable_name} (in {units}): ')
            variable = nonfloat_error_handling(variable, variable_name, units)
        return variable
    
    def outofrange_error_handling(variable, variable_name, units, lower_bound, upper_bound):
        while not lower_bound < variable < upper_bound:
            print(f'INVALID ENTRY: {variable_name} must be between '
                  f'{lower_bound} and {upper_bound} {units} (non-inclusive).\n')
            variable = input(f'Re-enter the {variable_name} (in {units}): ')
            variable = nonfloat_error_handling(variable, variable_name, units)
        return variable
    
    # pendulum ODE function
    def pendulum_ODE(t, y):
        return (y[1], -g*np.sin(y[0])/length)
    
    # pendulum position function
    def pendulum_position(angle):
        return (length*np.sin(angle), -length*np.cos(angle))
    
    # length input and error-handling
    length = input('Enter the length of your pendulum (in meters): ')
    length = nonfloat_error_handling(length, 'pendulum length', 'meters')
    length = nonpositive_error_handling(length, 'pendulum length', 'meters')
    
    # initial angle input and error-handling
    theta0 = input('\nEnter the initial angular position of your pendulum (in degrees): ')
    theta0 = nonfloat_error_handling(theta0, 'initial angular position', 'degrees')
    theta0 = outofrange_error_handling(theta0, 'initial angular position', 'degrees', 0, 90)
    
    # initial angular velocity input and error-handling
    theta_dot0 = input('\nEnter the initial angular velocity of your pendulum (in degrees per second): ')
    theta_dot0 = nonfloat_error_handling(theta_dot0, 'initial angular velocity', 'degrees per second')
        
    # convert deg to rad
    theta0 = np.deg2rad(theta0)
    theta_dot0 = np.deg2rad(theta_dot0)
    
    # solve pendulum ODE
    time_span = (0, 10)
    y0 = (theta0, theta_dot0)
    solution = solve_ivp(pendulum_ODE, time_span, y0,
                         t_eval=np.linspace(time_span[0],time_span[1],60*time_span[1]))
    
    '''
        Input ODE function, 1D independent variable range ("time_span"), 
    initial conditions vector ("y0": initial angular position and velocity).
    
        "solve_ivp" links "t_span" to "t" and "y0" to "y", setting the return vector
    of the ODE function (y[1], -g*np.sin(y[0])/length) to be the derivative of the
    vector (y[0], y[1]) that corresponds to "y0", all with respect to "t".
    
        Python now knows that "y[1]" is the derivative of "y[0]", and that
    "-g*np.sin(y[0])/length" is the derivative of "y[1]". It thus deduces that
    "-g*np.sin(y[0])/length" is the second derivative of "y[0]", and knows how to
    evolve the system of ODEs from the initial conditions (theta0, theta_dot0).
    
        "np.linspace(a,b,c)" tells "solve_ivp" to generate c equally spaced intervals
    from a to b, so c = fps*(b-a). As "np.linspace" is passed into the "t_eval" parameter,
    the system of ODEs iterates 600 times from 0 to 10 (60 fps). Each interval is
    appended into the array "t". Each iterated angular position and velocity is respectively
    appended into columns "y[0]" and "y[1]" of the array "y". Arrays are appended into by row.
    '''
    
    # store solution results
    time = solution.t
    theta = solution.y[0]
    theta_dot = solution.y[1]
    
    # convert rad to degrees
    theta = np.rad2deg(solution.y[0])
    theta_dot = np.rad2deg(solution.y[1])
    
    # save results to csv file
    with open('Simple_Pendulum.csv', 'w') as results_file:
        pass
        results_file.write('Angle (deg) vs Angular Velocity (deg/sec)\n\n')
        np.savetxt(results_file, np.transpose([theta, theta_dot]), delimiter=', ')
    
    # plot time domain of angular position and velocity
    plt.plot(time, theta, 'r', lw=2, label=r'$\theta$')
    plt.plot(time, theta_dot, 'b', lw=2, label=r'$\dot \theta$')
    plt.title('Simple Pendulum: Time Domain')
    plt.xlim(time_span)
    plt.ylim(-100, 100)
    plt.xlabel('Time [s]')
    plt.ylabel(r'$\theta$ [deg], $\dot \theta$ [deg/s]')
    plt.legend()
    plt.grid()
    plt.show()
    
    # plot phase diagram of angular position and velocity
    plt.plot(theta, theta_dot, 'b')
    plt.title('Simple Pendulum: Phase Diagram')
    plt.xlim(min(theta)-5, max(theta)+5)
    plt.ylim(-100, 100)
    plt.xlabel(r'$\theta$ [deg]')
    plt.ylabel(r'$\dot \theta$ [deg/s]')
    plt.grid()
    plt.show()
    
    # animate time domain of angular position and velocity
    fig, ax = plt.subplots()
    
    theta_curve, = ax.plot(time[0], theta[0], 'r')
    theta_dot_curve, = ax.plot(time[0], theta_dot[0], 'b')
    
    ax.set_title('Simple Pendulum: Time Domain')
    ax.set_xlim(time_span)
    ax.set_ylim(-100, 100)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel(r'$\theta$ [deg], $\dot \theta$ [deg/s]')
    ax.legend([r'$\theta$ [deg]', r'$\dot \theta$ [deg/s]'])
    ax.grid()
    
    def animate_time_domain(i):
        theta_curve.set_data(time[:i+1], theta[:i+1])
        theta_dot_curve.set_data(time[:i+1], theta_dot[:i+1])
    
    animation_time_domain = animation.FuncAnimation(fig, animate_time_domain,
                                                    frames=len(time))
    ffmpeg_writer = animation.FFMpegWriter(fps=60)
    animation_time_domain.save('SP_time_domain.mp4', writer=ffmpeg_writer)
    
    # animate phase diagram of angular position and velocity
    fig, ax = plt.subplots()
    
    phase_curve, = ax.plot(theta[0], theta_dot[0], 'b')
    phase_point, = ax.plot(theta[0], theta_dot[0], 'ro')
    
    ax.set_title('Simple Pendulum: Phase Diagram')
    ax.set_xlim(min(theta)-5, max(theta)+5)
    ax.set_ylim(-100, 100)
    ax.set_xlabel(r'$\theta$ [deg]')
    ax.set_ylabel(r'$\dot \theta$ [deg/s]')
    ax.grid()
    
    def animate_phase_diagram(i):
        phase_curve.set_data(theta[:i+1], theta_dot[:i+1])
        phase_point.set_data([theta[i]], [theta_dot[i]])
    
    animation_phase_diagram = animation.FuncAnimation(fig, animate_phase_diagram,
                                                      frames=len(time))
    ffmpeg_writer = animation.FFMpegWriter(fps=60)
    animation_phase_diagram.save('SP_phase_diagram.mp4', writer=ffmpeg_writer)
    
    # animate pendulum motion
    fig = plt.figure()
    ax = fig.add_subplot(aspect='equal')
    ax.set_xlim(-length, length)
    ax.set_ylim(-length*1.25, length*0.25)
    ax.grid()
    
    x0, y0 = pendulum_position(theta0)
    line, = ax.plot([0, x0], [0, y0], lw=2, c='k')
    circle = ax.add_patch(plt.Circle(pendulum_position(theta0), 0.05, fc='r', zorder=3))
    
    def animate_pendulum(i):
        x, y = pendulum_position(theta[i])
        line.set_data([0, x], [0, y])
        circle.set_center((x, y))
    
    animation_pendulum = animation.FuncAnimation(fig, animate_pendulum,
                                                 frames=len(time))
    ffmpeg_writer = animation.FFMpegWriter(fps=60)
    animation_pendulum.save('SP_animation.mp4', writer=ffmpeg_writer)