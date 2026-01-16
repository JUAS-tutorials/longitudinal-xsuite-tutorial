import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani

from scipy.constants import m_p, c, e

from xpart.longitudinal.rf_bucket import RFBucket

def plot_initial_distribution_and_rf_bucket(line, rfbucket, particle_distribution):
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    # Separatrix data
    xx = np.linspace(rfbucket.z_left, rfbucket.z_right, 1000)
    yy = rfbucket.separatrix(xx)

    # particle data
    particle_p0c = particle_distribution.p0c[0]

    # RF frequency
    frequency_rf = line.elements[0].frequency_rf
    
    # Particle beta
    particle_beta = line.particle_ref.beta0

    offset_x_separatrix =  -rfbucket.z_sfp
    
    # The x and y-axis factor are used to convert the data from (z, P) to (Phi, E) scales
    factor_x_axis = 1 / (particle_beta * c) * (frequency_rf * 2 * np.pi)
    (xlabel, xlim) = (r'$\phi$ [rad]', (-np.pi, 2*np.pi))
    additional_phi_offset = np.deg2rad(line.elements[0].lag_rf) # it's in degrees in Xsuite
    
    # factor_x_axis = 1
    # (xlabel, xlim) = ('z [m]', (-100, 100))

    factor_y_axis = particle_p0c/1e9
    (ylabel, ylim) = (r'$\Delta E$ [GeV]', (-2*yy.max()*factor_y_axis, 2*yy.max()*factor_y_axis))
    
    # factor_y_axis = 1
    # (ylabel, ylim) = (r'\Delta p / p_0', (-2*yy.max(), 2*yy.max()))

    # Plot the bucket separatrix (upper and lower)
    ax.plot((xx+offset_x_separatrix) * factor_x_axis + additional_phi_offset,  yy * factor_y_axis, linestyle='-', color='C1')
    ax.plot((xx+offset_x_separatrix) * factor_x_axis + additional_phi_offset, -yy * factor_y_axis, linestyle='-', color='C1')

    # Plot the particles
    ax.scatter((particle_distribution.zeta)*factor_x_axis + additional_phi_offset,
               particle_distribution.delta*factor_y_axis)

    # Set the labels and limits
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)


def plot_rf_bucket_frame(line, ii_frame, gamma_value_list, p0c_value_list, pzeta_value_list, rf_bucket_properties_dict,
zeta_value_list, delta_value_list, turn_number_list):

    gamma = gamma_value_list[ii_frame]
    energy = gamma * m_p * c**2 / (e*1e9) 
    particle_p0c = p0c_value_list[ii_frame]

    rf_bucket_properties_dict.update({'gamma': gamma})
    rfbucket = RFBucket(**rf_bucket_properties_dict)
    zeta_range = np.linspace(rfbucket.z_left, rfbucket.z_right, 2000)
    upper_separatrix = rfbucket.separatrix(zeta_range)
    lower_separatrix = -rfbucket.separatrix(zeta_range)



    # factor_x_axis = 1 / c * 1e9
    factor_x_axis = 1
    (xlabel, xlim) = ('z [m]', (-100, 100)) if factor_x_axis == 1 else ('t [ns]', (-0.2, 0.2))
    factor_y_axis = particle_p0c/1e9
    (ylabel, ylim) = (r'\Delta p / p_0', (-3*upper_separatrix.max(), 3*upper_separatrix.max())) if factor_y_axis == 1 else (r'$\Delta E$ [GeV]', (-2*upper_separatrix.max()*factor_y_axis, 2*upper_separatrix.max()*factor_y_axis))
    # ylim = (-0.2, 0.2)

    fig, ax = plt.subplots()
    offset_x_separatrix = - rfbucket.z_sfp

    upper_separatrix_plot, = ax.plot((zeta_range+offset_x_separatrix) * factor_x_axis,   upper_separatrix * factor_y_axis, linestyle='-', color='C1')
    lower_separatrix_plot, = ax.plot((zeta_range+offset_x_separatrix) * factor_x_axis,  -upper_separatrix * factor_y_axis, linestyle='-', color='C1')

    offset_x_mesh = 0

    z_mesh, dp_mesh = np.meshgrid(np.linspace(-100-offset_x_mesh, 100-offset_x_mesh, num=100), np.linspace(-2*0.005, 2*0.005, num=100))
    hamiltonian_contour = rfbucket.hamiltonian(z_mesh, dp_mesh, make_convex=False)

    while len(ax.collections) > 1:
        ax.collections[1].remove()

    # contour_phase_space = ax.contour((z_mesh+offset_x_mesh)*factor_x_axis, dp_mesh*factor_y_axis, hamiltonian_contour, 15, cmap=plt.cm.viridis_r)

    ax.set_title(f'Turn {turn_number_list[ii_frame]}\n $\gamma$ {rfbucket.gamma:.1f} $E_{{kin}}$ {(rfbucket.gamma-1)*m_p*c**2/(e*1e9):.2f} GeV')

    offset_x_particles = 0
    scatter = ax.scatter((zeta_value_list[ii_frame]+offset_x_particles)*factor_x_axis, delta_value_list[ii_frame] * factor_y_axis, c="b", s=20)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)


def plot_rf_bucket_frame_acceleration(line, ii_frame, gamma_value_list, p0c_value_list, pzeta_value_list, rf_bucket_properties_dict_list,
zeta_value_list, delta_value_list, turn_number_list):

    gamma = gamma_value_list[ii_frame]
    energy = gamma * m_p * c**2 / (e*1e9) 
    particle_p0c = p0c_value_list[ii_frame]

    rf_bucket_properties_dict = rf_bucket_properties_dict_list[ii_frame]

    rfbucket = RFBucket(**rf_bucket_properties_dict)
    zeta_range = np.linspace(rfbucket.z_left, rfbucket.z_right, 2000)
    upper_separatrix = rfbucket.separatrix(zeta_range)
    lower_separatrix = -rfbucket.separatrix(zeta_range)



    # factor_x_axis = 1 / c * 1e9
    factor_x_axis = 1
    (xlabel, xlim) = ('z [m]', (-100, 100)) if factor_x_axis == 1 else ('t [ns]', (-0.2, 0.2))
    factor_y_axis = particle_p0c/1e9
    (ylabel, ylim) = (r'\Delta p / p_0', (-3*upper_separatrix.max(), 3*upper_separatrix.max())) if factor_y_axis == 1 else (r'$\Delta E$ [GeV]', (-2*upper_separatrix.max()*factor_y_axis, 2*upper_separatrix.max()*factor_y_axis))
    ylim = (-0.2, 0.2)

    fig, ax = plt.subplots()
    offset_x_separatrix = - rfbucket.z_sfp

    upper_separatrix_plot, = ax.plot((zeta_range+offset_x_separatrix) * factor_x_axis,   upper_separatrix * factor_y_axis, linestyle='-', color='C1')
    lower_separatrix_plot, = ax.plot((zeta_range+offset_x_separatrix) * factor_x_axis,  -upper_separatrix * factor_y_axis, linestyle='-', color='C1')

    offset_x_mesh = 0

    z_mesh, dp_mesh = np.meshgrid(np.linspace(-100-offset_x_mesh, 100-offset_x_mesh, num=100), np.linspace(-2*0.005, 2*0.005, num=100))
    hamiltonian_contour = rfbucket.hamiltonian(z_mesh, dp_mesh, make_convex=False)

    while len(ax.collections) > 1:
        ax.collections[1].remove()

    # contour_phase_space = ax.contour((z_mesh+offset_x_mesh)*factor_x_axis, dp_mesh*factor_y_axis, hamiltonian_contour, 15, cmap=plt.cm.viridis_r)

    ax.set_title(f'Turn {turn_number_list[ii_frame]}\n $\gamma$ {rfbucket.gamma:.1f} $E_{{kin}}$ {(rfbucket.gamma-1)*m_p*c**2/(e*1e9):.2f} GeV')

    offset_x_particles = 0
    scatter = ax.scatter((zeta_value_list[ii_frame]+offset_x_particles)*factor_x_axis, delta_value_list[ii_frame] * factor_y_axis, c="b", s=20)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)


def plot_rf_bucket_animation_persistent(line, gamma_value_list, p0c_value_list, pzeta_value_list, rf_bucket_properties_dict,
zeta_value_list, delta_value_list, turn_number_list):

    gamma = gamma_value_list[0]
    energy = gamma * m_p * c**2 / (e*1e9) 
    particle_p0c = p0c_value_list[0]
    particle_beta = np.sqrt(1-1/(gamma**2))
    frequency_rf = line.elements[0].frequency_rf

    rf_bucket_properties_dict.update({'gamma': gamma}) 

    rfbucket = RFBucket(**rf_bucket_properties_dict)
    zeta_range = np.linspace(rfbucket.z_left, rfbucket.z_right, 1000)
    upper_separatrix = rfbucket.separatrix(zeta_range)
    lower_separatrix = -rfbucket.separatrix(zeta_range)

    factor_x_axis = (1 / (particle_beta * c)) * (2*np.pi*frequency_rf)
    (xlabel, xlim) = ('Phi [rad]', (-np.pi, 2*np.pi))
    additional_phi_offset = np.deg2rad(line.elements[0].lag_rf) # it's in degrees in Xsuite

    # factor_x_axis = 1
    # (xlabel, xlim) = ('z [m]', (-100, 100))
    # additional_phi_offset = 0
    
    factor_y_axis = particle_beta*particle_p0c/1e9
    (ylabel, ylim) = (r'$\Delta E$ [GeV]', (-2*upper_separatrix.max()*factor_y_axis, 2*upper_separatrix.max()*factor_y_axis))
    
    # factor_y_axis = 1
    # (ylabel, ylim) = (r'\Delta p / p_0', (-3*upper_separatrix.max(), 3*upper_separatrix.max()))

    fig, ax = plt.subplots()
    offset_x_separatrix = -rfbucket.z_sfp

    upper_separatrix_plot, = ax.plot((zeta_range+offset_x_separatrix) * factor_x_axis + additional_phi_offset,   upper_separatrix * factor_y_axis, linestyle='-', color='C1')
    lower_separatrix_plot, = ax.plot((zeta_range+offset_x_separatrix) * factor_x_axis + additional_phi_offset,  -upper_separatrix * factor_y_axis, linestyle='-', color='C1')

    scatter = ax.scatter((zeta_value_list[0])*factor_x_axis + additional_phi_offset, delta_value_list[0] * factor_y_axis, c="b", s=20)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plot_title = ax.set_title('Turn 0')

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)


    def update_frame_number():
        frame_number = 0
        frame_number_max = len(zeta_value_list)
        while frame_number < frame_number_max:
            frame_number += anim.direction
            yield frame_number

    def update_plot(frame_number):
        particle_p0c = p0c_value_list[frame_number]
        particle_beta = np.sqrt(1-1/(gamma_value_list[frame_number]**2))
        factor_x_axis = 1 / (particle_beta * c) * (2*np.pi*frequency_rf)
        additional_phi_offset = np.deg2rad(line.elements[0].lag_rf) # it's in degrees in Xsuite
        factor_y_axis = particle_p0c/1e9

        rf_bucket_properties_dict.update({'gamma': gamma_value_list[frame_number]}) 
        rfbucket = RFBucket(**rf_bucket_properties_dict)
        upper_separatrix = rfbucket.separatrix(zeta_range)
        lower_separatrix = -rfbucket.separatrix(zeta_range)

        upper_separatrix_plot.set_ydata(upper_separatrix*factor_y_axis)
        lower_separatrix_plot.set_ydata(lower_separatrix*factor_y_axis)

        x = np.hstack(zeta_value_list[:frame_number]) * factor_x_axis + additional_phi_offset
        y = np.hstack(pzeta_value_list[:frame_number]) * factor_y_axis
        data = np.stack([x, y]).T
        scatter.set_offsets(data)

        offset_x_mesh = 0

        z_mesh, dp_mesh = np.meshgrid(np.linspace(-100-offset_x_mesh, 100-offset_x_mesh, num=100), np.linspace(-2*0.005, 2*0.005, num=100))
        hamiltonian_contour = rfbucket.hamiltonian(z_mesh, dp_mesh, make_convex=False)

        while len(ax.collections) > 1:
            ax.collections[1].remove()

        # contour_phase_space = ax.contour((z_mesh+offset_x_mesh)*factor_x_axis, dp_mesh*factor_y_axis, hamiltonian_contour, 15, cmap=plt.cm.viridis_r)

        plot_title.set_text(f'Turn {turn_number_list[frame_number]}\n $\gamma$ {rfbucket.gamma:.1f} $E_{{kin}}$ {(rfbucket.gamma-1)*m_p*c**2/(e*1e9):.2f} GeV')

    def on_press(event):
        if event.key.isspace():
            if anim.running:
                anim.event_source.stop()
            else:
                anim.event_source.start()
            anim.running = not anim.running
        elif event.key == 'left':
            anim.direction = -1
        elif event.key == 'right':
            anim.direction = +1
        if event.key in ['left', 'right']:
            t = anim.frame_seq.__next__()
            update_plot(t)
            plt.draw()

    fig.canvas.mpl_connect('key_press_event', on_press)
    anim = ani.FuncAnimation(fig, update_plot, frames=update_frame_number, interval=200, repeat=True)
    anim.running = True
    anim.direction = +1
    plt.show()


def plot_rf_bucket_animation(line, gamma_value_list, p0c_value_list, pzeta_value_list, rf_bucket_properties_dict,
zeta_value_list, delta_value_list, turn_number_list):

    gamma = gamma_value_list[0]
    energy = gamma * m_p * c**2 / (e*1e9) 
    particle_p0c = p0c_value_list[0]
    particle_beta = np.sqrt(1-1/(gamma**2))
    frequency_rf = line.elements[0].frequency_rf

    rf_bucket_properties_dict.update({'gamma': gamma})
    rfbucket = RFBucket(**rf_bucket_properties_dict)
    zeta_range = np.linspace(rfbucket.z_left, rfbucket.z_right, 1000)
    upper_separatrix = rfbucket.separatrix(zeta_range)
    lower_separatrix = -rfbucket.separatrix(zeta_range)
    
    factor_x_axis = 1 / (particle_beta * c) * (2*np.pi*frequency_rf)
    (xlabel, xlim) = ('Phi [rad]', (-np.pi, 2*np.pi))
    additional_phi_offset = np.deg2rad(line.elements[0].lag_rf) # it's in degrees in Xsuite
    
    # factor_x_axis = 1
    # (xlabel, xlim) = ('z [m]', (-100, 100))
    # additional_phi_offset = 0
    
    factor_y_axis = particle_beta * particle_p0c/1e9
    (ylabel, ylim) = (r'\Delta p / p_0', (-3*upper_separatrix.max(), 3*upper_separatrix.max())) if factor_y_axis == 1 else (r'$\Delta E$ [GeV]', (-2*upper_separatrix.max()*factor_y_axis, 2*upper_separatrix.max()*factor_y_axis))
    ylim = (-0.1, 0.1)

    fig, ax = plt.subplots()
    offset_x_separatrix = -rfbucket.z_sfp

    upper_separatrix_plot, = ax.plot((zeta_range+offset_x_separatrix) * factor_x_axis + additional_phi_offset,   upper_separatrix * factor_y_axis, linestyle='-', color='C1')
    lower_separatrix_plot, = ax.plot((zeta_range+offset_x_separatrix) * factor_x_axis + additional_phi_offset,  -upper_separatrix * factor_y_axis, linestyle='-', color='C1')

    scatter = ax.scatter((zeta_value_list[0])*factor_x_axis + additional_phi_offset, delta_value_list[0] * factor_y_axis, c="b", s=20)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plot_title = ax.set_title('Turn 0')

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    def update_frame_number():
        frame_number = 0
        frame_number_max = len(zeta_value_list)-1
        while frame_number < frame_number_max:
            frame_number += anim.direction
            yield frame_number

    def update_plot(frame_number):
        particle_p0c = p0c_value_list[frame_number]
        particle_beta = np.sqrt(1-1/(gamma_value_list[frame_number]**2))
        factor_x_axis = 1 / (particle_beta * c) * (2*np.pi*frequency_rf)
        additional_phi_offset = np.deg2rad(line.elements[0].lag_rf) # it's in degrees in Xsuite
        factor_y_axis = particle_beta*particle_p0c/1e9

        rf_bucket_properties_dict.update({'gamma': gamma_value_list[frame_number]})
        rfbucket = RFBucket(**rf_bucket_properties_dict)

        zeta_range = np.linspace(rfbucket.z_left, rfbucket.z_right, 1000)
        upper_separatrix = rfbucket.separatrix(zeta_range)
        lower_separatrix = -rfbucket.separatrix(zeta_range)

        upper_separatrix_plot.set_xdata((zeta_range+offset_x_separatrix)*factor_x_axis + additional_phi_offset)
        lower_separatrix_plot.set_xdata((zeta_range+offset_x_separatrix)*factor_x_axis + additional_phi_offset)

        upper_separatrix_plot.set_ydata(upper_separatrix*factor_y_axis)
        lower_separatrix_plot.set_ydata(lower_separatrix*factor_y_axis)

        x = zeta_value_list[frame_number] * factor_x_axis + additional_phi_offset
        y = pzeta_value_list[frame_number] * factor_y_axis
        data = np.stack([x, y]).T
        scatter.set_offsets(data)

        offset_x_mesh = 0

        z_mesh, dp_mesh = np.meshgrid(np.linspace(-100-offset_x_mesh, 100-offset_x_mesh, num=100), np.linspace(-2*0.005, 2*0.005, num=100))
        hamiltonian_contour = rfbucket.hamiltonian(z_mesh, dp_mesh, make_convex=False)

        while len(ax.collections) > 1:
            ax.collections[1].remove()

        # contour_phase_space = ax.contour((z_mesh+offset_x_mesh)*factor_x_axis, dp_mesh*factor_y_axis, hamiltonian_contour, 15, cmap=plt.cm.viridis_r)

        plot_title.set_text(f'Turn {turn_number_list[frame_number]}\n $\gamma$ {rfbucket.gamma:.1f} $E_{{kin}}$ {(rfbucket.gamma-1)*m_p*c**2/(e*1e9):.2f} GeV')

    def on_press(event):
        if event.key.isspace():
            if anim.running:
                anim.event_source.stop()
            else:
                anim.event_source.start()
            anim.running = not anim.running
        elif event.key == 'left':
            anim.direction = -1
        elif event.key == 'right':
            anim.direction = +1
        if event.key in ['left', 'right']:
            t = anim.frame_seq.__next__()
            update_plot(t)
            plt.draw()

    fig.canvas.mpl_connect('key_press_event', on_press)
    anim = ani.FuncAnimation(fig, update_plot, frames=update_frame_number, interval=200, repeat=True)
    anim.running = True
    anim.direction = +1
    plt.show()


def plot_rf_bucket_animation_acceleration(line, gamma_value_list, p0c_value_list, pzeta_value_list, rf_bucket_properties_dict_list,
zeta_value_list, delta_value_list, turn_number_list):

    gamma = gamma_value_list[0]
    energy = gamma * m_p * c**2 / (e*1e9) 
    particle_p0c = p0c_value_list[0]
    particle_beta = np.sqrt(1-1/(gamma**2))
    frequency_rf = line.elements[0].frequency_rf

    rf_bucket_properties_dict = rf_bucket_properties_dict_list[0]

    rfbucket = RFBucket(**rf_bucket_properties_dict)
    zeta_range = np.linspace(rfbucket.z_left, rfbucket.z_right, 1000)
    upper_separatrix = rfbucket.separatrix(zeta_range)
    lower_separatrix = -rfbucket.separatrix(zeta_range)

    factor_x_axis = 1 / (particle_beta * c) * (2*np.pi*frequency_rf)
    (xlabel, xlim) = ('Phi [rad]', (-np.pi, 2*np.pi))
    additional_phi_offset = np.deg2rad(rf_bucket_properties_dict['phi_offset_list'])[0] # it's in degrees in Xsuite
    
    # factor_x_axis = 1
    # (xlabel, xlim) = ('z [m]', (-100, 100)) 
    
    factor_y_axis = particle_p0c/1e9
    (ylabel, ylim) = (r'$\Delta E$ [GeV]', (-2*upper_separatrix.max()*factor_y_axis, 2*upper_separatrix.max()*factor_y_axis))
    ylim = (-0.1, 0.1)
    
    # factor_y_axis = 1
    # (ylabel, ylim) = (r'\Delta p / p_0', (-3*upper_separatrix.max(), 3*upper_separatrix.max()))

    fig, ax = plt.subplots()
    offset_x_separatrix = -rfbucket.z_sfp

    upper_separatrix_plot, = ax.plot((zeta_range+offset_x_separatrix) * factor_x_axis+additional_phi_offset,   upper_separatrix * factor_y_axis, linestyle='-', color='C1')
    lower_separatrix_plot, = ax.plot((zeta_range+offset_x_separatrix) * factor_x_axis+additional_phi_offset,  -upper_separatrix * factor_y_axis, linestyle='-', color='C1')

    scatter = ax.scatter((zeta_value_list[0])*factor_x_axis+additional_phi_offset, delta_value_list[0] * factor_y_axis, c="b", s=20)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plot_title = ax.set_title('Turn 0')

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)


    def update_frame_number():
        frame_number = 0
        frame_number_max = len(zeta_value_list)-1
        while frame_number < frame_number_max:
            frame_number += anim.direction
            yield frame_number

    def update_plot(frame_number):
        particle_p0c = p0c_value_list[frame_number]
        particle_beta = np.sqrt(1-1/(gamma_value_list[frame_number]**2))
        
        factor_x_axis = 1 / (particle_beta * c) * (2*np.pi*frequency_rf)
        additional_phi_offset = np.deg2rad(rf_bucket_properties_dict_list[frame_number]['phi_offset_list'])[0] # it's in degrees in Xsuite
        factor_y_axis = particle_beta*particle_p0c/1e9

        rf_bucket_properties_dict = rf_bucket_properties_dict_list[frame_number]
        rfbucket = RFBucket(**rf_bucket_properties_dict)

        zeta_range = np.linspace(rfbucket.z_left, rfbucket.z_right, 1000)
        upper_separatrix = rfbucket.separatrix(zeta_range)
        lower_separatrix = -rfbucket.separatrix(zeta_range)

        offset_x_separatrix = -rfbucket.z_sfp

        upper_separatrix_plot.set_xdata((zeta_range+offset_x_separatrix)*factor_x_axis + additional_phi_offset)
        lower_separatrix_plot.set_xdata((zeta_range+offset_x_separatrix)*factor_x_axis + additional_phi_offset)

        upper_separatrix_plot.set_ydata(upper_separatrix*factor_y_axis)
        lower_separatrix_plot.set_ydata(lower_separatrix*factor_y_axis)

        x = zeta_value_list[frame_number] * factor_x_axis + additional_phi_offset
        y = pzeta_value_list[frame_number] * factor_y_axis
        data = np.stack([x, y]).T
        scatter.set_offsets(data)

        offset_x_mesh = 0

        z_mesh, dp_mesh = np.meshgrid(np.linspace(-100-offset_x_mesh, 100-offset_x_mesh, num=100), np.linspace(-2*0.005, 2*0.005, num=100))
        hamiltonian_contour = rfbucket.hamiltonian(z_mesh, dp_mesh, make_convex=False)

        while len(ax.collections) > 1:
            ax.collections[1].remove()

        # contour_phase_space = ax.contour((z_mesh+offset_x_mesh)*factor_x_axis, dp_mesh*factor_y_axis, hamiltonian_contour, 15, cmap=plt.cm.viridis_r)

        plot_title.set_text(f'Turn {turn_number_list[frame_number]}\n $\gamma$ {rfbucket.gamma:.1f} $E_{{kin}}$ {(rfbucket.gamma-1)*m_p*c**2/(e*1e9):.2f} GeV')

    def on_press(event):
        if event.key.isspace():
            if anim.running:
                anim.event_source.stop()
            else:
                anim.event_source.start()
            anim.running = not anim.running
        elif event.key == 'left':
            anim.direction = -1
        elif event.key == 'right':
            anim.direction = +1
        if event.key in ['left', 'right']:
            t = anim.frame_seq.__next__()
            update_plot(t)
            plt.draw()

    fig.canvas.mpl_connect('key_press_event', on_press)
    anim = ani.FuncAnimation(fig, update_plot, frames=update_frame_number, interval=200, repeat=True)
    anim.running = True
    anim.direction = +1
    plt.show()
