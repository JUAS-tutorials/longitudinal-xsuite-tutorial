import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani

from xpart.longitudinal.rf_bucket import RFBucket

def plot_initial_distribution_and_rf_bucket(rfbucket, particle_distribution):
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    # Separatrix data
    xx = np.linspace(rfbucket.z_left, rfbucket.z_right, 1000)
    yy = rfbucket.separatrix(xx)

    # particle data
    particle_p0c = particle_distribution.p0c[0]

    offset_x_separatrix =  -rfbucket.z_sfp

    # Conversion factor the particle data
    # The x and y-axis factor are used to convert the data from (time, Energy) to (length, dp/p) scales
    # factor_x_axis = 1 / c * 1e9
    factor_x_axis = 1
    (xlabel, xlim) = ('z [m]', (-100, 100)) if factor_x_axis == 1 else ('t [ns]', (-0.2, 0.2))

    factor_y_axis = particle_p0c/1e9
    # factor_y_axis = 1
    (ylabel, ylim) = (r'\Delta p / p_0', (-2*yy.max(), 2*yy.max())) if factor_y_axis == 1 else (r'$\Delta E$ [GeV]', (-2*yy.max()*factor_y_axis, 2*yy.max()*factor_y_axis))

    # Plot the bucket separatrix (upper and lower)
    ax.plot((xx+offset_x_separatrix) * factor_x_axis,  yy * factor_y_axis, linestyle='-', color='C1')
    ax.plot((xx+offset_x_separatrix) * factor_x_axis, -yy * factor_y_axis, linestyle='-', color='C1')

    # Show lines with the RF bucket limits and center
    ax.vlines((np.array([rfbucket.z_left, rfbucket.z_right, rfbucket.z_sfp_extr])+offset_x_separatrix) * factor_x_axis,
        ylim[0], 0,
        linestyle=':', color='xkcd:dark grey')

    # Plot the particles
    ax.scatter((particle_distribution.zeta)*factor_x_axis,
               particle_distribution.delta*factor_y_axis)

    # Set the labels and limits
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)


def plot_rf_bucket_frame(line, ii_frame, gamma_value_list, p0c_value_list, pzeta_value_list, rf_bucket_properties_dict,
m_p, e, c, zeta_value_list, delta_value_list, turn_number_list):

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
m_p, e, c, zeta_value_list, delta_value_list, turn_number_list):

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
m_p, e, c, zeta_value_list, delta_value_list, turn_number_list):

    gamma = gamma_value_list[0]
    energy = gamma * m_p * c**2 / (e*1e9) 
    particle_p0c = p0c_value_list[0]

    rf_bucket_properties_dict.update({'gamma': gamma}) 

    rfbucket = RFBucket(**rf_bucket_properties_dict)
    zeta_range = np.linspace(rfbucket.z_left, rfbucket.z_right, 1000)
    upper_separatrix = rfbucket.separatrix(zeta_range)
    lower_separatrix = -rfbucket.separatrix(zeta_range)



    # factor_x_axis = 1 / c * 1e9
    factor_x_axis = 1
    (xlabel, xlim) = ('z [m]', (-100, 100)) if factor_x_axis == 1 else ('t [ns]', (-0.2, 0.2))
    factor_y_axis = particle_p0c/1e9
    (ylabel, ylim) = (r'\Delta p / p_0', (-3*upper_separatrix.max(), 3*upper_separatrix.max())) if factor_y_axis == 1 else (r'$\Delta E$ [GeV]', (-2*upper_separatrix.max()*factor_y_axis, 2*upper_separatrix.max()*factor_y_axis))

    fig, ax = plt.subplots()
    offset_x_separatrix = -rfbucket.z_sfp

    upper_separatrix_plot, = ax.plot((zeta_range+offset_x_separatrix) * factor_x_axis,   upper_separatrix * factor_y_axis, linestyle='-', color='C1')
    lower_separatrix_plot, = ax.plot((zeta_range+offset_x_separatrix) * factor_x_axis,  -upper_separatrix * factor_y_axis, linestyle='-', color='C1')

    scatter = ax.scatter((zeta_value_list[0])*factor_x_axis, delta_value_list[0] * factor_y_axis, c="b", s=20)
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
        # factor_x_axis = 1 / c * 1e9
        factor_x_axis = 1
        factor_y_axis = particle_p0c/1e9

        rf_bucket_properties_dict.update({'gamma': gamma_value_list[frame_number]}) 
        rfbucket = RFBucket(**rf_bucket_properties_dict)
        upper_separatrix = rfbucket.separatrix(zeta_range)
        lower_separatrix = -rfbucket.separatrix(zeta_range)

        upper_separatrix_plot.set_ydata(upper_separatrix*factor_y_axis)
        lower_separatrix_plot.set_ydata(lower_separatrix*factor_y_axis)

        x = np.hstack(zeta_value_list[:frame_number]) * factor_x_axis
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
m_p, e, c, zeta_value_list, delta_value_list, turn_number_list):

    gamma = gamma_value_list[0]
    energy = gamma * m_p * c**2 / (e*1e9) 
    particle_p0c = p0c_value_list[0]

    rf_bucket_properties_dict.update({'gamma': gamma})
    rfbucket = RFBucket(**rf_bucket_properties_dict)
    zeta_range = np.linspace(rfbucket.z_left, rfbucket.z_right, 1000)
    upper_separatrix = rfbucket.separatrix(zeta_range)
    lower_separatrix = -rfbucket.separatrix(zeta_range)



    # factor_x_axis = 1 / c * 1e9
    factor_x_axis = 1
    (xlabel, xlim) = ('z [m]', (-100, 100)) if factor_x_axis == 1 else ('t [ns]', (-0.2, 0.2))
    factor_y_axis = particle_p0c/1e9
    (ylabel, ylim) = (r'\Delta p / p_0', (-3*upper_separatrix.max(), 3*upper_separatrix.max())) if factor_y_axis == 1 else (r'$\Delta E$ [GeV]', (-2*upper_separatrix.max()*factor_y_axis, 2*upper_separatrix.max()*factor_y_axis))
    ylim = (-0.1, 0.1)

    fig, ax = plt.subplots()
    offset_x_separatrix = -rfbucket.z_sfp

    upper_separatrix_plot, = ax.plot((zeta_range+offset_x_separatrix) * factor_x_axis,   upper_separatrix * factor_y_axis, linestyle='-', color='C1')
    lower_separatrix_plot, = ax.plot((zeta_range+offset_x_separatrix) * factor_x_axis,  -upper_separatrix * factor_y_axis, linestyle='-', color='C1')

    scatter = ax.scatter((zeta_value_list[0])*factor_x_axis, delta_value_list[0] * factor_y_axis, c="b", s=20)
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
        # factor_x_axis = 1 / c * 1e9
        factor_x_axis = 1
        factor_y_axis = particle_p0c/1e9

        rf_bucket_properties_dict.update({'gamma': gamma_value_list[frame_number]})
        rfbucket = RFBucket(**rf_bucket_properties_dict)

        zeta_range = np.linspace(rfbucket.z_left, rfbucket.z_right, 1000)
        upper_separatrix = rfbucket.separatrix(zeta_range)
        lower_separatrix = -rfbucket.separatrix(zeta_range)

        upper_separatrix_plot.set_xdata((zeta_range+offset_x_separatrix))
        lower_separatrix_plot.set_xdata((zeta_range+offset_x_separatrix))

        upper_separatrix_plot.set_ydata(upper_separatrix*factor_y_axis)
        lower_separatrix_plot.set_ydata(lower_separatrix*factor_y_axis)

        x = zeta_value_list[frame_number] * factor_x_axis
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
m_p, e, c, zeta_value_list, delta_value_list, turn_number_list):

    gamma = gamma_value_list[0]
    energy = gamma * m_p * c**2 / (e*1e9) 
    particle_p0c = p0c_value_list[0]

    rf_bucket_properties_dict = rf_bucket_properties_dict_list[0]

    rfbucket = RFBucket(**rf_bucket_properties_dict)
    zeta_range = np.linspace(rfbucket.z_left, rfbucket.z_right, 1000)
    upper_separatrix = rfbucket.separatrix(zeta_range)
    lower_separatrix = -rfbucket.separatrix(zeta_range)



    # factor_x_axis = 1 / c * 1e9
    factor_x_axis = 1
    (xlabel, xlim) = ('z [m]', (-100, 100)) if factor_x_axis == 1 else ('t [ns]', (-0.2, 0.2))
    factor_y_axis = particle_p0c/1e9
    (ylabel, ylim) = (r'\Delta p / p_0', (-3*upper_separatrix.max(), 3*upper_separatrix.max())) if factor_y_axis == 1 else (r'$\Delta E$ [GeV]', (-2*upper_separatrix.max()*factor_y_axis, 2*upper_separatrix.max()*factor_y_axis))
    ylim = (-0.1, 0.1)

    fig, ax = plt.subplots()
    offset_x_separatrix = -rfbucket.z_sfp

    upper_separatrix_plot, = ax.plot((zeta_range+offset_x_separatrix) * factor_x_axis,   upper_separatrix * factor_y_axis, linestyle='-', color='C1')
    lower_separatrix_plot, = ax.plot((zeta_range+offset_x_separatrix) * factor_x_axis,  -upper_separatrix * factor_y_axis, linestyle='-', color='C1')

    scatter = ax.scatter((zeta_value_list[0])*factor_x_axis, delta_value_list[0] * factor_y_axis, c="b", s=20)
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
        # factor_x_axis = 1 / c * 1e9
        factor_x_axis = 1
        factor_y_axis = particle_p0c/1e9

        rf_bucket_properties_dict = rf_bucket_properties_dict_list[frame_number]
        rfbucket = RFBucket(**rf_bucket_properties_dict)

        zeta_range = np.linspace(rfbucket.z_left, rfbucket.z_right, 1000)
        upper_separatrix = rfbucket.separatrix(zeta_range)
        lower_separatrix = -rfbucket.separatrix(zeta_range)

        offset_x_separatrix = -rfbucket.z_sfp

        upper_separatrix_plot.set_xdata((zeta_range+offset_x_separatrix))
        lower_separatrix_plot.set_xdata((zeta_range+offset_x_separatrix))

        upper_separatrix_plot.set_ydata(upper_separatrix*factor_y_axis)
        lower_separatrix_plot.set_ydata(lower_separatrix*factor_y_axis)

        x = zeta_value_list[frame_number] * factor_x_axis
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