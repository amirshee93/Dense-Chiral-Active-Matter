import matplotlib.pyplot as plt
import numpy as np

# Set Parameters
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'Helvetica'

font_size = 16
cmap = plt.cm.plasma # plt.cm.viridis

fig = plt.figure(figsize=(12, 3))
ax1 = plt.subplot2grid((1, 3), (0, 0))
ax2 = plt.subplot2grid((1, 3), (0, 1))
ax3 = plt.subplot2grid((1, 3), (0, 2))



# Adjust the left, right, bottom, and top margins of the figure
left_margin = 0.06   # Adjust this value as needed
right_margin = 0.90  # Adjust this value as needed
bottom_margin = 0.15 # Adjust this value as needed
top_margin = 0.95    # Adjust this value as needed
plt.subplots_adjust(left=left_margin, right=right_margin, bottom=bottom_margin, top=top_margin, wspace=0.80, hspace=0.0)

color1 = ['r']
marker1 = ['o']
color2 = ['b']
marker2 = ['^']


v0 = 0.01
def read_two_column_file(file_name):
    with open(file_name, 'r') as data:
        x1, x2 = [], []
        for line in data:
            p = line.split()
            x1.append(float(p[0]))
            x2.append(float(p[1]))
    return x1, x2


def potential_energy(la_nu, dr, omega):
    potential_energy_function = (v0**2)*(dr+la_nu)/(4.0*((dr+la_nu)**2 +omega**2))
    return potential_energy_function


# Figure (a)
omega = 0.1
# Define the time values
la_nu_values = np.linspace(0,1,1000)
dr_values = np.logspace(np.log10(0.001), np.log10(100), 1000) #np.linspace(0,1,100)
y_values = np.zeros((len(dr_values), len(la_nu_values)))
# Perform the integration for each time value
for j, dr in enumerate(dr_values):
        for i, la_nu in enumerate(la_nu_values):
            y_values[j, i] = potential_energy(la_nu,dr,omega)/v0**2

cmesh = ax1.pcolormesh(la_nu_values, dr_values, y_values, shading='auto', cmap=cmap)
cbar_ax = fig.add_axes([0.25, 0.30, 0.02, 0.50])   # Adjust the position and size as needed
cbar = fig.colorbar(cmesh, ax=ax1, cax=cbar_ax, orientation='vertical', pad=0.01)
cbar.ax.tick_params(direction="in", which='major', labelsize=font_size, size=4)  # Colorbar ticks
cbar.ax.tick_params(direction="in", which='minor', labelsize=font_size, size=3)  # Colorbar ticks
cbar.set_label('$E_{\\nu}/v_0^2$', fontsize=font_size, labelpad=0, color='k')


# Set properties of ax1
ax1.tick_params(which='major', direction='in', bottom=True, top=True, left=True, right=True)
ax1.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True)
ax1.set_yscale('log')
ax1.set_xlim([0.0, 1])
ax1.set_ylim([0.001, 100])
ax1.tick_params(axis='both', which='major', labelsize=font_size, size=3)  # Major ticks
ax1.tick_params(axis='both', which='minor', labelsize=font_size - 2, size=2)  # Minor ticks
ax1.set_xlabel('$\\lambda_{\\nu}$', fontsize=14, labelpad=0)
ax1.set_ylabel('$D_r$', fontsize=14, labelpad=0)
ax1.text(-0.30, 1.05, '(a)', transform=ax1.transAxes, fontsize=16, fontweight='bold', va='top')




# Figure (b)
omega = 1
# Define the time values
la_nu_values = np.linspace(0,3,1000)
dr_values = np.logspace(np.log10(0.001), np.log10(100), 1000) #np.linspace(0,1,100)
y_values = np.zeros((len(dr_values), len(la_nu_values)))
# Perform the integration for each time value
for j, dr in enumerate(dr_values):
        for i, la_nu in enumerate(la_nu_values):
            y_values[j, i] = potential_energy(la_nu,dr,omega)/v0**2

cmesh = ax2.pcolormesh(la_nu_values, dr_values, y_values, shading='auto', cmap=cmap)
cbar_ax = fig.add_axes([0.58, 0.30, 0.02, 0.50])   # Adjust the position and size as needed
cbar = fig.colorbar(cmesh, ax=ax1, cax=cbar_ax, orientation='vertical', pad=0.01)
cbar.ax.tick_params(direction="in", which='major', labelsize=font_size, size=4)  # Colorbar ticks
cbar.ax.tick_params(direction="in", which='minor', labelsize=font_size, size=3)  # Colorbar ticks
cbar.set_label('$E_{\\nu}/v_0^2$', fontsize=font_size, labelpad=0, color='k')


# Set properties of ax2
ax2.tick_params(which='major', direction='in', bottom=True, top=True, left=True, right=True)
ax2.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True)
ax2.set_yscale('log')
ax2.set_xlim([0.0, 3])
ax2.set_ylim([0.001, 100])
ax2.tick_params(axis='both', which='major', labelsize=font_size, size=3)  # Major ticks
ax2.tick_params(axis='both', which='minor', labelsize=font_size - 2, size=2)  # Minor ticks
ax2.set_xlabel('$\\lambda_{\\nu}$', fontsize=14, labelpad=0)
ax2.set_ylabel('$D_r$', fontsize=14, labelpad=0)
ax2.text(-0.30, 1.05, '(b)', transform=ax2.transAxes, fontsize=16, fontweight='bold', va='top')


# Figure (a)
omega = 10
# Define the time values
la_nu_values = np.linspace(0,6,1000)
dr_values = np.logspace(np.log10(0.001), np.log10(100), 1000) #np.linspace(0,1,100)
y_values = np.zeros((len(la_nu_values), len(dr_values)))
# Perform the integration for each time value
for j, dr in enumerate(dr_values):
        for i, la_nu in enumerate(la_nu_values):
            y_values[j, i] = potential_energy(la_nu,dr,omega)/v0**2

cmesh = ax3.pcolormesh(la_nu_values, dr_values, y_values, shading='auto', cmap=cmap)

cbar_ax = fig.add_axes([0.91, 0.30, 0.02, 0.50])   # Adjust the position and size as needed
cbar = fig.colorbar(cmesh, ax=ax1, cax=cbar_ax, orientation='vertical', pad=0.01)
cbar.ax.tick_params(direction="in", which='major', labelsize=font_size, size=4)  # Colorbar ticks
cbar.ax.tick_params(direction="in", which='minor', labelsize=font_size, size=3)  # Colorbar ticks
cbar.set_label('$E_{\\nu}/v_0^2$', fontsize=font_size, labelpad=0, color='k')

# Set properties of ax1
ax3.tick_params(which='major', direction='in', bottom=True, top=True, left=True, right=True)
ax3.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True)
ax3.set_yscale('log')
ax3.set_xlim([0.0, 6])
ax3.set_ylim([0.001, 100])
ax3.tick_params(axis='both', which='major', labelsize=font_size, size=3)  # Major ticks
ax3.tick_params(axis='both', which='minor', labelsize=font_size - 2, size=2)  # Minor ticks
ax3.set_xlabel('$\\lambda_{\\nu}$', fontsize=14, labelpad=0)
ax3.set_ylabel('$D_r$', fontsize=14, labelpad=0)
ax3.text(-0.30, 1.05, '(c)', transform=ax3.transAxes, fontsize=16, fontweight='bold', va='top')




plt.show()

fig.savefig('SM_fig2.png', dpi=300)
