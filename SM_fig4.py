import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from matplotlib.ticker import AutoMinorLocator
import matplotlib.colors as colors


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
right_margin = 0.92  # Adjust this value as needed
bottom_margin = 0.15 # Adjust this value as needed
top_margin = 0.95    # Adjust this value as needed
plt.subplots_adjust(left=left_margin, right=right_margin, bottom=bottom_margin, top=top_margin, wspace=0.25, hspace=0.30)

color1 = ['r']
marker1 = ['o']
color2 = ['b']
marker2 = ['^']

ngrid = 2000

G = 0.66
ratio = 4.3
B = G*(ratio-1)
a = 2.0
v0 = 0.01
def read_two_column_file(file_name):
    with open(file_name, 'r') as data:
        x1, x2 = [], []
        for line in data:
            p = line.split()
            x1.append(float(p[0]))
            x2.append(float(p[1]))
    return x1, x2



def integrand_chiral(q, t):
    n1 = (((((B+G)*(q**2))**2)*(omega**2-dr**2)+(omega**2+dr**2)**2)*np.cos(omega*t)-2.0*omega*dr*(((B+G)*(q**2))**2)*np.sin(omega*t))*np.exp(-dr*t)
    d1 = (((B+G)*(q**2))**2 + omega**2 - dr**2)**2 + 4.0 * (dr * omega)**2
    f1 = n1/d1
    n2 = ((((G*(q**2))**2)*(omega**2-dr**2)+(omega**2+dr**2)**2)*np.cos(omega*t)-2.0*omega*dr*((G*(q**2))**2)*np.sin(omega*t))*np.exp(-dr*t)
    d2 = ((G*(q**2))**2 + omega**2 - dr**2)**2 + 4.0 * (dr * omega)**2
    f2 = n2/d2
    n3 = (B+G)*(q**2)*dr*(dr**2 + omega**2 - ((B+G)*(q**2))**2 )*np.exp(-(B+G)*(q**2)*t)
    d3 = (dr**2 + omega**2 - ((B+G)*(q**2))**2)**2 + 4.0 * ((B+G)*(q**2) * omega)**2
    f3 = n3/d3
    n4 = G*(q**2)*dr*(dr**2 + omega**2 - (G*(q**2))**2 )*np.exp(-G*(q**2)*t)
    d4 = (dr**2 + omega**2 - (G*(q**2))**2)**2 + 4.0 * (G*(q**2) * omega)**2
    f4 = n4/d4

    integral_function = (((a**2)*v0**2)/(4.0*np.pi)) * q * (f1 + f2 - f3 - f4)
    return integral_function 

def integrand_achiral(q, t):
    n1 = (((((B+G)*(q**2))**2)*(-dr**2)+(dr**2)**2))*np.exp(-dr*t)
    d1 = (((B+G)*(q**2))**2 - dr**2)**2
    f1 = n1/d1
    n2 = ((((G*(q**2))**2)*(-dr**2)+(dr**2)**2))*np.exp(-dr*t)
    d2 = ((G*(q**2))**2 - dr**2)**2 
    f2 = n2/d2
    n3 = (B+G)*(q**2)*dr*(dr**2 - ((B+G)*(q**2))**2 )*np.exp(-(B+G)*(q**2)*t)
    d3 = (dr**2 - ((B+G)*(q**2))**2)**2
    f3 = n3/d3
    n4 = G*(q**2)*dr*(dr**2 - (G*(q**2))**2 )*np.exp(-G*(q**2)*t)
    d4 = (dr**2 - (G*(q**2))**2)**2
    f4 = n4/d4

    integral_function = (((a**2)*v0**2)/(4.0*np.pi)) * q * (f1 + f2 - f3 - f4)
    return integral_function 


# Figure (a)
dr = 0.01
omega = 0.1
# Define the time values
t_values = np.linspace(0,200,ngrid)
phia_values = np.linspace(0,1,ngrid) 
y_values = np.zeros((len(phia_values), len(t_values)))
# Perform the integration for each time value
for j, phia in enumerate(phia_values):
    initial_integral_chiral, _ = quad(integrand_chiral, 0, np.pi, args=(0,))
    initial_integral_achiral, _ = quad(integrand_achiral, 0, np.pi, args=(0,))
    for i, t in enumerate(t_values):
        result_chiral, _ = quad(integrand_chiral, 0, np.pi, args=(t,))
        result_achiral, _ = quad(integrand_achiral, 0, np.pi, args=(t,))
        y_values[j, i] = (phia * result_chiral + (1-phia) * result_achiral) / (phia * initial_integral_chiral + (1-phia) * initial_integral_achiral)

norm = colors.Normalize(vmin=-1, vmax=1)

cmesh = ax1.pcolormesh(t_values, phia_values, y_values, shading='auto', cmap=cmap, norm=norm)


# Set properties of ax1
ax1.tick_params(which='major', direction='in', bottom=True, top=True, left=True, right=True)
ax1.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True)

# Set minor tick locator
ax1.xaxis.set_minor_locator(AutoMinorLocator())
ax1.yaxis.set_minor_locator(AutoMinorLocator())

ax1.set_xlim([0.0, 200])
ax1.set_ylim([0, 1])
ax1.tick_params(axis='both', which='major', labelsize=font_size, size=3)  # Major ticks
ax1.tick_params(axis='both', which='minor', labelsize=font_size - 2, size=2)  # Minor ticks
ax1.set_xlabel('$t$', fontsize=14, labelpad=0)
ax1.set_ylabel('$\\phi_{A}/\\phi$', fontsize=14, labelpad=5)
ax1.text(-0.22, 1.05, '(a)', transform=ax1.transAxes, fontsize=16, fontweight='bold', va='top')



def integrand_chiral1(q, t):
    n1 = (((((B+G)*(q**2))**2)*(omega1**2-dr**2)+(omega1**2+dr**2)**2)*np.cos(omega1*t)-2.0*omega1*dr*(((B+G)*(q**2))**2)*np.sin(omega1*t))*np.exp(-dr*t)
    d1 = (((B+G)*(q**2))**2 + omega1**2 - dr**2)**2 + 4.0 * (dr * omega1)**2
    f1 = n1/d1
    n2 = ((((G*(q**2))**2)*(omega1**2-dr**2)+(omega1**2+dr**2)**2)*np.cos(omega1*t)-2.0*omega1*dr*((G*(q**2))**2)*np.sin(omega1*t))*np.exp(-dr*t)
    d2 = ((G*(q**2))**2 + omega1**2 - dr**2)**2 + 4.0 * (dr * omega1)**2
    f2 = n2/d2
    n3 = (B+G)*(q**2)*dr*(dr**2 + omega1**2 - ((B+G)*(q**2))**2 )*np.exp(-(B+G)*(q**2)*t)
    d3 = (dr**2 + omega1**2 - ((B+G)*(q**2))**2)**2 + 4.0 * ((B+G)*(q**2) * omega1)**2
    f3 = n3/d3
    n4 = G*(q**2)*dr*(dr**2 + omega1**2 - (G*(q**2))**2 )*np.exp(-G*(q**2)*t)
    d4 = (dr**2 + omega1**2 - (G*(q**2))**2)**2 + 4.0 * (G*(q**2) * omega1)**2
    f4 = n4/d4

    integral_function = (((a**2)*v0**2)/(4.0*np.pi)) * q * (f1 + f2 - f3 - f4)
    return integral_function 

def integrand_chiral2(q, t):
    n1 = (((((B+G)*(q**2))**2)*(omega2**2-dr**2)+(omega2**2+dr**2)**2)*np.cos(omega2*t)-2.0*omega2*dr*(((B+G)*(q**2))**2)*np.sin(omega2*t))*np.exp(-dr*t)
    d1 = (((B+G)*(q**2))**2 + omega2**2 - dr**2)**2 + 4.0 * (dr * omega2)**2
    f1 = n1/d1
    n2 = ((((G*(q**2))**2)*(omega2**2-dr**2)+(omega2**2+dr**2)**2)*np.cos(omega2*t)-2.0*omega2*dr*((G*(q**2))**2)*np.sin(omega2*t))*np.exp(-dr*t)
    d2 = ((G*(q**2))**2 + omega2**2 - dr**2)**2 + 4.0 * (dr * omega2)**2
    f2 = n2/d2
    n3 = (B+G)*(q**2)*dr*(dr**2 + omega2**2 - ((B+G)*(q**2))**2 )*np.exp(-(B+G)*(q**2)*t)
    d3 = (dr**2 + omega2**2 - ((B+G)*(q**2))**2)**2 + 4.0 * ((B+G)*(q**2) * omega2)**2
    f3 = n3/d3
    n4 = G*(q**2)*dr*(dr**2 + omega2**2 - (G*(q**2))**2 )*np.exp(-G*(q**2)*t)
    d4 = (dr**2 + omega2**2 - (G*(q**2))**2)**2 + 4.0 * (G*(q**2) * omega2)**2
    f4 = n4/d4

    integral_function = (((a**2)*v0**2)/(4.0*np.pi)) * q * (f1 + f2 - f3 - f4)
    return integral_function 


# Figure (b)
dr = 0.01
omega1 = 0.1
omega2 = 2 * omega1
# Define the time values
t_values = np.linspace(0,200,ngrid)
phia_values = np.linspace(0,1,ngrid) 
y_values = np.zeros((len(phia_values), len(t_values)))
# Perform the integration for each time value
for j, phia in enumerate(phia_values):
    initial_integral_chiral1, _ = quad(integrand_chiral1, 0, np.pi, args=(0,))
    initial_integral_chiral2, _ = quad(integrand_chiral2, 0, np.pi, args=(0,))
    for i, t in enumerate(t_values):
        result_chiral1, _ = quad(integrand_chiral1, 0, np.pi, args=(t,))
        result_chiral2, _ = quad(integrand_chiral2, 0, np.pi, args=(t,))
        y_values[j, i] = (phia * result_chiral1 + (1-phia) * result_chiral2) / (phia * initial_integral_chiral1 + (1-phia) * initial_integral_chiral2)

norm = colors.Normalize(vmin=-1, vmax=1)

cmesh = ax2.pcolormesh(t_values, phia_values, y_values, shading='auto', cmap=cmap, norm=norm)


# Set properties of ax1
ax2.tick_params(which='major', direction='in', bottom=True, top=True, left=True, right=True)
ax2.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True)

# Set minor tick locator
ax2.xaxis.set_minor_locator(AutoMinorLocator())
ax2.yaxis.set_minor_locator(AutoMinorLocator())

ax2.set_xlim([0.0, 200])
ax2.set_ylim([0, 1])
ax2.tick_params(axis='both', which='major', labelsize=font_size, size=3)  # Major ticks
ax2.tick_params(axis='both', which='minor', labelsize=font_size - 2, size=2)  # Minor ticks
ax2.set_xlabel('$t$', fontsize=14, labelpad=0)
ax2.set_ylabel('$\\phi_{A}/\\phi$', fontsize=14, labelpad=5)
ax2.text(-0.22, 1.05, '(b)', transform=ax2.transAxes, fontsize=16, fontweight='bold', va='top')


def integrand_dr1(q, t):
    n1 = (((((B+G)*(q**2))**2)*(omega**2-dr1**2)+(omega**2+dr1**2)**2)*np.cos(omega*t)-2.0*omega*dr1*(((B+G)*(q**2))**2)*np.sin(omega*t))*np.exp(-dr1*t)
    d1 = (((B+G)*(q**2))**2 + omega**2 - dr1**2)**2 + 4.0 * (dr1 * omega)**2
    f1 = n1/d1
    n2 = ((((G*(q**2))**2)*(omega**2-dr1**2)+(omega**2+dr1**2)**2)*np.cos(omega*t)-2.0*omega*dr1*((G*(q**2))**2)*np.sin(omega*t))*np.exp(-dr1*t)
    d2 = ((G*(q**2))**2 + omega**2 - dr1**2)**2 + 4.0 * (dr1 * omega)**2
    f2 = n2/d2
    n3 = (B+G)*(q**2)*dr1*(dr1**2 + omega**2 - ((B+G)*(q**2))**2 )*np.exp(-(B+G)*(q**2)*t)
    d3 = (dr1**2 + omega**2 - ((B+G)*(q**2))**2)**2 + 4.0 * ((B+G)*(q**2) * omega)**2
    f3 = n3/d3
    n4 = G*(q**2)*dr1*(dr1**2 + omega**2 - (G*(q**2))**2 )*np.exp(-G*(q**2)*t)
    d4 = (dr1**2 + omega**2 - (G*(q**2))**2)**2 + 4.0 * (G*(q**2) * omega)**2
    f4 = n4/d4

    integral_function = (((a**2)*v0**2)/(4.0*np.pi)) * q * (f1 + f2 - f3 - f4)
    return integral_function 

def integrand_dr2(q, t):
    n1 = (((((B+G)*(q**2))**2)*(omega**2-dr2**2)+(omega**2+dr2**2)**2)*np.cos(omega*t)-2.0*omega*dr2*(((B+G)*(q**2))**2)*np.sin(omega*t))*np.exp(-dr2*t)
    d1 = (((B+G)*(q**2))**2 + omega**2 - dr2**2)**2 + 4.0 * (dr2 * omega)**2
    f1 = n1/d1
    n2 = ((((G*(q**2))**2)*(omega**2-dr2**2)+(omega**2+dr2**2)**2)*np.cos(omega*t)-2.0*omega*dr2*((G*(q**2))**2)*np.sin(omega*t))*np.exp(-dr2*t)
    d2 = ((G*(q**2))**2 + omega**2 - dr2**2)**2 + 4.0 * (dr2 * omega)**2
    f2 = n2/d2
    n3 = (B+G)*(q**2)*dr2*(dr2**2 + omega**2 - ((B+G)*(q**2))**2 )*np.exp(-(B+G)*(q**2)*t)
    d3 = (dr2**2 + omega**2 - ((B+G)*(q**2))**2)**2 + 4.0 * ((B+G)*(q**2) * omega)**2
    f3 = n3/d3
    n4 = G*(q**2)*dr2*(dr2**2 + omega**2 - (G*(q**2))**2 )*np.exp(-G*(q**2)*t)
    d4 = (dr2**2 + omega**2 - (G*(q**2))**2)**2 + 4.0 * (G*(q**2) * omega)**2
    f4 = n4/d4

    integral_function = (((a**2)*v0**2)/(4.0*np.pi)) * q * (f1 + f2 - f3 - f4)
    return integral_function 


# Figure (c)
omega = 0.1
dr1 = 0.01
dr2 = 10 * dr1
# Define the time values
t_values = np.linspace(0,200,ngrid)
phia_values = np.linspace(0,1,ngrid) 
y_values = np.zeros((len(phia_values), len(t_values)))
# Perform the integration for each time value
for j, phia in enumerate(phia_values):
    initial_integral_dr1, _ = quad(integrand_chiral1, 0, np.pi, args=(0,))
    initial_integral_dr2, _ = quad(integrand_chiral2, 0, np.pi, args=(0,))
    for i, t in enumerate(t_values):
        result_dr1, _ = quad(integrand_dr1, 0, np.pi, args=(t,))
        result_dr2, _ = quad(integrand_dr2, 0, np.pi, args=(t,))
        y_values[j, i] = (phia * result_dr1 + (1-phia) * result_dr2) / (phia * initial_integral_dr1 + (1-phia) * initial_integral_dr2)

norm = colors.Normalize(vmin=-1, vmax=1)

cmesh = ax3.pcolormesh(t_values, phia_values, y_values, shading='auto', cmap=cmap, norm=norm)

cbar_ax = fig.add_axes([0.93, 0.30, 0.02, 0.50])   # Adjust the position and size as needed
cbar = fig.colorbar(cmesh, ax=ax2, cax=cbar_ax, orientation='vertical', pad=0.01)


# Set properties of ax1
ax3.tick_params(which='major', direction='in', bottom=True, top=True, left=True, right=True)
ax3.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True)

# Set minor tick locator
ax3.xaxis.set_minor_locator(AutoMinorLocator())
ax3.yaxis.set_minor_locator(AutoMinorLocator())

# Add minor ticks to the colorbar for ax2
cbar_ax3 = cbar.ax  # 'cbar' is the colorbar object for ax2
minor_locator_cbar_ax4 = AutoMinorLocator(5)
cbar_ax3.yaxis.set_minor_locator(minor_locator_cbar_ax4)
cbar.ax.tick_params(axis='y', which='minor', size=2)  # Customize minor tick parameters
cbar.set_ticks([-1, 0, 1])
cbar.ax.tick_params(direction="in", which='major', labelsize=font_size, size=4)  # Colorbar ticks
cbar.ax.tick_params(direction="in", which='minor', labelsize=font_size, size=3)  # Colorbar ticks
cbar.set_label('$\\langle \\mathbf{v}(t)\\cdot\\mathbf{v}(0) \\rangle$', fontsize=font_size, labelpad=1)

ax3.set_xlim([0.0, 200])
ax3.set_ylim([0, 1])
ax3.tick_params(axis='both', which='major', labelsize=font_size, size=3)  # Major ticks
ax3.tick_params(axis='both', which='minor', labelsize=font_size - 2, size=2)  # Minor ticks
ax3.set_xlabel('$t$', fontsize=14, labelpad=0)
ax3.set_ylabel('$\\phi_{A}/\\phi$', fontsize=14, labelpad=5)
ax3.text(-0.22, 1.05, '(c)', transform=ax3.transAxes, fontsize=16, fontweight='bold', va='top')




#plt.show()
fig.savefig("SM_fig4.png", dpi=600)   # PNG format
