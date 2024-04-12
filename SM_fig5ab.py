import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import quad
from matplotlib.ticker import AutoMinorLocator

fontsize = 18

# Set Parameters
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'Helvetica'

fig = plt.figure(figsize=(12, 3))
ax1 = plt.subplot2grid((1, 2), (0, 0))
ax2 = plt.subplot2grid((1, 2), (0, 1))
# Set tick locations and label sizes for each axis
def set_ticks_and_labelsize(ax, interval=0.5, labelsize=14):
    # Set ticks using MultipleLocator
    #ax.xaxis.set_major_locator(plt.MultipleLocator(interval))
    ax.yaxis.set_major_locator(plt.MultipleLocator(interval))
    
    # Set tick label size
    ax.tick_params(labelsize=labelsize)

# Assuming ax1, ax2, ..., ax5 are your axes objects
set_ticks_and_labelsize(ax1)
set_ticks_and_labelsize(ax2)
for ax in [ax1, ax2]:
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))



# Adjust the left, right, bottom, and top margins of the figure
left_margin = 0.1   # Adjust this value as needed
right_margin = 0.85  # Adjust this value as needed
bottom_margin = 0.08 # Adjust this value as needed
top_margin = 0.98    # Adjust this value as needed
plt.subplots_adjust(left=left_margin, right=right_margin, bottom=bottom_margin, top=top_margin, wspace=0.0, hspace=0.0)

def set_ticks_and_labelsize(ax, major_interval=100, minor_interval=10, labelsize=14):
    # Set major ticks using MultipleLocator
    major_locator = plt.MultipleLocator(major_interval)
    ax.xaxis.set_major_locator(major_locator)
    
    # Set minor ticks using MultipleLocator
    minor_locator = plt.MultipleLocator(minor_interval)
    ax.xaxis.set_minor_locator(minor_locator)

    # Set tick label size
    ax.tick_params(labelsize=labelsize)


a = 2.0
v0 = 0.01
G = 0.61
ratio = 4.33
B = G*(ratio-1)
def read_two_column_file(file_name):
    with open(file_name, 'r') as data:
        x1, x2 = [], []
        for line in data:
            p = line.split()
            x1.append(float(p[0]))
            x2.append(float(p[1]))
    return x1, x2


def integrand(q, t, dr, omega):
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

#############
#############
#############
# Figure (a)
dr = 0.01
omega = 0.1


file_path0 = f"../../uniform_chirality/autocorr/vel_autocorr/smooth/chi_0.1/dr_0.01/"
x1, x2 = read_two_column_file(file_path0 + "velocity_autocorrelation.txt")
ax1.plot(x1[::10], x2[::10], marker = '^', color='k', markersize=6, linestyle='None', label='$D_r\pm 0\% D_r$')

t_values = np.linspace(0, 200, 500) 
integral_values = []
initial_integral, _ = quad(integrand, 0, np.pi, args=(0,dr,omega))
for t in t_values:
    result, _ = quad(integrand, 0, np.pi, args=(t,dr,omega))
    normalized_result = result / initial_integral
    integral_values.append(normalized_result)
y_values = np.array(integral_values)
ax1.plot(t_values, y_values, linestyle='-', color='k')


file_path0 = f"../../binary_mixtures/heterogeneous/dr/vel_autocorr_sim/poly_0.2/"
x1, x2 = read_two_column_file(file_path0 + "velocity_autocorrelation.txt")
ax1.plot(x1[::10], x2[::10], marker = 's', color='b', markersize=6, linestyle='None', label='$D_r\pm 20\% D_r$')



df = pd.read_csv("../../binary_mixtures/heterogeneous/dr/dynamics/poly_0.2/drs_data.csv")
drs = df['Drs'].values

t_values = np.linspace(0, 200, 500) 
# Initialize an array to store summed integral values
summed_integral_values = np.zeros_like(t_values)

# Iterate over each (dr, omega) pair
for dr in drs:
    # Compute the initial integral for the (dr, omega) pair
    initial_integral, _ = quad(integrand, 0, np.pi, args=(0, dr, omega))

    # Compute the integral for each t value
    for i, t in enumerate(t_values):
        result, _ = quad(integrand, 0, np.pi, args=(t, dr, omega))
        normalized_result = result / initial_integral
        summed_integral_values[i] += normalized_result/3183.0

y_values = np.array(summed_integral_values)
np.savetxt("anlytic_vel_autocorr_poly_0.2.txt", np.column_stack((t_values, y_values)))
ax1.plot(t_values, y_values, linestyle='-', color='b')




file_path0 = f"../../binary_mixtures/heterogeneous/dr/vel_autocorr_sim/poly_0.4/"
x1, x2 = read_two_column_file(file_path0 + "velocity_autocorrelation.txt")
ax1.plot(x1[::10], x2[::10], marker = 'o', color='r', markersize=6, linestyle='None', label='$D_r\pm 40\% D_r$')


df = pd.read_csv("../../binary_mixtures/heterogeneous/dr/dynamics/poly_0.4/drs_data.csv")
drs = df['Drs'].values

t_values = np.linspace(0, 200, 500) 
# Initialize an array to store summed integral values
summed_integral_values = np.zeros_like(t_values)

# Iterate over each (dr, omega) pair
for dr in drs:
    # Compute the initial integral for the (dr, omega) pair
    initial_integral, _ = quad(integrand, 0, np.pi, args=(0, dr, omega))

    # Compute the integral for each t value
    for i, t in enumerate(t_values):
        result, _ = quad(integrand, 0, np.pi, args=(t, dr, omega))
        normalized_result = result / initial_integral
        summed_integral_values[i] += normalized_result/3183.0

y_values = np.array(summed_integral_values)
np.savetxt("anlytic_vel_autocorr_poly_0.4.txt", np.column_stack((t_values, y_values)))
ax1.plot(t_values, y_values, linestyle='-', color='r')





ax1.tick_params(which='major', direction='in', bottom=True, top=True, left=True, right=True)
ax1.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True)
ax1.set_xlim([0.0, 200.0])
ax1.set_ylim([-0.7, 1.1])
set_ticks_and_labelsize(ax1, major_interval=100, minor_interval=10)
ax1.set_xlabel('$t$', fontsize=fontsize, labelpad=-4)
ax1.set_ylabel('$\\langle\\mathbf{v}(t)\\cdot\\mathbf{v}(0)\\rangle/\\langle\\mathbf{v}(0)^2\\rangle$', fontsize=fontsize-2, labelpad=0)
ax1.legend(loc='upper right', fontsize=fontsize-7, markerfirst=True, labelspacing=0.5)
ax1.text(-0.15, 1.02, '(a)', transform=ax1.transAxes, fontsize=fontsize, fontweight='bold', va='top')




#############
#############
#############
# Figure (b)

def vel_qcorr(q, B, G, dr, omega):
    chi = dr / np.sqrt(dr**2 + omega **2)
    xil = np.sqrt((B+G)/np.sqrt(dr**2+omega**2))
    xit = np.sqrt(G/np.sqrt(dr**2+omega**2))
    A = 1.0+chi*(xil*q)**2
    B = 1.0 + 2.0 * chi*(xil*q)**2 + (xil*q)**4
    C = 1.0+chi*(xit*q)**2
    D = 1.0 + 2.0 * chi*(xit*q)**2 + (xit*q)**4
    vel_qcorr = (((A/B)+(C/D)) / 2)
    return vel_qcorr


file_path0 = f"../../uniform_chirality/vel_corr_sim/dr_0.01/chi_0.1/"
x1, x2 = read_two_column_file(file_path0 + "vel_corr.txt")
x1 = np.array(x1)
x2 = np.array(x2)
sort_indices = np.argsort(x1)
x1_sorted = x1[sort_indices]
x2_sorted = x2[sort_indices]
log_indices = np.unique(np.logspace(0, np.log10(len(x1_sorted)-1), 30, dtype=int))
x1_subsampled = x1_sorted[log_indices]
x2_subsampled = x2_sorted[log_indices]
ax2.plot(x1[::2], x2[::2], marker='^', color='k', markersize=6, linestyle='None')

q_values = np.linspace(0.01, np.pi, 500) 

y_values = np.array([vel_qcorr(q, B, G, dr, omega) for q in q_values]) 
ax2.plot(q_values, y_values, linestyle='--', color='k')

file_path1 = f"../../uniform_chirality/vel_corr_eqm/dr_0.01/chi_0.1/output_data/"
x1, x2 = read_two_column_file(file_path1 + "corr_data.txt")
ax2.plot(x1, x2, color='k', linewidth=2)


file_path0 = f"../../binary_mixtures/heterogeneous/dr/vel_corr_sim/poly_0.2/"
x1, x2 = read_two_column_file(file_path0 + "vel_corr.txt")
x1 = np.array(x1)
x2 = np.array(x2)
sort_indices = np.argsort(x1)
x1_sorted = x1[sort_indices]
x2_sorted = x2[sort_indices]
log_indices = np.unique(np.logspace(0, np.log10(len(x1_sorted)-1), 30, dtype=int))
x1_subsampled = x1_sorted[log_indices]
x2_subsampled = x2_sorted[log_indices]
ax2.plot(x1[::2], x2[::2], marker = 's', color='b', markersize=6, linestyle='None')

df = pd.read_csv("../../binary_mixtures/heterogeneous/dr/dynamics/poly_0.2/drs_data.csv")
drs = df['Drs'].values
q_values = np.linspace(0.01, np.pi, 500) 
summed_values = np.zeros_like(q_values)
for dr in drs:
    for i, q in enumerate(q_values):
            summed_values[i] += vel_qcorr(q, B, G, dr, omega)/3183.0
y_values = np.array(summed_values)
np.savetxt("anlytic_vel_corr_poly_0.2.txt", np.column_stack((q_values, y_values)))
#ax2.plot(q_values, y_values, linestyle='--', color='b')

file_path0 = f"../../binary_mixtures/heterogeneous/dr/vel_corr_sim/poly_0.4/"
x1, x2 = read_two_column_file(file_path0 + "vel_corr.txt")
x1 = np.array(x1)
x2 = np.array(x2)
sort_indices = np.argsort(x1)
x1_sorted = x1[sort_indices]
x2_sorted = x2[sort_indices]
log_indices = np.unique(np.logspace(0, np.log10(len(x1_sorted)-1), 30, dtype=int))
x1_subsampled = x1_sorted[log_indices]
x2_subsampled = x2_sorted[log_indices]
ax2.plot(x1[::2], x2[::2], marker = 'o', color='r', markersize=6, linestyle='None')

df = pd.read_csv("../../binary_mixtures/heterogeneous/dr/dynamics/poly_0.4/drs_data.csv")
drs = df['Drs'].values
q_values = np.linspace(0.01, np.pi, 500) 
summed_values = np.zeros_like(q_values)
for dr in drs:
    for i, q in enumerate(q_values):
            summed_values[i] += vel_qcorr(q, B, G, dr, omega)/3183.0
y_values = np.array(summed_values)
np.savetxt("anlytic_vel_corr_poly_0.4.txt", np.column_stack((q_values, y_values)))
#ax2.plot(q_values, y_values, linestyle='--', color='r')



#file_path0 = f"../../binary_mixtures/mixing_chiral/vel_corr_eqm/phi_{phichiral}/output_data/"
#x1, x2 = read_two_column_file(file_path0 + "corr_data.txt")
#ax2.plot(x1[::1], x2[::1], color=color, linestyle='-')

#q_values = np.linspace(min(x1), max(x1), 500) 
# Set properties of ax2
ax2.tick_params(which='major', direction='in', bottom=True, top=True, left=True, right=True)
ax2.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True)
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlim([0.1, 1.5*np.pi])
ax2.set_ylim([0.003, 1.1])
ax2.set_xlabel('$|\\mathbf{q}|$', fontsize=fontsize, labelpad=-2)
ax2.set_ylabel('$\\langle|\\mathbf{v}(\\mathbf{q})^2|\\rangle/N v_0^2$', fontsize=fontsize, labelpad=4)
#ax2.legend(loc='upper right', fontsize=12, markerfirst=True, labelspacing=0.2)
ax2.text(-0.15, 1.02, '(b)', transform=ax2.transAxes, fontsize=fontsize, fontweight='bold', va='top')


# Adjust spacing and show the plot
plt.tight_layout()
plt.show()

fig.savefig('SM_fig5ab.pdf', dpi=600)
