import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle

v0 =0.01
# Function to convert an angle to a color using HSV colormap for subfigure (a)
def vector_to_rgb_a(angle):
    # Ensure the angle is within 0 to 2π
    angle = angle % (2 * np.pi)
    
    # Normalize the angle from 0 to 2π to a range between 0 and 1
    normalized_angle = angle / (2 * np.pi)
    return plt.cm.twilight(normalized_angle)

# Function to convert an angle to a color using Viridis colormap for subfigure (b)
def vector_to_rgb_b(angle):
    return plt.cm.viridis(angle)

plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'Helvetica'

# Master figure
fig = plt.figure(figsize=(15, 5))

# Subfigure (a)
left_margin_a = 0.04   # Adjust this value as needed
right_margin_a = 0.435  # Adjust this value as needed
bottom_margin_a = 0.1 # Adjust this value as needed
top_margin_a = 0.98    # Adjust this value as needed
gs_a = plt.GridSpec(3, 4, left=left_margin_a, right=right_margin_a, bottom=bottom_margin_a, top=top_margin_a, wspace=0.0, hspace=0.04)
all_axes_a = []
omega_values_a = [10, 5, 1]
dr_values_a = [0.1, 1.0, 10.0, 100.0]
font_size = 18
num_rows_a = 3
num_cols_a = 4

all_axes_a = []
for i in range(num_rows_a):
    for j in range(num_cols_a):
        ax = plt.subplot(gs_a[i, j])
        ax.set_aspect('equal', 'box')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlim([-52, 52])
        ax.set_ylim([-52, 52])
        all_axes_a.append(ax)

for ax in all_axes_a:
    for spine in ax.spines.values():
        spine.set_linewidth(0)


# Annotate the top row
for i, ax in enumerate(all_axes_a[-num_cols_a:]):
    ax.text(0.5, -0.1, str(dr_values_a[i]), transform=ax.transAxes,
            ha='center', va='center', fontsize=font_size)
# Annotate the left column
for i, ax in enumerate(all_axes_a[::num_cols_a]):
    ax.text(-0.05, 0.5, str(omega_values_a[i]), transform=ax.transAxes, ha='right', va='center', fontsize=font_size)



center_subplot_col = all_axes_a[num_cols_a // 2]
center_subplot_row = all_axes_a[num_rows_a // 2]
center_subplot_row.text(-1.20, -0.50, "$\\Omega$", transform=center_subplot_row.transAxes,
                    ha='right', va='center', fontsize=font_size, rotation='vertical')

center_subplot_col.text(0.0, -2.30, "$D_r$", transform=center_subplot_col.transAxes,
                    ha='right', va='center', fontsize=font_size)

center_subplot_row.text(-1.150, 0.95, "(a)", transform=center_subplot_row.transAxes,
                    ha='right', va='center', fontsize=font_size)


# Create a colorbar reference image
data = np.linspace(0, 1, 256).reshape(1, -1)
data = np.vstack((data, data))

cbar_ax_position = [right_margin_a+0.005, bottom_margin_a+0.2, 0.02, top_margin_a-bottom_margin_a-0.4]
ax = fig.add_axes([0, 0, 0.01, 0.01], visible=False)
cax = ax.imshow(data, cmap=plt.cm.twilight)
cbar_ax = fig.add_axes(cbar_ax_position)
cbar = fig.colorbar(cax, cax=cbar_ax, orientation='vertical', ticks=[0, 0.5, 1])
cbar.ax.set_yticklabels(['$0$', '$\\pi$', '$2\\pi$'], fontsize=font_size)
cbar.ax.tick_params(direction='in', length=6)  # Set tick direction inward and length to 6
cbar.set_label(r'$\phi_{\mathbf{v}}$', fontsize=font_size, labelpad=-10)


for omega in omega_values_a:
    for dr in dr_values_a:
        subdir = f'omega_{omega}/dr_{dr}/dynamics_data/positions_1999.csv'
        df = pd.read_csv(subdir)
        xval = df['xval'].values
        yval = df['yval'].values
        vxval = df['vxval'].values
        vyval = df['vyval'].values
        angles = np.arctan2(vyval, vxval)

        col = dr_values_a.index(dr)
        row = omega_values_a.index(omega)
        ax = all_axes_a[row * num_cols_a + col]

        for i in range(len(xval)):
            ax.add_patch(Circle((xval[i], yval[i]), 1.0, color=vector_to_rgb_a(angles[i]), fill=True, lw=0.0))

# Subfigure (b)
left_margin_b = 0.53   # Adjust this value as needed
right_margin_b = 0.925  # Adjust this value as needed
bottom_margin_b = 0.1 # Adjust this value as needed
top_margin_b = 0.98    # Adjust this value as needed
gs_b = plt.GridSpec(3, 4, left=left_margin_b, right=right_margin_b, bottom=bottom_margin_b, top=top_margin_b, wspace=0.0, hspace=0.04)
all_axes_b = []
omega_values_b = [10, 5, 1]
dr_values_b = [0.1, 1.0, 10.0, 100.0]

num_rows_b = 3
num_cols_b = 4
font_size = 18
all_axes_b = []

for i in range(num_rows_b):
    for j in range(num_cols_b):
        ax = plt.subplot(gs_b[i, j])
        ax.set_aspect('equal', 'box')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlim([-52, 52])
        ax.set_ylim([-52, 52])
        all_axes_b.append(ax)

for ax in all_axes_b:
    for spine in ax.spines.values():
        spine.set_linewidth(0)



# Annotate the top row
for i, ax in enumerate(all_axes_b[-num_cols_b:]):
    ax.text(0.5, -0.1, str(dr_values_b[i]), transform=ax.transAxes,
            ha='center', va='center', fontsize=font_size)
# Annotate the left column
for i, ax in enumerate(all_axes_b[::num_cols_b]):
    ax.text(-0.05, 0.5, str(omega_values_b[i]), transform=ax.transAxes, ha='right', va='center', fontsize=font_size)



center_subplot_col = all_axes_b[num_cols_b // 2]
center_subplot_row = all_axes_b[num_rows_b // 2]
center_subplot_row.text(-1.20, -0.50, "$\\Omega$", transform=center_subplot_row.transAxes,
                    ha='right', va='center', fontsize=font_size, rotation='vertical')

center_subplot_col.text(0.0, -2.30, "$D_r$", transform=center_subplot_col.transAxes,
                    ha='right', va='center', fontsize=font_size)

center_subplot_row.text(-1.150, 0.95, "(b)", transform=center_subplot_row.transAxes,
                    ha='right', va='center', fontsize=font_size)


cbar_ax_position = [right_margin_b+0.005, bottom_margin_b+0.2, 0.02, top_margin_b-bottom_margin_b-0.4]
ax = fig.add_axes([0, 0, 0.01, 0.01], visible=False)
cax = ax.imshow(data, cmap=plt.cm.viridis)
cbar_ax = fig.add_axes(cbar_ax_position)
cbar = fig.colorbar(cax, cax=cbar_ax, orientation='vertical', ticks=[0, 0.5, 1])
cbar.ax.set_yticklabels([r'$0$', r'$0.5$', r'$1$'], fontsize=font_size)
cbar.ax.tick_params(direction='in', length=6)  # Set tick direction inward and length to 6
cbar.set_label(r'$|\mathbf{v}|/v_0$', fontsize=font_size, labelpad=0)


for omega in omega_values_b:
    for dr in dr_values_b:
        subdir = f'omega_{omega}/dr_{dr}/dynamics_data/positions_1999.csv'
        df = pd.read_csv(subdir)
        xval = df['xval'].values
        yval = df['yval'].values
        vxval = df['vxval'].values
        vyval = df['vyval'].values
        angles = np.sqrt(vxval**2 + vyval**2)/v0

        col = dr_values_b.index(dr)
        row = omega_values_b.index(omega)
        ax = all_axes_b[row * num_cols_b + col]

        for i in range(len(xval)):
            ax.add_patch(Circle((xval[i], yval[i]), 1.0, color=vector_to_rgb_b(angles[i]), fill=True, lw=0.0))

plt.show()

fig.savefig('SM_fig1.png', dpi=300)
