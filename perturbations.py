import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os

# Planck-like colormap
from matplotlib.colors import ListedColormap
colombi1_cmap = ListedColormap(np.loadtxt("Planck_Parchment_RGB.txt")/255.)
colombi1_cmap.set_bad("gray") # color of missing pixels
colombi1_cmap.set_under("white") 
cmap = colombi1_cmap

planck_colors = [cmap(i) for i in np.linspace(0,1,100)]
col = {
    'dark blue'     : planck_colors[0],
    'blue'          : planck_colors[10],
    'light blue'    : planck_colors[25],
    'beige'         : planck_colors[50],
    'yellow'        : planck_colors[60],
    'orange'        : planck_colors[68],
    'red'           : planck_colors[88],
    'black'         : 'black'
}

linestyles = ['-', 'dashed', ':']

k_array = [0.1, 0.01, 0.001]    # k-values considered
data = {}                       # Dictionary containing sim data


for k in k_array:
    # Filename for specific k-value
    filename = f'perturbation_data/perturbations_k{k}.txt' 

    # Unpack data
    x, Theta0, Theta1, Theta2, ThetaP0, ThetaP1, ThetaP2, Nu0, Nu1, Nu2, Phi, Psi, Pi, delta_b, delta_cdm, v_b, v_cdm = np.loadtxt(filename, unpack=True, skiprows=1)

    with open(filename, 'r') as infile:
        horizon_enter = float(infile.readline())

    # Store data in dictionary
    data[k] = {
        'x'             : x,
        'horizon_enter' : horizon_enter,
        'Theta0'        : Theta0,
        'Theta1'        : Theta1,
        'Theta2'        : Theta2,
        'ThetaP0'       : ThetaP0,
        'ThetaP1'       : ThetaP1,
        'ThetaP2'       : ThetaP2, 
        'Nu0'           : Nu0,
        'Nu1'           : Nu1,
        'Nu2'           : Nu2,
        'Phi'           : Phi,
        'Psi'           : Psi,
        'Pi'            : Pi,
        'delta_b'       : delta_b,
        'delta_cdm'     : delta_cdm,
        'v_b'           : v_b,
        'v_cdm'         : v_cdm
    }

###########################
# DENSITY PERTURBATIONS   #
###########################
fig, axs = plt.subplots(2, 2, figsize=(10,7), sharex=True)

cnames = ['dark blue', 'orange', 'red']
for k, cname in zip(k_array, cnames):
    x               = data[k]['x']
    horizon_enter   = data[k]['horizon_enter']

    # Perturbations
    delta_gamma = 4.0*data[k]['Theta0']
    delta_Nu    = 4.0*data[k]['Nu0']
    v_gamma     = -3.0*data[k]['Theta1']
    v_Nu        = -3.0*data[k]['Nu1']
    delta_b     = data[k]['delta_b']
    delta_cdm   = data[k]['delta_cdm']
    v_b         = data[k]['v_b']
    v_cdm       = data[k]['v_cdm']

    c = col[cname]

    # Density perturbations: radiation
    axs[0,0].plot(x, delta_gamma, color=c, alpha=0.5)
    axs[0,0].plot(x, delta_Nu, color=c, linestyle='dashed')

    # Velocity perturbations: radiation
    axs[1,0].plot(x, v_gamma, color=c, alpha=0.5)
    axs[1,0].plot(x, v_Nu, color=c, linestyle='dashed')

    # Density perturbations: matter
    axs[0,1].plot(x, abs(delta_cdm), color=c, alpha=0.5)
    axs[0,1].plot(x, abs(delta_b), color=c, linestyle='dashed')

    # Velocity perturbations: matter
    axs[1,1].plot(x, abs(v_cdm), color=c, alpha=0.5)
    axs[1,1].plot(x, abs(v_b), color=c, linestyle='dashed')

    if horizon_enter != float('nan'):
        cross_idx = np.argmin(abs(x-horizon_enter))
        for ax in axs.flat:
            ax.axvline(x[cross_idx], color=c, alpha=0.3, linestyle='dotted')

x_text, y_text = 0.9, 0.07
for i, ax in enumerate(axs.flat):
    ax.text(x_text, y_text, f'({chr(i+97)})', fontweight='bold', transform=ax.transAxes)


for ax in [axs[0,0], axs[1,0]]:
    ax.yaxis.tick_left()

axs[0,0].set_xlim(x[0], x[-1])
axs[1,0].set_xlabel(r'$x=\ln a$')
axs[0,0].set_ylabel(r'$\delta_\gamma, \;\delta_\nu$')
axs[1,0].set_ylabel(r'$v_\gamma, \;v_\nu$')

for ax in [axs[0,1], axs[1,1]]:
    ax.set_yscale('log')
    ax.yaxis.set_label_position('right')
    ax.yaxis.tick_right()
    
axs[0,1].set_xlim(x[0], x[-1])
axs[0,1].set_ylabel(r'$|\delta_b|,\;|\delta_\mathrm{CDM}|$')
axs[1,1].set_ylabel(r'$|v_b|,\;|v_\mathrm{CDM}|$')
axs[1,1].set_xlabel(r'$x=\ln a$')

line_solid = Line2D([], [], color='black', linestyle='-')
line_solid_faded = Line2D([], [], color='black', linestyle='-', alpha=0.5)
line_dashed = Line2D([], [], color='black', linestyle='dashed')
fs = 13

fig.legend(
    loc='upper center', 
    labels=[r'$k=$'+f' {k}/Mpc' for k in k_array], 
    handles=[Line2D([], [], color=col[cname], linestyle='-') for cname in cnames],
    ncols=3
)
axs[0,0].legend(
    loc='upper left', 
    labels=[r'$\delta_\gamma=4\Theta_0$', r'$\delta_\nu=4\mathcal{N}_0$'], 
    handles=[line_solid_faded, line_dashed],
    fontsize=fs
)
axs[1,0].legend(
    loc='upper left', 
    labels=[r'$v_\gamma=-3\Theta_1$', r'$v_\nu = -3\mathcal{N}_1$'], 
    handles=[line_solid_faded, line_dashed],
    fontsize=fs
)
axs[0,1].legend(
    loc='upper left', 
    labels=[r'$\delta_\mathrm{CDM}$', r'$\delta_b$'], 
    handles=[line_solid_faded, line_dashed],
    fontsize=fs
)
axs[1,1].legend(
    loc='upper left', 
    labels=[r'$v_\mathrm{CDM}$', r'$v_b$'], 
    handles=[line_solid_faded, line_dashed],
    fontsize=fs
)

fig.subplots_adjust(hspace=0, wspace=0, top=0.92, bottom=0.15)
fig.savefig('figures/matter_pert.pdf')




###########################
# PHI, PSI                #
###########################
fig, axs = plt.subplots(2, 1, figsize=(5,5), sharex=True, height_ratios=[0.6, 0.4])
cnames = ['dark blue', 'orange', 'red']

for k, cname in zip(k_array, cnames):
    x               = data[k]['x']
    horizon_enter   = data[k]['horizon_enter']
    Phi             = data[k]['Phi']
    Psi             = data[k]['Psi']

    c = col[cname]

    axs[0].plot(x, Phi, color=c, alpha=0.5)
    axs[0].plot(x, Psi, color=c, linestyle='dashed')
    axs[1].plot(x, Phi+Psi, color=c, label=r'$k=$'+f' {k}'+r'/Mpc')

    if horizon_enter != float('nan'):
        cross_idx = np.argmin(abs(x-horizon_enter))
        axs[0].axvline(x[cross_idx], color=c, alpha=0.3, linestyle='dotted')
        axs[1].axvline(x[cross_idx], color=c, alpha=0.3, linestyle='dotted')

axs[0].legend(
    loc='center left',
    labels=[r'$\Phi$', r'$\Psi$'], 
    handles=[line_solid_faded, line_dashed]
)
axs[0].set_xlim(x[0], x[-1])
axs[1].set_xlabel(r'$x=\ln a$')
axs[1].legend(loc='lower left', fontsize=11)
axs[0].set_ylabel(r'$\Phi,\Psi$')
axs[1].set_ylabel(r'$\Phi+\Psi$')

fig.subplots_adjust(hspace=0, bottom=0.15, left=0.15)
fig.savefig('figures/phi_psi_pert.pdf')





###########################
# PHOTON POLARISATION     #
###########################
fig, axs = plt.subplots(3, 1, figsize=(6,6), sharex=True, )
cnames = ['dark blue', 'orange', 'red']

for k, cname in zip(k_array, cnames):
    x           = data[k]['x']
    ThetaP0     = data[k]['ThetaP0']
    ThetaP1     = data[k]['ThetaP1']
    ThetaP2     = data[k]['ThetaP2']

    c = col[cname]

    axs[0].plot(x, ThetaP0, color=c, alpha=0.85, label=r'$k=$'+f' {k}'+r'/Mpc')
    axs[1].plot(x, ThetaP1, color=c, alpha=0.85)
    axs[2].plot(x, ThetaP2, color=c, alpha=0.85)

for ax in axs:
    ax.set_xlim(-9.5, 0)
    ax.set_xlabel(r'$x=\ln a$')
    ax.axvline(-6.97013, color='black', linestyle='dotted', alpha=1)

axs[0].text(0.80, 0.92, r'$\Theta^P_0(x,k)$', transform=axs[0].transAxes, 
            verticalalignment='top', color='black')
axs[1].text(0.80, 0.92, r'$\Theta^P_1(x,k)$', transform=axs[1].transAxes, 
            verticalalignment='top', color='black')
axs[2].text(0.80, 0.92, r'$\Theta^P_2(x,k)$', transform=axs[2].transAxes, 
            verticalalignment='top', color='black')

fig.legend(loc='upper center', ncols=3, fontsize=14)

fig.subplots_adjust(hspace=0, top=0.90, bottom=0.12)

fig.savefig('figures/polarization_pert.pdf')


###########################
# QUADRUPOLE              #
###########################
fig, ax = plt.subplots(1,1)

cnames = ['dark blue', 'orange', 'red']

for k, cname in zip(k_array, cnames):
    x       = data[k]['x']
    Theta2  = data[k]['Theta2']
    Nu2     = data[k]['Nu2']

    ax.plot(x, Theta2*np.exp(-2*x), color=col[cname], alpha=0.5)
    ax.plot(x, Nu2*np.exp(-2*x), color=col[cname], linestyle='dashed')

    if horizon_enter != float('nan'):
        cross_idx = np.argmin(abs(x-horizon_enter))
        axs[0].axvline(x[cross_idx], color=c, alpha=0.3, linestyle='dotted')
        axs[1].axvline(x[cross_idx], color=c, alpha=0.3, linestyle='dotted')

ax.set_xlabel(r'$x=\ln a$')
ax.set_ylabel(r'$\Theta_2/a^2,\;\mathcal{N}_2/a^2$')
ax.set_yscale('symlog', linthresh=1e2, linscale=1)


ax.legend(
    loc='upper left',
    labels=[r'$k=$'+f' {k}/Mpc' for k in k_array] + [r'$\Theta_2$', r'$\mathcal{N}_2$'], 
    handles=[Line2D([], [], color=col[cname], linestyle='-') for cname in cnames] + [line_solid_faded, line_dashed]
)

fig.savefig('figures/quadrupoles.pdf')

##########################
# PHOTON-BARYON COUPLING #
##########################
fig, ax = plt.subplots(1,1)
k = 0.1

x               = data[k]['x']
delta_b         = data[k]['delta_b']
delta_gamma     = 4.0*data[k]['Theta0']

ax.axvline(data[k]['horizon_enter'], color=col['dark blue'], alpha=0.3, linestyle='dotted')
ax.axvline(-6.97013, color='black', linestyle='dotted')
ax.plot(x, delta_b, color=col['dark blue'], linestyle='dashed')
ax.plot(x, delta_gamma, color=col['dark blue'], alpha=0.5)

ax.text(0.25, 1.03, r'HE', transform=ax.transAxes)
ax.text(0.71, 1.03, r'Rec.', transform=ax.transAxes)

ax.set_yscale('symlog', linthresh=5, linscale=1)
ax.set_xlabel(r'$x=\ln a$')
ax.set_ylabel(r'$\delta_b,\; \delta_\gamma$')
ax.set_ylim(-4, 1e2)
ax.set_xlim(-13, -5)

ax.legend(
    loc='upper left',
    labels=[r'$k=0.1/\mathrm{Mpc}$', r'$\delta_\gamma = 4\Theta_0$', r'$\delta_b$'], 
    handles=[Line2D([], [], color=col['dark blue'], linestyle='-'), line_solid_faded, line_dashed]
)

fig.savefig('figures/photon_baryon_coupling.pdf')

os.system(
"""
cd paper
pdflatex -interaction=batchmode ast5220_paper
cd ..
"""
)
