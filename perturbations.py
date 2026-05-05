import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

# Plot colors and linestyles
custom_col = ["#476CFF", "#DE4E61", "#32A639"]
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

for k, c in zip(k_array, custom_col):
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

    # Density perturbations: radiation
    axs[0,0].plot(x, abs(delta_gamma), color=c, alpha=0.5)
    axs[0,0].plot(x, abs(delta_Nu), color=c, linestyle='dashed')

    # Velocity perturbations: radiation
    axs[1,0].plot(x, abs(v_gamma), color=c, alpha=0.5)
    axs[1,0].plot(x, abs(v_Nu), color=c, linestyle='dashed')

    # Density perturbations: matter
    axs[0,1].plot(x, abs(delta_cdm), color=c, alpha=0.5)
    axs[0,1].plot(x, abs(delta_b), color=c, linestyle='dashed')

    # Velocity perturbations: matter
    axs[1,1].plot(x, abs(v_cdm), color=c, alpha=0.5)
    axs[1,1].plot(x, abs(v_b), color=c, linestyle='dashed')

    if horizon_enter != float('nan'):
        cross_idx = np.argmin(abs(x-horizon_enter))
        axs[0,1].axvline(x[cross_idx], color=c, alpha=0.3, linestyle='dotted')
        axs[1,1].axvline(x[cross_idx], color=c, alpha=0.3, linestyle='dotted')

for ax in [axs[0,0], axs[1,0]]:
    ax.set_yscale('log')
    ax.yaxis.tick_left()

axs[0,0].set_xlim(x[0], x[-1])
axs[1,0].set_xlabel(r'$x=\ln a$')
axs[0,0].set_ylabel(r'$|\delta_\gamma|, \;|\delta_\nu|$')
axs[1,0].set_ylabel(r'$|v_\gamma|, \;|v_\nu|$')

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
    loc='lower center', 
    labels=[r'$k=$'+f' {k}/Mpc' for k in k_array], 
    handles=[Line2D([], [], color=col, linestyle='-') for col in custom_col[0:3]],
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

fig.subplots_adjust(hspace=0, wspace=0, bottom=0.15)
fig.savefig('figures/matter_pert.pdf')




###########################
# PHI, PSI                #
###########################
fig, axs = plt.subplots(2, 1, figsize=(5.5,5.5), sharex=True, height_ratios=[0.6, 0.4])

for k, c in zip(k_array, custom_col):
    x               = data[k]['x']
    horizon_enter   = data[k]['horizon_enter']
    Phi             = data[k]['Phi']
    Psi             = data[k]['Psi']

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
    handles=[line_solid_faded, line_dashed],
    fontsize=13
)
axs[0].set_xlim(x[0], x[-1])
axs[1].set_xlabel(r'$x=\ln a$')
axs[1].legend(loc='lower left', fontsize=11)
axs[0].set_ylabel(r'$\Phi,\Psi$')
axs[1].set_ylabel(r'$\Phi+\Psi$')

fig.subplots_adjust(hspace=0, bottom=0.15, left=0.15)
fig.savefig('figures/phi_psi_pert.pdf')




###########################
# PHOTON POLARIZATION     #
###########################
fig, axs = plt.subplots(3, 1, figsize=(6,9), sharex=True, )

for k, c in zip(k_array[::-1], custom_col[2::-1]):
    x           = data[k]['x']
    ThetaP0     = data[k]['ThetaP0']
    ThetaP1     = data[k]['ThetaP1']
    ThetaP2     = data[k]['ThetaP2']

    axs[0].plot(x, ThetaP0, color=c, alpha=0.85, label=r'$k=$'+f' {k}'+r'/Mpc')
    axs[1].plot(x, ThetaP1, color=c, alpha=0.85)
    axs[2].plot(x, ThetaP2, color=c, alpha=0.85)

for ax in axs:
    ax.set_xlim(-9.5, 0)
    ax.set_xlabel(r'$x=\ln a$')

axs[0].text(0.80, 0.92, r'$\Theta^P_0(x,k)$', transform=axs[0].transAxes, 
            verticalalignment='top', color='black')
axs[1].text(0.80, 0.92, r'$\Theta^P_1(x,k)$', transform=axs[1].transAxes, 
            verticalalignment='top', color='black')
axs[2].text(0.80, 0.92, r'$\Theta^P_2(x,k)$', transform=axs[2].transAxes, 
            verticalalignment='top', color='black')

fig.legend(loc='lower center', ncols=3, frameon=True, shadow=True)

fig.subplots_adjust(hspace=0, bottom=0.12)

# ax.legend(loc='upper left')
# ax.set_yscale('log')
# ax.set_xlim(-10.0, x[-1])
# # ax.set_ylim(1e-15, 1e-2)
# ax.set_xlabel(r'$x=\ln a$')
# ax.set_ylabel(r'$|\delta_b|,\;|\delta_\mathrm{CDM}|$')
fig.savefig('figures/polarization_pert.pdf')
