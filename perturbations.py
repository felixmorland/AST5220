import matplotlib.pyplot as plt
import numpy as np

# Plot colors and linestyles
custom_col = ["#476CFF", "#DE4E61", "#32A639", "#ECC723", "#EC7023"]
linestyles = ['-', 'dashed', ':']

k_array = [0.001, 0.01, 0.1]    # k-values considered
data = {}                       # Dictionary containing sim data

for k in k_array:
    # Filename for specific k-value
    filename = f'perturbation_data/perturbations_k{k}.txt' 

    # Unpack data
    x, Theta0, Theta1, Theta2, ThetaP0, ThetaP1, ThetaP2, Nu0, Nu1, Nu2, Phi, Psi, Pi, delta_b, delta_cdm, v_b, v_cdm = np.loadtxt(filename, unpack=True)

    # Store data in dictionary
    data[k] = {
        'x': x,
        'Theta0'    : Theta0,
        'Theta1'    : Theta1,
        'Theta2'    : Theta2,
        'ThetaP0'   : ThetaP0,
        'ThetaP1'   : ThetaP1,
        'ThetaP2'   : ThetaP2, 
        'Nu0'       : Nu0,
        'Nu1'       : Nu1,
        'Nu2'       : Nu2,
        'Phi'       : Phi,
        'Psi'       : Psi,
        'Pi'        : Pi,
        'delta_b'   : delta_b,
        'delta_cdm' : delta_cdm,
        'v_b'       : v_b,
        'v_cdm'     : v_cdm
    }

###########################
# THETA0                  #
###########################
fig, axs = plt.subplots(2, 1, figsize=(6,8))

for k, c in zip(k_array, custom_col):
    x       = data[k]['x']
    delta_gamma = 4.0*data[k]['Theta0']
    delta_Nu    = 4.0*data[k]['Nu0']
    v_gamma     = -3.0*data[k]['Theta1']
    v_Nu        = -3.0*data[k]['Nu1']

    # Density perturbations
    axs[0].plot(x, delta_gamma, color=c, alpha=0.5, label=r'$k=$'+f' {k}'+r'/Mpc')
    axs[0].plot(x, delta_Nu, color=c, linestyle='dashed')

    # Velocity perturbations
    axs[1].plot(x, v_gamma, color=c, alpha=0.5, label=r'$k=$'+f' {k}'+r'/Mpc')
    axs[1].plot(x, v_Nu, color=c, linestyle='dashed')

axs[0].legend(loc='lower left')
axs[0].set_xlim(x[0], x[-1])
axs[0].set_xlabel(r'$x=\ln a$')
axs[0].set_ylabel(r'$\Theta_0(x,k)$')
fig.savefig('figures/radiation_pert.pdf')

###########################
# DENSITY PERT            #
###########################
fig, axs = plt.subplots(2, 1, figsize=(6,8))

for k, c in zip(k_array, custom_col):
    x           = data[k]['x']
    delta_b     = data[k]['delta_b']
    delta_cdm   = data[k]['delta_cdm']
    v_b         = data[k]['v_b']
    v_cdm       = data[k]['v_cdm']

    axs[0].plot(x, abs(delta_cdm), color=c, alpha=0.5, label=r'$k=$'+f' {k}'+r'/Mpc')
    axs[0].plot(x, abs(delta_b), color=c, linestyle='dashed')

    axs[1].plot(x, abs(v_cdm), color=c, alpha=0.5, label=r'$k=$'+f' {k}'+r'/Mpc')
    axs[1].plot(x, abs(v_b), color=c, linestyle='dashed')
    
axs[0].legend(loc='upper left')
axs[0].set_yscale('log')
axs[1].set_yscale('log')
axs[0].set_xlim(x[0], x[-1])
axs[0].set_xlabel(r'$x=\ln a$')
axs[0].set_ylabel(r'$|\delta_b|,\;|\delta_\mathrm{CDM}|$')
fig.savefig('figures/matter_pert.pdf')

###########################
# PHI, PSI                #
###########################
fig, axs = plt.subplots(2,1)

for k, c in zip(k_array, custom_col):
    x           = data[k]['x']
    Phi         = data[k]['Phi']
    Psi         = data[k]['Psi']
    axs[0].plot(x, Phi, color=c, alpha=0.5, label=r'$k=$'+f' {k}'+r'/Mpc')
    axs[0].plot(x, Psi, color=c, linestyle='dashed')
    axs[1].plot(x, Phi+Psi, color=c)

axs[0].legend(loc='upper left')
axs[0].set_xlim(x[0], x[-1])
axs[0].set_xlabel(r'$x=\ln a$')
axs[0].set_ylabel(r'$|\delta_b|,\;|\delta_\mathrm{CDM}|$')
fig.savefig('figures/phi_psi_pert.pdf')

###########################
# PHOTON POLARIZATION     #
###########################
fig, axs = plt.subplots(3, 1, figsize=(6,10), sharex=True, )

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