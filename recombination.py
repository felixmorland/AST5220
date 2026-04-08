import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import norm

# COLORS
cblue = "#6D79F8"
cred = "#DE4E61"
cgreen = "#32A639"

# Recombination data: H + He, reionisation
lna, Xe, ne, tau, dtaudx, ddtauddx, g_tilde, dgdx, ddgddx, sound_hor, Xe_saha = np.loadtxt('recombination_data/recombination_He_reion.txt', unpack=True)
decoupling = -6.97282

# Recombination data: H + He, NO reionisation
lna_He_noreion, Xe_He_noreion, ne_He_noreion, tau_He_noreion, dtaudx_He_noreion, ddtauddx_He_noreion, g_tilde_He_noreion, dgdx_He_noreion, ddgddx_He_noreion, sound_hor_He_noreion, Xe_saha_He_noreion = np.loadtxt('recombination_data/recombination_He_noreion.txt', unpack=True)

# Recombination data: H, NO reionisation
lna_H_noreion, Xe_H_noreion, ne_H_noreion, tau_H_noreion, dtaudx_H_noreion, ddtauddx_H_noreion, g_tilde_H_noreion, dgdx_H_noreion, ddgddx_H_noreion, sound_hor_H_noreion, Xe_saha_H_noreion = np.loadtxt('recombination_data/recombination_H_noreion.txt', unpack=True)

# Recombination data: H, reionisation
lna_H_reion, Xe_H_reion, ne_H_reion, tau_H_reion, dtaudx_H_reion, ddtauddx_H_reion, g_tilde_H_reion, dgdx_H_reion, ddgddx_H_reion, sound_hor_H_reion, Xe_saha_H_reion = np.loadtxt('recombination_data/recombination_H_reion.txt', unpack=True)


#===============
# OPTICAL DEPTH
#===============
fig, axs = plt.subplots(2, 2, figsize=(10,7))

axs[0,0].plot(lna_H_noreion, tau_H_noreion, color='black', linestyle='-')
axs[0,0].plot(lna_H_noreion, -dtaudx_H_noreion, color=cblue, linestyle='-')
axs[0,0].plot(lna_H_noreion, ddtauddx_H_noreion, color=cred, linestyle='-')
axs[0,0].text(0.6, 0.95, r'H no reionisation', transform=axs[0,0].transAxes, 
            verticalalignment='top', fontsize=13)
axs[0,0].text(0.81, 0.87, r'$Y_p=0$', transform=axs[0,0].transAxes, 
            verticalalignment='top', fontsize=13)
axs[0,0].text(0.06, 0.15, r'$\mathbf{(a)}$', transform=axs[0,0].transAxes, 
            verticalalignment='top', fontsize=16)

axs[0,1].plot(lna_H_reion, tau_H_reion, color='black', linestyle='-')
axs[0,1].plot(lna_H_reion, -dtaudx_H_reion, color=cblue, linestyle='-')
axs[0,1].plot(lna_H_reion, ddtauddx_H_reion, color=cred, linestyle='-')
axs[0,1].text(0.56, 0.95, r'H with reionisation', transform=axs[0,1].transAxes, 
            verticalalignment='top', fontsize=13)
axs[0,1].text(0.81, 0.87, r'$Y_p=0$', transform=axs[0,1].transAxes, 
            verticalalignment='top', fontsize=13)
axs[0,1].text(0.06, 0.15, r'$\mathbf{(b)}$', transform=axs[0,1].transAxes, 
            verticalalignment='top', fontsize=16)

axs[1,0].plot(lna_He_noreion, tau_He_noreion, color='black', linestyle='-')
axs[1,0].plot(lna_He_noreion, -dtaudx_He_noreion, color=cblue, linestyle='-')
axs[1,0].plot(lna_He_noreion, ddtauddx_He_noreion, color=cred, linestyle='-')
axs[1,0].text(0.5, 0.95, r'H + He no reionisation', transform=axs[1,0].transAxes, 
            verticalalignment='top', fontsize=13)
axs[1,0].text(0.75, 0.87, r'$Y_p=0.245$', transform=axs[1,0].transAxes, 
            verticalalignment='top', fontsize=13)
axs[1,0].text(0.06, 0.15, r'$\mathbf{(c)}$', transform=axs[1,0].transAxes, 
            verticalalignment='top', fontsize=16)

axs[1,1].plot(lna, tau, label=r'$\tau(x)$', color='black', linestyle='-')
axs[1,1].plot(lna, -dtaudx, label=r"$-\tau'(x)$", color=cblue, linestyle='-')
axs[1,1].plot(lna, ddtauddx, label=r"$\tau''(x)$", color=cred, linestyle='-')
axs[1,1].text(0.45, 0.95, r'H + He with reionisation', transform=axs[1,1].transAxes, 
            verticalalignment='top', fontsize=13)
axs[1,1].text(0.75, 0.87, r'$Y_p=0.245$', transform=axs[1,1].transAxes, 
            verticalalignment='top', fontsize=13)
axs[1,1].text(0.06, 0.15, r'$\mathbf{(d)}$', transform=axs[1,1].transAxes, 
            verticalalignment='top', fontsize=16)

for ax in axs.flat:
    ax.set_yscale('log')
    ax.set_xlim(-10.5,0)
    ax.set_ylim(0.9e-8, 1e7)
    ax.tick_params(axis='x', top=True)
    ax.set_xlabel(r'$x=\ln a$')
    ax.label_outer()

fig.legend(loc='upper center', frameon=True, shadow=True, ncols=3)
fig.subplots_adjust(hspace=0.0, wspace=0.0)
plt.savefig('figures/optical_depth_derivatives.pdf')


#========================
# FREE ELECTRON FRACTION
#========================
fig, axs = plt.subplots(2, 1, figsize=(6,6), gridspec_kw={'height_ratios': [2, 1], 'hspace': 0.0})
axs[0].plot(lna, Xe, color='black', label=r'H + He with reionisation')
axs[0].plot(lna, Xe_H_noreion, color='black', linestyle='dashed', label=r'H no reionisation')
axs[0].plot(lna, Xe_saha, color='black', linestyle=':', label=r'Saha (H + He)')
axs[0].axvspan(-7.5, -6.5, alpha=0.2, color='gray')
# ax.set_yscale('log')

axs[0].set_xlim(-12, 0)
axs[0].set_ylim(-0.05, 1.5)
axs[0].set_xlabel(r'$x=\ln a$')
axs[0].set_ylabel(r'$X_e(x)$')
axs[0].tick_params(axis='x', bottom=False, labeltop=True, labelbottom=False)
axs[0].xaxis.set_label_position('top')
axs[0].legend(loc='upper right', fontsize=13)


axs[1].plot(lna, Xe, color='black')
axs[1].plot(lna, Xe_H_noreion, color='black', linestyle='dashed')
axs[1].plot(lna, Xe_saha_H_noreion, color='black', linestyle=':')
axs[1].tick_params(axis='x', top=False)
axs[1].axvline(-7.14035, color=cblue, label=r'Saha $x_\mathrm{rec}=-7.14$')
axs[1].axvline(-6.97013, color=cred, label=r'Peebles $x_\mathrm{rec} = -6.97$')
axs[1].legend(loc='upper right', fontsize=13)
axs[1].set_xlabel(r'$x=\ln a$')
axs[1].set_ylabel(r'$X_e(x)$')
axs[1].set_xlim(-7.5, -6.5)
axs[1].set_ylim(1e-4, 5)
axs[1].set_yscale('log')

plt.savefig('figures/free_electron_frac.pdf')


#========================
# VISIBILITY FUNCTION
#========================

# three visibility-related quantities in a row
fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(6, 5),
                                gridspec_kw={'width_ratios': [2, 1], 'wspace': 0.08})

plot_data = [(g_tilde,    r"$\tilde{g}(x)$",      '-'),
             (dgdx/10,    r"$\tilde{g}'/10$",      'dashed'),
             (ddgddx/300, r"$\tilde{g}''/300$",    ':')]

for data, name, linestyle in plot_data:
    ax1.plot(lna, data, linestyle=linestyle, color='black', label=name)
    ax2.plot(lna, data, linestyle=linestyle, color='black')

ax1.set_xlim(-7.4, -6.3)
ax2.set_xlim(-2.75, -0.5)

for ax in (ax1, ax2):
    ax.tick_params(direction='in')

ax1.axvline(decoupling, color=cblue, label=r'$x_\mathrm{dec}=$'+f' {decoupling:.3f}')

# Hide the spines between the two axes
ax1.spines['right'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax1.yaxis.tick_left()
ax2.yaxis.tick_right()
ax2.tick_params(left=False)

# Draw diagonal break marks on both sides of the cut
d = 0.015  # size of diagonal lines
kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False, linewidth=1)
ax1.plot([1-d, 1+d], [-d, +d], **kwargs)
ax1.plot([1-d, 1+d], [1-d, 1+d], **kwargs)

kwargs.update(transform=ax2.transAxes)
ax2.plot([-2*d, +2*d], [-d, +d], **kwargs)
ax2.plot([-2*d, +2*d], [1-d, 1+d], **kwargs)

# Single shared x-label
fig.text(0.5, 0.02, r'$x = \ln a$', ha='center')
handles, labels = ax1.get_legend_handles_labels()
ax2.legend(handles=handles, labels=labels, loc='upper right')

fig.savefig('figures/visibility_functions.pdf')