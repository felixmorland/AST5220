import matplotlib.pyplot as plt
import numpy as np
import healpy as hp

# Planck-like colormap
from matplotlib.colors import ListedColormap
colombi1_cmap = ListedColormap(np.loadtxt("../data/Planck_Parchment_RGB.txt")/255.)
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

def CMB_map():
    #====================================================
    # Load data
    #====================================================
    data  = np.loadtxt('../results/powerspectrum/cells.txt')
    ells  = data[:, 0].astype(int)
    cl_TT = data[:, 1]   # ell(ell+1)/2pi * C_ell [muK^2]
    cl_EE = data[:, 2]   # (if polarization)
    cl_TE = data[:, 3]   # (if polarization)

    #====================================================
    # Parameters
    #====================================================
    TCMB      = 2.7255e6   # muK
    lmax      = int(ells[-1])
    nside     = 2048        # lmax ~ 2*nside
    fwhm_rad  = np.radians(0.0)
    seed      = 5

    #====================================================
    # Remove ell(ell+1)/2pi normalization
    #====================================================
    norm  = ells * (ells + 1) / (2 * np.pi)
    cl_TT /= norm
    cl_EE /= norm
    cl_TE /= norm

    #====================================================
    # Build full Cl arrays from ell=0 to ell=lmax
    #====================================================
    cl_TT_full = np.zeros(lmax + 1)
    cl_EE_full = np.zeros(lmax + 1)
    cl_BB_full = np.zeros(lmax + 1)   # assume no B-modes
    cl_TE_full = np.zeros(lmax + 1)

    cl_TT_full[ells] = cl_TT
    cl_EE_full[ells] = cl_EE
    cl_TE_full[ells] = cl_TE

    #====================================================
    # Generate maps
    #====================================================
    np.random.seed(seed)
    map_T, map_Q, map_U = hp.synfast(
        [cl_TT_full, cl_EE_full, cl_BB_full, cl_TE_full],
        nside   = nside,
        lmax    = lmax,
        fwhm    = fwhm_rad,
        new     = True,
    )

    #====================================================
    # Plot
    #====================================================
    plot_kwargs = dict(cmap=cmap, unit=r"$\mu K$", xsize=3000)
    hp.mollview(map_T, title="Simulated CMB map (TT)", min=-500, max=500, **plot_kwargs)
    # hp.mollview(map_Q, title="CMB E-polarisation (EE)", **plot_kwargs)
    # hp.mollview(map_U, title="CMB Cross (TE)", **plot_kwargs)

    plt.savefig('../figures/simulated_CMB_map.pdf')

# Power spectra
filenames = {
    'CMB power-spectrum'    : '../results/powerspectrum/cells.txt',
    'Master level PS'       : '../results/powerspectrum/cells_master.txt',
    'POL'                   : '../results/powerspectrum/Polarization_ps.txt',
    'Doppler'               : '../results/powerspectrum/Doppler_ps.txt',
    'ISW'                   : '../results/powerspectrum/ISW_ps.txt',
    'SW'                    : '../results/powerspectrum/SW_ps.txt'
}
spectra = {}
for file in filenames.keys():
    try:
        ells, cell_TT, cell_EE, cell_TE = np.loadtxt(filenames[file], unpack=True)
        spectra[file] = {
            'ells'  : ells,
            'TT'    : cell_TT,
            'EE'    : cell_EE,
            'TE'    : cell_TE
        }
    except:
        ells, cell_TT  = np.loadtxt(filenames[file], unpack=True)
        spectra[file] = {
            'ells'  : ells,
            'TT'    : cell_TT
        }
k_array, matter_ps = np.loadtxt('../results/powerspectrum/matter_ps.txt', unpack=True, skiprows=1) 
with open('../results/powerspectrum/matter_ps.txt', 'r') as infile:
    k_eq = float(infile.readline())

# Observational data
low_ell_TT_data     = np.loadtxt('../data/low_ell_TT.txt', skiprows=1)
high_ell_TT_data    = np.loadtxt('../data/high_ell_TT.txt', skiprows=1)
high_ell_EE_data    = np.loadtxt('../data/high_ell_EE.txt', skiprows=1)
high_ell_TE_data    = np.loadtxt('../data/high_ell_TE.txt', skiprows=1)

SDSS_gal            = np.loadtxt('../data/SDSS_DR7_LRG.txt', skiprows=1)
WMAP_ACT            = np.loadtxt('../data/WMAP_ACT.txt', skiprows=1)
Ly_alpha            = np.loadtxt('../data/lyalpha.txt', skiprows=1)



# Transfer functions
tf_dataT = np.genfromtxt('../results/powerspectrum/transfer_funcT.txt')
transfer_funcsT = {'k':tf_dataT[1:,0]}
for ell in range(1, tf_dataT.shape[1]):
    transfer_funcsT[int(tf_dataT[0,ell])] = tf_dataT[1:,ell]

tf_dataE = np.genfromtxt('../results/powerspectrum/transfer_funcE.txt')
transfer_funcsE = {'k':tf_dataE[1:,0]}
for ell in range(1, tf_dataE.shape[1]):
    transfer_funcsE[int(tf_dataE[0,ell])] = tf_dataE[1:,ell]

# Source functions
source_funcs = {}


###############################
# MATTER POWER SPECTRUM       #
###############################
fig, ax = plt.subplots(1, 1, figsize=(6,5))

ax.axvline(k_eq, linestyle='dashed', color='black', alpha=0.5)
fig.text(0.31, 1.03, r'$k_\mathrm{eq}$', ha='center', transform=ax.transAxes)
ax.plot(k_array, matter_ps, color=col['red'], label=r'Theory')

ax.errorbar(WMAP_ACT[:,0], WMAP_ACT[:,1], 
            yerr=WMAP_ACT[:,2]-WMAP_ACT[:,1],
            fmt='o', markerfacecolor='none', color=col['dark blue'], linewidth=1,
            markeredgecolor=col['dark blue'], capsize=3, label=r'CMB (WMAP + ACT)')

ax.errorbar(SDSS_gal[:,0], SDSS_gal[:,1], 
            yerr=SDSS_gal[:,2], alpha=0.8,
            fmt='x', markerfacecolor='none', color='black', linewidth=1,
            markeredgecolor='black', capsize=2, label=r'SDSS Galaxies (DR7 LRG)')

ax.errorbar(Ly_alpha[:,0], Ly_alpha[:,1], 
            yerr=Ly_alpha[:,2]-Ly_alpha[:,1],
            fmt='^', markerfacecolor='none', color=col['orange'], linewidth=1,
            markeredgecolor=col['orange'], capsize=3, label=r'Ly$\alpha$ forest (BOSS)')

axins = ax.inset_axes([0.15, 0.32, 0.4, 0.4], transform=ax.transAxes) 
axins.plot(k_array, matter_ps, color=col['red'])
axins.errorbar(SDSS_gal[:,0], SDSS_gal[:,1], 
            yerr=SDSS_gal[:,2], alpha=0.6,
            fmt='x', markerfacecolor='none', color='black', linewidth=1,
            markeredgecolor='black', capsize=2, label=r'SDSS Galaxies (DR7 LRG)')
axins.set_yscale('log')
axins.set_xscale('log')

x1, x2, y1, y2 = 1.5e-1, 2.9e-1, 8e2, 3e3
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
axins.set_xticklabels([])
axins.set_yticklabels([])
axins.set_xticklabels([], minor=True)
axins.set_yticklabels([], minor=True)

ax.indicate_inset_zoom(axins, edgecolor="black")

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(1,1e5)
ax.set_xlim(2e-3, k_array[-1])
ax.legend(loc='lower right', fontsize=12)

ax.set_xlabel(r'$k\;[h/\mathrm{Mpc}]$')
ax.set_ylabel(r'$\mathcal{P}(k)\;[(\mathrm{Mpc}/h)^3]$')

fig.savefig('../figures/matter_power_spectrum.pdf')

###############################
# CMB POWER SPECTRUM          #
###############################
fig, axs = plt.subplots(2, 1, figsize=(9,7), height_ratios=[0.6, 0.4], sharex=True)
ells = spectra['POL']['ells']

s_ls    = 144.413           # Mpc
c       = 2.99792458e5  # km/s
H0      = 67            # km/s/Mpc

ell_peak = 2*np.pi*c / (H0*s_ls)

# axs[0].axvline(ell_peak, linestyle='dashed', color='black', alpha=0.5)
axs[0].plot(ells, spectra['CMB power-spectrum']['TT'], color=col['red'], label=r'CMB power spectrum')
axs[0].plot(ells, spectra['Master level PS']['TT'], color=col['orange'], linestyle='dashed',  label=r'Simplified model$^\dagger$')

axs[0].errorbar(low_ell_TT_data[:,0], low_ell_TT_data[:,1], 
            yerr=[low_ell_TT_data[:,2], low_ell_TT_data[:,3]], markersize=5,
            fmt='x', markerfacecolor='none', color=col['dark blue'], linewidth=1,
            markeredgecolor=col['dark blue'], capsize=3, label=r'Low-$\ell$, Planck 2018')

axs[0].errorbar(high_ell_TT_data[:,0], high_ell_TT_data[:,1], 
            yerr=[high_ell_TT_data[:,3], high_ell_TT_data[:,2]], markersize=5,
            fmt='o', markerfacecolor='none', color='black', linewidth=1,
            markeredgecolor='black', capsize=3, label=r'High-$\ell$, Planck 2018')


cnames = ['light blue', 'orange', 'blue', 'red']
for i, ps in enumerate(spectra.keys()):
    if i == 0:
        axs[1].plot(ells, spectra[ps]['TT'], label='Total', color='black')
        continue
    if i == 1:
        continue
    axs[1].plot(ells, spectra[ps]['TT'], label=ps, color=col[cnames[i-2]], linewidth=1.5, linestyle='dotted')


for ax in axs:
    ax.set_xlim(2, 3e3)
    ax.set_xscale('log')
    ax.set_xlabel(r'Multipole $\ell$')
    ax.set_ylabel(r'$\ell(\ell+1)C_\ell/2\pi\;[\mu\mathrm{K}^2]$')

axs[0].set_ylim(-2e2, 6.5e3)

axins = axs[0].inset_axes([0.80, 0.6, 0.2, 0.4], transform=axs[0].transAxes) 
axins.plot(ells, spectra['CMB power-spectrum']['TT'], color=col['red'])
axins.plot(ells, spectra['Master level PS']['TT'], color=col['orange'], linestyle='dashed')
axins.errorbar(high_ell_TT_data[:,0], high_ell_TT_data[:,1], 
            yerr=[high_ell_TT_data[:,3], high_ell_TT_data[:,2]], markersize=5,
            fmt='o', markerfacecolor='none', color='black', linewidth=1, alpha=0.5,
            markeredgecolor='black', capsize=3, label=r'High-$\ell$ TT (Planck 2018)')
axins.set_xscale('log')

x1, x2, y1, y2 = 1500, 2200, 100, 700
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
# axins.set_xticklabels([])
# axins.set_yticklabels([])
axins.set_xticklabels([], minor=True)
axins.set_yticklabels([], minor=True)

axs[0].indicate_inset_zoom(axins, edgecolor="black")

axs[0].legend(loc='upper left', fontsize=13)
axs[1].legend(loc='lower left', fontsize=13, ncols=1)
axs[1].set_yscale('log')

plt.figtext(0.12, 0.02, r'$\dagger$ Simulation run \textit{without} helium, reionisation, neutrinos, nor polarisation.', 
            horizontalalignment='left', fontsize=10, color='gray')
fig.subplots_adjust(hspace=0)
fig.savefig('../figures/CMB_power_spectrum.pdf')




##########################
# TRANSFER FUNCTIONS (T) #
##########################
ells = [5, 50, 100, 200, 500, 1000]
k_axis = transfer_funcsT['k']

fig, axs = plt.subplots(2, 1, figsize=(6,6.5), sharex=True)
colors = [cmap(i) for i in np.linspace(1, 0.7, len(ells))]
for i, ell in enumerate(ells):
    axs[0].plot(k_axis*(c/H0), np.sqrt(ell*(ell+1))*transfer_funcsT[ell], color=colors[i], alpha=0.65)
    axs[1].plot(k_axis*(c/H0), ell*(ell+1)*abs(transfer_funcsT[ell])**2/k_axis, color=colors[i], alpha=0.65, label=r'$\ell=$' + f' {ell}')

x_text, y_text = 0.9, 0.85
for i, ax in enumerate(axs.flat):
    ax.text(x_text, y_text, f'({chr(i+97)})', fontweight='bold', transform=ax.transAxes)

fig.legend(fontsize=14, loc='lower center', ncols=3)
axs[0].set_xlim(0, 600)
axs[1].set_yscale('symlog', linthresh=3, linscale=1)
axs[1].set_ylim(-0.02, 200)

axs[1].set_xlabel(r'$k\eta_0$')
axs[0].set_ylabel(r'$\sqrt{\ell(\ell+1)}\Theta_\ell$')
axs[1].set_ylabel(r'$\ell(\ell+1)|\Theta_\ell|^2/k$')

fig.subplots_adjust(hspace=0, bottom=0.2)
fig.savefig('../figures/transfer_funcsT.pdf')



###############################
# POLARISATION POWER SPECTRA #
###############################
fig, axs = plt.subplots(2, 1, figsize=(9,6), sharex=True)

full_ps = 'CMB power-spectrum'
ells = spectra[full_ps]['ells']
axs[0].plot(ells, spectra[full_ps]['TE'], color=col['red'], label=r'Cross PS (theory)')

axs[0].errorbar(high_ell_TE_data[:,0], high_ell_TE_data[:,1], 
            yerr=[high_ell_TE_data[:,3], high_ell_TE_data[:,2]], markersize=4,
            fmt='x', markerfacecolor='none', markeredgewidth=0.5, color='black', linewidth=1,
            markeredgecolor='black', capsize=3, label=r'Planck 2018 (TE)')

axs[1].plot(ells, spectra[full_ps]['EE'], color=col['red'], label=r'E-mode PS (theory)')

# Fix normalisation to match data
ells = high_ell_EE_data[:,0]
high_ell_EE_data[:,1]  *= 2.0*np.pi*1e5/(ells*(ells+1.0))
high_ell_EE_data[:,2]  *= 2.0*np.pi*1e5/(ells*(ells+1.0))
high_ell_EE_data[:,3]  *= 2.0*np.pi*1e5/(ells*(ells+1.0))

axs[1].errorbar(high_ell_EE_data[:,0], high_ell_EE_data[:,1], 
            yerr=[high_ell_EE_data[:,3], high_ell_EE_data[:,2]], markersize=4,
            fmt='x', markerfacecolor='none', markeredgewidth=0.5, color='black', linewidth=1,
            markeredgecolor='black', capsize=3, label=r'Planck 2018 (EE)')

axs[0].set_ylabel(r'$\ell(\ell+1)C^\mathrm{TE}_\ell/2\pi\;[\mu\mathrm{K}^2]$')
axs[1].set_ylabel(r'$C^\mathrm{EE}_\ell\;[10^{-5}\mu\mathrm{K}^2]$')
axs[1].set_xlabel(r'Multipole $\ell$')
axs[1].set_xlim(0,2000)

for ax in axs.flat:
    ax.legend()

fig.subplots_adjust(hspace=0)
fig.savefig('../figures/polarisation_spectra.pdf')