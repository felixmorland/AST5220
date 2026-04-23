import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import healpy as hp

custom_col = ["#476CFF", "#32A639", "#ECC723", "#EC7023", "#DE4E61"]

# Planck-like colormap
from matplotlib.colors import ListedColormap
colombi1_cmap = ListedColormap(np.loadtxt("Planck_Parchment_RGB.txt")/255.)
colombi1_cmap.set_bad("gray") # color of missing pixels
colombi1_cmap.set_under("white") 
cmap = colombi1_cmap

def CMB_map():
    #====================================================
    # Load data
    #====================================================
    data  = np.loadtxt('power_spectrum_data/cells.txt')
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
    seed      = 42

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
    hp.mollview(map_T, title="CMB Temperature (TT)", min=-500, max=500, **plot_kwargs)
    # hp.mollview(map_Q, title="CMB E-polarisation (EE)", **plot_kwargs)
    # hp.mollview(map_U, title="CMB Cross (TE)", **plot_kwargs)

    plt.savefig('figures/simulated_CMB_map.pdf')

# Power spectra
filenames = {
    'CMB power-spectrum'    : 'power_spectrum_data/cells.txt',
    'POL'                   : 'power_spectrum_data/Polarization_ps.txt',
    'Doppler'               : 'power_spectrum_data/Doppler_ps.txt',
    'ISW'                   : 'power_spectrum_data/ISW_ps.txt',
    'SW'                    : 'power_spectrum_data/SW_ps.txt'
}
spectra = {}
for file in filenames.keys():
    ells, cell_TT, cell_EE, cell_TE = np.loadtxt(filenames[file], unpack=True)
    spectra[file] = {
        'ells'  : ells,
        'TT'    : cell_TT,
        'EE'    : cell_EE,
        'TE'    : cell_TE
    }
k_array, matter_ps = np.loadtxt('power_spectrum_data/matter_ps.txt', unpack=True)

# Observational data
low_ell_TT_data = np.loadtxt('data/low_ell_TT.txt', skiprows=1)
high_ell_TT_data = np.loadtxt('data/high_ell_TT.txt', skiprows=1)
high_ell_EE_data = np.loadtxt('data/high_ell_EE.txt', skiprows=1)
high_ell_TE_data = np.loadtxt('data/high_ell_TE.txt', skiprows=1)

SDSS_gal = np.loadtxt('data/SDSS_DR7_LRG.txt', skiprows=1)
WMAP_ACT = np.loadtxt('data/WMAP_ACT.txt', skiprows=1)

# Transfer functions
tf_data = np.genfromtxt('power_spectrum_data/transfer_func.txt')
transfer_funcs = {'k':tf_data[1:,0]}
for ell in range(1, tf_data.shape[1]):
    transfer_funcs[int(tf_data[0,ell])] = tf_data[1:,ell]


###############################
# MATTER POWER SPECTRUM       #
###############################
fig, ax = plt.subplots(1, 1, figsize=(8,6))

h = 0.67
ax.plot(k_array/h, matter_ps/(h**3), label=r'Theory')

ax.errorbar(SDSS_gal[:,0], SDSS_gal[:,1], 
            yerr=SDSS_gal[:,2], markersize=3,
            fmt='o', markerfacecolor='none', color='black', linewidth=1,
            markeredgecolor='black', capsize=0, label=r'SDSS Galaxies (DR7 LRG)')

ax.errorbar(WMAP_ACT[:,0], WMAP_ACT[:,1], 
            yerr=WMAP_ACT[:,2], markersize=3,
            fmt='^', markerfacecolor='none', color='black', linewidth=1,
            markeredgecolor='black', capsize=3, label=r'CMB (WMAP + ACT)')


ax.set_xscale('log')
ax.set_yscale('log')
ax.legend()

fig.savefig('figures/matter_power_spectrum.pdf')

###############################
# CMB POWER SPECTRUM          #
###############################
fig, axs = plt.subplots(2, 1, figsize=(7,8), height_ratios=[0.7, 0.3], sharex=True)
colors = [custom_col[1]] + 4*['black']
widths = [2] + 4*[1]
linestyles = ['-', (0,(1,1)), (0,(5,5)), (0,(1,5)), (0,(3,5,1,5))]
markers = ['none', 'none', 's', 'o', '^']
ells = spectra['POL']['ells']

axs[0].plot(ells, spectra['CMB power-spectrum']['TT'], color=custom_col[4], label=r'CMB power-spectrum (TT)')

axs[0].errorbar(low_ell_TT_data[:,0], low_ell_TT_data[:,1], 
            yerr=[low_ell_TT_data[:,3], low_ell_TT_data[:,2]], markersize=5,
            fmt='o', markerfacecolor='none', color='blue', linewidth=1,
            markeredgecolor='blue', capsize=3, label=r'Low-$\ell$ TT')

axs[0].errorbar(high_ell_TT_data[:,0], high_ell_TT_data[:,1], 
            yerr=[high_ell_TT_data[:,3], high_ell_TT_data[:,2]], markersize=5,
            fmt='^', markerfacecolor='none', color='black', linewidth=1,
            markeredgecolor='black', capsize=0, label=r'High-$\ell$ TT')

indices = np.logspace(0, np.log10(len(ells)-1), 10).astype(int)
for i, ps in enumerate(spectra.keys()):
    if i == 0:
        continue
    axs[1].plot(ells, spectra[ps]['TT'], label=ps, color=colors[i], linestyle=linestyles[i], linewidth=widths[i], marker=markers[i], markevery=indices, markersize=3)

for ax in axs:
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel(r'$\ell$')
    ax.set_ylabel(r'$\ell(\ell+1)C_\ell/2\pi$')

axs[0].legend(loc='lower left', fontsize=13)
axs[1].legend(fontsize=11.5)
fig.subplots_adjust(hspace=0)
fig.savefig('figures/CMB_power_spectrum.pdf')




##########################
# TRANSFER FUNCTIONS     #
##########################
c   = 2.99792458e5  # km/s
H0  = 67            # km/s/Mpc
ells = [6, 100, 200, 500, 1000]
k_axis = transfer_funcs['k']*(c/H0)

fig, ax = plt.subplots(1,1)

for i, ell in enumerate(ells):
    ax.plot(k_axis, transfer_funcs[ell], label=r'$\ell=$' + f' {ell}', color=custom_col[i])

ax.legend()
ax.set_xlim(0, 500)
ax.set_ylim(-0.02, 0.02)
fig.savefig('figures/transfer_funcs.pdf')

###############################
# POLARISATION POWER SPECTRA #
###############################
fig, axs = plt.subplots(2,1)

full_ps = 'CMB power-spectrum'
ells = spectra[full_ps]['ells']
axs[0].plot(ells, spectra[full_ps]['TE'])
axs[1].plot(ells, spectra[full_ps]['EE'])

fig.savefig('figures/polarisation_spectra.pdf')