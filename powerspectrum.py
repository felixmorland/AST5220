import matplotlib.pyplot as plt
import numpy as np
import healpy as hp

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
    nside     = 512        # lmax ~ 2*nside
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
        verbose = False
    )

    #====================================================
    # Plot
    #====================================================
    plot_kwargs = dict(cmap="RdBu_r", unit=r"$\mu K$")

    hp.mollview(map_T, title="CMB Temperature (T)",       **plot_kwargs)
    hp.mollview(map_Q, title="CMB Polarization (Q)",      **plot_kwargs)
    hp.mollview(map_U, title="CMB Polarization (U)",      **plot_kwargs)

    plt.show()

# Power spectra
filenames = {
    'full'              : 'power_spectrum_data/cells.txt',
    'only SW'           : 'power_spectrum_data/SW_ps.txt',
    'only ISW'          : 'power_spectrum_data/ISW_ps.txt',
    'only Doppler'      : 'power_spectrum_data/Doppler_ps.txt',
    'only Polarization' : 'power_spectrum_data/Polarization_ps.txt'
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

# Transfer functions
tf_data = np.genfromtxt('power_spectrum_data/transfer_func.txt')
transfer_funcs = {'k':tf_data[1:,0]}
for ell in range(1, tf_data.shape[1]):
    transfer_funcs[int(tf_data[0,ell])] = tf_data[1:,ell]


###############################
# MATTER POWER SPECTRUM       #
###############################
fig, ax = plt.subplots(1,1)
ax.plot(k_array, matter_ps)

fig.savefig('figures/matter_power_spectrum.pdf')

###############################
# CMB POWER SPECTRUM          #
###############################
fig, ax = plt.subplots(1,1, figsize=(10,8))
for ps in spectra.keys():
    ax.plot(spectra[ps]['ells'], spectra[ps]['TT'], label=ps)


ax.errorbar(low_ell_TT_data[:,0], low_ell_TT_data[:,1], 
            yerr=[low_ell_TT_data[:,3], low_ell_TT_data[:,2]], markersize=3,
            fmt='o', markerfacecolor='none', color='black', linewidth=1,
            markeredgecolor='black', capsize=3, label=r'Low-$\ell$ TT')
ax.errorbar(high_ell_TT_data[:,0], high_ell_TT_data[:,1], 
            yerr=[high_ell_TT_data[:,3], high_ell_TT_data[:,2]], markersize=3,
            fmt='o', markerfacecolor='none', color='black', linewidth=1,
            markeredgecolor='black', capsize=2, label=r'High-$\ell$ TT')
# ax.errorbar(high_ell_EE_data[:,0], high_ell_EE_data[:,1], 
#             yerr=[high_ell_EE_data[:,3], high_ell_EE_data[:,2]], markersize=3,
#             fmt='o', markerfacecolor='none', color='black', linewidth=1,
#             markeredgecolor='black', capsize=2, label=r'High-$\ell$ EE')
ax.set_yscale('log')
ax.set_xscale('log')
ax.legend()

fig.savefig('figures/CMB_power_spectrum.pdf')

##########################
# TRANSFER FUNCTIONS     #
##########################
c   = 2.99792458e5  # km/s
H0  = 67            # km/s/Mpc
ells = [6, 100, 200, 500, 1000]
k_axis = transfer_funcs['k']*(c/H0)

fig, ax = plt.subplots(1,1)

for ell in ells:
    ax.plot(k_axis, transfer_funcs[ell], label=r'$\ell=$'+f' {ell}')

ax.legend()
ax.set_xlim(0, 500)
ax.set_ylim(-0.02, 0.02)
fig.savefig('figures/transfer_funcs.pdf')