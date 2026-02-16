import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import norm

# COLORS
cblue = "#6D79F8"
cred = "#DE4E61"
cgreen = "#32A639"

# Observational data
chi2, h, OmegaM, OmegaK = np.loadtxt('results_supernovafitting.txt', skiprows=1000, unpack=True)
OmegaLambda = 1 - OmegaM - OmegaK

# Data from best fit run
lna, eta, Hp, dHpdx, ddHpddx, OmegaB, OmegaCDM, OmegaL, OmegaR, OmegaNu, OmegaK, lumdist = np.loadtxt('cosmology.txt', unpack=True) 

# Initial density parameters
OmegaM0 = 0.05 + 0.205027
OmegaR0 = 3.47427e-05 + 5.0223e-05
OmegaL0 = 0.665937

# Equality times
x_MR_eq = np.log(OmegaR0/OmegaM0)
x_LM_eq = np.log(OmegaM0/OmegaL0)/3

# Data from fiducial cosmological parameters run
fidlna, fideta, fidHp, fiddHpdx, fidddHpddx, fidOmegaB, fidOmegaCDM, fidOmegaL, fidOmegaR, fidOmegaNu, fidOmegaK, fidlumdist = np.loadtxt('fiducial_cosmology.txt', unpack=True) 

# Minimal chi^2
min_index = np.argmin(chi2)
onesigma = chi2 < chi2[min_index] + 3.53
twosigma = chi2 < chi2[min_index] + 8.02

# Scatter plot OmegaM vs. OmegaLambda
fig, ax = plt.subplots(1,1,figsize=(6,5))
ax.scatter(OmegaM[twosigma], OmegaLambda[twosigma], color=cblue, label=r'$2\sigma$ constraint')
ax.scatter(OmegaM[onesigma], OmegaLambda[onesigma], color=cred, label=r'$1\sigma$ constraint')
ax.scatter(OmegaM[min_index], OmegaLambda[min_index], edgecolor='black', marker='s', facecolor='none', label=r'Best fit')
x = np.linspace(0,1,100)
ax.plot(x, 1-x, linestyle='--', color='black', label='Flat Universe')
ax.set_xlim(0, 1)
ax.set_ylim(0, 1.2)
ax.set_xlabel(r'$\Omega_{M0}$')
ax.set_ylabel(r'$\Omega_{\Lambda0}$')
ax.legend()

plt.savefig('figures/supernova_fitting.pdf')


# Histogram
fig, axs = plt.subplots(2,1,figsize=(7,10))
counts, bins, patches = axs[0].hist(OmegaLambda, bins=30, color='black', alpha=0.3, density=True)
axs[0].set_xlabel(r'$\Omega_\Lambda$')
axs[0].set_ylabel(r'Probability density')
axs[0].set_title(r'Posterior for $\Omega_\Lambda$')

# Fit Gaussian
mu = np.mean(OmegaLambda)
sigma = np.std(OmegaLambda)
x_fit = np.linspace(OmegaLambda.min(), OmegaLambda.max(), 200)
axs[0].plot(x_fit, norm.pdf(x_fit, mu, sigma), color=cred, linewidth=2)
axs[0].axvline(mu, linestyle='--', color='black', label=r'$\langle\Omega_\Lambda\rangle=$'+f' {mu:.3f}')
axs[0].legend()

counts, bins, patches = axs[1].hist(h, bins=30, color='black', alpha=0.3, density=True)
axs[1].set_xlabel(r'$h$')
axs[1].set_ylabel(r'Probability density')
axs[1].set_title(r'Posterior for $h$')

# Fit Gaussian
mu = np.mean(h)
sigma = np.std(h)
x_fit = np.linspace(h.min(), h.max(), 200)
axs[1].plot(x_fit, norm.pdf(x_fit, mu, sigma), color=cred, linewidth=2)
axs[1].axvline(mu, linestyle='--', color='black', label=r'$\langle h\rangle=$'+f' {mu:.3f}')
axs[1].axvline(0.67, linestyle='-', color=cblue, label=r'Fiducial $h=0.67$')
axs[1].legend(loc='upper right')
axs[1].set_xlim(0.668, 0.735)

fig.subplots_adjust(hspace=0.3)
plt.savefig('figures/bestfit_OmegaLambda.pdf')


# Testing the model
fig, axs = plt.subplots(2,2,figsize=(12,8))
axs[0,0].semilogy(lna, Hp*3.08567758e17, color='black')
axs[0,0].set_title(r'$\mathcal{H}(x)\;\;\bigg(\frac{100\;\mathrm{km/s}}{\mathrm{Mpc}}\bigg)$')
axs[0,0].set_xlim(-20,5)
# axs[0,0].set_ylim(1e-1,2e3)
axs[0,0].set_xlabel(r'$x= \ln a$')

axs[0,1].semilogy(lna, eta/3.08567758e22, color='black')
axs[0,1].set_title(r'$\eta(x)\; (\mathrm{Mpc})$')
axs[0,1].set_xlim(-12,0)
axs[0,1].set_ylim(1e0,1e5)
axs[0,1].set_xlabel(r'$x=\ln a$')

axs[1,0].plot(lna, eta*Hp/2.99792458e8, color='black')
axs[1,0].set_title(r'$\eta\mathcal{H}/c$')
axs[1,0].set_xlim(-15,0)
axs[1,0].set_ylim(0.75,3)
axs[1,0].set_xlabel(r'$x=\ln a$')

axs[1,1].plot(lna, OmegaR + OmegaNu, color=cblue, label=r'$\Omega_R = \Omega_\gamma + \Omega_\nu$')
axs[1,1].plot(lna, OmegaB + OmegaCDM, color=cred, label=r'$\Omega_M = \Omega_b + \Omega_\mathrm{CDM}$')
axs[1,1].plot(lna, OmegaL, color=cgreen, label=r'$\Omega_\Lambda$')
axs[1,1].axvline(x_MR_eq, color='black', linestyle=':')
axs[1,1].axvline(x_LM_eq, color='black', linestyle=':')
axs[1,1].set_title(r'$\Omega_i(x)$')
axs[1,1].set_xlim(-20,5)
axs[1,1].set_xlabel(r'$x=\ln a$')
axs[1,1].legend(loc='center left', fontsize=11.5)

fig.subplots_adjust(hspace=0.35, right=0.95, left=0.05)

plt.savefig('figures/test_background_cosmology.pdf')


# Derivatives of the conformal Hubble function
fig, axs = plt.subplots(2,1,figsize=(6,5))
axs[0].plot(lna, dHpdx/Hp, color=cblue)
axs[0].tick_params(which='both', bottom=True, left=True, top=True, right=True, direction='in')
axs[0].axvline(x_MR_eq, color='black', linestyle=':')
axs[0].axvline(x_LM_eq, color='black', linestyle=':')
axs[0].set_ylabel(r"$\mathcal{H}'/\mathcal{H}$")
axs[0].set_xlim(-20,5)
axs[0].label_outer()
axs[0].text(0.13, 0.6, 'Radiation', transform=axs[0].transAxes, 
            verticalalignment='top')
axs[0].text(0.56, 0.6, 'Matter', transform=axs[0].transAxes, 
            verticalalignment='top')
axs[0].text(0.88, 0.6, r'$\Lambda$', transform=axs[0].transAxes, 
            verticalalignment='top')

axs[1].plot(lna, ddHpddx/Hp, color=cred)
axs[1].tick_params(which='both', bottom=True, left=True, top=True, right=True, direction='in')
axs[1].axvline(x_MR_eq, color='black', linestyle=':')
axs[1].axvline(x_LM_eq, color='black', linestyle=':')
axs[1].set_ylabel(r"$\mathcal{H}''/\mathcal{H}$")
axs[1].set_xlabel(r'$x=\ln a$')
axs[1].set_xlim(-20,5)
fig.subplots_adjust(hspace=0)

plt.savefig('figures/Hp_derivatives.pdf')


# Luminosity distance
fig, ax = plt.subplots(1,1,figsize=(6,5))
z_obs, dL_obs, sigma = np.loadtxt('data/supernovadata.txt', unpack=True)

z = np.exp(-lna)-1
ax.plot(z, lumdist/3.08567758e25/z, color=cblue, label=r'Best fit')

fidz = np.exp(-fidlna)-1
ax.plot(fidz, fidlumdist/3.08567758e25/fidz, color=cred, label=r'Fiducial')

ax.errorbar(z_obs, dL_obs/z_obs, yerr=sigma/z_obs, fmt='o', 
            markerfacecolor='none', color='black', linewidth=0.3,
            markeredgecolor='black', capsize=3, label=r'Supernovae')
ax.set_ylabel(r'Luminosity distance $d_L(z)/z\;(\mathrm{Gpc})$')
ax.set_xlabel(r'Redshift $z$')

ax.set_xscale('log')
ax.set_xlim(9e-3, 1.4)
ax.set_ylim(3, 8)
ax.legend(loc='upper left')
fig.savefig('figures/dL_obs_vs_pred.pdf')