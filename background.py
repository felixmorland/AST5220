import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import norm
from matplotlib.lines import Line2D
from scipy.spatial import ConvexHull
from matplotlib.ticker import LogLocator

# COLORS
cblue = "#476CFF"
cred = "#DE4E61"
cgreen = "#32A639" 

# Observational data
chi2, h, OmegaM, OmegaK = np.loadtxt('results_supernovafitting.txt', skiprows=1000, unpack=True)
OmegaLambda = 1 - OmegaM - OmegaK

# Data from best fit run
bflna, bfeta, bfHp, bfdHpdx, bfddHpddx, bfOmegaB, bfOmegaCDM, bfOmegaL, bfOmegaR, bfOmegaNu, bfOmegaK, bflumdist, bfangdist, bft = np.loadtxt('best_fit_cosmology.txt', unpack=True) 

# Initial density parameters
OmegaM0 = 0.05 + 0.267
OmegaR0 = 3.81093e-05 + 5.50896e-05
OmegaL0 = 0.682907

# Equality times
x_MR_eq = np.log(OmegaR0/OmegaM0)
x_LM_eq = np.log(OmegaM0/OmegaL0)/3
x_acc = np.log(OmegaM0/(2*OmegaL0))/3

# Data from fiducial cosmological parameters run
lna, eta, Hp, dHpdx, ddHpddx, OmegaB, OmegaCDM, OmegaL, OmegaR, OmegaNu, OmegaK, lumdist, angdist, t = np.loadtxt('fiducial_cosmology.txt', unpack=True) 

# Minimal chi^2
min_index = np.argmin(chi2)
onesigma = chi2 < chi2[min_index] + 3.53
twosigma = chi2 < chi2[min_index] + 8.02

# Scatter plot OmegaM vs. OmegaLambda

# ax.scatter(OmegaM[twosigma], OmegaLambda[twosigma], color=cblue, label=r'$2\sigma$ constraint')
# ax.scatter(OmegaM[onesigma], OmegaLambda[onesigma], color=cred, label=r'$1\sigma$ constraint')

fig, ax = plt.subplots(1, 1, figsize=(6, 5))

# Create convex hulls from the scatter points
if len(OmegaM[twosigma]) > 2:
    points_2sigma = np.column_stack([OmegaM[twosigma], OmegaLambda[twosigma]])
    hull_2sigma = ConvexHull(points_2sigma)
    
    # Plot filled region for 2-sigma
    ax.fill(points_2sigma[hull_2sigma.vertices, 0], 
            points_2sigma[hull_2sigma.vertices, 1],
            color='blue', alpha=0.3, label=r'$2\sigma$ constraint')

if len(OmegaM[onesigma]) > 2:
    points_1sigma = np.column_stack([OmegaM[onesigma], OmegaLambda[onesigma]])
    hull_1sigma = ConvexHull(points_1sigma)
    
    # Plot filled region for 1-sigma
    ax.fill(points_1sigma[hull_1sigma.vertices, 0], 
            points_1sigma[hull_1sigma.vertices, 1],
            color=cred, alpha=0.8, label=r'$1\sigma$ constraint')

x = np.linspace(0,1,100)
ax.plot(x, 1-x, linestyle='--', color='black', label='Flat Universe')
ax.scatter(OmegaM[min_index], OmegaLambda[min_index], edgecolor='black', marker='s', facecolor='none', label=r'Best-fit')
ax.scatter(OmegaM0, OmegaL0, edgecolor='black', marker='o', s=70, facecolor='none', label=r'Fiducial')
ax.set_xlim(0, 1)
ax.set_ylim(0, 1.2)
ax.set_xlabel(r'$\Omega_{M0}$')
ax.set_ylabel(r'$\Omega_{\Lambda0}$')
ax.legend()

plt.savefig('figures/supernova_fitting.pdf')


# Histogram
fig, axs = plt.subplots(2,1,figsize=(7,10))
counts, bins, patches = axs[0].hist(OmegaLambda, bins=30, color='black', alpha=0.3, density=True)
axs[0].set_xlabel(r'$\Omega_{\Lambda0}$')
axs[0].set_ylabel(r'Probability density')
axs[0].set_title(r'Posterior for $\Omega_{\Lambda0}$')
axs[0].set_xlim(0.1, 1.3)

# Fit Gaussian
mu = np.mean(OmegaLambda)
sigma = np.std(OmegaLambda)
x_fit = np.linspace(OmegaLambda.min(), OmegaLambda.max(), 200)
axs[0].plot(x_fit, norm.pdf(x_fit, mu, sigma), color=cred, linewidth=2)
axs[0].axvline(OmegaLambda[min_index], linestyle='--', color='black', label=r'Best-fit $\Omega_{\Lambda0}=$'+f' {OmegaLambda[min_index]:.3f}')
axs[0].axvline(OmegaL0, linestyle=':', color='black', label=r'Fiducial $\Omega_{\Lambda0}=0.683$')
axs[0].legend(loc='upper right')

counts, bins, patches = axs[1].hist(h, bins=30, color='black', alpha=0.3, density=True)
axs[1].set_xlabel(r'$h$')
axs[1].set_ylabel(r'Probability density')
axs[1].set_title(r'Posterior for $h$')

# Fit Gaussian
mu = np.mean(h)
sigma = np.std(h)
x_fit = np.linspace(h.min(), h.max(), 200)
axs[1].plot(x_fit, norm.pdf(x_fit, mu, sigma), color=cred, linewidth=2)
axs[1].axvline(0.701711, linestyle='--', color='black', label=r'Best-fit $ h=$'+f' {0.701711:.3f}')
axs[1].axvline(0.67, linestyle=':', color='black', label=r'Fiducial $h=0.67$')
axs[1].legend(loc='upper right')
axs[1].set_xlim(0.65, 0.75)

fig.subplots_adjust(hspace=0.3)
plt.savefig('figures/bestfit_OmegaLambda.pdf')



##################
# Conformal time #
##################
fig, axs = plt.subplots(2,1,figsize=(6,6))

axs[0].semilogy(lna, eta/3.08567758e22, color='black')
axs[0].set_ylabel(r'$\eta(x)\; [\mathrm{Mpc}]$')
axs[0].set_xlim(-18,0)
axs[0].set_ylim(2e-3,1e5)
axs[0].label_outer()

axs[1].axhline(1, linestyle='-.', color='black')
axs[1].plot(lna, eta*Hp/2.99792458e8, color=cred)
axs[1].set_ylabel(r'$\eta\mathcal{H}/c$')
axs[1].set_xlim(-18,0)
axs[1].set_ylim(0.75,3.3)
axs[1].set_xlabel(r'$x=\ln a$')

fig.subplots_adjust(hspace=0)

plt.savefig('figures/conformal_time.pdf')


################################################
# Derivatives of the conformal Hubble function #
################################################

fig, axs = plt.subplots(2,1,figsize=(6,6))
axs[0].semilogy(lna, Hp*3.08567758e17, color='black')
axs[0].set_ylabel(r'$\mathcal{H}(x)\;\;\tiny[100\;\mathrm{km/s}/\mathrm{Mpc}]$', fontsize=12)
axs[0].set_xlim(-15,5)
axs[0].set_ylim(3e-1, 1e5)
axs[0].yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=10))
axs[0].axvline(x_MR_eq, color='black', linestyle=(0, (3, 1, 1, 1, 1, 1)), linewidth=1)
axs[0].axvline(x_LM_eq, color='black', linestyle=':', linewidth=1)
axs[0].axvline(x_acc, color='black', linestyle='dashed', linewidth=1)
axs[0].label_outer()

axs[1].plot(lna, dHpdx/Hp, color=cblue, label=r"$\mathcal{H}'/\mathcal{H}$")
axs[1].plot(lna, ddHpddx/Hp, color=cred, label=r"$\mathcal{H}''/\mathcal{H}$")
axs[1].axvline(x_MR_eq, color='black', linestyle=(0, (3, 1, 1, 1, 1, 1)), linewidth=1)
axs[1].axvline(x_LM_eq, color='black', linestyle=':', linewidth=1)
axs[1].axvline(x_acc, color='black', linestyle='dashed', linewidth=1)
axs[1].legend(loc='lower right')
axs[1].set_xlim(-15,5)
axs[1].set_xlabel(r'$x=\ln a$')


axs[0].text(0.08, 1.15, 'Radiation', transform=axs[0].transAxes, 
            verticalalignment='top')
axs[0].text(0.47, 1.15, 'Matter', transform=axs[0].transAxes, 
            verticalalignment='top')
axs[0].text(0.85, 1.15, r'$\Lambda$', transform=axs[0].transAxes, 
            verticalalignment='top')

line1 = Line2D([], [], color='black', linestyle=(0, (3, 1, 1, 1, 1, 1)), label=r'$x_\mathrm{eq}^{MR}=-8.132$')
line2 = Line2D([], [], color='black', linestyle='dashed', label=r'$x_\mathrm{acc}=-0.487$')
line3 = Line2D([], [], color='black', linestyle=':', label=r'$x_\mathrm{eq}^{\Lambda R}=-0.256$')

fig.legend(ncols=3, loc='lower center', handles=[line1, line2, line3], fontsize=14, frameon=True)
fig.subplots_adjust(hspace=0, bottom=0.20)

plt.savefig('figures/Hp_derivatives.pdf')


######################
# Density parameters #
######################
fig, ax = plt.subplots(1, 1, figsize=(6,5))
ax.plot(lna, OmegaR + OmegaNu, color=cblue, label=r'$\Omega_R$')
ax.plot(lna, OmegaB + OmegaCDM, color=cred, label=r'$\Omega_M$')
ax.plot(lna, OmegaL, color=cgreen, label=r'$\Omega_\Lambda$')
ax.axvline(x_MR_eq, color='black', linestyle=(0, (3, 1, 1, 1, 1, 1)))
ax.axvline(x_LM_eq, color='black', linestyle=':')

handles, labels = ax.get_legend_handles_labels()

line1 = Line2D([], [], color='black', linestyle=(0, (3, 1, 1, 1, 1, 1)), label=r'$x_\mathrm{eq}^{MR}$')
line3 = Line2D([], [], color='black', linestyle=':', label=r'$x_\mathrm{eq}^{\Lambda R}$')

handles += [line1, line3]
labels += [r'$x_\mathrm{eq}^{MR}$', r'$x_\mathrm{eq}^{\Lambda R}$']

ax.legend(loc='center left', handles=handles, labels=labels, bbox_to_anchor=(0.05, 0.5))
ax.set_xlim(-20,5)
ax.set_xlabel(r'$x=\ln a$')
ax.set_ylabel(r'$\Omega_i(x)$')

fig.subplots_adjust(bottom=0.2)

plt.savefig('figures/density_parameters.pdf')



# Luminosity distance
fig, ax = plt.subplots(1,1,figsize=(6,5))
z_obs, dL_obs, sigma = np.loadtxt('data/supernovadata.txt', unpack=True)

bfz = np.exp(-bflna)-1
ax.plot(bfz, bflumdist/3.08567758e25/bfz, color=cblue, label=r'Best-fit')

z = np.exp(-lna)-1
ax.plot(z, lumdist/3.08567758e25/z, color=cred, label=r'Fiducial')

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

# Time evolution
fig, ax = plt.subplots(1, 1, figsize=(6,5))

ratex = np.linspace(-6.5,-3, 2)
ax.plot(ratex, np.exp(1.5*ratex + 5), color='gray', linewidth=0.8)
ax.text(0.11, 0.74, r'$\ln(t)\propto 3x/2$', transform=ax.transAxes, 
            verticalalignment='top', rotation=45, color='gray')

ax.semilogy(lna, t/(1e9*365.*24.*60.*60.), color='black')
ax.set_xlabel(r'$x = \ln a$')
ax.set_ylabel(r'$t(x)$ [Gyr]')
ax.set_xlim(-8,5)
ax.set_ylim(1e-5,100)
# ax.axvline(x_acc, color=cred, linestyle='dashed', label=r'$x_\mathrm{acc}=$'+f' {x_acc:.3f}')
ax.axvline(x_LM_eq, color='black', linestyle='dashed', label=r'$x_\mathrm{eq}^{\Lambda M}=$'+f' {x_LM_eq:.3f}', linewidth=0.8)
ax.tick_params(which='both', bottom=True, left=True, top=True, right=True, direction='in')
ax.legend(loc='lower right')

fig.savefig('figures/time_func_of_x.pdf')