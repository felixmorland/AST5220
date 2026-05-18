# Simulating cosmic perturbations and the CMB power spectrum from inflationary initial conditions in $\Lambda$ CDM
### A numerical project in Cosmology II, AST5220
by Felix Morland

MSc Astrophysics @ UiO, Spring 2026

## Description
This is a simple Einstein-Boltzmann solver that I developped as a project in the cosmology course AST5220 at the University of Oslo during the Spring semester 2026. It is based upon the [template made by Hans A. Winther](https://github.com/HAWinther/AST5220-Cosmology/), the guide made by [Petter Callin](https://arxiv.org/abs/astro-ph/0606683), as well as the lectures by [Winther](https://cmb.wintherscoming.no). 

The code assumes a flat FLRW background and adiabatic inflationary initial condition and evolves the linear scalar perturbation equations from deep in the radiation dominated era until today. To get around the stiffness problem of the perturbation system in the Thomson scattering-dominated era, we solve it using the tight coupling approximation, and to recover the higher order moments of temperature multipoles we use the *line-of-sight* integral,

$$\Theta_\ell(k, x=0) = \int_{-\infty}^{0} \tilde{S}(k,x) j_\ell[k(\eta_0-\eta)] \text{d} x,$$

with the source function

$$\tilde{S}(k,x) = \tilde{g}\left[ \Theta_0 + \Psi + \frac{1}{4}\Pi\right] + e^{-\tau} \left[\Psi^\prime-\Phi^\prime\right] \\ 
-\frac{1}{ck}\frac{\text{d}}{\text{d}x}(\mathcal{H}\tilde{g}v_b) + \frac{3}{4c^2k^2} \frac{\text{d}}{\text{d}x}
\left[\mathcal{H}\frac{\text{d}}{\text{d}x} (\mathcal{H}\tilde{g}\Pi)\right],$$

taken from [Seljak and Zaldarriaga (1996)](https://arxiv.org/abs/astro-ph/9603033).

The solver includes helium, neutrinos and polarisation, with the possibility of turning each one of these parameters on or off, and can model the late-stage reionisation of hydrogen and helium. It solves the recombination history of the Universe using the Saha approximation and Peebles ODE method, and the background is evolved according to the regular Friedmann equations

$$H^2(a) \equiv \bigg(\frac{\dot{a}}{a}\bigg)^2 = \frac{8\pi G}{3}\sum_i \rho_i,$$
$$\frac{\ddot{a}}{a} = -\frac{4\pi G}{3}\sum_i (\rho_i + 3p_i).$$

Despite its simplifications, the solver is built on first principles with most of its governing equation accesible to master students. I have thus found the solver's merits to lie in the intuition it builds through first principles, while still producing results in impressive agreement with Planck 2018, SDSS, and Ly $\alpha$ observation (see figures below).

<img src="figures/CMB_power_spectrum.pdf" width="300" height="200" alt="CMB power spectrum">
<img src="figures/matter_power_spectrum.pdf" width="300" height="200" alt="Matter power spectrum">

## How to set up and run the simulation
Other than the standard `C++` library, the code has *one* main dependency; the GSL scientific library. Make sure to have it installed and change the `Makefile` parameters specifying the location of your GSL library:
```
# Paths to GSL library
INC = -I$(HOME)/local/include 
LIBS = -L$(HOME)/local/lib -lgsl -lgslcblas
```

Next you need to specify the cosmology for your simulation. This can be done in `src/Main.cpp` where you set the simulation parameters by passing them to the `SimParams` struct;
+ `h` sets the Hubble factor $H_0 = 100h\text{ km/s/Mpc}$;
+ `OmegaCDM` is the density parameter for cold dark matter $\Omega_\text{CDM}$;
+ `OmegaB` is the baryon density parameter $\Omega_b$;
+ `OmegaK` controls the curvature $\Omega_K$;
+ `Neff` is the effective number of neutrino species $N_\text{eff}$;
+ `TCMB` sets the normalisation of the CMB temperature today $T_\text{CMB0}$;
+ `Yp` controls the helium mass fraction of the Universe $Y_p$;
+ `A_s`, `n_s` and `kpivot_mpc` determine the *primordial power spectrum*

$$\frac{k^3}{2\pi^2}\mathcal{P}_\text{prim.}=A_s\left(\frac{k [\text{Mpc}]}{k_\text{pivot Mpc}}\right)^{n_s-1};$$

+ and `reionisation`, `polarisation` and `neutrinos` are booleans that respectively decide whether or not to include reionisation, polarisation or neutrinos in the model. Note that if you turn off neutrinos, you should also remove them at the background level too, i.e. setting `Neff=0.0`, to omit strange relics in the results.

"Out of the box", the program has all these parameters set to their corresponding $\Lambda$ CDM fiducial values from Planck 2018. You can also set the ranges for the perturbation scales and the time variable $x=\ln a$ by changing `k_min` and `k_max`, as well as `x_start` and `x_end`. 

Once everything is set up and the cosmological parameters are set to your liking, the code is ready to be compiled. Simply run
```
make
``` 
to compile, and can then run the simulation by executing
```
./cmb
```

The results of the simulation should appear in the `results` directory. Some examples on how to plot and present the results are given in the Python-script available in the `python` directory.
