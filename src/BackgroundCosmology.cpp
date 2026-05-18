#include "BackgroundCosmology.h"
    
BackgroundCosmology::BackgroundCosmology(
    double h, 
    double OmegaB, 
    double OmegaCDM, 
    double OmegaK,
    double Neff, 
    double TCMB) :
  h(h),
  OmegaB(OmegaB),
  OmegaCDM(OmegaCDM),
  OmegaK(OmegaK),
  Neff(Neff), 
  TCMB(TCMB)
{
  // Derived parameters
  H0 = h * Constants.H0_over_h;
  OmegaR = (
    (M_PI*M_PI/15.) * std::pow(Constants.k_b*TCMB, 4.) / std::pow(Constants.hbar, 3.) / std::pow(Constants.c, 5.) * (8.*M_PI*Constants.G) / (3.*H0*H0)
  );
  OmegaNu = Neff*(7./8.) * std::pow(4./11., 4./3.) * OmegaR;
  OmegaLambda = 1. - OmegaB - OmegaCDM - OmegaR - OmegaNu - OmegaK;
}

/**
 * @brief Solves the background cosmology
 */
void BackgroundCosmology::solve(){
  // x-array
  const double xmin = -20.0;
  const double xmax = 5.0;
  const int npts = 1000;

  Vector x_array = Utils::linspace(xmin, xmax, npts);
  
  // Indicate new milestone
  std::cout << "\n";
  std::cout << "/===========================\\\n";
  std::cout << "|  I. BACKGROUND COSMOLOGY  |\n";
  std::cout << "\\===========================/\n\n";

  std::cout << "Solving...\n";
  // Solve for conformal time eta(x)
  Utils::StartTiming("Eta");
  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){
    detadx[0] = Constants.c / Hp_of_x(x);
    return GSL_SUCCESS;
  };

  // Initial condition for eta
  double eta_ini = Constants.c / Hp_of_x(xmin);
  Vector eta_ic{eta_ini};
  
  // Solve the ODE
  ODESolver ode;
  ode.solve(detadx, x_array, eta_ic);

  // Get and spline the result
  auto eta_array = ode.get_data_by_component(0);
  eta_of_x_spline.create(x_array, eta_array, "Conformal time eta(x)");

  Utils::EndTiming("Eta");

  // Solve for cosmic time
  Utils::StartTiming("Cosmic time");
  // The ODE for deta/dx
  ODEFunction dtdx = [&](double x, const double *t, double *dtdx){
    dtdx[0] = 1 / H_of_x(x);
    return GSL_SUCCESS;
  };

  // Initial condition for eta
  double t_ini = 1 / (2 * H_of_x(xmin));
  Vector t_ic{t_ini};
  
  // Solve the ODE
  ode.solve(dtdx, x_array, t_ic);

  // Get and spline the result
  auto t_array = ode.get_data_by_component(0);
  t_of_x_spline.create(x_array, t_array, "Cosmic time t(x)");

  Utils::EndTiming("Cosmic time");
}


/**
 * @brief Hubble function H(x)
 * @param x the log-scale factor
 */
double BackgroundCosmology::H_of_x(double x) const{

  // Scale factor from x
  double a = exp(x);
  // Hubble function H(a)
  double H = H0 * sqrt(
    (OmegaB + OmegaCDM)*pow(a, -3.) + (OmegaR + OmegaNu)*pow(a, -4.)
    + OmegaK*pow(a, -2.) + OmegaLambda
  );

  return H;
}

/**
 * @brief Conformal Hubble function \mathcal{H}(x)
 * @param x the log-scale factor
 */
double BackgroundCosmology::Hp_of_x(double x) const{
  // Conformal Hubble function
  return exp(x)*H_of_x(x);
}

/**
 * @brief First derivative of the conformal Hubble function
 * @param x the log-scale factor
 */
double BackgroundCosmology::dHpdx_of_x(double x) const{
  // Density parameters
  double M = get_OmegaB(x) + get_OmegaCDM(x);
  double R = get_OmegaR(x) + get_OmegaNu(x);
  double K = get_OmegaK(x);

  double dHpdx = Hp_of_x(x) * (1. - 3.*M/2. - 2.*R - K);

  return dHpdx;
}

/**
 * @brief Second derivative of the conformal Hubble function
 * @param x the log-scale factor
 */
double BackgroundCosmology::ddHpddx_of_x(double x) const{
  // Density parameters
  double M = get_OmegaB(x) + get_OmegaCDM(x);
  double R = get_OmegaR(x) + get_OmegaNu(x);
  double K = get_OmegaK(x);

  double ddHpddx = Hp_of_x(x) * (
    pow(1. - 3.*M/2. - 2.*R - K, 2) - pow(3.*M + 4.*R + 2.*K, 2)/2. 
    + (9.*M + 16.*R + 4.*K)/2.
  );

  return ddHpddx;
}

/**
 * @brief Baryon density parameter
 * @param x the log-scale factor
 * @return OmegaB
 */
double BackgroundCosmology::get_OmegaB(double x) const{ 
  // Baryon density parameter at x
  return OmegaB * H0*H0 / (pow(exp(x), 3) * pow(H_of_x(x), 2));
}
/**
 * @brief Photon density parameter
 * @param x the log-scale factor
 * @return OmegaR (Omega_gamma)
 */
double BackgroundCosmology::get_OmegaR(double x) const{ 
  // Photon density parameter at x
  return OmegaR * H0*H0 / (pow(exp(x), 4) * pow(H_of_x(x), 2));
}
/**
 * @brief Neutrino density parameter
 * @param x the log-scale factor
 * @return OmegaNu
 */
double BackgroundCosmology::get_OmegaNu(double x) const{ 
  // Neutrino density parameter at x
  return OmegaNu * H0*H0 / (pow(exp(x), 4) * pow(H_of_x(x), 2));
}
/**
 * @brief CDM density parameter
 * @param x the log-scale factor
 * @return OmegaCDM
 */
double BackgroundCosmology::get_OmegaCDM(double x) const{ 
  // CDM density parameter at x
  return OmegaCDM * H0*H0 / (pow(exp(x), 3) * pow(H_of_x(x), 2));
}
/**
 * @brief DE density parameter
 * @param x the log-scale factor
 * @return OmegaLambda
 */
double BackgroundCosmology::get_OmegaLambda(double x) const{ 
  // Cosmological constant density parameter at x
  return OmegaLambda * H0*H0 /  pow(H_of_x(x), 2);
}
/**
 * @brief Curvature density parameter
 * @param x the log-scale factor
 * @return OmegaK
 */
double BackgroundCosmology::get_OmegaK(double x) const{ 
  // Curvature density parameter at x
  return OmegaK * H0*H0 / (pow(exp(x), 2) * pow(H_of_x(x), 2));
}
/**
 * @brief Total non-relativistic matter density parameter
 * @param x the log-scale factor
 * @return OmegaM = OmegaB + OmegaCDM
 */
double BackgroundCosmology::get_OmegaM(double x) const {
  // Total non-relativistic matter density parameter at x
  return (OmegaB + OmegaCDM) * H0*H0 / (pow(exp(x), 3) * pow(H_of_x(x), 2));
}
/**
 * @brief Total radiation density parameter
 * @param x the log-scale factor
 * @return OmegaRtot = OmegaR + OmegaNu
 */
double BackgroundCosmology::get_OmegaRtot(double x) const {
  // Total relativistic species density parameter at x
  return (OmegaR + OmegaNu) * H0*H0 / (pow(exp(x), 4) * pow(H_of_x(x), 2));
}
/**
 * @brief Matter and neutrino density parameter
 * @param x the log-scale factor
 * @return OmegaM + OmegaNu
 */
double BackgroundCosmology::get_OmegaMnu(double x) const {
  // Matter + neutrinos density parameter at x
  return get_OmegaM(x) + get_OmegaNu(x);
}
/**
 * @brief Comoving distance
 * @param x the log-scale factor
 * @return \chi
 */   
double BackgroundCosmology::get_comoving_distance_of_x(double x) const{
  // Comoving distance
  return eta_of_x(0.0) - eta_of_x(x);
}

/**
 * @brief FLRW r-coordinate
 * @param x the log-scale factor
 * @return r(x)
 */
double BackgroundCosmology::get_r_coordinate(double x) const {
  // Comoving distance
  double chi = get_comoving_distance_of_x(x);
  // r coordinate
  double r;

  if(OmegaK == 0) {
    r = chi;
  } else if(OmegaK > 0) {
    r = chi * sinh(sqrt(abs(OmegaK)) * H0 * chi / Constants.c) / (sqrt(abs(OmegaK)) * H0 * chi / Constants.c);
  } else {
    r = chi * sin(sqrt(abs(OmegaK)) * H0 * chi / Constants.c) / (sqrt(abs(OmegaK)) * H0 * chi / Constants.c);
  };
  return r;
}

/**
 * @brief Luminosity distance
 * @param x the log-scale factor
 * @return d_L(x)
 */
double BackgroundCosmology::get_luminosity_distance_of_x(double x) const{
  // d_L = r/a
  return get_r_coordinate(x) * exp(-x);
}

/**
 * @brief Angular diameter distance
 * @param x the log-scale factor
 * @return d_A(x)
 */
double BackgroundCosmology::get_angular_diameter_distance_of_x(double x) const{
  // d_A = r*a
  return get_r_coordinate(x) * exp(x);
}

/**
 * @brief Conformal time
 * @param x the log-scale factor
 * @return \eta(x)
 */
double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
}
/**
 * @brief Cosmic time
 * @param x the log-scale factor
 * @return t(x)
 */
double BackgroundCosmology::get_t_of_x(double x) const{
  return t_of_x_spline(x);
}

/**
 * @brief Hubble parameter H0
 */
double BackgroundCosmology::get_H0() const{ 
  return H0; 
}
/**
 * @brief Reduced Hubble parameter h
 */
double BackgroundCosmology::get_h() const{ 
  return h; 
}
/**
 * @brief Effective neutrino number Neff
 */
double BackgroundCosmology::get_Neff() const{ 
  return Neff; 
}
/**
 * @brief CMB temperature
 * @param x the log-scale factor
 * @return T_CMB(x)
 */
double BackgroundCosmology::get_TCMB(double x) const{ 
  if(x == 0.0) return TCMB;
  return TCMB * exp(-x); 
}

/**
 * @brief Print information about the class
 */
void BackgroundCosmology::info() const{ 
  // Event times
  double x_MR_eq = log((OmegaR + OmegaNu)/(OmegaB + OmegaCDM));
  double x_LM_eq = log((OmegaB + OmegaCDM)/OmegaLambda)/3.;
  double x_acc = log((OmegaB + OmegaCDM)/(2.*OmegaLambda))/3.;

  std::cout << "\nParameters...\n";
  std::cout << "OmegaB:      " << OmegaB << "\n";
  std::cout << "OmegaCDM:    " << OmegaCDM << "\n";
  std::cout << "OmegaLambda: " << OmegaLambda << "\n";
  std::cout << "OmegaK:      " << OmegaK << "\n";
  std::cout << "OmegaNu:     " << OmegaNu << "\n";
  std::cout << "OmegaR:      " << OmegaR << "\n";
  std::cout << "Neff:        " << Neff << "\n";
  std::cout << "h:           " << h << "\n";
  std::cout << "TCMB:        " << TCMB << "\n";
  std::cout << "t0:          " << get_t_of_x(0) / Constants.Gyr << "\n";
  std::cout << "eta0:        " << eta_of_x(0) / Constants.c / Constants.Gyr << "\n\n";

  std::cout << "Events timeline...\n";
  std::cout << "x_MR_eq:     " << x_MR_eq << "\n";
  std::cout << "t_MR_eq:     " << get_t_of_x(x_MR_eq) / Constants.Gyr     << "\n";
  std::cout << "z_MR_eq:     " << exp(-x_MR_eq) - 1. << "\n\n";
  
  std::cout << "x_acc:       " << x_acc << "\n";
  std::cout << "t_acc:       " << get_t_of_x(x_acc) / Constants.Gyr << "\n";
  std::cout << "z_acc:       " << exp(-x_acc) - 1. << "\n\n";

  std::cout << "x_LM_eq:     " << x_LM_eq << "\n";
  std::cout << "t_LM_eq:     " << get_t_of_x(x_LM_eq) / Constants.Gyr << "\n";
  std::cout << "z_LM_eq:     " << exp(-x_LM_eq) - 1. << "\n";
  std::cout << std::endl;
} 

/**
 * @brief Write out results to file
 * @param filename location of the outfile
 */
void BackgroundCosmology::output(const std::string filename) const{
  const double x_min = -20.0;
  const double x_max =  5.0;
  const int    n_pts =  1000;
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                                          << " ";
    fp << eta_of_x(x)                                << " ";
    fp << Hp_of_x(x)                                 << " ";
    fp << dHpdx_of_x(x)                              << " ";
    fp << ddHpddx_of_x(x)                            << " ";
    fp << get_OmegaB(x)                              << " ";
    fp << get_OmegaCDM(x)                            << " ";
    fp << get_OmegaLambda(x)                         << " ";
    fp << get_OmegaR(x)                              << " ";
    fp << get_OmegaNu(x)                             << " ";
    fp << get_OmegaK(x)                              << " ";
    fp << get_luminosity_distance_of_x(x)            << " ";
    fp << get_angular_diameter_distance_of_x(x)      << " ";
    fp << get_t_of_x(x)                              << " ";
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

