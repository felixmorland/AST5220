#include"RecombinationHistory.h"

//====================================================
// Constructors
//====================================================
   
RecombinationHistory::RecombinationHistory(
    BackgroundCosmology *cosmo, 
    double Yp, bool reionisation) :
  cosmo(cosmo),
  Yp(Yp),
  reionisation(reionisation)
{}

//====================================================
// Do all the solving we need to do
//====================================================

void RecombinationHistory::solve(){
    
  // Compute and spline Xe, ne
  solve_number_density_electrons();
   
  // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx, ...
  solve_for_optical_depth_tau();

  // Sound horizon
  solve_for_sound_horizon();
}

// ====================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
// ====================================================
void RecombinationHistory::solve_number_density_electrons(){

  Utils::StartTiming("Xe");
  
  Vector x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);
  Vector Xe_arr(npts_rec_arrays);
  Vector Xe_saha_arr(npts_rec_arrays);
  Vector Xe_saha_arr_with_He(npts_rec_arrays);
  Vector Xe_peebles_arr(npts_rec_arrays);
  Vector ne_arr(npts_rec_arrays);

  // Reionisation parameters
  double f_He       = Yp / (4.0 - 4.0*Yp);
  double z_reion    = 8.0;
  double dz_reion   = 0.5;
  double z_Hereion  = 3.5;
  double dz_Hereion = 0.5; 
  double y_reion    = pow(1+z_reion, 1.5);
  double dy_reion   = 1.5 * sqrt(1+z_reion) * dz_reion;

  // Background cosmology
  const double H0           = cosmo->get_H0();
  const double OmegaB0      = cosmo->get_OmegaB();       
  const double rho_crit0    = 3.0 * H0*H0 / (8.0*M_PI*Constants.G);

  int switch_idx = npts_rec_arrays;
  for(int i = 0; i < npts_rec_arrays; i++){
    Utils::progressbar(double(i) / double(npts_rec_arrays));

    auto Xe_ne_data_H = electron_fraction_from_saha_equation(x_array[i]);
    Xe_saha_arr[i] = Xe_ne_data_H.first;

    if (Yp != 0) {
      auto Xe_ne_data_He     = electron_fraction_from_saha_equation_with_He(x_array[i]);
      Xe_saha_arr_with_He[i] = Xe_ne_data_He.first;
    }

    if(switch_idx == npts_rec_arrays && Xe_saha_arr[i] < Xe_saha_limit)
      switch_idx = i;
  }

  for(int i = 0; i < switch_idx; i++){
    Xe_arr[i]         = (Yp != 0) ? Xe_saha_arr_with_He[i] : Xe_saha_arr[i];
    Xe_peebles_arr[i] = Xe_arr[i];

    double nb  = OmegaB0 * rho_crit0 / (Constants.m_H * exp(3.0*x_array[i]));
    ne_arr[i]  = Xe_arr[i] * nb;
  }

  if(switch_idx < npts_rec_arrays){
    ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
      return this->rhs_peebles_ode(x, Xe, dXedx);
    };

    double Xe_initial = std::min(
      (Yp != 0) ? Xe_saha_arr_with_He[switch_idx-1] : Xe_saha_arr[switch_idx-1],
      Xe_saha_limit
    );
    Vector Xe_ic{ Xe_initial };

    // Use exact slice of x_array starting at switch_idx-1
    Vector x_peebles(x_array.begin() + switch_idx - 1, x_array.end());

    ODESolver peebles_Xe_ode;
    peebles_Xe_ode.solve(dXedx, x_peebles, Xe_ic);
    auto peebles_Xe_array = peebles_Xe_ode.get_data_by_component(0);

    // Write back, skipping index 0 (ic at switch_idx-1, already filled by Saha)
    for(int i = switch_idx; i < npts_rec_arrays; i++){
      int peebles_idx   = i - (switch_idx - 1);
      Xe_peebles_arr[i] = peebles_Xe_array[peebles_idx];
      Xe_arr[i]         = Xe_peebles_arr[i];

      double nb  = OmegaB0 * rho_crit0 / (Constants.m_H * exp(3.0*x_array[i]));
      ne_arr[i]  = Xe_arr[i] * nb;
    }
  }

  if (reionisation) {
    for(int i = 0; i < npts_rec_arrays; i++){
      if(x_array[i] < -4.0) continue;

      double y       = exp(-1.5 * x_array[i]);
      double z       = exp(-x_array[i]) - 1;
      double Xe_base = Xe_peebles_arr[i];   // always use clean Peebles/Saha value

      Xe_arr[i] = Xe_base
                + (1.0+f_He)/2.0 * (1.0 + tanh((y_reion  - y) / dy_reion ))
                +       f_He/2.0 * (1.0 + tanh((z_Hereion - z) / dz_Hereion));

      double nb  = OmegaB0 * rho_crit0 / (Constants.m_H * exp(3.0*x_array[i]));
      ne_arr[i]  = Xe_arr[i] * nb;
    }
  }

  Vector log_ne_arr(npts_rec_arrays);
  Vector log_Xe_arr(npts_rec_arrays);
  for(int i = 0; i < npts_rec_arrays; i++){
    log_Xe_arr[i] = (Xe_arr[i] > 0) ? log(Xe_arr[i]) : -1e100;
    log_ne_arr[i] = (ne_arr[i] > 0) ? log(ne_arr[i]) : -1e100;
  }

  log_ne_of_x_spline.create(x_array, log_ne_arr, "Spline log(ne(x))");
  log_Xe_of_x_spline.create(x_array, log_Xe_arr, "Spline log(Xe(x))");
  Xe_saha_of_x_spline.create(x_array, Xe_saha_arr, "Spline Saha Xe(x)");

  Utils::EndTiming("Xe");
}



//====================================================
// Solve the Saha equation to get ne and Xe
//====================================================
std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
  const double a           = exp(x);
 
  // Physical constants
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double epsilon_0   = Constants.epsilon_0;
  const double H0_over_h   = Constants.H0_over_h;

  // Fetch cosmological parameters
  const double OmegaB      = cosmo->get_OmegaB();
  const double TCMB        = cosmo->get_TCMB();
  const double H0          = cosmo->get_H0();

  // Baryon temperature
  const double Tb          = TCMB / a;
  
  // Critical density and number density of baryons
  const double rho_crit    = 3.0 * H0*H0 / (8.0*M_PI*G);
  const double nb          = OmegaB * rho_crit / (m_H * pow(a, 3));

  // RHS of the Saha equation
  double A = (1.0/nb) * pow(m_e*k_b*Tb / (2.0*M_PI*hbar*hbar), 1.5) * exp(-epsilon_0/(k_b*Tb));
  
  // Electron fraction and number density
  double Xe = (A < 400.0) ? 0.5 * A * (sqrt(1 + 4.0/A) - 1) : 1;
  double ne = Xe * nb;
  
  return std::pair<double,double>(Xe, ne);
}

//====================================================
// Solve the Saha equation WITH HELIUM to get ne and Xe
//====================================================
std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation_with_He(double x) const{
  const double a           = exp(x);
 
  // Physical constants
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double epsilon_0   = Constants.epsilon_0;
  const double H0_over_h   = Constants.H0_over_h;

  // Fetch cosmological parameters
  const double OmegaB      = cosmo->get_OmegaB();
  const double TCMB        = cosmo->get_TCMB();
  const double H0          = cosmo->get_H0();

  // Baryon temperature
  const double Tb          = TCMB / a;

  // Critical density and number density of baryons
  const double rho_crit    = 3.0 * H0*H0 / (8.0*M_PI*G);
  const double nb          = OmegaB * rho_crit / (m_H * pow(a, 3));

  // Thermal factor (no exponential)
  double A_base = pow(m_e*k_b*Tb / (2.0*M_PI*hbar*hbar), 1.5);

  // Saha RHS factors (already divided by nb)
  double r0 = 2.0 * A_base * exp(-Constants.xhi0 / (k_b*Tb)) / nb;  // He -> He+
  double r1 = 4.0 * A_base * exp(-Constants.xhi1 / (k_b*Tb)) / nb;  // He+ -> He2+
  double r2 =       A_base * exp(-epsilon_0       / (k_b*Tb)) / nb;  // H -> H+

  // Ion fractions as functions of fe
  // He+ = r0/fe / (1 + r0/fe + r0*r1/fe^2)
  // He2+ = (r1/fe) * He+
  // H+  = r2/fe / (1 + r2/fe)  =  r2 / (fe + r2)
  double fe = 1.0;
  for (int iter = 0; iter < 200; iter++) {
      double fe2    = fe * fe;
      double denom_He = fe2 + r0*fe + r0*r1;         // fe^2 * denom of He+
      double He_p   = r0*fe  / denom_He;              // He+  (as fraction of nHe)
      double He_2p  = r0*r1  / denom_He;              // He2+ (as fraction of nHe)
      double H_p    = r2     / (fe + r2);             // H+   (as fraction of nH)

      double fe_new = (2.0*He_2p + He_p) * Yp/4.0 + H_p * (1.0 - Yp);
      
      double diff = std::abs(fe_new - fe);
      fe = fe_new;
      if (diff < 1e-10) break;
  }

  double Xe = fe / (1-Yp);
  double ne = Xe * nb;
  
  return std::pair<double,double>(Xe, ne);
}

//====================================================
// The right hand side of the dXedx Peebles ODE
//====================================================
int RecombinationHistory::rhs_peebles_ode(double x, const double *Xe, double *dXedx){

  // Current value of a and X_e
  const double X_e         = Xe[0];
  const double a           = exp(x);

  // Physical constants in SI units
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double c           = Constants.c;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double sigma_T     = Constants.sigma_T;
  const double lambda_2s1s = Constants.lambda_2s1s;
  const double epsilon_0   = Constants.epsilon_0;

  // Fetch cosmological parameters
  const double OmegaB      = cosmo->get_OmegaB();
  const double TCMB        = cosmo->get_TCMB();
  const double H0          = cosmo->get_H0();
  const double H           = cosmo->H_of_x(x);

  // Baryon temperature
  const double Tb          = TCMB / a;
  
  // Critical density and number density of baryons
  const double rho_crit    = 3.0 * H0*H0 / (8.0*M_PI*G);
  const double nb          = OmegaB * rho_crit / (m_H * pow(a, 3));
  const double nH          = (1-Yp)*nb;

  // Lambda alpha
  const double n1s            = (1-X_e) * nH;
  const double lambda_alpha   = H * pow(3.0*epsilon_0, 3) / (pow(8.0*M_PI, 2)*pow(hbar*c, 3)*n1s);

  // Beta 2
  const double eps0_over_Tb = epsilon_0/(k_b*Tb);
  const double phi2     = 0.448 * log(eps0_over_Tb);
  const double alpha2   = 8.0*sigma_T*c / sqrt(3.0*M_PI) * sqrt(eps0_over_Tb) * phi2;
  const double beta     = alpha2 * pow(m_e*k_b*Tb/(2.0*M_PI*hbar*hbar), 1.5) * exp(-eps0_over_Tb);
  const double beta2    = (eps0_over_Tb > 200) ? 0.0 : beta * exp(3.0*epsilon_0 / (4.0*k_b*Tb));
  
  // Cr constant
  const double Cr       = (lambda_2s1s + lambda_alpha) / (lambda_2s1s + lambda_alpha + beta2);
  
  // The right hand side
  dXedx[0] = (Cr/H) * (beta*(1-X_e) - nH*alpha2*X_e*X_e);

  return GSL_SUCCESS;
}

//====================================================
// Solve for the optical depth tau, compute the 
// visibility function and spline the result
//====================================================

void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("opticaldepth");

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){

    double H  = cosmo->H_of_x(x);
    dtaudx[0] = -Constants.c * ne_of_x(x) * Constants.sigma_T / H;   

    return GSL_SUCCESS;
  };

  ODESolver tau_ode;
  Vector x_array = Utils::linspace(x_end, x_start, npts_rec_arrays);
  Vector tau_ic{0.0};
  tau_ode.solve(dtaudx, x_array, tau_ic);
  auto tau_array = tau_ode.get_data_by_component(0);


  Vector g_tilde_array(npts_rec_arrays);
  Vector dtaudx_array(npts_rec_arrays);
  double H;
  for (int i = 0; i < npts_rec_arrays; i++) {
    H                = cosmo->H_of_x(x_array[i]);
    dtaudx_array[i]  = -Constants.c * ne_of_x(x_array[i]) * Constants.sigma_T / H;  
    g_tilde_array[i] = -dtaudx_array[i] * exp(-tau_array[i]);
  }

  tau_of_x_spline.create(x_array, tau_array, "Spline optical depth tau(x)");
  dtaudx_of_x_spline.create(x_array, dtaudx_array, "Spline dtau/dx");
  g_tilde_of_x_spline.create(x_array, g_tilde_array, "Spline visibility function g_tilde(x)");

  Utils::EndTiming("opticaldepth");
}

void RecombinationHistory::solve_for_sound_horizon(){
  Utils::StartTiming("soundhorizon");

  // The ODE system ds/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dsdx = [&](double x, const double *tau, double *dtaudx){
    double OmegaR     = cosmo -> get_OmegaR();
    double OmegaB     = cosmo -> get_OmegaB();
    double R          = 4*OmegaR / (3*OmegaB) * exp(-x);

    double cs         = Constants.c * sqrt(R / (3 + 3*R));
    double Hp         = cosmo -> Hp_of_x(x);

    dtaudx[0]         = cs / Hp;   

    return GSL_SUCCESS;
  };

  ODESolver sound_horizon_ode;
  Vector x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);
  double OmegaR     = cosmo -> get_OmegaR();
  double OmegaB     = cosmo -> get_OmegaB();
  double R          = 4*OmegaR / (3*OmegaB) * exp(-x_start);

  double cs_init    = Constants.c * sqrt(R / (3 + 3*R));
  double Hp_init    = cosmo -> Hp_of_x(x_start);

  Vector sound_horizon_ic{cs_init / Hp_init};
  sound_horizon_ode.solve(dsdx, x_array, sound_horizon_ic);

  auto sound_horizon_array = sound_horizon_ode.get_data_by_component(0);
  sound_horizon_spline.create(x_array, sound_horizon_array, "Spline sound horizon");
}

//====================================================
// Get methods
//====================================================

double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const{
  return dtaudx_of_x_spline(x);
}

double RecombinationHistory::ddtauddx_of_x(double x) const{
return dtaudx_of_x_spline.deriv_x(x);
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{
  return g_tilde_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{
  return g_tilde_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::Xe_of_x(double x) const{
  return exp(log_Xe_of_x_spline(x));
}

double RecombinationHistory::ne_of_x(double x) const{
  return exp(log_ne_of_x_spline(x));
}

double RecombinationHistory::get_sound_horizon(double x) const{
  return sound_horizon_spline(x);
}

double RecombinationHistory::get_Yp() const{
  return Yp;
}

//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const{
  std::cout << "\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "Yp:             " << Yp                      << "\n";
  std::cout << "Reionisation:   " << bool(reionisation)      << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts       = 2000;
  const double x_min   = x_start;
  const double x_max   = x_end;

  Vector x_array = Utils::linspace(x_min, x_max, npts);
  auto print_data = [&] (const double x) {
    fp << x                    << " ";
    fp << Xe_of_x(x)           << " ";
    fp << ne_of_x(x)           << " ";
    fp << tau_of_x(x)          << " ";
    fp << dtaudx_of_x(x)       << " ";
    fp << ddtauddx_of_x(x)     << " ";
    fp << g_tilde_of_x(x)      << " ";
    fp << dgdx_tilde_of_x(x)   << " ";
    fp << ddgddx_tilde_of_x(x) << " ";
    fp << get_sound_horizon(x) << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

