#include"PowerSpectrum.h"

//====================================================
// Constructors
//====================================================

PowerSpectrum::PowerSpectrum(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec, 
    Perturbations *pert,
    double A_s,
    double n_s,
    double kpivot_mpc) : 
  cosmo(cosmo), 
  rec(rec), 
  pert(pert),
  A_s(A_s),
  n_s(n_s),
  kpivot_mpc(kpivot_mpc)
{}

//====================================================
// Do all the solving
//====================================================
void PowerSpectrum::solve(){

  generate_bessel_function_splines();

  //=========================================================================
  // TODO: Line of sight integration to get Theta_ell(k)
  // Implement line_of_sight_integration
  //=========================================================================
  line_of_sight_integration();

  //=========================================================================
  // TODO: Integration to get Cell by solving dCell^f/dlogk = Delta(k) * f_ell(k)^2
  // Implement solve_for_cell
  //=========================================================================
  auto cell_TT = solve_for_cell(thetaT_ell_of_k_spline, thetaT_ell_of_k_spline);
  cell_TT_spline.create(ells, cell_TT, "Cell_TT_of_ell");
  
  auto cell_TE = solve_for_cell(thetaT_ell_of_k_spline, thetaE_ell_of_k_spline);
  cell_TE_spline.create(ells, cell_TE, "Cell_TE_of_ell");

  auto cell_EE = solve_for_cell(thetaE_ell_of_k_spline, thetaE_ell_of_k_spline);
  cell_EE_spline.create(ells, cell_EE, "Cell_EE_of_ell");
}

//====================================================
// Generate splines of j_ell(z) needed for LOS integration
//====================================================

void PowerSpectrum::generate_bessel_function_splines(){
  Utils::StartTiming("besselspline");
  
  // Make storage for the splines
  j_ell_splines = std::vector<Spline>(ells.size());
    
  //=================================================================
  // Compute splines for bessel functions j_ell(z)
  //=================================================================
  const double z_max   = k_max * cosmo->eta_of_x(0.0);
  const double delta_z = M_PI / 8.0;
  const int n_z        = int(z_max / delta_z);

  Vector z_array = Utils::linspace(0.0, z_max, n_z);

  std::cout << "\nSplining Bessel functions...\n";
  for(size_t i = 0; i < ells.size(); i++){
    Utils::progressbar(double(i) / double(ells.size()));
    const int ell = ells[i];

    // Sample j_ell
    Vector j_ell_sample(n_z);
    for (int iz = 0; iz < n_z; iz++) {
      j_ell_sample[iz] = Utils::j_ell(ell, z_array[iz]);
    }

    // Make the j_ell_splines[i] spline
    j_ell_splines[i].create(z_array, j_ell_sample, "j_" + std::to_string(ell) + "_spline");
  }
  std::cout << "\n";
  Utils::EndTiming("besselspline");
}

//====================================================
// Do the line of sight integration for a single
// source function
//====================================================

Vector2D PowerSpectrum::line_of_sight_integration_single(
    Vector & k_array,
    Vector & x_array,
    std::function<double(double,double)> &source_function,
    std::string & func_name
){
  Utils::StartTiming(func_name);
    
  // Make storage for the results
  Vector2D result = Vector2D(ells.size(), Vector(k_array.size()));

  // Cosmological parameters
  const double eta0 = cosmo->eta_of_x(0.0);

  // Set timestep
  const double dx = x_array[1] - x_array[0];
  
  // Loop over all ell
  // for (int ell = 0; ell < ells.size(); ell++) {
  //   Utils::progressbar(double(ell) / double(ells.size()));
  //   double los_integral = 0.0;

  //   for(int ik = 0; ik < k_array.size(); ik++){
  //     const double k = k_array[ik];


  //     // Trapezoidal integral
  //     for (int ix = 0; ix < n_x; ix++) {

  //       double S_tilde = source_function(x_array[ix], k);
  //       double eta     = cosmo->eta_of_x(x_array[ix]);
  //       double j_ell   = j_ell_splines[ell](k*(eta0 - eta)); 

  //       if (ix == 0 || ix == n_x - 1) {
  //         los_integral += S_tilde*j_ell/2.0;
  //       }
  //       else {
  //         los_integral += S_tilde*j_ell;
  //       }
  //     }
  //     // Store the result for Source_ell(k) in results[ell][ik]
  //     result[ell][ik] = los_integral*dx;
  //   }
  // }
  // Loop over all k
  std::cout << "\nLine of sight integration: " << func_name << "...\n";
  for(int ik = 0; ik < k_array.size(); ik++){
    Utils::progressbar(double(ik) / double(k_array.size()));
    const double k = k_array[ik];

    Vector integral(ells.size(), 0.0);
    for (int ix = 0; ix < n_x; ix++) {
      double S_tilde = source_function(x_array[ix], k);
      double eta     = cosmo->eta_of_x(x_array[ix]);
      for (int ell = 0; ell < ells.size(); ell++) {
        double j_ell   = j_ell_splines[ell](k*(eta0 - eta));
        integral[ell] += S_tilde*j_ell;
      }
    }

    for (int ell = 0; ell < ells.size(); ell++) {
      result[ell][ik] = integral[ell]*dx;
    }
  }
  std::cout << "\n";
  Utils::EndTiming(func_name);
  return result;
}

//====================================================
// Do the line of sight integration
//====================================================
void PowerSpectrum::line_of_sight_integration(){
  const int nells = ells.size();
  
  //============================================================================
  // Make linear spaced k-array with dk ~ 2pi/eta0/N where N >~ 6
  //============================================================================
  const double eta0 = cosmo->eta_of_x(0.0);
  const double dk   = M_PI / (4.0*eta0);
  const int n_k_int = int((k_max - k_min) / dk);

  Vector k_array    = Utils::linspace(k_min, k_max, n_k_int);
  Vector x_array    = Utils::linspace(x_start, x_end, n_x);

  // Make storage for the splines we are to create
  thetaT_ell_of_k_spline      = std::vector<Spline>(nells);

  //============================================================================
  // Solve for Theta_ell(k) and spline the result
  //============================================================================

  // Make a function returning the source function
  std::function<double(double,double)> source_function_T = [&](double x, double k){
    return pert->get_Source_T(x,k);
  };

  // Do the line of sight integration
  std::string func_name = "thetaT";
  Vector2D thetaT_ell_of_k = line_of_sight_integration_single(k_array, x_array, source_function_T, func_name);

  // Spline the result and store it in thetaT_ell_of_k_spline
  for (int i = 0; i < ells.size(); i++) {
    thetaT_ell_of_k_spline[i].create(k_array, thetaT_ell_of_k[i], "thetaT_ell_of_k_spline_ell=" + std::to_string(ells[i]));
  }

  //============================================================================
  // Solve for ThetaE_ell(k) and spline
  //============================================================================
  if(Constants.polarization){

    // Make storage for the splines we are to create
    thetaE_ell_of_k_spline = std::vector<Spline>(nells);

    // Make a function returning the source function
    std::function<double(double,double)> source_function_E = [&](double x, double k){
      return pert->get_Source_E(x,k);
    };

    // Do the line of sight integration
    std::string func_name = "thetaE";
    Vector2D thetaE_ell_of_k = line_of_sight_integration_single(k_array, x_array, source_function_E, func_name);

    // Spline the result and store it in thetaE_ell_of_k_spline
    for (int i = 0; i < ells.size(); i++) {
      thetaE_ell_of_k_spline[i].create(k_array, thetaE_ell_of_k[i], "thetaE_ell_of_k_spline_ell=" + std::to_string(ells[i]));
    }

  }
}

//====================================================
// Compute Cell (could be TT or TE or EE) 
// Cell = Int_0^inf 4 * pi * P(k) f_ell g_ell dk/k
//====================================================
Vector PowerSpectrum::solve_for_cell(
    std::vector<Spline> & f_ell_spline,
    std::vector<Spline> & g_ell_spline
  ){
  const int nells       = ells.size();
  const double eta0     = cosmo->eta_of_x(0.0);
  const double dk       = M_PI / eta0 / 64.0;
  const int n           = int((k_max - k_min) / dk);

  Vector log_k_array    = Utils::linspace(log(k_min), log(k_max), n); 
  double dlog_k         = log_k_array[1] - log_k_array[0];
  Vector result(nells);

  //============================================================================
  // Integrate Cell = Int 4 * pi * P(k) f_ell g_ell dlog_k
  // Trapezoidal rule
  //============================================================================
  for (int i = 0; i < nells; i++) {
    int ell = ells[i];
    double integral = 0.0;

    for (int ik = 0; ik < n; ik++) {
      double k        = exp(log_k_array[ik]);
      double Pk       = primordial_power_spectrum(k);
      double f_ell    = f_ell_spline[i](k);
      double g_ell    = g_ell_spline[i](k);

      if (ik == 0 || ik == n - 1) {
        integral += 2.0 * M_PI * Pk * f_ell * g_ell * dlog_k;
      }
      else {
        integral += 4.0 * M_PI * Pk * f_ell * g_ell * dlog_k;
      }
    }
    result[i] = integral;
  }
  return result;
}

//====================================================
// The primordial power-spectrum
//====================================================

double PowerSpectrum::primordial_power_spectrum(const double k) const{
  return A_s * pow( Constants.Mpc * k / kpivot_mpc , n_s - 1.0);
}

//====================================================
// P(k) in units of (Mpc)^3
//====================================================

double PowerSpectrum::get_matter_power_spectrum(const double x, const double k) const{
  double Hp   = cosmo->Hp_of_x(x);
  double Phi  = pert->get_Phi(x,k);
  double Delta_M = 2.0 / 3.0 * pow(Constants.c*k / Hp, 2.0) * Phi;

  double pofk = 2.0 * pow(M_PI * Delta_M, 2.0) * pow(k * Constants.Mpc, -3.0) * primordial_power_spectrum(k);

  return pofk;
}

//====================================================
// Get methods
//====================================================
double PowerSpectrum::get_cell_TT(const double ell) const{
  return cell_TT_spline(ell);
}
double PowerSpectrum::get_cell_TE(const double ell) const{
  return cell_TE_spline(ell);
}
double PowerSpectrum::get_cell_EE(const double ell) const{
  return cell_EE_spline(ell);
}

//====================================================
// Output the cells to file
//====================================================

void PowerSpectrum::output_CMB_spectrum(std::string filename) const{
  // Output in standard units of muK^2
  std::ofstream fp(filename.c_str());
  const int ellmax = int(ells[ells.size()-1]);
  auto ellvalues = Utils::linspace(2, ellmax, ellmax-1);
  auto print_data = [&] (const double ell) {
    double normfactor  = (ell * (ell+1)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(), 2);
    double normfactorN = (ell * (ell+1)) / (2.0 * M_PI) 
      * pow(1e6 * cosmo->get_TCMB() *  pow(4.0/11.0, 1.0/3.0), 2);
    double normfactorL = (ell * (ell+1)) * (ell * (ell+1)) / (2.0 * M_PI);
    fp << ell                                 << " ";
    fp << cell_TT_spline( ell ) * normfactor  << " ";
    if(Constants.polarization){
      fp << cell_EE_spline( ell ) * normfactor  << " ";
      fp << cell_TE_spline( ell ) * normfactor  << " ";
    }
    fp << "\n";
  };
  std::for_each(ellvalues.begin(), ellvalues.end(), print_data);
}
void PowerSpectrum::output_matter_power_spectrum(std::string filename) const{

  std::ofstream fp(filename.c_str());
  auto k_array = exp(Utils::linspace(log(k_min), log(k_max), 2000));
  auto print_data = [&] (const double k) {
    fp << k * Constants.Mpc                             << " ";
    fp << get_matter_power_spectrum(0.0, k)             << " "; 
    fp << "\n";
  };
  std::for_each(k_array.begin(), k_array.end(), print_data);
}
void PowerSpectrum::output_transfer_func(std::string filename) const{

  std::ofstream fp(filename.c_str());
  fp << "k ";
  for (int i = 0; i < ells.size(); i++){
      int ell = ells[i];
      fp << ell << " ";
  }
  fp << "\n";

  auto k_array = exp(Utils::linspace(log(k_min), log(k_max), 2000));
  auto print_data = [&] (const double k) {
    fp << k * Constants.Mpc                        << " ";
    for (int i = 0; i < ells.size(); i++){
      fp << thetaT_ell_of_k_spline[i](k)           << " ";
    } 
    fp << "\n";
  };
  std::for_each(k_array.begin(), k_array.end(), print_data);
}

