#include"Perturbations.h"

//====================================================
// Constructors
//====================================================

Perturbations::Perturbations(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec) : 
  cosmo(cosmo), 
  rec(rec)
{
  // Neutrino fraction
  double Omega_gamma0 = cosmo->get_OmegaR();
  double Omega_nu0    = cosmo->get_OmegaNu();
  f_nu = Omega_nu0 / (Omega_gamma0 + Omega_nu0) * Constants.neutrinos;
}

//====================================================
// Do all the solving
//====================================================

void Perturbations::solve(){

  // Integrate all the perturbation equation and spline the result
  integrate_perturbations();

  // Compute source functions and spline the result
  // SW, ISW, Doppler, Polarization
  compute_source_functions(true, true, true, true);
}

//====================================================
// The main work: integrate all the perturbations
// and spline the results
//====================================================

void Perturbations::integrate_perturbations(){
  Utils::StartTiming("integrateperturbation");

  //===================================================================
  // TODO: Set up the k-array for the k's we are going to integrate over
  // Start at k_min end at k_max with n_k points with either a
  // quadratic or a logarithmic spacing
  //===================================================================
  Vector k_array = exp(Utils::linspace(log(k_min), log(k_max), n_k));
  Vector x_array = Utils::linspace(x_start, x_end, n_x);

  // 2D for full y-solution to spline afterwards
  int n_ell_tot_full = Constants.n_ell_tot_full;
  Vector2D y_array(n_ell_tot_full);
  for (int i = 0; i < n_ell_tot_full; i++)
    y_array[i] = Vector(n_x * n_k);

  // Physical constants
  double c = Constants.c;

  // Loop over all wavenumbers
  std::cout << "Integrating perturbations...\n";
  for(int ik = 0; ik < n_k; ik++){

    // Update progressbar
    Utils::progressbar(double(ik) / double(n_k));
 
    // Current value of k
    double k = k_array[ik];

    // Find value to integrate to
    auto end_tight     = get_tight_coupling_time(k, x_array);
    double x_end_tight = end_tight.first;
    int idx_end        = end_tight.second;

    //===================================================================
    // Tight coupling integration
    // set_ic : The IC at the start
    // rhs_tight_coupling_ode : The dydx for our coupled ODE system
    //===================================================================

    // Set up initial conditions in the tight coupling regime
    auto y_tight_coupling_ini = set_ic(x_start, k);

    // The tight coupling ODE system
    ODEFunction dydx_tight_coupling = [&](double x, const double *y, double *dydx){
      return rhs_tight_coupling_ode(x, k, y, dydx);
    };

    ODESolver y_tight_ode;
    Vector x_array_tc        = Vector(x_array.begin(), x_array.begin() + idx_end + 1);
    Vector x_array_after_tc  = Vector(x_array.begin() + idx_end, x_array.end());
    y_tight_ode.solve(dydx_tight_coupling, x_array_tc, y_tight_coupling_ini);

    // TC solution array
    int n_ell_tot_tc = Constants.n_ell_tot_tc;
    Vector y_tight_coupling(n_ell_tot_tc);

    //======================================
    // Fill in scalars
    //======================================
    for (int j = 0; j < Constants.n_scalars; j++) {
      auto y_tight_array = y_tight_ode.get_data_by_component(j);
      for (int ix = 0; ix < idx_end; ix++) {
        y_array[j][ix + n_x*ik] = y_tight_array[ix];
      }
      y_tight_coupling[j] = y_tight_array[idx_end];
    }

    //======================================
    // Fill in Theta
    //======================================
    // Theta 0,1
    for (int ell = 0; ell < Constants.n_ell_theta_tc; ell++) {
      auto y_tight_array = y_tight_ode.get_data_by_component(Constants.ind_start_theta_tc + ell);
      for (int ix = 0; ix < idx_end; ix++) {
        y_array[Constants.ind_start_theta + ell][ix + n_x*ik] = y_tight_array[ix];
      }
      y_tight_coupling[Constants.ind_start_theta_tc + ell] = y_tight_array[idx_end];
    }
    for (int ix = 0; ix < idx_end; ix++) {
      double Hp     = cosmo->Hp_of_x(x_array[ix]);
      double dtaudx = rec  ->dtaudx_of_x(x_array[ix]);
      // Theta 2
      y_array[Constants.ind_start_theta + 2][ix + n_x*ik] = Constants.polarization ?
                          -8.0*c*k/(15.0*Hp*dtaudx) * y_array[Constants.ind_start_theta+1][ix + n_x*ik] :
                          -20.0*c*k/(45.0*Hp*dtaudx) * y_array[Constants.ind_start_theta+1][ix + n_x*ik];
      
       // Theta ell > 2
      for (int ell = 3; ell < Constants.n_ell_theta; ell++) {
        y_array[Constants.ind_start_theta + ell][ix + n_x*ik] = -ell/(2.0*ell + 1.0) * (c*k/(Hp*dtaudx)) * y_array[Constants.ind_start_theta + ell - 1][ix + n_x*ik];
      }
    }

    //==========================================
    // Fill in ThetaP
    //==========================================
    if (Constants.polarization) {
      for (int ix = 0; ix < idx_end; ix++) {
        double Hp     = cosmo->Hp_of_x(x_array[ix]);
        double dtaudx = rec  ->dtaudx_of_x(x_array[ix]);
        // ThetaP0
        y_array[Constants.ind_start_thetap][ix + n_x*ik] = 1.25*y_array[Constants.ind_start_theta+2][ix + n_x*ik];
        // ThetaP1
        y_array[Constants.ind_start_thetap+1][ix + n_x*ik] = -c*k/(4.0*Hp*dtaudx)*y_array[Constants.ind_start_theta+2][ix + n_x*ik];
        // ThetaP2
        y_array[Constants.ind_start_thetap+2][ix + n_x*ik] = 0.25*y_array[Constants.ind_start_theta+2][ix + n_x*ik];

        // ell > 2
        for (int ell = 3; ell < Constants.n_ell_thetap; ell++) {
          y_array[Constants.ind_start_thetap + ell][ix + n_x*ik] = -ell/(2.0*ell + 1.0)*(c*k/(Hp*dtaudx))*y_array[Constants.ind_start_thetap + ell - 1][ix + n_x*ik];
        }
      }
    }
    //==========================================
    // Fill in Nu
    //==========================================
    if (Constants.neutrinos) {
      for (int ell = 0; ell < Constants.n_ell_neutrinos_tc; ell++) {
        auto y_tight_array = y_tight_ode.get_data_by_component(Constants.ind_start_nu_tc+ell);
        for (int ix = 0; ix < idx_end; ix++) {
          y_array[Constants.ind_start_nu+ell][ix + n_x*ik] = y_tight_array[ix];
        }
        y_tight_coupling[Constants.ind_start_nu_tc+ell] = y_tight_array[idx_end];
      }
    }

    //====================================================================
    // Full equation integration
    // set_ic_after_tight_coupling : The IC after tight coupling ends
    // rhs_full_ode : The dydx for our coupled ODE system
    //===================================================================

    // Set up initial conditions (y_tight_coupling is the solution at the end of tight coupling)
    auto y_full_ini = set_ic_after_tight_coupling(y_tight_coupling, x_end_tight, k);

    // The full ODE system
    ODEFunction dydx_full = [&](double x, const double *y, double *dydx){
      return rhs_full_ode(x, k, y, dydx);
    };

    ODESolver y_full_ode;
    y_full_ode.solve(dydx_full, x_array_after_tc, y_full_ini);

    for (int i = 0; i < Constants.n_ell_tot_full; i++) {
      auto y_full_array = y_full_ode.get_data_by_component(i);
      for (int ix = idx_end; ix < n_x; ix++) {
        y_array[i][ix + n_x*ik] = y_full_array[ix - idx_end];
      }
    }
  }
  std::cout << "\n";
  Utils::EndTiming("integrateperturbation");
  
  //===========================
  // Making splines
  //===========================
  delta_cdm_spline.create(x_array, k_array, y_array[Constants.ind_deltacdm]);
  delta_b_spline.create(x_array, k_array, y_array[Constants.ind_deltab]);
  v_cdm_spline.create(x_array, k_array, y_array[Constants.ind_vcdm]);
  v_b_spline.create(x_array, k_array, y_array[Constants.ind_vb]);
  Phi_spline.create(x_array, k_array, y_array[Constants.ind_Phi]);

  Theta_spline = std::vector<Spline2D>(Constants.n_ell_theta);
  for (int ell = 0; ell < Constants.n_ell_theta; ell++) {
    Theta_spline[ell].create(x_array, k_array, y_array[Constants.ind_start_theta+ell]);
  } 
  if (Constants.polarization) {
    ThetaP_spline = std::vector<Spline2D>(Constants.n_ell_thetap);
    for (int ell = 0; ell < Constants.n_ell_thetap; ell++) {
      ThetaP_spline[ell].create(x_array, k_array, y_array[Constants.ind_start_thetap+ell]);
    }
  }
  if (Constants.neutrinos) {
    Nu_spline = std::vector<Spline2D>(Constants.n_ell_neutrinos);
    for (int ell = 0; ell < Constants.n_ell_neutrinos; ell++) {
      Nu_spline[ell].create(x_array, k_array, y_array[Constants.ind_start_nu+ell]);
    }
  }

  Vector Psi_array(n_x*n_k);
  Vector Pi_array(n_x*n_k);

  // Constants
  const double H0       = cosmo->get_H0();
  const double OmegaR0  = cosmo->get_OmegaR();
  const double OmegaNu0 = cosmo->get_OmegaNu();

  for (int ix = 0; ix < n_x; ix++) {
    double x = x_array[ix];

    for (int ik = 0; ik < n_k; ik++) {
      double k = k_array[ik];

      Psi_array[ix + n_x*ik] = Constants.neutrinos ?
                    -get_Phi(x,k) - 12.0*pow(H0/(c*k), 2.0)*exp(-2.0*x) * (OmegaR0*get_Theta(x,k,2) + OmegaNu0*get_Nu(x,k,2)) :
                    -get_Phi(x,k) - 12.0*pow(H0/(c*k), 2.0)*exp(-2.0*x) * (OmegaR0*get_Theta(x,k,2));
      
      Pi_array[ix + n_x*ik] = Constants.polarization ?
                    get_Theta(x,k,2) + get_Theta_p(x,k,0) + get_Theta_p(x,k,2) :
                    get_Theta(x,k,2);
    }
  }
  Psi_spline.create(x_array, k_array, Psi_array, "Psi_spline");
  Pi_spline.create(x_array, k_array, Pi_array, "Pi_spline");
}

//====================================================
// Set IC at the start of the run (this is in the
// tight coupling regime)
//====================================================
Vector Perturbations::set_ic(const double x, const double k) const{

  // The vector we are going to fill
  Vector y_tc(Constants.n_ell_tot_tc);

  //=============================================================================
  // Compute where in the y_tc array each component belongs
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const int n_ell_tot_tc        = Constants.n_ell_tot_tc;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // References to the tight coupling quantities
  double &delta_cdm    =  y_tc[Constants.ind_deltacdm_tc];
  double &delta_b      =  y_tc[Constants.ind_deltab_tc];
  double &v_cdm        =  y_tc[Constants.ind_vcdm_tc];
  double &v_b          =  y_tc[Constants.ind_vb_tc];
  double &Phi          =  y_tc[Constants.ind_Phi_tc];
  double *Theta        = &y_tc[Constants.ind_start_theta_tc];
  double *Nu           = neutrinos ? &y_tc[Constants.ind_start_nu_tc] : nullptr;

  // Physical constants
  double c = Constants.c;

  // Background quantities
  double Hp       = cosmo->Hp_of_x(x);
  double H0       = cosmo->get_H0();
  double OmegaNu0 = cosmo->get_OmegaNu();

  // Scalar quantities (Gravitational potential, baryons and CDM)
  double Psi  = -1.0 / (3.0/2.0 + 2.0*f_nu/5.0);
  Phi         = -(1.0 + 2.0*f_nu/5.0) * Psi;

  delta_cdm   = -1.5 * Psi;
  delta_b     = delta_cdm;
  v_cdm       = -c*k/(2.0*Hp) * Psi;
  v_b         = v_cdm;

  // Photon temperature perturbations (monopole and dipole)
  *Theta         = - 1.0/2.0 * Psi;
  *(Theta+1)     = c*k/(6.0*Hp) * Psi;

  // Neutrino perturbations
  if(neutrinos){
    *Nu     = *Theta;
    *(Nu+1) = *(Theta+1);
    *(Nu+2) = - pow(c*k/H0, 2.0) * exp(2.0*x) * (Phi + Psi) / (12.0*OmegaNu0);
    for (int ell = 3; ell < n_ell_neutrinos_tc; ell++) {
      *(Nu+ell) = *(Nu+ell-1) * (c*k / ((2.0*ell + 1.0)*Hp));
    }
  }

  return y_tc;
}

//====================================================
// Set IC for the full ODE system after tight coupling 
// regime ends
//====================================================

Vector Perturbations::set_ic_after_tight_coupling(
    const Vector &y_tc, 
    const double x, 
    const double k) const{

  // Make the vector we are going to fill
  Vector y(Constants.n_ell_tot_full);
  
  //=============================================================================
  // Compute where in the y array each component belongs and where corresponding
  // components are located in the y_tc array
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Number of multipoles we have in the full regime
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // Number of multipoles we have in the tight coupling regime
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;

  // References to the tight coupling quantities
  const double &delta_cdm_tc    =  y_tc[Constants.ind_deltacdm_tc];
  const double &delta_b_tc      =  y_tc[Constants.ind_deltab_tc];
  const double &v_cdm_tc        =  y_tc[Constants.ind_vcdm_tc];
  const double &v_b_tc          =  y_tc[Constants.ind_vb_tc];
  const double &Phi_tc          =  y_tc[Constants.ind_Phi_tc];
  const double *Theta_tc        = &y_tc[Constants.ind_start_theta_tc];
  const double *Nu_tc           = neutrinos ? &y_tc[Constants.ind_start_nu_tc] : nullptr;

  // References to the quantities we are going to set
  double &delta_cdm       =  y[Constants.ind_deltacdm];
  double &delta_b         =  y[Constants.ind_deltab];
  double &v_cdm           =  y[Constants.ind_vcdm];
  double &v_b             =  y[Constants.ind_vb];
  double &Phi             =  y[Constants.ind_Phi];
  double *Theta           = &y[Constants.ind_start_theta];
  double *ThetaP          = polarization ? &y[Constants.ind_start_thetap] : nullptr;
  double *Nu              = neutrinos ? &y[Constants.ind_start_nu] : nullptr;

  // Physical constants
  double c = Constants.c;

  // Background quantities
  double Hp       = cosmo->Hp_of_x(x);
  double dtaudx   = rec  ->dtaudx_of_x(x);
  double H0       = cosmo->get_H0();
  double OmegaNu0 = cosmo->get_OmegaNu();

  // Scalar quantities (Gravitational potental, baryons and CDM)
  // These are taken from the TC solution
  Phi       = Phi_tc;
  delta_cdm = delta_cdm_tc;
  delta_b   = delta_b_tc;
  v_cdm     = v_cdm_tc;
  v_b       = v_b_tc;

  // Photon temperature perturbations (Theta_ell)
  // Monopole and dipole taken from TC solution
  *Theta     = *Theta_tc;
  *(Theta+1) = *(Theta_tc+1);
  *(Theta+2) = polarization ? 
              -8.0*c*k / (15.0*Hp*dtaudx) * *(Theta+1) : 
              -20.0*c*k / (45.0*Hp*dtaudx) * *(Theta+1);

  for (int ell = 3; ell < n_ell_theta; ell++) {
      *(Theta+ell) = -*(Theta+ell-1) * (ell*c*k / ((2.0*ell + 1.0)*Hp*dtaudx));
  }

  // Photon polarization perturbations (Theta_p_ell)
  if(polarization){
    *ThetaP     = 1.25 * *(Theta+2);
    *(ThetaP+1) = -c*k / (4*Hp*dtaudx) * *(Theta+2);
    *(ThetaP+2) = 0.25 * *(Theta+2);

    for (int ell = 3; ell < n_ell_thetap; ell++) {
      *(ThetaP+ell) = -*(ThetaP+ell-1) * (ell*c*k / ((2.0*ell + 1.0)*Hp*dtaudx));
    }
  }

  // Neutrino perturbations (N_ell)
  // Same as in TC regime
  if(neutrinos){
    for (int ell = 0; ell < n_ell_neutrinos; ell++) {
      *(Nu+ell) = *(Nu_tc+ell);
    }
  }

  return y;
}

//====================================================
// The time when tight coupling end
//====================================================

std::pair<double,int> Perturbations::get_tight_coupling_time(const double k, Vector x_array) const{
  double x_tight_coupling_end = 0.0;
  int idx_end = 0;

  for (int i = 0; i < n_x; i++) {
    double dtaudx = rec   ->  dtaudx_of_x(x_array[i]);
    double Hp     = cosmo ->  Hp_of_x(x_array[i]);
    // Check if tight coupling criterion breaks
    if (x_array[i] > -8.3 || abs(dtaudx) < 10*std::max(1.0, Constants.c*k/Hp)) {
      x_tight_coupling_end = x_array[i];
      break;
    }
    idx_end += 1;
  }
  return std::pair<double,int>(x_tight_coupling_end, idx_end);
}

//====================================================
// After integrsating the perturbation compute the
// source function(s)
//====================================================
void Perturbations::compute_source_functions(
  bool SW_on,
  bool ISW_on,
  bool Dop_on,
  bool Pol_on
){
  Utils::StartTiming("source");

  Vector k_array = exp(Utils::linspace(log(k_min), log(k_max), n_k));
  Vector x_array = Utils::linspace(x_start, x_end, n_x);

  // Make storage for the source functions (in 1D array to be able to pass it to the spline)
  Vector ST_array(k_array.size() * x_array.size());   // Temperature
  Vector SE_array(k_array.size() * x_array.size());   // Polarization

  for(auto ix = 0; ix < x_array.size(); ix++){
    const double x = x_array[ix];

    // From background
    double tau          = rec->tau_of_x(x);
    double g_tilde      = rec->g_tilde_of_x(x);
    double dgdx_tilde   = rec->dgdx_tilde_of_x(x);
    double ddgddx_tilde = rec->ddgddx_tilde_of_x(x);

    double Hp           = cosmo->Hp_of_x(x);
    double dHpdx        = cosmo->dHpdx_of_x(x);
    double ddHpddx      = cosmo->ddHpddx_of_x(x);
    double chi          = cosmo->get_comoving_distance_of_x(x);

    for(auto ik = 0; ik < k_array.size(); ik++){
      const double k = k_array[ik];
      const int index = ix + n_x * ik;

      // Perturbations
      double v_b          = get_v_b(x,k);
      double dv_bdx       = get_dv_bdx(x,k);
      double Theta0       = get_Theta(x,k,0);
      double Phi          = get_Phi(x,k);
      double dPhidx       = get_dPhidx(x,k);
      double Psi          = get_Psi(x,k);
      double dPsidx       = get_dPsidx(x,k);
      double ddPsiddx     = get_ddPsiddx(x,k);
      double Pi           = get_Pi(x,k);
      double dPidx        = get_dPidx(x,k);
      double ddPiddx      = get_ddPiddx(x,k);

      // The terms in the temperature source function
      double SW_term            = SW_on ?
                                  g_tilde * (Theta0 + Psi + Pi / 4.0) :
                                  0.0;

      double ISW_term           = ISW_on ?
                                  (dPsidx - dPhidx) * exp(-tau) :
                                  0.0;

      double Doppler_term       = Dop_on ?
                                  - (dHpdx*g_tilde*v_b + Hp*dgdx_tilde*v_b 
                                  + Hp*g_tilde*dv_bdx) / (Constants.c * k) :
                                  0.0;
      
      double Polarization_term  = Pol_on ?
                                  (3.0 / (4.0*pow(Constants.c*k, 2.0)))
                            * (dHpdx * (dHpdx*g_tilde*Pi + Hp*dgdx_tilde*Pi + Hp*g_tilde*dPidx)
                            + Hp * (ddHpddx*g_tilde*Pi + Hp*ddgddx_tilde*Pi + Hp*g_tilde*ddPiddx)
                            + 2.0*Hp * (dHpdx*dgdx_tilde*Pi + dHpdx*g_tilde*dPidx + Hp*dgdx_tilde*dPidx)) :
                                  0.0;
      
      ST_array[index] = SW_term + ISW_term + Doppler_term + Polarization_term;

      if(Constants.polarization){
        SE_array[index] = (x <= -3.0) ?
                          3.0*g_tilde*Pi / (4.0*Hp*pow(k*chi, 2.0)) :
                          SE_array[index - 1];
      }
    }
  }

  // Spline the source functions
  ST_spline.create(x_array, k_array, ST_array, "Source_Temp_x_k");
  if(Constants.polarization){
    SE_spline.create(x_array, k_array, SE_array, "Source_Pol_x_k");
  }

  Utils::EndTiming("source");
}

//====================================================
// The right hand side of the perturbations ODE
// in the tight coupling regime
//====================================================

// Derivatives in the tight coupling regime
int Perturbations::rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx){

  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const bool neutrinos          = Constants.neutrinos;
  const bool polarization       = Constants.polarization;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  const double &delta_b         =  y[Constants.ind_deltab_tc];
  const double &v_cdm           =  y[Constants.ind_vcdm_tc];
  const double &v_b             =  y[Constants.ind_vb_tc];
  const double &Phi             =  y[Constants.ind_Phi_tc];
  const double *Theta           = &y[Constants.ind_start_theta_tc];
  const double *Nu              = neutrinos ? &y[Constants.ind_start_nu_tc] : nullptr;

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm_tc];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab_tc];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm_tc];
  double &dv_bdx          =  dydx[Constants.ind_vb_tc];
  double &dPhidx          =  dydx[Constants.ind_Phi_tc];
  double *dThetadx        = &dydx[Constants.ind_start_theta_tc];
  double *dNudx           = neutrinos ? &dydx[Constants.ind_start_nu_tc] : nullptr;

  // Physical constants
  double c = Constants.c;

  // Background quantities
  double Hp       = cosmo->Hp_of_x(x);
  double dHp_dx   = cosmo->dHpdx_of_x(x);
  double dtaudx   = rec  ->dtaudx_of_x(x);
  double ddtauddx = rec  ->ddtauddx_of_x(x);
  double H0       = cosmo->get_H0();
  double OmegaNu0 = cosmo->get_OmegaNu();
  double OmegaB0  = cosmo->get_OmegaB();
  double OmegaR0  = cosmo->get_OmegaR();
  double OmegaCDM0 = cosmo->get_OmegaCDM();
  double R        = 4.0*OmegaR0 / (3.0*OmegaB0) * exp(-x);
  double eta      = cosmo->eta_of_x(x);

  // Scalar quantities (Phi, delta, v, ...)
  double Theta2 = polarization ?
                  -8.0*c*k / (15.0*Hp*dtaudx) * *(Theta+1) : 
                  -20.0*c*k / (45.0*Hp*dtaudx) * *(Theta+1);

  double Psi = neutrinos ?
              -Phi - 12.0*pow(H0/(c*k), 2.0)*exp(-2.0*x) * (OmegaR0*Theta2 + OmegaNu0* *(Nu+2)) :
              -Phi - 12.0*pow(H0/(c*k), 2.0)*exp(-2.0*x) * (OmegaR0*Theta2);

  dPhidx = neutrinos ?
          // With neutrinos
          Psi - pow(c*k/Hp, 2.0)*Phi/3.0 + pow(H0/Hp, 2.0)/2.0 * (OmegaCDM0*exp(-x)*delta_cdm + OmegaB0*exp(-x)*delta_b + 4.0*OmegaR0*exp(-2.0*x)* *Theta + 4.0*OmegaNu0*exp(-2.0*x)* *Nu) :
          // Without neutrinos
          Psi - pow(c*k/Hp, 2.0)*Phi/3.0 + pow(H0/Hp, 2.0)/2.0 * (OmegaCDM0*exp(-x)*delta_cdm + OmegaB0*exp(-x)*delta_b + 4.0*OmegaR0*exp(-2.0*x)* *Theta);

  
  
  ddelta_cdmdx = c*k/Hp * v_cdm - 3.0*dPhidx;
  dv_cdmdx     = -v_cdm - c*k/Hp * Psi;
  ddelta_bdx   = c*k/Hp * v_b - 3.0*dPhidx;

  // Photon monopole
  *dThetadx = -c*k/Hp * *(Theta+1) - dPhidx;

  double q      = (-((1.0-R)*dtaudx + (1.0+R)*ddtauddx)*(3.0* *(Theta+1) + v_b) - (c*k/Hp)*Psi 
                  + (1.0 - dHp_dx/Hp)*c*k/Hp*(-*Theta + 2.0*Theta2) - c*k/Hp* *dThetadx) / ((1.0+R)*dtaudx + dHp_dx/Hp - 1.0);
  
  // Baryon perturbation velocity
  dv_bdx        = (-v_b - c*k/Hp * Psi + R*(q + c*k/Hp*(-*Theta + 2.0*Theta2) 
                  - (c*k/Hp)*Psi)) / (1.0+R);
  // Photon dipole
  *(dThetadx+1) = (q - dv_bdx) / 3.0;

  // Neutrino mutlipoles (Nu_ell)
  if(neutrinos){
    // Monopole and dipole
    *dNudx = -(c*k/Hp) * *(Nu+1) - dPhidx;
    *(dNudx+1) = c*k/(3.0*Hp) * (*Nu - 2.0* *(Nu+2) + Psi);

    // 2 <= ell < ell_max
    for (int ell = 2; ell < n_ell_neutrinos_tc-1; ell++) {
      *(dNudx+ell) = c*k / ((2.0*ell + 1.0)*Hp) * (ell* *(Nu+ell-1) - (ell+1.0)* *(Nu+ell+1));
    }

    // ell_max
    *(dNudx+n_ell_neutrinos_tc-1) = (c*k/Hp)* (*(Nu+n_ell_neutrinos_tc-2)
                          - (n_ell_neutrinos_tc)/(k*eta)* *(Nu+n_ell_neutrinos_tc-1));
  }

  return GSL_SUCCESS;
}

//====================================================
// The right hand side of the full ODE
//====================================================

int Perturbations::rhs_full_ode(double x, double k, const double *y, double *dydx){
  
  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Index and number of the different quantities
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm];
  const double &delta_b         =  y[Constants.ind_deltab];
  const double &v_cdm           =  y[Constants.ind_vcdm];
  const double &v_b             =  y[Constants.ind_vb];
  const double &Phi             =  y[Constants.ind_Phi];
  const double *Theta           = &y[Constants.ind_start_theta];
  const double *ThetaP          = polarization ? &y[Constants.ind_start_thetap] : nullptr;
  const double *Nu              = neutrinos ? &y[Constants.ind_start_nu] : nullptr;

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm];
  double &dv_bdx          =  dydx[Constants.ind_vb];
  double &dPhidx          =  dydx[Constants.ind_Phi];
  double *dThetadx        = &dydx[Constants.ind_start_theta];
  double *dThetaPdx       = polarization ? &dydx[Constants.ind_start_thetap] : nullptr;
  double *dNudx           = neutrinos ? &dydx[Constants.ind_start_nu] : nullptr;

  // Physical constants
  double c = Constants.c;

  // Background quantities
  double Hp       = cosmo->Hp_of_x(x);
  double dHp_dx   = cosmo->dHpdx_of_x(x);
  double dtaudx   = rec  ->dtaudx_of_x(x);
  double ddtauddx = rec  ->ddtauddx_of_x(x);
  double H0       = cosmo->get_H0();
  double OmegaNu0 = cosmo->get_OmegaNu();
  double OmegaB0  = cosmo->get_OmegaB();
  double OmegaR0  = cosmo->get_OmegaR();
  double OmegaCDM0 = cosmo->get_OmegaCDM();
  double R        = 4.0*OmegaR0 / (3.0*OmegaB0) * exp(-x);
  double eta      = cosmo->eta_of_x(x);

  //===========================================
  // Scalar quantities (Phi, delta, v, ...)
  //===========================================
  double Psi = neutrinos ?
              -Phi - 12.0*pow(H0/(c*k), 2.0)*exp(-2.0*x) * (OmegaR0* *(Theta+2) + OmegaNu0* *(Nu+2)) :
              -Phi - 12.0*pow(H0/(c*k), 2.0)*exp(-2.0*x) * (OmegaR0* *(Theta+2));
  
  dPhidx = neutrinos ?
          // With neutrinos
          Psi - pow(c*k/Hp, 2.0)*Phi/3.0 + pow(H0/Hp, 2.0)/2.0 * (OmegaCDM0*exp(-x)*delta_cdm + OmegaB0*exp(-x)*delta_b + 4.0*OmegaR0*exp(-2.0*x)* *Theta + 4.0*OmegaNu0*exp(-2.0*x)* *Nu) :
          // Without neutrinos
          Psi - pow(c*k/Hp, 2.0)*Phi/3.0 + pow(H0/Hp, 2.0)/2.0 * (OmegaCDM0*exp(-x)*delta_cdm + OmegaB0*exp(-x)*delta_b + 4.0*OmegaR0*exp(-2.0*x)* *Theta);
  
  // Scalar perturbations
  ddelta_cdmdx = (c*k/Hp)*v_cdm - 3.0*dPhidx;
  dv_cdmdx     = -v_cdm - (c*k/Hp)*Psi;
  ddelta_bdx   = (c*k/Hp)*v_b - 3.0*dPhidx;
  dv_bdx       = -v_b - (c*k/Hp)*Psi + dtaudx*R*(*(Theta+1) + v_b/3.0);

  //====================================
  // Photon multipoles (Theta_ell)
  //====================================
  // Monopole an dipole
  *dThetadx     = -(c*k/Hp)* *(Theta+1) - dPhidx;
  *(dThetadx+1) = (c*k/(3.0*Hp)) * (*Theta - 2.0* *(Theta+2) + Psi) + dtaudx*(*(Theta+1) + v_b/3.0);

  double Pi = polarization ?
              *(Theta+2) + *(ThetaP) + *(ThetaP+2) :
              *(Theta+2);

  // 2 <= ell < ell_max
  for (int ell = 2; ell < n_ell_theta-1; ell++) {
    *(dThetadx+ell) = polarization ? 
                      (c*k/((2.0*ell + 1.0)*Hp)) * (ell* *(Theta+ell-1) - (ell+1)* *(Theta+ell+1)) + dtaudx*(*(Theta+ell) - Pi/10.0*(ell==2)) :
                      (c*k/((2.0*ell + 1.0)*Hp)) * (ell* *(Theta+ell-1) - (ell+1)* *(Theta+ell+1)) + dtaudx*(*(Theta+ell));
  }

  // ell_max
  *(dThetadx+n_ell_theta-1) = (c*k/Hp) * (*(Theta+n_ell_theta-2) - (n_ell_theta/(k*eta))* *(Theta+n_ell_theta-1)) + dtaudx* *(Theta+n_ell_theta-1);

  //=================================================
  // Photon polarization multipoles (Theta_p_ell)
  //=================================================
  if(polarization){
    // Monopole
    *(dThetaPdx) = -(c*k/Hp)* *(ThetaP+1) + dtaudx*(*ThetaP - Pi/2.0);

    // 1 <= ell < ell_max
    for (int ell = 1; ell < n_ell_thetap-1; ell++) {
      *(dThetaPdx+ell) = (c*k/((2.0*ell + 1.0)*Hp)) * (ell* *(ThetaP+ell-1) - (ell+1)* *(ThetaP+ell+1)) + dtaudx*(*(ThetaP+ell) - Pi/10.0*(ell==2));
    }

    // ell_max
    *(dThetaPdx+n_ell_thetap-1) = (c*k/Hp) * (*(ThetaP+n_ell_thetap-2) - (n_ell_thetap/(k*eta))* *(ThetaP+n_ell_thetap-1)) + dtaudx* *(ThetaP+n_ell_thetap-1);
  }

  //================================
  // Neutrino mutlipoles (Nu_ell)
  //================================
  if(neutrinos){
    // Monopole and dipole
    *dNudx = -(c*k/Hp) * *(Nu+1) - dPhidx;
    *(dNudx+1) = c*k/(3.0*Hp) * (*Nu - 2.0* *(Nu+2) + Psi);

    // 2 <= ell < ell_max
    for (int ell = 2; ell < n_ell_neutrinos-1; ell++) {
      *(dNudx+ell) = c*k / ((2.0*ell + 1.0)*Hp) * (ell* *(Nu+ell-1) - (ell+1.0)* *(Nu+ell+1));
    }

    // ell_max
    *(dNudx+n_ell_neutrinos-1) = (c*k/Hp)* (*(Nu+n_ell_neutrinos-2)
                          - (n_ell_neutrinos)/(k*eta)* *(Nu+n_ell_neutrinos-1));
  }

  return GSL_SUCCESS;
}

//====================================================
// Get methods
//====================================================

double Perturbations::get_delta_cdm(const double x, const double k) const{
  return delta_cdm_spline(x,k);
}
double Perturbations::get_delta_b(const double x, const double k) const{
  return delta_b_spline(x,k);
}
double Perturbations::get_v_cdm(const double x, const double k) const{
  return v_cdm_spline(x,k);
}
double Perturbations::get_v_b(const double x, const double k) const{
  return v_b_spline(x,k);
}
double Perturbations::get_dv_bdx(const double x, const double k) const{
  return v_b_spline.deriv_x(x,k);
}
double Perturbations::get_Phi(const double x, const double k) const{
  return Phi_spline(x,k);
}
double Perturbations::get_dPhidx(const double x, const double k) const{
  return Phi_spline.deriv_x(x,k);
}
double Perturbations::get_Psi(const double x, const double k) const{
  return Psi_spline(x,k);
}
double Perturbations::get_dPsidx(const double x, const double k) const{
  return Psi_spline.deriv_x(x,k);
}
double Perturbations::get_ddPsiddx(const double x, const double k) const{
  return Psi_spline.deriv_xx(x,k);
}
double Perturbations::get_Pi(const double x, const double k) const{
  return Pi_spline(x,k);
}
double Perturbations::get_dPidx(const double x, const double k) const{
  return Pi_spline.deriv_x(x,k);
}
double Perturbations::get_ddPiddx(const double x, const double k) const{
  return Pi_spline.deriv_xx(x,k);
}
double Perturbations::get_Source_T(const double x, const double k) const{
  return ST_spline(x,k);
}
double Perturbations::get_Source_E(const double x, const double k) const{
  return SE_spline(x,k);
}
double Perturbations::get_Theta(const double x, const double k, const int ell) const{
  return Theta_spline[ell](x,k);
}
double Perturbations::get_Theta_p(const double x, const double k, const int ell) const{
  return ThetaP_spline[ell](x,k);
}
double Perturbations::get_Nu(const double x, const double k, const int ell) const{
  return Nu_spline[ell](x,k);
}

//====================================================
// Print some useful info about the class
//====================================================

void Perturbations::info() const{
  std::cout << "\n";
  std::cout << "Info about perturbations class:\n";
  std::cout << "x_start:       " << x_start                << "\n";
  std::cout << "x_end:         " << x_end                  << "\n";
  std::cout << "n_x:     " << n_x              << "\n";
  std::cout << "k_min (1/Mpc): " << k_min * Constants.Mpc  << "\n";
  std::cout << "k_max (1/Mpc): " << k_max * Constants.Mpc  << "\n";
  std::cout << "n_k:     " << n_k              << "\n";
  if(Constants.polarization)
    std::cout << "We include polarization\n";
  else
    std::cout << "We do not include polarization\n";
  if(Constants.neutrinos)
    std::cout << "We include neutrinos\n";
  else
    std::cout << "We do not include neutrinos\n";

  std::cout << "Information about the perturbation system:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm         << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab           << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm             << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb               << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi              << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta      << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta          << "\n";
  if(Constants.polarization){
    std::cout << "ind_start_thetap:   " << Constants.ind_start_thetap   << "\n";
    std::cout << "n_ell_thetap:       " << Constants.n_ell_thetap       << "\n";
  }
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu       << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos    << "\n";
  }
  std::cout << "n_ell_tot_full:     " << Constants.n_ell_tot_full       << "\n";

  std::cout << "Information about the perturbation system in tight coupling:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm_tc      << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab_tc        << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm_tc          << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb_tc            << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi_tc           << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta_tc   << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta_tc       << "\n";
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu_tc    << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos_tc << "\n";
  }
  std::cout << "n_ell_tot_tc:       " << Constants.n_ell_tot_tc         << "\n";
  std::cout << std::endl;
}

//====================================================
// Output some results to file for a given value of k
//====================================================

void Perturbations::pert_output(const double k, const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts = 5000;
  auto x_array = Utils::linspace(x_start, x_end, npts);
  auto print_data = [&] (const double x) {
    fp << x                     << " ";
    fp << get_Theta(x,k,0)      << " ";
    fp << get_Theta(x,k,1)      << " ";
    fp << get_Theta(x,k,2)      << " ";
    fp << get_Theta_p(x,k,0)    << " ";
    fp << get_Theta_p(x,k,1)    << " ";
    fp << get_Theta_p(x,k,2)    << " ";
    fp << get_Nu(x,k,0)         << " ";
    fp << get_Nu(x,k,1)         << " ";
    fp << get_Nu(x,k,2)         << " ";
    fp << get_Phi(x,k)          << " ";
    fp << get_Psi(x,k)          << " ";
    fp << get_Pi(x,k)           << " ";
    fp << get_delta_b(x,k)      << " ";
    fp << get_delta_cdm(x,k)    << " ";
    fp << get_v_b(x,k)          << " ";
    fp << get_v_cdm(x,k)        << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

void Perturbations::source_func_output(const double k, const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts = 5000;
  auto x_array = Utils::linspace(x_start, x_end, npts);
  auto print_data = [&] (const double x) {
    double arg = k * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));
    fp << x                     << " ";
    fp << get_Source_T(x,k)  << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(5,   arg)           << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(50,  arg)           << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(500, arg)           << " ";
    fp << get_Source_E(x,k)  << " ";
    fp << get_Source_E(x,k) * Utils::j_ell(5,   arg)           << " ";
    fp << get_Source_E(x,k) * Utils::j_ell(50,  arg)           << " ";
    fp << get_Source_E(x,k) * Utils::j_ell(500, arg)           << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}