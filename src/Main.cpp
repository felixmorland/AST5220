#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"
#include "Perturbations.h"
#include "PowerSpectrum.h"
#include "SupernovaFitting.h"
#include <iomanip>
#include <sstream>

int main(int argc, char **argv){
  Utils::StartTiming("Everything");
  //==================================================================
  // Set parameters
  //==================================================================

  // Background parameters
  SimParams.h           = 0.67;
  SimParams.OmegaCDM    = 0.267;
  SimParams.OmegaK      = 0.0;
  SimParams.OmegaB      = 0.05;
  SimParams.Neff        = 3.046;
  SimParams.TCMB        = 2.7255;

  // Include neutrinos and/or polarisation ?
  SimParams.neutrinos     = true;
  SimParams.polarisation  = true;

  // Recombination parameters
  SimParams.Yp          = 0.245;
  SimParams.reionisation  = true;

  // Power-spectrum parameters
  SimParams.A_s         = 2.1e-9;
  SimParams.n_s         = 0.965;
  SimParams.kpivot_mpc  = 0.05;

  // Min and max k-value
  const double Mpc = 3.08567758e22;
  SimParams.k_min = 0.00005 / Mpc;
  SimParams.k_max = 1.0     / Mpc;
  
  // Min and max x-value
  SimParams.x_start = log(1e-8);
  SimParams.x_end   = 0.0;
  
  double h           = SimParams.h;
  double OmegaCDM    = SimParams.OmegaCDM;
  double OmegaK      = SimParams.OmegaK;
  double OmegaB      = SimParams.OmegaB;
  double Neff        = SimParams.Neff;
  double TCMB        = SimParams.TCMB;
  
  double Yp          = SimParams.Yp;
  bool reionisation  = SimParams.reionisation;
  
  bool neutrinos     = SimParams.neutrinos;
  bool polarization  = SimParams.polarisation;
  
  double A_s         = SimParams.A_s;
  double n_s         = SimParams.n_s;
  double kpivot_mpc  = SimParams.kpivot_mpc;

  Utils::print_startup();
  
  //==================================================================
  // Milestone I
  //==================================================================

  // Set up and solve the background
  BackgroundCosmology cosmo(h, OmegaB, OmegaCDM, OmegaK, Neff, TCMB);
  cosmo.solve();
  cosmo.info();
  
  // Output background evolution quantities
  cosmo.output("results/background/fiducial_cosmology.txt");

  // MCMC algorithm 
  //  mcmc_fit_to_supernova_data("data/supernovadata.txt", "results_supernovafitting.txt");

  //==================================================================
  // Milestone II
  //==================================================================
  
  // Solve the recombination history
  RecombinationHistory rec(&cosmo, Yp, reionisation);
  rec.solve();
  rec.info();

  // Output recombination quantities
  rec.output("results/recombination/recombination_He_reion.txt");

  //==================================================================
  // Milestone III
  //==================================================================
 
  // Solve the perturbations
  Perturbations pert(&cosmo, &rec);
  pert.solve();
  pert.info();
  
  // Output perturbation quantities
  Vector kvalues{0.001, 0.01, 0.05, 0.1, 0.2};
  
  for (int i = 0; i < kvalues.size(); i++) {
    double kvalue = kvalues[i] / Constants.Mpc;
    std::string filename = "results/perturbations/perturbations_k" + Utils::format_k(kvalues[i]) + ".txt";
    std::string source_func_filename = "results/perturbations/source_k" + Utils::format_k(kvalues[i]) + ".txt";
    pert.pert_output(kvalue, filename);
    pert.source_func_output(kvalue, source_func_filename);
  }
  
  //==================================================================
  // Milestone IV
  //==================================================================

  PowerSpectrum power(&cosmo, &rec, &pert, A_s, n_s, kpivot_mpc);
  power.solve();
  power.output_CMB_spectrum("results/powerspectrum/cells.txt");
  power.output_matter_power_spectrum("results/powerspectrum/matter_ps.txt");
  power.output_transfer_func("results/powerspectrum/transfer_funcT.txt", power.get_thetaT_ell_of_k_spline());
  power.output_transfer_func("results/powerspectrum/transfer_funcE.txt", power.get_thetaE_ell_of_k_spline());

  std::cout << "\n";
  Utils::EndTiming("Everything");
  std::cout << "\nSimulation complete!";
  std::cout << std::endl;
  return 0;
}
