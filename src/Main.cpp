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

  //=========================================================================
  // Parameters
  //=========================================================================

  // Background parameters
  double OmegaB      = 0.05;
  double Neff        = 3.046;
  double TCMB        = 2.7255;

  // Fiducial parameters
  double h           = 0.67;
  double OmegaCDM    = 0.267;
  double OmegaK      = 0.0;

  // Best-fit parameters
  // double h           = 0.701711;
  // double OmegaCDM    = 0.205027;
  // double OmegaK      = 0.0789514;

  // Recombination parameters
  double Yp          = 0.245;
  bool reionisation  = true;

  // Power-spectrum parameters
  double A_s         = 2.1e-9;
  double n_s         = 0.965;
  double kpivot_mpc  = 0.05;
  
  //=========================================================================
  // Module I
  //=========================================================================

  // Set up and solve the background
  BackgroundCosmology cosmo(h, OmegaB, OmegaCDM, OmegaK, Neff, TCMB);
  cosmo.solve();
  cosmo.info();
  
  // Output background evolution quantities
  // cosmo.output("fiducial_cosmology.txt");
  // cosmo.output("best_fit_cosmology.txt");

  //  mcmc_fit_to_supernova_data("data/supernovadata.txt", "results_supernovafitting.txt");

  //=========================================================================
  // Module II
  //=========================================================================
  
  // Solve the recombination history
  RecombinationHistory rec(&cosmo, Yp, reionisation);
  rec.solve();
  rec.info();

  // Output recombination quantities
  rec.output("recombination_data/recombination_He_noreion.txt");

  //=========================================================================
  // Module III
  //=========================================================================
 
  // Solve the perturbations
  Perturbations pert(&cosmo, &rec);
  pert.solve();
  pert.info();
  
  // Output perturbation quantities
  Vector kvalues{0.001, 0.01, 0.05, 0.1, 0.2};
  
  for (int i = 0; i < kvalues.size(); i++) {
    double kvalue = kvalues[i] / Constants.Mpc;
    std::string filename = "perturbation_data/perturbations_k" + Utils::format_k(kvalues[i]) + ".txt";
    std::string source_func_filename = "perturbation_data/source_k" + Utils::format_k(kvalues[i]) + ".txt";
    pert.pert_output(kvalue, filename);
    pert.source_func_output(kvalue, source_func_filename);
  }
  
  //=========================================================================
  // Module IV
  //=========================================================================

  PowerSpectrum power(&cosmo, &rec, &pert, A_s, n_s, kpivot_mpc);
  power.solve();
  power.output_CMB_spectrum("power_spectrum_data/cells.txt");
  power.output_matter_power_spectrum("power_spectrum_data/matter_ps.txt");
  power.output_transfer_func("power_spectrum_data/transfer_func.txt");

  Utils::EndTiming("Everything");
  return 0;
}
