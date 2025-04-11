/*
 * File: TWorldParameters.cpp
 *
 * Institute of Biomedical Engineering, 
 * Karlsruhe Institute of Technology (KIT)
 * https://www.ibt.kit.edu
 * 
 * Repository: https://github.com/KIT-IBT/CardioMechanics
 *
 * License: GPL-3.0 (See accompanying file LICENSE or visit https://www.gnu.org/licenses/gpl-3.0.html)
 *
 */


#include <TWorldParameters.h>
TWorldParameters::TWorldParameters(const char *initFile, ML_CalcType tinc) {
  P = new Parameter[vtLast];
  Init(initFile, tinc);
}

TWorldParameters::~TWorldParameters() {}

void TWorldParameters::PrintParameters() {
  cout << "TWorldParameters:"<<endl;
  for (int i = vtFirst; i < vtLast; i++) {
    cout << "\t" << P[i].name << "\t= "     << P[i].value <<endl;
  }
}

void TWorldParameters::Init(const char *initFile, ML_CalcType tinc) {
#if KADEBUG
  cerr << "Loading the TWorld parameter from " << initFile << " ...\n";
#endif  // if KADEBUG

//  P[VT_amplitude].name        = "amplitude";
//  P[VT_duration].name         = "duration";
//  P[VT_celltype].name         = "celltype";
//  P[VT_Na_o].name             = "Na_o";
//  P[VT_Ca_o].name             = "Ca_o";
//  P[VT_K_o].name              = "K_o";
//  P[VT_Cl_o].name             = "Cl_o";
//    
//  P[VT_INaF_Multiplier].name  = "INaF_Multiplier";
//  P[VT_INaL_Multiplier].name  = "INaL_Multiplier";
//  P[VT_Ito_Multiplier].name   = "Ito_Multiplier";
//  P[VT_PCa_Multiplier].name   = "PCa_Multiplier";
//  P[VT_IKr_Multiplier].name   = "IKr_Multiplier";
//  P[VT_IKs_Multiplier].name   = "IKs_Multiplier";
//  P[VT_IK1_Multiplier].name   = "IK1_Multiplier";
//  P[VT_INaCa_Multiplier].name = "INaCa_Multiplier";
//  P[VT_INaK_Multiplier].name  = "INaK_Multiplier";
//  P[VT_IKb_Multiplier].name   = "IKb_Multiplier";
//  P[VT_INab_Multiplier].name  = "INab_Multiplier";
//  P[VT_ICab_Multiplier].name  = "ICab_Multiplier";
//  P[VT_IpCa_Multiplier].name  = "IpCa_Multiplier";
//  P[VT_Ito_slow_Multiplier].name = "Ito_slow_Multiplier";
//  P[VT_Ito_fast_Multiplier].name = "Ito_fast_Multiplier";
//  P[VT_ICaL_fractionSS].name = "ICaL_fractionSS";
//  P[VT_CMDN_Multiplier].name = "CMDN_Multiplier";
//  P[VT_INaCa_fractionSS].name = "INaCa_fractionSS";
//  P[VT_Jup_Multiplier].name = "Jup_Multiplier";
//  P[VT_ICaCl_Multiplier].name = "ICaCl_Multiplier";
//  P[VT_IClb_Multiplier].name = "IClb_Multiplier";
//  P[VT_Jrel_Multiplier].name = "Jrel_Multiplier";
    

  /// state variables
    // Sodium current (INa, INaL)
    P[VT_Init_m].name = "Init_m";
    P[VT_Init_A_h].name = "Init_A_h";
    P[VT_Init_B_h].name = "Init_B_h";
    P[VT_Init_h].name = "Init_h";
    P[VT_Init_A_j].name = "Init_A_j";
    P[VT_Init_B_j].name = "Init_B_j";
    P[VT_Init_j].name = "Init_j";
    P[VT_Init_h_p].name = "Init_h_p";
    P[VT_Init_j_p].name = "Init_j_p";
    P[VT_Init_m_PKA].name = "Init_m_PKA";
    P[VT_Init_h_PKA].name = "Init_h_PKA";
    P[VT_Init_j_PKA].name = "Init_j_PKA";
    P[VT_Init_h_both].name = "Init_h_both";
    P[VT_Init_j_both].name = "Init_j_both";
    P[VT_Init_m_L].name = "Init_m_L";
    P[VT_Init_h_L].name = "Init_h_L";
    P[VT_Init_h_L_p].name = "Init_h_L_p";
    
    // L-type calcium current (I_CaL, I_CaNa, I_CaK)
    
    
    
//  P[VT_Init_Vm].name             = "Init_Vm";
//  P[VT_Init_Na_i].name           = "Init_Na_i";
//  P[VT_Init_Na_ss].name          = "Init_Na_ss";
//  P[VT_Init_K_i].name            = "Init_K_i";
//  P[VT_Init_K_ss].name           = "Init_K_ss";
//  P[VT_Init_Ca_i].name           = "Init_Ca_i";
//  P[VT_Init_Ca_ss].name          = "Init_Ca_ss";
//  P[VT_Init_Ca_nsr].name         = "Init_Ca_nsr";
//  P[VT_Init_Ca_jsr].name         = "Init_Ca_jsr";
//  P[VT_Init_Cl_i].name           = "Init_Cl_i";

//  P[VT_Init_a].name              = "Init_a";
//  P[VT_Init_i_fast].name         = "Init_i_fast";
//  P[VT_Init_i_slow].name         = "Init_i_slow";
//  P[VT_Init_a_CaMK].name         = "Init_a_CaMK";
//  P[VT_Init_i_CaMK_fast].name    = "Init_i_CaMK_fast";
//  P[VT_Init_i_CaMK_slow].name    = "Init_i_CaMK_slow";
//  P[VT_Init_d].name              = "Init_d";
//  P[VT_Init_f_fast].name         = "Init_f_fast";
//  P[VT_Init_f_slow].name         = "Init_f_slow";
//  P[VT_Init_f_Ca_fast].name      = "Init_f_Ca_fast";
//  P[VT_Init_f_Ca_slow].name      = "Init_f_Ca_slow";
//  P[VT_Init_j_Ca].name           = "Init_j_Ca";
//  P[VT_Init_n_ss].name           = "Init_n_ss";
//  P[VT_Init_n_i].name            = "Init_n_i";
//  P[VT_Init_f_CaMK_fast].name    = "Init_f_CaMK_fast";
//  P[VT_Init_f_Ca_CaMK_fast].name = "Init_f_Ca_CaMK_fast";
//  P[VT_Init_C_0].name            = "Init_C_0";
//  P[VT_Init_C_1].name            = "Init_C_1";
//  P[VT_Init_C_2].name            = "Init_C_2";
//  P[VT_Init_O].name              = "Init_O";
//  P[VT_Init_I].name              = "Init_I";
//  P[VT_Init_x_s1].name           = "Init_x_s1";
//  P[VT_Init_x_s2].name           = "Init_x_s2";
//  P[VT_Init_J_rel_NP].name       = "Init_J_rel_NP";
//  P[VT_Init_J_rel_CaMK].name     = "Init_J_rel_CaMK";
//  P[VT_Init_CaMK_trap].name      = "Init_CaMK_trap";



   
    
    /// PKA Fractions
    P[VT_fINa_PKA].name = "fINa_PKA";
    P[VT_fICaL_PKA].name = "fICaL_PKA";
    P[VT_fINaK_PKA].name = "fINaK_PKA";
    P[VT_fIKs_PKA].name = "fIKs_PKA";
    P[VT_fIfPLB_PKA].name = "fIfPLB_PKA";
    P[VT_fTnI_PKA].name = "fTnI_PKA";
    P[VT_fMyBPC_PKA].name = "fMyBPC_PKA";
    
  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Land-Niederer model of contraction
  ////////////////////////////////////////////////////////////////////////////////////////
  P[VT_perm50].name      = "perm50";
  P[VT_TRPN_n].name      = "TRPN_n";
  P[VT_koff].name        = "koff";
  P[VT_dr].name          = "dr";
  P[VT_wfrac].name       = "wfrac";
  P[VT_A_eff].name       = "A_eff";
  P[VT_ktm_unblock].name = "ktm_unblock";
  P[VT_beta_1].name      = "beta_1";
  P[VT_beta_0].name      = "beta_0";
  P[VT_gamma].name       = "gamma";
  P[VT_gamma_wu].name    = "gamma_wu";
  P[VT_phi].name         = "phi";
  P[VT_nperm].name       = "nperm";
  P[VT_ca50].name        = "ca50";
  P[VT_Tref].name        = "Tref";
  P[VT_nu].name          = "nu";
  P[VT_mu].name          = "mu";
  P[VT_xi].name          = "xi";
  P[VT_a].name           = "a";
  P[VT_k].name           = "k";
  P[VT_eta_l].name       = "eta_l";
  P[VT_eta_s].name       = "eta_s";
  P[VT_k_ws].name           = "VT_k_ws";
  P[VT_k_uw].name           = "VT_k_uw";
  P[VT_cdw].name            = "VT_cdw";
  P[VT_cds].name            = "VT_cds";
  P[VT_k_wu].name           = "VT_k_wu";
  P[VT_k_su].name           = "VT_k_su";
  P[VT_A].name              = "VT_A";
  P[VT_XS_init].name        = "XS_init";
  P[VT_XW_init].name        = "XW_init";
  P[VT_TRPN_init].name      = "TRPN_init";
  P[VT_TmBlocked_init].name = "TmBlocked_init";
  P[VT_ZETAS_init].name     = "ZETAS_init";
  P[VT_ZETAW_init].name     = "ZETAW_init";
  P[VT_Cd_init].name        = "Cd_init";
    
  P[VT_k_ws].readFromFile      = false;
  P[VT_k_uw].readFromFile      = false;
  P[VT_cdw].readFromFile       = false;
  P[VT_cds].readFromFile       = false;
  P[VT_k_wu].readFromFile      = false;
  P[VT_k_su].readFromFile      = false;
  P[VT_A].readFromFile         = false;
  P[VT_XSSS].readFromFile      = false;
  P[VT_XWSS].readFromFile      = false;
  P[VT_ktm_block].readFromFile = false;
    

  ParameterLoader EPL(initFile, EMT_TWorld);
  for (int x = 0; x < vtLast; x++)
    if (P[x].readFromFile)
      P[x].value = EPL.getParameterValue(P[x].name.c_str(), P[x].readFromFile);
  Calculate();
  InitTable(tinc);
#if KADEBUG
  cerr << "#Init() done ... \n";
#endif  // if KADEBUG
}  // TWorldParameters::Init

void TWorldParameters::Calculate() {
#if KADEBUG
  cerr << "#TWorldParameters - Calculate ..." << endl;
#endif  // if KADEBUG

  P[VT_k_ws].value = 0.004 * P[VT_mu].value * P[VT_xi].value;
  P[VT_k_uw].value = 0.026  * P[VT_nu].value;

  P[VT_cdw].value = P[VT_phi].value * P[VT_k_uw].value * (1.0-P[VT_dr].value)*(1.0-P[VT_wfrac].value) / ((1.0-P[VT_dr].value)*P[VT_wfrac].value);
  P[VT_cds].value = P[VT_phi].value * P[VT_k_ws].value * ((1.0-P[VT_dr].value)*P[VT_wfrac].value) / P[VT_dr].value;

  P[VT_k_wu].value = P[VT_k_uw].value * (1.0/P[VT_wfrac].value - 1.0) - P[VT_k_ws].value;
  P[VT_k_su].value = P[VT_k_ws].value * (1.0/P[VT_dr].value - 1.0) * P[VT_wfrac].value;

  P[VT_A].value = (P[VT_A_eff].value) * (P[VT_dr].value) / ((1.0-P[VT_dr].value) * P[VT_wfrac].value + P[VT_dr].value);

  P[VT_XSSS].value      = P[VT_dr].value * 0.5;
  P[VT_XWSS].value      = (1.0 - P[VT_dr].value) * P[VT_wfrac].value * 0.5;
  P[VT_ktm_block].value = P[VT_ktm_unblock].value * pow(P[VT_perm50].value, P[VT_nperm].value) * 0.5 / (0.5 - P[VT_XSSS].value - P[VT_XWSS].value);

}  // TWorldParameters::Calculate

void TWorldParameters::InitTable(ML_CalcType tinc) {
#if KADEBUG
  cerr << "#TWorldParameters - InitTable ..." << endl;
#endif  // if KADEBUG

  tinc *= 1000.0;  // second to millisecond conversion

  for (double V = -RangeTabhalf+.0001; V < RangeTabhalf; V += dDivisionTab) {
    // V in mV
    int Vi = (int)(DivisionTab*(RangeTabhalf+V)+.5);
  }
}  // TWorldParameters::InitTable
