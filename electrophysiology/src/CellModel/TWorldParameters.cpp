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

  P[VT_amplitude].name        = "amplitude";
  P[VT_duration].name         = "duration";
  P[VT_celltype].name         = "celltype";
  P[VT_Na_o].name             = "Na_o";
  P[VT_Ca_o].name             = "Ca_o";
  P[VT_K_o].name              = "K_o";
  P[VT_Cl_o].name             = "Cl_o";
  P[VT_INaF_Multiplier].name  = "INaF_Multiplier";
  P[VT_INaL_Multiplier].name  = "INaL_Multiplier";
  P[VT_Ito_Multiplier].name   = "Ito_Multiplier";
  P[VT_PCa_Multiplier].name   = "PCa_Multiplier";
  P[VT_IKr_Multiplier].name   = "IKr_Multiplier";
  P[VT_IKs_Multiplier].name   = "IKs_Multiplier";
  P[VT_IK1_Multiplier].name   = "IK1_Multiplier";
  P[VT_INaCa_Multiplier].name = "INaCa_Multiplier";
  P[VT_INaK_Multiplier].name  = "INaK_Multiplier";
  P[VT_IKb_Multiplier].name   = "IKb_Multiplier";
  P[VT_INab_Multiplier].name  = "INab_Multiplier";
  P[VT_ICab_Multiplier].name  = "ICab_Multiplier";
  P[VT_IpCa_Multiplier].name  = "IpCa_Multiplier";
  P[VT_Ito_slow_Multiplier].name = "Ito_slow_Multiplier";
  P[VT_Ito_fast_Multiplier].name = "Ito_fast_Multiplier";
  P[VT_ICaL_fractionSS].name = "ICaL_fractionSS";
  P[VT_CMDN_Multiplier].name = "CMDN_Multiplier";
  P[VT_INaCa_fractionSS].name = "INaCa_fractionSS";
    
  /// state variables
  P[VT_Init_Vm].name             = "Init_Vm";
  P[VT_Init_Na_i].name           = "Init_Na_i";
  P[VT_Init_Na_ss].name          = "Init_Na_ss";
  P[VT_Init_K_i].name            = "Init_K_i";
  P[VT_Init_K_ss].name           = "Init_K_ss";
  P[VT_Init_Ca_i].name           = "Init_Ca_i";
  P[VT_Init_Ca_ss].name          = "Init_Ca_ss";
  P[VT_Init_Ca_nsr].name         = "Init_Ca_nsr";
  P[VT_Init_Ca_jsr].name         = "Init_Ca_jsr";
  P[VT_Init_Cl_i].name           = "Init_Cl_i";
  P[VT_Init_m].name              = "Init_m";
  P[VT_Init_A_h].name            = "Init_A_h";
  P[VT_Init_B_h].name            = "Init_B_h";  
  P[VT_Init_h].name              = "Init_h";  
  P[VT_Init_A_j].name            = "Init_A_j";
  P[VT_Init_B_j].name            = "Init_B_j"; 
  P[VT_Init_j].name              = "Init_j";
  P[VT_Init_h_p].name            = "Init_h_p";
  P[VT_Init_j_p].name            = "Init_j_p";
  P[VT_Init_m_L].name            = "Init_m_L";
  P[VT_Init_h_L].name            = "Init_h_L";
  P[VT_Init_h_L_CaMK].name       = "Init_h_L_CaMK";
  P[VT_Init_a].name              = "Init_a";
  P[VT_Init_i_fast].name         = "Init_i_fast";
  P[VT_Init_i_slow].name         = "Init_i_slow";
  P[VT_Init_a_CaMK].name         = "Init_a_CaMK";
  P[VT_Init_i_CaMK_fast].name    = "Init_i_CaMK_fast";
  P[VT_Init_i_CaMK_slow].name    = "Init_i_CaMK_slow";
  P[VT_Init_d].name              = "Init_d";
  P[VT_Init_f_fast].name         = "Init_f_fast";
  P[VT_Init_f_slow].name         = "Init_f_slow";
  P[VT_Init_f_Ca_fast].name      = "Init_f_Ca_fast";
  P[VT_Init_f_Ca_slow].name      = "Init_f_Ca_slow";
  P[VT_Init_j_Ca].name           = "Init_j_Ca";
  P[VT_Init_n_ss].name           = "Init_n_ss";
  P[VT_Init_n_i].name            = "Init_n_i";
  P[VT_Init_f_CaMK_fast].name    = "Init_f_CaMK_fast";
  P[VT_Init_f_Ca_CaMK_fast].name = "Init_f_Ca_CaMK_fast";
  P[VT_Init_C_0].name            = "Init_C_0";
  P[VT_Init_C_1].name            = "Init_C_1";
  P[VT_Init_C_2].name            = "Init_C_2";
  P[VT_Init_O].name              = "Init_O";
  P[VT_Init_I].name              = "Init_I";
  P[VT_Init_x_s1].name           = "Init_x_s1";
  P[VT_Init_x_s2].name           = "Init_x_s2";
  P[VT_Init_J_rel_NP].name       = "Init_J_rel_NP";
  P[VT_Init_J_rel_CaMK].name     = "Init_J_rel_CaMK";
  P[VT_Init_CaMK_trap].name      = "Init_CaMK_trap";

  /// Parameters that are not read from .ev file
  P[VT_R].name              = "R";
  P[VT_T].name              = "T";
  P[VT_F].name              = "F";
  P[VT_L].name              = "L";
  P[VT_r].name              = "r";
  P[VT_tau_h_L].name        = "tau_h_L";
  P[VT_A_h].name            = "A_h";
  P[VT_v_cell].name         = "v_cell";
  P[VT_A_geo].name          = "A_geo";
  P[VT_A_cap].name          = "A_cap";
  P[VT_v_myo].name          = "v_myo";
  P[VT_v_nsr].name          = "v_nsr";
  P[VT_v_jsr].name          = "v_jsr";
  P[VT_v_ss].name           = "v_ss";
  P[VT_P_CaNa].name         = "P_CaNa";
  P[VT_P_CaK].name          = "P_CaK";
  P[VT_P_Ca_CaMK].name      = "P_Ca_CaMK";
  P[VT_P_CaNa_CaMK].name    = "P_CaNa_CaMK";
  P[VT_P_CaK_CaMK].name     = "P_CaK_CaMK";
  P[VT_A_f_fast].name       = "A_f_fast";
  P[VT_A_f_slow].name       = "A_f_slow";
  P[VT_k_Na1].name          = "k_Na1";
  P[VT_k_Na2].name          = "k_Na2";
  P[VT_k_asymm].name        = "k_asymm";
  P[VT_k_Ca_on].name        = "k_Ca_on";
  P[VT_k_Ca_off].name       = "k_Ca_off";
  P[VT_k_m_1].name          = "k_m_1";
  P[VT_k_p_2].name          = "k_p_2";
  P[VT_k_p_4].name          = "k_p_4";
  P[VT_MgADP].name          = "MgADP";
  P[VT_MgATP].name          = "MgATP";
  P[VT_K_MgATP].name        = "K_MgATP";
  P[VT_beta_tau].name       = "beta_tau";
  P[VT_P_Ca].name           = "P_Ca";
  P[VT_h_10].name           = "h_10";
  P[VT_h_11].name           = "h_11";
  P[VT_h_12].name           = "h_12";
  P[VT_k_1].name            = "k_1";
  P[VT_k_2].name            = "k_2";
  P[VT_k_5].name            = "k_5";
  P[VT_beta_1].name         = "beta_1";
  P[VT_alpha_2].name        = "alpha_2";
  P[VT_alpha_4].name        = "alpha_4";
  P[VT_alpha_rel].name      = "alpha_rel";
  P[VT_beta_tau_CaMK].name  = "beta_tau_CaMK";
  P[VT_alpha_rel_CaMK].name = "alpha_rel_CaMK";
  P[VT_tau_h_LCaMK].name    = "tau_h_LCaMK";
  P[VT_RToverF].name        = "RToverF";

  P[VT_R].readFromFile              = false;
  P[VT_T].readFromFile              = false;
  P[VT_F].readFromFile              = false;
  P[VT_L].readFromFile              = false;
  P[VT_r].readFromFile              = false;
  P[VT_tau_h_L].readFromFile        = false;
  P[VT_A_h].readFromFile            = false;
  P[VT_v_cell].readFromFile         = false;
  P[VT_A_geo].readFromFile          = false;
  P[VT_A_cap].readFromFile          = false;
  P[VT_v_myo].readFromFile          = false;
  P[VT_v_nsr].readFromFile          = false;
  P[VT_v_jsr].readFromFile          = false;
  P[VT_v_ss].readFromFile           = false;
  P[VT_P_CaNa].readFromFile         = false;
  P[VT_P_CaK].readFromFile          = false;
  P[VT_P_Ca_CaMK].readFromFile      = false;
  P[VT_P_CaNa_CaMK].readFromFile    = false;
  P[VT_P_CaK_CaMK].readFromFile     = false;
  P[VT_A_f_fast].readFromFile       = false;
  P[VT_A_f_slow].readFromFile       = false;
  P[VT_k_Na1].readFromFile          = false;
  P[VT_k_Na2].readFromFile          = false;
  P[VT_k_asymm].readFromFile        = false;
  P[VT_k_Ca_on].readFromFile        = false;
  P[VT_k_Ca_off].readFromFile       = false;
  P[VT_k_m_1].readFromFile          = false;
  P[VT_k_p_2].readFromFile          = false;
  P[VT_k_p_4].readFromFile          = false;
  P[VT_MgADP].readFromFile          = false;
  P[VT_MgATP].readFromFile          = false;
  P[VT_K_MgATP].readFromFile        = false;
  P[VT_beta_tau].readFromFile       = false;
  P[VT_P_Ca].readFromFile           = false;
  P[VT_h_10].readFromFile           = false;
  P[VT_h_11].readFromFile           = false;
  P[VT_h_12].readFromFile           = false;
  P[VT_k_1].readFromFile            = false;
  P[VT_k_2].readFromFile            = false;
  P[VT_k_5].readFromFile            = false;
  P[VT_beta_1].readFromFile         = false;
  P[VT_alpha_2].readFromFile        = false;
  P[VT_alpha_4].readFromFile        = false;
  P[VT_alpha_rel].readFromFile      = false;
  P[VT_beta_tau_CaMK].readFromFile  = false;
  P[VT_alpha_rel_CaMK].readFromFile = false;
  P[VT_tau_h_LCaMK].readFromFile    = false;
  P[VT_RToverF].readFromFile        = false;

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

//  P[VT_R].value        = 8314.0;
//  P[VT_v_cell].value =  1000.0 * 3.14 * P[VT_r].value * P[VT_r].value * P[VT_L].value;

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
