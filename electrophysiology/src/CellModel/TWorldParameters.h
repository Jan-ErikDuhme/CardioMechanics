/*
 * File: TWorldParameters.h
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


#ifndef TWORLDPARAMETERS_H
#define TWORLDPARAMETERS_H

/// TWorld model variations.

#include <ParameterLoader.h>

namespace NS_TWorldParameters {
enum varType {
  VT_amplitude = vtFirst,
  VT_duration,
  VT_celltype,
  VT_Na_o,
  VT_Ca_o,
  VT_K_o,
  VT_Cl_o,
  VT_INaF_Multiplier,
  VT_INaL_Multiplier,
  VT_Ito_Multiplier,
  VT_PCa_Multiplier,
  VT_IKr_Multiplier,
  VT_IKs_Multiplier,
  VT_IK1_Multiplier,
  VT_INaCa_Multiplier,
  VT_INaK_Multiplier,
  VT_IKb_Multiplier,
  VT_INab_Multiplier,
  VT_ICab_Multiplier,
  VT_IpCa_Multiplier,
  VT_Ito_slow_Multiplier,
  VT_Ito_fast_Multiplier,
  VT_ICaL_fractionSS,
  VT_CMDN_Multiplier,
  VT_INaCa_fractionSS,
  VT_Jup_Multiplier,
  VT_ICaCl_Multiplier,
  VT_IClb_Multiplier,
  VT_Jrel_Multiplier,
    
    
  /// PKA fractions
  VT_fINa_PKA,
  VT_fICaL_PKA,
  VT_fINaK_PKA,
  VT_fIKs_PKA,
  VT_fIfPLB_PKA,
  VT_fTnI_PKA,
  VT_fMyBPC_PKA,

  /// state variables
  VT_Init_Vm,
  VT_Init_Na_i,
  VT_Init_Na_ss,
  VT_Init_K_i,
  VT_Init_K_ss,
  VT_Init_Ca_i,
  VT_Init_Ca_ss,
  VT_Init_Ca_nsr,
  VT_Init_Ca_jsr,
  VT_Init_Cl_i,
    
  // Sodium current (INa, INaL)
  VT_Init_m,
  VT_Init_A_h,
  VT_Init_B_h,
  VT_Init_h,
  VT_Init_A_j,
  VT_Init_B_j,  
  VT_Init_j,
  VT_Init_h_p,
  VT_Init_j_p,
  VT_Init_m_PKA,
  VT_Init_h_PKA,
  VT_Init_j_PKA,
  VT_Init_h_both,
  VT_Init_j_both,
  VT_Init_m_L,
  VT_Init_h_L,
  VT_Init_h_L_p,
    
  // L-type calcium current (I_CaL, I_CaNa, I_CaK)
    
    
  VT_Init_a,
  VT_Init_i_fast,
  VT_Init_i_slow,
  VT_Init_a_CaMK,
  VT_Init_i_CaMK_fast,
  VT_Init_i_CaMK_slow,
  VT_Init_d,
  VT_Init_f_fast,
  VT_Init_f_slow,
  VT_Init_f_Ca_fast,
  VT_Init_f_Ca_slow,
  VT_Init_j_Ca,
  VT_Init_n_ss,
  VT_Init_n_i,
  VT_Init_f_CaMK_fast,
  VT_Init_f_Ca_CaMK_fast,
  VT_Init_C_0,
  VT_Init_C_1,
  VT_Init_C_2,
  VT_Init_O,
  VT_Init_I,
  VT_Init_x_s1,
  VT_Init_x_s2,
  VT_Init_J_rel_NP,
  VT_Init_J_rel_CaMK,
  VT_Init_CaMK_trap,
    
    
  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Land-Niederer model of contraction
  ////////////////////////////////////////////////////////////////////////////////////////
  VT_perm50,
  VT_TRPN_n,
  VT_koff,
  VT_dr,
  VT_wfrac,
  VT_A_eff,
  VT_ktm_unblock,
  VT_beta_1,
  VT_beta_0,
  VT_gamma,
  VT_gamma_wu,
  VT_phi,
  VT_nperm,
  VT_ca50,
  VT_Tref,
  VT_nu,
  VT_mu,
  VT_xi,
  VT_a,
  VT_k,
  VT_eta_l,
  VT_eta_s,
    
  VT_k_ws,
  VT_k_uw,
  VT_cdw,
  VT_cds,
  VT_k_wu,
  VT_k_su,
  VT_A,
  VT_XSSS,
  VT_XWSS,
  VT_ktm_block,
    
  VT_XS_init,
  VT_XW_init,
  VT_TRPN_init,
  VT_TmBlocked_init,
  VT_ZETAS_init,
  VT_ZETAW_init,
  VT_Cd_init,
    
  vtLast
};
}  // namespace NS_TWorldParameters

using namespace NS_TWorldParameters;

class TWorldParameters : public vbNewElphyParameters {
 public:
  TWorldParameters(const char *, ML_CalcType);
  ~TWorldParameters();
  void PrintParameters();
  void Calculate();
  void InitTable(ML_CalcType);
  void Init(const char *, ML_CalcType);

  ML_CalcType m_inf[RTDT];
  ML_CalcType exptau_m[RTDT];
  ML_CalcType h_inf[RTDT];
  ML_CalcType exptau_h[RTDT];
  ML_CalcType j_inf[RTDT];
  ML_CalcType exptau_j[RTDT];
  ML_CalcType h_p_inf[RTDT];
  ML_CalcType exptau_h_p[RTDT];
  ML_CalcType j_p_inf[RTDT];
  ML_CalcType exptau_j_p[RTDT];
  ML_CalcType m_L_inf[RTDT];
  ML_CalcType exptau_m_L[RTDT];
  ML_CalcType h_L_inf[RTDT];
  ML_CalcType exptau_h_L[RTDT];
  ML_CalcType h_L_CaMK_inf[RTDT];
  ML_CalcType exptau_h_L_CaMK[RTDT];
  ML_CalcType a_inf[RTDT];
  ML_CalcType exptau_a[RTDT];
  ML_CalcType i_inf[RTDT];
  ML_CalcType exptau_i_fast[RTDT];
  ML_CalcType exptau_i_slow[RTDT];
  ML_CalcType A_i_fast[RTDT];
  ML_CalcType A_i_slow[RTDT];
  ML_CalcType a_CaMK_inf[RTDT];
  ML_CalcType exptau_a_CaMK[RTDT];
  ML_CalcType i_CaMK_inf[RTDT];
  ML_CalcType exptau_i_CaMK_fast[RTDT];
  ML_CalcType exptau_i_CaMK_slow[RTDT];
  ML_CalcType d_inf[RTDT];
  ML_CalcType exptau_d[RTDT];
  ML_CalcType f_inf[RTDT];
  ML_CalcType exptau_f_fast[RTDT];
  ML_CalcType exptau_f_slow[RTDT];
  ML_CalcType f_Ca_inf[RTDT];
  ML_CalcType exptau_f_Ca_fast[RTDT];
  ML_CalcType exptau_f_Ca_slow[RTDT];
  ML_CalcType A_f_Ca_fast[RTDT];
  ML_CalcType A_f_Ca_slow[RTDT];
  ML_CalcType j_Ca_inf[RTDT];
  ML_CalcType exptau_j_Ca[RTDT];
  ML_CalcType f_CaMK_inf[RTDT];
  ML_CalcType exptau_f_CaMK_fast[RTDT];
  ML_CalcType f_Ca_CaMK_inf[RTDT];
  ML_CalcType exptau_f_Ca_CaMK_fast[RTDT];
  ML_CalcType x_r_inf[RTDT];
  ML_CalcType exptau_x_r_fast[RTDT];
  ML_CalcType exptau_x_r_slow[RTDT];
  ML_CalcType A_x_r_fast[RTDT];
  ML_CalcType A_x_r_slow[RTDT];
  ML_CalcType C_0[RTDT];
  ML_CalcType C_1[RTDT];
  ML_CalcType C_2[RTDT];
  ML_CalcType O[RTDT];
  ML_CalcType I[RTDT];
  ML_CalcType R_Kr[RTDT];
  ML_CalcType x_s1_inf[RTDT];
  ML_CalcType exptau_x_s1[RTDT];
  ML_CalcType x_s2_inf[RTDT];
  ML_CalcType exptau_x_s2[RTDT];
  ML_CalcType x_K1_inf[RTDT];
  ML_CalcType exptau_x_K1[RTDT];
  ML_CalcType R_K1[RTDT];
  ML_CalcType x_Kb[RTDT];
};  // class TWorldParameters

#endif  // ifndef TWORLDPARAMETERS_H
