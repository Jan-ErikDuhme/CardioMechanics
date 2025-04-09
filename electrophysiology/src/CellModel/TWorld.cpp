/*
 * File: TWorld.cpp
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


#include <TWorld.h>

TWorld::TWorld(TWorldParameters *pp) {
  ptTeaP = pp;
#ifdef HETERO
  PS = new ParameterSwitch(ptTeaP, NS_TWorldParameters::vtLast);
#endif  // ifdef HETERO
  Init();
}

TWorld::~TWorld() {}

#ifdef HETERO

inline bool TWorld::AddHeteroValue(string desc, double val) {
  Parameter EP(desc, val);

  return PS->addDynamicParameter(EP);
}

#else  // ifdef HETERO

inline bool TWorld::AddHeteroValue(string desc, double val) {
  throw kaBaseException("compile with HETERO to use this feature!\n");
}

#endif  // ifdef HETERO

inline int TWorld::GetSize(void) {
  return sizeof(TWorld)-sizeof(vbElphyModel<ML_CalcType>)-sizeof(TWorldParameters *)
#ifdef HETERO
         -sizeof(ParameterSwitch *)
#endif  // ifdef HETERO
  ;
}

inline unsigned char TWorld::getSpeed(ML_CalcType adVm) {
  return (unsigned char)5;
}

void TWorld::Init() {
#if KADEBUG
  cerr << "#initializing Class: TWorld ... " << endl;
        #endif  // if KADEBUG

//  m              = v(VT_Init_m);
//  A_h            = v(VT_Init_A_h);
//  B_h            = v(VT_Init_B_h);
//  h              = v(VT_Init_h);
//  j              = v(VT_Init_j);
//  h_p            = v(VT_Init_h_p);
//  j_p            = v(VT_Init_j_p);
//  m_L            = v(VT_Init_m_L);
//  h_L            = v(VT_Init_h_L);
//  h_L_CaMK       = v(VT_Init_h_L_CaMK);
//  a              = v(VT_Init_a);
//  i_fast         = v(VT_Init_i_fast);
//  i_slow         = v(VT_Init_i_slow);
//  a_CaMK         = v(VT_Init_a_CaMK);
//  i_CaMK_fast    = v(VT_Init_i_CaMK_fast);
//  i_CaMK_slow    = v(VT_Init_i_CaMK_slow);
//  d              = v(VT_Init_d);
//  f_fast         = v(VT_Init_f_fast);
//  f_slow         = v(VT_Init_f_slow);
//  f_Ca_fast      = v(VT_Init_f_Ca_fast);
//  f_Ca_slow      = v(VT_Init_f_Ca_slow);
//  j_Ca           = v(VT_Init_j_Ca);
//  n_ss           = v(VT_Init_n_ss);
//  n_i            = v(VT_Init_n_i);
//  f_CaMK_fast    = v(VT_Init_f_CaMK_fast);
//  f_Ca_CaMK_fast = v(VT_Init_f_Ca_CaMK_fast);
//  C_0            = v(VT_Init_C_0);
//  C_1            = v(VT_Init_C_1);
//  C_2            = v(VT_Init_C_2);
//  O              = v(VT_Init_O);
//  I              = v(VT_Init_I);
//  x_s1           = v(VT_Init_x_s1);
//  x_s2           = v(VT_Init_x_s2);
//  Na_i           = v(VT_Init_Na_i);
//  Na_ss          = v(VT_Init_Na_ss);
//  K_i            = v(VT_Init_K_i);
//  K_ss           = v(VT_Init_K_ss);
//  Ca_i           = v(VT_Init_Ca_i);
//  Ca_ss          = v(VT_Init_Ca_ss);
//  Ca_nsr         = v(VT_Init_Ca_nsr);
//  Ca_jsr         = v(VT_Init_Ca_jsr);
//  Cl_i           = v(VT_Init_Cl_i);
//  CaMK_trap      = v(VT_Init_CaMK_trap);
//  J_rel_NP       = v(VT_Init_J_rel_NP);
//  J_rel_CaMK     = v(VT_Init_J_rel_CaMK);
    
  
    fINa_PKA = v(VT_Init_fINa_PKA);
    fICaL_PKA = v(VT_Init_fICaL_PKA);
    fINaK_PKA = v(VT_Init_fINaK_PKA);
    fIKs_PKA = v(VT_Init_fIKs_PKA);
    fPLB_PKA = v(VT_Init_fPLB_PKA);
    fTnI_PKA = v(VT_Init_fTnI_PKA);
    fMyBPC_PKA = v(VT_Init_fMyBPC_PKA);
    ICaL_fractionSS = v(VT_Init_ICaL_fractionSS)
    
}  // TWorld::Init

ML_CalcType TWorld::Calc(double tinc,  ML_CalcType V,  ML_CalcType i_external,  ML_CalcType stretch, int euler) {
  tinc *= 1000.0;  // second to millisecond conversion
  ML_CalcType V_m = V * 1000.0;

  const int Vi = (int)(DivisionTab*(RangeTabhalf+V_m)+.5);  // array position

  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Model parameters
  ////////////////////////////////////////////////////////////////////////////////////////
  // Constants
  double R = 8314.0; // [J/kmol*K]
  double F = 96485.0; // [C/mol]
  double T = 310.0; // [K]
  double FoRT = F / R / T;
  double Cmem = 1.3810e-10; // [F] membrane capacitance
  double Qpow = (T - 310.0) / 10.0;
    
  // Cell Tometry
  double cellLength = 100.0; // cell length [um]
  double cellRadius = 10.25; // cell radius [um]
  double Vcell = 3.14 * pow(cellRadius, 2) * cellLength * 1e-15; // [L]
  double Vmyo = 0.65 * Vcell;
  double Vsr = 0.035 * Vcell;
  double Vsl = 0.02 * Vcell;
  double Vjunc = 0.0539 * 0.01 * Vcell;
    
  // Fractional currents in compartments
  double Fdyad = 0.11;
  double Fsl = 1 - Fdyad;
    
  // Fixed ion concentrations
  double Cl_o = 150.0; // Extracellular Cl  [mM]
  double Mg_i = 0.5; // Intracellular Mg  [mM]
    
  // Nerst Potentials
  const ML_CalcType E_Na_dyad = (1.0 / FoRT) * log(Na_o / Na_dyad);
  const ML_CalcType E_Na_sl = (1.0 / FoRT) * log(Na_o / Na_sl);
  const ML_CalcType E_K = (1.0 / FoRT) * log(K_o / K_myo);
  const ML_CalcType E_Ca_dyad = (1.0 / FoRT / 2.0) * log(Ca_o / Ca_dyad);
  const ML_CalcType E_Ca_sl = (1.0 / FoRT / 2.0) * log(Ca_o / Ca_sl);
  const ML_CalcType E_Cl = (1.0 / FoRT) * log(Cl_myo / Cl_o);

    
  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Buffering parameters
  ////////////////////////////////////////////////////////////////////////////////////////
  double Bmax_Naj = 7.561;
  double Bmax_Nasl = 1.65;
  double koff_na = 1e-3;
  double kon_na = 0.1e-3;
  double Bmax_TnClow = 70.0e-3;
  double Bmax_TnChigh = 140.0e-3;
  double koff_tnchca = 0.032e-3;
  double kon_tnchca = 2.37;
  double koff_tnchmg = 3.33e-3;
  double kon_tnchmg = 3.0e-3;
  double Bmax_CaM = 24.0e-3 * CMDN_Multiplier; //TODO: Add Multiplier to v()
  double koff_cam = 238.0e-3;
  double kon_cam = 34.0;
  double Bmax_myosin = 140.0e-3;
  double koff_myoca = 0.46e-3;
  double kon_myoca = 13.8;
  double koff_myomg = 0.057e-3;
  double kon_myomg = 0.0157;
  double Bmax_SR = 17.85854e-3;
  double koff_sr = 60.0e-3;
  double kon_sr = 100.0;
  double Bmax_SLlowsl = 33.923e-3 * Vmyo / Vsl;
  double Bmax_SLlowj = 4.89983e-4 * Vmyo / Vjunc;
  double koff_sll = 1300.0e-3;
  double kon_sll = 100.0;
  double Bmax_SLhighsl = 12.15423e-3 * Vmyo / Vsl;
  double Bmax_SLhighj = 1.75755e-4 * Vmyo / Vjunc;
  double koff_slh = 30.0e-3;
  double kon_slh = 100.0;
  double Bmax_Csqn = 136.55214e-3 * Vmyo / Vsr;
  double koff_csqn = 65.0;
  double kon_csqn = 100.0;
    
    
  /////////////////////////////////////////////////////////////////////////////////////////
  ///        CaMK and Ca signalling
  ////////////////////////////////////////////////////////////////////////////////////////
  double PP1_tot = 0.13698; //TODO: Change to parameter using v(), needs to be changed if PKA signalling is simulted
  double CaMK0  = 2.0 * 0.05; // Equilibrium fraction of active CaMKII binding sites
  double Km_CaMK_Ca = 5.0 * 0.0015; //[mmol/L] CaMKII affinity for Ca2+/CaM activation %Adjusted because of smaller cleft space
    
  const ML_CalcType CaMK_bound = CaMK0 * (1.0 - CaMK_trap) / (1 + Km_CaMK_Ca / Ca_dyad);
  const ML_CalyType CaMK_active = CaMK_bound + CaMK_trap; // Fraction of active CaMKII
    
  double alpha = 0.05;
  double beta  = 6.8e-4;
  const ML_CalcType dCaMK_trap = alpha * CaMK_bound * CaMK_active - beta * CaMK_trap * (0.1 + 0.9 * PP1_tot / 0.1371); //dCaMK_Trap/dt
  CaMK_trap += tinc * dCaMK_trap;
    
  double tau_plb = 100000.0; // [ms] Time constant of CaMKII PLB phosphorylation
  double tau_ryr = 10000.0; // [ms] Time constant of CaMKII RyR phosphorylation
  double tau_cal = tau_ryr; // Time constant of CaMKII ICaL phosphorylation
    
  double K_Phos_CaMK = 0.35;  // Affinity of PLB, ICaL etc for CaMKII
  const ML_CalcType CaMK_Phos_ss_ICaL =  CaMK_active / (CaMK_active + K_Phos_CaMK);
  const ML_CalcType CaMK_Phos_ss_RyR =  CaMK_active / (CaMK_active + 1.0);
  const ML_CalcType CaMK_Phos_ss_PLB =  CaMK_active / (CaMK_active + 10.0);

  const ML_CalcType dCaMK_f_ICaL = (CaMK_Phos_ss_ICaL - CaMK_f_ICaL) / tau_cal;
  CaMK_f_ICaL += tinc * dCaMK_f_ICaL;
  const ML_CalcType dCaMK_f_RyR = (CaMK_Phos_ss_RyR - CaMK_f_RyR)  / tau_ryr;
  CaMK_f_RyR += tinc * dCaMK_f_RyR;
  const ML_CalcType dCaMK_f_PLB = (CaMK_Phos_ss_PLB - CaMK_f_PLB)  / tau_plb;
  CaMK_f_PLB += tinc * dCaMK_f_PLB;
    
  double alpha_serca = 0.05;
  const ML_CalcType bound_serca = CaMK0 * (1.0 - casig_serca_trap) / (1.0 + Km_CaMK_Ca / Ca_dyad);
  const ML_CalcType casig_SERCA_act = bound_serca + casig_serca_trap; // Fraction of active CaMKII
  const ML_CalcType dcasig_serca_trap = alpha_serca * bound_serca * casig_SERCA_act - beta * casig_serca_trap * (0.1 + 0.9 * PP1_tot / 0.1371); //dCaMK_Trap/dt
  casig_serca_trap += tinc * dcasig_serca_trap;
    
    
  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Sodium current (INa, INaL)
  ////////////////////////////////////////////////////////////////////////////////////////
  ///////////// calulate I_Na //////////
  // m gate
  const ML_CalcType m_inf = 1.0 / ((1.0 + exp(-(V_m + 56.86)/9.03))*(1.0 + exp(-(V_m + 56.86)/9.03)));
  const ML_CalcType tau_m = 0.1292 * exp(-(((V_m + 45.79)/15.54) * ((V_m + 45.79)/15.54))) + 0.06487 * exp(-(((V_m - 4.823)/51.12) * ((V_m - 4.823)/51.12)));
  m = m_inf - (m_inf - m) * exp(-tinc / tau_m);

  // h and j gate
  if (V_m < -40.0) {
      A_h = (0.057*exp((-(V_m+80.0)/6.8)));
      B_h = ((2.7*exp((0.079*V_m)))+(3.1e5*exp((0.3485*V_m))));
      A_j = ((((-2.5428e4*exp((0.2444*V_m))) - (6.948e-6*exp((-0.04391*V_m))))*(V_m+37.78))/(1.0+exp((0.311*(V_m+79.23)))));
      B_j = ((0.02424*exp((-0.01052*V_m)))/(1.0+exp((-0.1378*(V_m+40.14)))));
  } else {
      A_h = 0.0;
      B_h = (0.77/(0.13*(1.+exp((-(V_m+10.66)/11.1)))));
      A_j = 0.0;
      B_j =((0.6*exp((0.057*V_m)))/(1.0+exp((-0.1*(V_m+32.0))))) ;
  }
  const ML_CalcType tau_h = 1.0/(A_h+B_h);
  const ML_CalcType h_inf = (1.0/((1.0+exp(((V_m+71.55)/7.43)))*(1.0+exp(((V_m+71.55)/7.43)))));
  h = h_inf - (h_inf - h) * exp(-tinc / tau_h);
    
  const ML_CalcType tau_j = 1.0/(A_j+B_j);
  const ML_CalcType j_inf = (1.0/((1.0+exp(((V_m+71.55)/7.43)))*(1.0+exp(((V_m+71.55)/7.43)))));;
  j = j_inf - (j_inf - j) * exp(-tinc / tau_j);

  // gating CaMK-P
  const ML_CalcType h_p_inf = (1.0/((1.0+exp(((V_m+77.55)/7.43)))*(1.0+exp(((V_m+77.55)/7.43)))));
  h_p = h_p_inf - (h_p_inf - h_p) * exp(-tinc / tau_h);
  const ML_CalcType tau_j_p = 1.46 * tau_j;
  j_p = j_inf - (j_inf - j_p) * exp(-tinc / tau_j_p);
    
  // gating PKA
  const ML_CalcType m_PKA_inf = 1.0 / ((1.0 + exp(-(V_m + 56.86)/9.03))*(1.0 + exp(-(V_m + 56.86)/9.03)));
  m_PKA = m_PKA_inf - (m_PKA_inf - m_PKA) * exp(-tinc / tau_m);
      
  const ML_CalcType h_PKA_inf = (1.0/((1.0+exp(((V_m+76.55)/7.43)))*(1.0+exp(((V_m+76.55)/7.43)))));
  h_PKA = h_PKA_inf - (h_PKA_inf - h_PKA) * exp(-tinc / tau_h);
    
  const ML_CalcType j_PKA_inf = h_PKA_inf;
  j_PKA = j_PKA_inf - (j_PKA_inf - j_PKA) * exp(-tinc / tau_j_p);
      
  // Both Phosphorylated
  const ML_CalcType h_both_inf = (1.0/((1.0+exp(((V_m+82.55)/7.43)))*(1.0+exp(((V_m+82.55)/7.43)))));
  h_both = h_both_inf - (h_both_inf - h_both) * exp(-tinc / tau_h);
  const ML_CalcType j_both_inf = h_both_inf;
  j_both = j_both_inf - (j_both_inf - j_both) * exp(-tinc / tau_j_p);
      
  // Putting together the channels behavior and fraction
  double G_Na = 22.08788 * v(VT_INaF_Multiplier);
  double G_Na_PKA = G_Na * 1.25;
      
  const ML_CalcType phi_INa_CaMK = CaMK_f_RyR;
  const ML_CalcType phi_INa_PKA = fINa_PKA;
  const ML_CalcType phi_INa_Both = phi_INa_CaMK * phi_INa_PKA;
  const ML_CalcType phi_INa_CaMKonly = phi_INa_CaMK - phi_INa_Both;
  const ML_CalcType phi_INa_PKAonly = phi_INa_PKA - phi_INa_Both;

  I_Na_Base_NP = G_Na * m*m*m * h * j; //Non-Phosphorylated
  I_Na_Base_CaMK = G_Na * m*m*m * h_p * j_p;
  I_Na_Base_PKA = G_Na_PKA * m_PKA*m_PKA*m_PKA * h_PKA * j_PKA;
  I_Na_Base_Both = G_Na_PKA * m_PKA*m_PKA*m_PKA * h_both * j_both;
      
  // 4 population
  I_Na_Base = ((1-phi_INa_CaMKonly-phi_INa_PKAonly-phi_INa_Both)*I_Na_Base_NP + phi_INa_CaMKonly*I_Na_Base_CaMK + phi_INa_PKAonly*I_Na_Base_PKA + phi_INa_Both*I_Na_Base_Both);
  I_NaFast_junc = Fdyad * I_Na_Base*(V_m - E_Na_dyad);
  I_NaFast_sl = (1 - dyad) * I_Na_Base * (V_m - E_Na_sl);

  ///////////// calulate I_NaL //////////
  // m gate
  const ML_CalcType m_L_inf = 1.0 / (1.0 + exp(-(V_m + 42.85) / 5.264));
  const ML_CalcType tau_m_L = tau_m;
  m_L = m_L_inf - (m_L_inf - m_L) * exp(-tinc / tau_m_L);
    
  // h gate
  const ML_CalcType h_L_inf = 1.0 / (1.0 + exp((V_m + 87.61) / 7.488));
  h_L = h_L_inf - (h_L_inf - h_L) * exp(-tinc / v(VT_tau_h_L));
    
  // gating CaMK-P
  const ML_CalcType h_L_p_inf = 1.0 / (1.0 + exp((V_m + 93.81) / 7.488));
  const ML_CalcType tau_h_L_p = 3.0 * v(VT_tau_h_L);
  h_L_p = h_L_p_inf - (h_L_p_inf - h_L_p) * exp(-tinc / tau_h_L_p);
    
  // Putting together the channels behavior and fraction
  const ML_CalcType phi_INaL_CaMK = phi_INa_CaMK;
  double G_Na_L = 0.04229 * v(VT_INaL_Multiplier) * (1.0 + phi_INaL_CaMK); //TODO: Ist noch falsch
  I_NaL_junc = Fjunc  * G_Na_L * (V_m - E_Na_junc) * m_L * ((1.0 - phi_INaL_CaMK) * h_L + phi_INaL_CaMK * h_L_p);;
  I_NaL_sl = Fsl * G_Na_L * (V_m - E_Na_sl) * m_L * ((1.0 - phi_INaL_CaMK) * h_L + phi_INaL_CaMK * h_L_p);;

  ///////////// Combinations I_Na and I_NaL //////////
  I_Na_junc = I_NaFast_junc + I_NaL_junc ;
  I_Na_sl = I_NaFast_sl + I_NaL_sl ;
  I_NaFast = I_NaFast_junc+I_NaFast_sl;
  I_NaL = I_NaL_junc+I_NaL_sl;
    
    
  /////////////////////////////////////////////////////////////////////////////////////////
  ///        L-type calcium current (I_CaL, I_CaNa, I_CaK)
  ////////////////////////////////////////////////////////////////////////////////////////

  
    // d gate
    const ML_CalcType d_inf = min(1.0763*exp((-1.007*exp((-0.0829*(V_m+3.62483))))), 1.0);
    const ML_CalcType tau_d = 1.5 + 1.0 / (exp(-0.05 * (V_m + 6.0)) + exp(0.09 * (V_m + 14.0)));
    d = d_inf - (d_inf -d) * exp(-tinc / tau_d);
    
    // f gate
    const ML_CalcType f_inf = 1.0 / (1.0 + exp((V_m + 19.58) / 3.696));
    const ML_CalcType tau_f_fast = 6.17111 + 1.0 / (0.00126 * exp(-(V_m + 26.63596) / 9.69961) + 0.00126 * exp((V_m + 26.63596) / 9.69961));
    const ML_CalcType tau_f_slow = 2719.22489 + 1.0 / (7.19411e-05 * exp(-(V_m + 5.74631) / 10.87690) + 7.19411e-05 * exp((V_m + 5.74631) / 16.31535));
    
    const ML_CalcType A_f_fast = 0.52477;
    const ML_CalcType A_f_slow = 1.0 - A_f_fast;
    f_fast = f_inf - (f_inf - f_fast) * exp(-tinc / tau_f_fast);
    f_slow = f_inf - (f_inf - f_slow) * exp(-tinc / tau_f_slow);
    const ML_CalcType f = A_f_fast * f_fast + A_f_slow * f_slow;
    
    const ML_CalcType f_Ca_inf      = f_inf;
    const ML_CalcType tau_f_Ca_fast = 13.50673 + 1.0 / (0.15420 * exp(-(V_m - 1.31611) / 11.33960) + 0.15420 * exp((V_m - 1.31611) / 11.33960));
    const ML_CalcType tau_f_Ca_slow = 177.95813 + 1.0 / (4.73955e-04 * exp((-V_m+0.79049) / 0.81777) + 4.73955e-04 * exp((V_m+2.40474) /1.90812));
    const ML_CalcType A_f_Ca_fast = 0.3 + 0.6 / (1.0 + exp((V_m - 9.24247) / 27.96201));
    const ML_CalcType A_f_Ca_slow = 1.0 - A_f_Ca_fast;
    f_Ca_fast = f_Ca_inf - (f_Ca_inf - f_Ca_fast) * exp(-tinc / tau_f_Ca_fast);
    f_Ca_slow = f_Ca_inf - (f_Ca_inf - f_Ca_slow) * exp(-tinc / tau_f_Ca_slow);
    const ML_CalcType f_Ca = A_f_Ca_fast * f_Ca_fast + A_f_Ca_slow * f_Ca_slow;
    
    // j gate
    double tau_j_Ca = 66.0;
    const ML_CalcType j_Ca_inf = 1. / (1 + exp((V_m + 17.66945) / (3.21501)));
    j_Ca = j_Ca_inf - (j_Ca_inf - j_Ca) * exp(-tinc / tau_j_Ca);
    
    // gating CaMK-P
    const ML_CalcType tau_f_p_fast = 2.5 * tau_f_fast;
    const ML_CalcType f_p_inf = f_inf;
    f_p_fast = f_p_inf - (f_p_inf - f_p_fast) * exp(-tinc / tau_f_p_fast);
    const ML_CalcType f_p = A_f_fast * f_p_fast + A_f_slow * f_slow;
    
    const ML_CalcType tau_f_Ca_p_fast = 2.5 * tau_f_Ca_fast;
    const ML_CalcType f_Ca_p_inf = f_inf;
    f_Ca_p_fast = f_Ca_p_inf - (f_Ca_p_inf - f_Ca_p_fast) * exp(-tinc / tau_f_Ca_p_fast);
    const ML_CalcType f_Ca_p = A_f_Ca_fast * f_Ca_p_fast + A_f_Ca_slow * f_Ca_slow;
    
    // SS nca
    double Kmn = 0.00222;
    double k2n = 957.85903;
    const ML_CalcType k_m2_n = j_Ca * 0.84191;
    const ML_CalcType alpha_n_Ca_dyad   = 1.0 / ((k2n / k_m2_n) + pow((1.0 + (Kmn / Ca_dyad)), 3.80763));
    const ML_CalcType dn_Ca_dyad = alpha_n_Ca_dyad * k2n - n_Ca_dyad * k_m2_n;
    n_Ca_dyad += tinc * dn_Ca_dyad;
    
    // myoplasmic nca
    const ML_CalcType alpha_n_Ca_sl   = 1.0 / ((k2n / k_m2_n) + pow((1.0 + (Kmn / Ca_sl)), 3.80763));
    const ML_CalcType dn_Ca_sl = alpha_n_Ca_sl * k2n - n_Ca_sl * k_m2_n;
    n_Ca_sl += tinc * dn_Ca_sl;
    
    // SS driving force
    const ML_CalcType I_o = (0.5 * (v(VT_Na_o) + v(VT_K_o) + v(VT_Cl_o) + (4.0 * v(VT_Ca_o))) / 1000.0);
    const ML_CalcType I_sl = ((0.5 * (Na_sl + K_sl + Cl_i + (4.* Ca_ss))) / 1000.0);
    
    double constA = (1.82e6*pow((74.*310.),-1.5)); // Diel constant and temperature as constants
    
   
    
    
   
    
    
    
    
    

    
    
    
    
    
    // Channel conductances
    double P_Ca = 1.5768e-04 * v(VT_IpCa_Multiplier);
    double P_Ca_PKA = P_Ca * 1.9; // BetaAdrenergic;  (%Gong--> 1.9)
    if (v(VT_celltype) == 1.0) { // epi
        P_Ca = P_Ca * 1.025;
        P_Ca_PKA = P_Ca_P * 1.025;
    } else if (v(VT_celltype) == 2.0) { // mid
        P_Ca = P_Ca * 1.1;
        P_Ca_PKA = P_Ca_P * 1.1;
    }
    double P_CaNa = 1.1737 / 1.8969 * 0.00125 * P_Ca;
    double P_CaK = 1.1737 / 1.8969 * 3.574e-4 * P_Ca;
    double P_Ca_p = 1.1 * P_Ca;
    double P_CaNa_p = 1.1737 / 1.8969 * 0.00125 * P_Ca_p;
    double P_CaK_p = 1.1737 / 1.8969 * 3.574e-4 * P_Ca_p;
    double P_CaNa_PKA = 1.1737 / 1.8969 * 0.00125 * P_Ca_PKA;
    double P_CaK_PKA = 1.1737 / 1.8969 * 3.574e-4 * P_Ca_PKA;
      
    
    // Putting together the channels behavior and fraction
    const ML_CalcType phi_ICaL_CaMK = CaMK_f_ICaL;
    const ML_CalcType phi_ICaL_PKA = fICaL_PKA;
    const ML_CalcType phi_ICaL_Both = phi_ICaL_CaMK * phi_ICaL_PKA;
    const ML_CalcType phi_ICaL_CaMKonly = phi_ICaL_CaMK - phi_ICaL_Both;
    const ML_CalcType phi_ICaL_PKAonly = phi_ICaL_PKA - phi_ICaL_Both;
    
    I_CaL_dyad_NP = ;
    I_CaL_dyad_CaMK = ;
    I_CaL_dyad_PKA = ;
    I_CaL_dyad_Both = ;
    I_CaL_sl_NP = ;
    I_CaL_sl_CaMK = ;
    I_CaL_sl_PKA = ;
    I_CaL_sl_Both = ;
    
    I_CaNa_dyad_NP = ;
    I_CaNa_dyad_CaMK = ;
    I_CaNa_dyad_PKA = ;
    I_CaNa_dyad_Both = ;
    I_CaNa_sl_NP = ;
    I_CaNa_sl_CaMK = ;
    I_CaNa_sl_PKA = ;
    I_CaNa_sl_Both = ;
    
    I_CaK_dyad_NP = ;
    I_CaK_dyad_CaMK = ;
    I_CaK_dyad_PKA = ;
    I_CaK_dyad_Both = ;
    I_CaK_sl_NP = ;
    I_CaK_sl_CaMK = ;
    I_CaK_sl_PKA = ;
    I_CaK_sl_Both = ;
    
    //4 population combination
    I_CaL_dyad = ;
    I_CaNa_dyad = ;
    I_CaK_dyad = ;
    I_CaL_sl = ;
    I_CaNa_sl = ;
    I_CaK_sl = ;
    
    // Weigh I_CaL in sl and dyad
    I_CaL_sl = I_CaL_sl * (1 - ICaL_fractionSS);
    I_CaNa_sl = I_CaNa_sl * (1 - ICaL_fractionSS);
    I_CaK_sl = I_CaK_sl * (1 - ICaL_fractionSS);
    I_CaL_dyad = I_CaL_dyad * ICaL_fractionSS;
    I_CaNa_dyad = I_CaNa_dyad * ICaL_fractionSS;;
    I_CaK_dyad = I_CaK_dyad * ICaL_fractionSS;;
    
    
    
    
    
//    /// calculate I_CaL, I_CaNa, I_CaK

//

//
//
//
//
//
//    const ML_CalcType gamma_Ca_ss = exp((-constA*4.*((sqrt(I_ss)/(1.+sqrt(I_ss))) - (0.3*I_ss))));
//    const ML_CalcType gamma_Ca_o = exp((-constA*4.*((sqrt(I_o)/(1.+sqrt(I_o))) - (0.3*I_o))));
//    const ML_CalcType gamma_Na_ss = exp((-constA*1.*((sqrt(I_ss)/(1.+sqrt(I_ss))) - (0.3*I_ss))));
//    const ML_CalcType gamma_Na_o = exp((-constA*1.*((sqrt(I_o)/(1.+sqrt(I_o))) - (0.3*I_o))));
//    const ML_CalcType gamma_K_ss = exp((-constA*1.*((sqrt(I_ss)/(1.+sqrt(I_ss))) - (0.3*I_ss))));
//    const ML_CalcType gamma_K_o = exp((-constA*1.*((sqrt(I_o)/(1.+sqrt(I_o))) - (0.3*I_o))));
//    const ML_CalcType Psi_Ca_ss = ((4.*VFFoverRT*((gamma_Ca_ss*Ca_ss*exp((2.*VFoverRT))) - (gamma_Ca_o*v(VT_Ca_o))))/(exp((2.*VFoverRT)) - 1.));
//    const ML_CalcType Psi_CaNa_ss = ((1.*VFFoverRT*((gamma_Na_ss*Na_ss*exp((1.*VFoverRT))) - (gamma_Na_o)*v(VT_Na_o)))/(exp((1.*VFoverRT)) - 1.));
//    const ML_CalcType Psi_CaK_ss = ((1.*VFFoverRT*((gamma_K_ss*K_ss*exp((1.*VFoverRT))) - (gamma_K_o*v(VT_K_o))))/(exp((1.*VFoverRT)) - 1.));
//
//    double PCa_b = 8.3757e-05;
//    const ML_CalcType PCa = v(VT_IpCa_Multiplier) * PCa_b;
//    const ML_CalcType PCaNa = (0.00125*PCa);
//    const ML_CalcType PCaK = (3.574e-4*PCa);
//    const ML_CalcType PCap = (1.1*PCa);
//    const ML_CalcType PCaNap = (0.00125*PCap);
//    const ML_CalcType PCaKp = (3.574e-4*PCap);
//
//    double ICaL_fractionSS = 0.8;
//    const ML_CalcType phi_ICaL_CaMK = 1.0 / (1.0 + KmCaMK / CaMK_active);
//
//    I_CaL_ss = (ICaL_fractionSS*(((1. - phi_ICaL_CaMK)*PCa*Psi_Ca_ss*d*((f*(1. - n_ss))+(j_Ca*f_Ca*n_ss)))+(phi_ICaL_CaMK*PCap*Psi_Ca_ss*d*((f_CaMK*(1. - n_ss))+(j_Ca*f_Ca_CaMK*n_ss)))));
//    I_CaNa_ss = (ICaL_fractionSS*(((1. - phi_ICaL_CaMK)*PCaNa*Psi_CaNa_ss*d*((f*(1. - n_ss))+(j_Ca*f_Ca*n_ss)))+(phi_ICaL_CaMK*PCaNap*Psi_CaNa_ss*d*((f_CaMK*(1. - n_ss))+(j_Ca*f_Ca_CaMK*n_ss)))));
//    I_CaK_ss = (ICaL_fractionSS*(((1. - phi_ICaL_CaMK)*PCaK*Psi_CaK_ss*d*((f*(1. - n_ss))+(j_Ca*f_Ca*n_ss)))+(phi_ICaL_CaMK*PCaKp*Psi_CaK_ss*d*((f_CaMK*(1. - n_ss))+(j_Ca*f_Ca_CaMK*n_ss)))));
//
//
//    const ML_CalcType alpha_n_Ca_i   = 1.0 / ((k2n / k_m2_n) + ((1. + (Kmn / Ca_i))*(1. + (Kmn / Ca_i))*(1. + (Kmn / Ca_i))*(1. + (Kmn / Ca_i))));
//    n_i = alpha_n_Ca_i * (k2n / k_m2_n) - (alpha_n_Ca_i * (k2n / k_m2_n) - n_i) * exp(-k_m2_n * tinc);
//
//
//    const ML_CalcType I_i = ((0.5*(Na_i+K_i+Cl_i+(4.*Ca_i)))/1000.);
//    const ML_CalcType gamma_Ca_i = exp((-constA*4.*((sqrt(I_i)/(1.+sqrt(I_i))) - (0.3*I_i))));
//    const ML_CalcType gamma_Na_i = exp((-constA*1.*((sqrt(I_i)/(1.+sqrt(I_i))) - (0.3*I_i))));
//    const ML_CalcType gamma_K_i = exp((-constA*1.*((sqrt(I_i)/(1.+sqrt(I_i))) - (0.3*I_i))));
//    const ML_CalcType Psi_Ca_i = ((4.*VFFoverRT*((gamma_Ca_i*Ca_i*exp((2.*VFoverRT))) - (gamma_Ca_o*v(VT_Ca_o))))/(exp((2.*VFoverRT)) - 1.));
//    const ML_CalcType Psi_CaNa_i = (1.*VFFoverRT*((gamma_Na_i*Na_i*exp((1.*VFoverRT))) - (gamma_Na_o*v(VT_Na_o)))/(exp((1.*VFoverRT)) - 1.));
//    const ML_CalcType Psi_CaK_i = ((1.*VFFoverRT*((gamma_K_i*K_i*exp((1.*VFoverRT))) - (gamma_K_o*v(VT_K_o))))/(exp((1.*VFoverRT)) - 1.));
//
//    I_CaL_i = ((1. - ICaL_fractionSS)*(((1. - phi_ICaL_CaMK)*PCa*Psi_Ca_i*d*((f*(1. - n_i))+(j_Ca*f_Ca*n_i)))+(phi_ICaL_CaMK*PCap*Psi_Ca_i*d*((f_CaMK*(1. - n_i))+(j_Ca*f_Ca_CaMK*n_i)))));
//    I_CaNa_i = ((1. - ICaL_fractionSS)*(((1. - phi_ICaL_CaMK)*PCaNa*Psi_CaNa_i*d*((f*(1. - n_i))+(j_Ca*f_Ca*n_i)))+(phi_ICaL_CaMK*PCaNap*Psi_CaNa_i*d*((f_CaMK*(1. - n_i))+(j_Ca*f_Ca_CaMK*n_i)))));
//    I_CaK_i = ((1. - ICaL_fractionSS)*(((1. - phi_ICaL_CaMK)*PCaK*Psi_CaK_i*d*((f*(1. - n_i))+(j_Ca*f_Ca*n_i)))+(phi_ICaL_CaMK*PCaKp*Psi_CaK_i*d*((f_CaMK*(1. - n_i))+(j_Ca*f_Ca_CaMK*n_i)))));
//
//
//    I_CaL = I_CaL_ss + I_CaL_i;
//    I_CaNa = I_CaNa_ss + I_CaNa_i;
//    I_CaK = I_CaK_ss + I_CaK_i;
  
    
  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Transient outward current (Ito)
  ////////////////////////////////////////////////////////////////////////////////////////
  // activation (a gate)
  const ML_CalcType a_inf = 1.0 / (1.0 + exp(-(V_m - 19.0) / 13.0));; //Info: changed x -> a and y -> i to better match Tomek and OHara syntax
  const ML_CalcType tau_a_slow = 9.0 / (1.0 + exp((V_m + 3.0) / 15.0)) + 0.5;
  a_slow = a_inf - (a_inf - a_slow) * exp(-tinc / tau_a_slow);
    
  const ML_CalcType tau_a_fast = 8.5 * exp(-pow(((V_m + 45.0) / 50.0), 2)) + 0.5;
  a_fast = a_inf - (a_inf - a_fast) * exp(-tinc / tau_a_fast);
    
  // inactivation (i gate)
  const ML_CalcType i_inf = 1.0 / (1.0 + exp((V_m + 19.5) / 5.0));
  const ML_CalcType tau_i_slow = 800.0 / (1.0 +exp((V_m + 60.0) / 10.0)) + 30.0;
  i_slow = i_inf - (i_inf - i_slow) * exp(-tinc / tau_i_slow);
    
  const ML_CalcType tau_i_fast = 85.0 * exp(-pow((V_m+40), 2/220)) + 7.0;
  i_fast = i_inf - (i_inf - i_fast) * exp(-tinc / tau_i_fast);
    
  // gating CaMK-P
  const ML_CalcType a_p_inf = 1.0 / (1.0 + exp(-(V_m - 29.0) / 13.0));
  a_p_slow = a_p_inf - (a_p_inf - a_p_slow) * exp(-tinc / tau_a_slow);
  a_p_fast = a_p_inf - (a_p_inf - a_p_fast) * exp(-tinc / tau_a_fast);
    
  const ML_CalcType delta_p_develop = 1.354 + 1.0e-04 / (exp((V_m - 167.4) / 15.89) + exp(-(V_m - 12.23) / 0.2154));
  const ML_CalcType delta_p_recover = 1.0 - 0.5 / (1.0 + exp((V_m + 70.0) / 20.0));
  const ML_CalcType tau_i_p_slow = tau_i_slow * delta_p_develop * delta_p_recover;
  const ML_CalcType tau_i_p_fast = tau_i_fast * delta_p_develop * delta_p_recover;
  i_p_slow = i_inf - (i_inf - i_p_slow) * exp(-tinc / tau_i_p_slow);
  i_p_fast = i_inf - (i_inf - i_p_fast) * exp(-tinc / tau_i_p_fast);
    
  // Putting together the channels behavior and fraction
  if (v(VT_celltype) == 1.0) { // epi
      double G_to_slow = 0.02036 * v(VT_Ito_slow_Multiplier);
      double G_to_fast = 0.29856 * v(VT_Ito_fast_Multiplier);
  } else if (v(VT_celltype) == 2.0) { // mid
      double G_to_slow = 0.04632 * v(VT_Ito_slow_Multiplier);
      double G_to_fast = 0.14928 * v(VT_Ito_fast_Multiplier);
  } else if (v(VT_celltype) == 0.0) { //endo
      double G_to_slow = 0.07210 * v(VT_Ito_slow_Multiplier);
      double G_to_fast = 0.01276 * v(VT_Ito_fast_Multiplier);
  }
  const ML_CalcType phi_Ito_CaMK = phi_INa_CaMK; //TODO: Hierbei bin ich mir unsicher aus Tomek übernommen (MATLAB -> fItop = camk_f_RyR)
  I_to_slow = v(VT_Ito_Multiplier) * G_to_slow * (V_m - E_K) * ((1.0 - phi_Ito_CaMK) * a_slow * i_slow + phi_Ito_CaMK * a_p_slow * i_p_slow);
  I_to_fast = v(VT_Ito_Multiplier) * G_to_fast * (V_m - E_K) * ((1.0 - phi_Ito_CaMK) * a_fast * i_fast + phi_Ito_CaMK * a_p_fast * i_p_fast);
  I_to = I_to_slow + I_to_fast;

  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Rapid delayed rectifier current (IKr)
  ////////////////////////////////////////////////////////////////////////////////////////
  const ML_CalcType alpha_Kr = 0.1161 * exp(0.299 * VFoverRT);
  const ML_CalcType beta_Kr = 0.2442 * exp(-1.604 * VFoverRT);
  double alpha_Kr_1 = 1.25 * 0.1235;
  double beta_Kr_1 = 0.1911;
  const ML_CalcType alpha_Kr_2 = 0.0578 * exp(0.971 * VFoverRT);
  const ML_CalcType beta_Kr_2 = 0.349e-3 * exp(-1.062 * VFoverRT);
  const ML_CalcType alpha_Kr_i = 0.2533 * exp(0.5953 * VFoverRT);
  const ML_CalcType beta_Kr_i = 0.04568 * exp(-0.8209 * VFoverRT);
  const ML_CalcType alpha_C2_to_I = 0.52e-4 * exp(1.525 * VFoverRT);
  const ML_CalcType beta_I_to_C2 = (beta_Kr_2 * beta_Kr_i * alpha_C2_to_I) / (alpha_Kr_2 * alpha_Kr_i);
  const ML_CalcType dC_0 = C_1 * beta_Kr - C_0 * alpha_Kr;
  C_0 += tinc * dC_0;
  const ML_CalcType dC_1 = C_0 * alpha_Kr + C_2 * beta_Kr_1 - C_1 * (beta_Kr + alpha_Kr_1);
  C_1 += tinc * dC_1;
  const ML_CalcType dC_2 = C_1 * alpha_Kr_1 + O * beta_Kr_2 + I * beta_I_to_C2 - C_2 * (beta_Kr_1 + alpha_Kr_2 + alpha_C2_to_I);
  C_2 += tinc * dC_2;
  const ML_CalcType dO = (((alpha_Kr_2 * C_2) + (beta_Kr_i * I)) - ((beta_Kr_2 + alpha_Kr_i) * O));
  O += tinc * dO;
  const ML_CalcType dI = (((alpha_C2_to_I * C_2) + (alpha_Kr_i * O)) - ((beta_I_to_C2 + beta_Kr_i) * I));
  I += tinc * dI;
    
  // Putting together the channels behavior and fraction
  double G_Kr = 0.043 * sqrt(v(VT_K_o) / 5.0) * v(VT_IKr_Multiplier);
  if (v(VT_celltype) == 1.0) { // epi
      G_Kr = G_Kr * 1.25;
  } else if (v(VT_celltype) == 2.0) { // mid
      G_Kr = G_Kr * 0.7;
  }
  I_Kr = G_Kr * O * (V_m - E_K);
    
  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Slow delayed rectifier current (IKs)
  ////////////////////////////////////////////////////////////////////////////////////////
  // gating
  const ML_CalcType k_PKA = phi_IKs_PKA; //TODO: No idea what that is
  const ML_CalcType V_h_0 = -1.0 - 10.0 * k_PKA;
  const ML_CalcType V_h_max = -12.0 - 9.0 * k_PKA;
  const ML_CalcType tau_0 = 26.0 + 9.0 * k_PKA;
  const ML_CalcType tau_max = 40.0 + 4.0 * k_PKA;
    
  const ML_CalcType V_junc = V_h_0 + (V_h_max - V_h_0) / (1.0 + pow((350e-06 / Ca_junc), 4.0)); // Regulated by PKA
  const ML_CalcType xs_junc_inf = 1.0 / (1.0 + exp(-(V_m - V_junc) / 25.0));
  const ML_CalcType V_tau_junc = tau_0 + (tau_max - tau_0) / (1.0 + pow((150e-06/Ca_junc), 3.0)); // Regulated by PKA
  const ML_CalcType tau_xs_junc = 2.0 * (50.0 + (50.0 + 350.0 * exp(-(pow((V_m + 30.0), 2.0)) / 4000.0)) * 1.0 / (1.0 + exp(-(V_m + V_tau_junc) / 10.0)));
  xs_junc = xs_junc_inf - (xs_junc_inf - xs_junc) * exp(-tinc / tau_xs_junc);
    
  const ML_CalcType V_sl = V_h_0 + (V_h_max - V_h_0) / (1.0 + pow((350e-06 / Ca_sl), 4.0)); // Regulated by PKA
  const ML_CalcType xs_sl_inf = 1.0 / (1.0 + exp(-(V_m - V_sl) / 25.0));
  const ML_CalcType V_tau_sl = tau_0 + (tau_max - tau_0) / (1.0 + pow((150e-06/Ca_sl), 3.0)); // Regulated by PKA
  const ML_CalcType tau_xs_sl = 2.0 * (50.0 + (50.0 + 350.0 * exp(-(pow((V_m + 30.0), 2.0)) / 4000.0)) * 1.0 / (1.0 + exp(-(V_m + V_tau_sl) / 10.0)));
  xs_sl = xs_sl_inf - (xs_sl_inf - xs_sl) * exp(-tinc / tau_xs_sl);
    
  // conductances
  double G_Ks_factor_SA = 2.97002 * v(VT_IKs_Multiplier);
  if (v(VT_celltype) == 2.0) { // mid
      G_Ks_factor_SA = 0.5 * G_Ks_factor_SA;
  }
  double G_Ks_factor = 0.01;
  const ML_CalcType G_Ks_0 = G_Ks_factor * (0.2 + 0.2 * k_PKA);
  const ML_CalcType G_Ks_max = G_Ks_factor * (0.8 + 7.0 * k_PKA);
  const ML_CalcType G_Ks_junc = G_Ks_0 + (G_Ks_max - G_Ks_0) / (1.0 + pow((150e-06 / Ca_junc), 1.3)); // Regulated by PKA
  const ML_CalcType G_Ks_sl = G_Ks_0 + (G_Ks_max - G_Ks_0) / (1.0 + pow((150e-06 / Ca_sl), 1.3)); // Regulated by PKA
    
  // Nernst potential
  double pNaK = 0.01833;
  const ML_CalcType E_Ks = v(VT_RToverF) * log((v(VT_K_o) + pNaK * v(VT_Na_o)) / (K_i + pNaK * Na_i));
    
  // Putting together the channels behavior and fraction
  I_Ks_junc = Fjunc * G_Ks_factor_SA * G_Ks_junc * xs_junc * xs_junc * (V_m - E_Ks); //TODO: Fjunc was?
  I_Ks_sl = Fsl * G_Ks_factor_SA * G_Ks_sl * xs_sl * xs_sl * (V_m - E_Ks); //TODO: Fsl was?
  I_Ks = I_Ks_junc + I_Ks_sl;
    
  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Inward rectifier current (IK1)
  ////////////////////////////////////////////////////////////////////////////////////////
  // gating
  const ML_CalcType a_K1 = 4.094 / (1.0 + exp(0.1217 * (V_m - E_K - 49.934)));
  const ML_CalcType b_K1 = (15.72 * exp(0.0674 * (V_m - E_K - 3.257)) + exp(0.0618 * (V_m - E_K - 594.31))) / (1.0 + exp(-0.1629 * (V_m - E_K + 14.207)));
  const ML_CalcType K1_SS = (a_K1) / (a_K1 + b_K1);
    
  // Putting together the channels behavior and fraction
  double G_K1 = 0.6992 + v(VT_IK1_Multiplier);
  if (v(VT_celltype) == 1.0) { // epi
      G_K1 = G_K1 * 1.1;
  } else if (v(VT_celltype) == 2.0) { // mid
      G_K1 = G_K1 * 1.3;
  }
  I_K1 = G_K1 * sqrt(v(VT_K_o)/5.0) * K1_SS * (V_m - E_K);
    
  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Sodium-calcium exchanger (INaCa)
  ////////////////////////////////////////////////////////////////////////////////////////
  double z_Ca = 2.0;
  double kna1 = 11.9712;
  double kna2 = 2.76;
  double kna3 = 88.767;
  double kasymm = 19.4258;
  double wna  = 3.2978e+04;
  double wca = 5.1756e+04;
  double wnaca  = 2.7763e+03;
  double kcaon = 3.4164e+06;
  double kcaoff = 3.8532e+03;
  double qna = 0.6718;
  double qca  = 0.0955;
  const ML_CalcType h_Ca = exp(qca * (V_m - 8.3117) * VFoverRT);
  const ML_CalcType h_Na = exp(qna * (V_m - 8.3117) * VFoverRT);
    
  // calculate I_NaCa_i
  const ML_CalcType h_1_i = 1.0 + (Na_i * (1.0 + h_Na)) / kna3;
  const ML_CalcType h_2_i = (Na_i * h_Na) / (kna3 * h_1_i);
  const ML_CalcType h_3_i = 1.0 / h_1_i;
  const ML_CalcType h_4_i = 1.0 + (Na_i * (1.0 + Na_i / kna2)) / kna1;
  const ML_CalcType h_5_i = (Na_i * Na_i) / (h_4_i * kna1 * kna2);
  const ML_CalcType h_6_i = 1.0 / h_4_i;
  const ML_CalcType h_7_i = 1.0 + (v(VT_Na_o) * (1.0 + 1.0 / h_Na)) / kna3;
  const ML_CalcType h_8_i = v(VT_Na_o) / (kna3 * h_Na * h_7_i);
  const ML_CalcType h_9_i = 1.0 / h_7_i;
  const ML_CalcType h_10_i = kasymm + 1.0 + v(VT_Na_o) / kna1 * (1.0 + v(VT_Na_o) / kna2);
  const ML_CalcType h_11_i = v(VT_Na_o) * v(VT_Na_o) / (h_10_i * kna1 * kna2);
  const ML_CalcType h_12_i = 1.0 / h_10_i;
    
  const ML_CalcType k_1_i = h_12_i * v(VT_Ca_o) * kcaon;
  const ML_CalcType k_2_i = kcaoff;
  const ML_CalcType k_3_d_i = h_9_i * wca;
  const ML_CalcType k_3_dd_i = h_8_i * wnaca;
  const ML_CalcType k_3_i = k_3_d_i + k_3_dd_i;
  const ML_CalcType k_4_d_i = (h_3_i * wca) / h_Ca;
  const ML_CalcType k_4_dd_i = h_2_i * wnaca;
  const ML_CalcType k_4_i = k_4_d_i + k_4_dd_i;
  const ML_CalcType k_5_i = kcaoff;
  const ML_CalcType k_6_i = h_6_i * Ca_i * kcaon;
  const ML_CalcType k_7_i = h_5_i * h_2_i * wna;
  const ML_CalcType k_8_i = h_8_i * h_11_i * wna;
 
  const ML_CalcType x_1_i = k_2_i * k_4_i * (k_7_i + k_6_i) + k_5_i * k_7_i * (k_2_i + k_3_i);
  const ML_CalcType x_2_i = k_1_i * k_7_i * (k_4_i + k_5_i) + k_4_i * k_6_i * (k_1_i + k_8_i);
  const ML_CalcType x_3_i = k_1_i * k_3_i * (k_7_i + k_6_i) + k_8_i * k_6_i * (k_2_i + k_3_i);
  const ML_CalcType x_4_i = k_2_i * k_8_i * (k_4_i + k_5_i) + k_3_i * k_5_i * (k_1_i + k_8_i);
    
  const ML_CalcType E_1_i = x_1_i / (x_1_i + x_2_i + x_3_i + x_4_i);
  const ML_CalcType E_2_i = x_2_i / (x_1_i + x_2_i + x_3_i + x_4_i);
  const ML_CalcType E_3_i = x_3_i / (x_1_i + x_2_i + x_3_i + x_4_i);
  const ML_CalcType E_4_i = x_4_i / (x_1_i + x_2_i + x_3_i + x_4_i);
    
  double KmCaAct = 150.0e-6;
  const ML_CalcType allo_i = 1.0 / (1.0 + ((KmCaAct / Ca_i) * (KmCaAct / Ca_i)));
  double z_Na = 1.0;
  const ML_CalcType J_NaCa_Na_i = 3.0 * (E_4_i * k_7_i - E_1_i * k_8_i) + E_3_i * k_4_dd_i - E_2_i * k_3_dd_i;
  const ML_CalcType J_NaCa_Ca_i = E_2_i * k_2_i - E_1_i * k_1_i;
  double G_NaCa = 0.00179 * v(VT_INaCa_Multiplier);
    
  I_NaCa_i = G_NaCa * allo_i * (z_Na * J_NaCa_Na_i + z_Ca * J_NaCa_Ca_i) * (1-INaCa_fractionSS); //TODO: Was ist INaCa_fractionSS
    
  // calculate I_NaCa_ss
  const ML_CalcType h_1_ss = 1.0 + (Na_sl * (1.0 + h_Na)) / kna3;
  const ML_CalcType h_2_ss = (Na_sl * h_Na) / (kna3 * h_1_ss);
  const ML_CalcType h_3_ss = 1.0 / h_1_ss;
  const ML_CalcType h_4_ss = 1.0 + (Na_sl * (1.0 + Na_sl / kna2)) / kna1;
  const ML_CalcType h_5_ss = (Na_sl * Na_sl) / (h_4_ss * kna1 * kna2);
  const ML_CalcType h_6_ss = 1.0 / h_4_ss;
  const ML_CalcType h_7_ss = 1.0 + v(VT_Na_o) / kna3 * (1.0 + 1.0 / h_Na);
  const ML_CalcType h_8_ss = v(VT_Na_o) / (kna3 * h_Na * h_7_ss);
  const ML_CalcType h_9_ss = 1.0 / h_7_ss;
  const ML_CalcType h_10_ss = kasymm + 1.0 + v(VT_Na_o) / kna1 * (1 + v(VT_Na_o) / kna2);
  const ML_CalcType h_11_ss = v(VT_Na_o) * v(VT_Na_o) / (h_10_ss * kna1 * kna2);
  const ML_CalcType h_12_ss = 1.0 / h_10_ss;
    
  const ML_CalcType k_1_ss = h_12_ss * v(VT_Ca_o) * kcaon;
  const ML_CalcType k_2_ss = kcaoff;
  const ML_CalcType k_3_d_ss = h_9_ss * wca;
  const ML_CalcType k_3_dd_ss = h_8_ss * wnaca;
  const ML_CalcType k_3_ss = k_3_d_ss + k_3_dd_ss;
  const ML_CalcType k_4_d_ss = (h_3_ss * wca) / h_Ca;
  const ML_CalcType k_4_dd_ss = h_2_ss * wnaca;
  const ML_CalcType k_4_ss = k_4_d_ss + k_4_dd_ss;
  const ML_CalcType k_5_ss = kcaoff;
  const ML_CalcType k_6_ss = h_6_ss * Ca_sl * kcaon;
  const ML_CalcType k_7_ss = h_5_ss * h_2_ss * wna;
  const ML_CalcType k_8_ss = h_8_ss * h_11_ss * wna;
    
  const ML_CalcType x_1_ss = k_2_ss * k_4_ss * (k_7_ss + k_6_ss) + k_5_ss * k_7_ss * (k_2_ss + k_3_ss);
  const ML_CalcType x_2_ss = k_1_ss * k_7_ss * (k_4_ss + k_5_ss) + k_4_ss * k_6_ss * (k_1_ss + k_8_ss);
  const ML_CalcType x_3_ss = k_1_ss * k_3_ss * (k_7_ss + k_6_ss) + k_8_ss * k_6_ss * (k_2_ss + k_3_ss);
  const ML_CalcType x_4_ss = k_2_ss * k_8_ss * (k_4_ss + k_5_ss) + k_3_ss * k_5_ss * (k_1_ss + k_8_ss);

  const ML_CalcType E_1_ss = x_1_ss / (x_1_ss + x_2_ss + x_3_ss + x_4_ss);
  const ML_CalcType E_2_ss = x_2_ss / (x_1_ss + x_2_ss + x_3_ss + x_4_ss);
  const ML_CalcType E_3_ss = x_3_ss / (x_1_ss + x_2_ss + x_3_ss + x_4_ss);
  const ML_CalcType E_4_ss = x_4_ss / (x_1_ss + x_2_ss + x_3_ss + x_4_ss);
    
  const ML_CalcType allo_ss = 1.0 / (1.0 + (KmCaAct / Ca_ss) * (KmCaAct / Ca_ss)); //TODO: Unklar ob Ca_ss oder Ca_sl
  const ML_CalcType J_NaCa_Na_ss = 3.0 * (E_4_ss * k_7_ss - E_1_ss * k_8_ss) + E_3_ss * k_4_dd_ss - E_2_ss * k_3_dd_ss;
  const ML_CalcType J_NaCa_Ca_ss = E_2_ss * k_2_ss - E_1_ss * k_1_ss;
    
  I_NaCa_ss = INaCa_fractionSS * G_NaCa * allo_ss * (z_Na * J_NaCa_Na_ss + z_Ca * J_NaCa_Ca_ss); //TODO: Was ist INaCa_fractionSS

  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Sodium-potassium pump (INaK)
  ////////////////////////////////////////////////////////////////////////////////////////
  double I_bar_NaK = 2.10774 * v(VT_INaK_Multiplier);
  double KmNaip = 11;
  double KmNaip_PKA = 8.4615;
  double KmKo = 1.5;
    
  const ML_CalcType f_NaK = 0.75 + (0.00375 - ((140 - v(VT_Na_o)) / 50) * 0.001) * V_m; //Varying the slope mainly based on https://rupress.org/jgp/article-pdf/94/3/539/1814046/539.pdf
    
  I_NaK_junc_noPKA = Fjunc * I_bar_NaK * f_NaK * v(VT_K_o) / (1.0 + pow((KmNaip / Na_junc), 4.0)) /(v(VT_K_o) + KmKo);
  I_NaK_junc_PKA = Fjunc * I_bar_NaK * f_NaK * v(VT_K_o) / (1.0 + pow((KmNaip_PKA / Na_junc),4.0)) /(v(VT_K_o) + KmKo);
  I_NaK_sl_noPKA = Fsl * I_bar_NaK * f_NaK * v(VT_K_o) / (1.0 + pow((KmNaip / Na_sl), 4.0)) / (v(VT_K_o) + KmKo);
  I_NaK_sl_PKA = Fsl * I_bar_NaK * f_NaK * v(VT_K_o) / (1.0 + pow((KmNaip_PKA / Na_sl), 4.0)) / (v(VT_K_o) + KmKo);
    
  I_NaK_junc = (1.0 - phi_INaK_PKA) * I_NaK_junc_noPKA + phi_INaK_PKA * I_NaK_junc_PKA; //TODO: Wo kommt phi_INaK_PKA her???
  I_NaK_sl = (1.0 - phi_INaK_PKA) * I_NaK_sl_noPKA + phi_INaK_PKA * I_NaK_sl_PKA;
  I_NaK = I_Nak_junc + I_NaK_sl;
    
  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Chloride currents (ICaCl, IClb)
  ////////////////////////////////////////////////////////////////////////////////////////
  // calculate I_CaCl
  double G_CaCl = 0.01615 * v(VT_ICaCl_Multiplier);
  double Kd_CaCl = 100e-03;
    
  I_CaCl_junc = 0.5 * Fjunc * G_CaCl / (1.0 + Kd_CaCl / Ca_junc)*(V_m - E_Cl);
  I_CaCl_sl = 0.5 * Fsl * G_CaCl / (1.0 + Kd_CaCl / Ca_sl)*(V_m - E_Cl);
  I_ClCa = I_ClCa_junc+I_ClCa_sl;
    
  // calculate I_Clb
  double G_Clb = 0.00241 * v(VT_IClb_Multiplier);
  I_Clb = G_Clb * (V_m - E_Cl);
  
  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Background currents (INab, ICab, IKb)
  ////////////////////////////////////////////////////////////////////////////////////////
  // calculate I_Nab
  double G_Nab = 2.0 * 0.297e-03 * v(VT_INab_Multiplier);
  I_Nab_junc = Fjunc * G_Nab * (V_m - E_Na_junc);
  I_Nab_sl = Fsl * G_Nab * (V_m - E_Na_sl);
  I_Nab = I_Nab_junc + I_Nab_sl;
    
  // calculate I_Kb
  double G_Kb = 0.010879 * v(VT_IKb_Multiplier);
  const ML_CalcType x_Kb = 1.0 / (1.0 + exp(-(V_m - 10.8968) / (23.9871)));
  I_Kb = G_Kb * x_Kb * (V_m - E_K);
    
  // calculate I_Cab
  double G_Cab = 5.15575e-04 * v(VT_ICab_Multiplier);
  I_Cab_junc = Fjunc * G_Cab * (V_m - E_ca_junc);
  I_Cab_sl = Fsl * G_Cab * (V_m - E_ca_sl);
  I_Cab = I_Cab_junc + I_Cab_sl;
    
  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Calcium release from SR (Jrel, Jleak)
  ////////////////////////////////////////////////////////////////////////////////////////
    
    double ks = 26.6 * v(VT_Jrel_Multiplier);
    double koCa = 23.87221;
    double kom = 0.16219;
    double kiCa = 0.39871;
    double kim = 0.04311;
    double ec50SR = 0.75385;
    double steepnessCaSR = 5.09473;
    double caExpFactor = 0.68655 ;
    double caTransFactor = 0.94428 ;
    double caExpFactor2 = 2.06273 ;
    double caTransFactor2 = 0.52967 ;

    ///////////// calculate directly I_CaL-coupled RyRy //////////
    ///
    
    // inactivation
    
    // slower inactivation
    
    
    
    ///////////// Main Ca-sensitive RyRs //////////
    
    
    // not phosphorylated by CaMKII
    
    
    // And also a version of phosphorylated
    
    
    // Total release
    
    
    // Additional leak
    
    
    
    
    
//    /// SR Calcuim release flux J_rel
//    double K_inf_rel = 1.7;
//    ML_CalcType J_rel_NP_inf = (v(VT_alpha_rel) * (-I_CaL)) / (1.0 + pow((K_inf_rel/ Ca_jsr), 8.0));
//    ML_CalcType J_rel_CaMK_inf = (v(VT_alpha_rel_CaMK) * (-I_CaL)) / (1.0 + pow((K_inf_rel/ Ca_jsr), 8.0));
//    if (v(VT_celltype) == 2.0) {
//      J_rel_NP_inf   *= 1.7;
//      J_rel_CaMK_inf *= 1.7;
//    }
//    const ML_CalcType tau_rel_NP_b = v(VT_beta_tau) / (1.0 + (0.0123 / Ca_jsr));
//    const ML_CalcType tau_rel_NP   = tau_rel_NP_b < 0.001 ? 0.001 : tau_rel_NP_b;
//    J_rel_NP = J_rel_NP_inf - (J_rel_NP_inf - J_rel_NP) * exp(-tinc / tau_rel_NP);
//    
//    const ML_CalcType tau_rel_CaMK_b = v(VT_beta_tau_CaMK) / (1.0 + (0.0123 / Ca_jsr));
//    const ML_CalcType tau_rel_CaMK   = tau_rel_CaMK_b < 0.001 ? 0.001 : tau_rel_CaMK_b;
//    J_rel_CaMK = J_rel_CaMK_inf - (J_rel_CaMK_inf - J_rel_CaMK) * exp(-tinc / tau_rel_CaMK);
//    
//    const ML_CalcType phi_rel_CaMK = phi_INa_CaMK;
//    J_rel        = 1.5378 * (1.0 - phi_rel_CaMK) * J_rel_NP + phi_rel_CaMK * J_rel_CaMK;
//
//    /// J_up
//    ML_CalcType J_up_NP   = (0.005425 * Ca_i) / (0.00092 + Ca_i);
//    ML_CalcType J_up_CaMK = (2.75 * 0.005425 * Ca_i) / (0.00092 - 0.00017 + Ca_i);
//    if (v(VT_celltype) == 1.0) {
//      J_up_NP   *= 1.3;
//      J_up_CaMK *= 1.3;
//    }
//    const ML_CalcType phi_up_CaMK = (1.0 /(1.0 + (KmCaMK / CaMK_active)));
//    J_leak      = (0.0048825 * Ca_nsr) / 15.0;
//

    
    
    
  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Calcium reuptake to the SR (Jup)
  ////////////////////////////////////////////////////////////////////////////////////////
    //    J_up = ((1.0 - phi_up_CaMK) * J_up_NP + (phi_up_CaMK * J_up_CaMK) - J_leak);
    
    
    
  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Sarcolemmal calcium pump (pCa)
  ////////////////////////////////////////////////////////////////////////////////////////
//    /// calculate I_pCa
//    double G_pCa = 0.0005;
//    I_pCa = v(VT_IpCa_Multiplier) * ((G_pCa * Ca_i) / (0.0005 + Ca_i));
    
    
    
    

//    double kmcmdn  = 0.00238;
//    double trpnmax = 0.07;
//    double kmtrpn  = 0.0005;
//    double BSRmax  = 0.047;
//    double KmBSR   = 0.00087;
//    double BSLmax  = 1.124;
//    double KmBSL   = 0.0087;
//    double csqnmax = 10.0;
//    double kmcsqn  = 0.8;
    
    
    kon_na
    Bmax_Naj
    koff_na
    Bmax_Nasl
    
    Bmax_TnClow
    
    
  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Buffering
  ////////////////////////////////////////////////////////////////////////////////////////
    const ML_CalcType dBuffer_NaBj = ; //TODO: Put equation
    Buffer_NaBj += tinc * dBuffer_NaBj;
    const ML_CalcType dBuffer_NaBsl = ; //TODO: Put equation
    Buffer_NaBsl += tinc * dBuffer_NaBsl;
    const ML_CalcType dBuffer_TnClow = ; //TODO: Put equation
    Buffer_TnClow += tinc * dBuffer_TnClow;
    const ML_CalcType dBuffer_TnCHc = ; //TODO: Put equation
    Buffer_TnCHc += tinc * dBuffer_TnCHc;
    const ML_CalcType dBuffer_TnCHm = ; //TODO: Put equation
    Buffer_TnCHm += tinc * dBuffer_TnCHm;
    const ML_CalcType dBuffer_CaM = ; //TODO: Put equation
    Buffer_CaM += tinc * dBuffer_CaM;
    const ML_CalcType dBuffer_Myosin_ca = ; //TODO: Put equation
    Buffer_Myosin_ca += tinc * dBuffer_Myosin_ca;
    const ML_CalcType dBuffer_Myosin_mg = ; //TODO: Put equation
    Buffer_Myosin_mg += tinc * dBuffer_Myosin_mg;
    const ML_CalcType dBuffer_SRB = ; //TODO: Put equation
    Buffer_SRB += tinc * dBuffer_SRB;
    const ML_CalcType dBuffer_SLLj = ; //TODO: Put equation
    Buffer_SLLj += tinc * dBuffer_SLLj;
    const ML_CalcType dBuffer_SLLsl = ; //TODO: Put equation
    Buffer_SLLsl += tinc * dBuffer_SLLsl;
    const ML_CalcType dBuffer_SLHj = ; //TODO: Put equation
    Buffer_SLHj += tinc * dBuffer_SLHj;
    const ML_CalcType dBuffer_SLHsl = ; //TODO: Put equation
    Buffer_SLHsl += tinc * dBuffer_SLHsl;
    const ML_CalcType dBuffer_Csqn = ; //TODO: Put equation
    Buffer_Csqn += tinc * dBuffer_Csqn;
    
    
  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Diffusion
  ////////////////////////////////////////////////////////////////////////////////////////
    double J_Ca_dyad_sl = 1.0 / 3.06685e12;
    double J_Ca_sl_myo = 1.0 / 0.74556e11;
    double J_Na_dyad_sl = 1.0 / (1.6382e12 / 3.0 * 100.0);
    double J_Na_sl_myo = 1.0 / (1.8308e10 / 3.0 * 100.0);
    
    const ML_CalcType J_CaBuffer_myo = ; //TODO: Put equation
    const ML_CalcType J_CaBuffer_dyad = ; //TODO: Put equation
    const ML_CalcType J_CaBuffer_sl = ; //TODO: Put equation

    
  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Total ion currents and concentration changes
  ////////////////////////////////////////////////////////////////////////////////////////
    // Sodium Concentration
    I_Na_tot_dyad = I_Na_dyad + I_Nab_dyad + 3 * I_NaCx_dyad + 3 * I_NaK_dyad + I_CaNa_dayd;
    I_Na_tot_sl = I_Na_sl + I_Nab_sl + 3 * I_NaCx_sl + 3 * I_NaK_sl + I_CaNa_sl;
    I_Na_tot = I_Na_tot_dyad + I_Na_tot_sl;
    const ML_CalcType dNa_dyad = ; //TODO: Put equation
    Na_dyad += tinc * dNa_dyad;
    const ML_CalcType dNa_sl = ; //TODO: Put equation
    Na_sl += tinc * dNa_sl;
    const ML_CalcType dNa_myo = ; //TODO: Put equation
    Na_myo += tinc * dNa_myo;
    
    // Potassium Concentration
    I_K_tot = I_to + I_Kr + I_Ks + I_K1 - (2 * I_NaK) + I_CaK + I_Kb + i_external;
    const ML_CalcType dK_myo = ; //TODO: Put equation
    K_myo += tinc * dK_myo;
    
    // Cloride Concentration
    I_Cl_tot = I_CaCl + I_Clb;
    const ML_CalcType dCl_myo = ; //TODO: Put equation
    Cl_myo += tinc * dCl_myo;
    
    // Calcium Concentration
    I_Ca_tot_dyad = I_Ca_dyad + I_Cab_dyad + I_pCa_dyad - (2 * I_NaCx_dyad);
    I_Ca_tot_sl = I_Ca_sl + I_Cab_sl + I_pCa_sl - (2 * I_NaCx_sl);
    I_Ca_tot = I_Ca_tot_dyad + I_Ca_tot_sl;
    const ML_CalcType dCa_dyad = ; //TODO: Put equation
    Ca_dyad += tinc * dCa_dyad;
    const ML_CalcType dCa_sl = ; //TODO: Put equation
    Ca_sl += tinc * dCa_sl;
    const ML_CalcType dCa_myo = ; //TODO: Put equation
    Ca_myo += tinc * dCa_myo;
    const ML_CalcType dCa_SR = ; //TODO: Put equation
    Ca_SR += tinc * dCa_SR;

    
  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Change Membrane Potential (I_tot)
  ////////////////////////////////////////////////////////////////////////////////////////
    I_tot = -(I_Na_tot + I_Cl_tot + I_Ca_tot + I_K_tot + i_external);

    
    
    
    
    
  return 0.001 * tinc * I_tot;
}  // TWorld::Calc

void TWorld::Print(ostream &tempstr, double tArg,  ML_CalcType V) {
  tempstr << tArg << ' ' << V << ' ' <<
    m << ' '  << A_h << ' '  << B_h << ' '  << h << ' ' << A_j << ' '  << B_j << ' '<< j << ' ' << h_p << ' ' << j_p << ' ' << m_L << ' ' <<
    h_L << ' ' << h_L_CaMK << ' ' << a << ' ' <<
    i_fast << ' ' << i_slow << ' ' << a_CaMK << ' ' << i_CaMK_fast << ' ' << i_CaMK_slow << ' ' << d << ' ' << f_fast <<
    ' ' << f_slow << ' ' << f_Ca_fast << ' ' << f_Ca_slow << ' ' <<
    j_Ca << ' ' << f_CaMK_fast << ' ' << f_Ca_CaMK_fast << ' ' << n_ss << ' ' << n_i << ' ' << C_0 << ' ' << C_1 << ' ' << C_2 << ' ' << O << ' ' << I << ' ' <<
    x_s1 << ' ' << x_s2 <<' ' <<
    Na_i << ' ' << Na_ss << ' ' << K_i << ' ' << K_ss << ' ' << Ca_i << ' ' << Ca_ss << ' ' << Ca_nsr << ' ' <<
    Ca_jsr << ' ' << Cl_i << ' ' << CaMK_trap << ' ' << J_rel_NP << ' ' << J_rel_CaMK << ' '
  ;
}

void TWorld::LongPrint(ostream &tempstr, double tArg,  ML_CalcType V) {
  Print(tempstr, tArg, V);

  tempstr << I_Na << ' ' << I_Na_late << ' ' << I_to << ' ' << I_CaL_i << ' ' << I_CaL_ss << ' ' << I_CaNa_i << ' ' << I_CaNa_ss << ' ' << I_CaK_i << ' ' << I_CaK_ss << ' ' << I_CaL << ' ' << I_CaNa << ' ' << I_CaK << ' ' <<
    I_Kr << ' ' << I_Ks << ' ' << I_K1 << ' ' << I_NaCa_i << ' ' << I_NaCa_ss << ' ' << I_NaK << ' ' << I_CaCl << ' ' << I_Nab << ' ' <<
    I_Cab << ' ' << I_Kb << ' ' << I_Clb << ' '<< I_pCa << ' '<< J_diff_Na << ' '<< J_diff_Ca << ' '<< J_diff_K << ' '<< J_leak << ' '<< J_rel << ' '<< J_tr << ' '<< J_up << ' '<< cur_I_tot << ' '
  ;
}  // TWorld::LongPrint

void TWorld::GetParameterNames(vector<string> &getpara) {
  const string ParaNames[] =
  {
    "m", "A_h", "B_h",          "h",  "A_j", "B_j",                "j",                            "h_p",
    "j_p",
    "m_L",         "h_L",
    "h_L_CaMK",    "a",                  "i_fast",                  "i_slow",                       "a_CaMK",
    "i_CaMK_fast", "i_CaMK_slow",        "d",                       "f_fast",                       "f_slow",
    "f_Ca_fast",
    "f_Ca_slow",
    "j_Ca",        "f_CaMK_fast",        "f_Ca_CaMK_fast",          "n_ss",          "n_i", "C_0", "C_1", "C_2", "O", "I",
    "x_s1",
    "x_s2",               "Na_i",                    "Na_ss",                        "K_i",
    "K_ss",
    "Ca_i",        "Ca_ss",
    "Ca_nsr",      "Ca_jsr", "Cl_i",             "CaMK_trap",               "J_rel_NP",
    "J_rel_CaMK"
  };

  for (int i = 0; i < sizeof(ParaNames)/sizeof(ParaNames[0]); i++)
    getpara.push_back(ParaNames[i]);
}  // TWorld::GetParameterNames

void TWorld::GetLongParameterNames(vector<string> &getpara) {
  GetParameterNames(getpara);
  const string ParaNames[] =
  {
    "I_Na", "I_Na_late", "I_to", "I_CaL_i", "I_CaL_ss", "I_CaNa_i", "I_CaNa_ss", "I_CaK_i", "I_CaK_ss",      "I_CaL",      "I_CaNa",      "I_CaK",      "I_Kr",      "I_Ks",
    "I_K1",
    "I_NaCa_i",
    "I_NaCa_ss",
    "I_NaK",     "I_CaCl",     "I_Nab",     "I_Cab",     "I_Kb",     "I_Clb",       "I_pCa", "J_diff_Na", "J_diff_Ca", "J_diff_K", "J_leak", "J_rel", "J_tr", "J_up", "cur_I_tot"
  };
  for (int i = 0; i < sizeof(ParaNames)/sizeof(ParaNames[0]); i++)
    getpara.push_back(ParaNames[i]);
}
