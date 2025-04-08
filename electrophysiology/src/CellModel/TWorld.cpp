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

  m              = v(VT_Init_m);
  A_h            = v(VT_Init_A_h);
  B_h            = v(VT_Init_B_h);
  h              = v(VT_Init_h);
  j              = v(VT_Init_j);
  h_p            = v(VT_Init_h_p);
  j_p            = v(VT_Init_j_p);
  m_L            = v(VT_Init_m_L);
  h_L            = v(VT_Init_h_L);
  h_L_CaMK       = v(VT_Init_h_L_CaMK);
  a              = v(VT_Init_a);
  i_fast         = v(VT_Init_i_fast);
  i_slow         = v(VT_Init_i_slow);
  a_CaMK         = v(VT_Init_a_CaMK);
  i_CaMK_fast    = v(VT_Init_i_CaMK_fast);
  i_CaMK_slow    = v(VT_Init_i_CaMK_slow);
  d              = v(VT_Init_d);
  f_fast         = v(VT_Init_f_fast);
  f_slow         = v(VT_Init_f_slow);
  f_Ca_fast      = v(VT_Init_f_Ca_fast);
  f_Ca_slow      = v(VT_Init_f_Ca_slow);
  j_Ca           = v(VT_Init_j_Ca);
  n_ss           = v(VT_Init_n_ss);
  n_i            = v(VT_Init_n_i);
  f_CaMK_fast    = v(VT_Init_f_CaMK_fast);
  f_Ca_CaMK_fast = v(VT_Init_f_Ca_CaMK_fast);
  C_0            = v(VT_Init_C_0);
  C_1            = v(VT_Init_C_1);
  C_2            = v(VT_Init_C_2);
  O              = v(VT_Init_O);
  I              = v(VT_Init_I);
  x_s1           = v(VT_Init_x_s1);
  x_s2           = v(VT_Init_x_s2);
  Na_i           = v(VT_Init_Na_i);
  Na_ss          = v(VT_Init_Na_ss);
  K_i            = v(VT_Init_K_i);
  K_ss           = v(VT_Init_K_ss);
  Ca_i           = v(VT_Init_Ca_i);
  Ca_ss          = v(VT_Init_Ca_ss);
  Ca_nsr         = v(VT_Init_Ca_nsr);
  Ca_jsr         = v(VT_Init_Ca_jsr);
  Cl_i           = v(VT_Init_Cl_i);
  CaMK_trap      = v(VT_Init_CaMK_trap);
  J_rel_NP       = v(VT_Init_J_rel_NP);
  J_rel_CaMK     = v(VT_Init_J_rel_CaMK);
  #ifdef TRPN
  Ca_TRPN = v(VT_Init_Ca_TRPN)/0.07;
  #endif  // ifdef TRPN
}  // TWorld::Init

ML_CalcType TWorld::Calc(double tinc,  ML_CalcType V,  ML_CalcType i_external,  ML_CalcType stretch, int euler) {
  tinc *= 1000.0;  // second to millisecond conversion
  ML_CalcType V_m = V * 1000.0;
  
  #if defined(TRPN) || defined(ISAC)
  ML_CalcType lambda_m = std::min(1.2, stretch);
  #endif  // if defined(TRPN) || defined(ISAC)

  const int Vi = (int)(DivisionTab*(RangeTabhalf+V_m)+.5);  // array position

  /// CaMK constants
  double KmCaMK = 0.15;
  double aCaMK  = 0.05;
  double bCaMK  = 0.00068;
  double CaMKo  = 0.05;
  double KmCaM  = 0.0015;

  /// update CaMK
  const ML_CalcType CaMK_bound  = ((CaMKo * (1.0 - CaMK_trap)) / (1.0 + (KmCaM / Ca_ss)));
  const ML_CalcType CaMK_active = (CaMK_bound + CaMK_trap);
  const ML_CalcType dCaMK_trap = ((aCaMK * CaMK_bound * (CaMK_bound + CaMK_trap)) -
    (bCaMK * CaMK_trap));
  CaMK_trap += tinc * dCaMK_trap;

  /// reversal potentials
  const ML_CalcType E_Na = v(VT_RToverF) * log(v(VT_Na_o) / Na_i);
  const ML_CalcType E_K  = v(VT_RToverF) * log(v(VT_K_o) / K_i);
  double PKNa            = 0.01833;
  const ML_CalcType E_Ks = v(VT_RToverF) * log((v(VT_K_o) + PKNa * v(VT_Na_o)) / (K_i + PKNa * Na_i));
  double zcl = 1.0;
  const ML_CalcType E_Cl = ((v(VT_R) * v(VT_T)) / (zcl * v(VT_F))) * log((Cl_i/v(VT_Cl_o)));

  // convenient shorthand calculations
  const ML_CalcType VFFoverRT = (V_m * v(VT_F) * v(VT_F)) / (v(VT_R) * v(VT_T));
  const ML_CalcType VFoverRT  = (V_m * v(VT_F)) / (v(VT_R) * v(VT_T));
  const ML_CalcType util_1    = v(VT_A_cap) / (v(VT_F) * v(VT_v_myo));
  const ML_CalcType util_2    = v(VT_v_ss) / v(VT_v_myo);
  const ML_CalcType util_3    = v(VT_A_cap) / (v(VT_F) * v(VT_v_ss));
    
  /////////////////////////////////////////////////////////////////////////////////////////
  ///        Sodium current (INa, INaL)
  ////////////////////////////////////////////////////////////////////////////////////////
  
  ///////////// calulate I_Na //////////
  // m gate
  const ML_CalcType m_inf = 1.0 / ((1.0 + exp(-(V_m + 56.86)/9.03))*(1.0 + exp(-(V_m + 56.86)/9.03))); //passt
  const ML_CalcType tau_m = 0.1292 * exp(-(((V_m + 45.79)/15.54) * ((V_m + 45.79)/15.54))) + 0.06487 * exp(-(((V_m - 4.823)/51.12) * ((V_m - 4.823)/51.12))); //passt
  m = m_inf - (m_inf - m) * exp(-tinc / tau_m); //passt

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
  // TODO: h scheint im MATLAB Code komplett zu fehlen??? hier von Tomek sollte eigentlich passen
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
      
  const ML_CalcType phi_INa_CaMK = ; //TODO: Formel unklar
  const ML_CalcType phi_INa_PKA = fINaP; //TODO: Wat is dat??? // PKA-P fraction as assigned as input, take the value 0 or 1
  const ML_CalcType phi_INa_Both = phi_INa_CaMK * phi_INa_PKA;
  const ML_CalcType phi_INa_CaMKonly = phi_INa_CaMK - phi_INa_Both;
  const ML_CalcType phi_INa_PKAonly = phi_INa_PKA - phi_INa_Both;

  I_Na_Base_NP = G_Na * m*m*m * h * j; //Non-Phosphorylated
  I_Na_Base_CaMK = G_Na * m*m*m * h_p * j_p;
  I_Na_Base_PKA = G_Na_PKA * m_PKA*m_PKA*m_PKA * h_PKA * j_PKA;
  I_Na_Base_Both = G_Na_PKA * m_PKA*m_PKA*m_PKA * h_both * j_both;
      
  // 4 population
  I_Na_Base = ((1-phi_INa_CaMKonly-phi_INa_PKAonly-phi_INa_Both)*I_Na_Base_NP + phi_INa_CaMKonly*I_Na_Base_CaMK + phi_INa_PKAonly*I_Na_Base_PKA + phi_INa_Both*I_Na_Base_Both);
  I_NaFast_junc = Fjunc * I_Na_Base*(V_m - E_Na_junc); // TODO: Was ist Fjunc und Nernst junc
  I_NaFast_sl = (1 - Fjunc) * I_Na_Base * (V_m - E_Na_sl); //TODO: Was is Fjunc und Nerst sl

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
    const ML_CalcType phi_ICaL_CaMK = ; //TODO: unklar was genau hier hin
    const ML_CalcType phi_ICaL_PKA = ; //TODO: unklar was genau hier hin
  
    // d gate
    const ML_CalcType d_inf = min(1.0763*exp((-1.007*exp((-0.0829*(V_m+3.62483))))), 1.0);
    const ML_CalcType tau_d = 1.5 + 1.0 / (exp(-0.05 * (V_m + 6.0)) + exp(0.09 * (V_m + 14.0)));
    d = d_inf - (d_inf -d) * exp(-tinc / tau_d);
    
    // f gate
    const ML_CalcType f_inf = 1.0 / (1.0 + exp((V_m + 19.58) / 3.696));
    const ML_CalcType tau_f_fast = 6.17111 + 1.0 / (0.00126 * exp(-(V_m + 26.63596) / 9.69961) + 0.00126 * exp((V_m + 26.63596) / 9.69961));
    const ML_CalcType tau_f_slow = 2719.22489 + 1.0 / (7.19411e-05 * exp(-(V_m + 5.74631) / 10.87690) + 7.19411e-05 * exp((V_m + 5.74631) / 16.31535));
    
    const ML_CalcType A_f_fast = 0.52477; //TODO: can be replaced by using v()
    const ML_CalcType A_f_slow = 1.0 - A_f_fast; //TODO: can be replaced by using v()
    f_fast = f_inf - (f_inf - f_fast) * exp(-tinc / tau_f_fast);
    f_slow = f_inf - (f_inf - f_slow) * exp(-tinc / tau_f_slow);
    const ML_CalcType f = A_f_fast * f_fast + A_f_slow * f_slow;
   //TODO: other option -> const ML_CalcType f = v(VT_A_f_fast) * f_fast + v(VT_A_f_slow) * f_slow
    
    const ML_CalcType f_Ca_inf      = f_inf;
    const ML_CalcType tau_f_Ca_fast = 13.50673 + 1.0 / (0.15420 * exp(-(V_m - 1.31611) / 11.33960) + 0.15420 * exp((V_m - 1.31611) / 11.33960));
    const ML_CalcType tau_f_Ca_slow = 177.95813 + 1.0 / (4.73955e-04 * exp((-V_m+0.79049) / 0.81777) + 4.73955e-04 * exp((V_m+2.40474) /1.90812));
    const ML_CalcType A_f_Ca_fast = 0.3 + 0.6 / (1.0 + exp((V_m - 9.24247) / 27.96201));
    const ML_CalcType A_f_Ca_slow = 1.0 - A_f_Ca_fast;
    f_Ca_fast = f_Ca_inf - (f_Ca_inf - f_Ca_fast) * exp(-tinc / tau_f_Ca_fast);
    f_Ca_slow = f_Ca_inf - (f_Ca_inf - f_Ca_slow) * exp(-tinc / tau_f_Ca_slow);
    const ML_CalcType f_Ca        = A_f_Ca_fast * f_Ca_fast + A_f_Ca_slow * f_Ca_slow;
    
    // j gate
    double tau_j_Ca = 66.0;
    const ML_CalcType j_Ca_inf = 1. / (1 + exp((V_m + 17.66945) / (3.21501)));
    j_Ca = j_Ca_inf - (j_Ca_inf - j_Ca) * exp(-tinc / tau_j_Ca);
    
    // gating CaMK-P
    const ML_CalcType tau_f_p_fast = 2.5 * tau_f_fast;
    const ML_CalcType f_p_inf = f_inf;
    f_p_fast = f_p_inf - (f_p_inf - f_p_fast) * exp(-tinc / tau_f_p_fast);
    const ML_CalcType f_p = A_f_fast * f_p_fast + A_f_slow * f_slow;
    //TODO: Alternative -> const ML_CalcType f_p = v(VT_A_f_fast) * f_p_fast + v(VT_A_f_slow) * f_slow;
    
    const ML_CalcType tau_f_Ca_p_fast = 2.5 * tau_f_Ca_fast;
    const ML_CalcType f_Ca_p_inf = f_inf;
    f_Ca_p_fast = f_Ca_p_inf - (f_Ca_p_inf - f_Ca_p_fast) * exp(-tinc / tau_f_Ca_p_fast);
    const ML_CalcType f_Ca_p = A_f_Ca_fast * f_Ca_p_fast + A_f_Ca_slow * f_Ca_slow;
    
    // SS nca
    double Kmn = 0.00222;
    double k2n = 957.85903;
    const ML_CalcType k_m2_n = j_Ca * 0.84191;
    const ML_CalcType alpha_n_Ca_junc   = 1.0 / ((k2n / k_m2_n) + pow((1.0 + (Kmn / Ca_junc)), 3.80763)); //TODO: Ca_junc does not exist yet
    n_junc = alpha_n_Ca_junc * (k2n / k_m2_n) - (alpha_n_Ca_junc * (k2n / k_m2_n) - n_ss) * exp(-k_m2_n * tinc); //TODO: Ist noch falsch
    
    // myoplasmic nca
    const ML_CalcType alpha_n_Ca_sl   = 1.0 / ((k2n / k_m2_n) + pow((1.0 + (Kmn / Ca_sl)), 3.80763)); //TODO: Ca_sl does not exist yet
    n_sl = alpha_n_Ca_sl * (k2n / k_m2_n) - (alpha_n_Ca_junc * (k2n / k_m2_n) - n_ss) * exp(-k_m2_n * tinc); //TODO: Ist noch falsch
    
    // SS driving force
    
    
//    /// calculate I_CaL, I_CaNa, I_CaK

//

//
//    const ML_CalcType I_o = (0.5*(v(VT_Na_o)+v(VT_K_o)+v(VT_Cl_o)+(4.*v(VT_Ca_o)))/1000.);
//    const ML_CalcType I_ss = ((0.5*(Na_ss+K_ss+Cl_i+(4.*Ca_ss)))/1000.);
//    double constA = (1.82e6*pow((74.*310.),-1.5)); // Diel constant and temperature as constants
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
    
    
//    /// Calculate I_CaCl
//    I_CaCl_junc = (0.2843 * (V_m - E_Cl)) / (1. + (0.1/Ca_ss));
//    I_CaCl_sl = (0 * 0.2843 * (V_m - E_Cl)) / (1. + (0.1/Ca_i));
//    I_CaCl = I_CaCl_junc + I_CaCl_sl;
//    
//    /// calculate I_Clb
//    double P_Clb = 1.98e-3;
//    I_Clb = P_Clb * (V_m - E_Cl);
    
  //////////////////////// Background currents (INab, ICab, IKb) ////////////////////////
    
//    /// calculate I_Kb
//    const ML_CalcType x_Kb = 1.0 / (1.0 + exp(-(V_m - 10.8968) / (23.9871)));
//    double G_Kb            = 0.0189;
//    I_Kb = v(VT_IKb_Multiplier) * G_Kb * x_Kb * (V_m - E_K);
//
//    /// calculate I_Nab
//    double P_Nab = 1.9239e-9;
//    I_Nab = v(VT_INab_Multiplier) * P_Nab * z_Na * z_Na * VFFoverRT *
//      ((Na_i * exp_z_Na - v(VT_Na_o)) / (exp_z_Na - 1.0));
//
//    /// calculate I_Cab
//    double P_Cab = 5.9194e-8;
//    I_Cab = v(VT_ICab_Multiplier) * P_Cab * z_Ca * z_Ca * VFFoverRT *
//      ((gamma_Ca_i * Ca_i * exp_z_Ca - gamma_Ca_o * v(VT_Ca_o)) / (exp_z_Ca - 1.0));
    
  //////////////////////// Calcium release from SR (Jrel, Jleak) ////////////////////////
    
  //////////////////////// Calcium reuptake to the SR (Jup) ////////////////////////
    
  //////////////////////// Sarcolemmal calcium pump (pCa) ////////////////////////
    
//    /// calculate I_pCa
//    double G_pCa = 0.0005;
//    I_pCa = v(VT_IpCa_Multiplier) * ((G_pCa * Ca_i) / (0.0005 + Ca_i));
    

















  /// stretch activated channel
 #ifdef ISAC
  I_SAC = v(VT_G_sac) * (V_m - v(VT_E_sac)) / (1. + v(VT_K_sac) * exp(-v(VT_alpha) * (stretch - 1.)));
  I_SAC *= v(VT_ISAC_SWITCH);
 #endif  // ifdef ISAC

  /// I_tot
   cur_I_tot =
    I_Na + I_Na_late + I_to + I_CaL + I_CaNa + I_CaK + I_Kr + I_Ks + I_K1 + I_NaCa_i + I_NaCa_ss + I_NaK +
    I_Nab + I_Cab + I_Kb + I_pCa + I_Clb + I_CaCl 
  #ifdef ISAC
    + I_SAC
  #endif  // ifdef ISAC
    + i_external;
  I_tot = -cur_I_tot;

  /// diffusion fluxes J_diff_Na, J_diff_Ca, J_diff_K
  double tau_diff_Na          = 2.0;
  double tau_diff_Ca          = 0.2;
  double tau_diff_K           = 2.0;
  J_diff_Na = (Na_ss - Na_i) / tau_diff_Na;
  J_diff_Ca = (Ca_ss - Ca_i) / tau_diff_Ca;
  J_diff_K  = (K_ss - K_i) / tau_diff_K;

  /// SR Calcuim release flux J_rel
  double K_inf_rel = 1.7;
  ML_CalcType J_rel_NP_inf = (v(VT_alpha_rel) * (-I_CaL)) / (1.0 + pow((K_inf_rel/ Ca_jsr), 8.0));
  ML_CalcType J_rel_CaMK_inf = (v(VT_alpha_rel_CaMK) * (-I_CaL)) / (1.0 + pow((K_inf_rel/ Ca_jsr), 8.0));
  if (v(VT_celltype) == 2.0) {
    J_rel_NP_inf   *= 1.7;
    J_rel_CaMK_inf *= 1.7;
  }
  const ML_CalcType tau_rel_NP_b = v(VT_beta_tau) / (1.0 + (0.0123 / Ca_jsr));
  const ML_CalcType tau_rel_NP   = tau_rel_NP_b < 0.001 ? 0.001 : tau_rel_NP_b;
  J_rel_NP = J_rel_NP_inf - (J_rel_NP_inf - J_rel_NP) * exp(-tinc / tau_rel_NP);
  
  const ML_CalcType tau_rel_CaMK_b = v(VT_beta_tau_CaMK) / (1.0 + (0.0123 / Ca_jsr));
  const ML_CalcType tau_rel_CaMK   = tau_rel_CaMK_b < 0.001 ? 0.001 : tau_rel_CaMK_b;
  J_rel_CaMK = J_rel_CaMK_inf - (J_rel_CaMK_inf - J_rel_CaMK) * exp(-tinc / tau_rel_CaMK);
  
  const ML_CalcType phi_rel_CaMK = phi_INa_CaMK;
  J_rel        = 1.5378 * (1.0 - phi_rel_CaMK) * J_rel_NP + phi_rel_CaMK * J_rel_CaMK;

  /// J_up
  ML_CalcType J_up_NP   = (0.005425 * Ca_i) / (0.00092 + Ca_i);
  ML_CalcType J_up_CaMK = (2.75 * 0.005425 * Ca_i) / (0.00092 - 0.00017 + Ca_i);
  if (v(VT_celltype) == 1.0) {
    J_up_NP   *= 1.3;
    J_up_CaMK *= 1.3;
  }
  const ML_CalcType phi_up_CaMK = (1.0 /(1.0 + (KmCaMK / CaMK_active)));
  J_leak      = (0.0048825 * Ca_nsr) / 15.0;

  J_up = ((1.0 - phi_up_CaMK) * J_up_NP + (phi_up_CaMK * J_up_CaMK) - J_leak);


  /// J_tr
  double tau_tr          = 60.0;
  J_tr = (Ca_nsr - Ca_jsr) / tau_tr;

  /// Concentrations
  const ML_CalcType cur_Nai = I_Na + I_Na_late + 3.0 * I_NaCa_i + I_CaNa_i + 3.0 * I_NaK + I_Nab;
  const ML_CalcType dNa_i = -(cur_Nai) * util_1 + J_diff_Na * util_2;
  Na_i += tinc * dNa_i;
  const ML_CalcType dNa_ss = -(I_CaNa_ss + 3.0 * I_NaCa_ss) * util_3 - J_diff_Na;
  Na_ss += tinc * dNa_ss;

  //const ML_CalcType cur_Ki = I_to + I_Kr + I_Ks + I_K1 + I_Kb - 2.0 * I_NaK + I_CaK_i;
  const ML_CalcType cur_Ki = I_to + I_Kr + I_Ks + I_K1 + I_Kb + i_external - 2.0 * I_NaK + I_CaK_i;
  const ML_CalcType dK_i = -(cur_Ki) * util_1 + J_diff_K * util_2;
  K_i += tinc * dK_i;
  const ML_CalcType dK_ss = -(I_CaK_ss) * util_3 - J_diff_K;
  K_ss += tinc * dK_ss;

  /// Ca buffer constants
  double cmdnmax = 0.05;
  if (v(VT_celltype) == 1) {
    cmdnmax *= 1.3;
  }
  double kmcmdn  = 0.00238;
  double trpnmax = 0.07;
  double kmtrpn  = 0.0005;
  double BSRmax  = 0.047;
  double KmBSR   = 0.00087;
  double BSLmax  = 1.124;
  double KmBSL   = 0.0087;
  double csqnmax = 10.0;
  double kmcsqn  = 0.8;
  #ifdef TRPN
  const ML_CalcType Ca_T50 = (v(VT_Ca_T50) + v(VT_beta1) * (lambda_m - 1.0));
  const ML_CalcType dCaTRPN = v(VT_k_TRPN)*(pow(1000*Ca_i/Ca_T50, v(VT_n_TRPN))*(1.-Ca_TRPN)-Ca_TRPN);// needs muM with Land parameters
  Ca_TRPN += tinc * dCaTRPN;
  const ML_CalcType I_Trpn = trpnmax * dCaTRPN;
  #endif  // ifdef TRPN
  
  const ML_CalcType beta_Cai = 1.0 / (1.0 + ((cmdnmax * kmcmdn) / pow((kmcmdn + Ca_i), 2.0))
  #ifndef TRPN
                             + ((trpnmax * kmtrpn) / pow((kmtrpn + Ca_i), 2.0))
  #endif  // ifndef TRPN
                                      );
  const ML_CalcType cur_Cai = I_CaL_i + I_pCa + I_Cab - 2.0 * I_NaCa_i;
  

  //const ML_CalcType dCa_i = beta_Cai * (-(cur_Cai) * 0.5 * util_1 - J_up * (v(VT_v_nsr) / v(VT_v_myo)) + J_diff_Ca * util_2
  const ML_CalcType dCa_i = beta_Cai * (-(cur_Cai) * 0.5 * util_1 - J_up * (v(VT_v_nsr) / v(VT_v_myo)) + J_diff_Ca * util_2
  #ifdef TRPN
      - I_Trpn
  #endif  // ifdef TRPN
      );
  Ca_i += tinc * dCa_i;
  
  const ML_CalcType beta_Cass = 1.0 /
    (1.0 +
     ((BSRmax * KmBSR) /
      pow((KmBSR + Ca_ss), 2.0)) + ((BSLmax * KmBSL) / pow((KmBSL + Ca_ss), 2.0)));
  const ML_CalcType dCa_ss = beta_Cass *
    (-(I_CaL_ss - 2.0 * I_NaCa_ss) * 0.5 * util_3 + J_rel * (v(VT_v_jsr) / v(VT_v_ss)) - J_diff_Ca);
  Ca_ss += tinc * dCa_ss;
  const ML_CalcType dCa_nsr = J_up - ((J_tr * v(VT_v_jsr)) / v(VT_v_nsr));
  Ca_nsr += tinc * dCa_nsr;
  const ML_CalcType beta_Cajsr = 1.0 / (1.0 + ((csqnmax * kmcsqn) / pow((kmcsqn + Ca_jsr), 2.0)));
  const ML_CalcType dCa_jsr    = beta_Cajsr * (J_tr - J_rel);
  Ca_jsr += tinc * dCa_jsr;


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
    #ifdef TRPN
    << 0.07*Ca_TRPN << ' ' // export in mM
    #endif  // ifdef TRPN
  ;
}

void TWorld::LongPrint(ostream &tempstr, double tArg,  ML_CalcType V) {
  Print(tempstr, tArg, V);

  tempstr << I_Na << ' ' << I_Na_late << ' ' << I_to << ' ' << I_CaL_i << ' ' << I_CaL_ss << ' ' << I_CaNa_i << ' ' << I_CaNa_ss << ' ' << I_CaK_i << ' ' << I_CaK_ss << ' ' << I_CaL << ' ' << I_CaNa << ' ' << I_CaK << ' ' <<
    I_Kr << ' ' << I_Ks << ' ' << I_K1 << ' ' << I_NaCa_i << ' ' << I_NaCa_ss << ' ' << I_NaK << ' ' << I_CaCl << ' ' << I_Nab << ' ' <<
    I_Cab << ' ' << I_Kb << ' ' << I_Clb << ' '<< I_pCa << ' '<< J_diff_Na << ' '<< J_diff_Ca << ' '<< J_diff_K << ' '<< J_leak << ' '<< J_rel << ' '<< J_tr << ' '<< J_up << ' '<< cur_I_tot << ' '
  #ifdef ISAC
    << I_SAC << ' '
  #endif // ifdef ISAC
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
   #ifdef TRPN
    ,              "Ca_TRPN"
   #endif  // ifdef TRPN
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
#ifdef ISAC
    , "I_SAC"
#endif // ifdef ISAC
  };
  for (int i = 0; i < sizeof(ParaNames)/sizeof(ParaNames[0]); i++)
    getpara.push_back(ParaNames[i]);
}
