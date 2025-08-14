#ifndef __Spectral__
#define __Spectral__

#include "numerics.h"
#include "highPrec.h"
#include "stat.h"

using namespace std;

class CHLT {

public:
  CHLT()  {syst=0;}
  
 ////////////////////////////////////

  PrecVect g, g_2;
  distr_t rho;
  double syst;
};


void Get_bt(PrecVect &bt, const PrecFloat E, const function<PrecFloat(int, PrecFloat)> &Bn,const int tmin,const int tmax);  

void Get_Atr(PrecMatr& Atr, int tmin, int tmax,const function<PrecFloat(int, int)> &Atr_lambda ) ;

void Get_ft(PrecVect& ft, int tmin, int tmax,  const function<PrecFloat(int)> &func_f);


void Compute_covariance_matrix(PrecMatr &B, int tmin, int tmax, Vfloat &covariance, const distr_t_list &corr);
  
void Get_optimal_lambda(const PrecMatr &Atr, const PrecMatr &B,const PrecVect &ft, const PrecFloat & M2, double& lambda_opt, double& lambda_opt_2,  const distr_t_list & corr, int tmin, int tmax,const double mult, double mult2, double Ag_ov_A0_tg,  string path );

CHLT Get_INVLT(int tmin, int tmax, const PrecFloat M2, const function<PrecFloat(int,int)> &Atr_lambda, const function<PrecFloat(int)> &func_f,  Vfloat &covariance, const distr_t_list& corr, double mult, double mult2, double Ag_ov_A0_tg,  string out_path, bool INCLUDE_ERRORS, int prec);


#endif
