#include "../include/Spectral.h"
#include "highPrec.h"

const int verbosity_lev=1;
const bool DEACTIVATE_FULL_COV=true;


using namespace std;


void Get_bt(PrecVect &bt, const PrecFloat E, const function<PrecFloat(int, PrecFloat)> &Bn,const int tmin,const int tmax) {

  bt.resize(tmax-tmin+1);

  for(int t=tmin;t<=tmax;t++) bt(t-tmin) = Bn(t, E);

};



void Get_Atr(PrecMatr& Atr, int tmin, int tmax,const function<PrecFloat( int, int)> &Atr_lambda )  {

  Atr.resize(tmax-tmin+1, tmax-tmin+1);
#pragma omp parallel for
  for(int t=tmin;t<= tmax; t++) {
    for(int r=t; r<= tmax; r++)
      {
	Atr(t-tmin, r-tmin) = Atr_lambda( t,r);
	Atr(r-tmin, t-tmin) = Atr(t-tmin,r-tmin);
      }
  }
  return;
}


void Get_ft(PrecVect& ft, int tmin, int tmax,  const function<PrecFloat(int)> &func_f) {

  ft.resize(tmax-tmin+1);
#pragma omp parallel for
  for(int t=tmin;t<=tmax;t++) {

    if(verbosity_lev>=2) cout<<"Computing time: "<<t<<"..."<<flush;

    ft(t-tmin) = func_f(t);
  }     
  return;
}


void Compute_covariance_matrix(PrecMatr &B, int tmin, int tmax, Vfloat &covariance, const distr_t_list &corr) {

  
  PrecFloat norm= sqr(PrecFloat(corr.ave(0)));

  
  B.resize(tmax-tmin+1, tmax-tmin+1);


  for(int t=tmin;t<=tmax;t++) {
    for(int r=t;r<=tmax;r++) {
      
	B(t-tmin,r-tmin) = covariance[t*corr.size()+ r];
	B(r-tmin, t-tmin) = B(t-tmin,r-tmin);
	if(r != t && DEACTIVATE_FULL_COV) { B(t-tmin,r-tmin) = 0; B(r-tmin,t-tmin) = 0; }
    }
  }

  B= (1.0/norm)*B;
 
  return;
}

void automated_plateaux_search(const PrecMatr &Atr, const PrecMatr &Btr,const PrecVect &ft, const PrecFloat & M2, double& lambda_opt, double& lambda_opt_2, double &ch2_ret,  const distr_t_list & corr, int tmin, int tmax, string path ) {

  bool UseJack= corr.UseJack;

  ofstream Print_R_at_lambda(path);

  //Print header
  Print_R_at_lambda<<"# $1=lambda,  $2= A[g]/A[0], $3=B[g], $4= val, $5=err, $6=S_FLAG $7=mult"<<endl; 
  
  Print_R_at_lambda.precision(10);
    
  if(verbosity_lev) cout<<"Automated plateaux search ..."<<flush;

 
  const auto A=
    [&M2, &Atr, &ft](const PrecVect& gmin) -> PrecFloat
    {
      PrecVect Atr_g = Atr*gmin;
      PrecFloat g_Atr_g = gmin.transpose()*Atr_g;
      PrecFloat ft_g = ft.transpose()*gmin;
      return (M2 + g_Atr_g -2*ft_g)/M2  ;
    };

  const auto B=
    [&Btr](const PrecVect& gmin) -> PrecFloat
    {
      PrecVect B_g = Btr*gmin;
      PrecFloat g_B_g = gmin.transpose()*B_g;
      return g_B_g ;
    };


  double r= 0.7; //0.5;
  //compute minimum possible A[g]/A[0]
  PrecMatr C_min = Atr/M2;
  PrecMatr C_inv_min = C_min.inverse();
  PrecVect ft_l_min = ft/M2;
  PrecVect gm_min = C_inv_min*ft_l_min;     
  PrecFloat A_val_min = A(gm_min);
  int Npoints= 36;
  double As=0.3;
  int Npoints_max =  (int)(  (log(A_val_min.get()) - log(As))/log(r) + 1 );
  Npoints = min(36, Npoints_max);
  Vfloat Ags;
  vector<PrecFloat> lambdas;
  vector<PrecFloat> A_tofit;
  vector<PrecFloat> B_tofit;
  vector<distr_t> R_tofit;
  for(int i=0;i<Npoints;i++) Ags.push_back( As*pow(r,i));  
  
  bool JUST_STARTED=true;
  PrecFloat l_start=1;
  PrecFloat l_low =0;
  PrecFloat l_up = 1;
  int Nit_Ag0=0;
  
  for(auto & Ag_T: Ags) {

    
    bool lambda_found_Ag_A0=false;
    Nit_Ag0=0;
    l_low=0;
    
    //bisection search for given A[g]/A[0]
    while( !lambda_found_Ag_A0 ) {

     
      PrecFloat lambda_mid  =  (Nit_Ag0==0 && JUST_STARTED)?l_start:(l_up+l_low)/2;
      PrecMatr C = Atr*(1-lambda_mid)/M2 + Btr*lambda_mid;
      PrecMatr C_inv = C.inverse();
      PrecVect ft_l = ft*(1-lambda_mid)/M2;
      PrecVect gm = C_inv*ft_l;

      
      PrecFloat A_val = A(gm);
      PrecFloat B_val = B(gm);
      PrecFloat W_val = (1-lambda_mid)*A_val + lambda_mid*B_val;
      
      
      cout.precision(10);

      if(A_val > Ag_T) { // lambda_mid is new l_low
	l_up =lambda_mid;
      }
      else { //lambda_mid is new l_up
	l_low = lambda_mid;
      }

      
      
      Nit_Ag0++;
      if( (A_val > 0.95*Ag_T && A_val < 1.05*Ag_T)) lambda_found_Ag_A0=true;


      //##########################################################################################
      //compute anti-Laplace transform corresponding to lambda_mid:
      distr_t R_E_lambda(UseJack);
      
      for(int ijack=0; ijack< corr.distr_list[1].size(); ijack++) {

	PrecFloat spec_lambda_d_jack=0;

	for(int t=tmin;t<=tmax;t++) spec_lambda_d_jack = spec_lambda_d_jack + gm(t-tmin)*corr.distr_list[t].distr[ijack]; 

	R_E_lambda.distr.push_back( spec_lambda_d_jack.get());
      }

      PrecFloat mult_est = A_val/B_val;


      //cout<<"A: "<<A_val<<endl<<"lambda: "<<lambda_mid<<endl<<flush;
     
      Print_R_at_lambda<<lambda_mid<<" "<<A_val<<" "<<B_val<<" "<<R_E_lambda.ave()<<" "<<R_E_lambda.err();
      Print_R_at_lambda<<" "<<lambda_found_Ag_A0<<" "<<mult_est<<endl;

      //##########################################################################################

      if(lambda_found_Ag_A0) {  lambdas.push_back( lambda_mid) ; R_tofit.push_back(R_E_lambda); A_tofit.push_back(A_val); B_tofit.push_back(B_val); }

     
    }
    JUST_STARTED=false;
  }

  

  //we can now find the interval where ch2/dof is O(1)
  int Npoints_fit = 7; //must be and odd number
  assert( (Npoints_fit%2) != 0);
  assert(Npoints_fit < Npoints);
  int fit_start=0;
  bool found_interval=false;
  double ch2_thresh=1;

  auto constant_fit = [&Npoints_fit](const vector<distr_t> &A ) -> distr_t {
    int N= A.size();
    assert(N==Npoints_fit);
    double tot_w=0;
    distr_t ret = 0.0*A[0];
    for(int i=0;i<N;i++) { tot_w += 1.0/pow(A[i].err(),2); ret = ret + A[i]/pow(A[i].err(),2); }
    return ret/tot_w;
  };

  auto ch2_func = [&Npoints_fit](const vector<distr_t> &A,const distr_t &res) -> double {
    int N=A.size();
    assert(N==Npoints_fit);
    double ch2=0;
    for(int i=0;i<N;i++) ch2 += pow( (A[i].ave() - res.ave())/A[i].err(),2);
    return ch2;
  };
  
  while(!found_interval) {

    vector<distr_t> R;
    
    for(int ip=0;ip<Npoints_fit;ip++) R.push_back( R_tofit[fit_start+ip]);
    if( ch2_func(R, constant_fit(R))/4 < ch2_thresh )   found_interval=true;
    if(!found_interval) fit_start++;
    if( (fit_start >= ((signed)R_tofit.size() - Npoints_fit +1)) && !found_interval) { fit_start= 0; ch2_thresh += 0.5; }
  }

  int id_lambda= fit_start + Npoints_fit/2;
  int id_lambda2 = fit_start + Npoints_fit -1;

  
  lambda_opt=lambdas[id_lambda].get();
  lambda_opt_2 = lambdas[id_lambda2].get();
  ch2_ret=ch2_thresh;

  //print reconstruction corresponding to lambda and lambda2

  //lambda
  Print_R_at_lambda<<lambda_opt<<" "<<A_tofit[id_lambda]<<" "<<B_tofit[id_lambda]<<" "<<R_tofit[id_lambda].ave()<<" "<<R_tofit[id_lambda].err();
  Print_R_at_lambda<<" "<<2<<" "<<A_tofit[id_lambda]/B_tofit[id_lambda]<<endl;
  //lambda2
  Print_R_at_lambda<<lambda_opt_2<<" "<<A_tofit[id_lambda2]<<" "<<B_tofit[id_lambda2]<<" "<<R_tofit[id_lambda2].ave()<<" "<<R_tofit[id_lambda2].err();
  Print_R_at_lambda<<" "<<3<<" "<<A_tofit[id_lambda2]/B_tofit[id_lambda2]<<endl;


  Print_R_at_lambda.close();

  cout<<"done!"<<endl<<flush;
  
  return;

}


void Get_optimal_lambda(const PrecMatr &Atr, const PrecMatr &Btr,const PrecVect &ft, const PrecFloat & M2, double& lambda_opt, double& lambda_opt_2,  const distr_t_list & corr, int tmin, int tmax,const double mult, double mult2, double Ag_ov_A0_tg,  string path ) {

  bool UseJack= corr.UseJack;
 

  int MAX_Iters = 1000;

  int Global_id=0;

  ofstream Print_R_at_lambda(path);

  //Print header
  Print_R_at_lambda<<"# $1=lambda,  $2= A[g]/A[0], $3=B[g], $4= val, $5=err, $6=S_FLAG $7=mult"<<endl; 
  
  Print_R_at_lambda.precision(10);
    
  if(verbosity_lev) cout<<"Finding optimal lambda* ..."<<endl;

 

 
   
  const auto A=
    [&M2, &Atr, &ft](const PrecVect& gmin) -> PrecFloat
    {
      PrecVect Atr_g = Atr*gmin;
      PrecFloat g_Atr_g = gmin.transpose()*Atr_g;
      PrecFloat ft_g = ft.transpose()*gmin;
      return (M2 + g_Atr_g -2*ft_g)/M2  ;
    };

  const auto B=
    [&Btr](const PrecVect& gmin) -> PrecFloat
    {
      PrecVect B_g = Btr*gmin;
      PrecFloat g_B_g = gmin.transpose()*B_g;
      return g_B_g ;
    };

 
  //bisection search
  PrecFloat l_low =0;
  PrecFloat l_up = 1;
  PrecFloat l_start=1.0;
  int Nit=0;
  int Nit_Ag0=0;
  int Nit_2=0;
  int Nit_100=0;
  double lambda_balance;
  double lambda_balance_2;
  double lambda_balance_100;
  bool lambda_balance_found=false;
  bool lambda_balance_found_2=false;
  PrecFloat Ag_ov_A0_target=1e-3;
  if(Ag_ov_A0_tg > 0) Ag_ov_A0_target= PrecFloat(Ag_ov_A0_tg);
  Vfloat Ags_mult ={10.0, 1.0}; 
 
  if(verbosity_lev) cout<<"Finding lambdas corresponding to A[g]/A[0] in  "<<Ag_ov_A0_target<<" * { 1, 10}"<<endl;

  //###############################################################################################
  

  
  

  //################################################################################################
  
  int counter_Ag_m=0;
 
  for(auto & Ag_m: Ags_mult) {
  
    bool lambda_found_Ag_A0=false;
    Nit_Ag0=0;
    l_low=0;
    
    //bisection search for given A[g]/A[0]
    while( !lambda_found_Ag_A0 ) {

     
      PrecFloat lambda_mid  =  (Nit_Ag0==0 && counter_Ag_m==0)?l_start:(l_up+l_low)/2;
      PrecMatr C = Atr*(1-lambda_mid)/M2 + Btr*lambda_mid;
      PrecMatr C_inv = C.inverse();
      PrecVect ft_l = ft*(1-lambda_mid)/M2;
     
      PrecVect gm = C_inv*ft_l;
         
      PrecFloat A_val = A(gm);
      PrecFloat B_val = B(gm);
      PrecFloat W_val = (1-lambda_mid)*A_val + lambda_mid*B_val;

     
      cout.precision(10);

      //##########################################################################################
      //compute anti-Laplace transform corresponding to lambda_mid:
      distr_t R_E_lambda(UseJack);
      
      for(int ijack=0; ijack< corr.distr_list[1].size(); ijack++) {

	PrecFloat spec_lambda_d_jack=0;

	for(int t=tmin;t<=tmax;t++) spec_lambda_d_jack = spec_lambda_d_jack + gm(t-tmin)*corr.distr_list[t].distr[ijack]; 

	R_E_lambda.distr.push_back( spec_lambda_d_jack.get());
      }

      PrecFloat mult_est = A_val/B_val;
     
      Print_R_at_lambda<<lambda_mid<<" "<<A_val<<" "<<B_val<<" "<<R_E_lambda.ave()<<" "<<R_E_lambda.err();
      Print_R_at_lambda<<" 0 "<<mult_est<<endl;

      //##########################################################################################

     
   

      if(A_val > Ag_ov_A0_target*Ag_m) { // lambda_mid is new l_low
	l_up =lambda_mid;
      }
      else { //lambda_mid is new l_up
	l_low = lambda_mid;
      }

      
      
      Nit_Ag0++;
      if( (A_val < 1.5*Ag_ov_A0_target*Ag_m && A_val > 0.7*Ag_ov_A0_target*Ag_m)) lambda_found_Ag_A0=true;

      if(Nit_Ag0 >= 60) { 
	lambda_found_Ag_A0=true;
      }
      Global_id++;
    }
   
    counter_Ag_m++;
  }

  if(verbosity_lev) cout<<"A[g]/A[0] scan completed!"<<endl;
  cout.precision(10);
 


  //#############################################################################################################################################
  //#############################################################################################################################################
  //#############################################################################################################################################
  //#############################################################################################################################################
  
  if(verbosity_lev) cout<<"Finding lambda from balance condition A=mult*B, mult= "<<mult<<endl;
  l_up=1.0;
  l_low =0.0;
  
 
  //bisection search for condition A = mult*B
  while( !lambda_balance_found ) {

   
    //evaluate the minimum at midpoint
    PrecFloat lambda_mid = (Nit==0)?l_start:(l_up+l_low)/2;
    PrecMatr C = Atr*(1-lambda_mid)/M2 + Btr*lambda_mid;
    PrecMatr C_inv = C.inverse();
    PrecVect ft_l = ft*(1-lambda_mid)/M2;
      
    PrecVect gm = C_inv*ft_l;

    PrecFloat A_val = A(gm);
    PrecFloat B_val = B(gm);
    PrecFloat W_val = (1-lambda_mid)*A_val + lambda_mid*B_val;
    
          
    double mult_est = (A_val/B_val).get();
         
    if(mult> mult_est) { // lambda_mid is new l_low
      l_low =lambda_mid;
    }
    else { //lambda_mid is new l_up
      l_up = lambda_mid;
    }
    
    lambda_balance= lambda_mid.get();
    
    Nit++;
    if(fabs(mult_est - mult)/(mult) < 0.01) lambda_balance_found=true;

    
    if(Nit > MAX_Iters) {
      cout<<"###### FAILED CONVERGENCE #########"<<endl;
      cout<<"###### INFO #######################"<<endl;
      cout<<"-----------------------------"<<endl;
      cout<<"Current values of g[t] & f[t]:"<<endl;
      for(int t=tmin;t<=tmax;t++) cout<<"t: "<<t<<"   "<<gm(t-tmin).get()<<"      "<<ft(t-tmin)<<endl;
      cout<<"-----------------------------"<<endl;
      cout<<"M2: "<<M2.get()<<endl;
      cout<<"A[g]/A[0]: "<<A_val<<", B[g]: "<<B_val<<endl;
      cout<<"Printing inverse W_tr = A_tr*(1-lambda) + B_tr*lambda matrix: "<<endl;
      cout<<"-------------------------------------------------------------"<<endl;
      cout<<C_inv<<endl;
      cout<<"-------------------------------------------------------------"<<endl;
      cout<<"Printing inverse A_tr matrix: "<<endl;
      cout<<"-------------------------------------------------------------"<<endl;
      cout<<Atr.inverse()<<endl;
      cout<<"-------------------------------------------------------------"<<endl;
      cout<<"lambda_low: "<<l_low<<" lambda_up: "<<l_up<<endl;
      cout<<"####################################"<<endl;
      crash("After "+to_string(Nit)+" iterations, balance condition A = mult*B cannot be obtained");
    
      
    }


    //##########################################################################################
    //compute anti-Laplace transform corresponding to lambda_mid:
    distr_t R_E_lambda(UseJack);
      
    for(int ijack=0; ijack< corr.distr_list[1].size(); ijack++) {

      PrecFloat spec_lambda_d_jack=0;

      for(int t=tmin;t<=tmax;t++) spec_lambda_d_jack = spec_lambda_d_jack + gm(t-tmin)*corr.distr_list[t].distr[ijack]; 

      R_E_lambda.distr.push_back( spec_lambda_d_jack.get());
    }

      

    Print_R_at_lambda<<lambda_mid<<" "<<A_val<<" "<<B_val<<" "<<R_E_lambda.ave()<<" "<<R_E_lambda.err();
    Print_R_at_lambda<<" "<<lambda_balance_found<<" "<<mult_est<<endl;

    //##########################################################################################

     
    
    Global_id++;
  }


  if(verbosity_lev) cout<<"lambda_opt = "<<lambda_balance<<endl;

  //#############################################################################################################################################
  //#############################################################################################################################################
  //#############################################################################################################################################
  //#############################################################################################################################################







  //#############################################################################################################################################
  //#############################################################################################################################################
  //#############################################################################################################################################
  //#############################################################################################################################################

  double k=mult2/mult;
  l_low =0.0;
  if(verbosity_lev) cout<<"Finding lambda from balance condition A= mult2*B, mult2= "<<mult2<<endl;
 
 
  //bisection search for condition A = mult*B
  while( !lambda_balance_found_2 ) {

   

    //evaluate the minimum at midpoint
    PrecFloat lambda_mid = (Nit_2==0)?l_up:(l_up+l_low)/2;
    PrecMatr C = Atr*(1-lambda_mid)/M2 + Btr*lambda_mid;
    PrecMatr C_inv = C.inverse();
    PrecVect ft_l = ft*(1-lambda_mid)/M2;
      
    PrecVect gm = C_inv*ft_l;

  
    PrecFloat A_val = A(gm);
    PrecFloat B_val = B(gm);
    PrecFloat W_val = (1-lambda_mid)*A_val + lambda_mid*B_val;
  
      
    double mult_est = (A_val/B_val).get() ;
         
    if(k*mult> mult_est) { // lambda_mid is new l_low
      l_low =lambda_mid;
    }
    else { //lambda_mid is new l_up
      l_up = lambda_mid;
    }

    lambda_balance_2= lambda_mid.get();
      
    Nit_2++;
    if(fabs(mult_est - k*mult)/(k*mult) < 0.01) lambda_balance_found_2=true;

      

    if(Nit_2 > MAX_Iters) {
      cout<<"###### FAILED CONVERGENCE #########"<<endl;
      cout<<"###### INFO #######################"<<endl;
      cout<<"-----------------------------"<<endl;
      cout<<"Current values of g[t] & f[t]:"<<endl;
      for(int t=tmin;t<=tmax;t++) cout<<"t: "<<t<<"   "<<gm(t-tmin).get()<<"      "<<ft(t-tmin)<<endl;
      cout<<"-----------------------------"<<endl;
      cout<<"M2: "<<M2.get()<<endl;
      cout<<"A[g]/A[0]: "<<A_val<<", B[g]: "<<B_val<<endl;
      cout<<"Printing inverse W_tr = A_tr*(1-lambda) + B_tr*lambda matrix: "<<endl;
      cout<<"-------------------------------------------------------------"<<endl;
      cout<<C_inv<<endl;
      cout<<"-------------------------------------------------------------"<<endl;
      cout<<"Printing inverse A_tr matrix: "<<endl;
      cout<<"-------------------------------------------------------------"<<endl;
      cout<<Atr.inverse()<<endl;
      cout<<"-------------------------------------------------------------"<<endl;
      cout<<"lambda_low: "<<l_low<<" lambda_up: "<<l_up<<endl;
      cout<<"####################################"<<endl;
      crash("After "+to_string(Nit)+" iterations, balance condition A = mult*B cannot be obtained");
    }


    //##########################################################################################
    //compute anti-Laplace transform corresponding to lambda_mid:
    distr_t R_E_lambda(UseJack);
      
    for(int ijack=0; ijack< corr.distr_list[1].size(); ijack++) {

      PrecFloat spec_lambda_d_jack=0;

      for(int t=tmin;t<=tmax;t++) spec_lambda_d_jack = spec_lambda_d_jack + gm(t-tmin)*corr.distr_list[t].distr[ijack]; 

      R_E_lambda.distr.push_back( spec_lambda_d_jack.get());
    }

    Print_R_at_lambda<<lambda_mid<<" "<<A_val<<" "<<B_val<<" "<<R_E_lambda.ave()<<" "<<R_E_lambda.err();
    Print_R_at_lambda<<" "<<2*lambda_balance_found_2<<" "<<mult_est<<endl;

    //##########################################################################################

    Global_id++;
    
   
  }


  if(verbosity_lev) cout<<"lambda_opt_2 = "<<lambda_balance_2<<endl;

  //#############################################################################################################################################
  //#############################################################################################################################################
  //#############################################################################################################################################
  //#############################################################################################################################################



  k=0.01;
  l_low =0.0;
  if(verbosity_lev) cout<<"Finding lambda from balance condition A="<<to_string_with_precision(k,2)<<"*mult*B, mult: "<<mult<<endl;
  bool lambda_balance_found_100=false;
 
 
  //bisection search for condition A = mult*B
  while( !lambda_balance_found_100 ) {

   

    //evaluate the minimum at midpoint
    PrecFloat lambda_mid = (Nit_100==0)?l_up:(l_up+l_low)/2;
    PrecMatr C = Atr*(1-lambda_mid)/M2 + Btr*lambda_mid;
    PrecMatr C_inv = C.inverse();
    PrecVect ft_l = ft*(1-lambda_mid)/M2;
       
    PrecVect gm = C_inv*ft_l;

    PrecFloat A_val = A(gm);
    PrecFloat B_val = B(gm);
    PrecFloat W_val = (1-lambda_mid)*A_val + lambda_mid*B_val;
      
    double mult_est = (A_val/B_val).get();
         
    if(k*mult> mult_est) { // lambda_mid is new l_low
      l_low =lambda_mid;
    }
    else { //lambda_mid is new l_up
      l_up = lambda_mid;
    }

    lambda_balance_100= lambda_mid.get();
      
    Nit_100++;
    if(fabs(mult_est - k*mult)/(k*mult) < 0.01) lambda_balance_found_100=true;


    bool skipping_N_100=false;

    if(Nit_100 > 20) {
     
      lambda_balance_found_100=true;
      skipping_N_100=true;
    }


    //##########################################################################################
    //compute anti-Laplace transform corresponding to lambda_mid:
    distr_t R_E_lambda(UseJack);
      
    for(int ijack=0; ijack< corr.distr_list[1].size(); ijack++) {

      PrecFloat spec_lambda_d_jack=0;

      for(int t=tmin;t<=tmax;t++) spec_lambda_d_jack = spec_lambda_d_jack + gm(t-tmin)*corr.distr_list[t].distr[ijack]; 

      R_E_lambda.distr.push_back( spec_lambda_d_jack.get());
    }

    Print_R_at_lambda<<lambda_mid<<" "<<A_val<<" "<<B_val<<" "<<R_E_lambda.ave()<<" "<<R_E_lambda.err();
    Print_R_at_lambda<<" "<<3*(lambda_balance_found_100==1 && skipping_N_100==0)<<" "<<mult_est<<endl;

    //##########################################################################################

     
    Global_id++;
   
  }

  if(verbosity_lev) cout<<"lambda_opt_100 = "<<lambda_balance_100<<endl;


 
  //#############################################################################################################################################
  //#############################################################################################################################################
  //#############################################################################################################################################
  //#############################################################################################################################################

  

  lambda_opt_2 = lambda_balance_2;
  lambda_opt =  lambda_balance;
 
  Print_R_at_lambda.close();

  
  return;

}


CHLT Get_INVLT(int tmin, int tmax, const PrecFloat M2, const function<PrecFloat(int,int)> &Atr_lambda, const function<PrecFloat(int)> &func_f,  Vfloat &covariance, const distr_t_list& corr, double mult, double mult2, double Ag_ov_A0_tg,  string out_path, bool INCLUDE_ERRORS, int prec, string AUTO) {


  CHLT CC;

  double ch2=0.0;


  cout<<"precision used: "<<prec<<endl<<flush;

  
  bool UseJack= corr.distr_list[0].UseJack;
  int Njacks = corr.distr_list[0].size();

  
   
  PrecFloat::setDefaultPrecision(prec);

 
  PrecMatr Atr;
  PrecVect ft;
  PrecMatr B;


 
  
  if(verbosity_lev) cout<<"computing f(t)..."<<flush;
  Get_ft(ft, tmin, tmax,  func_f);
  if(verbosity_lev) { cout<<"done!"<<endl<<flush;}
  Get_Atr(Atr, tmin, tmax, Atr_lambda );

  PrecMatr Atr_2=Atr;
  PrecVect ft_2= ft;

 
  cout.precision(10);
 
  if(INCLUDE_ERRORS) Compute_covariance_matrix(B,tmin,tmax,covariance, corr);

  double lambda_opt= INCLUDE_ERRORS?0.9:0.0;
  double lambda_opt_2=INCLUDE_ERRORS?0.9:0.0;

  if(INCLUDE_ERRORS) {
    if(AUTO == "AUTO") automated_plateaux_search(Atr, B,ft, M2, lambda_opt,lambda_opt_2,ch2, corr, tmin, tmax,out_path+".auto_search");
    else Get_optimal_lambda(Atr, B,ft, M2, lambda_opt,lambda_opt_2, corr, tmin, tmax,mult, mult2, Ag_ov_A0_tg,out_path+".stab_analysis");
  }

  if(verbosity_lev>=2) cout<<"Stability analysis completed!"<<endl;
    							         
  if(INCLUDE_ERRORS) {
    Atr_2 = Atr*(1-lambda_opt_2)/M2 + B*lambda_opt_2;
    Atr = Atr*(1-lambda_opt)/M2 + B*lambda_opt;
    ft_2 = ft*(1-lambda_opt_2)/M2;
    ft= ft*(1-lambda_opt)/M2;
  }


  //invert Atr

  const PrecMatr Atr_inv = Atr.inverse();
  const PrecMatr Atr_inv_2= Atr_2.inverse();

  
  //get g(t) 
  

  PrecVect g = Atr_inv*ft;
  PrecVect g_2= Atr_inv_2*ft_2;

  CC.g = g;
  CC.g_2 = g_2;

  //compute rho

  distr_t rho = 0.0*Get_id_distr(Njacks,UseJack);
  distr_t rho_2 = 0.0*Get_id_distr(Njacks,UseJack);
  
  for(int t=tmin;t<=tmax;t++)  {
    rho = rho + g(t-tmin).get()*corr.distr_list[t];
    rho_2 = rho_2 + g_2(t-tmin).get()*corr.distr_list[t];
  }
  
  CC.rho = rho;

  if(INCLUDE_ERRORS)  {

    if(ch2 <= 1.0) CC.syst = erf(fabs( (rho - rho_2).ave()/(sqrt(2)*rho_2.err())))*fabs( (rho - rho_2).ave());
    else {
      CC.syst = erf(fabs( (rho - rho_2).ave()/(sqrt(2*ch2)*rho_2.err())))*fabs( (rho - rho_2).ave());
      double x= pow(CC.rho.err(),2)*(ch2-1.0); 
      CC.syst = sqrt( pow(CC.syst,2) + x);
      
    }
  }

  return CC;

  
}



