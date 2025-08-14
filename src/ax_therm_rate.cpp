#include "../include/ax_therm_rate.h"

using namespace std;

const int prec = 2*256;  //double precision in this scale is 52


void get_axion_therm_rate(double m_T, double s_T, int NT, int Nboots,distr_t_list &C, string out_path, int INCLUDE_ERRORS) {

  //set precision for precfloat
  PrecFloat::setDefaultPrecision(prec);
  const PrecFloat PI= precPi();  // define precision PI = 3.14159.....


  //compute covariance matrix
  Vfloat Cov, Corr_matrix;   
  Vfloat TT, RR; 
  for(int tt=0; tt <NT; tt++)
    for(int rr=0;rr <NT;rr++) {
      TT.push_back(tt);
      RR.push_back(rr);
      Cov.push_back(  C.distr_list[tt]%C.distr_list[rr] );
      Corr_matrix.push_back(  C.distr_list[tt]%C.distr_list[rr]/(C.err(tt)*C.err(rr)));
    }
  if(INCLUDE_ERRORS != -1) Print_To_File({}, {TT,RR,Corr_matrix, Cov}, out_path+".Cov", "", "");
  else  INCLUDE_ERRORS=0;

  
  

  //##########################################//
  //       SET  HLT PARAMETERS         //
  double E0=0.0;
  double Emax= 3.0;   // minimizes the norm in the interval [E0, Emax]. If Emax < E0 the interval becomes [E0,\infty]
  PrecFloat alpha=1.99;  // <f|f> = \int_E0^Emax dx e^alpha*x f^2(x) 
  double Ag_target= 1e-3;  //in the stability analysis it will perform a scan in lambda searching for values of A[g]/A[0] around Ag_target
  double mult= 1e-8; // optimal lambda found from condition A[g]/A[0] = mult*B[g]
  double mult2 = 1e-9; //systematic error estimated from difference between results obtained using A[g]/A[0]= mult*B[g] and A[g]/A[0] = mult2*B[g], via erf(). 
  int tmin=1;
  int tmax= NT/2;  
  double T= 1.0/NT; //temperature in lattice units
  double s= s_T*T;  //sigma in lattice units
  double m= m_T*T;  //center of the regularized delta in lattice units
  //##########################################//

  cout<<"Printing info: "<<endl;
  cout<<"[tmin,tmax]:  ["<<tmin<<","<<tmax<<endl;
  cout<<"s/T: "<<s_T<<" as: "<<s<<endl;
  cout<<"m/T: "<<m_T<<" am: "<<m<<endl;
  cout<<"[E0,Emax]: ["<<E0<<","<<Emax<<"] alpha: "<<alpha<<endl;
  cout<<"NT: "<<NT<<endl;
  cout<<"Nboots: "<<Nboots<<endl;
  cout<<"INCLUDE ERRORS: "<<INCLUDE_ERRORS<<endl;
 
 


  
  //######################################################################
  //######################################################################
  //############## LAMBDA FUNCTIONS TO BE USED IN HLT-RECO  ##############
  //######################################################################
  //######################################################################

  // B(n,E) is the basis function to be used. The standard one used at zero-temperature is B(n,E) = exp(-nE);
  // at m/T=0 (sphaleron rate) the basis used is B(n,E) = (1/2*PI)*E*NT*cosh( E*( n - NT/2))/sinh(E*NT/2)  indicated with [1] below
   
  auto Bn = [&](int n, PrecFloat E) -> PrecFloat {
    //return (1.0/(2.0*PI))*exp(-E*n)*(1.0+ exp(-E*(NT-2.0*n)))/safe_1m_expmx_ov_x(E*NT);   [1]
    //return (1.0/(2.0*PI))*exp(-E*n)*(1+ exp(-E*(NT-2*n)))*E*NT/(1-exp(-E*NT));
    return (1.0/(2.0*PI))*(exp(-E*n) + exp(-E*(NT-n)));                              // cosh(x*(n-NT/2))*e^-x*NT/2
  };

  //set the function G we want to approximate
  auto G = [&](PrecFloat E) -> PrecFloat {
    PrecFloat x = (E-m)/s;
    return (1.0/s_T)*pow(2/PI,2)/safe_sinhx_ov_x(x);   //[ x/sinhx  ]
    //return (1.0/(s_T*sqrt(2*PI)))*exp( -x*x/2);      //[ gaussian ]
  };

  // G*G
  auto G2 = [&](PrecFloat E) -> PrecFloat {
    return G(E)*G(E)*exp(alpha*E);
  };

  // Atr(n1,n2) = < B(n1) | B(n2) >
  auto Atr = [&](int n1, int n2) -> PrecFloat {
    auto func = [&](PrecFloat E) -> PrecFloat {
      return Bn(n1,E)*Bn(n2,E)*exp(alpha*E);
    };
    return  (Emax > E0)?integrateUpToXmax(func, E0,Emax, false):integrateUpToInfinite(func,E0,false);   //called with true sets verbosity level to 1
  };

  // ft(n) = < B(n) | G > ,  where G is the function we want to approximate
  auto ft = [&](int n) -> PrecFloat {
    auto func = [&](PrecFloat E) -> PrecFloat {
      PrecFloat ret= Bn(n,E)*G(E)*exp(alpha*E);
      return ret;
    };
    return   (Emax > E0)?integrateUpToXmax(func, E0,Emax, false):integrateUpToInfinite(func,E0,false);   //called with true sets verbosity level to 1
  };

  //M2 = || G ||^2  
  PrecFloat M2 =  (Emax > E0)?integrateUpToXmax(G2, E0,Emax, false):integrateUpToInfinite(G2,E0,false);   //called with true sets verbosity level to 1

  //######################################################################
  //######################################################################
  //######################################################################
  //######################################################################

  
  
  //Do HLT analysis
  CHLT ax_therm_rate =  Get_INVLT(tmin,tmax, M2, Atr, ft,  Cov, C, mult, mult2,  Ag_target, out_path, INCLUDE_ERRORS, prec);  //in this example tmin = 1 , tmax = NT/2
  


  //print bootstrap samples of the smeared axion therm rate

  ofstream print_res(out_path+".res");
  print_res<<"res= "<<ax_therm_rate.rho.ave()<<" +- "<<ax_therm_rate.rho.err()<<" +-  "<<ax_therm_rate.syst<<" (syst)"<<endl;
  for(int iboot=0;iboot<Nboots;iboot++) {
    print_res<<ax_therm_rate.rho.distr[iboot]<<endl;
  }
  print_res.close();


  //print reconstructed and exact kernel function
  vector<double> Reco;  //reconstructed function
  vector<double> Exact; //exact function
  vector<double> Erg;  //sampled energies
   
  double step_size= 0.01;
  long int Npoints= (int)((m+20*s-E0)/(s*step_size)) ;
  cout<<"Npoints: "<<Npoints<<endl;
  for(int ip=0; ip<Npoints;ip++) {
    PrecVect bt;
    PrecFloat E;
    E = E0 + (ip*s)*step_size;
    Get_bt(bt, E, Bn, tmin,tmax);  
    PrecFloat reco_result  = (ax_therm_rate.g).transpose()*bt;
    Reco.push_back(reco_result.get());
    Exact.push_back(G(E).get());
    Erg.push_back(E.get()*NT);
  }

  Print_To_File({} , {Erg, Reco, Exact} , out_path+".reco", "", "");
  


  return;
}
