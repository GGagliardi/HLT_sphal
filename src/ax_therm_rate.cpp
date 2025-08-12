#include "../include/ax_therm_rate.h"

using namespace std;

const int prec = 2*256;  //double precision in this scale is 52


void get_axion_therm_rate(double m_T, double s_T, int NT, int Nboots,const distr_t_list &C, string out_path, int INCLUDE_ERRORS) {

  //set precision for precfloat
  PrecFloat::setDefaultPrecision(prec);

  
  

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
  if(INCLUDE_ERRORS == -1) INCLUDE_ERRORS=0;

  
  

  //##########################################//
  //       SET  HLT PARAMETERS         //


  double E0=0.0;
  double Ag_target= 1e-3;  //in the stability analysis it will perform a scan in lambda searching for values of A[g]/A[0] around Ag_target
  double mult= 1e1; // optimal lambda found from condition A[g]/A[0] = mult*B[g]
  int tmin=1;
  int tmax= NT/2;
  double T= 1.0/NT;
  double s= s_T*T;
  double m= m_T*T;

  //##########################################//

  cout<<"Printing info: "<<endl;
  cout<<"tmin: "<<tmin<<endl;
  cout<<"tmax: "<<tmax<<endl;
  cout<<"s/T: "<<s_T<<endl;
  cout<<"m/T: "<<m_T<<endl;
  cout<<"NT: "<<NT<<endl;
  cout<<"Nboots: "<<Nboots<<endl;
  cout<<"a*s: "<<s<<endl;
  cout<<"a*m: "<<m<<endl;
  cout<<"INCLUDE ERRORS: "<<INCLUDE_ERRORS<<endl;
 
 


  
  //######################################################################
  //######################################################################
  //############## LAMBDA FUNCTIONS TO BE USED IN HLT-RECO  ##############
  //######################################################################
  //######################################################################

  // B(n,E) is the basis function to be used. The standard one used at zero-temperature is B(n,E) = exp(-nE);
  // at m/T=0 (sphaleron rate) the basis used is B(n,E) = (1/2*PI)*E*NT*cosh( E*( n - NT/2))/sinh(E*NT/2)
  // at m/T != 0 (axion therm rate) the basis used is B(n,E) = (1/2*PI)*(1 - exp(-E*NT))*cosh( E*( n - NT/2))/sinh(E*NT/2)
  
  auto Bn = [&NT](int n, PrecFloat E) -> PrecFloat {
    //if(E==0) return 1.0/precPi();

    return  (1.0/(2.0*precPi()))*exp(-E*PrecFloat(1.0*n))*(1.0+ exp(-E*PrecFloat((1.0*NT-2.0*n))))/safe_1m_expmx_ov_x(E*NT);
    
    //return (1.0/(2.0*precPi()))*exp(-E*n)*(1+ exp(-E*(NT-2*n)))*E*NT/(1-exp(-E*NT));
    
    //return (1.0/(2.0*precPi()))*( exp(-E*PrecFloat(1.0*n)) + exp(-E*PrecFloat(1.0*(NT-n)))); //exp(-E*n)*(1+ exp(-E*(NT-2*n))); 
  };

  //set the function G we want to approximate
  auto G = [&NT, &m, &s_T, &s](PrecFloat E) -> PrecFloat {
    PrecFloat x = (E-m)/s;
    //if(x==0) return  PrecFloat((1.0/s_T))*pow(2/precPi(),2);   
    return (1.0/s_T)*pow(2/precPi(),2)/safe_sinhx_ov_x(x);

    //return (1.0/s_T)*exp( -x*x/2); 
  };

  // | G*G |
  auto G2 = [&G](PrecFloat E) -> PrecFloat {
    return G(E)*G(E); //*exp(PrecFloat(1.99)*E);
  };

  // Atr(n1,n2) = < B(n1) | B(n2) >
  auto Atr = [&NT, &Bn, &E0](int n1, int n2) -> PrecFloat {
    auto func = [&NT, &Bn, &n1, &n2](PrecFloat E) -> PrecFloat {
      return Bn(n1,E)*Bn(n2,E); //*exp(PrecFloat(1.99)*E);
    };
    return  integrateUpToInfinite(func, E0, false);   //called with true sets verbosity level to 1
  };

  // ft(n) = < B(n) | f > ,  where f is the function we want to approximate
  auto ft = [&NT, &Bn, &G, &E0](int n) -> PrecFloat {
    auto func = [&Bn, &G, &n](PrecFloat E) -> PrecFloat {
      PrecFloat ret= Bn(n,E)*G(E); //*exp(PrecFloat(1.99)*E);
      return ret;
    };
    return  integrateUpToInfinite(func, E0, false);   //called with true sets verbosity level to 1
  };

  //M2 = || f ||^2  
  PrecFloat M2 =  integrateUpToInfinite(G2, E0, false);   //called with true sets verbosity level to 1

  //######################################################################
  //######################################################################
  //######################################################################
  //######################################################################

  
  
  //Do HLT analysis
  CHLT ax_therm_rate =  Get_INVLT(tmin,tmax, M2, Atr, ft,  Cov, C, mult, Ag_target, out_path, INCLUDE_ERRORS, prec);  //in this example tmin = 1 , tmax = NT/2
  


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
    PrecFloat reco_result  = 0;// (ax_therm_rate.g).transpose()*bt;
    for(int t=tmin;t<=tmax;t++) reco_result += (ax_therm_rate.g(t-tmin))*Bn(t,E);
    Reco.push_back( reco_result.get());
    Exact.push_back( G(E).get());
    Erg.push_back(E.get()*NT);
  }

  Print_To_File({} , {Erg, Reco, Exact} , out_path+".reco", "", "");
  


  return;
}
