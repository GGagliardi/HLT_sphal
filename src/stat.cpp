#include "../include/stat.h"

using namespace std;


Pfloat JackAve(const Vfloat &JackDistr) {

  int Nclusters = JackDistr.size();
  Pfloat res= make_pair(0.0,0.0);
  for(auto & data: JackDistr) res.first += data/Nclusters;
  for(auto & data: JackDistr) res.second += pow((data - res.first),2);

  res.second = sqrt( ((Nclusters-1.0)/(double)Nclusters)*res.second);

 
  return res;
}

double Compute_jack_cov(const Vfloat& A,const Vfloat& B) {

  if(A.size() != B.size()) crash(" In Compute_jack_cov: Number of clusters for the two observables A,B don't match");

  int Nclusters = A.size();

  double barA = JackAve(A).first;
  double barB = JackAve(B).first;

  double result=0;

  for(int ijack=0;ijack<Nclusters;ijack++) result += (A[ijack]-barA)*(B[ijack]-barB);

  //Vfloat C = Multiply_vectors(A,B);
  //double AB = JackAve(C).first;
  return (Nclusters-1.0)*result/Nclusters;
  
}

void Compute_covariance_matrix(bool UseJack, Eigen::MatrixXd& Cov, int nargs,...) {

  
  Cov.resize(nargs,nargs);
  va_list args;
  va_start(args, nargs);
  VVfloat A(nargs);
  for(int i=0; i<nargs;i++) A[i] = va_arg(args,Vfloat);
  va_end(args);

  for(int i=0; i<nargs;i++)
    for(int j=i;j<nargs;j++) {
      double res;
      if(UseJack) res= Compute_jack_cov(A[i], A[j]);
      else res = Compute_boot_cov(A[i], A[j]);
      Cov(i,j) = res;
      Cov(j,i) = res;
    }

  return;
}


void Compute_correlation_matrix(bool UseJack, Eigen::MatrixXd& Corr, int nargs,...) {

  
  Corr.resize(nargs,nargs);
  va_list args;
  va_start(args, nargs);
  VVfloat A(nargs);
  for(int i=0; i<nargs;i++) A[i] = va_arg(args,Vfloat);
  va_end(args);

  for(int i=0; i<nargs;i++)
    for(int j=i;j<nargs;j++) {
      double num,den;
      if(UseJack) {
	num= Compute_jack_cov(A[i], A[j]);
	den = sqrt( Compute_jack_cov(A[i], A[i]) * Compute_jack_cov(A[j], A[j]));
      }
      else {
	num = Compute_boot_cov(A[i], A[j]);
	den = sqrt( Compute_boot_cov(A[i], A[i]) * Compute_boot_cov(A[j], A[j]));
      }
      Corr(i,j) = num/den;
      Corr(j,i) = num/den;
    }

  return;
}


void Compute_autocorrelation_time(const Vfloat &data, string Path, string Tag) {

  //compute empirical autocorrelation function


  auto GT = [](const Vfloat & data,Vfloat &rho, Vfloat  &tau_int, int tmax)  {

    int N= data.size();

    rho.clear();
    tau_int.clear();

    double t_accumulated = -1.0;
    
    for(int t=0; t< tmax;t++) {

      double bar_x_t=0;
      double bar_y_t=0;
      double num=0;
      double den=0;
      double den1=0;
      double den2=0;
      for(int i=0;i< N-t;i++) bar_x_t += (1.0/((double)(N-t)))*data[i];
      for(int i=0;i< N-t;i++) bar_y_t += (1.0/((double)(N-t)))*data[t+i];
      for(int i=0;i< N-t;i++) num += (data[i] - bar_x_t)*(data[t+i]-bar_y_t);
      for(int i=0;i< N-t;i++) {den1 += pow( data[i] - bar_x_t, 2); den2 += pow( data[t+i]-bar_y_t,2);}
      den = sqrt( den1*den2);
      rho.push_back(num/den);
      t_accumulated += 2.0*num/den;
      tau_int.push_back(t_accumulated);
      
    }

    return;
    

  };

  int tmax= data.size()/8 - 20;

  Vfloat rho, tau_int;
  GT(data, rho, tau_int, tmax);

  int N=data.size();
  int N4= N/4;
  int N8= N/8;

  //divide sample in 4 parts and compute autocorr on each part
  VVfloat data4;
  VVfloat rho4(4), tau_int4(4);
  for(int i=0;i<4;i++) {
    data4.emplace_back( data.begin() + i*N4, data.begin() + (i+1)*N4);
    assert((signed)data4[i].size() == N4);
    GT(data4[i], rho4[i], tau_int4[i], tmax);
  }



  //divide sample in 8 parts and compute autocorr on each part
  VVfloat data8;
  VVfloat rho8(8), tau_int8(8);
  for(int i=0;i<8;i++) {
    data8.emplace_back( data.begin() + i*N8, data.begin() + (i+1)*N8);
    assert((signed)data8[i].size() == N8);
    GT(data8[i], rho8[i], tau_int8[i], tmax);
  }

 
  //get errors

  Vfloat rho_err_4, tau_err_4;
  Vfloat rho_err_8, tau_err_8;

  for(int t=0;t<tmax;t++) {

    double err_rho=0; double err_tau=0;

    for(int i=0;i<4;i++) err_rho += pow(rho4[i][t] - rho[t],2)/3.0;
    err_rho = sqrt(err_rho);
    for(int i=0;i<4;i++) err_tau += pow(tau_int4[i][t] - tau_int[t],2)/3.0;
    err_tau = sqrt(err_tau);
    
    rho_err_4.push_back( err_rho);
    tau_err_4.push_back( err_tau);


    err_rho=0; err_tau=0;

    for(int i=0;i<8;i++) err_rho += pow(rho8[i][t] - rho[t],2)/7.0;
    err_rho = sqrt(err_rho);
    for(int i=0;i<8;i++) err_tau += pow(tau_int8[i][t] - tau_int[t],2)/7.0;
    err_tau = sqrt(err_tau);
    
    rho_err_8.push_back( err_rho);
    tau_err_8.push_back( err_tau);

  }

  
 
  
  

  
  
  Print_To_File({}, {rho,rho_err_4, rho_err_8, tau_int, tau_err_4, tau_err_8}, Path+"/autocorr_"+Tag+".dat", "", "");



  return;
}








Pfloat BootAve(const Vfloat &BootDistr) {

  Pfloat ave_err = make_pair(0.0,0.0);


  int N= BootDistr.size();
  for(int i=0; i<N; i++)  ave_err.first += BootDistr[i]/N;
  for(int i=0; i<N; i++) ave_err.second += pow( BootDistr[i]-ave_err.first,2)/(1.0*N);

  ave_err.second = sqrt(ave_err.second);
  return ave_err;
}



double Compute_boot_cov(const Vfloat& A,const Vfloat& B) {

    double barA = BootAve(A).first;
    double barB = BootAve(B).first;

  
    int N=A.size();
    double res=0.0;
    if(N != (signed)B.size()) crash("In Compute_boot_cov: A.size() != B.size()");
    for(int i=0;i<N;i++) res += (1.0/(N))*(A[i]-barA)*(B[i] - barB); 
    return res;
}




double distr_t::ave() const {
  
   if(UseJack) return JackAve(this->distr).first;

   return BootAve(distr).first;
 }
double distr_t::err() const {
  if(UseJack) {
    return JackAve(this->distr).second;
  }
  else return BootAve(this->distr).second;
}

Pfloat distr_t::ave_err() const {
   if(UseJack) return JackAve(this->distr);
   else return BootAve(distr);
 }
 
int distr_t::size() const {return distr.size();}

Vfloat distr_t_list::ave() const {
    Vfloat res;
    for(int i=0; i < this->size();i++) res.push_back(distr_list[i].ave());
    return res;
  }

double distr_t_list::ave(int i_distr) const {
    if(i_distr >= (signed)distr_list.size()) crash("In distr_t_list function ave called with positional argument greater than distr_list size");
    return distr_list[i_distr].ave();
  }
Vfloat distr_t_list::err()  const{
    Vfloat res;
    for(int i=0; i < this->size();i++) res.push_back(this->distr_list[i].err());
    return res;
  }
double distr_t_list::err(int i_distr) const {
    if(i_distr >=  (signed)distr_list.size()) crash("In distr_t_list function ave called with positional argument greater than distr_list size");
    return this->distr_list[i_distr].err();
  }

vector<Pfloat> distr_t_list::ave_err() const {
    vector<Pfloat> res;
    for(int i=0; i < this->size();i++) res.push_back(this->distr_list[i].ave_err());
    return res;
  }


  
Pfloat distr_t_list::ave_err(int i_distr) const{
     if(i_distr >=  (signed)distr_list.size()) crash("In distr_t_list function ave called with positional argument greater than distr_list size");
     return distr_list[i_distr].ave_err();
  }

int distr_t_list::size() const { return distr_list.size();}

Vfloat distr_t_list::Get_distr_index(int k) const {
    Vfloat res;
    for(int i=0; i < this->size(); i++) {
      if(this->distr_list[i].size() <= k ) crash("In distr_t_list, call to I_list(int k) invalid. Positional argument k is too large");
      res.push_back(this->distr_list[i].distr[k]);
    }
    return res;
}

//operator overloading

distr_t operator*(const distr_t& A, const distr_t& B) {

  if(A.size() != B.size()) crash("In distr_t, call to A*B is invalid. A and B have different sizes");
   distr_t res(A.UseJack,A.size());
         
   for(int i=0; i < A.size(); i++) res.distr[i] = A.distr[i]*B.distr[i];
	 
   return res;
 }
 
 distr_t operator+(const distr_t& A,const  distr_t& B) {
   if(A.size() != B.size()) crash("In distr_t, call to A+B is invalid. A and B have different sizes");
   distr_t res(A.UseJack, A.size());
   for(int i=0; i < A.size(); i++) res.distr[i] = A.distr[i]+B.distr[i];
   return res;
 }
 
 distr_t operator-(const distr_t& A, const distr_t& B) {
   if(A.size() != B.size()) crash("In distr_t, call to A-B is invalid. A and B have different sizes");
   distr_t res(A.UseJack,A.size());
   for(int i=0; i < A.size(); i++) res.distr[i] = A.distr[i]-B.distr[i];
   return res;
 }

 distr_t operator/(const distr_t& A, const distr_t& B) {
   if(A.size() != B.size()) crash("In distr_t, call to A/B is invalid. A and B have different sizes");
   distr_t res(A.UseJack,A.size());
   for(int i=0; i < A.size(); i++) {
     //if(fabs(B.distr[i]) < eps(16)) crash("In dist_t, call to A/B is invalid. B has a zero element: "+to_string_with_precision(B.distr[i], 16));
     res.distr[i] = A.distr[i]/B.distr[i];
   }
   return res;
 }

 double operator%(const distr_t& A,const  distr_t& B) {

   if(A.UseJack != B.UseJack) crash("Error in distr_t%distr_t, the two distribution are not using the same method");
   if(A.UseJack) return Compute_jack_cov(A.distr, B.distr);
   else return Compute_boot_cov(A.distr, B.distr);
 }

  

 distr_t operator*(const distr_t& A, double B) {
   distr_t res(A.UseJack,A.size());
   for(int i=0; i < A.size(); i++) res.distr[i] = A.distr[i]*B;
   return res;
 }

 distr_t operator*(double B,  const distr_t& A) {
   return A*B;
 }

 distr_t operator+(const distr_t& A, double B) {
  
   distr_t res(A.UseJack,A.size());
   for(int i=0; i < A.size(); i++) res.distr[i] = A.distr[i]+B;
   return res;
 }

 distr_t operator+(double B,const distr_t& A) {
   return A+B;
 }

 distr_t operator-(const distr_t& A, double B) {
   distr_t res(A.UseJack,A.size());
   
   for( int i=0; i < A.size(); i++) res.distr[i] = A.distr[i]-B;
   return res;
 }

 distr_t operator-(double B,const distr_t& A) {
   distr_t res(A.UseJack,A.size());
   for(int i=0; i < A.size(); i++) res.distr[i] = B-A.distr[i];
   return res;
 }

 distr_t operator/(const distr_t& A, double B) {
   distr_t res(A.UseJack,A.size());
   //if(fabs(B)<eps(16)) crash("cannot call operator distr_t/double with double=0");
   for( int i=0; i < A.size(); i++) res.distr[i] = A.distr[i]/B;
   return res;
 }

 distr_t operator/(double B,const distr_t& A) {
   distr_t res(A.UseJack,A.size());
   for( int i=0; i < A.size(); i++) {
     //if(fabs(A.distr[i])<eps(16)) crash("Cannot call double/distr_t with distr_t containing zeros"); 
     res.distr[i] = B/A.distr[i];
   }
     return res;
 }

 distr_t operator*(const distr_t& A,const Vfloat& B) {
   if(A.size() != (signed)B.size()) crash("In distr_t, call to A*B is invalid. A and B have different sizes");
   distr_t res(A.UseJack,A.size());
   for(int i=0; i < A.size(); i++) res.distr[i] = A.distr[i]*B[i];
   return res;
 }
 
 distr_t operator+(const distr_t& A, const Vfloat& B) {
   if(A.size() != (signed)B.size()) crash("In distr_t, call to A+B is invalid. A and B have different sizes");
   distr_t res(A.UseJack, A.size());
   for( int i=0; i < A.size(); i++) res.distr[i] = A.distr[i]+B[i];
   return res;
 }
 
 distr_t operator-(const distr_t& A, const Vfloat& B) {
   
   if(A.size() != (signed)B.size()) crash("In distr_t, call to A-B is invalid. A and B have different sizes");
   distr_t res(A.UseJack,A.size());
   for( int i=0; i < A.size(); i++) res.distr[i] = A.distr[i]-B[i];
   return res;
 }
 
 distr_t operator/(const distr_t& A,const Vfloat& B) {
   if(A.size() != (signed)B.size()) crash("In distr_t, call to A/B is invalid. A and B have different sizes");
   distr_t res(A.UseJack,A.size());
   for( int i=0; i < A.size(); i++) {
     //if(fabs(B[i]) < eps(16)) crash("In dist_t, call to A/B is invalid. B has a zero element");
     res.distr[i] = A.distr[i]/B[i];
   }
   return res;
 }

 distr_t operator*(const Vfloat &B,const distr_t& A) {
   return A*B;
 }
 distr_t operator+(const Vfloat &B, const distr_t& A) {
   return A+B;
 }
 distr_t operator-(const Vfloat& B, const distr_t& A) {
   return -1.0*A + B;
 }
 distr_t operator/(const Vfloat& B, const distr_t& A) {
   if(A.size() != (signed)B.size()) crash("In distr_t, call to A/B is invalid. A and B have different sizes");
   distr_t res(A.UseJack,A.size());
   for(int i=0; i < A.size(); i++) {
     //if(fabs(A.distr[i]) < eps(16)) crash("In dist_t, call to A/B is invalid. B has a zero element");
     res.distr[i] = B[i]/A.distr[i];
   }
   return res;
 }


  
distr_t_list operator*(const distr_t_list& A, const distr_t_list& B) {
  if(A.size() != B.size()) crash("In distr_t_list, call to A*B is invalid. A and B have different sizes");
  distr_t_list res(A.UseJack,A.size());
  for(int i=0; i < A.size(); i++) res.distr_list[i] =A.distr_list[i]*B.distr_list[i];
    return res;
}

distr_t_list operator+(const distr_t_list& A,const distr_t_list& B) {

  if(A.size() != B.size()) crash("In distr_t_list, call to A+B is invalid. A and B have different sizes");
  distr_t_list res(A.UseJack,A.size());
  for(int i=0; i < A.size(); i++) res.distr_list[i] = A.distr_list[i]+B.distr_list[i];
  return res;
}

distr_t_list operator-(const distr_t_list& A, const distr_t_list& B) {

  if(A.size() != B.size()) crash("In distr_t_list, call to A-B is invalid. A and B have different sizes");
  distr_t_list res(A.UseJack,A.size());
  for(int i=0; i < A.size(); i++) res.distr_list[i] = A.distr_list[i]-B.distr_list[i];
  return res;
}

distr_t_list operator/(const distr_t_list& A,const distr_t_list& B) {
  if(A.size() != B.size()) crash("In distr_t_list, call to A/B is invalid. A and B have different sizes");
  distr_t_list res(A.UseJack,A.size());
  for(int i=0; i < A.size(); i++) res.distr_list[i] = A.distr_list[i]/B.distr_list[i];
    
  return res;
}

Vfloat operator%(const distr_t_list& A,const distr_t_list& B) {
  if(A.size() != B.size()) crash("In distr_t_list, call to A%B is invalid. A and B have different sizes");
  Vfloat cov(A.size());
  for(int i=0;i<A.size();i++) {
    if(A.distr_list[i].size() != B.distr_list[i].size()) crash("In distr_t_list, call to A%B invalid, the two operands have different distr_t sizes");
    cov[i] = A.distr_list[i]%B.distr_list[i];
  }
  return cov;
}

////////////////////////////////

distr_t_list operator*(const distr_t_list& A, const distr_t& B) {

  distr_t_list res(A.UseJack,A.size());
  for(int i=0; i < A.size(); i++) res.distr_list[i] =A.distr_list[i]*B;
  return res;
}
distr_t_list operator*(const distr_t& B, const distr_t_list& A) {return A*B;}

distr_t_list operator+(const distr_t_list& A,const distr_t& B) {
  distr_t_list res(A.UseJack,A.size());
  for( int i=0; i < A.size(); i++) res.distr_list[i] = A.distr_list[i]+B;
  return res;
}
distr_t_list operator+(const distr_t& B, const distr_t_list& A) { return A+B;}

distr_t_list operator-(const distr_t_list& A,const distr_t& B) {
  distr_t_list res(A.UseJack,A.size());
  for( int i=0; i < A.size(); i++) res.distr_list[i] = A.distr_list[i]-B;
  return res;
}

distr_t_list operator-(const distr_t& B,const distr_t_list& A) {
  distr_t_list res(A.UseJack,A.size());
  for(int i=0; i < A.size(); i++) res.distr_list[i] = B- A.distr_list[i];
  return res;
}

distr_t_list operator/(const distr_t_list& A,const distr_t& B) {
  distr_t_list res(A.UseJack,A.size());
  for(int i=0; i < A.size(); i++) res.distr_list[i] = A.distr_list[i]/B;
  return res;
}

distr_t_list operator/(const distr_t& B,const distr_t_list& A) {
  distr_t_list res(A.UseJack,A.size());
  for( int i=0; i < A.size(); i++) res.distr_list[i] = B/A.distr_list[i];
  return res;
}
  

Vfloat operator%(const distr_t_list& A,const distr_t& B) {
  Vfloat cov(A.size());
  for( int i=0;i<A.size();i++) {
    if(A.distr_list[i].size() != B.size()) crash("In distr_t_list, call to A%B<distr_t_list, distr_t>  invalid, the two operands have different distr_t sizes");
    cov[i] = A.distr_list[i]%B;
  }
  return cov;
}

Vfloat operator%(const distr_t& B, const distr_t_list& A) {
  return A%B;
}

///////////////////////////////


distr_t_list operator*(const distr_t_list& A, double B) {
  distr_t_list res(A.UseJack,A.size());
    for( int i=0; i < A.size(); i++) res.distr_list[i] = A.distr_list[i]*B;
    return res;
}

distr_t_list operator*(double B, const distr_t_list& A) {return A*B;}

distr_t_list operator+(const distr_t_list& A, double B) {
  distr_t_list res(A.UseJack,A.size());
  for(int i=0; i < A.size(); i++) res.distr_list[i] = A.distr_list[i]+B;
     return res;
}

distr_t_list operator+(double B, const distr_t_list& A) { return A+B;}

distr_t_list operator-(const distr_t_list& A, double B) {
  distr_t_list res(A.UseJack,A.size());
  for( int i=0; i < A.size(); i++) res.distr_list[i] = A.distr_list[i]-B;
  return res;
  }

distr_t_list operator-(double B, const distr_t_list& A) {
  distr_t_list res(A.UseJack,A.size());
  for(int i=0; i < A.size(); i++) res.distr_list[i] = B- A.distr_list[i];
  return res;
}

distr_t_list operator/(const distr_t_list& A, double B) {
  distr_t_list res(A.UseJack,A.size());
  for( int i=0; i < A.size(); i++) res.distr_list[i] = A.distr_list[i]/B;
  return res;
}
  
distr_t_list operator/(double B, const distr_t_list& A) {
  distr_t_list res(A.UseJack,A.size());
  for(int i=0; i < A.size(); i++) res.distr_list[i] = B/A.distr_list[i];
  return res;
}

distr_t_list operator+(const distr_t_list &A, const Vfloat& B) {
  distr_t_list res(A.UseJack, A.size());
  if(A.size() != (signed)B.size()) crash("Call to operator distr_t_list*Vfloat is invalid, sizeof(Vfloat) and size(distr_t_list) do not coincide");

  for(int i=0; i<A.size();i++) res.distr_list[i] = A.distr_list[i]+B[i];

  return res;
}

distr_t_list operator-(const distr_t_list &A, const Vfloat& B) {
  distr_t_list res(A.UseJack, A.size());
  if(A.size() != (signed)B.size()) crash("Call to operator distr_t_list*Vfloat is invalid, sizeof(Vfloat) and size(distr_t_list) do not coincide");

  for(int i=0; i<A.size();i++) res.distr_list[i] = A.distr_list[i]-B[i];

  return res;
}

distr_t_list operator+(const Vfloat &B, const distr_t_list &A) { return A + B; }
distr_t_list operator-(const Vfloat& B, const distr_t_list& A) { return -1.0*(A-B); }


distr_t_list operator*(const distr_t_list &A, const Vfloat& B) {
  distr_t_list res(A.UseJack, A.size());
  if(A.size() != (signed)B.size()) crash("Call to operator distr_t_list*Vfloat is invalid, sizeof(Vfloat) and size(distr_t_list) do not coincide");

  for(int i=0; i<A.size();i++) res.distr_list[i] = A.distr_list[i]*B[i];

  return res;
}

distr_t_list operator*(const Vfloat& B, const distr_t_list& A) { return A*B;}


distr_t Get_id_jack_distr(int N) {


  distr_t id_jack_distr;

  for(int i=0;i<N;i++) id_jack_distr.distr.push_back(1);

  return id_jack_distr;
}

distr_t Get_id_distr(int N, bool UseJack) {


  distr_t id_jack_distr(UseJack);

  for(int i=0;i<N;i++) id_jack_distr.distr.push_back(1);

  return id_jack_distr;
}

distr_t_list Get_id_jack_distr_list(int size, int N) {
  distr_t_list return_distr(1);
  for(int i=0; i<size;i++) return_distr.distr_list.push_back( Get_id_jack_distr(N));
  return return_distr;
}

distr_t_list Get_id_distr_list(int size, int N, bool UseJack) {
  distr_t_list return_distr(UseJack);
  for(int i=0; i<size;i++) return_distr.distr_list.push_back( Get_id_distr(N,UseJack));
  return return_distr;
}




distr_t_list EXP_DL (const distr_t_list &A) {
  distr_t_list B = A;
  for (int t = 0; t < A.size(); t++)
    for (int i = 0; i < A.distr_list[t].size(); i++)
      B.distr_list[t].distr[i] = exp(A.distr_list[t].distr[i]);
  return B;
}

distr_t_list EXPT_DL(const distr_t_list &A) {
  distr_t_list B = A;
  for (int t = 0; t < A.size(); t++)
    for (int i = 0; i < A.distr_list[t].size(); i++)
      B.distr_list[t].distr[i] = exp(A.distr_list[t].distr[i]*t);
  return B;
}


distr_t_list COSH_DL (const distr_t_list &A) {
  distr_t_list B = A;
  for (int t = 0; t < A.size(); t++)
    for (int i = 0; i < A.distr_list[t].size(); i++)
      B.distr_list[t].distr[i] = cosh(A.distr_list[t].distr[i]);
  return B;
}

distr_t_list SINH_DL (const distr_t_list & A) { distr_t_list B=A;
      for(int t=0;t<A.size();t++)
	for(int i=0;i<A.distr_list[t].size();i++) B.distr_list[t].distr[i] = sinh(A.distr_list[t].distr[i]);
      return B;
 }

distr_t_list SQRT_DL (const distr_t_list &A) {
  distr_t_list B = A;
  for(int t=0;t<A.size();t++)
    for(int i=0;i<A.distr_list[t].size();i++)
      B.distr_list[t].distr[i] = sqrt(A.distr_list[t].distr[i]);
  return B;
}

distr_t_list LOG_DL (const distr_t_list &A) {
  distr_t_list B = A;
  for(int t=0;t<A.size();t++)
    for(int i=0;i<A.distr_list[t].size();i++)
      B.distr_list[t].distr[i] = log(A.distr_list[t].distr[i]);
  return B;
}

distr_t_list POW_DL (const distr_t_list &A, int n) {
  distr_t_list B = A;
  for(int t=0;t<A.size();t++)
    for(int i=0;i<A.distr_list[t].size();i++)
      B.distr_list[t].distr[i] = pow(A.distr_list[t].distr[i],n);
  return B;
}


distr_t_list EXPT_D(const distr_t &A, int T) {
  distr_t_list B(A.UseJack, T, A.size()); 
  for (int t = 0; t < T; t++)
    for (int i = 0; i < A.size(); i++)
      B.distr_list[t].distr[i] = exp(A.distr[i]*t);
  return B;
}


distr_t EXP_D (const distr_t &A) {
  distr_t B = A;
  for (int i = 0; i < A.size(); i++)
    B.distr[i] = exp(A.distr[i]);
  return B;
}

distr_t COSH_D (const distr_t &A) {
  distr_t B = A;
  for (int i = 0; i < A.size(); i++)
    B.distr[i] = cosh(A.distr[i]);
  return B;
}

distr_t SINH_D (const distr_t &A) {
  distr_t B = A;
  for (int i = 0; i < A.size(); i++)
    B.distr[i] = sinh(A.distr[i]);
  return B;
}

distr_t SINHM_D (const distr_t &A) {
  distr_t B = A;
  for (int i = 0; i < A.size(); i++)
    B.distr[i] = asinh(A.distr[i]);
  return B;
}

distr_t ASIN_D (const distr_t &A) {
  distr_t B = A;
  for (int i = 0; i < A.size(); i++)
    B.distr[i] = asin(A.distr[i]);
  return B;
}

distr_t SIN_D (const distr_t &A) {
  distr_t B = A;
  for (int i = 0; i < A.size(); i++)
    B.distr[i] = sin(A.distr[i]);
  return B;
}

distr_t TANH_D (const distr_t &A) {
  distr_t B = A;
  for (int i = 0; i < A.size(); i++)
    B.distr[i] = tanh(A.distr[i]);
  return B;
}

distr_t LOG_D (const distr_t &A) {
  distr_t B = A;
  for (int i = 0; i < A.size(); i++)
    B.distr[i] = log(A.distr[i]);
  return B;
}

distr_t SQRT_D (const distr_t &A) {
  distr_t B = A;
  for (int i = 0; i < A.size(); i++)
    B.distr[i] = sqrt(A.distr[i]);
  return B;
}

distr_t POW_D (const distr_t  &A,int n) {
  distr_t B = A;
  for (int i = 0; i < A.size(); i++)
    B.distr[i] = pow(A.distr[i],n);
  return B;
}



