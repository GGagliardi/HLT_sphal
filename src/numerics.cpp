#include "../include/numerics.h"

using namespace std;


void D(int k) { 
  cout<<"D("<<k<<")"<<endl;
}

double eps(int k) {
  double epsilon= 1;
  for(int i=0;i<k;i++) epsilon/=10;

  return epsilon;
}
long long int fact(int n) {
  if(n==0) return 1;
  else return n*fact(n-1);
}


long long int ipow(int a,int n) {
  long long int result= a;
  if(a < 0 || n < 0) crash("ERR0R: ipow called with negative arguments");
  if(n==0) return 1;
  else return result*ipow(a,n-1);
}

double fpow(double a, int n) {
  double result = a;
  if(n<0) crash("ERROR: fpow called with negative arguments");
  if(n==0) return 1;
  else return result*fpow(a,n-1);
}

long long int BinomialCoeff(int n,int k) {

  vector<long long int> C(k+1,0);
   
    C[0] = 1;  // nC0 is 1 
  
    for (int i = 1; i <= n; i++) 
    { 
        // Compute next row of pascal triangle using 
        // the previous row 
        for (int j = min(i, k); j > 0; j--) 
            C[j] = C[j] + C[j-1]; 
    } 
    return C[k];   
}


double RatioPol(Vfloat &x_num, Vfloat &x_den, Vfloat &Num, Vfloat &Den, Vint &NumPow, Vint &DenPow, string MODE_POL) {

  if(x_num.size() == 0 && x_den.size() == 0) crash("RatioPol called with zero NUM and DEN");

  double NUM_VAL=0;
  double DEN_VAL=0;
  double VAL=0;

  if(MODE_POL=="SINGLE_RATIO") {

  if(x_num.size()==0) NUM_VAL=1.0;
  else {
    for(int i=0;i<(signed)x_num.size();i++) NUM_VAL += pow(x_num[i],NumPow[i])*Num[i];
  }
  if(x_den.size() ==0) DEN_VAL=1.0;
  else {

    for(int j=0; j<(signed)x_den.size();j++) DEN_VAL += pow(x_den[j],DenPow[j])*Den[j];
  }
  VAL= NUM_VAL/DEN_VAL;

  }
  else if(MODE_POL=="MULTIPLE_RATIO") {
    if((signed)x_num.size() != (signed)x_den.size()) crash("RatioPol called in MULTIPLE_RATIO MODE but x_num.size and y_num.size do not coincide");
    for(int i=0; i<(signed)x_num.size();i++) VAL+= pow(x_num[i],NumPow[i])*Num[i]/(pow(x_den[i],DenPow[i])*Den[i]); 
  }

  else crash("MODE: "+MODE_POL+" in RatioPol not yet implemented");


  return VAL;
}

void crash(string Message) {
  cout << Message <<endl;
  exit(-1);
  return ;
}


double F(double x, double a, int alpha) {
  
  if(alpha==0) return a-1.0/cosh((alpha+1)*x);
  else if(alpha==-1) return a -cosh(alpha*x);
  
  double return_value= a -cosh(alpha*x)/cosh((alpha+1)*x);
  return return_value;
}


void derivative(Vfloat &RES, Vfloat &INPUT, string MODE) {

  int size= INPUT.size();

  for(int t=0; t<size;t++) {
    int t1,t2;
    if(MODE=="BACKWARD") {t2=t; t1= (t2-1)%size;}
    if(MODE=="FORWARD") {t1=t; t2= (t1+1)%size;}
    if(MODE=="SYMMETRIC") {t1= (t-1)%size; t2=(t+1)%size;}

    double temp_res= INPUT[t2]-INPUT[t1];
    if(MODE=="SYMMETRIC") temp_res /= 2.0;
    RES.push_back(temp_res);


  }

  return;
}


double FTAN(double x1,  int t, int NT) {

  double dNT = (double)NT/2.0;
  double dt = (double)t;
  
  double RES=  +1.0*(dNT -dt +1.0)*tanh(x1*(dNT-dt+1.0)) -1.0*(dNT -dt)*tanh(x1*(dNT-dt));

  // cout<<"FTAN: "<<RES<<endl<<flush;

  return RES;

}

double FTAN_SYMM(double x1,  int t, int NT) {

  double dNT = (double)NT/2.0;
  double dt = (double)t;
  
  double RES=  +1.0*(dNT -dt +1.0)*tanh(x1*(dNT-dt+1.0)) -1.0*(dNT -dt-1.0)*tanh(x1*(dNT-dt-1.0));

  // cout<<"FTAN: "<<RES<<endl<<flush;

  return RES;


}



double Root_Brent(double R, int nt, int NT) {


  
  
   
  //cout<<"R: "<<R<<" t:"<<nt<<" T: "<<NT<<endl;
  if( R < 1 && nt < NT/2) return log(fabs(R)) ;
  if( R > 1 && nt >= NT/2) return log(1.0/R);
  
  if(NT%2 != 0) crash("Temporal lattice extent is not divisible by two!");
  NT= NT/2; //it is what enter the expression of the effective mass
  int alpha = nt-NT;
  if(alpha==0 && R<= 1) return acosh(1.0/R);
  double Precision = 1e-9;
  double delta = 0.001;
  //solve the equation cosh(meff(alpha+1))*corr_ratio = cosh(meff*alpha) using the Brent method!
  //initialize iteration
  double b=0;
  double a;
  a= min(fabs(1./(alpha+1)),fabs(1./alpha));

  while (F(a,R,alpha)*F(b,R,alpha) > 0) {a  +=  min(fabs(1./(alpha+1)),fabs(1./alpha));}
 
 



  if(fabs(F(a,R,alpha)) < fabs(F(b,R,alpha))) {double atemp=a; a=b; b=atemp;}

  double c=a;
  bool FLAG = true;
  double s=b;
  double d=0;

  
  
  while(F(s,R,alpha) !=0 && fabs(b-a)>= Precision*fabs((double)(b+a)/2.0) ) {

  
    if((F(a,R,alpha) != F(c,R,alpha)) && (F(b,R,alpha) != F(c,R,alpha))) {//inverse quadratic interpolation
      s= a*F(b,R,alpha)*F(c,R,alpha)/((F(a,R,alpha)-F(b,R,alpha))*(F(a,R,alpha)-F(c,R,alpha))) + b*F(a,R,alpha)*F(c,R,alpha)/((F(b,R,alpha)-F(a,R,alpha))*(F(b,R,alpha)-F(c,R,alpha))) + c*F(a,R,alpha)*F(b,R,alpha)/((F(c,R,alpha)-F(a,R,alpha))*(F(c,R,alpha)-F(b,R,alpha)));
    }
    else s= b-(F(b,R,alpha)*(b-a)/(F(b,R,alpha)-F(a,R,alpha)));

    double s1= (double)(3*a+b/4.0);
    if( (s < s1 || s> b) || (FLAG==true && fabs(s-b) >= (double)fabs(b-c)/2) || (FLAG==false && fabs(s-b) >= (double)fabs(c-d)/2) || (FLAG==true && fabs(b-c) < delta) || (FLAG==false && fabs(c-d) < delta)) {
      
      FLAG=true;
      s= (a+b)/2.0;
      
    }
    
    else FLAG= false;

    d= c;
    c= b;
    if (F(a,R,alpha)*F(s,R,alpha)<0) b=s;
    else a=s;

    if(fabs(F(a,R,alpha)) < fabs(F(b,R,alpha))) {double atemp=a; a=b; b=atemp;}

   
  }

  
  return s;
}



double Root_Brent_sinh(double R, int nt, int NT) {

  if(nt==NT/2 -1 || nt==NT/2) return 0.0;
 
  if( R < 1 && nt < NT/2) return log(fabs(R)) ;
  if( R > 1 && nt >= NT/2) return log(1.0/R);
   
  if(NT%2 != 0) crash("Temporal lattice extent is not divisible by two!");
  NT= NT/2; //it is what enter the expression of the effective mass
  int alpha = nt-NT;
  
 
  

  auto F=[&](double x) {return R- sinh(x*alpha)/sinh(x*(alpha+1)); };
  double Precision = 1e-9;
  double delta = 0.001;
  //solve the equation sinh(meff(alpha+1))*corr_ratio = sinh(meff*alpha) using the Brent method!
  //initialize iteration
  double b=1e-8;
  double a;
  a= min(fabs(1./(alpha+1)),fabs(1./alpha));

  while (F(a)*F(b) > 0) {a  +=  min(fabs(1./(alpha+1)),fabs(1./alpha));}
 
 



  if(fabs(F(a)) < fabs(F(b))) {double atemp=a; a=b; b=atemp;}

  double c=a;
  bool FLAG = true;
  double s=b;
  double d=0;

  
  
  while(F(s) !=0 && fabs(b-a)>= Precision*fabs((double)(b+a)/2.0) ) {

  
    if((F(a) != F(c)) && (F(b) != F(c))) {//inverse quadratic interpolation
      s= a*F(b)*F(c)/((F(a)-F(b))*(F(a)-F(c))) + b*F(a)*F(c)/((F(b)-F(a))*(F(b)-F(c))) + c*F(a)*F(b)/((F(c)-F(a))*(F(c)-F(b)));
    }
    else s= b-(F(b)*(b-a)/(F(b)-F(a)));

    double s1= (double)(3*a+b/4.0);
    if( (s < s1 || s> b) || (FLAG==true && fabs(s-b) >= (double)fabs(b-c)/2) || (FLAG==false && fabs(s-b) >= (double)fabs(c-d)/2) || (FLAG==true && fabs(b-c) < delta) || (FLAG==false && fabs(c-d) < delta)) {
      
      FLAG=true;
      s= (a+b)/2.0;
      
    }
    
    else FLAG= false;

    d= c;
    c= b;
    if (F(a)*F(s)<0) b=s;
    else a=s;

    if(fabs(F(a)) < fabs(F(b))) {double atemp=a; a=b; b=atemp;}

   
  }

  
  return s;
}


double DoConstantFit(Vfloat &data, Vfloat &err) {

  double result=0;
  double weight=0;


  

  if (data.size() != err.size()) {
    cout<<"data size and err size do not match in double DoConstantFit(Vfloat&, Vfloat&)"<<endl;;
    cout<<"Printing info..."<<endl;
    printV(data, "data", 1);
    printV(err, "err", 1);
    crash("Abort.");
  }
  
  for(unsigned int k=0; k<data.size(); k++) {result += pow(1.0/err[k],2)*data[k]; weight += pow(1.0/err[k],2);}

  if(isnan(result) || isnan(weight)) { //fit without errors
    //cout<<"Switching from chi^2 minimization to least square minimization......"<<endl;
    result=0;
    weight=0;
    for(auto &meas: data) {result += meas; weight++;}
  }

  //normalize
  result /= weight;

  return result;
}

double lin_interpolator(double y1, double y2, double Dx1, double Dx2,double Dx) {

  double n1= y1*Dx2;
  double n2= y2*Dx1;

  double a0= (n1-n2)/(Dx2-Dx1);
  double a1 = (y1-a0)/Dx1;

  return a0+a1*Dx;
}


double quad_interpolator(double y1, double y2, double y3, double Dx1, double Dx2, double Dx3, double Dx) {


  
	double n1 = y1*Dx2*Dx3/( (Dx1-Dx2)*(Dx1-Dx3));
	double n2 = -y2*Dx1*Dx3/( (Dx1-Dx2)*(Dx2-Dx3));
	double n3 = y3*Dx1*Dx2/( (Dx1-Dx3)*(Dx2-Dx3));

	double prod= (Dx1-Dx2)*(Dx1-Dx3)*(Dx2-Dx3);

	double a0 = n1+n2+n3;
	double a1 = (pow(Dx3,2)*(y1-y2) + pow(Dx1,2)*(y2-y3) + pow(Dx2,2)*(y3-y1))/prod ;
	double a2 = (Dx3*(y2-y1) + Dx2*(y1-y3) + Dx1*(y3-y2))/prod;

	return a0 + a1*Dx + a2*pow(Dx,2);
  



}


void Print_To_File(const vector<string>& row_id, const VVfloat &data, string Path, string MODE, string Header) {

  ofstream Print;
  if(MODE=="APP") Print.open(Path, ofstream::app);
  else Print.open(Path, ofstream::out);

  Print.precision(16);

  if(data.size() == 0) crash("In Print_To_File an empty vectors has been provided");

  unsigned int rows= data[0].size();

  if((signed)row_id.size() != 0 && row_id.size() != rows) crash(" In Print_To_File data size and data_id size do not match"); 

  Print<<Header<<endl;
  for(unsigned int j=0; j < rows; j++) {

    Print<<j;
    if(row_id.size() != 0) Print<<setw(20)<<row_id[j];
    
    for(unsigned int i=0; i<data.size(); i++) {
      if(data[i].size() != rows) crash("In Print_To_File, number of measurements is non-constant over columns. Exiting...");
      Print<<setw(30)<<data[i][j];
    }
    Print<<endl;
  }
  Print.close();
  return;
}





void cascade_resize( vector<vector<double>>& arr, const Vint& A ) {

  if((signed)A.size() != 2) crash("In cascade_resize cannot resize vv<double> with a Vint of size != 2");
  arr.resize(A[0]);
  for( auto &a: arr) a.resize(A[1]);

  return;
  
}

void cascade_resize( vector<vector<vector<double>>>& arr,const Vint &A) {
 
  if((signed)A.size() != 3) crash("In cascade_resize cannot resize vvv<double> with a Vint of size != 3");
  arr.resize(A[0]);
  for(auto & a: arr) cascade_resize( a, slicing(A,1, A.size() -1));


  return;
}

void cascade_resize( vector<vector<vector<vector<double>>>>& arr,const Vint &A) {
  
  if((signed)A.size() != 4) crash("In cascade_resize cannot resize vvvv<double> with a Vint of size != 4");
  arr.resize(A[0]);
  for(auto & a: arr) cascade_resize( a, slicing(A,1, A.size() -1));


  return;
}


void cascade_resize( vector<vector<vector<vector<vector<double>>>>>& arr,const Vint &A) {

  if((signed)A.size() != 5) crash("In cascade_resize cannot resize vvvvv<double> with a Vint of size != 5");
  arr.resize(A[0]);
  for(auto & a: arr) cascade_resize( a, slicing(A,1, A.size() -1));


  return;
}
