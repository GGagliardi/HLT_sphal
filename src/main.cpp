#include "../include/Spectral.h"
#include "../include/ax_therm_rate.h"

using namespace std;


int main(int narg, char** argv)
{
  
  if(narg != 8 && narg != 5) {
    cout<<"Use ./HLT_spallate  m_T s_T NT Nboots correlator_path out_path INCLUDE_ERRORS"<<endl;
    cout<<"or use ./HLT_spallate m_T s_T NT out_path"<<endl;
    exit(-1);
  }
   
  if(narg==8) {

    double m_T = stod(argv[1]);
    double s_T = stod(argv[2]);
    int NT = stoi(argv[3]);
    int Nboots= stoi(argv[4]);
    string InputFile = argv[5];
    string out_path = argv[6];
    int INCLUDE_ERRORS= stoi(argv[7]);

    //read correlator
    int NThalf_p1= NT/2+1 ;
     



  
    //read correlator
    distr_t_list C(0, NT, Nboots);
    ifstream read(InputFile);
    int Nrows= Nboots*NThalf_p1;
    
    for(int irow=0;irow<Nrows;irow++) {
      
      int dummy;
      double dummy2,dummy3;
      double corr;
      
      read>>dummy>>dummy2>>corr>>dummy3;
      
      int t= irow%NThalf_p1;
      int boot= irow/NThalf_p1;
       
      C.distr_list[t].distr[boot] = corr;
      C.distr_list[(NT-t)%NT].distr[boot] = corr;
    }
    
    read.close();
    
    //change sign to C
    C=-1.0*C;
  

     //print correlator
     Print_To_File({}, {C.ave(), C.err()}, out_path+".corr", "", "");
    
    
     get_axion_therm_rate(m_T, s_T, NT, Nboots, C, out_path, INCLUDE_ERRORS);
  }

  else { //nargs=5

    double m_T = stod(argv[1]);
    double s_T = stod(argv[2]);
    int NT = stoi(argv[3]);
    string out_path = argv[4];
    int Nboots=1000;
    //generate fake correlator
    distr_t_list C = Get_id_distr_list(NT,Nboots,0);
    get_axion_therm_rate(m_T, s_T, NT, Nboots, C, out_path, -1);
    
    
  }
  
  return 0;
}











