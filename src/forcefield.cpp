#include <fstream>
#include <iostream>
#include "forcefield.h"
#include <cmath>

using namespace std;
int combination(const int a_int);
/*
 * Default constructor for forcefield
 */
CForcefield::CForcefield(const char* a_prm){
  char      tmp[5000];

  ifstream  file;

  file.open(a_prm, ios::in);
  if( !file ) {
    cout << endl << "Error opening file " << a_prm << "." << endl;
    exit(1);
  }

  cout << endl << " Reading forcefield file " << a_prm << "..." ;

  file >> tmp;

  int size=combination(atypes);
  m_mass.resize(atypes,0.0); m_charge.resize(atypes,0.0); m_sigma2.resize(size,0.0); m_eps.resize(size,0.0);

  /*this part needs to be changed if more atoms are included later*/
  while( !file.eof() ) {
      if( tmp[0] == '#' ) {
        file.getline(tmp, 300, '\n');
         cout << tmp <<"\t";
      }
      /*
       * Reads the parameters for sodium
       */
      else if( !strcmp(tmp, "Na") ) {
        file >> tmp;
        m_mass[Na]=atof(tmp); //in 10-26
        file >> tmp;
        m_sigma2[Na*atypes]=atof(tmp);  // in nm2
        file >> tmp;
        m_eps[Na*atypes]=atof(tmp);     // in kJ/mol
        file >> tmp;
        m_charge[Na]=atof(tmp);           //unit charge
      }
      /*
       * Reads the parameters for chloride
       */
      else if( !strcmp(tmp, "Cl") ) {
        file >> tmp;
        m_mass[Cl]=atof(tmp);
        file >> tmp;
        m_sigma2[Cl*atypes]=atof(tmp);
        file >> tmp;
        m_eps[Cl*atypes]=atof(tmp);
        file >> tmp;
        m_charge[Cl]=atof(tmp);
     }
      file >> tmp;
  }

  file.close();

  /*
   * Computes the combined sigma and epsilon for atom pairs
   */
  for(int i=0; i<atypes; i++){
     for(int j=i+1; j<atypes; j++){
          m_eps[i*atypes-i+j]=sqrt(m_eps[i*atypes]* m_eps[j*atypes]);
          m_sigma2[i*atypes-i+j]=0.5*(m_sigma2[i*atypes] +  m_sigma2[j*atypes]);
     }
  }
  for(int i=0; i<size; i++){
      m_eps[i]/=(Navo*1000);                    // This factor of 1000 comes from converting from kJ to J
      m_sigma2[i]*=m_sigma2[i];
  }

  cout << "  done initializing forcefield" << endl;

};

/*
 * Writes out the forcefield
 */
void CForcefield::display(){

  cout << "forcefield parameters:" << endl;
  cout << "atomtype\tmass\t\teps\t\tsigma\t\tcharge\n";
  for(int i=0; i<atypes; i++){
    if (i==0){cout <<"Na"<<"\t\t";}
    else if (i==1){cout <<"Cl"<<"\t\t";}
    cout  << m_mass[i] << "\t\t"<<sqrt(m_eps[i*atypes])<<"\t\t" <<  m_sigma2[i*atypes] <<"\t\t" <<m_charge[i]<<"\n";
  }

};

/*
 * Returns the combination value
 */
int combination(const int a_int){

   if (a_int==1) { return 1;}
   else{
      return (a_int+combination(a_int-1));
   }
};
