#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "parameters.h"

using namespace std;

/* Constructors */
/*
 * Default constructor
 */
CParam::CParam(){

  //Periodic boundary parameters
  m_pbcs=true;
  m_cell_x=0;m_cell_y=0;m_cell_z=0;

  /*simulation initialization*/
  m_restartvel=false; m_readpos=false; m_movieflag=false;
  m_kinstep=0; m_equilstep=0; m_moviefilefreq=0; m_restartfreq=0;m_mdoutfreq=0;
  m_natom_eachside=0;

  m_cutoff=100;

  m_dt=0; m_T=0;

  strcpy(m_method,"PAIRWISE");
}

/* initialization from a given control filename*/
CParam::CParam(char* a_fname){

  //Periodic boundary parameters
  m_pbcs= true;
  m_cell_x=0;m_cell_y=0;m_cell_z=0;

  /*simulation initialization*/
  m_restartvel=false; m_readpos=false; m_movieflag=false;
  m_kinstep=0; m_equilstep=0; m_moviefilefreq=0; m_restartfreq=0;m_mdoutfreq=0;
  m_natom_eachside=0;

  m_cutoff=100;

  m_dt=0; m_T=0;

  strcpy(m_method,"PAIRWISE");

  char      tmp[5000];

  *this = CParam(); //Basic Initialization
  ifstream  file;

  file.open(a_fname, ios::in);
  if( !file ) {
    cout << endl << "Error opening file " << a_fname << "." << endl;
    exit(1);
  }

  cout << endl << " Reading parameter file " << a_fname << "..." ;

  file >> tmp;
  while( !file.eof() ) {
    if( tmp[0] == '#' ) {
      file.getline(tmp, 300, '\n');
    }
    else if( !strcmp(tmp, "TEMP") ) {
      file >> tmp;
      m_T = atof(tmp);
    }
    else if( !strcmp(tmp, "PBCON") ) {
      m_pbcs = true;
    }
    else if( !strcmp(tmp, "PBCOFF") ) {
      m_pbcs = false;
    }
    else if( !strcmp(tmp, "CELLX") ) {
      file >> tmp;
      m_cell_x = atof(tmp);
    }
    else if( !strcmp(tmp, "CELLY") ) {
      file >> tmp;
      m_cell_y = atof(tmp);
    }
    else if( !strcmp(tmp, "CELLZ") ) {
      file >> tmp;
      m_cell_z = atof(tmp);
    }
    else if( !strcmp(tmp, "TIMESTEP") ) {
      file >> tmp;
      m_dt = atof(tmp);
    }
    else if( !strcmp(tmp, "KINSTEP") ) {
      file >> tmp;
      m_kinstep = atoi(tmp);
    }
    else if( !strcmp(tmp, "EQUILSTEP") ) {
      file >> tmp;
      m_equilstep = atoi(tmp);
    }
    else if( !strcmp(tmp, "RESTARTVELON") ) {
      m_restartvel=true;
    }
    else if( !strcmp(tmp, "VELFNAME") ) {
      file >> tmp;
      strcpy(m_readvelfname, tmp);
    }
    else if( !strcmp(tmp, "READPOSON") ) {
      m_readpos=true;
    }
    else if( !strcmp(tmp, "NATOMEACHSIDE") ) {
      file >> tmp;
      m_natom_eachside = atoi(tmp);
    }
    else if( !strcmp(tmp, "POSFNAME") ) {
      file >> tmp;
      strcpy(m_readposfname, tmp);
    }
    else if( !strcmp(tmp, "MOVIEFLAGON") ) {
      m_movieflag=true;
    }
    else if( !strcmp(tmp, "MOVIEFNAME") ) {
      file >> tmp;
      strcpy(m_moviefname, tmp);
    }
    else if( !strcmp(tmp, "MOVIEFILEFREQ") ) {
      file >> tmp;
      m_moviefilefreq = atoi(tmp);
    }
    else if( !strcmp(tmp, "MDOUTFNAME") ) {
      file >> tmp;
      strcpy(m_mdoutfname, tmp);
    }
    else if( !strcmp(tmp, "MDOUTFREQ") ) {
      file >> tmp;
      m_mdoutfreq = atoi(tmp);
    }
    else if( !strcmp(tmp, "RESTARTVELFNAME") ) {
      file >> tmp;
      strcpy(m_restartvelfname, tmp);
    }
    else if( !strcmp(tmp, "RESTARTPOSFNAME") ) {
      file >> tmp;
      strcpy(m_restartposfname, tmp);
    }
    else if( !strcmp(tmp, "RESTARTFREQ") ) {
      file >> tmp;
      m_restartfreq = atoi(tmp);
    }
    else if( !strcmp(tmp, "PARMFILE") ) {
      file >> tmp;
      strcpy(m_parmfname, tmp);
    }
    else if( !strcmp(tmp, "CUTOFF") ) {
      file >> tmp;
      m_cutoff=atof(tmp);
    }
    else if( !strcmp(tmp, "METHOD") ) {
      file >> tmp;
      strcpy(m_method, tmp);
    }
    else {
      cout << endl << " Unrecognized parameter '" << tmp
           << "' in input file " << a_fname << "."
           << "  Skipping..." << endl;
      file.getline(tmp, 300, '\n');
    }
    file >> tmp;
  }

  file.close();

  cout << "  done reading control file" << endl;

};

/*
 * Writes out the parameters
 */
void CParam::Display(){

  cout << endl << "--------------------------------------------------------"
    << endl;
  /* Move information */
  cout << "will add more information here later"<<endl;
  printf("Periodic boundary conditions: ");
  if(m_pbcs) {
    printf("ON\n");
    printf("Cell lengths (a,b,c): %8.5f  %8.5f  %8.5f\n",m_cell_x,m_cell_y,m_cell_z);
  } else {
    printf("OFF\n");
  }
  fflush(stdout);

  if (!strcmp(m_method,"PME")){
   cout <<"Particle mesh ewald method is being used.\n";
   cout << "PME cutoff:\t" << m_cutoff ;
  }
  else{
   cout <<"Particle mesh ewald method is not used.\n";
  }

  cout << "--------------------------------------------------------"
    << endl << endl;
};


/*
 * This initializes the direction unit vector for calculating the neigbouring cubes for the nblist
 */
void CParam::createDirectionUnitVector(){
     int i = 0;
     for (int x = -1; x < 2; x++){
         for (int y = -1; y < 2; y++){
             for (int z = -1; z < 2; z++){
                 int temp[DIM];
                 temp[0] = x; temp[1] = y; temp[2] = z;
                 Point P(temp);
                 directionUnitVector.push_back(P);
                 i++;
                 }
             }
         }
     };

