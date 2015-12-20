
#include <iostream>
#include <fstream>
//#include <cmath>
#include <stdio.h>
#include <assert.h>
#include "config.h"
#include "system.h"
#include <random>
#include <chrono>
using namespace std;

/*
 * Default constructor
 */
CSystem::CSystem(){
m_nmols=0; m_totalatom=0;
};

/*
 * Constructor given the parameter and force field files
 */
CSystem::CSystem(CParam& a_param, CForcefield& a_ff){

  //reading initial position file
  if (a_param.getReadpos()){
    readPosition(a_param.getReadPosfname());
  }
  //creating own lattice
  else{
    if (a_param.getNatomEachside() < 1){
      cout << "Tag NATOMEACHSIDE needs to be mentioned in control parameter file.\n"; exit (-1);
    };
    createLattice(a_param.getCELLX(),a_param.getCELLY(),a_param.getCELLZ(), a_param.getNatomEachside());
  }
  //reading initial velocity from file
  if (a_param.getRestartVel()){
    readVelocity(a_param.getReadVelfname());
  }
  //initializing the system based on given temperature
  else{
    initializeVelocity(a_param.getTemp(), a_ff.getMass());
  }
};

/*
 * Read positions from input xyz file
 */
void CSystem::readPosition(const char* a_posfname){
    ifstream posfile;
    char      tmp[1000];
    vector<double> temp_pos(DIM,0.0);
    /* Open the system configuration file */
    posfile.open(a_posfname, ios::in);
    if(!posfile) {
      cout << endl << " Error: cannot open file '" << a_posfname << "'" << endl;
      exit(-1);
    }

    printf("\nReading system configuration from '%s'............\n",a_posfname);
    fflush(stdout);

    posfile >> tmp ;
    m_totalatom=atoi(tmp);
    m_nmols=m_totalatom; //later on it ca be extended so that one molecule can contain more than one atom
    posfile >> tmp ;
    int atomno=0;
    /* Reading atom positions*/
    posfile >> tmp;
    while( !posfile.eof() ) {

       if (!strcmp(tmp,"Na")){
         posfile >> tmp; temp_pos[0]=atof(tmp);
         posfile >> tmp; temp_pos[1]=atof(tmp);
         posfile >> tmp; temp_pos[2]=atof(tmp);
         vector<int> atomindex(1, atomno);
         vector<int> atomname(1, Na);
         m_mollist.push_back(CMolecule(1, temp_pos, atomindex, atomno, atomname));
         atomno++;
       }
       else if (!strcmp(tmp,"Cl")){
         posfile >> tmp; temp_pos[0]=atof(tmp);
         posfile >> tmp; temp_pos[1]=atof(tmp);
         posfile >> tmp; temp_pos[2]=atof(tmp);
         vector<int> atomindex(1, atomno);
         vector<int> atomname(1, Cl);
         m_mollist.push_back(CMolecule(1, temp_pos, atomindex, atomno, atomname));
         atomno++;

       }
       else{
         printf("Atom type '%s' not found.\n Currently Code supports only Na and Cl.", tmp);
       }
       posfile >> tmp;
    };

    posfile.close();

    /* checking whether the number of atoms matches the XYZ header line*/
    if ( atomno <  m_totalatom) {
      cout << "XYZ file does not have all the atoms mentioned in the first line of XYZ file"  << endl; exit(-1);
    }
    else if ( atomno >  m_totalatom) {
      cout << "XYZ file has more atoms than the number mentioned in the first line of XYZ file"  << endl; exit(-1);
    }

    cout << "  done reading initial position file" << endl;

};

/*
 * Creates lattice if input xyz file is not given
 */
void CSystem::createLattice(const double a_cellx, const double a_celly, const double a_cellz, const int a_natom_eachside){

    double xSpacing = a_cellx / a_natom_eachside ; // x spacing between two adjacent atoms.
    double ySpacing = a_celly / a_natom_eachside ; // y spacing between two adjacent atoms.
    double zSpacing = a_cellz / a_natom_eachside ; // z spacing between two adjacent atoms.
    vector<double> temp_pos(DIM,0.0);
    int atomno=0; //running index of atom numbers
    printf("Creating lattice.............\n");
    for (int i = 0; i < a_natom_eachside; i++)
        {
        for (int j = 0; j < a_natom_eachside; j++)
            {
            for (int k = 0; k < a_natom_eachside; k++)
                {
                temp_pos[0] = i * xSpacing;
                temp_pos[1] = j * ySpacing;
                temp_pos[2] = k * zSpacing;
                vector<int> atomindex(1, atomno);
                if ((i+j+k)%2 == 0)
                    {
                       vector<int> atomname(1, Na);
                       m_mollist.push_back(CMolecule(1, temp_pos, atomindex, atomno, atomname));
                    }
                else
                    {
                       vector<int> atomname(1, Cl);
                       m_mollist.push_back(CMolecule(1, temp_pos, atomindex, atomno, atomname));
                    }
                atomno++;
                }
            }
        }
    m_totalatom=atomno;
    m_nmols=m_totalatom;
    printf("Created lattice containing %d atoms.\n", atomno);

};

/*
 * Read in velocity from input file
 */
void CSystem::readVelocity(const char* a_velfname){

    ifstream velfile;
    char      tmp[1000];
    CVector3<double> temp_vel(0.0);
    /* Open the system configuration file */
    velfile.open(a_velfname, ios::in);
    if(!velfile) {
      cout << endl << " Error: cannot open file '" << a_velfname << "'" << endl;
      exit(-1);
    }

    printf("Reading initlal velocity from '%s'............\n",a_velfname);
    fflush(stdout);

    int atomno=0;
    for(int i=0; i<m_nmols; i++){

      for(int j=0; j<m_mollist[i].getNatom(); j++){
    /* Reading initial velocity*/
           if ( !velfile.eof() ) {
              for (int k=0; k<DIM ; k++){
                velfile >> tmp; temp_vel[k]=atof(tmp);
              }
              m_mollist[i].getAtom(j).setVelocity(temp_vel);
              atomno++;
           }
           else {
              cout << "Total number of velocity vectors in the restart file ( "<< atomno <<") does not match the number of atoms in the system (" << m_totalatom << " )\n" ;
              exit(-1);
           }
      };

    };

    velfile.close();
    if (atomno!=m_totalatom){
       cout << "Total number of velocity vectors in the restart file ( "<< atomno <<") does not match the number of atoms in the system (" << m_totalatom << " )\n" ;
    }

};

/*
 * Initialize velocity with random numbers from Gaussian distribution at input temperature
 */
void CSystem::initializeVelocity(const double a_Temp, vector<double>& a_mass){

   CVector3<double> temp_vel(0.0);
   cout << "\nInitializing velocities at temp " <<a_Temp << ".........\n";
   vector<double> scalingfactor(atypes, 0.0);
   for (int i =0; i<atypes; i++){
       scalingfactor[i]=sqrt(kB*0.001*a_Temp/a_mass[i]);
   }
   cout << scalingfactor[0] <<"\t"<<scalingfactor[1]<<"\n";

   // obtain a seed from the timer
   double mean =0;
   double stddev=1;

   typedef std::chrono::high_resolution_clock myclock;
   myclock::time_point beginning = myclock::now();

   myclock::duration d = myclock::now() - beginning;
   unsigned seed2 = d.count();
//   unsigned seed2 = 1;

   default_random_engine generator(seed2);
   normal_distribution<double> distribution(mean, stddev);

   CVector3<double> velSum(0.0);
   for(int i=0; i<m_nmols; i++){
      for(int j=0; j<m_mollist[i].getNatom(); j++){
            /* Reading initial velocity*/

            for (int k=0; k<DIM ; k++){
              temp_vel[k]=distribution(generator);
            }
            temp_vel*=scalingfactor[m_mollist[i].getAtom(j).getAtomtype()];
            m_mollist[i].getAtom(j).setVelocity(temp_vel);
            velSum+=temp_vel;
      }
   }
   velSum/=1.0*m_totalatom;

   /*Substracting vel sum, so that the net vel is zero*/
   for(int i=0; i<m_nmols; i++){
      for(int j=0; j<m_mollist[i].getNatom(); j++){
          m_mollist[i].getAtom(j).getVelocity()-=velSum;
          m_mollist[i].getAtom(j).getVelocity()/=1000;   //converting to nm/fs
      }
   }

   cout << "done\n";
};

/*
 * Sets all forces to zero
 */
void CSystem::setForcesZero(){
    for(int i=0; i<m_nmols; i++){
      for(int j=0; j<m_mollist[i].getNatom(); j++){
         m_mollist[i].getAtom(j).setForces(0.0);
      }
    }
}

/*
 * Creates cubes for neighbor list implementation
 */
void CSystem::createCubes(CParam& a_params){

    double lx = a_params.getCELLX();
    double ly = a_params.getCELLY();
    double lz = a_params.getCELLZ();

    double cutoff = a_params.getCutoff();

    int Nx = ceil(lx/(cutoff/cubeFactor));
    int Ny = ceil(ly/(cutoff/cubeFactor));
    int Nz = ceil(lz/(cutoff/cubeFactor));

    m_ncubes = Nx*Ny*Nz;
    cubeSet.resize(m_ncubes);

    int pointIndices[DIM];
    pointIndices[0] = 0;
    pointIndices[1] = 0;
    pointIndices[2] = 0;
    Point m_lowCorner(pointIndices);
    pointIndices[0] = Nx-1;
    pointIndices[1] = Ny-1;
    pointIndices[2] = Nz-1;
    Point m_highCorner(pointIndices);

    for (int i = 0; i < m_ncubes; i++){
        cubeSet[i].getCubeCenter() = getPoint(i, m_lowCorner, m_highCorner);
        cubeSet[i].setCubeSize(cutoff/cubeFactor);
    }

};

/*
 * Assigns particles to cubes
 */
void CSystem::assignParticlesToCubes(CParam& a_params){

    double lx = a_params.getCELLX();
    double ly = a_params.getCELLY();
    double lz = a_params.getCELLZ();

    double cutoff = a_params.getCutoff();

    int Nx = ceil(lx/(cutoff/cubeFactor));
    int Ny = ceil(ly/(cutoff/cubeFactor));
    int Nz = ceil(lz/(cutoff/cubeFactor));

    int pointIndices[DIM];
    pointIndices[0] = 0;
    pointIndices[1] = 0;
    pointIndices[2] = 0;
    Point m_lowCorner(pointIndices);
    pointIndices[0] = Nx-1;
    pointIndices[1] = Ny-1;
    pointIndices[2] = Nz-1;
    Point m_highCorner(pointIndices);

    double lby2[DIM];
    lby2[0] = a_params.getCELLX()/2.0;
    lby2[1] = a_params.getCELLY()/2.0;
    lby2[2] = a_params.getCELLZ()/2.0;

    double tempCubeSize = cutoff/cubeFactor;

    for (int i = 0; i < m_nmols; i++){
        for (int j = 0; j < m_mollist[i].getNatom(); j++){
            for (int k = 0; k < DIM; k++){
                pointIndices[k] = floor((m_mollist[i].getAtom(j).getCoord()[k] + lby2[k])/tempCubeSize);
            }
        Point tempPoint(pointIndices);
        cubeSet[getindex(tempPoint, m_lowCorner, m_highCorner)].getInteriorMolecule().push_back(m_mollist[i]);
        }
    }

    for (int i = 0; i < cubeSet.size(); i++){
        for (int j =0; j < cubeSet[i].getInteriorMolecule().size(); j++){
            }
        }
}

/*
 * Gets point from index
 */
Point CSystem::getPoint(int k, Point& m_lowCorner, Point& m_highCorner) const
{
  assert(k >= 0);
  int tuple[DIM];
  for (int i = 0;i < DIM; i++)
    {
      int factor = (m_highCorner[i] - m_lowCorner[i] + 1);
      int kred = k%factor;
      tuple[i] = kred + m_lowCorner[i];
      k = (k - kred)/factor;
    }
  Point pt(tuple);
  return pt;
};

/*
 * Gets index from point
 */
int CSystem::getindex(const Point& a_pt, Point& m_lowCorner, Point& m_highCorner) const
{
  int factor = 1;
  int linindex = a_pt[0] - m_lowCorner[0];
  for (int i = 1;i < DIM;i++)
    {
      factor = factor*(m_highCorner[i-1] - m_lowCorner[i-1]+1);
      linindex = linindex + (a_pt[i] - m_lowCorner[i])*factor;
    }
  return linindex;
};

void CSystem::display(){

    cout << "System Information:\n" <<"Total number of molecules: \t" << m_nmols <<"\nTotal number of atoms:  \t" <<m_totalatom <<endl;


};

/*
 * Writes out the whole system
 */
void CSystem::fullDisplay(){

    cout << "System Information:\n" <<"Total number of molecules: \t" << m_nmols <<"\nTotal number of atoms:  \t" <<m_totalatom <<endl;
    for (int i=0; i<m_nmols; i++){
      cout << "Molecule " <<i <<" :\n";
      m_mollist[i].Display();
    }
}

/*
 * Shifts coordinates of atoms
 */
void CSystem::ShiftCoord(CForcefield& a_ff)
{
  CVector3<double>   center(0.0);

  GetCM(center, a_ff.getMass());

  for(int i=0; i<m_nmols; i++){
    m_mollist[i].ShiftCoord(center);
  }
  GetCM(center, a_ff.getMass());
}

/*
 * Gets center of mass of system
 */
void CSystem::GetCM(CVector3<double>& a_center, vector<double>& a_mass){

    double mass_sum = 0.0;
    CVector3<double> dummy(0.0);
    a_center*=0;
    for(int i=0; i<m_nmols; i++){
       for(int j=0; j<m_mollist[i].getNatom(); j++){
           //a_center+=m_mollist[i].getAtom(j).getCoord() * a_mass[m_mollist[i].getAtom(j).getAtomtype()];
           dummy=m_mollist[i].getAtom(j).getCoord() * a_mass[m_mollist[i].getAtom(j).getAtomtype()];
           a_center+=dummy;
           mass_sum+=a_mass[m_mollist[i].getAtom(j).getAtomtype()];
      }
    }
    a_center/=mass_sum;
    return ;
}

/*
 * Writes the velocities for restarting
 */
void CSystem::writeRestartvel(const char* a_velfname){

    ofstream velfile;
    velfile.open(a_velfname, ios::out);

    for(int i=0; i<m_nmols; i++){
       for(int j=0; j<m_mollist[i].getNatom(); j++){

         for (int dim=0; dim<DIM; dim++){
            velfile << m_mollist[i].getAtom(j).getVelocity()[dim] <<"\t\t";
         }
         velfile <<endl;
       }
    }
   velfile.close();
};

/*
 * Writes out coordinates for restarting
 */
void CSystem::writeRestartpos(const char* a_posfname){

    ofstream posfile;
    posfile.open(a_posfname, ios::out);

    posfile << "\t" <<m_totalatom <<"\n";
    posfile <<"\tXYZ\n";
    for(int i=0; i<m_nmols; i++){
       for(int j=0; j<m_mollist[i].getNatom(); j++){
           if (m_mollist[i].getAtom(j).getAtomtype() == 0)
            {
            posfile << "Na" << "\t";
            }
           else
            {
            posfile << "Cl" << "\t";
            }
         for (int dim=0; dim<DIM; dim++){
            posfile << m_mollist[i].getAtom(j).getCoord()[dim] <<"\t\t";
         }
         posfile <<endl;
       }
    }
   posfile.close();
};

/*
 * Writes out the trajectory
 */
void CSystem::writeMovie(ofstream& a_file, int a_time){

    /* Make sure the file is open*/
    if(!a_file) {
         cout << endl << " passed .xyz file is not open" << endl;
         exit(1);
    }
    a_file << "\t" <<m_totalatom <<"\n";
    a_file <<"\ttime" <<a_time<<"\n";
    for(int i=0; i<m_nmols; i++){
       for(int j=0; j<m_mollist[i].getNatom(); j++){
           if (m_mollist[i].getAtom(j).getAtomtype() == 0)
            {
            a_file << "Na" << "\t";
            }
           else
            {
            a_file << "Cl" << "\t";
            }
         for (int dim=0; dim<DIM; dim++){
            a_file << m_mollist[i].getAtom(j).getCoord()[dim] <<"\t\t";
         }
         a_file <<endl;
       }
    }
}

