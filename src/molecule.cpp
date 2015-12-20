#include "molecule.h"
#include <assert.h>

/*
 * Constructor given input values
 */
CMolecule::CMolecule(const int a_natoms, vector<double> a_tuple, vector<int> a_atomindex, int a_moleculeIndex, vector<int> a_atomtype){

   m_molIndex = a_moleculeIndex;
   m_natoms=a_natoms;
   assert(a_tuple.size() == a_natoms*DIM);
   m_atomlist.resize(m_natoms);
   for (int i =0; i<a_natoms; i++){
      CVector3<double> a_position(0.0);
      for (int j =0; j<DIM; j++){
         a_position[j]=a_tuple[i*DIM+j];
      }
      m_atomlist[i]=CAtom( a_position, a_atomindex[i], a_atomtype[i]);
   }
};

/*
 * Write out the molecule
 */
void CMolecule::Display(){

  cout << "Number of atoms:\t" << m_natoms <<"\n";
  for(int i =0; i<m_natoms;i++){

     m_atomlist[i].Display();
  }
}

/*
 * Shift coordinates
 */
void CMolecule::ShiftCoord(CVector3<double>& a_center){

   for (int i=0; i<m_natoms;i++){
      m_atomlist[i].getCoord()-=a_center;
   }
}
