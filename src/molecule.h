/*!
 * The molecule class holds all atoms within a molecule.
 */

#ifndef _CMOLECULE_H_
#define _CMOLECULE_H_

#include <vector>
#include "atom.h"
class CMolecule{

private:
   int m_natoms; /*!< Number of atoms in the molecule */
   int m_molIndex; /*!< Global index of the molecule */
   std::vector<CAtom> m_atomlist; /*!< Vector of atom objects */

public:
    /* Constructors */
    CMolecule(){}; /*!<  */
    CMolecule(const int a_natoms, vector<double> a_tuple, vector<int> a_atomindex, int a_moleculeIndex, vector<int> a_atomtype); /*!<  */

    /*Functions*/
    void Display(); /*!< Writes out the molecule */
    void ShiftCoord(CVector3<double>&); /*!< Shifts the coordinates of the molecules by the vector */

    CAtom& operator[](const int a_index){ return m_atomlist[a_index];} /*!< Returns a reference to an atom within the molecule by index */
    CAtom& getAtom(const int a_index){ return m_atomlist[a_index];} /*!< Returns a reference to an atom within the molecule by index */
    CAtom getActualAtom(const int a_index){ return m_atomlist[a_index];} /*!< Returns the actual atom object within the molecule by index  */
    const int getNatom(){ return m_natoms;} /*!< Returns the number of atoms within the molecule */
    const int getMolIndex() {return m_molIndex;}; /*!< Returns the global index of the molecule */
    void setMolIndex(int index) { m_molIndex = index;}; /*!< Sets the global index of the molecule */

};


#endif
