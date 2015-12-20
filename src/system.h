/*!
 * This is the entire system class, which contains all the molecules.
 */

#ifndef _CSYSTEM_H_
#define _CSYSTEM_H_

#include "parameters.h"
#include "forcefield.h"
#include "molecule.h"
#include "Point.H"
#include "Cubes.h"

class CSystem {

  private:
    int m_nmols/*!< Number of molecules */, m_totalatom /*!< Number of atoms in the system */;
    vector<CMolecule> m_mollist /*!< Vector of molecules in the system */;
    vector<CCubes> cubeSet;
    int m_ncubes;
  public:
  /* Constructors */
    CSystem(); /*!< Default constructor */
    CSystem(CParam&, CForcefield& ); /*!< Constructor based on the parameter and forcefield files */

  /*System initialization*/
    void readPosition(const char*); /*!< Reads positions from the input coordinate file */
    void createLattice(const double, const double, const double, const int); /*!< Creates a default lattice if no input file is specified */
    void readVelocity(const char*); /*!< Reads velocities from the input velocities file */
    void initializeVelocity(const double, vector<double>&); /*!< Initializes random velocities from a Gaussian if no input velocities file is specified */
    void setForcesZero(); /*!< Sets the forces on all atom to zero */
    void createCubes(CParam& a_params); /*!< Creates cubes for neighbor list assignment */
    void assignParticlesToCubes(CParam& a_params); /*!< Assigns particles/atoms to cubes */
    Point getPoint(int k, Point& m_lowCorner, Point& m_highCorner) const; /*!< Returns point by index */
    int getindex(const Point& a_pt, Point& m_lowCorner, Point& m_highCorner) const; /*!< Returns index by point */
    vector<CCubes>& getCubeSet() {return cubeSet;}; /*!< Returns the entire cubeset */


  /*Functions*/
    void fullDisplay(); /*!< Writes out the entire system */
    void display(); /*!< Writes out a note about the system information */
    void ShiftCoord(CForcefield&); /*!< Shifts the coordinates based to center the system */
    void GetCM(CVector3<double>&, vector<double>& a_mass); /*!< Calculates the coordinates of the system center of mass */

   /*writing to files*/
   void writeRestartvel(const char*); /*!< Writes velocities for restarting the system */
   void writeRestartpos(const char*); /*!< Writes coordinates for restarting the system */
   void writeMovie(ofstream&, int); /*!< Writes the trajectory for visualization */

  /* Destructor */
    ~CSystem(){}; /*!< Default destructor */

  /*get values*/
  const int getNmols(){return m_nmols;} /*!< Returns the number of molecules in the system */
  CMolecule& getMol(const int index){return m_mollist[index];} /*!< Returns a molecule based on index */
};

#endif
