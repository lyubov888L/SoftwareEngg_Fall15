/*!
 * The atom class holds all the properties of an atom.
 */

#ifndef _CATOM_H_
#define _CATOM_H_

#include "vector3.h"

class CAtom{

private:
   int m_atomindex/*!< Index of the atom */, m_atomtype/*!< Atom type */;
   CVector3<double> m_coord; /*!< Vector of atomic coordinates */
   CVector3<double> m_vel; /*!< Vector of atomic velocities */
   CVector3<double> m_force; /*!< Vector of atomic forces */


public:
    /* Constructors */
    CAtom(){}; /*!< Default constructor for the atom */
    CAtom(CVector3<double> a_position, const int a_atomindex, const int a_atomtype); /*!< Atom constructor given the positions and atom information. Velocity will be randomized. */
    CAtom(CVector3<double> a_position, CVector3<double> a_vel,  const int a_atomindex, const int a_atomtype);  /*!< Atom constructor given both the positions and velocities, as well as the atom information */

    /*Functions*/
    void Display(); /*!< Writes out all the properties of the atom, i.e. position, velocity, forces*/
    void setVelocity(CVector3<double>& a_vel); /*!< Sets the velocity of the atom to the the input velocity vector */
    void setForces(double a_force); /*!< Sets the forves in all directions to the input value */
    CVector3<double>& getVelocity(){ return m_vel;} /*!< Returns the velocity vector of the atom */
    CVector3<double>& getCoord(){ return m_coord;} /*!< Returns the coordinates vector of the atom */
    CVector3<double>& getForces(){ return m_force;} /*!< Returns the forces vector of the atom */

    /* get types*/
    const int getAtomtype(){ return m_atomtype;} /*!< Returns the atom type */
    const int getAtomIndex(){ return m_atomindex;} /*!< Returns the atom index */

};

#endif
