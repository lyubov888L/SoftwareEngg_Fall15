#include "atom.h"

/*
 * Constructor given the positions and atom information. Velocities will be randomized
 */
CAtom::CAtom(CVector3<double> a_position, const int a_atomindex, const int a_atomtype){
   m_atomindex=a_atomindex; m_atomtype=a_atomtype;
   m_coord=a_position;
   m_force=0.0;  //set forces to 0.0
};

/*
 * Constructor given the positions and velocities and atom information
 */
CAtom::CAtom(CVector3<double> a_position, CVector3<double> a_vel,  const int a_atomindex, const int a_atomtype){
   m_atomindex=a_atomindex; m_atomtype=a_atomtype;
   m_coord=a_position;
   m_vel=a_vel;
   m_force=0.0;  //set forces to 0.0
};

/*
 * Velocities set to input velocity values
 */
void CAtom::setVelocity(CVector3<double>& a_vel){
   m_vel=a_vel;
}

/*
 * Forces set to input force values
 */
void CAtom::setForces(double a_force){
    m_force=a_force;
}

/*
 * Writes out all the properties of the atom, i.e. position, velocity, forces
 */
void CAtom::Display(){

   cout << "Atomtype:\t" << m_atomtype;
   cout <<"\nPosition:\t";
   m_coord.Display();
   cout <<"Velocity:\t";
   m_vel.Display();
   cout << "Force:\t";
   m_force.Display();
}
