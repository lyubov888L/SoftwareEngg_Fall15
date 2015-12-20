/*!
 * This class performs the dynamic update step, and calculates the updated positions and velocities.
 */

#ifndef __CMOVES_H
#define __CMOVES_H

#include "Nonbonded.h"

class CMoves{
public:

   CMoves(){}; /*!< Default constructor for moves class object */

   ~CMoves(){};  /*!< Default destructor for moves class object */

   void propagate(CSystem&, CForcefield&, CParam&); /*!< This class performs the dynamic update using the Velocity Verlet Algorithm */

   void calForce(CSystem&, CForcefield&, CParam&, int timeStep); /*!< This class calls the appropriate force calculation method based on the keyfile argument */

   void assignCubes(CSystem&){}; /*!< This class assigns cubes for the neighbor list calculation */
};

#endif
