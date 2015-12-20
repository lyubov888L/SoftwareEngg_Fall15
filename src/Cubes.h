/*!
 * This class is an element to store atoms to be able to use neighborlists to calculate interactions within cutoffs.
 */
#ifndef _CCUBES_H_
#define _CCUBES_H_

#include <iostream>
#include "Point.H"
#include "atom.h"

class CCubes{

    private:

        float m_cubeSize; /*!< Dimension of the cube */
        Point m_cubeCenter; /*!< Center of the cube */
        vector<CMolecule> m_interiorMolecules; /*!< Vector of interior molecules contained within a cube */

    public:

        double getCubeSize() {return m_cubeSize;}; /*!< Return size of the cube */
        void setCubeSize(double a_cubeSize) {m_cubeSize = a_cubeSize;}; /*!< Set the size of the cube */
        Point& getCubeCenter() {return m_cubeCenter;}; /*!< Return the center of the cube */
        vector<CMolecule>& getInteriorMolecule() {return m_interiorMolecules;}; /*!< Returns the vector of interior molecules contained in the cube */

};

#endif
