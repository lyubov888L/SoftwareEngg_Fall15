/*!
 * This class calculates the non bonded forces on the atoms.
 */

#ifndef __VDW_H
#define __VDW_H

#include "system.h"
#include <complex>
#include <vector>
#include <memory>
#include "Box.H"
#include "RectMDArray.H"
#include "FFTMD.H"
#include "FFTW1D.H"
#include "FFT1D.H"
#include "FFT1DBRI.H"
#include "fftw3.h"

class CNonbonded{

public:

   /*Functions*/
   void calInterNonbondedNaive(CSystem&, CForcefield&, CParam&); /*!< Calculates the non bonded forces using direct pairwise computation */
   void calInterNonbondedNblist(CSystem&, CForcefield&, CParam&); /*!< Calculates the non bonded forces using a neighbor list approach, has to be used in conjunction with the PME method */
   void calInterNonbondedPME(CSystem&, CForcefield&, CParam&); /*!< Calculates the non bonded forces using the PME method. The real space forces are calculated via a neighbor list */

};

#endif
