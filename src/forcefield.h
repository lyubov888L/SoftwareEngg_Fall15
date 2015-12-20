/*!
 * The forcefield class holds all the force field parameters. 
 */


#ifndef __FORCEFIELD_H
#define __FORCEFIELD_H
#include <vector>
#include "config.h"
using namespace std;
class CForcefield {
 private:

     std::vector<double> m_mass/*!< Vector of masses for each atom type */, m_charge/*!< Vector of charges for each atom type */, m_sigma2/*!< Vector of Van der Waals sigma parameters for each atom type */, m_eps/*!< Vector of Van der Waals epsilon parameters for each atom type */;


 public:

    /* Initializer*/
    CForcefield(){}; /*!< Default constructor for force field class object */
    //Initializing from a parameter file
    CForcefield(const char*); /*!< Constructor for force field class object from input configuration file */
    void display(); /*!< Function to write out all the parameters */

    vector<double>& getMass(){return m_mass;} /*!< Returns vector of masses */
    const double getMassAtom(const int i) { return m_mass[i];} /*!< Returns the mass of specific atom type */

    const double getEps(const int a_value1, const int a_value2){
       if (a_value1 <a_value2){ return m_eps[a_value1*(atypes-1)+ a_value2]; }
       else { return m_eps[a_value2*(atypes-1)+ a_value1];}
    } /*!< Returns the combined epsilon value for a particular pair of atom types */

    const double getSigma2(const int a_value1, const int a_value2){
       if (a_value1 <a_value2){ return m_sigma2[a_value1*(atypes-1)+ a_value2]; }
       else { return m_sigma2[a_value2*(atypes-1)+ a_value1];}
    } /*!< Returns the combined sigma value for a particular pair of atom types */

    const double getCharges(const int i) {return m_charge[i];}  /*!< Returns the charge for a particular atom type */
};
#endif
