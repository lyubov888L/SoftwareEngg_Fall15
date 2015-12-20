/*!
 * This class contains the parameters requiered for the simulation
 */

#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

#include <iostream>
#include <vector>
#include "config.h"
#include "Point.H"

using namespace std;

class CParam {
  private:


  /* simulation cell parameters */
  bool     m_pbcs;  /*!< Boolean flag for periodic boundary conditions */
  double   m_cell_x,m_cell_y,m_cell_z; /*!< Box dimensions */

  /*simulation initialization*/
  bool m_restartvel /*!< Boolean flag to read in velocities */, m_readpos/*!< Boolean flag to read in positions */, m_movieflag/*!< Boolean flag to write out movie trajectory */;
  char m_readposfname[50]/*!< Name of the initial coordinates file */, m_readvelfname[50]/*!< Name of the initial velocities file */, m_restartposfname[50] /*!< Name of the restarts coordinates file */, m_restartvelfname[50] /*!< Name of the restart velocities file */, m_moviefname[50] /*!< Name of the trajectory file */, m_mdoutfname[50]/*!< Name of the mdout data file */;
  int m_kinstep/*!< Number of dynamic steps */, m_equilstep/*!< Number of equilibration steps */, m_moviefilefreq/*!< Frequency of writing out trajectory file */, m_mdoutfreq/*!< Frequency of writing out the mdout data file */, m_restartfreq/*!< Frequency of writing out restart file */, m_natom_eachside/*!< Number of atoms on each side of the box for initialization */;


  /* forcefield parameters */
  char m_parmfname[50]; /*!< Name of the parameters file */
  double m_cutoff; /*!< Cutoff for VdW and Ewald calculations */
  char m_method[50]; /*!< Name of the method for electrostatic computations */


  /* move parameters */
  double   m_dt/*!< Time step for dynamics */, m_T/*!< Temperature */;

  vector<Point> directionUnitVector; /*!< Direction unit vector for calculating the neigbouring cubes for the nblist */

  public:
  /* Constructors */
  //Basic Initialization
  CParam(); /*!< Default constructor */
  //Initialization from a source file
  CParam(char*); /*!< Constructor from input parameter file */
  /* Destructor */
  ~CParam() {}; /*!< Destructor */

  /*Display the values*/
  void Display(); /*!< Writes out the parameters */

  /*return values*/
  const int getKinstep(){return m_kinstep;} /*!< Return number of dynamics steps */
  const double getDt(){ return m_dt;} /*!< Return timestep */
  const bool getReadpos(){return m_readpos;} /*!< Return read positions flag */
  const char* getReadPosfname(){return m_readposfname;} /*!< Return positions filename */
  const bool getRestartVel(){return m_restartvel;} /*!< Return read velocities flag */
  const char* getReadVelfname(){return m_readvelfname;} /*!< Return velocities filename */
  const double getCELLX(){return m_cell_x;} /*!< Return cell x-dimension */
  const double getCELLY(){return m_cell_y;} /*!< Return cell x-dimension */
  const double getCELLZ(){return m_cell_z;} /*!< Return cell x-dimension */
  const int getNatomEachside(){return m_natom_eachside;} /*!< Return number of atoms along each side of box */
  const double getTemp() {return m_T;} /*!< Return simulation temperature */
  const char* getParmfile(){return m_parmfname;} /*!< Return parameters file name */
  const double getCutoff(){return m_cutoff;} /*!< Return Van der Waal's and Ewald cutoff */
  const bool getPBCS() { return m_pbcs;} /*!< Return periodic boundary conditions flag */
  const int getRestartfreq() { return m_restartfreq;} /*!< Returns frequency for writing restart file */
  const int getMdoutfreq() { return m_mdoutfreq;} /*!< Returns frequency for mdout data file */
  const int getMoviefilefreq() { return m_moviefilefreq;} /*!< Returns frequency for writing out movie file */
  const bool getMovieflag() { return m_movieflag;} /*!< Returns flag for writing out movie file */
  const char* getMoviefname() { return m_moviefname;} /*!< Returns name of movie file */
  const char* getRestartPosfname() { return m_restartposfname;} /*!< Returns name of restart positions file */
  const char* getRestartVelfname() { return m_restartvelfname;} /*!< Returns name of restart velocities file */
  const char* getMethod() { return m_method;} /*!< Returns electrostatic computation method cutoff */

  void createDirectionUnitVector(); /*!< Initializes the direction unit vector for calculating the neigbouring cubes for the nblist */
  vector<Point>& getDirectionUnitVector() {return directionUnitVector;}; /*!< Returns direction unit vector for calculating neigbbouring cubes for the nblist */
};

#endif
