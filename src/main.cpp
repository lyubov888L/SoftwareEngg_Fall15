#include <iostream>
#include <fstream>
#include <cmath>
#include "moves.h"

using namespace std;
void parseflags(int, char**, char*);
void flaglist();
int main(int argc, char *argv[])
    {
    char controlfile[50];

    // Reading the control file name and movie flag
    parseflags(argc, argv, controlfile);

    /* Read control file and write parameters feedback to screen */
    CParam runparams(controlfile);
    runparams.createDirectionUnitVector();

    /*Initialize the forcefield*/
    CForcefield ff(runparams.getParmfile());

    /* Initialize the system configuration */
    CSystem sys(runparams, ff);
    sys.ShiftCoord(ff);
    sys.createCubes(runparams);
    sys.assignParticlesToCubes(runparams);

    // Set up the dynamics loop here.
    CMoves move;
    move.propagate ( sys , ff, runparams);
//    sys.fullDisplay();

    // Exit
    cout << "\n\nSimulation completed successfully!! \n\n" ;
    return 0;
    }

void parseflags(int argc, char* argv[], char* controlfile){

  /* this function is there so that later on more flags can be added if needed*/
  if( argc > 2 ) {
   for(int i=1; i<argc; i++) {
      if( *(argv[i]) == '-' ) {
         if ( *(argv[i]+1) =='p') {
            strcpy(controlfile, argv[i+1]);
         }
         else {
             cerr << " Unknown flag " << argv[i] << endl << endl;
             flaglist();
         }
      }
    }
  }
  else { flaglist();};
  return;

};

void flaglist() {
  cout << endl
       << " Usage: "
       << " [-p parameter file]"
       << endl << endl;
  cout << " Supported flags: " << endl << endl;
  cout << "  p - ";
  cout << "read specified parameter file \n";
  cout << endl;

  exit(1);
}

