#include "moves.h"
#include <cmath>
#include <fstream>

using namespace std;

/*
 * Propagte performs the actual simulation, updating the positions, velocities and forces
 */
void CMoves::propagate(CSystem& a_sys, CForcefield& a_ff, CParam& a_prm){

   int kinstep=a_prm.getKinstep(); // number of steps
   double dt=a_prm.getDt(); // timestep
   double dtby2=dt*0.5; // half timestep
   ofstream movie; // stream for writing out trajectory

   bool pbcs_on=a_prm.getPBCS(); // periodic boundary conditions flag
   CVector3<double> boxsize(0.0), halfbox(0.0); // box dimensions and half box dimensions
   // get box dimensions
   if (pbcs_on){
      boxsize[0]=a_prm.getCELLX();
      boxsize[1]=a_prm.getCELLY();
      boxsize[2]=a_prm.getCELLZ();
      halfbox=boxsize/2.0;
   }

   /*Setting up movie file*/
   if (a_prm.getMovieflag()){
        movie.open(a_prm.getMoviefname(), ios::out);
        cout << "\nOpening " <<a_prm.getMoviefname()<< " to save movie.\n";
        if (!movie){
            cout << "Error opening .xyz movie file " << a_prm.getMoviefname() << "." << endl;
            exit(-1);
        }
   }
   printf( "\nDoing dynamics run for %d timesteps using a step size %+.0e seconds....", kinstep, dt) ;
   CVector3<double> vel_halfStep(0.0), coord_halfStep(0.0);
   calForce(a_sys, a_ff, a_prm, 0); // Initialize force
   for (int i =0; i<kinstep; i++){
          cout << "Step number = " << i << endl;
          /* Get new positions, and half-step velocities */
          for(int mol=0; mol<a_sys.getNmols(); mol++){
              for(int atm=0; atm<a_sys.getMol(mol).getNatom(); atm++){
                     vel_halfStep=a_sys.getMol(mol).getAtom(atm).getForces()*(dtby2/a_ff.getMassAtom(a_sys.getMol(mol).getAtom(atm).getAtomtype()));

                     coord_halfStep=vel_halfStep+a_sys.getMol(mol).getAtom(atm).getVelocity();
                     coord_halfStep*=dt;
                     /*update half step vel*/
                     a_sys.getMol(mol).getAtom(atm).getVelocity()+=vel_halfStep;

                     /*update new position*/
                     a_sys.getMol(mol).getAtom(atm).getCoord()+=coord_halfStep;
                     if (pbcs_on){
                         a_sys.getMol(mol).getAtom(atm).getCoord().minImage(halfbox, boxsize);
		          }
              }
          }
          /*get new forces at t+dt*/
          calForce(a_sys, a_ff, a_prm, i);

          /* finish updating half-step velocities */
          for(int mol=0; mol<a_sys.getNmols(); mol++){
              for(int atm=0; atm<a_sys.getMol(mol).getNatom(); atm++){
                     vel_halfStep=a_sys.getMol(mol).getAtom(atm).getForces()*(dtby2/a_ff.getMassAtom(a_sys.getMol(mol).getAtom(atm).getAtomtype()));
                     a_sys.getMol(mol).getAtom(atm).getVelocity()+=vel_halfStep;
              }
          }

          /* write restart file if needed*/
          if (i%a_prm.getRestartfreq()==0){
             a_sys.writeRestartvel(a_prm.getRestartVelfname());
             a_sys.writeRestartpos(a_prm.getRestartPosfname());
          /* write movie file if needed */
             if (a_prm.getMovieflag() && i%a_prm.getMoviefilefreq() == 0){
                a_sys.writeMovie(movie, i);
             }
          }
   }

   if (a_prm.getMovieflag()){
      movie.close();
   }

   cout << "done :)\n" << endl;
   cout << "Storing final configuration and velocities.\n";
   if (kinstep%a_prm.getRestartfreq()==0){
      a_sys.writeRestartvel(a_prm.getRestartVelfname());
      a_sys.writeRestartpos(a_prm.getRestartPosfname());
   }
};

/*
 * Calls the appropriate force calculation function based on the method
 */
void CMoves::calForce(CSystem& a_sys, CForcefield& a_ff, CParam& a_prm, int timeStep){
    /*
     * Particles are assigned to cubes (for neighborlist implementation) at this interval
     */
    if (timeStep%reassignCubesFreq == 0){
        a_sys.assignParticlesToCubes(a_prm);
        }

   CNonbonded nbded; // Initializing a nonbonded class object
   a_sys.setForcesZero(); // Initializing forces to zero before performing force calculation based on specified method
   if ( !strcmp(a_prm.getMethod(),"PME")){
        nbded.calInterNonbondedPME(a_sys, a_ff, a_prm); // Force calculation via PME method
    }
   else if ( !strcmp(a_prm.getMethod(),"NBLIST")) {
        nbded.calInterNonbondedNblist(a_sys, a_ff, a_prm); // Force calculation via NB list. This function call here is only for testing. The actual application of the NB list force call is from within the PME calculation.
   }
   else {
        nbded.calInterNonbondedNaive(a_sys, a_ff, a_prm); // Force calculation via naive pairwise calculation.
   }
};
