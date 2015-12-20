#include "Nonbonded.h"
#include <cmath>

/*
 * This function calculates the Lennard-Jones force and Coulomb force between two atoms.
 */

void calcForce(CAtom& atom1, CAtom& atom2, CVector3<double>& forcevectorVDW, CVector3<double>& forcevectorCoulomb, CVector3<double>& boxsize, CVector3<double>& halfbox, bool pbcs_on, double cutoff, CForcefield& a_ff, CParam& a_prm){

    CVector3<double> rij; // Vector between two atoms
    int index1, index2; // Atom types
    double r2, s2r2, s6r6, erfFactor; // distance between atoms, sigma by r square, (sigma by r)^6, error function factor

    rij=atom1.getCoord() - atom2.getCoord();
    // Apply periodic box conditions
    if (pbcs_on){
        rij.minImage(halfbox, boxsize);
    }
    r2=rij.SqLength(); // Distance between atoms
    if ( !strcmp(a_prm.getMethod(),"PME"))
        {
        erfFactor = (1.0 - erf(sqrt(alphaForGamma*r2))); // Reduction factor for real space electrostatics when PME is applied
        }
    else if ( !strcmp(a_prm.getMethod(),"NBLIST"))
        {
        erfFactor = 1.0; // Error function is only used when using PME, hence set to 1.0 otherwise
        }
    else
        {
        cutoff = 100.0; // Very large cutoff when all electrostatics are to be calculated in real space
        erfFactor = 1.0; // Error function is only used when using PME, hence set to 1.0 otherwise
        }
    // Forces calcualted only within cutoff
    if (sqrt(r2) < cutoff and r2 > 0.00001){
    index1=atom1.getAtomtype();
    index2=atom2.getAtomtype();
    s2r2=a_ff.getSigma2(index1, index2)/r2; // sigma by rsquare
    s6r6=s2r2*s2r2*s2r2; // (sigma by r)^6
    forcevectorVDW=rij*(24*a_ff.getEps(index1, index2)*s6r6*(2*s6r6 -1)/r2); // Calculate VdW force

    forcevectorCoulomb=(rij/sqrt(r2)) * k_e * (a_ff.getCharges(index1)*a_ff.getCharges(index2))/(r2)*erfFactor; // Calculate Coulomb force
    }
};

void CNonbonded::calInterNonbondedNaive(CSystem& a_sys, CForcefield& a_ff, CParam& a_prm){

   CVector3<double> fij, boxsize(0.0), halfbox, forcevectorVDW(0.0), forcevectorCoulomb(0.0);
   double cutoff=a_prm.getCutoff();
   bool pbcs_on=a_prm.getPBCS();
   if (pbcs_on){
      boxsize[0]=a_prm.getCELLX();
      boxsize[1]=a_prm.getCELLY();
      boxsize[2]=a_prm.getCELLZ();
      halfbox=boxsize/2.0;
   }

   for(int mol1=0; mol1<a_sys.getNmols(); mol1++){
      for(int atm1=0; atm1<a_sys.getMol(mol1).getNatom(); atm1++){

          for(int mol2=mol1+1; mol2<a_sys.getNmols(); mol2++){
              for(int atm2=0; atm2<a_sys.getMol(mol2).getNatom(); atm2++){
                // Atoms between which forces are calculated
                  CAtom atom1 = a_sys.getMol(mol1).getAtom(atm1);
                  CAtom atom2 = a_sys.getMol(mol2).getAtom(atm2);
                  calcForce(atom1, atom2, forcevectorVDW, forcevectorCoulomb, boxsize, halfbox, pbcs_on, cutoff, a_ff, a_prm);

                  a_sys.getMol(mol1).getAtom(atm1).getForces()-=forcevectorVDW;
                  a_sys.getMol(mol2).getAtom(atm2).getForces()+=forcevectorVDW;

                  a_sys.getMol(mol1).getAtom(atm1).getForces()-=forcevectorCoulomb;
                  a_sys.getMol(mol2).getAtom(atm2).getForces()+=forcevectorCoulomb;

                  }
              }
          }
   }
};

void CNonbonded::calInterNonbondedNblist(CSystem& a_system, CForcefield& a_ff, CParam& a_params)
    {
    double factor = 0.5;

    CVector3<double> fij, boxsize(0.0), halfbox, forcevectorVDW(0.0), forcevectorCoulomb(0.0);
    double cutoff=a_params.getCutoff();
    bool pbcs_on=a_params.getPBCS();
    if (pbcs_on)
        {
        boxsize[0]=a_params.getCELLX();
        boxsize[1]=a_params.getCELLY();
        boxsize[2]=a_params.getCELLZ();
        halfbox=boxsize/2.0;
        }

    int sizeOfCubeSet = a_system.getCubeSet().size();
    for (int i = 0; i < sizeOfCubeSet; i++)
        {
        /*
         * This part accounts for periodic boundary conditions. The neighboring points are calculated
         * by taking the minimum image into account. If the center point is on the boundary, it's
         * neighboring point is calculated by wrapping around.
         */
        Point center = a_system.getCubeSet()[i].getCubeCenter();
        int centerLinIndex = a_system.getindex(center, a_system.getCubeSet()[0].getCubeCenter(), a_system.getCubeSet()[sizeOfCubeSet-1].getCubeCenter());
        for (int j = 0; j < a_params.getDirectionUnitVector().size(); j++)
            {
            Point tempNbPoint = center + a_params.getDirectionUnitVector()[j];
            int tempTuple[DIM]; // Temp point to get image
            for(int l = 0; l < DIM; l++)
                {
                tempTuple[l] = 0;
                if (tempNbPoint[l] < a_system.getCubeSet()[0].getCubeCenter()[l])
                    {
                    tempTuple[l] = a_system.getCubeSet()[sizeOfCubeSet - 1].getCubeCenter()[l];
                    }
                else if (tempNbPoint[l] > a_system.getCubeSet()[sizeOfCubeSet - 1].getCubeCenter()[l])
                    {
                    tempTuple[l] = a_system.getCubeSet()[0].getCubeCenter()[l];
                    }
                else
                    {
                    tempTuple[l] = tempNbPoint[l];
                    }
                }
            Point nbPoint(tempTuple); // Point created using actual neighbor or image, whichever is closer
            int nbLinIndex = a_system.getindex(nbPoint, a_system.getCubeSet()[0].getCubeCenter(), a_system.getCubeSet()[sizeOfCubeSet-1].getCubeCenter());
            /*
             * Looping over all molecule pairs within neighboring cubes only
             */
            for(int mol1=0; mol1 < a_system.getCubeSet()[centerLinIndex].getInteriorMolecule().size(); mol1++)
                {
                for(int atm1=0; atm1 < a_system.getCubeSet()[centerLinIndex].getInteriorMolecule()[mol1].getNatom(); atm1++)
                    {
                    for(int mol2=0; mol2 < a_system.getCubeSet()[nbLinIndex].getInteriorMolecule().size(); mol2++)
                        {
                        for(int atm2=0; atm2 < a_system.getCubeSet()[nbLinIndex].getInteriorMolecule()[mol2].getNatom(); atm2++)
                            {
                            /*
                             * If the two atoms are the same, the force should not be calculated,
                             * i.e. factor = 1, else, the factor = 0.5 to account for double counting
                             * of atoms.
                             */
                            if (a_system.getCubeSet()[centerLinIndex].getInteriorMolecule()[mol1].getMolIndex() == a_system.getCubeSet()[nbLinIndex].getInteriorMolecule()[mol2].getMolIndex())
                                {
                                factor = 0.0;
                                }
                            else
                                {
                                factor = 0.5;
                                }
                            CAtom atom1 = a_system.getCubeSet()[centerLinIndex].getInteriorMolecule()[mol1].getAtom(atm1);
                            CAtom atom2 = a_system.getCubeSet()[nbLinIndex].getInteriorMolecule()[mol2].getAtom(atm2);
                            calcForce(atom1, atom2, forcevectorVDW, forcevectorCoulomb, boxsize, halfbox, pbcs_on, cutoff, a_ff, a_params);
                            // The force vectors are multiplied by the factor
                            forcevectorVDW*=factor;
                            forcevectorCoulomb*=factor;

                            a_system.getMol(a_system.getCubeSet()[centerLinIndex].getInteriorMolecule()[mol1].getMolIndex()).getAtom(atm1).getForces()-=forcevectorVDW;
                            a_system.getMol(a_system.getCubeSet()[nbLinIndex].getInteriorMolecule()[mol2].getMolIndex()).getAtom(atm2).getForces()+=forcevectorVDW;

                            a_system.getMol(a_system.getCubeSet()[centerLinIndex].getInteriorMolecule()[mol1].getMolIndex()).getAtom(atm1).getForces()-=forcevectorCoulomb;
                            a_system.getMol(a_system.getCubeSet()[nbLinIndex].getInteriorMolecule()[mol2].getMolIndex()).getAtom(atm2).getForces()+=forcevectorCoulomb;
                            }
                        }
                    }
                }
            }
        }
    };


/*
 * Calculating forces by the PME method
 */
void CNonbonded::calInterNonbondedPME(CSystem& a_system, CForcefield& a_ff, CParam& a_params)
    {
    // Calling the pairwise for the less than cutoff calculations.
    calInterNonbondedNblist(a_system, a_ff, a_params);

    CVector3<double> forcevectorCoulomb(0.0);

    // Declare required variables.
    double lx = a_params.getCELLX();
    double ly = a_params.getCELLY();
    double lz = a_params.getCELLZ();

    double cutoff = a_params.getCutoff(); // Cutoff

    int Nx = ceil(lx/(cutoff/gridFactor)); // High corner for grid
    int Ny = ceil(ly/(cutoff/gridFactor));
    int Nz = ceil(lz/(cutoff/gridFactor));

    /*
     * Setting high and low corner
     */
    int gridPoint[DIM];
    int pointIndices[DIM];
    pointIndices[0] = 0;
    pointIndices[1] = 0;
    pointIndices[2] = 0;
    Point m_lowCorner(pointIndices);
    pointIndices[0] = Nx-1;
    pointIndices[1] = Ny-1;
    pointIndices[2] = Nz-1;
    Point m_highCorner(pointIndices);
    double s[DIM];

    double lby2[DIM];
    lby2[0] = a_params.getCELLX()/2.0;
    lby2[1] = a_params.getCELLY()/2.0;
    lby2[2] = a_params.getCELLZ()/2.0;

    /*
     * Unit vectors for accessing neighboring grid points.
     */
    int v001[DIM];
    int v010[DIM];
    int v011[DIM];
    int v100[DIM];
    int v101[DIM];
    int v110[DIM];
    int v111[DIM];
    v001[0] = 0; v001[1] = 0, v001[2] = 1;
    v010[0] = 0; v010[1] = 1, v010[2] = 0;
    v011[0] = 0; v011[1] = 1, v011[2] = 1;
    v100[0] = 1; v100[1] = 0, v100[2] = 0;
    v101[0] = 1; v101[1] = 0, v101[2] = 1;
    v110[0] = 1; v110[1] = 1, v110[2] = 0;
    v111[0] = 1; v111[1] = 1, v111[2] = 1;
    Point e001(v001);
    Point e010(v010);
    Point e011(v011);
    Point e100(v100);
    Point e101(v101);
    Point e110(v110);
    Point e111(v111);
    vector<Point> unitPointsVector(2*DIM);
    vector<Point> imagePointsVector(2*DIM);
    vector<Point> unitPointsVector2(pow(2,DIM));
    vector<Point> imagePointsVector2(pow(2,DIM));
    unitPointsVector[0] = e100;
    unitPointsVector[1] = e100*-1;
    unitPointsVector[2] = e010;
    unitPointsVector[3] = e010*-1;
    unitPointsVector[4] = e001;
    unitPointsVector[5] = e001*-1;

    unitPointsVector2[0] = e001;
    unitPointsVector2[1] = e010;
    unitPointsVector2[2] = e011;
    unitPointsVector2[3] = e100;
    unitPointsVector2[4] = e101;
    unitPointsVector2[5] = e110;
    unitPointsVector2[6] = e111;

    double tempGridSize = cutoff/gridFactor;

    /// Declaring variables requred for the 4 step ewald protocol.
    Box omega_box(m_lowCorner, m_highCorner);
    RectMDArray<complex<double> > omega(omega_box);

    vector<vector<double> > Ugrid;
    vector<vector<double> > Uparticles;
    int numGridPoints = omega.getBox().sizeOf();
    int numParticles = a_system.getNmols();

    /// Initialize the omegaGrid values to zero.
    for (int i = 0; i < numGridPoints; i++)
        {
        omega[i] = complex<double>(0.0, 0.0);
        }

    /// Deposit the charges on the grid
    for (int i = 0; i < a_system.getNmols(); i++){
        for (int j = 0; j < a_system.getMol(i).getNatom(); j++){
            for (int k = 0; k < DIM; k++){
                pointIndices[k] = floor((a_system.getMol(i).getAtom(j).getCoord()[k] + lby2[k])/tempGridSize);
                s[k] = ((a_system.getMol(i).getAtom(j).getCoord()[k] - pointIndices[k] * tempGridSize) / tempGridSize);
            }
        Point tempPoint(pointIndices);
        omega[tempPoint] += complex<double>(a_ff.getCharges(a_system.getMol(i).getAtom(j).getAtomtype()) * (1 - s[0]) * (1 - s[1]) * (1 - s[2]), 0.0);
        omega[tempPoint+e001] += complex<double>(a_ff.getCharges(a_system.getMol(i).getAtom(j).getAtomtype()) * (1 - s[0]) * (1 - s[1]) * (s[2]), 0.0);
        omega[tempPoint+e010] += complex<double>(a_ff.getCharges(a_system.getMol(i).getAtom(j).getAtomtype()) * (1 - s[0]) * (s[1]) * (1 - s[2]), 0.0);
        omega[tempPoint+e011] += complex<double>(a_ff.getCharges(a_system.getMol(i).getAtom(j).getAtomtype()) * (1 - s[0]) * (s[1]) * (s[2]), 0.0);
        omega[tempPoint+e100] += complex<double>(a_ff.getCharges(a_system.getMol(i).getAtom(j).getAtomtype()) * (s[0]) * (1 - s[1]) * (1 - s[2]), 0.0);
        omega[tempPoint+e101] += complex<double>(a_ff.getCharges(a_system.getMol(i).getAtom(j).getAtomtype()) * (s[0]) * (1 - s[1]) * (s[2]), 0.0);
        omega[tempPoint+e110] += complex<double>(a_ff.getCharges(a_system.getMol(i).getAtom(j).getAtomtype()) * (s[0]) * (s[1]) * (1 - s[2]), 0.0);
        omega[tempPoint+e111] += complex<double>(a_ff.getCharges(a_system.getMol(i).getAtom(j).getAtomtype()) * (s[0]) * (s[1]) * (s[2]), 0.0);
        }
    }

    /// Performing a DIM dimensional fourier transform, to obtain the potential on the grid.
    shared_ptr<FFT1D> p_fft = dynamic_pointer_cast<FFT1D >(shared_ptr<FFTW1D>(new FFTW1D(numGridPoints)));
    FFTMD fftmd(p_fft);
    fftmd.forwardCC(omega);
    for (int k = 0; k < numGridPoints; k++)
        {
            if (k != 0)
                {
                omega[k] = (4*M_PI/(k*k)) * (exp((-(k*k))/(4*alphaForGamma))) * omega[k];
                }
            else
                {
                omega[k] = complex<double>(0.0,0.0);
                }
        }
    fftmd.inverseCC(omega);

    /// Create the grid field container (vector<vector<Real>> in this case) and initialize to zero.
    Ugrid.resize(numGridPoints);
    for (int i = 0; i < numGridPoints; i++)
        {
        Ugrid[i].resize(DIM);
        Ugrid[i][0] = 0.0;
        Ugrid[i][1] = 0.0;
        Ugrid[i][2] = 0.0;
        }

    /// Calculate the fields on the grid points via finite differences.
    Box box = omega.getBox();
    for (int i = 0; i < numGridPoints; i++)
        {
        gridPoint[0] = omega.getBox().getPoint(i)[0];
        gridPoint[1] = omega.getBox().getPoint(i)[1];
        gridPoint[2] = omega.getBox().getPoint(i)[2];
        Point p(gridPoint);

    for (int m = 0; m < unitPointsVector.size(); m++)
            {
            Point tempNbPoint = p + unitPointsVector[m];
            int tempTuple[DIM];
            for(int l = 0; l < DIM; l++)
                {
                tempTuple[l] = 0;
                if (tempNbPoint[l] < box.getLowCorner()[l])
                    {
                    tempTuple[l] = tempNbPoint[l] + (box.getHighCorner()[l] - box.getLowCorner()[l]);
                    }
                else if (tempNbPoint[l] > box.getHighCorner()[l])
                    {
                    tempTuple[l] = tempNbPoint[l] - (box.getHighCorner()[l] - box.getLowCorner()[l]);
                    }
                else
                    {
                    tempTuple[l] = tempNbPoint[l];
                    }
                }
            Point newNbPoint(tempTuple);
            imagePointsVector[m] = newNbPoint;
            }

        /// Compute field only for interior grid points.
        Ugrid[i][0] = -((real(omega[imagePointsVector[0]]) - real(omega[imagePointsVector[1]])) / (2*tempGridSize));
        Ugrid[i][1] = -((real(omega[imagePointsVector[2]]) - real(omega[imagePointsVector[3]])) / (2*tempGridSize));
        Ugrid[i][2] = -((real(omega[imagePointsVector[4]]) - real(omega[imagePointsVector[5]])) / (2*tempGridSize));
        }

    /// Create the particles field container (vector<vector<Real>> in this case) and initialize to zero.
    Uparticles.resize(numParticles);
    for (int i = 0; i < numParticles; i++)
        {
        Uparticles[i].resize(DIM);
        Uparticles[i][0] = 0.0;
        Uparticles[i][1] = 0.0;
        Uparticles[i][2] = 0.0;
        }

    for (int i = 0; i < a_system.getNmols(); i++)
        {
        for (int j = 0; j < a_system.getMol(i).getNatom(); j++)
            {
            for (int k = 0; k < DIM; k++)
                {
                pointIndices[k] = floor((a_system.getMol(i).getAtom(j).getCoord()[k] + lby2[k])/tempGridSize);
                s[k] = ((a_system.getMol(i).getAtom(j).getCoord()[k] - pointIndices[k] * tempGridSize) / tempGridSize);
                }
            Point p(pointIndices);
            /*
             * This part accounts for periodic boundary conditions. The neighboring points are calculated
             * by taking the minimum image into account. If the point p is on the boundary, it's
             * neighboring point is calculated by wrapping around.
             */
            for (int m = 0; m < unitPointsVector2.size(); m++)
                {
                Point tempNbPoint = p + unitPointsVector2[m];
                int tempTuple[DIM];
                for(int l = 0; l < DIM; l++)
                    {
                    tempTuple[l] = 0;
                    if (tempNbPoint[l] < box.getLowCorner()[l])
                        {
                        tempTuple[l] = tempNbPoint[l] + (box.getHighCorner()[l] - box.getLowCorner()[l]);
                        }
                    else if (tempNbPoint[l] > box.getHighCorner()[l])
                        {
                        tempTuple[l] = tempNbPoint[l] - (box.getHighCorner()[l] - box.getLowCorner()[l]);
                        }
                    else
                        {
                        tempTuple[l] = tempNbPoint[l];
                        }
                    }
                Point newNbPoint(tempTuple);
                imagePointsVector2[m] = newNbPoint; // Vector containing nearest image points
                }
            // Calculating the field on the particles from the field at the grid points
            Uparticles[i][0] = Ugrid[omega.getBox().getIndex(p)][0] * (1 - s[0]) * (1 - s[1]) * (1 - s[2])
                            + Ugrid[omega.getBox().getIndex(imagePointsVector2[0])][0] * (1 - s[0]) * (1 - s[1]) * (s[2])
                            + Ugrid[omega.getBox().getIndex(imagePointsVector2[1])][0] * (1 - s[0]) * (s[1]) * (1 - s[2])
                            + Ugrid[omega.getBox().getIndex(imagePointsVector2[2])][0] * (1 - s[0]) * (s[1]) * (s[2])
                            + Ugrid[omega.getBox().getIndex(imagePointsVector2[3])][0] * (s[0]) * (1 - s[1]) * (1 - s[2])
                            + Ugrid[omega.getBox().getIndex(imagePointsVector2[4])][0] * (s[0]) * (1 - s[1]) * (s[2])
                            + Ugrid[omega.getBox().getIndex(imagePointsVector2[5])][0] * (s[0]) * (s[1]) * (1 - s[2])
                            + Ugrid[omega.getBox().getIndex(imagePointsVector2[6])][0] * (s[0]) * (s[1]) * (s[2]);
            Uparticles[i][1] = Ugrid[omega.getBox().getIndex(p)][1] * (1 - s[0]) * (1 - s[1]) * (1 - s[2])
                            + Ugrid[omega.getBox().getIndex(imagePointsVector2[0])][1] * (1 - s[0]) * (1 - s[1]) * (s[2])
                            + Ugrid[omega.getBox().getIndex(imagePointsVector2[1])][1] * (1 - s[0]) * (s[1]) * (1 - s[2])
                            + Ugrid[omega.getBox().getIndex(imagePointsVector2[2])][1] * (1 - s[0]) * (s[1]) * (s[2])
                            + Ugrid[omega.getBox().getIndex(imagePointsVector2[3])][1] * (s[0]) * (1 - s[1]) * (1 - s[2])
                            + Ugrid[omega.getBox().getIndex(imagePointsVector2[4])][1] * (s[0]) * (1 - s[1]) * (s[2])
                            + Ugrid[omega.getBox().getIndex(imagePointsVector2[5])][1] * (s[0]) * (s[1]) * (1 - s[2])
                            + Ugrid[omega.getBox().getIndex(imagePointsVector2[6])][1] * (s[0]) * (s[1]) * (s[2]);
            Uparticles[i][2] = Ugrid[omega.getBox().getIndex(p)][2] * (1 - s[0]) * (1 - s[1]) * (1 - s[2])
                            + Ugrid[omega.getBox().getIndex(imagePointsVector2[0])][2] * (1 - s[0]) * (1 - s[1]) * (s[2])
                            + Ugrid[omega.getBox().getIndex(imagePointsVector2[1])][2] * (1 - s[0]) * (s[1]) * (1 - s[2])
                            + Ugrid[omega.getBox().getIndex(imagePointsVector2[2])][2] * (1 - s[0]) * (s[1]) * (s[2])
                            + Ugrid[omega.getBox().getIndex(imagePointsVector2[3])][2] * (s[0]) * (1 - s[1]) * (1 - s[2])
                            + Ugrid[omega.getBox().getIndex(imagePointsVector2[4])][2] * (s[0]) * (1 - s[1]) * (s[2])
                            + Ugrid[omega.getBox().getIndex(imagePointsVector2[5])][2] * (s[0]) * (s[1]) * (1 - s[2])
                            + Ugrid[omega.getBox().getIndex(imagePointsVector2[6])][2] * (s[0]) * (s[1]) * (s[2]);
            // Calculating the force on the particle from the field
            forcevectorCoulomb[0] = Uparticles[i][0]* k_e * a_ff.getCharges(a_system.getMol(i).getAtom(j).getAtomtype());
            forcevectorCoulomb[1] = Uparticles[i][1]* k_e * a_ff.getCharges(a_system.getMol(i).getAtom(j).getAtomtype());
            forcevectorCoulomb[2] = Uparticles[i][2]* k_e * a_ff.getCharges(a_system.getMol(i).getAtom(j).getAtomtype());
            //Updating the forces on each atom
            a_system.getMol(i).getAtom(j).getForces()-=forcevectorCoulomb;
            }
        }
    };
