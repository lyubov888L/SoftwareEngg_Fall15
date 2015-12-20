#include <cmath>
#include <complex>
#include <vector>
#include <cstdio>
#include <iostream>
using namespace std;
#include "PowerItoI.H"
#include "FFTCTRecursive.H"

void applyFFT(vector<complex<double> >& a_fhat, const vector<complex<double> >& a_f, int a_level)
{
    int N = a_f.size();
    a_fhat.resize(N);
    if (a_level > 2)
        {
        vector<complex<double> > fEven = evenPoints(a_f);
        vector<complex<double> > fOdd = oddPoints(a_f);
        N = a_f.size();
        a_fhat.resize(N);
        vector<complex<double> > fHatHalfEven, fHatHalfOdd;
        applyFFT(fHatHalfEven, fEven, a_level/2);
        applyFFT(fHatHalfOdd, fOdd, a_level/2);
        /// loop to compute fHatEven[k] + z^k*fHatOdd[k]
        for (int k = 0; k < a_level/2; k++)
            {
            complex<double> z(cos(-(2.0*M_PI*k)/float(a_level)), sin((-2.0*M_PI*k)/float(a_level)));
            a_fhat[k] = fHatHalfEven[k] + z*fHatHalfOdd[k];
            a_fhat[k+a_level/2] = fHatHalfEven[k] - z*fHatHalfOdd[k];
            }
       }
    else
        {
        a_fhat[1] = a_f[0] - a_f[1];
        a_fhat[0] = a_f[0] + a_f[1];
        }
    /// z = -1 at level 1
};

void applyFFTInv(vector<complex<double> >& a_f, const vector<complex<double> >& a_fhat, int a_level)
{
    int N = a_fhat.size();
    a_f.resize(N);
    if (a_level > 2)
        {
        vector<complex<double> > fHatEven = evenPoints(a_fhat);
        vector<complex<double> > fHatOdd = oddPoints(a_fhat);
        N = a_fhat.size();
        a_f.resize(N);
        vector<complex<double> > fHalfEven, fHalfOdd;
        applyFFTInv(fHalfEven, fHatEven, a_level/2);
        applyFFTInv(fHalfOdd, fHatOdd, a_level/2);
        /// loop to compute fHatEven[k] + z^k*fHatOdd[k]
        for (int k = 0; k < a_level/2; k++)
            {
            complex<double> z(cos((2.0*M_PI*k)/float(a_level)), sin((2.0*M_PI*k)/float(a_level)));
            a_f[k] = (fHalfEven[k] + z*fHalfOdd[k]);
            a_f[k+a_level/2] = (fHalfEven[k] - z*fHalfOdd[k]);
            }
       }
    else
        {
        a_f[1] = a_fhat[0] - a_fhat[1];
        a_f[0] = a_fhat[0] + a_fhat[1];
        }
    /// z = -1 at level 1
};



vector<complex<double> > evenPoints(const vector<complex<double> >& a_vec)
{
    vector<complex<double> > evenPointsVec;
    for(int i = 0; i < a_vec.size(); i++)
        {
        if (i%2 == 0)
            {
            evenPointsVec.push_back(a_vec[i]);
            }
        }
    return evenPointsVec;
};


vector<complex<double> > oddPoints(const vector<complex<double> >& a_vec)
{
    vector<complex<double> > oddPointsVec;
    for(int i = 0; i < a_vec.size(); i++)
        {
        if (i%2 != 0)
            {
            oddPointsVec.push_back(a_vec[i]);
            }
        }
    return oddPointsVec;
};
