#include <cmath>
#include <complex>
#include <vector>
#include <cstdio>
#include <iostream>
#include <cassert>
#include "PowerItoI.H"
#include "FFTCTRecursive.H"
#include "FFT1D.H"
#include "FFT1DRecursive.H"
FFT1DRecursive::FFT1DRecursive(const unsigned int& a_M):FFT1D(a_M)
{
};
FFT1DRecursive::~FFT1DRecursive(){};
void FFT1DRecursive::forwardFFTCC(vector<complex<double> > & a_fHat,
                            const vector<complex<double> >& a_f) const
{
    assert (a_fHat.size()==m_N);
    for (unsigned int j = 0; j < m_N; j++)
        {
        a_fHat[j]= a_f[j];
        }
    bool isForward = true;
    applyFFT(a_fHat, a_f, m_N);
}
void FFT1DRecursive::inverseFFTCC(vector<complex<double> >& a_f ,
                                   const vector<complex<double> > & a_fHat ) const
{
    assert (a_fHat.size()==m_N);
    for (unsigned int j = 0; j < m_N; j++)
        {
        a_f[j]=a_fHat[j];
        }
    bool isForward = false;
    applyFFTInv(a_f, a_fHat, m_N);
}
