#include <iostream>
using namespace std;

/*
 * Constructor with input value
 */
template <class X> CVector3<X>::CVector3(const X a_value){
    m_tuple.resize(DIM);
    for (int i = 0; i < DIM; i++){
        m_tuple[i] = a_value;
    }
};

/*
 * Access operator to access the desired element from the internal tuple
 */
template <class X> X& CVector3<X>::operator[](const int a_index) {
  return m_tuple[a_index];
};

/*
 * Writes out the object tuple
 */
template <class X> void CVector3<X>::Display(){
  for (int i=0; i<DIM; i++){
    cout << m_tuple[i] << "\t";
  }
  cout <<endl;
}

/*
 * Copy constructor
 */
template <class X> void CVector3<X>::copy(CVector3<X>& a_vector){
   for(int i =0; i<DIM; i++){
   m_tuple[i]=a_vector[i];
   };
};

/*
 * Internal Vector Addition Operator
 */
template <class X> void CVector3<X>::operator+=( CVector3& a_v){
   for(int i =0; i<DIM; i++){
   m_tuple[i]+=a_v[i];
   };
};

/*
 * Internal Vector Subtraction Operator
 */
template <class X> void CVector3<X>::operator-=( CVector3& a_v){
   for(int i =0; i<DIM; i++){
   m_tuple[i]-=a_v[i];
   };
};

/*
 * Internal Scalar Division Operator
 */
template <class X> void CVector3<X>::operator/=(const double a){
   for(int i =0; i<DIM; i++){
   m_tuple[i]/=a;
   };
};

/*
 * Internal Scalar Multiplication Operator
 */
template <class X> void CVector3<X>::operator*=(const double a){
   for(int i =0; i<DIM; i++){
   m_tuple[i]*=a;
   };
};

/*
 * Return Scalar Multiplication Operator
 */
template <class X> CVector3<X> CVector3<X>::operator*(double a_const){

    CVector3 returnvec(0.0);

    for(int i =0; i<DIM; i++){
      returnvec[i]=a_const*m_tuple[i];
    }
    return returnvec;
};

/*
 * Return Scalar Division Operator
 */
template <class X> CVector3<X> CVector3<X>::operator/(double a_const){

    CVector3 returnvec(0.0);

    for(int i =0; i<DIM; i++){
      returnvec[i]=m_tuple[i]/a_const;
    }
    return returnvec;
};

/*
 * Return Vector Subtraction Operator
 */
template <class X> CVector3<X> CVector3<X>::operator-(CVector3<X>& a_vec){

    CVector3 returnvec(0.0);

    for(int i =0; i<DIM; i++){
      returnvec[i]=m_tuple[i]-a_vec[i];
    }
    return returnvec;
};

/*
 * Return Vector Addition Operator
 */
template <class X> CVector3<X> CVector3<X>::operator+(CVector3<X>& a_vec){

    CVector3 returnvec(0.0);

    for(int i =0; i<DIM; i++){
      returnvec[i]=m_tuple[i]+a_vec[i];
    }
    return returnvec;
};

/*
 * Compute Minimum Image
 */
template <class X> void CVector3<X>::minImage(CVector3<X>& a_halfbox, CVector3<X>& a_boxsize){

    for(int i =0; i<DIM; i++){
        if ( m_tuple[i] > a_halfbox[i]){
             m_tuple[i]-=a_boxsize[i];
        }
        else if (m_tuple[i] < -a_halfbox[i]){
             m_tuple[i]+=a_boxsize[i];
        }
    }
};

/*
 * Calculate Magnitude of Vector
 */
template <class X> double CVector3<X>::SqLength() const{

   double sq=0;
   for(int i =0; i<DIM; i++){
      sq+=m_tuple[i]*m_tuple[i];
   }
   return sq;
};
