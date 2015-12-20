/*!
 * Cvector3 class is an extension of the standard vector class with additional functionality to be able to hold and manipulate coordinates and velocities for the atoms.
 */

#ifndef _CVECTOR3_H_
#define _CVECTOR3_H_
#include <vector>
#include "config.h"
using namespace std;
template <class X>
class CVector3
{
public:
   /*Constructors*/
   CVector3(){}; /*!< Default constructor */
   CVector3(const X a_value); /*!< Constructor with input value */

   X& operator[](const int); /*!< Access operator to access the desired element from the internal tuple */

   void Display(); /*!< Writes out the object tuple */

   /*Operator overloads*/
   void copy(CVector3&); /*!< Copy constructor */
   void operator+=( CVector3&); /*!< Internal Vector Addition Operator */
   void operator-=( CVector3&); /*!< Internal Vector Subtraction Operator */
   void operator/=(const double); /*!< Internal Scalar Division Operator */
   void operator*=(const double); /*!< Internal Scalar Multiplication Operator */

   CVector3 operator*(double a_const); /*!< Return Scalar Multiplication Operator */
   CVector3 operator/(double a_const); /*!< Return Scalar Division Operator */
   CVector3 operator-(CVector3& a_vec); /*!< Return Vector Subtraction Operator */
   CVector3 operator+(CVector3& a_vec); /*!< Return Vector Addition Operator */

   void minImage(CVector3& a_halfbox, CVector3& a_boxsize); /*!< Compute Minimum Image */
   double SqLength() const; /*!< Calculate Magnitude of Vector */

private:
   std::vector<X> m_tuple;  /*!< Vector of coordinates/velocities/forces */

};

#include "vector3implem.h"

#endif
