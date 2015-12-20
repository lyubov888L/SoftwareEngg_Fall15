
#include "include/molecule.h"
int main(){

//    CVector3<double> position(3.4);
//    position.Display();
//    position[1]=1;
//    position.Display();
//    CAtom test(position, 1, 0);
//    vector<double> a_tuple(6,2.2);
//    a_tuple[1]=1.2; a_tuple[5]=2.3;
//    vector<int> a_atomindex(2);
//    a_atomindex[0]=0;
//    a_atomindex[1]=1;
//    vector<int> a_atomtype(2);
//    a_atomtype[0]=0;
//    a_atomtype[1]=1;
//    CMolecule newmol(2, a_tuple, a_atomindex, a_atomtype);

   CVector3<double> position1(3.4);
   CVector3<double>& position2;
   CVector3<double> position3;

   position1-=position2;
   position1.Display();
   position2.Display();
   position1.copy(position2);
   position1.Display();

   //position2=position1*1.1;
   position2.Display();
   position3=position2-position1;
   position3.Display();

    return 1;
   
}
