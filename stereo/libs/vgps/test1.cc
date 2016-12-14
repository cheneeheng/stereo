#include <stdio.h>
#include "Vgps.h"
using namespace std;

#define Sqr(a) ((a)*(a))
static
void init_vec(XVMatrix &vec_1,XVMatrix &vec_2, XVColVector &D)
{
  FILE *fptr=fopen("test.txt","r");
  for(int i=0;i<30;i++)
  {
    fscanf(fptr,"%lf %lf %lf %lf %lf %lf %lf",
        &(vec_1[0][i]),&(vec_1[1][i]),&(vec_1[2][i]),&(vec_2[0][i]),
        &(vec_2[1][i]),&(vec_2[2][i]),&(D[i]));
  }
  fclose(fptr);
}
// R= 0.1969 0.1710 -0.9654 14.3
//   -0.9448 0.2962 -0.1401  2.2
//    0.2620 0.9397  0.2198 -3.4

//#define NEVER
int main(int argc, char ** argv)
{
  XVColVector D(30),T(3),T1(3);
  static XVMatrix vec_1(3,30),vec_2(3,30),vec_3(3,30),R(3,3),R1(3,3);
  Vgps vgps(R,T);
 R1=R; T1=T;
 init_vec(vec_1,vec_2,D);
 vgps.find_pose(vec_1,vec_2,D,R,T,800);
 //cout << "D_n=" << endl << D << endl;
 cout << "T=" << endl<< T << endl;
 cout << "R=" << endl << R << endl;
// cout << "D_n=" << endl << D << endl;
  return 0;
}
