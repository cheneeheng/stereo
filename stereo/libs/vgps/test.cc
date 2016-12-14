#include <stdio.h>
#include "Vgps.h"
using namespace std;

#define Sqr(a) ((a)*(a))
static
void init_vec(XVMatrix &vec_1,XVMatrix &vec_2,XVMatrix &vec_3,
              XVColVector &D)
{
  FILE *fptr=fopen("robot_input.txt","r");
  for(int i=0;i<1655;i++)
  {
    float n1;
    fscanf(fptr,"%lf %lf %lf %lf %lf %lf %lf",
        &(vec_1[0][i]),&(vec_1[1][i]),&(vec_2[0][i]),&(vec_2[1][i]),
        &(vec_3[0][i]),&(vec_3[1][i]),&(vec_3[2][i]));
    vec_1[2][i]=1.0,vec_2[2][i]=1.0;
    vec_3[2][i]+=5;
    vec_3[1][i]+=10;
    vec_1[0][i]=(vec_1[0][i]-512)/960;
    vec_1[1][i]=(vec_1[1][i]-384)/960;
    vec_2[0][i]=(vec_2[0][i]-512)/960;
    vec_2[1][i]=(vec_2[1][i]-384)/960;
    n1=sqrt(Sqr(vec_1[0][i])+Sqr(vec_1[1][i])+Sqr(vec_1[2][i]));
    vec_1[0][i]/=n1;
    vec_1[1][i]/=n1;
    vec_1[2][i]/=n1;
    n1=sqrt(Sqr(vec_2[0][i])+Sqr(vec_2[1][i])+Sqr(vec_2[2][i]));
    vec_2[0][i]/=n1;
    vec_2[1][i]/=n1;
    vec_2[2][i]/=n1;
    D[i]=sqrt(Sqr(vec_3[0][i])+Sqr(vec_3[1][i])+Sqr(vec_3[2][i]));
  }
  fclose(fptr);
}
// R= 0.1969 0.1710 -0.9654 14.3
//   -0.9448 0.2962 -0.1401  2.2
//    0.2620 0.9397  0.2198 -3.4

//#define NEVER
int main(int argc, char ** argv)
{
  XVColVector D(1655),T(3),T1(3);
  static XVMatrix vec_1(3,1655),vec_2(3,1655),vec_3(3,1655),R(3,3),R1(3,3);
  Vgps vgps(R,T);
 R1=R; T1=T;
 init_vec(vec_1,vec_2,vec_3,D);
 vgps.find_pose(vec_1,vec_3,D,R,T,10000);
 vgps.find_pose(vec_2,vec_3,D,R1,T1,10000);
 R=R1*R.t();
 T=T1-R*T;
 //cout << "D_n=" << endl << D << endl;
 cout << "T=" << endl<< T << endl;
 cout << "R=" << endl << R << endl;
 cout << "T1=" << endl<< T1 << endl;
 cout << "R1=" << endl << R1 << endl;
// cout << "D_n=" << endl << D << endl;
  return 0;
}
