#include "Vgps.h"
using namespace std;

static
void init_vec(XVMatrix &vec_1,XVMatrix &vec_2,XVColVector &D)
{
// reference points
  vec_1[0][0]=0.8165;
  vec_1[1][0]=0.4082;
  vec_1[2][0]=0.4082;
  
  vec_1[0][1]=0.9478;
  vec_1[1][1]=0.2962;
  vec_1[2][1]=0.1185;
  
  vec_1[0][2]=0.5185;
  vec_1[1][2]=0.8296;
  vec_1[2][2]=0.2074;
  
  vec_1[0][3]=0.5345;
  vec_1[1][3]=-0.2673;
  vec_1[2][3]=0.8018;
  D[0]=2.4495;D[1]=8.4410;D[2]=4.8218;
  D[3]=3.7417;

// observation
  vec_2[0][0]=0.9919;
  vec_2[1][0]=0.0333;
  vec_2[2][0]=-0.1225;
  
  vec_2[0][1]=0.9521;
  vec_2[1][1]=-0.2954;
  vec_2[2][1]=0.0785;
  
  vec_2[0][2]=0.9946;
  vec_2[1][2]=0.0605;
  vec_2[2][2]=0.0846;
  
  vec_2[0][3]=0.9645;
  vec_2[1][3]=-0.0337;
  vec_2[2][3]=-0.2618;
}

// R= 0.1969 0.1710 -0.9654 14.3
//   -0.9448 0.2962 -0.1401  2.2
//    0.2620 0.9397  0.2198 -3.4

#define NEVER
int main(int argc, char ** argv)
{
  XVColVector D(4),T(3);
  XVMatrix vec_1(3,4),vec_2(3,4),R(3,3);
  Vgps vgps(R,T);

 init_vec(vec_1,vec_2,D);
 for(int k=0;k<100;k++)
 {
  vgps.find_pose(vec_2,vec_1,D,R,T,400);
#ifdef NEVER
  cout << "T=" << endl<< T << endl;
  cout << "R=" << endl << R << endl;
  cout << "D_n=" << endl << D << endl;
#endif
  vgps.calc_3D(vec_2.Column(1),vec_1.Column(1),R,T);
#ifdef NEVER
  cerr << vgps.validate_set(vec_2,vec_1,R,T)<< endl;
  cerr << vec_2.n_of_cols() << endl;
#endif
}
  return 0;
}
