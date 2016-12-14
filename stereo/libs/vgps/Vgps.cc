// Licensed software... any distribution or replication
// without a written permission from D. Burschka burschka@ieee.org
// is prohibited
#include <unistd.h>
#include "Vgps.h"

#define Sqr(a) ((a)*(a))

using namespace std;

double 
inline Vgps::determinant(XVMatrix &x)
{

 return x[0][0]*x[1][1]*x[2][2]+
      	x[0][1]*x[1][2]*x[2][0]+
      	x[1][0]*x[2][1]*x[0][2]-
        x[0][2]*x[1][1]*x[2][0]-
	x[1][0]*x[0][1]*x[2][2]-
        x[1][2]*x[2][1]*x[0][0];
}

#define AddVector(K,a,b,v) \
	{XVColVector t=(v).Column(0);\
	K[a][b]=t[0],K[a][b+1]=t[1],K[a][b+2]=t[2];}
#define ZeroVector(K,a,b) \
	{K[a][b]=0.0,K[a][b+1]=0.0,K[a][b+2]=0.0;}
#define SetResidual(res,f0,f1,f2) {\
    AddVector(res,0,0,R*ref_set.Column(f0));\
    AddVector(res,1,0,-R*ref_set.Column(f1));\
    ZeroVector(res,2,0);\
    AddVector(res,3,0,point_set_new.Column(f0));\
    AddVector(res,4,0,-point_set_new.Column(f1));\
    ZeroVector(res,5,0);\
    AddVector(res,0,3,R*ref_set.Column(f0));\
    ZeroVector(res,1,3);\
    AddVector(res,2,3,-R*ref_set.Column(f2));\
    AddVector(res,3,3,point_set_new.Column(f0));\
    ZeroVector(res,4,3);\
    AddVector(res,5,3,-point_set_new.Column(f2));}

bool
Vgps::validate_set(XVMatrix &point_set_new,XVMatrix &ref_set,
                   const XVMatrix &R,const XVColVector &T)
{
  XVMatrix residual(6,6);
  XVMatrix V(6,6);
  XVColVector eig(6);
  int v0,v1,v2; // reference vectors
  int trial=0,rank,size_result;
  bool correct;

  if(T.ip()<1e-3) return false; // translation is necessary
  // both sets nee to be equal
  if(ref_set.n_of_cols()!=point_set_new.n_of_cols()) return false;
  v0=0,v1=1,v2=2;
  do{
    SetResidual(residual,v0,v1,v2);
    residual.SVDcmp(eig,V);
    rank=0;
    for(int i=0;i<6;i++) if(fabs(eig[i])>1e-3)  rank++;
    correct=(rank==5);
    if(!correct) {
      trial++;
      v0=trial;
      v1=(trial+1)%ref_set.n_of_cols();
      v2=(trial+2)%ref_set.n_of_cols();
    }
  }while(!correct && trial<ref_set.n_of_cols());
  // could not find basic triangle?
  if (!correct) return correct;

  size_result=0;
  for(int i=0;i<ref_set.n_of_cols();i++)
  {
    if(i==v0 || i==v1 || i==v2)
    {
        size_result++;
        continue;
    }
    else
    {
      SetResidual(residual,v0,v1,i);
      residual.SVDcmp(eig,V);
      rank=0;
      for(int j=0;j<6;j++) if(fabs(eig[j])>1e-3)  rank++;
      if(rank==5)
        size_result++;
      else
      {
	//mark as wrong
        ref_set[0][i]=1000;
      }
    }
  }
  XVMatrix new_ref(3,size_result),new_set(3,size_result);
  for(int i=0,w_index=0;i<ref_set.n_of_cols();i++)
  {
   // regular normalized vectors are smaller than 1 in length
   if(ref_set[0][i]<10)
   {
    new_set[0][w_index]=point_set_new[0][i],
    new_set[1][w_index]=point_set_new[1][i],
    new_set[2][w_index]=point_set_new[2][i];
    new_ref[0][w_index]=ref_set[0][i],
    new_ref[1][w_index]=ref_set[1][i],
    new_ref[2][w_index]=ref_set[2][i];
    w_index++;
   }
  }
  ref_set=new_ref;
  point_set_new=new_set;
  return correct;
}

bool 
Vgps::find_pose(XVMatrix &new_obs,XVMatrix &ref_obs,
                XVColVector &D, XVMatrix &R,XVColVector &T,
		int num_iterations)
{
   XVColVector  z(3),p1(3),p2(3),f(3),D_n=D;
   XVMatrix pp,V(3,3),old_T(3,3);
 
   // enough correspondences for computation?
   if(new_obs.n_of_cols()<3 ||
      new_obs.n_of_cols()!=ref_obs.n_of_cols()) return false;

   // zero mean vectors
   for(int i=0;i<3;i++) z[i]=0;
   old_T=z;
   p1=z;
   for (int i=0;i<ref_obs.n_of_cols();i++) 
   {
   	p1=p1+(ref_obs.Column(i)*D[i]);
	f=R*(ref_obs.Column(i)*D[i])+T;
        D_n[i]=sqrt(f.ip());
   }
   p1/=ref_obs.n_of_cols();
   for(int lo=0;lo<num_iterations;lo++)
   //do
   {
     // calculate means
     p2=z;
     for (int i=0;i<new_obs.n_of_cols();i++)
        p2=p2+(new_obs.Column(i)*D_n[i]);
     p2/=new_obs.n_of_cols();

     for(int i=0;i<3;i++) for(int j=0;j<3;j++) R[i][j]=0.0;
     // compute covariance matrix
     for (int k=0;k<new_obs.n_of_cols();k++)
      for(int i=0;i<3;i++)
       for(int j=0;j<3;j++)
         R[i][j]+=(ref_obs[i][k]*D[k]-p1[i])*(new_obs[j][k]*D_n[k]-p2[j]);
     R.SVDcmp(f,V);
     pp=V*R.t();
     if(determinant(pp)<0)
     {
       pp[0][2]=-pp[0][2];pp[1][2]=-pp[1][2];pp[2][2]=-pp[2][2];
     }
     R=pp;
     T=p2-R*p1;
     for(int i=0;i<new_obs.n_of_cols();i++)
     {
       f=R*(ref_obs.Column(i)*D[i])+T;
       D_n[i]=sqrt(Sqr(f[0])+Sqr(f[1])+Sqr(f[2]));
     }
     f=old_T-T;
     old_T=T;
   }//while(Sqr(f[0])+Sqr(f[1])+Sqr(f[2])<0.0001);
   return true;

}

XVColVector
Vgps::calc_3D(XVColVector new_v,XVColVector ref_v,
              const XVMatrix &R,const XVColVector &T)
{
   XVMatrix disp_matrix(2,3);
   XVColVector D(2);

   AddVector(disp_matrix,0,0,new_v);
   AddVector(disp_matrix,1,0,-R*ref_v);
   D=(disp_matrix*disp_matrix.t()).i()*disp_matrix*T;
   return D;
}

Vgps::Vgps(XVMatrix &R, XVColVector &T)
{
  // initialize with no motion at all
  R[0][0]=R[1][1]=R[2][2]=1.0;
  R[0][1]=R[1][0]=R[0][2]=R[2][0]=R[1][2]=R[2][1]=0.0;
  T[0]=T[1]=T[2]=0.0;
}

Vgps::~Vgps()
{
}
