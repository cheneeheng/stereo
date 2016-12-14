// Licensed software... any distribution or replication
// without a written permission from D. Burschka burschka@ieee.org
// is prohibited

#include <math.h>
#include "XVMatrix.h"


class Vgps
{
   protected:
     inline double 	determinant(XVMatrix &x);
   public:
     			Vgps(XVMatrix &R,XVColVector &T);
			~Vgps();
      bool		validate_set(XVMatrix &new_set,
      				     XVMatrix &ref_set,
				     const XVMatrix &R,
				     const XVColVector &T);
      bool 		find_pose(XVMatrix &new_obs, XVMatrix &ref_obs,
				  XVColVector &D, XVMatrix &R,
				  XVColVector &T,
				  int num_iterations=1000);
      XVColVector	calc_3D(XVColVector new_v,XVColVector ref_v,
      				const XVMatrix &R,const XVColVector &T);
};

