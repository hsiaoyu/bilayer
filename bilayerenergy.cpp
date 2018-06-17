#include "DDd01.h"
#include "DDd02.h"
#include "DDd03.h"
#include "DDd11.h"
#include "DDd12.h"
#include "DDd13.h"
#include "DDd21.h"
#include "DDd22.h"
#include "DDd23.h"
#include "DDdtop01.h"
#include "DDdtop02.h"
#include "DDdtop03.h"
#include "DDdtop11.h"
#include "DDdtop12.h"
#include "DDdtop13.h"
#include "DDdtop21.h"
#include "DDdtop22.h"
#include "DDdtop23.h"
#include "DDv01.h"
#include "DDv02.h"
#include "DDv03.h"
#include "DDv11.h"
#include "DDv12.h"
#include "DDv13.h"
#include "DDv21.h"
#include "DDv22.h"
#include "DDv23.h"
void hessiantri(double dAtop, double dA, double L1top, double L2top, double L1, double L2, double htop, double hbot, std::vector<Eigen::Vector3d> d, std::vector<Eigen::Vector3d> dtop, std::vector<Eigen::Vector3d> v, Eigen::Matrix3d ainvtop, Eigen::Matrix3d ainv, Eigen::MatrixXd & hessian)
{
   hessian.resize(27,27);
   Eigen::RowVectorXd tmp(27);

   DDv01(double dAtop, double dA, double L1top, double L2top, double L1, double L2, double htop, double hbot, std::vector<Eigen::Vector3d> d, std::vector<Eigen::Vector3d> dtop, std::vector<Eigen::Vector3d> v, Eigen::Matrix3d ainvtop, Eigen::Matrix3d ainv, Eigen::RowVectorXd & tmp);
   hessian.row(0) = tmp;
   DDv02(double dAtop, double dA, double L1top, double L2top, double L1, double L2, double htop, double hbot, std::vector<Eigen::Vector3d> d, std::vector<Eigen::Vector3d> dtop, std::vector<Eigen::Vector3d> v, Eigen::Matrix3d ainvtop, Eigen::Matrix3d ainv, Eigen::RowVectorXd & tmp);
   hessian.row(1) = tmp;
   DDv03(double dAtop, double dA, double L1top, double L2top, double L1, double L2, double htop, double hbot, std::vector<Eigen::Vector3d> d, std::vector<Eigen::Vector3d> dtop, std::vector<Eigen::Vector3d> v, Eigen::Matrix3d ainvtop, Eigen::Matrix3d ainv, Eigen::RowVectorXd & tmp);
   hessian.row(2) = tmp;
   DDd01(double dAtop, double dA, double L1top, double L2top, double L1, double L2, double htop, double hbot, std::vector<Eigen::Vector3d> d, std::vector<Eigen::Vector3d> dtop, std::vector<Eigen::Vector3d> v, Eigen::Matrix3d ainvtop, Eigen::Matrix3d ainv, Eigen::RowVectorXd & tmp);
   hessian.row(3) = tmp;
   DDd02(double dAtop, double dA, double L1top, double L2top, double L1, double L2, double htop, double hbot, std::vector<Eigen::Vector3d> d, std::vector<Eigen::Vector3d> dtop, std::vector<Eigen::Vector3d> v, Eigen::Matrix3d ainvtop, Eigen::Matrix3d ainv, Eigen::RowVectorXd & tmp);
   hessian.row(4) = tmp;
   DDd03(double dAtop, double dA, double L1top, double L2top, double L1, double L2, double htop, double hbot, std::vector<Eigen::Vector3d> d, std::vector<Eigen::Vector3d> dtop, std::vector<Eigen::Vector3d> v, Eigen::Matrix3d ainvtop, Eigen::Matrix3d ainv, Eigen::RowVectorXd & tmp);
   hessian.row(5) = tmp;
   DDdtop01(double dAtop, double dA, double L1top, double L2top, double L1, double L2, double htop, double hbot, std::vector<Eigen::Vector3d> d, std::vector<Eigen::Vector3d> dtop, std::vector<Eigen::Vector3d> v, Eigen::Matrix3d ainvtop, Eigen::Matrix3d ainv, Eigen::RowVectorXd & tmp);
   hessian.row(6) = tmp;
   DDdtop02(double dAtop, double dA, double L1top, double L2top, double L1, double L2, double htop, double hbot, std::vector<Eigen::Vector3d> d, std::vector<Eigen::Vector3d> dtop, std::vector<Eigen::Vector3d> v, Eigen::Matrix3d ainvtop, Eigen::Matrix3d ainv, Eigen::RowVectorXd & tmp);
   hessian.row(7) = tmp;
   DDdtop03(double dAtop, double dA, double L1top, double L2top, double L1, double L2, double htop, double hbot, std::vector<Eigen::Vector3d> d, std::vector<Eigen::Vector3d> dtop, std::vector<Eigen::Vector3d> v, Eigen::Matrix3d ainvtop, Eigen::Matrix3d ainv, Eigen::RowVectorXd & tmp);
   hessian.row(8) = tmp;

   DDv11(double dAtop, double dA, double L1top, double L2top, double L1, double L2, double htop, double hbot, std::vector<Eigen::Vector3d> d, std::vector<Eigen::Vector3d> dtop, std::vector<Eigen::Vector3d> v, Eigen::Matrix3d ainvtop, Eigen::Matrix3d ainv, Eigen::RowVectorXd & tmp);
   hessian.row(9) = tmp;
   DDv12(double dAtop, double dA, double L1top, double L2top, double L1, double L2, double htop, double hbot, std::vector<Eigen::Vector3d> d, std::vector<Eigen::Vector3d> dtop, std::vector<Eigen::Vector3d> v, Eigen::Matrix3d ainvtop, Eigen::Matrix3d ainv, Eigen::RowVectorXd & tmp);
   hessian.row(10) = tmp;
   DDv13(double dAtop, double dA, double L1top, double L2top, double L1, double L2, double htop, double hbot, std::vector<Eigen::Vector3d> d, std::vector<Eigen::Vector3d> dtop, std::vector<Eigen::Vector3d> v, Eigen::Matrix3d ainvtop, Eigen::Matrix3d ainv, Eigen::RowVectorXd & tmp);
   hessian.row(11) = tmp;
   DDd11(double dAtop, double dA, double L1top, double L2top, double L1, double L2, double htop, double hbot, std::vector<Eigen::Vector3d> d, std::vector<Eigen::Vector3d> dtop, std::vector<Eigen::Vector3d> v, Eigen::Matrix3d ainvtop, Eigen::Matrix3d ainv, Eigen::RowVectorXd & tmp);
   hessian.row(12) = tmp;
   DDd12(double dAtop, double dA, double L1top, double L2top, double L1, double L2, double htop, double hbot, std::vector<Eigen::Vector3d> d, std::vector<Eigen::Vector3d> dtop, std::vector<Eigen::Vector3d> v, Eigen::Matrix3d ainvtop, Eigen::Matrix3d ainv, Eigen::RowVectorXd & tmp);
   hessian.row(13) = tmp;
   DDd13(double dAtop, double dA, double L1top, double L2top, double L1, double L2, double htop, double hbot, std::vector<Eigen::Vector3d> d, std::vector<Eigen::Vector3d> dtop, std::vector<Eigen::Vector3d> v, Eigen::Matrix3d ainvtop, Eigen::Matrix3d ainv, Eigen::RowVectorXd & tmp);
   hessian.row(14) = tmp;
   DDdtop11(double dAtop, double dA, double L1top, double L2top, double L1, double L2, double htop, double hbot, std::vector<Eigen::Vector3d> d, std::vector<Eigen::Vector3d> dtop, std::vector<Eigen::Vector3d> v, Eigen::Matrix3d ainvtop, Eigen::Matrix3d ainv, Eigen::RowVectorXd & tmp);
   hessian.row(15) = tmp;
   DDdtop12(double dAtop, double dA, double L1top, double L2top, double L1, double L2, double htop, double hbot, std::vector<Eigen::Vector3d> d, std::vector<Eigen::Vector3d> dtop, std::vector<Eigen::Vector3d> v, Eigen::Matrix3d ainvtop, Eigen::Matrix3d ainv, Eigen::RowVectorXd & tmp);
   hessian.row(16) = tmp;
   DDdtop13(double dAtop, double dA, double L1top, double L2top, double L1, double L2, double htop, double hbot, std::vector<Eigen::Vector3d> d, std::vector<Eigen::Vector3d> dtop, std::vector<Eigen::Vector3d> v, Eigen::Matrix3d ainvtop, Eigen::Matrix3d ainv, Eigen::RowVectorXd & tmp);
   hessian.row(17) = tmp;

   DDv21(double dAtop, double dA, double L1top, double L2top, double L1, double L2, double htop, double hbot, std::vector<Eigen::Vector3d> d, std::vector<Eigen::Vector3d> dtop, std::vector<Eigen::Vector3d> v, Eigen::Matrix3d ainvtop, Eigen::Matrix3d ainv, Eigen::RowVectorXd & tmp);
   hessian.row(18) = tmp;
   DDv22(double dAtop, double dA, double L1top, double L2top, double L1, double L2, double htop, double hbot, std::vector<Eigen::Vector3d> d, std::vector<Eigen::Vector3d> dtop, std::vector<Eigen::Vector3d> v, Eigen::Matrix3d ainvtop, Eigen::Matrix3d ainv, Eigen::RowVectorXd & tmp);
   hessian.row(19) = tmp;
   DDv23(double dAtop, double dA, double L1top, double L2top, double L1, double L2, double htop, double hbot, std::vector<Eigen::Vector3d> d, std::vector<Eigen::Vector3d> dtop, std::vector<Eigen::Vector3d> v, Eigen::Matrix3d ainvtop, Eigen::Matrix3d ainv, Eigen::RowVectorXd & tmp);
   hessian.row(20) = tmp;
   DDd21(double dAtop, double dA, double L1top, double L2top, double L1, double L2, double htop, double hbot, std::vector<Eigen::Vector3d> d, std::vector<Eigen::Vector3d> dtop, std::vector<Eigen::Vector3d> v, Eigen::Matrix3d ainvtop, Eigen::Matrix3d ainv, Eigen::RowVectorXd & tmp);
   hessian.row(21) = tmp;
   DDd22(double dAtop, double dA, double L1top, double L2top, double L1, double L2, double htop, double hbot, std::vector<Eigen::Vector3d> d, std::vector<Eigen::Vector3d> dtop, std::vector<Eigen::Vector3d> v, Eigen::Matrix3d ainvtop, Eigen::Matrix3d ainv, Eigen::RowVectorXd & tmp);
   hessian.row(22) = tmp;
   DDd23(double dAtop, double dA, double L1top, double L2top, double L1, double L2, double htop, double hbot, std::vector<Eigen::Vector3d> d, std::vector<Eigen::Vector3d> dtop, std::vector<Eigen::Vector3d> v, Eigen::Matrix3d ainvtop, Eigen::Matrix3d ainv, Eigen::RowVectorXd & tmp);
   hessian.row(23) = tmp;
   DDdtop21(double dAtop, double dA, double L1top, double L2top, double L1, double L2, double htop, double hbot, std::vector<Eigen::Vector3d> d, std::vector<Eigen::Vector3d> dtop, std::vector<Eigen::Vector3d> v, Eigen::Matrix3d ainvtop, Eigen::Matrix3d ainv, Eigen::RowVectorXd & tmp);
   hessian.row(24) = tmp;
   DDdtop22(double dAtop, double dA, double L1top, double L2top, double L1, double L2, double htop, double hbot, std::vector<Eigen::Vector3d> d, std::vector<Eigen::Vector3d> dtop, std::vector<Eigen::Vector3d> v, Eigen::Matrix3d ainvtop, Eigen::Matrix3d ainv, Eigen::RowVectorXd & tmp);
   hessian.row(25) = tmp;
   DDdtop23(double dAtop, double dA, double L1top, double L2top, double L1, double L2, double htop, double hbot, std::vector<Eigen::Vector3d> d, std::vector<Eigen::Vector3d> dtop, std::vector<Eigen::Vector3d> v, Eigen::Matrix3d ainvtop, Eigen::Matrix3d ainv, Eigen::RowVectorXd & tmp);
   hessian.row(26) = tmp;
}

double energytot(const RestConfig & simRest, const Eigen::VectorXd & dofs, Eigen::VectorXd & dE. std::vector<Eigen::Triplet<double> > & hesslist){
   double e = 0;
   Eigen::VectorXd g(27);
   dE.setZero();
   hesslist.clear();
   Eigen::MatrixXd hess(27,27);
   for(int i = 0; i < simRest.F.rows(); i++){
       Eigen::Vector3d tmp;
       std::vector<Eigen::Vector3d> d;
       std::vector<Eigen::Vector3d> dtop;
       std::vector<Eigen::Vector3d> v;
       int FaceID[3];
       for(int j = 0; j < 3; j++){
          FaceID[j] = simRest.F(i,j);
          tmp = dofs.segment(9*FaceID[j], 3);
          v.push_back(tmp);
       //   std::cout << "v" << FaceID[j] << "  " << v[j](0) << " " << v[j](1) << " " << v[j](2) << std::endl; 
          tmp = dofs.segment(9*FaceID[j]+3, 3);
          d.push_back(tmp);
          tmp = dofs.segment(9*FaceID[j]+6, 3);
          dtop.push_back(tmp);
       }
       e += energytri(simRest.dAs[i], simRest.dAs[i],  simRest.L1top, simRest.L2top, simRest.L1, simRest.L2, simRest.htop, simRest.hbot, d, dtop, v, simRest.ainvs[i], simRest.ainvs[i]);
       gradtri(simRest.dAs[i], simRest.dAs[i], simRest.L1top, simRest.L2top, simRest.L1, simRest.L2, simRest.htop, simRest.hbot, d, dtop, v, simRest.ainvs[i], simRest.ainvs[i], g);
       hessiantri(simRest.dAs[i], simRest.dAs[i], simRest.L1top, simRest.L2top, simRest.L1, simRest.L2, simRest.htop, simRest.hbot, d, dtop, v, simRest.ainvs[i], simRest.ainvs[i], hess);
       for(int j = 0; j < 3; j++){
           for(int k = 0; k < 9; k++){
              dE(9*FaceID[j]+k) += (double) (simRest.fixDOF(9*FaceID[j]+k)) * g(9*j+k);
           }
       }
       for(int j = 0; j < 27; j++){
           for(int k = 0; k < 27; k++){
              hesslist.push_back(Eigen::Triplet<double>(9*FaceID[j/9]+(j%9), 9*FaceID[k/9]+(k%9)), (double)(simRest.fixDOF(9*FaceID[j/9]+(j%9)) * simRest.fixDOF(9*FaceID[k/9]+(k%9)))* hess(j,k));
           }
       }
   }
   
   return e;
}
