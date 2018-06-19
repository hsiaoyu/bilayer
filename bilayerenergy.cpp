#include "bilayerenergy.h"
//#include <iostream>
#include <Eigen/Sparse>
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
#include "energy1.h"
#include "energy2.h"
#include "energy3.h"
#include "energy4.h"
#include "energy5.h"
#include "Dd01.h"
#include "Dd02.h"
#include "Dd03.h"
#include "Dd11.h"
#include "Dd12.h"
#include "Dd13.h"
#include "Dd21.h"
#include "Dd22.h"
#include "Dd23.h"
#include "Ddtop01.h"
#include "Ddtop02.h"
#include "Ddtop03.h"
#include "Ddtop11.h"
#include "Ddtop12.h"
#include "Ddtop13.h"
#include "Ddtop21.h"
#include "Ddtop22.h"
#include "Ddtop23.h"
#include "Dv01.h"
#include "Dv02.h"
#include "Dv03.h"
#include "Dv11.h"
#include "Dv12.h"
#include "Dv13.h"
#include "Dv21.h"
#include "Dv22.h"
#include "Dv23.h"

double energytot(const RestConfig & simRest, const Eigen::VectorXd & dofs, Eigen::VectorXd * dE, std::vector<Eigen::Triplet<double> > * hesslist){
   double e = 0;
   Eigen::VectorXd g(27);
   Eigen::MatrixXd hess(27,27);
   if(dE)
      dE->setZero();
   if(hesslist)
      hesslist->clear();
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
          tmp = dofs.segment(9*FaceID[j]+3, 3);
          d.push_back(tmp);
          tmp = dofs.segment(9*FaceID[j]+6, 3);
          dtop.push_back(tmp);
       }
       e += energytri(simRest.dAs[i], simRest.dAs[i],  simRest.L1top, simRest.L2top, simRest.L1, simRest.L2, simRest.htop, simRest.hbot, d, dtop, v, simRest.ainvs[i], simRest.ainvs[i]);
       
       if(dE){
          gradtri(simRest.dAs[i], simRest.dAs[i], simRest.L1top, simRest.L2top, simRest.L1, simRest.L2, simRest.htop, simRest.hbot, d, dtop, v, simRest.ainvs[i], simRest.ainvs[i], g);
          for(int j = 0; j < 3; j++){
              for(int k = 0; k < 9; k++){
                 dE->coeffRef(9*FaceID[j]+k) += (double) (simRest.fixDOF(9*FaceID[j]+k)) * g(9*j+k);
              }
          }
       }
 
       if(hesslist){
          hessiantri(simRest.dAs[i], simRest.dAs[i], simRest.L1top, simRest.L2top, simRest.L1, simRest.L2, simRest.htop, simRest.hbot, d, dtop, v, simRest.ainvs[i], simRest.ainvs[i], hess);
          for(int j = 0; j < 27; j++){
              for(int k = 0; k < 27; k++){
                 hesslist->push_back(Eigen::Triplet<double>(9*FaceID[j/9]+(j%9), 9*FaceID[k/9]+(k%9), (double)(simRest.fixDOF(9*FaceID[j/9]+(j%9)) * simRest.fixDOF(9*FaceID[k/9]+(k%9)))* hess(j,k)));
              }
          }
       }
   }
   
   return e;
}

void hessiantri(double dAtop, double dA, double L1top, double L2top, double L1, double L2, double htop, double hbot, std::vector<Eigen::Vector3d> d, std::vector<Eigen::Vector3d> dtop, std::vector<Eigen::Vector3d> v, Eigen::Matrix3d ainvtop, Eigen::Matrix3d ainv, Eigen::MatrixXd & hessian)
{
   hessian.resize(27,27);
   Eigen::RowVectorXd tmp(27);

   DDv01(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv, tmp);
   hessian.row(0) = tmp;
   DDv02(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv, tmp);
   hessian.row(1) = tmp;
   DDv03(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv, tmp);
   hessian.row(2) = tmp;
   DDd01(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv, tmp);
   hessian.row(3) = tmp;
   DDd02(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv, tmp);
   hessian.row(4) = tmp;
   DDd03(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv, tmp);
   hessian.row(5) = tmp;
   DDdtop01(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv, tmp);
   hessian.row(6) = tmp;
   DDdtop02(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv, tmp);
   hessian.row(7) = tmp;
   DDdtop03(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv, tmp);
   hessian.row(8) = tmp;

   DDv11(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv, tmp);
   hessian.row(9) = tmp;
   DDv12(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv, tmp);
   hessian.row(10) = tmp;
   DDv13(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv, tmp);
   hessian.row(11) = tmp;
   DDd11(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv, tmp);
   hessian.row(12) = tmp;
   DDd12(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv, tmp);
   hessian.row(13) = tmp;
   DDd13(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv, tmp);
   hessian.row(14) = tmp;
   DDdtop11(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv, tmp);
   hessian.row(15) = tmp;
   DDdtop12(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv, tmp);
   hessian.row(16) = tmp;
   DDdtop13(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv, tmp);
   hessian.row(17) = tmp;

   DDv21(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv, tmp);
   hessian.row(18) = tmp;
   DDv22(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv, tmp);
   hessian.row(19) = tmp;
   DDv23(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv, tmp);
   hessian.row(20) = tmp;
   DDd21(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv, tmp);
   hessian.row(21) = tmp;
   DDd22(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv, tmp);
   hessian.row(22) = tmp;
   DDd23(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv, tmp);
   hessian.row(23) = tmp;
   DDdtop21(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv, tmp);
   hessian.row(24) = tmp;
   DDdtop22(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv, tmp);
   hessian.row(25) = tmp;
   DDdtop23(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv, tmp);
   hessian.row(26) = tmp;

   //std::cout << "hessian" << std::endl << hessian << std::endl;
}

double energytri(double dAtop, double dA, double L1top, double L2top, double L1, double L2, double htop, double hbot, std::vector<Eigen::Vector3d> d, std::vector<Eigen::Vector3d> dtop, std::vector<Eigen::Vector3d> v, Eigen::Matrix3d ainvtop, Eigen::Matrix3d ainv)
{
  double e, tmp;
 
  e = energy1(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv);
  e += energy2(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv);
  e += energy3(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv);
  e += energy4(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv);
  e += energy5(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv);
 
  return e;
}

void gradtri(double dAtop, double dA, double L1top, double L2top, double L1, double L2, double htop, double hbot, std::vector<Eigen::Vector3d> d, std::vector<Eigen::Vector3d> dtop, std::vector<Eigen::Vector3d> v, Eigen::Matrix3d ainvtop, Eigen::Matrix3d ainv, Eigen::VectorXd & gradient)
{
      gradient[0] = Dv01(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv);
      gradient[1] = Dv02(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv);
      gradient[2] = Dv03(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv);
      gradient[3] = Dd01(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv);
      gradient[4] = Dd02(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv);
      gradient[5] = Dd03(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv);
      gradient[6] = Ddtop01(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv);
      gradient[7] = Ddtop02(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv);
      gradient[8] = Ddtop03(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv);
      
      gradient[9] = Dv11(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv);
      gradient[10] = Dv12(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv);
      gradient[11] = Dv13(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv);
      gradient[12] = Dd11(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv);
      gradient[13] = Dd12(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv);
      gradient[14] = Dd13(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv);
      gradient[15] = Ddtop11(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv);
      gradient[16] = Ddtop12(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv);
      gradient[17] = Ddtop13(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv);
      
      gradient[18] = Dv21(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv);
      gradient[19] = Dv22(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv);
      gradient[20] = Dv23(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv);
      gradient[21] = Dd21(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv);
      gradient[22] = Dd22(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv);
      gradient[23] = Dd23(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv);
      gradient[24] = Ddtop21(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv);
      gradient[25] = Ddtop22(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv);
      gradient[26] = Ddtop23(dAtop, dA, L1top, L2top, L1, L2, htop, hbot, d, dtop, v, ainvtop, ainv);
}

