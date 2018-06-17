#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <fstream>
#include "simSetup.h"
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

double energytri(double dAtop, double dA, double L1top, double L2top, double L1, double L2, double htop, double hbot, std::vector<Eigen::Vector3d> d, std::vector<Eigen::Vector3d> dtop, std::vector<Eigen::Vector3d> v, Eigen::Matrix3d ainvtop, Eigen::Matrix3d ainv);
void gradtri(double dAtop, double dA, double L1top, double L2top, double L1, double L2, double htop, double hbot, std::vector<Eigen::Vector3d> d, std::vector<Eigen::Vector3d> dtop, std::vector<Eigen::Vector3d> v, Eigen::Matrix3d ainvtop, Eigen::Matrix3d ainv, Eigen::VectorXd & gradient);
double energytot(const RestConfig & simRest, const Eigen::VectorXd & dofs, Eigen::VectorXd & dE);

int main(int argc, char *argv[]){
   if(argc < 2){
       std::cout << "Need to pass in input file!" << std::endl;
       return 1;
    }
   
    std::ifstream infile(argv[1]);
    std::string str;
    getline(infile, str);
    std::string outName;
    outName = str.substr(0,str.size()-4);
    outName += "CPD.obj";  
    const char *cstr = str.c_str();
    RestConfig simRest;
    igl::readOBJ(cstr, simRest.V, simRest.F);
    getline(infile, str);
    simRest.fixDOF.resize(9 * simRest.V.rows());
    genConstraints(str, simRest.fixDOF, simRest.V.rows());
    //for(int i = 0; i < simRest.V.rows(); i++){
    //    for(int j = 0; j < 9; j++)
    //       std::cout << simRest.fixDOF(9*i+j) << " " ;
    //    std::cout << std::endl;
    //}

    double Youngs, poisson, Youngstop, poissontop, htop, hbot, compression, tol, stepsize, maxiter, lsalpha;
    infile >> Youngs;
    infile >> poisson;
    infile >> Youngstop;
    infile >> poissontop;
    infile >> simRest.hbot;
    infile >> simRest.htop;
    infile >> compression;
    infile >> stepsize;
    infile >> tol;
    infile >> maxiter;
    infile >> lsalpha;

    for(int i=0; i<simRest.F.rows(); i++){
       Eigen::Matrix3d tmp;
       Eigen::Vector3d v0, v1, v2;
       v0 = simRest.V.row(simRest.F(i,0));
       v1 = simRest.V.row(simRest.F(i,1));
       v2 = simRest.V.row(simRest.F(i,2));
       tmp = calcAbar(v0, v1, v2);
       simRest.dAs.push_back(sqrt(tmp.determinant()));
       simRest.ainvs.push_back(tmp.inverse());
    }
    
    simRest.L1top =  0.5 * LameAlpha(Youngstop,poissontop);
    simRest.L2top = LameBeta(Youngstop,poissontop);
    simRest.L1 =  0.5 * LameAlpha(Youngs,poisson);
    simRest.L2 = LameBeta(Youngs,poisson);
    
    Eigen::VectorXd simDof(9 * simRest.V.rows());
    setupDOF(simRest.V, compression, simDof); 
    
    double residual = 100;
    double Energy;
    int count = 0;
    double beta = 0.5;
    double linestep = stepsize;
    Eigen::VectorXd dE(9 * simRest.V.rows()), guess(9 * simRest.V.rows());
    while(residual > tol){
           double Energy = energytot(simRest, simDof, dE);
           //for(int i = 0; i < simRest.V.rows(); i ++){
           //    std::cout << "total grad  " << i <<  std::endl; 
           //    for(int j = 0; j < 9; j ++)
           //         std::cout <<  dE(9*i+j) << " ";
           //    std::cout << std::endl; 
           //}
           linestep = stepsize;
           guess = simDof - linestep * dE;
           Eigen::VectorXd dEguess(9 * simRest.V.rows());
           double guessE = energytot(simRest, guess, dEguess);
           while (guessE > Energy - lsalpha * linestep * dE.norm() * dE.norm()){
           //while (guessE > Energy){
                 linestep  = beta * linestep;
                 std::cout << "current linestep = " << linestep << std::endl;
                 guess = simDof - linestep * dE;
                 guessE = energytot(simRest, guess, dEguess);
           }
           residual = dEguess.norm();
           simDof = guess;
           count ++;
           std::cout << "Energy = " << guessE << " |dE| = " << residual <<std::endl;
    }
    
    Eigen::MatrixXd Vout(simRest.V.rows() * 3, 3);
    Eigen::MatrixXi Fout(simRest.F.rows() * 3 ,3);
    Fout.block(0,0,simRest.F.rows(),3)= simRest.F;   

    for(int i = 0; i < simRest.V.rows(); i++){
        for(int j = 0; j < 3; j++){
            Vout(i,j) = simDof(9*i+j);
            Vout(i+simRest.V.rows(),j) = simDof(9*i+j)+simRest.hbot*simDof(9*i+3+j);
            Vout(i+2*simRest.V.rows(),j) = simDof(9*i+j)+simRest.hbot*simDof(9*i+3+j)+simRest.htop*simDof(9*i+6+j);
        }
    }
   
    int nVerts = simRest.V.rows();
    for(int i = 0; i < simRest.F.rows(); i++){
        Eigen::RowVector3i nv(nVerts, nVerts, nVerts);
        Fout.row(i+simRest.F.rows()) = simRest.F.row(i) + nv;
        Fout.row(i+2*simRest.F.rows()) = simRest.F.row(i) + 2 * nv;
    } 
  
    igl::writeOBJ(outName, Vout, Fout); 
 
//    std::vector<Eigen::Vector3d> d;
//    Eigen::Vector3d tmp;
//    tmp << -1,-2,2;
//    d.push_back(tmp);
//    tmp << 1,3,7;
//    d.push_back(tmp);
//    tmp << 4,5,6;
//    d.push_back(tmp);
//    
//    Eigen::Matrix3d ainv;
//    ainv.setConstant(3);
//
//    double L, Ltop, dA, h, energy;
//    L = 3;
//    Ltop = 2;
//    dA = 0.5;
//    h = 0.1;
//    energy = energytri(dA,dA,Ltop, Ltop, L, L, h, h, d, d, d, ainv, ainv);
//    std::cout << "energy = " << energy << std::endl;
//    Eigen::VectorXd grads(27);
//    gradtri(dA,dA, Ltop, Ltop, L, L, h, h, d, d, d, ainv, ainv, grads);
//    std::cout << "grad " << std::endl << grads << std::endl;
    return 0;
}

double energytot(const RestConfig & simRest, const Eigen::VectorXd & dofs, Eigen::VectorXd & dE){
   double e = 0;
   Eigen::VectorXd g(27);
   dE.setZero();
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
       for(int j = 0; j < 3; j++){
           for(int k = 0; k < 9; k++){
              dE(9*FaceID[j]+k) += (double) (simRest.fixDOF(9*FaceID[j]+k)) * g(9*j+k);
       //       std::cout <<  g(9*j+k) << " ";
           }
       }
   }
   
   return e;
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

