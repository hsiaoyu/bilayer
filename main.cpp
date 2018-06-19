#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Sparse>
//#include <Eigen/Eigenvalues>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <fstream>
#include "simSetup.h"
#include "bilayerenergy.h"

void OutputSim(std::string outName, int outidx, int nVerts, double hbot, double htop, Eigen::MatrixXi & Fin, Eigen::VectorXd & simDof);
void zeroout(Eigen::MatrixXd & M, const Eigen::VectorXi & fixdof);

int main(int argc, char *argv[]){
   if(argc < 2){
       std::cout << "Need to pass in input file!" << std::endl;
       return 1;
    }
   
    std::ifstream infile(argv[1]);
    std::string str;
    getline(infile, str);
    RestConfig simRest;
    igl::readOBJ(str.c_str(), simRest.V, simRest.F);
    std::string outName;
    outName = str.substr(0,str.size()-4);
    outName += "CPD";  
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
    double linestep;
    std::vector<Eigen::Triplet<double> > hesslist;
    Eigen::VectorXd dE(9 * simRest.V.rows()), guess(9 * simRest.V.rows());
    while(residual > tol){
           /*
           //----------------Gradient Descent-------------------------
           double Energy = energytot(simRest, simDof, dE, NULL);
           //for(int i = 0; i < simRest.V.rows(); i ++){
           //    std::cout << "total grad  " << i <<  std::endl; 
           //    for(int j = 0; j < 9; j ++)
           //         std::cout <<  dE(9*i+j) << " ";
           //    std::cout << std::endl; 
           //}
           linestep = stepsize;
           guess = simDof - linestep * dE;
           Eigen::VectorXd dEguess(9 * simRest.V.rows());
           double guessE = energytot(simRest, guess, dEguess, NULL);
           while (guessE > Energy - lsalpha * linestep * dE.norm() * dE.norm()){
           //while (guessE > Energy){
                 linestep  = beta * linestep;
                 std::cout << "current linestep = " << linestep << std::endl;
                 guess = simDof - linestep * dE;
                 guessE = energytot(simRest, guess, dEguess, NULL);
           }
           residual = dEguess.norm();
           simDof = guess;
           //-------------------------------------------------------------
           */
           
           linestep = stepsize;
           double Energy = energytot(simRest, simDof, &dE, &hesslist);
           Eigen::SparseMatrix<double> hE(9*simRest.V.rows(), 9*simRest.V.rows());
           hE.setFromTriplets(hesslist.begin(), hesslist.end());

           //------------------Sparse Solve---------------------
           Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
           solver.compute(hE);
           //Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver(hE);
           Eigen::VectorXd upd = solver.solve(dE);

           guess = simDof - linestep * upd;
           double guessE = energytot(simRest, guess, NULL, NULL);
           while (guessE > Energy - lsalpha * linestep * dE.norm() * dE.norm()){
                 linestep  = beta * linestep;
                 std::cout << "current linestep = " << linestep << std::endl;
                 guess = simDof - linestep * upd;
                 guessE = energytot(simRest, guess, NULL, NULL);
           }

           guessE = energytot(simRest, guess, &dE, NULL);
           residual = dE.norm();
           simDof = guess;
           
           std::cout << "Energy = " << guessE << " |dE| = " << residual <<std::endl;
       
           if(residual < tol){
              OutputSim(outName, count, simRest.V.rows(), simRest.hbot, simRest.htop, simRest.F, simDof);
              Eigen::MatrixXd hEdense = hE;
              Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(hEdense);
              Eigen::MatrixXd D = (es.eigenvalues().cwiseAbs()).asDiagonal();
              Eigen::MatrixXd U = es.eigenvectors();
              Eigen::MatrixXd hEabs = U * D * U.inverse();
              zeroout(hEabs, simRest.fixDOF);
              upd = hEabs.ldlt().solve(dE);
              guess = simDof - linestep * upd;
              double guessE = energytot(simRest, guess, &dE, NULL);
              residual = dE.norm();
              simDof = guess;
              
              count ++;
              std::cout << "abs hessian convertion" << std::endl << "Energy = " << guessE << " |dE| = " << residual <<std::endl;
           }

    }
    
    //Eigen::MatrixXd Vout(simRest.V.rows() * 3, 3);
    //Eigen::MatrixXi Fout(simRest.F.rows() * 3 ,3);
    //Fout.block(0,0,simRest.F.rows(),3)= simRest.F;   

    //for(int i = 0; i < simRest.V.rows(); i++){
    //    for(int j = 0; j < 3; j++){
    //        Vout(i,j) = simDof(9*i+j);
    //        Vout(i+simRest.V.rows(),j) = simDof(9*i+j)+simRest.hbot*simDof(9*i+3+j);
    //        Vout(i+2*simRest.V.rows(),j) = simDof(9*i+j)+simRest.hbot*simDof(9*i+3+j)+simRest.htop*simDof(9*i+6+j);
    //    }
    //}
   
    //int nVerts = simRest.V.rows();
    //for(int i = 0; i < simRest.F.rows(); i++){
    //    Eigen::RowVector3i nv(nVerts, nVerts, nVerts);
    //    Fout.row(i+simRest.F.rows()) = simRest.F.row(i) + nv;
    //    Fout.row(i+2*simRest.F.rows()) = simRest.F.row(i) + 2 * nv;
    //} 
  
    //igl::writeOBJ(outName, Vout, Fout); 
 
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

void zeroout(Eigen::MatrixXd & M, const Eigen::VectorXi & fixdof){
     for (int i = 0; i < fixdof.size(); i++){
          if(fixdof(i) == 0){
             M.row(i).setZero();
             M.col(i).setZero();
          }
     }
}

void OutputSim(std::string outName, int outidx, int nVerts, double hbot, double htop, Eigen::MatrixXi & Fin, Eigen::VectorXd & simDof){
    Eigen::MatrixXd Vout(nVerts * 3, 3);
    Eigen::MatrixXi Fout(Fin.rows() * 3 ,3);
    Fout.block(0,0,Fin.rows(),3)= Fin;   

    for(int i = 0; i < nVerts; i++){
        for(int j = 0; j < 3; j++){
            Vout(i,j) = simDof(9*i+j);
            Vout(i+nVerts,j) = simDof(9*i+j) + hbot * simDof(9*i+3+j);
            Vout(i+2*nVerts,j) = simDof(9*i+j) + hbot * simDof(9*i+3+j) + htop * simDof(9*i+6+j);
        }
    }
   
    for(int i = 0; i < Fin.rows(); i++){
        Eigen::RowVector3i nv(nVerts, nVerts, nVerts);
        Fout.row(i+Fin.rows()) = Fin.row(i) + nv;
        Fout.row(i+2*Fin.rows()) = Fin.row(i) + 2 * nv;
    } 
  
    std::stringstream outname;
    outname << outName << "_" << outidx << ".obj"; 
    igl::writeOBJ(outname.str(), Vout, Fout); 
}

