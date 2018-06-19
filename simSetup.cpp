#include "simSetup.h"
#include <fstream>
#include <stdlib.h>

using namespace Eigen;

//Only valid for flat rest shape
Matrix3d calcAbar(const Vector3d &p0, const Vector3d &p1, const Vector3d &p2)
{
    Matrix3d result;
    result << (p1 - p0).dot(p1 - p0), (p1 - p0).dot(p2 - p0), 0,
        (p1 - p0).dot(p2 - p0), (p2 - p0).dot(p2 - p0),0,
        0, 0, 1;
    return result;
}

double LameAlpha(double YoungsModulus, double PoissonsRatio)
{
    //return YoungsModulus * PoissonsRatio / (1.0 + PoissonsRatio) / (1.0 - PoissonsRatio);
    return YoungsModulus * PoissonsRatio / (1.0 + PoissonsRatio) / (1.0 - 2.0 * PoissonsRatio);
}

double LameBeta(double YoungsModulus, double PoissonsRatio)
{
    //return YoungsModulus / (1.0 + PoissonsRatio);
    return  YoungsModulus / (2.0 * (1.0 + PoissonsRatio));
}

void setupDOF(const Eigen::MatrixXd & V, const double compression, Eigen::VectorXd & dof){
     dof.resize(V.rows()*9);
     dof.setZero();
     for(int i = 0; i < V.rows(); i++){
         for (int j = 0; j < 3; j++){
              dof(9*i+j) = V(i,j);
              if(j == 0)
                 dof(9*i+j) = compression * V(i,j);
         }
         //dof(9*i+5) = 1 ;
         //dof(9*i+8) = 1 ;
         dof(9*i+5) = 1 + 0.001 * rand()/(RAND_MAX + 1.0);
         dof(9*i+8) = 1 + 0.001 * rand()/(RAND_MAX + 1.0);
     }
}

void genConstraints(std::string str, Eigen::VectorXi & fixDOF, double nVerts){
    std::fstream infile;
    infile.open(str);
    fixDOF.setConstant(1);
    int nEdges, boundary, v1, v2;
    infile >> nEdges;     
    infile >> boundary;
    for(int i = 0; i < nEdges; i++){
       infile >> boundary;
       infile >> v1;
       infile >> v2;
       infile >> boundary;
       v1--;
       v2--;
       if(boundary > 0){
          fixDOF(9*v1+boundary-1) = 0;
          fixDOF(9*v1+3+boundary-1) = 0;
          fixDOF(9*v1+6+boundary-1) = 0;
          fixDOF(9*v2+boundary-1) = 0;
          fixDOF(9*v2+3+boundary-1) = 0;
          fixDOF(9*v2+6+boundary-1) = 0;
       }
    }
    
    for(int i =0; i <nVerts; i++)
       fixDOF(9*i+2) = 0;
}
