#ifndef BILAYERENERGY_H
#define BILAYERENERGY_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include "simSetup.h"

void hessiantri(double dAtop, double dA, double L1top, double L2top, double L1, double L2, double htop, double hbot, std::vector<Eigen::Vector3d> d, std::vector<Eigen::Vector3d> dtop, std::vector<Eigen::Vector3d> v, Eigen::Matrix3d ainvtop, Eigen::Matrix3d ainv, Eigen::MatrixXd & hessian);
double energytri(double dAtop, double dA, double L1top, double L2top, double L1, double L2, double htop, double hbot, std::vector<Eigen::Vector3d> d, std::vector<Eigen::Vector3d> dtop, std::vector<Eigen::Vector3d> v, Eigen::Matrix3d ainvtop, Eigen::Matrix3d ainv);
void gradtri(double dAtop, double dA, double L1top, double L2top, double L1, double L2, double htop, double hbot, std::vector<Eigen::Vector3d> d, std::vector<Eigen::Vector3d> dtop, std::vector<Eigen::Vector3d> v, Eigen::Matrix3d ainvtop, Eigen::Matrix3d ainv, Eigen::VectorXd & gradient);
double energytot(const RestConfig & simRest, const Eigen::VectorXd & dofs, Eigen::VectorXd * dE, std::vector<Eigen::Triplet<double> > * hesslist);

#endif
