#ifndef SIMSETUP_H
#define SIMSETUP_H
#include <Eigen/Dense>

struct RestConfig{
       Eigen::MatrixXd V;
       Eigen::MatrixXi F;
       Eigen::VectorXi fixDOF;
       std::vector<Eigen::Matrix3d> ainvs;
       std::vector<double> dAs;
       double L1top;
       double L2top;
       double L1;
       double L2;
       double hbot;
       double htop;
};

Eigen::Matrix3d calcAbar(const Eigen::Vector3d &p0, const Eigen::Vector3d &p1, const Eigen::Vector3d &p2);
double LameAlpha(double YoungsModulus, double PoissonsRatio);
double LameBeta(double YoungsModulus, double PoissonsRatio);
void setupDOF(const Eigen::MatrixXd & V, const double compression, Eigen::VectorXd & dof);
void genConstraints(std::string str, Eigen::VectorXi & fixDOF, double nVerts);

#endif
