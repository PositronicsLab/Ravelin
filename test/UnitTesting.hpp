#ifndef _UNITTESTING_HPP_
#define _UNITTESTING_HPP_

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <cassert>

//#define SINGLE_PRECISION

#ifdef SINGLE_PRECISION

   const static double ZERO_TOL = std::numeric_limits<float>::epsilon();

    #include <Ravelin/MatrixNf.h>
    #include <Ravelin/VectorNf.h>
    #include <Ravelin/Vector3f.h>
    #include <Ravelin/Matrix3f.h>

    typedef Eigen::MatrixXf MatE;
    typedef Eigen::VectorXf VecE;
    typedef Ravelin::MatrixNf MatR;
    typedef Ravelin::VectorNf VecR;
    const double TOL = 2e-6;
#else
    const static double ZERO_TOL = std::numeric_limits<double>::epsilon();

    #include <Ravelin/MatrixNd.h>
    #include <Ravelin/Vector3d.h>
    #include <Ravelin/VectorNd.h>
    #include <Ravelin/Matrix3d.h>

    typedef Eigen::MatrixXd MatE;
    typedef Eigen::VectorXd VecE;
    typedef Ravelin::MatrixNd MatR;
    typedef Ravelin::VectorNd VecR;
    const double TOL = 2e-16;
#endif

    const static double NEAR_ZERO = sqrt(ZERO_TOL);

    /// Test Utils
    bool checkError(MatE E, MatR R);
    bool checkError(VecE E, VecR R);
    bool checkError(VecR E, VecR R);
    bool checkError(MatR E, MatR R);
    MatR randM(unsigned r, unsigned c);
    VecR randV(unsigned r);
    MatR asRavelin(MatE E);
    VecR asRavelin(VecE E);

#endif //_UNITTESTING_HPP_
