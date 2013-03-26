#ifndef _UNITTESTING_HPP_
#define _UNITTESTING_HPP_
#include <Ravelin/MatrixNd.h>
#include <Ravelin/MatrixNf.h>
#include <Ravelin/VectorNd.h>
#include <Ravelin/VectorNf.h>
#include <Ravelin/Vector3d.h>
#include <Ravelin/Vector3f.h>
#include <Ravelin/Matrix3d.h>
#include <Ravelin/Matrix3f.h>
#include <Ravelin/Operators.h>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <cassert>

#define SINGLE_PRECISION

#ifdef SINGLE_PRECISION
    typedef Eigen::MatrixXf MatE;
    typedef Eigen::VectorXf VecE;
    typedef Ravelin::MatrixNf MatR;
    typedef Ravelin::VectorNf VecR;
    const double TOL = 2e-6;
#else
    typedef Eigen::MatrixXd MatE;
    typedef Eigen::VectorXd VecE;
    typedef Ravelin::MatrixNd MatR;
    typedef Ravelin::VectorNd VecR;
    const double TOL = 2e-16;
#endif

    /// Checks the error between the final value of the Eigen and Ravelin Arithmetic
    extern bool checkError(MatE E, MatR R);
    extern bool checkError(VecE E, VecR R);
    extern MatR asRavelin(MatE E);
    extern VecR asRavelin(VecE E);

#endif //_UNITTESTING_HPP_
