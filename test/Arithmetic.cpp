#include <Ravelin/MatrixNd.h>
#include <Ravelin/VectorNd.h>
#include <Ravelin/Vector3d.h>
#include <Ravelin/Matrix3d.h>
#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/Dense>
#include <cassert>

    bool checkError(Eigen::MatrixXd E, Ravelin::MatrixNd R){
        double TOL = 2e-17;
        double error = 0;
        for(int i=0;i<R.rows();i++){
            for(int j=0;j<R.columns();j++){
                double err = R(i,j) - E(i,j);
                error += err*err;
            }
        }
        return (TOL > error);
    }

    int main(){

        Eigen::MatrixXd E(10,10);
        Ravelin::MatrixNd R(10,10);
        E.setRandom();
        R.set_zero();

        // Check if the error between the matricies is small (fail)
        assert(checkError(E,  R));

        return 0;
    }

	// Run Unit Tests

