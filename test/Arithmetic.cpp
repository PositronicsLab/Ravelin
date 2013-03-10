#include <Ravelin/MatrixNd.h>
#include <Ravelin/VectorNd.h>
#include <Ravelin/Vector3d.h>
#include <Ravelin/Matrix3d.h>
#include <Ravelin/Operators.h>
#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/Dense>
#include <cassert>

    /// Checks the error between the final value of the Eigen and Ravelin Arithmetic
    bool checkError(Eigen::MatrixXd E, Ravelin::MatrixNd R){
        double TOL = 2e-17;
        double error = 0;
        for(int i=0;i<R.rows();i++){
            for(int j=0;j<R.columns();j++){
                double err = R(i,j) - E(i,j);
                error += err*err;
            }
        }
        std::cout << "l2-norm error = Ravelin-Eigen:\n" << error << std::endl;
        return (TOL > error);
    }

    /// Creates a Ravelin Matrix equal to your random Eigen3 Matrix
    Ravelin::MatrixNd asRavelin(Eigen::MatrixXd E){
        int rows = E.rows();
        int cols = E.cols();
        Ravelin::MatrixNd R(rows,cols);
        for(int i=0;i<rows;i++)
            for(int j=0;j<cols;j++)
                R(i,j) = E(i,j);
        return R;
    }

    /// Run Unit Tests
    int main(){

        // Loop Through all column and vector sizes
        for(int i=0;i<10;i++){
             for(int j=0;j<10;j++){
                 // Determine matrix size
                 Eigen::MatrixXd E1(i,j);
                 E1.setRandom();
                 Ravelin::MatrixNd R1 = asRavelin(E1);

                 Eigen::MatrixXd E2(j,i);
                 E2.setRandom();
                 Ravelin::MatrixNd R2 = asRavelin(E2);

                 // Perform operation on matricies
                 Eigen::MatrixXd resultE = E1 * E2;
                 Ravelin::MatrixNd resultR(i,i);
                 Ravelin::MatrixNd::mult(R1,R2,resultR);

                 // Check if the error between the matricies is small
                 assert(checkError(resultE,  resultR));
            }
        }
        std::cout << "Multiplication has passed unit testing..." << std::endl;

        return 0;
    }


