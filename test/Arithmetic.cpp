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
    bool checkError(Eigen::VectorXd E, Ravelin::VectorNd R){
        double TOL = 2e-17;
        double error = 0;
        for(int i=0;i<R.size();i++){
                double err = R[i] - E(i);
                error += err*err;
        }
        std::cout << "l2-norm error = Ravelin-Eigen:\n" << error << std::endl;
        return (TOL > error);
    }
    bool checkError(Eigen::MatrixXf E, Ravelin::MatrixNf R){
        double TOL = 1e-17;
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
    bool checkError(Eigen::VectorXf E, Ravelin::VectorNf R){
        double TOL = 1e-17;
        double error = 0;
        for(int i=0;i<R.size();i++){
                double err = R[i] - E(i);
                error += err*err;
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

    Ravelin::MatrixNf asRavelin(Eigen::MatrixXf E){
	int rows = E.rows();
        int cols = E.cols();
        Ravelin::MatrixNf R(rows,cols);
        for(int i=0;i<rows;i++)
            for(int j=0;j<cols;j++)
                R(i,j) = E(i,j);
        return R;
    }
    Ravelin::VectorNd asRavelin(Eigen::VectorXd E){
        int rows = E.size();
        Ravelin::VectorNd R(rows);
	for(int i = 0; i < rows; i++)
		R[i] = E(i); 
        return R;
    }

    Ravelin::VectorNf asRavelin(Eigen::VectorXf E){
        int rows = E.size();
        Ravelin::VectorNf R(rows);
        for(int i=0;i<rows;i++)
                R[i] = E(i);
        return R;
    }
    // MatrixNd*MatrixNd
    void MDMDMult() {
        // Loop Through all column and vector sizes
        for(int i=0;i<10;i++){
             for(int j=0;j<10;j++){
                 // Determine matrix size
                 Eigen::MatrixXd E1(std::pow(2,i),std::pow(2,j));
                 E1.setRandom();
                 Ravelin::MatrixNd R1 = asRavelin(E1);

                 for(int k = 0; k < 10; k++) {

                        Eigen::MatrixXd E2(std::pow(2,j),std::pow(2,k));
                        E2.setRandom();
                        Ravelin::MatrixNd R2 = asRavelin(E2);

                        // Perform multiplication on matricies
                        Eigen::MatrixXd resultE = E1 * E2;
                        Ravelin::MatrixNd resultR(std::pow(2,i),std::pow(2,k));
                        Ravelin::MatrixNd::mult(R1,R2, resultR);

                        // Check if the error between the matricies is small
                        assert(checkError(resultE,  resultR));
                }
            }
        }
        std::cout << "MatrixNd*MatrixNd has passed unit testing..." << std::endl;
    }
    // MatrixNd*VectorNd
    void MDVDMult() {
        for(int i=0;i<10;i++){
             for(int j=0;j<10;j++){
                 // Determine matrix size
                 Eigen::MatrixXd E1(std::pow(2,i),std::pow(2,j));
                 E1.setRandom();
                 Ravelin::MatrixNd R1 = asRavelin(E1);

                 Eigen::VectorXd E2(std::pow(2,j));
                 E2.setRandom();
                 Ravelin::VectorNd R2 = asRavelin(E2);

                 // Perform multiplication on matricies
                 Eigen::MatrixXd resultE = E1 * E2;
                 Ravelin::MatrixNd resultR(std::pow(2,i),1);
                 Ravelin::MatrixNd::mult(R1,R2, resultR);

                 // Check if the error between the matricies is small
                 assert(checkError(resultE,  resultR));
            }
        }
        std::cout << "MatrixNd*MatrixNd has passed unit testing..." << std::endl;
    }
    // MatrixNf*VectorNf
    void MFVFMult() {
        for(int i=0;i<10;i++){
             for(int j=0;j<10;j++){
                 // Determine matrix size
                 Eigen::MatrixXf E1(std::pow(2,i),std::pow(2,j));
                 E1.setRandom();
                 Ravelin::MatrixNf R1 = asRavelin(E1);

                 Eigen::VectorXf E2(std::pow(2,j));
                 E2.setRandom();
                 Ravelin::VectorNf R2 = asRavelin(E2);

                 // Perform multiplication on matricies
                 Eigen::MatrixXf resultE = E1 * E2;
                 Ravelin::MatrixNf resultR(std::pow(2,i),1);
                 Ravelin::MatrixNf::mult(R1,R2, resultR);
		 std::cout << "i=" << i << "\tj=" << j << std::endl;
                 // Check if the error between the matricies is small
                 assert(checkError(resultE,  resultR));
            }
        }
        std::cout << "MatrixNd*MatrixNd has passed unit testing..." << std::endl;
    }
    // VectorNd*VectorNd
    void VDVDDot() {
	for(int i=0;i<10;i++){
        // Determine matrix size
                Eigen::VectorXd E1(std::pow(2,i));
                E1.setRandom();
                Ravelin::VectorNd R1 = asRavelin(E1);
                Eigen::VectorXd E2(std::pow(2,i));
                E2.setRandom();
                Ravelin::VectorNd R2 = asRavelin(E2);
                // Perform multiplication on matricies
                Eigen::VectorXd resultE = E1.transpose()*E2;
                Ravelin::VectorNd resultR(1);
		resultR = R1.dot(R2);
                // Check if the error between the matricies is small
                assert(checkError(resultE,  resultR));
        }
        std::cout << "VectorNd+VectorNd has passed unit testing..." << std::endl;
    }
    void VFVFDot() {
	for(int i=0;i<10;i++){
        // Determine matrix size
                Eigen::VectorXf E1(std::pow(2,i));
                E1.setRandom();
                Ravelin::VectorNf R1 = asRavelin(E1);
                Eigen::VectorXf E2(std::pow(2,i));
                E2.setRandom();
                Ravelin::VectorNf R2 = asRavelin(E2);
                // Perform multiplication on matricies
                Eigen::VectorXf resultE = E1.transpose()*E2;
                Ravelin::VectorNf resultR(1);
                resultR = R1.dot(R2);
                // Check if the error between the matricies is small
                assert(checkError(resultE,  resultR));
        }
        std::cout << "VectorNd+VectorNd has passed unit testing..." << std::endl;
    }
    // MatrixNf*MatrixNf
    void MFMFMult() {
        for(int i=0;i<10;i++){
             for(int j=0;j<10;j++){
                 // Determine matrix size
                 Eigen::MatrixXf E1(std::pow(2,i),std::pow(2,j));
                 E1.setRandom();
                 Ravelin::MatrixNf R1 = asRavelin(E1);
                         for(int k = 0; k < 10; k++) {
                                Eigen::MatrixXf E2(std::pow(2,j),std::pow(2,k));
                        E2.setRandom();
                        Ravelin::MatrixNf R2 = asRavelin(E2);
                        // Perform multiplication on matricies
                        Eigen::MatrixXf resultE = E1 * E2;
                        Ravelin::MatrixNf resultR(std::pow(2,i),std::pow(2,k));
                        Ravelin::MatrixNf::mult(R1,R2, resultR);
                        std::cout << "i=" << i << "\tj=" << j << "\tk=" << k << std::endl;
                        // Check if the error between the matricies is small
                        assert(checkError(resultE,  resultR));
                }
            }
        }
        std::cout << "MatrixNf*MatrixNf has passed unit testing..." << std::endl;
    }
    void MDMDAdd() {
        for(int i = 0; i < 10; i++) {
                for(int j = 0; j < 10; j++) {
                        Eigen::MatrixXd E1(std::pow(2,i),std::pow(2,j));
                        E1.setRandom();
                        Ravelin::MatrixNd R1 = asRavelin(E1);

                        Eigen::MatrixXd E2(std::pow(2,i),std::pow(2,j));
                        E2.setRandom();
                        Ravelin::MatrixNd resultR = asRavelin(E2);

                        Eigen::MatrixXd resultE = E1 + E2;
                        resultR += R1;
                        assert(checkError(resultE, resultR));
                }
        }
        std::cout << "MatrixNd+MatrixNd has passed unit testing..." << std::endl;
    }
    // MatrixNf+MatrixNf
    void MFMFAdd() {
        for(int i = 0; i < 10; i++) {
                for(int j = 0; j < 10; j++) {
                        Eigen::MatrixXf E1(std::pow(2,i),std::pow(2,j));
                        E1.setRandom();
                        Ravelin::MatrixNf R1 = asRavelin(E1);

                        Eigen::MatrixXf E2(std::pow(2,i),std::pow(2,j));
                        E2.setRandom();
                        Ravelin::MatrixNf resultR = asRavelin(E2);

                        Eigen::MatrixXf resultE = E1 + E2;
                        resultR += R1;
                        assert(checkError(resultE, resultR));
                }
        }
        std::cout << "MatrixNf+MatrixNf has passed unit testing..." << std::endl;
    }
    // MatrixNd-MatrixNd
    void MDMDSub() {
        for(int i = 0; i < 10; i++) {
                for(int j = 0; j < 10; j++) {
                        Eigen::MatrixXd E1(std::pow(2,i),std::pow(2,j));
                        E1.setRandom();
                        Ravelin::MatrixNd R1 = asRavelin(E1);

                        Eigen::MatrixXd E2(std::pow(2,i),std::pow(2,j));
                        E2.setRandom();
                        Ravelin::MatrixNd resultR = asRavelin(E2);

                        Eigen::MatrixXd resultE = E2 - E1;
                        resultR -= R1;
                        assert(checkError(resultE, resultR));

                }
        }
        std::cout << "MatrixNd-MatrixNd has passed unit testing..." << std::endl;
    }
    // MatrixNf-MatrixNf
    void MFMFSub() {
        for(int i = 0; i < 10; i++) {
                for(int j = 0; j < 10; j++) {
                        Eigen::MatrixXf E1(std::pow(2,i),std::pow(2,j));
                        E1.setRandom();
                        Ravelin::MatrixNf R1 = asRavelin(E1);

                        Eigen::MatrixXf E2(std::pow(2,i),std::pow(2,j));
                        E2.setRandom();
                        Ravelin::MatrixNf resultR = asRavelin(E2);

                        Eigen::MatrixXf resultE = E2 - E1;
                        resultR -= R1;
                        assert(checkError(resultE, resultR));

                }
        }
        std::cout << "MatrixNf-MatrixNf has passed unit testing..." << std::endl;
    }
    // VectorNd+VectorNd
    void VDVDAdd() {
        for(int i=0;i<10;i++){
        // Determine matrix size
                Eigen::VectorXd E1(std::pow(2,i));
                E1.setRandom();
                Ravelin::VectorNd R1 = asRavelin(E1);
                Eigen::VectorXd E2(std::pow(2,i));
                E2.setRandom();
                Ravelin::VectorNd resultR = asRavelin(E2);
                // Perform multiplication on matricies
                Eigen::VectorXd resultE = E1 + E2;
                resultR += R1;
                // Check if the error between the matricies is small
                assert(checkError(resultE,  resultR));
        }
        std::cout << "VectorNd+VectorNd has passed unit testing..." << std::endl;
    }
    // VectorNf+VectorNf
    void VFVFAdd() {
        for(int i=0;i<10;i++){
        // Determine matrix size
                Eigen::VectorXf E1(std::pow(2,i));
                E1.setRandom();
                Ravelin::VectorNf R1 = asRavelin(E1);
                Eigen::VectorXf E2(std::pow(2,i));
                E2.setRandom();
                Ravelin::VectorNf resultR = asRavelin(E2);
                // Perform multiplication on matricies
                Eigen::VectorXf resultE = E1 + E2;
                resultR += R1;
                // Check if the error between the matricies is small
                assert(checkError(resultE,  resultR));
        }
        std::cout << "VectorNf+VectorNf has passed unit testing..." << std::endl;
    }
     // VectorNd-VectorNd
    void VDVDSub() {
	for(int i=0;i<10;i++){
        // Determine matrix size
                Eigen::VectorXd E1(std::pow(2,i));
                E1.setRandom();
                Ravelin::VectorNd R1 = asRavelin(E1);
                Eigen::VectorXd E2(std::pow(2,i));
                E2.setRandom();
                Ravelin::VectorNd resultR = asRavelin(E2);
                // Perform multiplication on matricies
                Eigen::VectorXd resultE = E1 - E2;
                resultR -= R1;
                // Check if the error between the matricies is small
                assert(checkError(resultE,  resultR));
        }
        std::cout << "VectorNd+VectorNd has passed unit testing..." << std::endl;
    }
    // VectorNf-VectorNf
    void VFVFSub() {
	for(int i=0;i<10;i++){
        // Determine matrix size
                Eigen::VectorXf E1(std::pow(2,i));
                E1.setRandom();
                Ravelin::VectorNf R1 = asRavelin(E1);
                Eigen::VectorXf E2(std::pow(2,i));
                E2.setRandom();
                Ravelin::VectorNf resultR = asRavelin(E2);
                // Perform multiplication on matricies
                Eigen::VectorXf resultE = E1 - E2;
                resultR -= R1;
                // Check if the error between the matricies is small
                assert(checkError(resultE,  resultR));
        }
        std::cout << "VectorNf+VectorNf has passed unit testing..." << std::endl;

    }
    
    /// Run Unit Tests
    int main() {
	MDMDMult();
	MDVDMult();
//	MFMFMult();
//	MFVFMult();
	VDVDDot();
	VFVFDot();
	VDVDAdd();
	VFVFAdd();
	MDMDAdd();
	MFMFAdd();
	VDVDSub();
	VFVFSub();
	MDMDSub();
	MFMFSub();
	return 0;
    }

