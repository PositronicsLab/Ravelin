#include <UnitTesting.hpp>
#include <iostream>

    /// Checks the error between the final value of the Eigen and Ravelin Arithmetic 
    double checkError(std::ostream& out, const std::string& str, const MatE& E, const MatR& R)
    {
        double error = 0; 
        for(int i=0;i<R.rows();i++){ 
            for(int j=0;j<R.columns();j++){ 
                double err = R(i,j) - E(i,j); 
                error = std::max(error, std::fabs(err));
            } 
        } 
        if (error > TOL*R.norm_inf()*std::max(R.rows(), R.columns()))
          out << "[FAIL] " << str << " [err=" << error << "]" << std::endl;
        else
          out << "[PASS] " << str << " [err=" << error << "]" << std::endl; 
        return error;
    } 

    double checkError(std::ostream& out, const std::string& str, const VecE& E, const VecR& R)
    {
        double error = 0; 
        for(int i=0;i<R.size();i++){ 
                double err = R[i] - E(i); 
                error = std::max(error, std::fabs(err));
        }
        if (error > TOL*R.norm_inf()*R.rows())
          out << "[FAIL] " << str << " [err=" << error << "]" << std::endl;
        else
          out << "[PASS] " << str << " [err=" << error << "]" << std::endl; 
        return error;
    } 
 
    /// Checks the error between the final value of Ravelin matrices 
    double checkError(std::ostream& out, const std::string& str, const MatR& E, const MatR& R)
    {
        double error = 0; 
        for(int i=0;i<R.rows();i++){ 
            for(int j=0;j<R.columns();j++){ 
                double err = R(i,j) - E(i,j); 
                error = std::max(error, std::fabs(err));
            } 
        } 
        if (error > TOL)
          out << "[FAIL] " << str << " [err=" << error << "]" << std::endl;
        else
          out << "[PASS] " << str << " [err=" << error << "]" << std::endl; 
        return error;
    } 

    /// Checks the error between the final value of Ravelin vectors 
    double checkError(std::ostream& out, const std::string& str, const VecR& E, const VecR& R)
    {
        double error = 0; 
        for(int i=0;i<R.size();i++){ 
                double err = R[i] - E[i]; 
                error += err*err; 
        }
        error = std::sqrt(error); 
        if (error > TOL)
          out << "[FAIL] " << str << " [err=" << error << "]" << std::endl;
        else
          out << "[PASS] " << str << " [err=" << error << "]" << std::endl; 
        return error;
    } 
 

    /// Creates a Ravelin Matrix from an Eigen Matrix
    MatR asRavelin(const MatE& E)
    {
        int rows = E.rows(); 
        int cols = E.cols(); 
        MatR R(rows,cols);
        for(int i=0;i<rows;i++) 
            for(int j=0;j<cols;j++) 
                R(i,j) = E(i,j); 
        return R; 
    } 

    /// A random Ravelin Matrix
    MatR randM(unsigned rows, unsigned cols)
    {
        MatR R(rows,cols);
        for(int i=0;i<rows;i++)
            for(int j=0;j<cols;j++)
                R(i,j) = (double) std::rand() / RAND_MAX;
        return R;
    }

    /// A random Ravelin Vector
    VecR randV(unsigned rows)
    {
        VecR R(rows);
        for(int i=0;i<rows;i++)
            R[i] = std::rand() % 10;
        return R;
    }

    /// Creates a Ravelin Vector from an Eigen Vector
    VecR asRavelin(const VecE& E)
    {
        int rows = E.size(); 
        VecR R(rows);
        for(int i = 0; i < rows; i++) 
                R[i] = E(i);  
        return R; 
    } 

    /// Creates a Ravelin Vector from an Eigen Vector
    VecR toVec(const MatR& E)
    {
        int rows = E.rows();
        VecR R(rows);
        for(int i = 0; i < rows; i++)
                R[i] = E(i,1);
        return R;
    }

    /// Creates a Ravelin Vector from an Eigen Vector
    MatR toMat(const VecR& E)
    {
        int rows = E.rows();
        MatR R(rows,1);
        for(int i = 0; i < rows; i++)
                R(i,1) = E[i];
        return R;
    }
 

