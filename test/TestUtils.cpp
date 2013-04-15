    #include <UnitTesting.hpp>



    /// Checks the error between the final value of the Eigen and Ravelin Arithmetic 
    bool checkError(MatE E, MatR R)
    {
        double error = 0; 
        for(int i=0;i<R.rows();i++){ 
            for(int j=0;j<R.columns();j++){ 
                double err = R(i,j) - E(i,j); 
                error += err*err; 
            } 
        } 
//        std::cout << "e: " << error << std::endl;
        return (TOL > error);
    } 

    bool checkError(VecE E, VecR R)
    {
        double error = 0; 
        for(int i=0;i<R.size();i++){ 
                double err = R[i] - E(i); 
                error += err*err; 
        } 
//        std::cout << "e: " << error << std::endl;
        return (TOL > error);
    } 

    bool checkError(MatR E, MatR R)
    {
        double error = 0;
        for(int i=0;i<R.rows();i++){
            for(int j=0;j<R.columns();j++){
                double err = R(i,j) - E(i,j);
                error += err*err;
            }
        }
//        std::cout << "e: " << error << std::endl;
        return (TOL > error);
    }

    bool checkError(VecR E, VecR R)
    {
        double error = 0;
        for(int i=0;i<R.size();i++){
                double err = R[i] - E[i];
                error += err*err;
        }
//        std::cout << "e: " << error << std::endl;
        return (TOL > error);
    }

    /// Creates a Ravelin Matrix from an Eigen Matrix
    MatR asRavelin(MatE E)
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
                R(i,j) = double(std::rand()) / double(RAND_MAX);
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
    VecR asRavelin(VecE E)
    {
        int rows = E.size(); 
        VecR R(rows);
        for(int i = 0; i < rows; i++) 
                R[i] = E(i);  
        return R; 
    } 

    /// Creates a Ravelin Vector from an Eigen Vector
    VecR toVec(MatR E)
    {
        int rows = E.rows();
        VecR R(rows);
        for(int i = 0; i < rows; i++)
                R[i] = E(i,1);
        return R;
    }

    /// Creates a Ravelin Vector from an Eigen Vector
    MatR toMat(VecR E)
    {
        int rows = E.rows();
        MatR R(rows,1);
        for(int i = 0; i < rows; i++)
                R(i,1) = E[i];
        return R;
    }
 

