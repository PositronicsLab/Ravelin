    #include <UnitTesting.hpp>



    /// Checks the error between the final value of the Eigen and Ravelin Arithmetic 
    bool checkError(MatE E, MatR R){
        double error = 0; 
        for(int i=0;i<R.rows();i++){ 
            for(int j=0;j<R.columns();j++){ 
                double err = R(i,j) - E(i,j); 
                error += err*err; 
            } 
        } 
        std::cout << "e: " << error << std::endl;
        return (TOL > error);
    } 
    bool checkError(VecE E, VecR R){
        double error = 0; 
        for(int i=0;i<R.size();i++){ 
                double err = R[i] - E(i); 
                error += err*err; 
        } 
        std::cout << "e: " << error << std::endl;
        return (TOL > error);
    } 

    /// Creates a Ravelin Matrix equal to your random Eigen3 Matrix 
    MatR asRavelin(MatE E){
        int rows = E.rows(); 
        int cols = E.cols(); 
        MatR R(rows,cols);
        for(int i=0;i<rows;i++) 
            for(int j=0;j<cols;j++) 
                R(i,j) = E(i,j); 
        return R; 
    } 
 
    VecR asRavelin(VecE E){
        int rows = E.size(); 
        VecR R(rows);
        for(int i = 0; i < rows; i++) 
                R[i] = E(i);  
        return R; 
    } 
 

