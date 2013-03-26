#include <UnitTesting.hpp>

const unsigned MIN_SIZE = 0, MAX_SIZE = 9;

    // MatrixNd*MatrixNd
    void MM_Mult()
    {
        std::cout << ">> MatrixNd*MatrixNd" << std::endl;
        // Loop Through all column and vector sizes
        for(int i=0;i<MAX_SIZE;i++){
             for(int j=0;j<MAX_SIZE;j++){
                 unsigned r = std::pow(2,i), c = std::pow(2,j);
                 // Determine matrix size
                 MatE E1(r,c);
                 E1.setRandom();
                 MatR R1 = asRavelin(E1);

                 for(int k = 0; k < 10; k++) {
                        unsigned ck = std::pow(2,k);

                        MatE E2(c,ck);
                        E2.setRandom();
                        MatR R2 = asRavelin(E2);

                        // Perform multiplication on matricies
                        MatE resultE = E1 * E2;
                        MatR resultR(r,ck);
                        MatR::mult(R1,R2, resultR);

                        // Check if the error between the matricies is small
                        std::cout << r << "x" << c << "*" << c << "x" << ck << std::endl;
                        if(!checkError(resultE,  resultR)){
                            assert(false);
                        }
                }
            }
        }
        std::cout << "<< MatrixNd*MatrixNd has passed unit testing..." << std::endl;
    }

    // MatrixNd'*MatrixNd
    void MM_TMult()
    {
        std::cout << ">> MatrixNd*MatrixNd\' ..." << std::endl;
        // Loop Through all column and vector sizes
        for(int i=0;i<MAX_SIZE;i++){
             for(int j=0;j<MAX_SIZE;j++){
                 unsigned r = std::pow(2,i), c = std::pow(2,j);
                 // Determine matrix size
                 MatE E1(r,c);
                 E1.setRandom();
                 MatR R1 = asRavelin(E1);

                 for(int k = 0; k < 10; k++) {
                     unsigned ck = std::pow(2,k);
                     MatE E2(r,ck);
                     E2.setRandom();
                     MatR R2 = asRavelin(E2);

                     // Perform multiplication on matricies
                     MatE resultE = E1.transpose() * E2;
                     MatR resultR(c,ck);
                     R1.transpose_mult(R2, resultR);

                     // Check if the error between the matricies is small
                     std::cout << r << "x" << c << "\' *" << r << "x" << ck << std::endl;
                     if(!checkError(resultE,  resultR)){
                         assert(false);
                     }
                 }
            }
        }
        std::cout << "<< MatrixNd*MatrixNd\' has passed unit testing..." << std::endl;
    }

    // MatrixNd*VectorNd
    void MV_Mult()
    {
        for(int i=0;i<MAX_SIZE;i++){
             for(int j=0;j<MAX_SIZE;j++){
                 // Determine matrix size
                 MatE E1(std::pow(2,i),std::pow(2,j));
                 E1.setRandom();
                 MatR R1 = asRavelin(E1);

                 VecE E2(std::pow(2,j));
                 E2.setRandom();
                 VecR R2 = asRavelin(E2);

                 // Perform multiplication on matricies
                 MatE resultE = E1 * E2;
                 MatR resultR(std::pow(2,i),1);
                 MatR::mult(R1,R2, resultR);

                 // Check if the error between the matricies is small
                 assert(checkError(resultE,  resultR));
            }
        }
        std::cout << "MatrixNd*VectorNd has passed unit testing..." << std::endl;
    }

    // VectorNd*VectorNd
    void VV_Dot()
    {
        for(int i=0;i<MAX_SIZE;i++){
            // Determine matrix size
            VecE E1(std::pow(2,i));
            E1.setRandom();
            VecR R1 = asRavelin(E1);
            VecE E2(std::pow(2,i));
            E2.setRandom();
            VecR R2 = asRavelin(E2);
            // Perform multiplication on matricies
            VecE resultE = E1.transpose()*E2;
            VecR resultR(1);
            resultR = R1.dot(R2);
            // Check if the error between the matricies is small
            assert(checkError(resultE,  resultR));
        }
        std::cout << "VectorNd dot VectorNd has passed unit testing..." << std::endl;
    }

    // MatrixNd*VectorNd
    void VM_Mult()
    {
        for(int i=0;i<MAX_SIZE;i++){
             for(int j=0;j<MAX_SIZE;j++){
                 // Determine matrix size
                 MatE E1(std::pow(2,j),std::pow(2,i));
                 E1.setRandom();
                 MatR R1 = asRavelin(E1);

                 VecE E2(std::pow(2,j));
                 E2.setRandom();
                 VecR R2 = asRavelin(E2);

                 // Perform multiplication on matricies
                 MatE resultE = E1.transpose() * E2;
                 MatR resultR(std::pow(2,i),1);
                 R1.transpose_mult(R2, resultR);

                 // Check if the error between the matricies is small
                 assert(checkError(resultE,  resultR));
            }
        }
        std::cout << "MatrixNd*VectorNd has passed unit testing..." << std::endl;
    }

    // VectorNd*VectorNd'
//    void VDVDOuterProd()
//    {
//        for(int i=0;i<MAX_SIZE;i++){
//            // Determine matrix size
//            VecE E1(std::pow(2,i));
//            E1.setRandom();
//            VecR R1 = asRavelin(E1);
//            VecE E2(std::pow(2,i));
//            E2.setRandom();
//            VecR R2 = asRavelin(E2);
//            // Perform multiplication on matricies
//            VecE resultE = E1*E2.transpose();
//            VecR resultR(1);
//            R1 R2,resultR);
//            // Check if the error between the matricies is small
//            assert(checkError(resultE,  resultR));
//        }
//        std::cout << "VectorNd dot VectorNd has passed unit testing..." << std::endl;
//    }

    // Addition Subtraction
    void MM_Add() {
        for(int i = 0; i < 10; i++) {
            for(int j = 0; j < 10; j++) {
                MatE E1(std::pow(2,i),std::pow(2,j));
                E1.setRandom();
                MatR R1 = asRavelin(E1);

                MatE E2(std::pow(2,i),std::pow(2,j));
                E2.setRandom();
                MatR resultR = asRavelin(E2);

                MatE resultE = E1 + E2;
                resultR += R1;
                assert(checkError(resultE, resultR));
            }
        }
        std::cout << "MatrixNd+MatrixNd has passed unit testing..." << std::endl;
    }

    // VectorNd+VectorNd
    void VV_Add()
    {
        for(int i=0;i<MAX_SIZE;i++){
        // Determine matrix size
                VecE E1(std::pow(2,i));
                E1.setRandom();
                VecR R1 = asRavelin(E1);
                VecE E2(std::pow(2,i));
                E2.setRandom();
                VecR resultR = asRavelin(E2);
                // Perform multiplication on matricies
                VecE resultE = E1 + E2;
                resultR += R1;
                // Check if the error between the matricies is small
                assert(checkError(resultE,  resultR));
        }
        std::cout << "VectorNd+VectorNd has passed unit testing..." << std::endl;
    }

     // VectorNd-VectorNd
    void VV_Sub() {
    for(int i=0;i<MAX_SIZE;i++){
        // Determine matrix size
            unsigned r = std::pow(2,i);
            VecE E1(r);
            E1.setRandom();
            VecR resultR = asRavelin(E1);
            VecE E2(r);
            E2.setRandom();
            VecR R2 = asRavelin(E2);
            // Perform multiplication on matricies
            VecE resultE = E1 - E2;
            resultR -= R2;
            // Check if the error between the matricies is small
            assert(checkError(resultE,  resultR));
        }
        std::cout << "VectorNd+VectorNd has passed unit testing..." << std::endl;
    }
    

    // MatrixNd-MatrixNd
    void MM_Sub()
    {
        for(int i = 0; i < 10; i++) {
                for(int j = 0; j < 10; j++) {
                        MatE E1(std::pow(2,i),std::pow(2,j));
                        E1.setRandom();
                        MatR R1 = asRavelin(E1);

                        MatE E2(std::pow(2,i),std::pow(2,j));
                        E2.setRandom();
                        MatR resultR = asRavelin(E2);

                        MatE resultE = E2 - E1;
                        resultR -= R1;
                        assert(checkError(resultE, resultR));

                }
        }
        std::cout << "MatrixNd-MatrixNd has passed unit testing..." << std::endl;
    }

    /// Run Unit Tests
int main() {
    MM_Mult();
    MM_TMult();
    MV_Mult();
    VM_Mult();
    VV_Dot();
    VV_Add();
    MM_Add();
    VV_Sub();
    MM_Sub();
	return 0;
}

