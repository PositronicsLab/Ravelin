#include <UnitTesting.hpp>
#include <gtest/gtest.h>

static const unsigned MIN_SIZE = 0, MAX_SIZE = 7;

    // MatrixNd*MatrixNd
    TEST(Arithmetic, MM_Mult)
    {
        // Loop Through all column and vector sizes
        for(int i=0;i<MAX_SIZE;i++){
             for(int j=0;j<MAX_SIZE;j++){
                 unsigned r = 1 << j, c = 1 << j;
                 // Determine matrix size
                 MatE E1(r,c);
                 E1.setRandom();
                 MatR R1 = asRavelin(E1);

                 for(int k = 0; k < 10; k++) {
                        unsigned ck = 1 << k;

                        MatE E2(c,ck);
                        E2.setRandom();
                        MatR R2 = asRavelin(E2);

                        // Perform multiplication on matricies
                        MatE resultE = E1 * E2;
                        MatR resultR(r,ck);
                        MatR::mult(R1,R2, resultR);

                        // Check if the error between the matricies is small
                        checkError(std::cerr, "MatrixNd*MatrixNd", resultE, resultR);
                }
            }
        }
    }

    // MatrixNd'*MatrixNd
    TEST(Arithmetic, MM_TMult)
    {
        // Loop Through all column and vector sizes
        for(int i=0;i<MAX_SIZE;i++){
             for(int j=0;j<MAX_SIZE;j++){
                 unsigned r = 1 << i, c = 1 << j;
                 // Determine matrix size
                 MatE E1(r,c);
                 E1.setRandom();
                 MatR R1 = asRavelin(E1);

                 for(int k = 0; k < 10; k++) {
                     unsigned ck = 1 << k;
                     MatE E2(r,ck);
                     E2.setRandom();
                     MatR R2 = asRavelin(E2);

                     // Perform multiplication on matricies
                     MatE resultE = E1.transpose() * E2;
                     MatR resultR(c,ck);
                     R1.transpose_mult(R2, resultR);

                     // Check if the error between the matricies is small
                     checkError(std::cerr, "MatrixNd\'*MatrixNd", resultE,  resultR);
                 }
            }
        }
    }

    // MatrixNd*VectorNd
    TEST(Arithmetic, MV_Mult)
    {
        for(int i=0;i<MAX_SIZE;i++){
             for(int j=0;j<MAX_SIZE;j++){
                 // Determine matrix size
                 MatE E1(1 << i,1 << j);
                 E1.setRandom();
                 MatR R1 = asRavelin(E1);

                 VecE E2(1 << j);
                 E2.setRandom();
                 VecR R2 = asRavelin(E2);

                 // Perform multiplication on matricies
                 MatE resultE = E1 * E2;
                 MatR resultR(1 << i,1);
                 MatR::mult(R1,R2, resultR);

                 // Check if the error between the matricies is small
                 checkError(std::cerr, "MatrixNd*VectorNd", resultE,  resultR);
            }
        }
    }

    // VectorNd'*VectorNd
    TEST(Arithmetic, VV_Dot)
    {
        for(int i=0;i<MAX_SIZE;i++){
            // Determine matrix size
            VecE E1(1 << i);
            E1.setRandom();
            VecR R1 = asRavelin(E1);
            VecE E2(1 << i);
            E2.setRandom();
            VecR R2 = asRavelin(E2);
            // Perform multiplication on matricies
            VecE resultE = E1.transpose()*E2;
            VecR resultR(1);
            resultR = R1.dot(R2);
            // Check if the error between the matricies is small
            checkError(std::cerr, "VectorNd\'*VectorNd", resultE,  resultR);
        }
    }

    // MatrixNd*VectorNd
    TEST(Arithmetic, VM_Mult)
    {
        for(int i=0;i<MAX_SIZE;i++){
             for(int j=0;j<MAX_SIZE;j++){
                 // Determine matrix size
                 MatE E1(1 << j,1 << i);
                 E1.setRandom();
                 MatR R1 = asRavelin(E1);

                 VecE E2(1 << j);
                 E2.setRandom();
                 VecR R2 = asRavelin(E2);

                 // Perform multiplication on matricies
                 MatE resultE = E1.transpose() * E2;
                 MatR resultR(1 << i,1);
                 R1.transpose_mult(R2, resultR);

                 // Check if the error between the matricies is small
                 checkError(std::cerr, "MatrixNd*VectorNd", resultE,  resultR);
             }
        }
    }

    // Addition Subtraction
    TEST(Arithmetic, MM_Add)
    {
        for(int i = 0; i < 10; i++) {
            for(int j = 0; j < 10; j++) {
                MatE E1(1 << i,1 << j);
                E1.setRandom();
                MatR R1 = asRavelin(E1);

                MatE E2(1 << i,1 << j);
                E2.setRandom();
                MatR resultR = asRavelin(E2);

                MatE resultE = E1 + E2;
                resultR += R1;
                checkError(std::cerr, "MatrixNd+MatrixNd", resultE,  resultR);
            }
        }
    }

    // VectorNd+VectorNd
    TEST(Arithmetic, VV_Add)
    {
        for(int i=0;i<MAX_SIZE;i++){
        // Determine matrix size
                VecE E1(1 << i);
                E1.setRandom();
                VecR R1 = asRavelin(E1);
                VecE E2(1 << i);
                E2.setRandom();
                VecR resultR = asRavelin(E2);
                // Perform multiplication on matricies
                VecE resultE = E1 + E2;
                resultR += R1;
                // Check if the error between the matricies is small
                checkError(std::cerr, "VectorNd+VectorNd", resultE,  resultR);
        }
    }

     // VectorNd-VectorNd
    TEST(Arithmetic, VV_Sub){
    for(int i=0;i<MAX_SIZE;i++){
        // Determine matrix size
            unsigned r = 1 << i;
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
            checkError(std::cerr, "VectorNd-VectorNd", resultE, resultR);
        }
    }
    

    // MatrixNd-MatrixNd
    TEST(Arithmetic, MM_Sub)
    {
        for(int i = 0; i < 10; i++) {
                for(int j = 0; j < 10; j++) {
                        MatE E1(1 << i,1 << j);
                        E1.setRandom();
                        MatR R1 = asRavelin(E1);

                        MatE E2(1 << i,1 << j);
                        E2.setRandom();
                        MatR resultR = asRavelin(E2);

                        MatE resultE = E2 - E1;
                        resultR -= R1;
                        checkError(std::cerr, "MatrixNd-MatrixNd", resultE,  resultR);

                }
        }
    }

    // Diag(VectorNd)*MatrixNd
    TEST(Arithmetic, Diag_Mult)
    {
        for(int i=0;i<MAX_SIZE;i++){
             int r = 1 << i;
             // Determine matrix size
//                 MatE E1;
             VecE e1(r);
             e1.setRandom();
//                 E1 = e1.asDiagonal();
             VecR r1 = asRavelin(e1);

             MatE E2(r,r);
             E2.setRandom();
             MatR R2 = asRavelin(E2);

             // Perform multiplication on matricies
             MatE resultE = e1.asDiagonal() * E2;
             MatR resultR(r,r);
             MatR::diag_mult(r1, R2, resultR);

             // Check if the error between the matricies is small
             checkError(std::cerr, "diag_mult()", resultE,  resultR);

             // Perform multiplication on matricies
             resultE = e1.asDiagonal() * E2.transpose();
             MatR::diag_mult_transpose(r1,R2, resultR);

             // Check if the error between the matricies is small
             checkError(std::cerr, "diag_mult_transpose()", resultE,  resultR);
        }
    }
