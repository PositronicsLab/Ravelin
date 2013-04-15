#include <UnitTesting.hpp>

static const unsigned MIN_SIZE = 0, MAX_SIZE = 7;

#ifdef SINGLE_PRECISION
    #include <Ravelin/LinAlgf.h>
    typedef Ravelin::LinAlgf LinAlg;
#else
    #include <Ravelin/LinAlgd.h>
    typedef Ravelin::LinAlgd LinAlg;
#endif

LinAlg * LA;
/*
// VectorNd*VectorNd'
void VV_OuterProd()
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
        MatE resultE = E1*E2.transpose();
        MatR resultR(std::pow(2,i),std::pow(2,i));
//        LA->outer_prod(,R2,resultR);
        // Check if the error between the matricies is small
        if(!checkError(resultE,  resultR)){
            std::cout << "<< [FAIL] VectorNd*VectorNd\'" << std::endl;
            return;
        }
    }
    std::cout << "<< [PASS] VectorNd*VectorNd\'" << std::endl;
}
*/
void testSVD(){
    std::cout << ">> testSVD: " << std::endl;

    for(int i=1;i<MAX_SIZE;i++){
        unsigned r = std::pow(2,i),c = 2;
        MatR A,B,AB(r,r);

        // Create PD matrix A
        A = randM(r,c);
        B = A;
        B.transpose();
        A.mult(B,AB);
        MatR U,V;
        VecR s;

        A = AB;
        /// Test the SVD1 Decomp
        LA->svd1(AB,U,s,V);
        MatR US(U.columns(),s.rows());

        MatR S(s.rows(),s.rows());
        S.set_zero();
        for(unsigned i=0;i<s.rows();i++)
            S(i,i) = s[i];

        U.mult(S,US);
        US.mult_transpose(V,AB);
        assert(checkError(A,AB));
        std::cout << "[PASS] svd1: " << std::endl;

        AB = A;
        /// Test the SVD2 Decomp
        LA->svd2(AB,U,s,V);
        S.set_zero();
        for(unsigned i=0;i<s.rows();i++)
            S(i,i) = s[i];

        U.mult(S,US);
        US.mult_transpose(V,AB);
        assert(checkError(A,AB));
        std::cout << "[PASS] svd2: " << std::endl;

        /// Test SVD Solve
        AB = A;
        LA->svd2(AB,U,s,V);

        MatR x = randM(r,1),b(r,1);
        //Ax1 = b
        A.mult(x,b);
        MatR xb = b;
        //(A^-1)b = x2
        LA->solve_LS_fast(U,s,V,xb);
        // x1 == x2 ?
        std::cout << x << std::endl;
        std::cout << xb << std::endl;
        assert(checkError(x,xb));
        std::cout << "[PASS] solve_LS_fast(U,S,V): " << std::endl;
    }
    std::cout << "<< [PASS] testSVD: " << std::endl;
}

/*
void testLU(){
    for(int i=1;i<MAX_SIZE;i++){
        unsigned s = std::pow(2,i),c = 2;
        MatR A,B,AB(s,s);

        // Create PD matrix A
        A = randM(s,c);
        B = A;
        B.transpose();
        A.mult(B,AB);
        A = AB;
        MatR U,V;
        VecR S;

        LA->factor_LU(AB);


        /// Test the SVD Decomp
        MatR SVT(S.rows(),V.rows());
        MatR::diag_mult_transpose(S,V,SVT);
        SVT.mult(U,AB);
        assert(checkError(A,AB));

        /// Test SVD Solve
        MatR x = randM(s,1),b(s,1);
        //Ax1 = b
        A.mult(x,b);
        MatR xb = b;
        //(A^-1)b = x2
        LA->solve_LS_fast(U,S,V,xb);
        // x1 == x2 ?
        assert(checkError(x,xb));
    }
}
*/

void testChol(){
    std::cout << ">> testChol: " << std::endl;

    for(int i=1;i<MAX_SIZE;i++){
        unsigned s = std::pow(2,i),c = s;
        MatR A,B,AB(s,s);

        std::cout << "SIZE: " << s << std::endl;

        /// Test SPSD
        bool PSD = false;
        while(!PSD){
            A = randM(s,c);
            B = A;
            B.transpose();
            A.mult(B,AB);
            A = AB;
            MatR U,V;
            VecR S;

            /// NOTE: TEST SVD
            LA->svd(AB,U,S,V);
            if((*std::min_element(S.begin(),S.end())) < (ZERO_TOL * S.rows() * (*std::max_element(S.begin(),S.end()))))
                continue;
            PSD = true;

            /// Test SPSD
                AB =A;
                LA->is_SPSD(AB,ZERO_TOL * S.rows() * (*std::max_element(S.begin(),S.end())));
                std::cout << " [PASS] is_SPSD(SPSD): " << std::endl;

            /// Test Eigen
                AB =A;
            VecR evals(S.rows());
            /// Test eigen decomp
                LA->eig_symm (AB,evals);
                // S should be same as eigenvalues
                checkError(evals,S);
                std::cout << " [PASS] eig_symm: " << std::endl;

            /// Test Eigen Vec...
//                eig_symm_plus (X &A_evecs, Y &evals);

        }


        // Create SPD matrix A
        bool PD = false;
        while(!PD){
            A = randM(s,c);
            B = A;
            B.transpose();
            A.mult(B,AB);
            A = AB;
            MatR U,V;
            VecR S;

            /// NOTE: TEST SVD
            LA->svd(AB,U,S,V);
            if((*std::max_element(S.begin(),S.end())) / (*std::min_element(S.begin(),S.end())) > 1e5)
                continue;
            PD = true;


            /// Test SPD
                AB =A;
                LA->is_SPD(AB,ZERO_TOL * S.rows() * (*std::max_element(S.begin(),S.end())));
                std::cout << " [PASS] is_SPD: " << std::endl;

            /// Test SPSD
                AB =A;
                LA->is_SPSD(AB,ZERO_TOL * S.rows() * (*std::max_element(S.begin(),S.end())));
                std::cout << " [PASS] is_SPSD(SPD): " << std::endl;

        }

        // Invert A
        B = A;
        LA->factor_chol(B);

        /// TEST CHOL SOLVE
            MatR x = randM(s,1),b(s,1);
            //Ax1 = b
            A.mult(x,b);
            MatR xb = b;
            //(A^-1)b = x2
            LA->solve_chol_fast(B,xb);
            // x1 == x2 ?
            assert(checkError(x,xb));
            std::cout << " [PASS] solve_chol_fast: " << std::endl;

        B = A;
        LA->factor_chol(B);

        /// TEST CHOL INVERSE
            // B = (A^-1)
            LA->inverse_chol(B);
            // A * B
            A.mult(B,AB);
            // I == A*B ?
            MatR I(s,s);
            I.set_identity();
            assert(checkError(AB,I));
            std::cout << " [PASS] inverse_chol: " << std::endl;

        B = A;

        /// TEST SPD INVERSE
            // B = (A^-1)
            LA->inverse_SPD(B);
            // A * B
            A.mult(B,AB);
            // I == A*B ?
            I.set_identity();
            assert(checkError(AB,I));
            std::cout << " [PASS] inverse_SPD: " << std::endl;

    }
    std::cout << "<< [PASS] testChol: " << std::endl;

}

void testLA(){
    std::cout << ">> testLA: " << std::endl;

    for(int i=1;i<MAX_SIZE;i++){

        unsigned s = std::pow(2,i),c = s;
        MatR A,B,AB(s,s);

        std::cout << "SIZE: " << s << std::endl;

        // Create PD matrix A
        bool PD = false;
        while(!PD){
            A = randM(s,c);
            B = A;
            B.transpose();
            A.mult(B,AB);
            A = AB;
            MatR U,V;
            VecR S;

            /// NOTE: TEST SVD
            LA->svd(AB,U,S,V);

            if((*std::max_element(S.begin(),S.end())) / (*std::min_element(S.begin(),S.end())) < 1e5)
                PD = true;
        }

        B = A;
        MatR x,b;
        b.resize(s,1);

        /// TEST solve_LS_fast1
            x = randM(s,1);
            //Ax1 = b
            A.mult(x,b);
            MatR xb = b;
            //(A^-1)b = x2
            LA->solve_LS_fast1(B,xb);
            // x1 == x2 ?
            assert(checkError(x,xb));
            std::cout << " [PASS] solve_LS_fast1: " << std::endl;

        B = A;

        /// TEST solve_LS_fast2
            x = randM(s,1);
            //Ax1 = b
            A.mult(x,b);
            xb = b;
            //(A^-1)b = x2
            LA->solve_LS_fast2(B,xb);
            std::cout << " [PASS] solve_LS_fast2: " << std::endl;
            // x1 == x2 ?
            assert(checkError(x,xb));

        B = A;

        /// TEST solve_fast
            x = randM(s,1);
            //Ax1 = b
            A.mult(x,b);
            xb = b;
            //(A^-1)b = x2
            LA->solve_fast(B,xb);
            std::cout << " [PASS] solve_fast: " << std::endl;
            // x1 == x2 ?
            assert(checkError(x,xb));

        B = A;

        /// TEST solve_symmetric_fast
            x = randM(s,1);
            //Ax1 = b
            A.mult(x,b);
            xb = b;
            //(A^-1)b = x2
            LA->solve_symmetric_fast(B,xb);
            std::cout << " [PASS] solve_symmetric_fast: " << std::endl;
            // x1 == x2 ?
            assert(checkError(x,xb));

    }
    std::cout << "<< [PASS] testLA: " << std::endl;

}

void TestLinearAlgebra(){
    LA = new LinAlg();
    testSVD();
    testChol();
    testLA();
}


