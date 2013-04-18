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

void testSVD(){
    std::cout << ">> testSVD: " << std::endl;

    for(int i=1;i<MAX_SIZE;i++){
        unsigned r = 1 << (i-1),c = r;
        std::cout << "SIZE: " << r << std::endl;
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
        checkError(std::cout, "svd1", A,AB);

        AB = A;
        /// Test the SVD2 Decomp
        LA->svd2(AB,U,s,V);
        S.set_zero();
        for(unsigned i=0;i<s.rows();i++)
            S(i,i) = s[i];

        U.mult(S,US);
        US.mult_transpose(V,AB);
        checkError(std::cout, "svd2", A,AB);

        /// Test SVD Solve
        AB = A;
//        LA->svd2(AB,U,s,V);

        MatR x = randM(r,1),b(r,1);
        //Ax1 = b
        A.mult(x,b);
        MatR xb = b;
        //(A^-1)b = x2
        LA->solve_LS_fast(U,s,V,xb);
        // x1 == x2 ?
        checkError(std::cout, "solve_LS_fast(U,S,V)", x,xb);
    }
}

void testLU(){
    std::cout << ">> testLU: " << std::endl;

    for(int i=1;i<MAX_SIZE;i++){
        unsigned s = 1 << (i-1),c = 2;
        MatR A,B,AB(s,s);

        // Create PD matrix A
        A = randM(s,s);
        B = A;
        B.transpose();
        A.mult(B,AB);
        A = AB;
        std::vector<int> P;

        MatR LU = A;
        LA->factor_LU(LU,P);

        MatR L(s,s), U(s,s);
        L.set_zero();
        U.set_zero();
        /// Test the SVD Decomp
        for(unsigned ii=0;ii<s;ii++){
            L(ii,ii) = 1;
            for(unsigned jj=0;jj<s;jj++){
                if(jj<ii){
                    L(ii,jj) = LU(ii,jj);
                } else {
                    U(ii,jj) = LU(ii,jj);
                }
            }
        }
        MatR LU2(s,s);
        L.mult(U,LU2);

        // Permutation is confusing... needs to be done back to front
        MatR A2(s,s);
        A2 = LU2;
        for(int ii = A2.rows()-1;ii>=0;ii--){
            VecR row1(A2.columns()), row2(A2.columns());
            A2.get_row(ii,row1);
            A2.get_row(P[ii]-1,row2);

            A2.set_row(ii,row2);
            A2.set_row(P[ii]-1,row1);
        }

        checkError(std::cout, "factor_LU", A,A2);




        /// Test SVD Solve
        MatR x = randM(s,1),b(s,1);
        //Ax1 = b
        A.mult(x,b);
        MatR xb = b;
//        for(int ii = A2.rows()-1;ii>=0;ii--)
//            std::cout << P[ii] << std::endl;
//            std::swap(xb(ii,0),xb(P[ii]-1,0));
//        for(int ii = A2.rows()-1;ii>=0;ii--)
//            std::swap(x(ii,0),x(P[ii]-1,0));
        //(A^-1)b = x2

        LA->solve_LU_fast(LU,false,P,xb);

        // x1 == x2 ?
        checkError(std::cout, "solve_LU_fast", x,xb);

    }
    std::cout << "<< [PASS] testLU: " << std::endl;
}

void testChol(){
    std::cout << ">> testChol: " << std::endl;

    for(int i=1;i<MAX_SIZE;i++){
        unsigned s = 1 << (i-1), c = s;
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
                checkError(std::cout, "eig_symm", evals,S);

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
            checkError(std::cout, "solve_chol_fast", x,xb);

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
            checkError(std::cout, "inverse_chol", AB,I);

        B = A;

        /// TEST SPD INVERSE
            // B = (A^-1)
            LA->inverse_SPD(B);
            // A * B
            A.mult(B,AB);
            // I == A*B ?
            I.set_identity();
            checkError(std::cout, "inverse_SPD", AB,I);

    }

}

void testLA(){
    std::cout << ">> testLA: " << std::endl;

    for(int i=1;i<MAX_SIZE;i++){

        unsigned s = 1 << (i-1), c = s;
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
            checkError(std::cout, "solve_LS_fast1", x,xb);

        B = A;

        /// TEST solve_LS_fast2
            x = randM(s,1);
            //Ax1 = b
            A.mult(x,b);
            xb = b;
            //(A^-1)b = x2
            LA->solve_LS_fast2(B,xb);
            // x1 == x2 ?
            checkError(std::cout, "solve_LS_fast2", x,xb);

        B = A;

        /// TEST solve_fast
            x = randM(s,1);
            //Ax1 = b
            A.mult(x,b);
            xb = b;
            //(A^-1)b = x2
            LA->solve_fast(B,xb);
            // x1 == x2 ?
            checkError(std::cout, "solve_fast", x,xb);

        B = A;

        /// TEST solve_symmetric_fast
            x = randM(s,1);
            //Ax1 = b
            A.mult(x,b);
            xb = b;
            //(A^-1)b = x2
            LA->solve_symmetric_fast(B,xb);
            // x1 == x2 ?
            checkError(std::cout, "solve_symmetric_fast", x,xb);

    }

}

void TestLinearAlgebra(){
    LA = new LinAlg();
    testSVD();
    testChol();
    testLA();
    testLU();
}


