#ifndef SOLVERS_H
#define SOLVERS_H

#include<bits/stdc++.h>
#include <Eigen/Dense>
//#include <cuda_runtime.h>

//#include <cublas_v2.h>
//#include <cusolverDn.h>

#include <complex>
//#include <cuComplex.h>
#include <chrono>
using namespace std;


void LinearSolver(	Eigen::MatrixXcd &A,
					Eigen::VectorXcd &X,
					Eigen::VectorXcd &B,
					double &relerr 			)
{
	X = A.householderQr().solve(B);

	relerr = (A*X - B).norm()/B.norm();
}

void LinearSolver_1(	Eigen::MatrixXcd &A,
					Eigen::VectorXcd &X,
					Eigen::VectorXcd &B,
					double &relerr 			)
{
	//X = A.BDCSVD().solve(B);

    cout<<"Determinant of the matrix A = "<<A.determinant()<<endl;
    Eigen::MatrixXcd ker = A.fullPivLu().kernel();

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver(A);
    if (eigensolver.info() != Eigen::Success) abort();
    std::cout << "The eigenvalues of A are:\n" <<setprecision(15)<< eigensolver.eigenvalues() << std::endl;
    cout<<"kernel of the matrix A ="<<endl;
    cout<<ker<<endl;
    std::cout << "Here's a matrix whose columns are eigenvectors of A \n"
        << "corresponding to these eigenvalues:\n"
        << eigensolver.eigenvectors() << std::endl;


	//relerr = (A*X - B).norm()/B.norm();
}

void LinearSolver_2(	Eigen::MatrixXcd &A,
					Eigen::VectorXcd &X,
					Eigen::VectorXcd &B,
					double &relerr 			)
{
    //this makes it run on a cpu 0 always
    //cpu_set_t set;
    //CPU_ZERO(&set);        // clear cpu mask
    //CPU_SET(0, &set);      // set cpu 0
    //sched_setaffinity(0, sizeof(cpu_set_t), &set);  // 0 is the calling process
	//----------------------------------------------//

    
    X = A.partialPivLu().solve(B);

	relerr = (A*X - B).norm()/B.norm();
}
/*void LinearSolver_2(	Eigen::MatrixXcd &A,
					Eigen::VectorXcd &X,
					Eigen::VectorXcd &B,
					double &relerr 			)
{
	
    Eigen::JacobiSVD<MatrixXf, ComputeThinU | ComputeThinV> svd(m);
    X = svd.solve(B);

	relerr = (A*X - B).norm()/B.norm();
}
*/
/*void cusolverComplex(   Eigen::MatrixXcd &A1,
                        Eigen::VectorXcd &B1,
                        Eigen::VectorXcd &X1,
                        int matsize          )
{
   //printf("lalala\n");
    cusolverDnHandle_t cusolverH = NULL;
    cublasHandle_t cublasH = NULL;
    cublasStatus_t cublas_status = CUBLAS_STATUS_SUCCESS;
    cusolverStatus_t cusolver_status = CUSOLVER_STATUS_SUCCESS;    
    cudaError_t cudaStat1 = cudaSuccess;
    cudaError_t cudaStat2 = cudaSuccess;
    cudaError_t cudaStat3 = cudaSuccess;
    cudaError_t cudaStat4 = cudaSuccess;
    const int m = matsize;
    const int lda = m;
    const int ldb = m;
    const int nrhs = 1; // number of right hand side vectors
//       | 1 2 3 |
//   A = | 4 5 6 |
//       | 2 1 1 |
//
//   x = (1 1 1)'
//   b = (6 15 4)'
 //
      //printf("line250\n");
    int somenum = 0;
    cuDoubleComplex* A = new cuDoubleComplex[lda*m];// = { 1.0, 4.0, 2.0, 2.0, 5.0, 1.0, 3.0, 6.0, 1.0}; 
//printf("line254\n");
    for(int i=0;i<m;i++)
    {
       for(int j=0;j<m;j++)
       {
          A[somenum] = make_cuDoubleComplex(A1(j,i).real(),A1(j,i).imag());
          somenum = somenum + 1;
       }
    }
    somenum = 0;
//    double X[ldb*nrhs] = { 1.0, 1.0, 1.0}; // exact solution
    cuDoubleComplex* B = new cuDoubleComplex[ldb*nrhs];// = { 6.0, 15.0, 4.0}; 
    for(int i=0;i<m;i++)
    {
       B[i] = make_cuDoubleComplex(B1(i).real(),B1(i).imag());
    }
    cuDoubleComplex* XC = new cuDoubleComplex[ldb*nrhs]; // solution matrix from GPU

    cuDoubleComplex *d_A = NULL; // linear memory of GPU  
    cuDoubleComplex *d_tau = NULL; // linear memory of GPU 
    cuDoubleComplex *d_B  = NULL; 
    int *devInfo = NULL; // info in gpu (device copy)
    cuDoubleComplex *d_work = NULL;
    int  lwork = 0; 
//printf("line277\n");
    int info_gpu = 0;

    const cuDoubleComplex one = make_cuDoubleComplex(1.0,0.0);
//    const double one = 1;

    //printf("A = (matlab base-1)\n");
    //printMatrixComplex(m, m, A, lda, "A");
    //printf("=====\n");
    //printf("B = (matlab base-1)\n");
    //printMatrixComplex(m, nrhs, B, ldb, "B");
    //printf("=====\n");
   //printf("line288");
// step 1: create cusolver/cublas handle
    cusolver_status = cusolverDnCreate(&cusolverH);
    assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);

    cublas_status = cublasCreate(&cublasH);
    assert(CUBLAS_STATUS_SUCCESS == cublas_status);
    
// step 2: copy A and B to device
    cudaStat1 = cudaMalloc ((void**)&d_A  , sizeof(cuDoubleComplex) * lda * m);
    cudaStat2 = cudaMalloc ((void**)&d_tau, sizeof(cuDoubleComplex) * m);
    cudaStat3 = cudaMalloc ((void**)&d_B  , sizeof(cuDoubleComplex) * ldb * nrhs);
    cudaStat4 = cudaMalloc ((void**)&devInfo, sizeof(int));
    assert(cudaSuccess == cudaStat1);
    assert(cudaSuccess == cudaStat2);
    assert(cudaSuccess == cudaStat3);
    assert(cudaSuccess == cudaStat4);
   //printf("line305");
    cudaStat1 = cudaMemcpy(d_A, A, sizeof(cuDoubleComplex) * lda * m   , cudaMemcpyHostToDevice);
    cudaStat2 = cudaMemcpy(d_B, B, sizeof(cuDoubleComplex) * ldb * nrhs, cudaMemcpyHostToDevice);
    assert(cudaSuccess == cudaStat1);
    assert(cudaSuccess == cudaStat2);


 
// step 3: query working space of geqrf and ormqr
    cusolver_status = cusolverDnZgeqrf_bufferSize(
        cusolverH, 
        m, 
        m, 
        d_A, 
        lda, 
        &lwork);
    assert (cusolver_status == CUSOLVER_STATUS_SUCCESS);
 
    //cudaStat1 = cudaMalloc((void**)&d_work, sizeof(double)*lwork);
    cudaStat1 = cudaMalloc((void**)&d_work, sizeof(cuDoubleComplex)*lwork);
    assert(cudaSuccess == cudaStat1);

// step 4: compute QR factorization
    cusolver_status = cusolverDnZgeqrf(
        cusolverH, 
        m, 
        m, 
        d_A, 
        lda, 
        d_tau, 
        d_work, 
        lwork, 
        devInfo);
    cudaStat1 = cudaDeviceSynchronize();
    assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);
    assert(cudaSuccess == cudaStat1);

    // check if QR is good or not
    cudaStat1 = cudaMemcpy(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost);
    assert(cudaSuccess == cudaStat1);

    //printf("after geqrf: info_gpu = %d\n", info_gpu);
    assert(0 == info_gpu);

// step 5: compute Q^T*B
    //cusolver_status = cusolverDnZunmqr_bufferSize(
    //cusolverH,
    //CUBLAS_SIDE_LEFT, 
    //    CUBLAS_OP_C,
    //m,
    //nrhs,
    //m,
    //d_A,
    //lda,
    //d_tau,
    //d_B,
    //ldb,
    //lwork);
   //cudaStat1 = cudaDeviceSynchronize();
    //assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);
    
     cusolver_status= cusolverDnZunmqr(
        cusolverH, 
        CUBLAS_SIDE_LEFT, 
        CUBLAS_OP_C,
        m, 
        nrhs, 
        m, 
        d_A, 
        lda,
        d_tau,
        d_B,
        ldb,
        d_work,
        lwork,
        devInfo);
    cudaStat1 = cudaDeviceSynchronize();
    assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);
    assert(cudaSuccess == cudaStat1);

    
 
    // check if QR is good or not
    cudaStat1 = cudaMemcpy(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost);
    assert(cudaSuccess == cudaStat1);

    //printf("after ormqr: info_gpu = %d\n", info_gpu);
    assert(0 == info_gpu);

// step 6: compute x = R \ Q^T*B

    cublas_status = cublasZtrsm(
         cublasH,
         CUBLAS_SIDE_LEFT,
         CUBLAS_FILL_MODE_UPPER,
         CUBLAS_OP_N, 
         CUBLAS_DIAG_NON_UNIT,
         m,
         nrhs,
         &one,
         d_A,
         lda,
         d_B,
         ldb);
    cudaStat1 = cudaDeviceSynchronize();
    assert(CUBLAS_STATUS_SUCCESS == cublas_status);
    assert(cudaSuccess == cudaStat1);

    cudaStat1 = cudaMemcpy(XC, d_B, sizeof(cuDoubleComplex)*ldb*nrhs, cudaMemcpyDeviceToHost);
    assert(cudaSuccess == cudaStat1);

    std::complex<double> ii = {0.0,1.0};
    for(int i=0;i<m;i++) X1(i) = cuCreal(XC[i]) + ii*cuCimag(XC[i]);
    //printf("X = (matlab base-1)\n");
    //printMatrixComplex(m, nrhs, XC, ldb, "X");

   

// free resources
    if (d_A    ) cudaFree(d_A);
    if (d_tau  ) cudaFree(d_tau);
    if (d_B    ) cudaFree(d_B);
    if (devInfo) cudaFree(devInfo);
    if (d_work ) cudaFree(d_work);


    if (cublasH ) cublasDestroy(cublasH);   
    if (cusolverH) cusolverDnDestroy(cusolverH);   

    cudaDeviceReset();
    
    delete[] A;
    delete[] B;
    delete[] XC;

    //eigencheck();

    //return 0;
}



*/
#endif
