#ifndef SOLVERS_GPU_H
#define SOLVERS_GPU_H

#include<bits/stdc++.h>
#include <Eigen/Dense>
#include <cuda_runtime.h>

#include <cublas_v2.h>
#include <cusolverDn.h>

#include <complex>
#include <cuComplex.h>
#include <chrono>
using namespace std;


void cusolverComplex(   Eigen::MatrixXcd &A1,
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
    //cout<<"cusolver status = "<<cusolver_status<<endl;
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



void cusolverComplex_mat(   Eigen::MatrixXcd &A1,
                            Eigen::MatrixXcd &B1,
                            Eigen::MatrixXcd &X1,
                            int mat_length,
                            int mat_width          )
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

    
    const int m = mat_length;
    const int lda = m;
    const int ldb = m;
    const int nrhs = mat_width; // number of right hand side vectors
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
    for(int i=0;i<nrhs;i++)
    {
        for(int j=0;j<ldb;j++)
        {
            B[somenum] = make_cuDoubleComplex(B1(j,i).real(),B1(j,i).imag());
            somenum = somenum + 1; 
            //std::cout<<"somenum = "<<somenum<<'\t'<<" i,j = "<<i<<'\t'<<j<<std::endl; 
        }
    }
    //std::cout<<"ami ekhane"<<std::endl; 

    somenum = 0; 

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
    //cout<<"cusolver status = "<<cusolver_status<<endl;
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
    somenum = 0;
    for(int i=0;i<nrhs;i++) 
    {
        for(int j=0;j<ldb;j++)
        {
            X1(j,i) = cuCreal(XC[somenum]) + ii*cuCimag(XC[somenum]);
            somenum = somenum + 1; 
        }
    }
    //std::cout<<"ami ekhane asi"<<std::endl;
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






void cusolverComplexAsync(  Eigen::MatrixXcd &A1,
                            Eigen::VectorXcd &B1,
                            Eigen::VectorXcd &X1,
                            int matsize,
                            cudaStream_t mystream          )
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
    cudaStat1 = cudaMemcpyAsync(d_A, A, sizeof(cuDoubleComplex) * lda * m   , cudaMemcpyHostToDevice, mystream);
    cudaStat2 = cudaMemcpyAsync(d_B, B, sizeof(cuDoubleComplex) * ldb * nrhs, cudaMemcpyHostToDevice, mystream);
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
    //cudaStat1 = cudaDeviceSynchronize();
    assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);
    //assert(cudaSuccess == cudaStat1);

    // check if QR is good or not
    cudaStat1 = cudaMemcpyAsync(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost, mystream);
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
    //cudaStat1 = cudaDeviceSynchronize();
    assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);
    //assert(cudaSuccess == cudaStat1);

    
 
    // check if QR is good or not
    cudaStat1 = cudaMemcpyAsync(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost, mystream);
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
    //cudaStat1 = cudaDeviceSynchronize();
    assert(CUBLAS_STATUS_SUCCESS == cublas_status);
    //assert(cudaSuccess == cudaStat1);

    cudaStat1 = cudaMemcpyAsync(XC, d_B, sizeof(cuDoubleComplex)*ldb*nrhs, cudaMemcpyDeviceToHost, mystream);
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

    //cudaDeviceReset();
    
    delete[] A;
    delete[] B;
    delete[] XC;

    //eigencheck();

    //return 0;
}




#endif
