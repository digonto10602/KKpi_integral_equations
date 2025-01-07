#ifndef LINEAR_SOLVERS_H
#define LINEAR_SOLVERS_H

#include<bits/stdc++.h>
#include <Eigen/Dense>
#include <complex>
#include <chrono>


//These solvers are based on matrices rather than vectors like 
//we did previously
void LinearSolver(	Eigen::MatrixXcd &A,
					Eigen::MatrixXcd &X,
					Eigen::MatrixXcd &B,
					double &relerr 			)
{
	X = A.householderQr().solve(B);

	relerr = (A*X - B).norm()/B.norm();
}

void LinearSolver_2(	Eigen::MatrixXcd &A,
					    Eigen::MatrixXcd &X,
					    Eigen::MatrixXcd &B,
					    double &relerr 			)
{
	X = A.partialPivLu().solve(B);

	relerr = (A*X - B).norm()/B.norm();
}

void LinearSolver_3(	Eigen::MatrixXcd &A, 
						Eigen::VectorXcd &X,
						Eigen::VectorXcd &B,
						double &relerr			)
{
	X = A.partialPivLu().solve(B);

	relerr = (A*X - B).norm()/B.norm();
}

#endif