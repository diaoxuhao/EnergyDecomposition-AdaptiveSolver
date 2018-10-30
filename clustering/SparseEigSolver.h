#pragma once
#include <iostream>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/Eigenvalues>
#include <Eigen/SparseCholesky>

#include <SparseSymShiftSolve.h>
#include <SymEigsShiftSolver.h>
#include <SymEigsSolver.h>
#include <SparseGenMatProd.h>

void eigs(Eigen::SparseMatrix<double> &A, int num, std::string sigma, Eigen::VectorXd &eigenvalue, Eigen::MatrixXd &eigenvector)
{
	int nrow = A.rows();
	int ncol = A.cols();
	if (nrow != ncol)
	{
		std::cout << "Error: The input matrix must be Symmetric!" << std::endl;
		abort();
	}

	eigenvalue.resize(ncol, 1);
	if ((double)num > nrow)
	{
		std::cout << "Error: Number of eigenvalues requested k > n!" << std::endl;
		abort();
	}


	if (nrow <= 6 || ncol <= 6)
	{
		Eigen::SelfAdjointEigenSolver< Eigen::SparseMatrix<double> > eigensolver(A);
		if (eigensolver.info() == Eigen::Success)
		{
			eigenvalue = eigensolver.eigenvalues();
			eigenvector = eigensolver.eigenvectors();
		}
		else
		{
			std::cout << "Error: Problems finding eienvalues!" << std::endl;
			abort();
		}
	}
	else
	{
		if (sigma == "sm")
		{
			Spectra::SparseSymShiftSolve<double> op(A);
			Spectra::SymEigsShiftSolver< double, Spectra::LARGEST_MAGN, Spectra::SparseSymShiftSolve<double> > eigensolver(&op, num, 2 * num, -0.001);
			eigensolver.init();

			int nconv = eigensolver.compute();
			if (eigensolver.info() == Spectra::SUCCESSFUL)
			{
				eigenvalue = (eigensolver.eigenvalues()).reverse();
				eigenvector = (eigensolver.eigenvectors()).rowwise().reverse();
			}
			else
			{
				std::cout << "Error: Problems finding eienvalues!" << std::endl;
				abort();
			}
		}
		else if (sigma == "lm")
		{
			Spectra::SparseGenMatProd<double> op(A);
			Spectra::SymEigsSolver< double, Spectra::LARGEST_MAGN, Spectra::SparseGenMatProd<double> > eigensolver(&op, num, 2 * num);
			eigensolver.init();

			int nconv = eigensolver.compute();
			if (eigensolver.info() == Spectra::SUCCESSFUL)
			{
				eigenvalue = (eigensolver.eigenvalues()).reverse();
				eigenvector = (eigensolver.eigenvectors()).rowwise().reverse();
			}
			else
			{
				std::cout << "Error: Problems finding eienvalues!" << std::endl;
				abort();
			}
		}
		else
		{
			std::cout << "Error: Wrong sigma input. Possible choices are 'lm' and 'sm'." << std::endl;
			abort();
		}
	}

}


double norm2(Eigen::SparseMatrix<double> &A)
{
	int nrow = A.rows();
	int ncol = A.cols();
	if (nrow != ncol)
	{
		std::cout << "Error: The input matrix must be Symmetric!" << std::endl;
		abort();
	}

	if (nrow <= 6 || ncol <= 6)
	{
		Eigen::VectorXd eigenvalue(ncol);
		Eigen::SelfAdjointEigenSolver< Eigen::SparseMatrix<double> > eigensolver(A);
		if (eigensolver.info() == Eigen::Success)
		{
			eigenvalue = eigensolver.eigenvalues();
			return eigenvalue(eigenvalue.size() - 1);
		}
		else
		{
			std::cout << "Error: Problems finding eienvalues!" << std::endl;
			abort();
		}
	}
	else
	{
		Eigen::VectorXd eigenvalue(1);
		Spectra::SparseGenMatProd<double> op(A);
		Spectra::SymEigsSolver< double, Spectra::LARGEST_MAGN, Spectra::SparseGenMatProd<double> > eigensolver(&op, 1, 4);
		eigensolver.init();

		int nconv = eigensolver.compute();
		if (eigensolver.info() == Spectra::SUCCESSFUL)
		{
			eigenvalue = eigensolver.eigenvalues();
			return eigenvalue(0);
		}
		else
		{
			std::cout << "Error: Problems finding eienvalues!" << std::endl;
			abort();
		}
	}
}


double cond(Eigen::SparseMatrix<double> &A)
{
	int nrow = A.rows();
	int ncol = A.cols();
	if (nrow != ncol)
	{
		std::cout << "Error: The input matrix must be Symmetric!" << std::endl;
		abort();
	}

	if (nrow <= 6 || ncol <= 6)
	{
		Eigen::VectorXd eigenvalue(ncol);
		Eigen::SelfAdjointEigenSolver< Eigen::SparseMatrix<double> > eigensolver(A);
		if (eigensolver.info() == Eigen::Success)
		{
			eigenvalue = eigensolver.eigenvalues();
			if (std::abs(eigenvalue(0)) < 1e-10)
			{
				std::cout << "Warning:Smallest eigenvalue is close to zero." << std::endl;
			}
			return eigenvalue(eigenvalue.size() - 1) / eigenvalue(0);
		}
		else
		{
			std::cout << "Error: Problems finding eienvalues!" << std::endl;
			abort();
		}
	}
	else
	{
		Eigen::VectorXd eigenvalue_sm(1);
		Eigen::VectorXd eigenvalue_lm(1);

		Spectra::SparseSymShiftSolve<double> opsm(A);
		Spectra::SymEigsShiftSolver< double, Spectra::LARGEST_MAGN, Spectra::SparseSymShiftSolve<double> > eigensolver_sm(&opsm, 1, 4, -0.001);
		eigensolver_sm.init();

		int nconv_sm = eigensolver_sm.compute();
		if (eigensolver_sm.info() == Spectra::SUCCESSFUL)
		{
			eigenvalue_sm = eigensolver_sm.eigenvalues();
		}
		else
		{
			std::cout << "Error: Problems finding eienvalues!" << std::endl;
			abort();
		}
		Spectra::SparseGenMatProd<double> oplm(A);
		Spectra::SymEigsSolver< double, Spectra::LARGEST_MAGN, Spectra::SparseGenMatProd<double> > eigensolver_lm(&oplm, 1, 4);
		eigensolver_lm.init();

		int nconv_lm = eigensolver_lm.compute();
		if (eigensolver_lm.info() == Spectra::SUCCESSFUL)
		{
			eigenvalue_lm = eigensolver_lm.eigenvalues();
		}
		else
		{
			std::cout << "Error: Problems finding eienvalues!" << std::endl;
			abort();
		}
		if (std::abs(eigenvalue_sm(0)) < 1e-10)
		{
			std::cout << "Warning:Smallest eigenvalue is close to zero." << std::endl;
		}
		return eigenvalue_lm(0) / eigenvalue_sm(0);
	}


}