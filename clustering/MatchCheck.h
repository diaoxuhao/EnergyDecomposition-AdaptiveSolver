#pragma once
#include <iostream>
#include <algorithm>
#include <vector>
#include <Eigen/Core>
#include "SparseEigSolver.h"

void Constructmatrix(Eigen::SparseMatrix<double> &subOED, Eigen::SparseMatrix<double> &supOED, Patch *P1, Patch *P2)
{
	PVedge *temp1;
	Vedge *temp2;
	temp1 = P1->Member;
	for (int i = 0; i< P1->Psize; ++i)
	{
		double sd = 0;
		double so = 0;
		temp2 = temp1->Outvertex->Neighbour;

		for (int j = 0; j<temp1->Outvertex->Degree; ++j)
		{
			if (temp2->Outvertex->Relativecor[0] == P1->Index)
			{
				sd += temp2->diagonal;
				so += std::abs(temp2->offdiagonal);
				subOED.coeffRef(temp1->Outvertex->Relativecor[1], temp2->Outvertex->Relativecor[1]) = temp2->offdiagonal;
				supOED.coeffRef(temp1->Outvertex->Relativecor[1], temp2->Outvertex->Relativecor[1]) = temp2->offdiagonal;
			}
			if (temp2->Outvertex->Relativecor[0] == P2->Index)
			{
				sd += temp2->diagonal;
				so += std::abs(temp2->offdiagonal);
				subOED.coeffRef(temp1->Outvertex->Relativecor[1], temp2->Outvertex->Relativecor[1] + P1->Psize) = temp2->offdiagonal;
				supOED.coeffRef(temp1->Outvertex->Relativecor[1], temp2->Outvertex->Relativecor[1] + P1->Psize) = temp2->offdiagonal;
			}
			temp2 = temp2->Next;
		}
		subOED.coeffRef(temp1->Outvertex->Relativecor[1], temp1->Outvertex->Relativecor[1]) = sd + temp1->Outvertex->Selfloop;
		//supOED.coeffRef(temp1->Outvertex->Relativecor[1], temp1->Outvertex->Relativecor[1]) = 2 * temp1->Outvertex->Wdegree - s + temp1->Outvertex->Selfloop;
		supOED.coeffRef(temp1->Outvertex->Relativecor[1], temp1->Outvertex->Relativecor[1]) = temp1->Outvertex->totaldiagonal + temp1->Outvertex->totaloffdiagonal - so + temp1->Outvertex->Selfloop;



		temp1 = temp1->Next;
	}

	temp1 = P2->Member;
	for (int i = 0; i< P2->Psize; ++i)
	{
		double sd = 0;
		double so = 0;
		temp2 = temp1->Outvertex->Neighbour;
		for (int j = 0; j<temp1->Outvertex->Degree; ++j)
		{
			if (temp2->Outvertex->Relativecor[0] == P1->Index)
			{
				sd += temp2->diagonal;
				so += std::abs(temp2->offdiagonal);
				subOED.coeffRef(temp1->Outvertex->Relativecor[1] + P1->Psize, temp2->Outvertex->Relativecor[1]) = temp2->offdiagonal;
				supOED.coeffRef(temp1->Outvertex->Relativecor[1] + P1->Psize, temp2->Outvertex->Relativecor[1]) = temp2->offdiagonal;
			}
			if (temp2->Outvertex->Relativecor[0] == P2->Index)
			{
				sd += temp2->diagonal;
				so += std::abs(temp2->offdiagonal);
				subOED.coeffRef(temp1->Outvertex->Relativecor[1] + P1->Psize, temp2->Outvertex->Relativecor[1] + P1->Psize) = temp2->offdiagonal;
				supOED.coeffRef(temp1->Outvertex->Relativecor[1] + P1->Psize, temp2->Outvertex->Relativecor[1] + P1->Psize) = temp2->offdiagonal;
			}
			temp2 = temp2->Next;
		}
		subOED.coeffRef(temp1->Outvertex->Relativecor[1] + P1->Psize, temp1->Outvertex->Relativecor[1] + P1->Psize) = sd + temp1->Outvertex->Selfloop;
		//supOED.coeffRef(temp1->Outvertex->Relativecor[1] + P1->Psize, temp1->Outvertex->Relativecor[1] + P1->Psize) = 2 * temp1->Outvertex->Wdegree - s + temp1->Outvertex->Selfloop;
		supOED.coeffRef(temp1->Outvertex->Relativecor[1] + P1->Psize, temp1->Outvertex->Relativecor[1] + P1->Psize) = temp1->Outvertex->totaldiagonal + temp1->Outvertex->totaloffdiagonal - so + temp1->Outvertex->Selfloop;

		temp1 = temp1->Next;
	}

}

bool Checkmatrix(Eigen::SparseMatrix<double> &subOED, Eigen::SparseMatrix<double> &supOED, Eigen::VectorXd *eigfunc, 
	double *delta, double *lambda, double epsilon, double ratio, int maxPsizeControl)
{

	int npart = subOED.cols();

	if (npart < maxPsizeControl)
	{
		Eigen::VectorXd evalue;
		Eigen::VectorXd temp;
		Eigen::VectorXd eigenvalue;
		Eigen::MatrixXd eigenvector;

		eigs(subOED, 2, "sm", eigenvalue, eigenvector);
		*lambda = eigenvalue(1);
		temp = eigenvector.col(0);
		*eigfunc = temp;

		Eigen::VectorXd x(npart);
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
		solver.compute(supOED);
		x = solver.solve(temp);
		*delta = temp.dot(temp) / x.dot(temp);
		return (*lambda * epsilon >= 1 && *delta / *lambda < ratio);
	}
	else
	{
		return 0;
	}
}

bool Check(Patch *P1, Patch *P2, double *delta, double *lambda, double epsilon, double ratio, int maxPsizeControl)
{
	Eigen::SparseMatrix<double> subOED(P1->Psize + P2->Psize, P1->Psize + P2->Psize);
	Eigen::SparseMatrix<double> supOED(P1->Psize + P2->Psize, P1->Psize + P2->Psize);
	Eigen::VectorXd eigfunction;
	Constructmatrix(subOED, supOED, P1, P2);

	int value = Checkmatrix(subOED, supOED, &eigfunction, delta, lambda, epsilon, ratio, maxPsizeControl);
	PVedge *temp1;

	if (value)
	{
		temp1 = P1->Member;
		for (int i = 0; i < P1->Psize; ++i)
		{
			temp1->Eigfun = eigfunction(i);
			temp1 = temp1->Next;
		}
		temp1 = P2->Member;
		for (int i = 0; i < P2->Psize; ++i)
		{
			temp1->Eigfun = eigfunction(i + P1->Psize);
			temp1 = temp1->Next;
		}
	}

	return value;
}