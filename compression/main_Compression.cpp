#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCore>
#include <iomanip>
#include "DataStructure_Compression.h"
#include "LocalBasis.h"
#include "ReduceGraph.h"


// arguments : 1. matrix A,   2. graph,   3. partition,
//             4. neighbour,  5. epsilon,  6. locPsi,
//             7. matrix Ast,  8. reduced graph 			  	

int main(int argv, char* argc[])
{
	typedef Eigen::Triplet<double> triplet;

	std::fstream Matrix_A, Partition, Patch_Neighbour;
	int n, m, M, new_m;
	int i, j, k=0, s, r;
	int temp_index;
	double skip;
	double epsilon;
	int Max_Psize=0, Max_Degree=0;
	double totaltime = 0;

	std::vector<triplet> tp_A;

	Matrix_A.open(argc[1]);
	Partition.open(argc[3]);
	Patch_Neighbour.open(argc[4]);
	epsilon=atof(argc[5]);

	Matrix_A >> n;
	Matrix_A >> m;
	Partition >> M;
	Patch_Neighbour >> skip;

	std::vector<Vertex> Vertex_list(n);
	std::vector<Vedge> Vedge_list(m);
	std::vector<Patch> Patch_list(M);
	std::vector<PVedge> PVedge_list(m);
	std::vector<Pedge> Pedge_list(m);

	for (i = 0; i < n; ++i)
	{
		Matrix_A >> Vertex_list[i].Index;
		Matrix_A >> Vertex_list[i].Diagonal;
		Matrix_A >> Vertex_list[i].Degree; 
		tp_A.push_back(triplet(i, i, Vertex_list[i].Diagonal));
		Vertex_list[i].Neighbour = &Vedge_list[k]; 
		for (j=0; j<Vertex_list[i].Degree; ++j)
		{
			Matrix_A >> Vedge_list[k].Index;
			Matrix_A >> Vedge_list[k].Offdiag;
			tp_A.push_back(triplet(i, Vedge_list[k].Index, Vedge_list[k].Offdiag));
			if (j<Vertex_list[i].Degree-1)
			{
				Vedge_list[k].Next = &Vedge_list[k+1];
			}
			else
			{
				Vedge_list[k].Next = NULL;
			}
			++k;
		}
	}

	int k1 = 0, k2 = 0;
	new_m = 0;
	for (i = 0; i < M; ++i)
	{
		Partition >> Patch_list[i].Index;
		Partition >> skip;
		Partition >> skip;
		Partition >> Patch_list[i].Psize;
		if (Patch_list[i].Psize>Max_Psize)
		{
			Max_Psize=Patch_list[i].Psize;
		}
		Patch_list[i].Mark=0;
		Patch_list[i].Member=&PVedge_list[k1];
		for (j=0; j<Patch_list[i].Psize; ++j)
		{
			Partition >> PVedge_list[k1].Index;
			Vertex_list[PVedge_list[k1].Index].Pindex=i;
			Partition >> Vertex_list[PVedge_list[k1].Index].Phi;
			Partition >> Vertex_list[PVedge_list[k1].Index].Hous;
			if (j<Patch_list[i].Psize-1)
			{
				PVedge_list[k1].Next = &PVedge_list[k1+1];
			}
			else
			{
				PVedge_list[k1].Next = NULL;
			}
			++k1;
		}

		Patch_Neighbour >> skip;
		Patch_Neighbour >> Patch_list[i].Degree;
		if (Patch_list[i].Degree>Max_Degree)
		{
			Max_Degree=Patch_list[i].Degree;
		}
		new_m += Patch_list[i].Degree;
		Patch_list[i].Neighbour=&Pedge_list[k2];
		for (j=0; j<Patch_list[i].Degree; ++j)
		{
			Patch_Neighbour >> Pedge_list[k2].Index;
			if (j<Patch_list[i].Degree-1)
			{
				Pedge_list[k2].Next = &Pedge_list[k2+1];
			}
			else
			{
				Pedge_list[k2].Next = NULL;
			}
			++k2;
		}
	}

	Matrix_A.close();
	Partition.close();
	Patch_Neighbour.close();

	std::clock_t start;
	double duration;
	//----------------Compute locPsi---------------------------- 
	s=500*Max_Psize;

	std::vector<int> supsize1(M);

	start = std::clock();
	std::vector<triplet> tp_Psi;
	LocalPsi(&Vertex_list, &Patch_list, &tp_Psi, &supsize1, epsilon, M, s);
	std::cout << tp_Psi.size() << std::endl;
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	totaltime += duration;
	std::cout << "Compute locPsi: " << duration << "s." << std::endl;

	Eigen::SparseMatrix<double> locPsi(n, M);
	locPsi.setFromTriplets(tp_Psi.begin(), tp_Psi.end());
	
	start = std::clock();
	std::ofstream localPsi_file;
	localPsi_file.open(argc[6]);
	localPsi_file << M << "\n";
	for (int i = 0; i < locPsi.outerSize(); ++i)
	{
		localPsi_file << i << " " << supsize1[i] << "\n";
		for (Eigen::SparseMatrix<double>::InnerIterator it(locPsi, i); it; ++it)
		{
			localPsi_file << it.row() << " " << std::setprecision(16) << it.value() << " ";
		}
		localPsi_file << "\n";
	}
	localPsi_file.close();
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "Store locPsi: " << duration << "s." << std::endl;

	//------------Compute Ast------------------------ 

	Eigen::SparseMatrix<double> A(n,n);
	Eigen::SparseMatrix<double> Ast(M,M);
	std::vector<triplet> tp_Ast;

	A.setFromTriplets(tp_A.begin(), tp_A.end());
	start = std::clock();	
	Ast=locPsi.transpose()*(A*locPsi);
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	totaltime += duration;
	std::cout << "Compute Ast: " << duration << "s." << std::endl;

	for (i = 0; i < M; ++i)
	{
		for (Eigen::SparseMatrix<double>::InnerIterator it(Ast, i); it; ++it)
		{
			if ( it.row()==i )
			{
				tp_Ast.push_back(triplet(it.row(), i, it.value()));
			}
			if ( it.row()>i && std::abs(it.value())>epsilon*0.1 )
			{
				tp_Ast.push_back(triplet(it.row(), i, it.value()));
				tp_Ast.push_back(triplet(i, it.row(), it.value()));
			}	
		}
	}
	Ast.setFromTriplets(tp_Ast.begin(), tp_Ast.end());

	start = std::clock();
	std::ofstream Ast_file;
	Ast_file.open(argc[7]);
	Ast_file << M << " " << Ast.nonZeros()-M << "\n";
	for (i=0; i<M; ++i)
	{
		Ast_file << i << " ";
		Ast_file << std::setprecision(16) << Ast.coeff(i,i) << " ";
		Ast_file << Ast.innerVector(i).nonZeros()-1 << "\n";
		for (Eigen::SparseMatrix<double>::InnerIterator it(Ast, i); it; ++it)
		{
			if ( it.row() != i )
			{
				Ast_file << it.row() << " " << std::setprecision(16) << it.value() << " ";
			}	
		}
		Ast_file << "\n";
	}
	Ast_file.close();
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "Store Ast: " << duration << "s." << std::endl;

	//------------Reduce Graph-----------------------

	start = std::clock();
	LocalPsi_0(&Vertex_list, &Patch_list, epsilon/10000, M, Max_Psize);
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	totaltime += duration;
	std::cout << "Compute locPsi_0: " << duration << "s." << std::endl;

	std::fstream Graph;
	std::ofstream Re_Graph;
	
	Pedge *tmpP;
	PVedge *tmpPV;

	Graph.open(argc[2]);

	Graph >> n;
	Graph >> m;

	Vedge_list.resize(m);

	k=0;
	for (i=0; i<n; ++i)
	{
		Graph >> skip;
		Graph >> Vertex_list[i].Selfloop;
		Graph >> Vertex_list[i].Degree;
		Vertex_list[i].Neighbour=&Vedge_list[k];
		for (j=0; j<Vertex_list[i].Degree; ++j)
		{
			Graph >> Vedge_list[k].Index;
			Graph >> Vedge_list[k].Diagonal;
			Graph >> Vedge_list[k].Offdiag;
			if (j<Vertex_list[i].Degree-1)
			{
				Vedge_list[k].Next = &Vedge_list[k+1];
			}
			else
			{
				Vedge_list[k].Next = NULL;
			}
			++k;
		}
	}

	Graph.close();

	start = std::clock();
	ReduceGraph(&Vertex_list, &Patch_list, M, Max_Degree);
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	totaltime += duration;
	std::cout << "Reduce Graph: " << duration << "s." << std::endl;

	Re_Graph.open(argc[8]);
	
	Re_Graph << M << " " << new_m << "\n";

	for (i=0; i<M; ++i)
	{
		Re_Graph << Patch_list[i].Index << " ";
		Re_Graph << std::setprecision(16) << Patch_list[i].Selfloop << " ";
		Re_Graph << Patch_list[i].Degree << "\n";
		tmpP=Patch_list[i].Neighbour;
		for (j=0; j<Patch_list[i].Degree; ++j)
		{	
			Re_Graph << tmpP->Index << " ";
			Re_Graph << std::setprecision(16) << tmpP->Diagonal << " ";
			Re_Graph << std::setprecision(16) << tmpP->Offdiag << " ";
			tmpP=tmpP->Next;
		}
		Re_Graph << "\n";

	}
	Re_Graph.close();

	std::cout << "Total time: " << totaltime << "s." << std::endl;

	return 0;
}

