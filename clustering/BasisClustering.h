#pragma once
#include <stdio.h>
#include <math.h>
#include <ctime>
#include <algorithm>
#include <fstream>
#include "PairingPackage.h"


void Clustering(std::vector<Vertex> *Vertex_list, int N, int M, double epsilon, double ratio, int maxPsizeControl, int level)
{
	std::vector<Patch> Patch_list(N);
	std::vector<Pedge> Pedge_list(M);
	std::vector<Patch *> Patch_key(N);
	std::vector<PVedge> Patch_vertex_list(N);
	
	Vedge *Ve_temp;
	Pedge *Pe_temp;
	PVedge *PVe_temp;

	// Initialize Patch_list
	int k = 0;
	for (int i = 0; i < N; ++i)
	{
		Patch_key[i] = &Patch_list[i];

		Patch_list[i].Index = i;
		Patch_list[i].Mark = 0;
		Patch_list[i].Delta = (*Vertex_list)[i].Selfloop + (*Vertex_list)[i].totaldiagonal + (*Vertex_list)[i].totaloffdiagonal;
		Patch_list[i].Cluster = &Patch_list[i];

		// Initialize members
		Patch_list[i].Psize = 1;
		Patch_vertex_list[i].Outvertex = &(*Vertex_list)[i];
		Patch_vertex_list[i].Next = NULL;
		Patch_list[i].Member = &Patch_vertex_list[i];
		Patch_list[i].Lastmember = &Patch_vertex_list[i];

		// Initialize neighbours
		Patch_list[i].Degree = (*Vertex_list)[i].Degree;
		Patch_list[i].Neighbour = &Pedge_list[k];
		Ve_temp = (*Vertex_list)[i].Neighbour;
		for (int j = 0; j < (*Vertex_list)[i].Degree; ++j)
		{
			//Pedge_list[k].Bond = Ve_temp->Weight;
			Pedge_list[k].Bond = std::abs(Ve_temp->offdiagonal);

			Pedge_list[k].Outpatch = &Patch_list[Ve_temp->Outvertex->Index];
			if (j<(*Vertex_list)[i].Degree - 1) Pedge_list[k].Next = &Pedge_list[k + 1];
			++k;
			Ve_temp = Ve_temp->Next;
		}
		Patch_list[i].Lastneighbour = &Pedge_list[k - 1];

	}

	std::cout << N << " " << M << std::endl;

	int Patch_N = 0, s = 0, Pedge_M = 0;

	std::clock_t start;
	double duration;
	start = std::clock();

	Patch_N = Parallel_Pairing(&Patch_key, N, epsilon, ratio, maxPsizeControl);

	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "Computation time: " << duration << "s." << std::endl;

	std::cout << std::endl;
	std::cout << Patch_N << std::endl;

	//duration;
	//start = std::clock();

	for (int i = 0; i < Patch_N; ++i)
	{
		if (Patch_key[i]->Psize > 1)
		{
			Eigen::VectorXd u = Eigen::VectorXd::Zero(Patch_key[i]->Psize);
			PVe_temp = Patch_key[i]->Member;
			for (int j = 0; j < Patch_key[i]->Psize; ++j)
			{
				u(j) = PVe_temp->Outvertex->Eigfun;
				PVe_temp = PVe_temp->Next;
			}

			double nu = u.norm();
			if (nu != 0)
			{
				u /= nu;
				u(0) += ((u(0) >= 0) - (u(0) < 0));
				u /= sqrt(std::abs(u(0)));
			}
			else
			{
				u(0) = sqrt(2);
			}
			u /= u.norm();
			PVe_temp = Patch_key[i]->Member;
			for (int j = 0; j < Patch_key[i]->Psize; ++j)
			{
				PVe_temp->Outvertex->U = u(j);
				PVe_temp = PVe_temp->Next;
			}
			//std::cout << "HV = " << u.transpose() << std::endl;
		}
		else
		{
			Patch_key[i]->Member->Outvertex->U = 1;
		}
	}

	//duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	//std::cout << "Computation time: " << duration << "s." << std::endl;


	std::ofstream ClusterOutput;
	std::ofstream NeighbourOutput;

	std::string clusterfile = "cluster" + std::to_string(level) + ".txt";
	std::string neighborfile = "neighbour" + std::to_string(level) + ".txt";
	ClusterOutput.open(clusterfile);
	NeighbourOutput.open(neighborfile);

	ClusterOutput << Patch_N << "\n";
	NeighbourOutput << Patch_N << "\n";

	double max_delta = 0;

	for (int i = 0; i < Patch_N; ++i)
	{
		Pedge_M += Patch_key[i]->Degree;
		s += Patch_key[i]->Psize;
		ClusterOutput << Patch_key[i]->Index << " ";
		ClusterOutput << Patch_key[i]->Delta << " ";

		NeighbourOutput << Patch_key[i]->Index << " ";

		if (max_delta<  Patch_key[i]->Delta) max_delta = Patch_key[i]->Delta;

		ClusterOutput << Patch_key[i]->Lambda << " ";
		ClusterOutput << "\n";
		ClusterOutput << Patch_key[i]->Psize << " ";
		PVe_temp = Patch_key[i]->Member;
		for (int j = 0; j<Patch_key[i]->Psize; ++j)
		{
			ClusterOutput << PVe_temp->Outvertex->Index << " ";
			ClusterOutput << PVe_temp->Outvertex->Eigfun << " ";
			ClusterOutput << PVe_temp->Outvertex->U << " ";
			PVe_temp = PVe_temp->Next;
		}
		ClusterOutput << "\n";

		NeighbourOutput << Patch_key[i]->Degree << "\n";
		Pe_temp = Patch_key[i]->Neighbour;
		for (int j = 0; j<Patch_key[i]->Degree; ++j)
		{
			NeighbourOutput << Pe_temp->Outpatch->Index << " ";
			Pe_temp = Pe_temp->Next;
		}
		NeighbourOutput << "\n";
	}

	ClusterOutput.close();
	NeighbourOutput.close();

	std::cout << max_delta << std::endl;
	std::cout << s << std::endl;



	// Construct U;
	//duration;
	//start = std::clock();
	
	//std::vector<Eigen::MatrixXd> U_list(Patch_N);
	//for (int i = 0; i < Patch_N; ++i)
	//{
	//	U_list[i] = Eigen::MatrixXd::Zero(Patch_key[i]->Psize, Patch_key[i]->Psize);
	//}

	//for (int i = 0; i < Patch_N; ++i)
	//{
	//	Eigen::MatrixXd R = Eigen::MatrixXd::Zero(Patch_key[i]->Psize, Patch_key[i]->Psize);
	//	
	//	for (int j = 1; j < Patch_key[i]->Psize; ++j)
	//	{
	//		for (int k = 1; k < Patch_key[i]->Psize; ++k)
	//		{
	//			R(j, k) = rand() / (RAND_MAX + 1.0);
	//		}
	//	}

	//	PVe_temp = Patch_key[i]->Member;
	//	for (int j = 0; j < Patch_key[i]->Psize; ++j)
	//	{
	//		R(j,0) = PVe_temp->Outvertex->Eigfun;
	//		PVe_temp = PVe_temp->Next;
	//	}
	//	
	//	for (int k = 0; k < Patch_key[i]->Psize; ++k)
	//	{
	//		Eigen::VectorXd u(Patch_key[i]->Psize - k);
	//		for (int t = 0; t < Patch_key[i]->Psize - k; ++t)
	//		{
	//			u(t) = R(t + k, k);
	//		}
	//		
	//		// Eigen::VectorXd u = R.col(k).tail(Patch_key[i]->Psize - k + 1);
	//		double nu = u.norm();
	//		if (nu != 0)
	//		{
	//			u /= nu;
	//			double sgn = ((u(0) >= 0) - (u(0) < 0));
	//			u(0) = u(0) + sgn;
	//			u = u / sqrt(abs(u(0)));
	//		}
	//		else
	//		{
	//			u(0) = sqrt(2);
	//		}
	//		for (int t = k; t < Patch_key[i]->Psize; ++t)
	//		{
	//			U_list[i].coeffRef(t, k) = u(t - k);
	//		}
	//		
	//		//for (int t = k; t < Patch_key[i]->Psize; ++t)
	//		//{
	//		//	U_list[i].coeffRef(t, k) = u(t - k);
	//		//}
	//		R.block(k, k, Patch_key[i]->Psize - k, Patch_key[i]->Psize - k) -= u* (u.transpose() * R.block(k, k, Patch_key[i]->Psize - k, Patch_key[i]->Psize - k));
	//	}
	//	//std::cout << "Householder vector for Patch" << i << ":" << std::endl;
	//	//std::cout << U_list[i] << std::endl;
	//}

	////duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	////std::cout << "Computation time: " << duration << "s." << std::endl;


	// for (i=0;i<Patch_N;++i)
	// {	
	// 	s=s+Patch_key[i]->Psize;
	// 	std::cout<< Patch_key[i]->Psize <<" ";
	// 	Ve_temp=Patch_key[i]->Member;
	// for (j=0;j<Patch_key[i]->Psize;++j)
	// {
	// 	std::cout<<Ve_temp->Outvertex->Index<<" ";
	// 	Ve_temp=Ve_temp->Next;
	// }
	// std::cout<<std::endl;
	// }

	// for (i=0;i<N;++i)
	// {
	// 	std::cout<< Patch_key[i]->Index<<" "<< Patch_list[i].Delta<< " ";
	// 	Pe_temp=Patch_list[i].Neighbour;
	// 	for (j=0;j<Vertex_list[i].Degree;++j)
	// 	{
	// 		std::cout<< Pe_temp->Bond<< " ";
	// 		Pe_temp=Pe_temp->Next;
	// 	}
	// 	std::cout<< std::endl;
	// }

}