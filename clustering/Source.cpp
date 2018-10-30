#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "DataStructure.h"
#include "BasisClustering.h"


int main(int argv, char* argc[])
{
	int k = 0;
	int temp_index;

	std::fstream fp;
	int N, M;
	fp.open(argc[1]);
	fp >> N;
	fp >> M;

	std::vector<Vertex> Vertex_list(N);
	std::vector<Vedge> Edge_list(M);

	double epsilon = atof(argc[2]);
	double ratio = atof(argc[3]);
	int maxPsizeControl = atoi(argc[4]);
	int level = atoi(argc[5]);

	for (int i = 0; i < N; ++i)
	{
		fp >> Vertex_list[i].Index;
		fp >> Vertex_list[i].Selfloop;
		fp >> Vertex_list[i].Degree;
		Vertex_list[i].Wdegree = 0;
		Vertex_list[i].Neighbour = &Edge_list[k];
		Vertex_list[i].Relativecor[0] = Vertex_list[i].Index;
		Vertex_list[i].Relativecor[1] = 0;
		Vertex_list[i].Eigfun = 1;
		for (int j = 0; j < Vertex_list[i].Degree; ++j)
		{
			fp >> temp_index;
			Edge_list[k].Outvertex = &Vertex_list[temp_index];

			fp >> Edge_list[k].diagonal;
			fp >> Edge_list[k].offdiagonal;

			Edge_list[k].diagonal = std::abs(Edge_list[k].diagonal);
			Edge_list[k].offdiagonal = Edge_list[k].offdiagonal;
			Vertex_list[i].totaldiagonal += Edge_list[k].diagonal;
			Vertex_list[i].totaloffdiagonal += std::abs(Edge_list[k].offdiagonal);

			Vertex_list[i].Wdegree += Edge_list[k].diagonal;

			if (j < Vertex_list[i].Degree - 1) Edge_list[k].Next = &Edge_list[k + 1];
			++k;
		}
	}

	Clustering(&Vertex_list, N, M, epsilon, ratio, maxPsizeControl, level);

	fp.close();

	return 0;
}