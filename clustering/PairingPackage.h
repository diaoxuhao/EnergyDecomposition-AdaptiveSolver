#pragma once
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include "MatchCheck.h"

// For updating neighbour patches
bool Compare_Patch_Index(Pedge *a, Pedge *b)
{
	return (a->Outpatch->Index)<(b->Outpatch->Index);
}

// For sort patches by Delta and Mark
bool Compare_Patch(Patch *a, Patch *b)
{
	if (a->Mark >= 0 && b->Mark >= 0)
		return ((a->Delta)>(b->Delta));
	else return (a->Mark>b->Mark);
}

void Combine_Patches(Patch *P1, Patch *P2, double delta, double lambda)
{
	P2->Index = P1->Index;

	P2->Cluster = P1;
	P1->Cluster = P1;

	P1->Delta = delta;
	P1->Lambda = lambda;

	P1->Mark = 1; // Absorb
	P2->Mark = 2; // Absorbed

				  // Combine members 
	P1->Psize = P1->Psize + P2->Psize;
	P1->Lastmember->Next = P2->Member;
	P1->Lastmember = P2->Lastmember;


	PVedge *PVe_temp;
	PVe_temp = P1->Member;
	for (int i = 0; i<P1->Psize; ++i)
	{
		PVe_temp->Outvertex->Eigfun = PVe_temp->Eigfun;
		PVe_temp->Outvertex->Relativecor[0] = P1->Index;
		PVe_temp->Outvertex->Relativecor[1] = i;
		PVe_temp = PVe_temp->Next;
	}

	// Combine neighbours
	P1->Degree = P1->Degree + P2->Degree;
	P1->Lastneighbour->Next = P2->Neighbour;
	P1->Lastneighbour = P2->Lastneighbour;


}

void Update_Neighbours(Patch *P)
{
	// Construct temp for sorting neighbour patches' index
	int t = 0, s = 0;
	int i;
	std::vector<Pedge *> temp(P->Degree);
	Pedge *q1, *q2;
	q1 = P->Neighbour;
	q2 = P->Neighbour;
	for (i = 0; i<P->Degree; ++i)
	{
		if (q1->Outpatch->Mark>0 || P->Mark>0)
		{
			temp[t] = q1;
			++t;
			q1 = q1->Next;
			if (s == 0)
			{
				q2 = q1;
				P->Neighbour = q1;
			}
			else
			{
				q2->Next = q1;
			}
		}
		else
		{
			if (s>0)
			{
				q2 = q1;
			}
			q1 = q1->Next;
			++s;
		}
	}

	if (t>0) std::sort(temp.begin(), temp.begin() + t, Compare_Patch_Index);

	// Compress neighbour patches
	P->Degree = s;
	for (i = 0; i<t; ++i)
	{
		if (temp[i]->Outpatch->Index == P->Index) continue;
		if (P->Degree == 0)
		{
			P->Neighbour = temp[i];
			++P->Degree;
			q2 = temp[i];
			q2->Outpatch = q2->Outpatch->Cluster;
			continue;
		}

		if (temp[i]->Outpatch->Index == q2->Outpatch->Index)
			q2->Bond += temp[i]->Bond;
		else
		{
			q2->Next = temp[i];
			q2 = temp[i];
			q2->Outpatch = q2->Outpatch->Cluster;
			++P->Degree;
		}
	}
	P->Lastneighbour = q2;
}

void Find_pair_small(Patch *P, double epsilon, double ratio, int maxPsizeControl)
{
	int i, s1 = 0, s2 = 0, s3 = 0;
	double delta = 0;
	double lambda = 0;
	double delta_match = 0;
	double lambda_match = 0;
	double Max_bond = 0;
	int in_exist = 1;

	Pedge *temp = P->Neighbour;
	Patch *Match_pair;

	for (i = 0; i<P->Degree; ++i)
	{
		if (temp->Outpatch->Mark>0)
		{
			s1 = 1;
		}

		if (temp->Outpatch->Mark <= 0 && (temp->Bond / temp->Outpatch->Psize) > Max_bond)
		{

			s2 = Check(P, temp->Outpatch, &delta, &lambda, epsilon, ratio, maxPsizeControl);
			//std::cout << s2 << std::endl;
			if (s2 == 1 && (P->Psize + temp->Outpatch->Psize <= 10*maxPsizeControl))
			{
				Max_bond = temp->Bond / temp->Outpatch->Psize;
				lambda_match = lambda;
				delta_match = delta;
				Match_pair = temp->Outpatch;
				s3 = 1;
			}
		}
		temp = temp->Next;
	}

	if (s3 == 1)
	{
		s2 = Check(P, Match_pair, &delta, &lambda, epsilon, ratio, maxPsizeControl);
		Combine_Patches(P, Match_pair, delta_match, lambda_match);
	}
	else if (s1 == 0)
	{
		P->Mark = -1;
	}

}


void Find_pair(Patch *P, double epsilon, double ratio, int maxPsizeControl)
{
	int i, s1 = 0, s2 = 0;
	double Max_bond = 0;
	double delta = 0;
	double lambda = 0;

	Pedge *temp = P->Neighbour;
	Patch *Match_pair;
	for (i = 0; i<P->Degree; ++i)
	{
		if (temp->Outpatch->Mark>0)
		{
			s1 = 1;
			//std::cout << temp->Outpatch->Mark << std::endl;
		}
		if (temp->Outpatch->Mark <= 0 && (temp->Bond / temp->Outpatch->Psize) > Max_bond)
		{
			//std::cout << "Patch ID: " << temp->Outpatch->Index << std::endl;
			Match_pair = temp->Outpatch;
			Max_bond = temp->Bond / temp->Outpatch->Psize;
		}
		temp = temp->Next;
	}
	if (Max_bond>0)
	{
		s2 = Check(P, Match_pair, &delta, &lambda, epsilon, ratio, maxPsizeControl);

	}
	if (s2 == 1)
	{	
		//std::cout << "here" << std::endl;
		//std::cout << "Combining Patch " << P->Index << " with Patch " << Match_pair->Index << std::endl;
		Combine_Patches(P, Match_pair, delta, lambda);
		//std::cout << "here out" << std::endl;
	}
	else if (s1 == 0)
	{
		P->Mark = -1;
	}
}

void Pairing_Once(std::vector<Patch *> *Patch_key, int Patch_num[], double epsilon, double ratio, int maxPsizeControl)
{
	int Newnum = 0, Newactivenum = 0;
	Patch *temp_pair; // for storing matched pair

	for (int i = 0; i < Patch_num[1]; ++i)
	{
		if ((*Patch_key)[i]->Mark == 0)
		{
			if ((*Patch_key)[i]->Psize <= 2)  // d=2
			{
				Find_pair_small((*Patch_key)[i], epsilon, ratio, maxPsizeControl);
			}
			else
			{
				Find_pair((*Patch_key)[i], epsilon, ratio, maxPsizeControl);
			}
		}
	}
	// Update neighbour patches, and compress Patch_key
	// Push the remaining patches pointers to the front
	for (int i = 0; i < Patch_num[0]; ++i)
	{
		if ((*Patch_key)[i]->Mark != 2)
		{
			if ((*Patch_key)[i]->Mark != -1)
			{
				++Newactivenum;
			}
			Update_Neighbours((*Patch_key)[i]);
			(*Patch_key)[Newnum] = (*Patch_key)[i];
			++Newnum;
		}
	}

	// Push the active patches to the front and sort them by Delta
	if (Newnum<Patch_num[1])
	{
		std::stable_sort((*Patch_key).begin(), (*Patch_key).begin() + Newnum, Compare_Patch);
	}
	else
	{
		std::stable_sort((*Patch_key).begin(), (*Patch_key).begin() + Patch_num[1], Compare_Patch);
	}

	for (int i = 0; i<Newactivenum; ++i)
	{
		(*Patch_key)[i]->Mark = 0;
	}

	Patch_num[0] = Newnum;
	Patch_num[1] = Newactivenum;

}

int Parallel_Pairing(std::vector<Patch *> *Patch_key, int initial_num, double epsilon, double ratio, int maxPsizeControl)
{
	int Patch_num[2] = { initial_num,initial_num };
	std::stable_sort((*Patch_key).begin(), (*Patch_key).begin() + initial_num, Compare_Patch);


	while (Patch_num[1]>0)
	{
		Pairing_Once(Patch_key, Patch_num, epsilon, ratio, maxPsizeControl);
		//system("Pause");
		std::cout << Patch_num[1] << " " << Patch_num[0] << std::endl;
	}
	


	for (int i = 0; i<Patch_num[0]; ++i)
		(*Patch_key)[i]->Index = i;

	return Patch_num[0];
}



