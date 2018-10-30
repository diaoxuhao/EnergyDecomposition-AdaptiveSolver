#include <stdio.h>
#include <math.h>
#include <vector>
#include <ctime>

void ReduceGraph(std::vector<Vertex> *Vertex_list, std::vector<Patch> *Patch_list, int M, int Max_Degree)
{
	int i,j,k,t1,t2;
	double psi;
	Pedge *tmpP;
	PVedge *tmpPV;
	Vedge *tmpV;
	std::vector<double> tmpDiagonal(Max_Degree);
	std::vector<double> tmpOffdiag(Max_Degree);


	for (i=0; i<M; ++i)
	{
		tmpP=(*Patch_list)[i].Neighbour;
		for (j=0; j<(*Patch_list)[i].Degree; ++j)
		{
			(*Patch_list)[tmpP->Index].Mark=j;
			tmpDiagonal[j]=0;
			tmpOffdiag[j]=0;
			tmpP=tmpP->Next;
		}
		tmpPV=(*Patch_list)[i].Member;
		for (j=0; j<(*Patch_list)[i].Psize; ++j)
		{
			t1=tmpPV->Index;
			psi=(*Vertex_list)[t1].Psi;
			(*Patch_list)[i].Selfloop += psi*psi*(*Vertex_list)[i].Selfloop;
			tmpV=(*Vertex_list)[t1].Neighbour;
			for (k=0; k<(*Vertex_list)[t1].Degree; ++k)
			{
				t2=tmpV->Index;
				if ((*Vertex_list)[t2].Pindex==i)
				{
					(*Patch_list)[i].Selfloop += psi*psi*tmpV->Diagonal + psi*(*Vertex_list)[t2].Psi*tmpV->Offdiag;
				}
				else
				{
					tmpDiagonal[(*Patch_list)[(*Vertex_list)[t2].Pindex].Mark] += psi*psi*tmpV->Diagonal;
					tmpOffdiag[(*Patch_list)[(*Vertex_list)[t2].Pindex].Mark] += psi*(*Vertex_list)[t2].Psi*tmpV->Offdiag;
				}
				tmpV=tmpV->Next;
			}
			tmpPV=tmpPV->Next;
		}
		tmpP=(*Patch_list)[i].Neighbour;
		for (j=0; j<(*Patch_list)[i].Degree; ++j)
		{
			tmpP->Diagonal=tmpDiagonal[j];
			tmpP->Offdiag=tmpOffdiag[j];
			tmpP=tmpP->Next;
		}
	}
}
