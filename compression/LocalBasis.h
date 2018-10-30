#include <stdio.h>
#include <math.h>
#include <vector>
#include <ctime>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCore>
#include <omp.h>

typedef Eigen::Triplet<double> triplet;

void Extend_Patch(std::vector<int> *tmpsup, std::vector<int> *seed, std::vector<int> *segment, std::vector<Vertex> *Vertex_list, std::vector<Patch> *Patch_list, int *seed_str, int *seed_end, int *sup_str, int *sup_end, int *Pnum)
{
	Pedge *tmpP;
	PVedge *tmpPV;
	Vedge* tmpV;
	int t1,t2,i,j,k;
	int new_seed_end=*seed_end;

	*sup_str=*sup_end;

	for (i=*seed_str; i<*seed_end; ++i)
	{
		t1=(*seed)[i];
		tmpP=(*Patch_list)[t1].Neighbour;
		for (j=0; j<(*Patch_list)[t1].Degree; ++j)
		{
			t2=tmpP->Index;
			if ((*Patch_list)[t2].Mark==0)
			{
				(*seed)[new_seed_end]=t2;
				++new_seed_end;
				(*Patch_list)[t2].Mark=1;
				tmpPV=(*Patch_list)[t2].Member;
				for (k=0; k<(*Patch_list)[t2].Psize; ++k)
				{
					(*tmpsup)[*sup_end]=tmpPV->Index;
					(*Vertex_list)[tmpPV->Index].Relative=*sup_end;
					++(*sup_end);
					tmpPV=tmpPV->Next;
				}			
				*Pnum=*Pnum+1;
				(*segment)[*Pnum]=*sup_end;
			}
			tmpP=tmpP->Next;
		}
	}
	*seed_str=*seed_end;
	*seed_end=new_seed_end;

}

void Extend_Patch(std::vector<int> *tmpsup, std::vector<int> *seed, std::vector<int> *segment, std::vector<Vertex> *Vertex_list, std::vector<Patch> *Patch_list, int *seed_str, int *seed_end, int *sup_str, int *sup_end, int *Pnum, std::vector<int> *private_mark, std::vector<int> *private_relative)
{
	Pedge *tmpP;
	PVedge *tmpPV;
	Vedge* tmpV;
	int t1, t2, i, j, k;
	int new_seed_end = *seed_end;

	*sup_str = *sup_end;

	for (i = *seed_str; i<*seed_end; ++i)
	{
		t1 = (*seed)[i];
		tmpP = (*Patch_list)[t1].Neighbour;
		for (j = 0; j<(*Patch_list)[t1].Degree; ++j)
		{
			t2 = tmpP->Index;
			if ((*private_mark)[t2] == 0)
			{
				(*seed)[new_seed_end] = t2;
				++new_seed_end;
				(*private_mark)[t2] = 1;
				tmpPV = (*Patch_list)[t2].Member;
				for (k = 0; k<(*Patch_list)[t2].Psize; ++k)
				{
					(*tmpsup)[*sup_end] = tmpPV->Index;
					(*private_relative)[tmpPV->Index] = *sup_end;
					++(*sup_end);
					tmpPV = tmpPV->Next;
				}
				*Pnum = *Pnum + 1;
				(*segment)[*Pnum] = *sup_end;
			}
			tmpP = tmpP->Next;
		}
	}
	*seed_str = *seed_end;
	*seed_end = new_seed_end;

}

void Extend_A(std::vector<Vertex> *Vertex_list, std::vector<Patch> *Patch_list, std::vector<triplet> *tp, std::vector<double> *tmpw, std::vector<int> *tmpsup, int sup_str, int sup_end)
{
	int i, j, t1, t2;
	Vedge *tmpV;

	for (i = sup_str; i<sup_end; ++i)
	{
		t1 = (*tmpsup)[i];
		(*tp).push_back(triplet(i, i, (*Vertex_list)[t1].Diagonal));

		(*tmpw)[i] = (*Vertex_list)[t1].Hous;
		tmpV = (*Vertex_list)[t1].Neighbour;
		for (j = 0; j<(*Vertex_list)[t1].Degree; ++j)
		{
			t2 = tmpV->Index;

			if ((*Patch_list)[(*Vertex_list)[t2].Pindex].Mark == 1 && (*Vertex_list)[t2].Relative < i)
			{
				(*tp).push_back(triplet(i, (*Vertex_list)[t2].Relative, tmpV->Offdiag));
				(*tp).push_back(triplet((*Vertex_list)[t2].Relative, i, tmpV->Offdiag));
			}
			tmpV = tmpV->Next;
		}
	}
}

void Extend_A(std::vector<Vertex> *Vertex_list, std::vector<Patch> *Patch_list, std::vector<triplet> *tp, std::vector<double> *tmpw, std::vector<int> *tmpsup, int sup_str, int sup_end, std::vector<int> *private_mark, std::vector<int> *private_relative, int recent_loop)
{
	int i, j, t1, t2;
	Vedge *tmpV;

	for (i = sup_str; i<sup_end; ++i)
	{
		t1 = (*tmpsup)[i];
		(*tp).push_back(triplet(i, i, (*Vertex_list)[t1].Diagonal));

		(*tmpw)[i] = (*Vertex_list)[t1].Hous;
		tmpV = (*Vertex_list)[t1].Neighbour;
		for (j = 0; j<(*Vertex_list)[t1].Degree; ++j)
		{
			t2 = tmpV->Index;

			if ((*private_mark)[(*Vertex_list)[t2].Pindex] == 1 && (*private_relative)[t2] < i)
			{
				(*tp).push_back(triplet(i, (*private_relative)[t2], tmpV->Offdiag));
				(*tp).push_back(triplet((*private_relative)[t2], i, tmpV->Offdiag));
			}
			tmpV = tmpV->Next;
		}
	}
}


void Extend_A(std::vector<Vertex> *Vertex_list, std::vector<Patch> *Patch_list, std::vector<triplet> *tp, std::vector<double> *tmpw, std::vector<int> *tmpsup, int sup_str, int sup_end, std::vector<int> private_mark, std::vector<int> private_relative)
{
	int i, j, t1, t2;
	Vedge *tmpV;

	for (i = sup_str; i<sup_end; ++i)
	{
		t1 = (*tmpsup)[i];
		(*tp).push_back(triplet(i, i, (*Vertex_list)[t1].Diagonal));

		(*tmpw)[i] = (*Vertex_list)[t1].Hous;
		tmpV = (*Vertex_list)[t1].Neighbour;
		for (j = 0; j<(*Vertex_list)[t1].Degree; ++j)
		{
			t2 = tmpV->Index;

			if (private_mark[(*Vertex_list)[t2].Pindex] == 1 && private_relative[t2] < i)
			{
				(*tp).push_back(triplet(i, private_relative[t2], tmpV->Offdiag));
				(*tp).push_back(triplet(private_relative[t2], i, tmpV->Offdiag));
			}
			tmpV = tmpV->Next;
		}
	}
}

void Multiply_Utrans(Eigen::VectorXd &u, Eigen::VectorXd &v, std::vector<double> *tmpw, std::vector<int> *segment, int Pnum)
{
	int i,j;
	double t;
	int k=0;

	for (i=0; i<Pnum; ++i)
	{
		t=0;
		for (j=(*segment)[i]; j<(*segment)[i+1]; ++j)
		{
			t=t+u(j)*(*tmpw)[j];
		}
		for (j=0; j<(*segment)[i+1]-(*segment)[i]-1; ++j)
		{
			v(k+j)=u((*segment)[i]+j+1)-2*t*(*tmpw)[(*segment)[i]+j+1];
		}
		k=k+(*segment)[i+1]-(*segment)[i]-1;
	}
}

void Multiply_U(Eigen::VectorXd &v, Eigen::VectorXd &u, std::vector<double> *tmpw, std::vector<int> *segment, int Pnum)
{
	int i,j;
	double t;
	int k=0;

	for (i=0; i<Pnum; ++i)
	{
		t=0;
		for (j=0; j<(*segment)[i+1]-(*segment)[i]-1; ++j)
		{
			t=t+v(k+j)*(*tmpw)[(*segment)[i]+j+1];
		}
		u((*segment)[i])=-2*t*(*tmpw)[(*segment)[i]];
		for (j=0; j<(*segment)[i+1]-(*segment)[i]-1; ++j)
		{
			u((*segment)[i]+j+1)=v(k+j)-2*t*(*tmpw)[(*segment)[i]+j+1];
		}
		k=k+(*segment)[i+1]-(*segment)[i]-1;

	}
}

double Extend_Psi(std::vector<Vertex> *Vertex_list, std::vector<triplet> *tp, int tp_length, std::vector<double> *tmppsi, std::vector<double> *tmpw, std::vector<int> *tmpsup, int sup_end, std::vector<int> *segment, int Pnum, double epsilon )
{
	Eigen::VectorXd r(sup_end-Pnum);
	Eigen::VectorXd v(sup_end);
	Eigen::VectorXd p = Eigen::VectorXd::Zero(sup_end-Pnum);
	Eigen::VectorXd Ap(sup_end-Pnum);
	Eigen::VectorXd x = Eigen::VectorXd::Zero(sup_end-Pnum);

	Eigen::SparseMatrix<double> Ar(sup_end,sup_end);
	int i,j;
	double alpha, beta, a0, a1, norm_b;

	Ar.setFromTriplets((*tp).begin(), (*tp).begin() + tp_length);

	for (i=0; i<sup_end; ++i)
	{
		v(i)=(*tmppsi)[i];
	}

	v=Ar*v;

	Multiply_Utrans(v,r,tmpw,segment,Pnum);

	norm_b=r.dot(r);
	a1=norm_b;
	beta=0;

	while ( a1>epsilon*norm_b )
	{
		p=r+beta*p;
		Multiply_U(p,v,tmpw,segment,Pnum);
		v=Ar*v;
		Multiply_Utrans(v,Ap,tmpw,segment,Pnum);
		alpha=a1/(p.dot(Ap));
		x=x+alpha*p;
		r=r-alpha*Ap;
		a0=a1;
		a1=r.dot(r);
		beta=a1/a0;
	}

	Multiply_U(x,v,tmpw,segment,Pnum);

	for (i=0; i<sup_end; ++i)
	{
		(*tmppsi)[i] = (*tmppsi)[i]-v[i];
	}

	return v.dot(Ar*v);

}

int LocalPsi(std::vector<Vertex> *Vertex_list, std::vector<Patch> *Patch_list, std::vector<triplet> *tp1, std::vector<int> *supsize, double epsilon, int M, int s)
{
	#pragma omp parallel
	{
		std::vector<triplet> private_tp;
		std::vector<int> private_seed;

		#pragma omp for
		for (int i = 0; i < M; ++i)
		{
			std::vector<triplet> tp;
			std::vector<double> tmpw(s);
			std::vector<double> tmppsi(s);
			std::vector<int> tmpsup(s);
			std::vector<int> seed(s);
			std::vector<int> segment(s);

			std::vector<int> private_mark((*Patch_list).size());
			std::vector<int> private_relative((*Vertex_list).size());

			for (int k = 0; k < private_mark.size(); ++k)
			{
				private_mark[k] = 0;
			}
			for (int k = 0; k < private_relative.size(); ++k)
			{
				private_relative[k] = (*Vertex_list)[k].Relative;
			}

			int Pnum;
			int seed_str, seed_end;
			int tp_length;
			int sup_str, sup_end;

			double energy_gap, new_energy_gap;
			double ratio;
			PVedge *tmpPV;

			Pnum = 0;
			seed_str = 0;
			seed_end = 1;
			sup_str = 0;
			segment[0] = 0;
			private_mark[i] = 1;
			sup_end = (*Patch_list)[i].Psize;
			++Pnum;
			segment[Pnum] = sup_end;
			tmpPV = (*Patch_list)[i].Member;
			seed[0] = i;
			for (int j = 0; j < sup_end; ++j)
			{
				tmpsup[j] = tmpPV->Index;
				tmppsi[j] = (*Vertex_list)[tmpPV->Index].Phi;
				private_relative[tmpPV->Index] = j;
				tmpPV = tmpPV->Next;
			}

			Extend_A(Vertex_list, Patch_list, &tp, &tmpw, &tmpsup, sup_str, sup_end, &private_mark, &private_relative, i);
			tp_length = tp.size();

			energy_gap = Extend_Psi(Vertex_list, &tp, tp_length, &tmppsi, &tmpw, &tmpsup, sup_end, &segment, Pnum, epsilon);
			
			Extend_Patch(&tmpsup, &seed, &segment, Vertex_list, Patch_list, &seed_str, &seed_end, &sup_str, &sup_end, &Pnum, &private_mark, &private_relative);
			
			Extend_A(Vertex_list, Patch_list, &tp, &tmpw, &tmpsup, sup_str, sup_end, &private_mark, &private_relative, i);
			tp_length = tp.size();

			for (int j = sup_str; j < sup_end; ++j)
			{
				tmppsi[j] = 0;
			}

			new_energy_gap = Extend_Psi(Vertex_list, &tp, tp_length, &tmppsi, &tmpw, &tmpsup, sup_end, &segment, Pnum, epsilon);
			
			ratio = std::abs(new_energy_gap / energy_gap);
			energy_gap = new_energy_gap;
			
			while (ratio != 1 && ratio*std::abs(energy_gap) > (1 - ratio)*epsilon)
			{
				Extend_Patch(&tmpsup, &seed, &segment, Vertex_list, Patch_list, &seed_str, &seed_end, &sup_str, &sup_end, &Pnum, &private_mark, &private_relative);

				Extend_A(Vertex_list, Patch_list, &tp, &tmpw, &tmpsup, sup_str, sup_end, &private_mark, &private_relative, i);
				tp_length = tp.size();

				for (int j = sup_str; j < sup_end; ++j)
				{	
					tmppsi[j] = 0;
				}

				new_energy_gap = Extend_Psi(Vertex_list, &tp, tp_length, &tmppsi, &tmpw, &tmpsup, sup_end, &segment, Pnum, epsilon);

				ratio = new_energy_gap / energy_gap;
				energy_gap = new_energy_gap;

			}

			private_tp.resize(0);
			for (int j = 0; j < sup_end; ++j)
			{
				private_tp.push_back(triplet(tmpsup[j], i, tmppsi[j]));
			}
			private_seed.resize(seed_end);
			for (int j = 0; j < seed_end; ++j)
			{
				private_seed[j] = seed[j];
			}
			(*supsize)[i] = sup_end;

			#pragma omp critical
			(*tp1).insert((*tp1).end(), private_tp.begin(), private_tp.end());

			#pragma omp critical
			for (int j = 0; j<private_seed.size(); ++j)
			{
				(*Patch_list)[private_seed[j]].Mark = 0;
			}
		}
	}

	return 1;
}

int LocalPsi_r(std::vector<Vertex> *Vertex_list, std::vector<Patch> *Patch_list, std::vector<triplet> *tp2, std::vector<int> *supsize, double epsilon, int M, int s, int radius  )
{
	
	int k,i,j;
	int Pnum;
	std::vector<int> segment(s);
	int seed_str, seed_end, sup_str, sup_end;
	Eigen::SparseMatrix<double> A(s,s);
	std::vector<double> tmppsi(s);
	std::vector<double> tmpw(s);
	double energy_gap, new_energy_gap;
	std::vector<int> tmpsup(s);
	std::vector<int> seed(s);
	double ratio;
	PVedge *tmpPV;

	for (i = 0; i < M; ++i)
	{
		std::vector<triplet> tp;
		int tp_length;

		Pnum=0;
		seed_str=0;
		seed_end=1;
		sup_str=0;
		segment[0]=0;
		(*Patch_list)[i].Mark=1;
		sup_end=(*Patch_list)[i].Psize;
		++Pnum;
		segment[Pnum]=sup_end;
		tmpPV=(*Patch_list)[i].Member;
		seed[0]=i;
		for (j=0; j<sup_end; ++j)
		{
			tmpsup[j]=tmpPV->Index;
			tmppsi[j]=(*Vertex_list)[tmpPV->Index].Phi;
			(*Vertex_list)[tmpPV->Index].Relative=j;
			tmpPV=tmpPV->Next;
		}
		Extend_A(Vertex_list, Patch_list, &tp, &tmpw, &tmpsup, sup_str, sup_end);
		tp_length = tp.size();

		Extend_Psi(Vertex_list, &tp, tp_length, &tmppsi, &tmpw, &tmpsup, sup_end, &segment, Pnum, epsilon);

		for (j=0; j<radius; ++j)
		{
			Extend_Patch(&tmpsup, &seed, &segment, Vertex_list, Patch_list, &seed_str, &seed_end, &sup_str, &sup_end, &Pnum);

			Extend_A(Vertex_list, Patch_list, &tp, &tmpw, &tmpsup, sup_str, sup_end);
			tp_length = tp.size();

			for (j=sup_str; j<sup_end; ++j) tmppsi[j]=0;

			Extend_Psi(Vertex_list, &tp, tp_length, &tmppsi, &tmpw, &tmpsup, sup_end, &segment, Pnum, epsilon);
		}

		for (j = 0; j<sup_end; ++j)
		{
			(*tp2).push_back(triplet(tmpsup[j], i, tmppsi[j]));
		}

		for (j=0; j<seed_end; ++j)
		{
			(*Patch_list)[seed[j]].Mark=0;
		}
		(*supsize)[i]=sup_end;
	}
	return 1;
}

void LocalPsi_0(std::vector<Vertex> *Vertex_list, std::vector<Patch> *Patch_list, double epsilon, int M, int s)
{
	#pragma omp parallel
	{
		#pragma omp for
		for (int i = 0; i < M; ++i)
		{
			std::vector<triplet> tp;
			std::vector<int> segment(s);
			Eigen::SparseMatrix<double> A(s, s);
			std::vector<double> tmppsi(s);
			std::vector<double> tmpw(s);
			std::vector<int> tmpsup(s);
			PVedge *tmpPV;

			std::vector<int> private_mark((*Patch_list).size());
			std::vector<int> private_relative((*Vertex_list).size());
			for (int k = 0; k < private_mark.size(); ++k)
			{
				private_mark[k] = 0;
			}
			for (int k = 0; k < private_relative.size(); ++k)
			{
				private_relative[k] = (*Vertex_list)[k].Relative;
			}

			int sup_str, sup_end;
			int tp_length;

			segment[0] = 0;
			sup_end = (*Patch_list)[i].Psize;
			segment[1] = sup_end;
			private_mark[i] = 1;
			tmpPV = (*Patch_list)[i].Member;
			for (int j = 0; j < sup_end; ++j)
			{
				tmpsup[j] = tmpPV->Index;
				tmppsi[j] = (*Vertex_list)[tmpPV->Index].Phi;
				private_relative[tmpPV->Index] = j;
				tmpPV = tmpPV->Next;
			}
			Extend_A(Vertex_list, Patch_list, &tp, &tmpw, &tmpsup, 0, sup_end, &private_mark, &private_relative, i);
			
			tp_length = tp.size();

			Extend_Psi(Vertex_list, &tp, tp_length, &tmppsi, &tmpw, &tmpsup, sup_end, &segment, 1, epsilon);

			#pragma omp critical
			for (int j = 0; j < sup_end; ++j)
			{
				(*Vertex_list)[tmpsup[j]].Psi = tmppsi[j];
			}

			#pragma omp critical
			(*Patch_list)[i].Mark = 0;
		}
	}
}

//void LocalPsi_0(std::vector<Vertex> *Vertex_list, std::vector<Patch> *Patch_list, double epsilon, int M, int s)
//{
//
//	int k, i, j;
//	std::vector<int> segment(s);
//	int sup_str, sup_end;
//	Eigen::SparseMatrix<double> A(s, s);
//	std::vector<double> tmppsi(s);
//	std::vector<double> tmpw(s);
//	std::vector<int> tmpsup(s);
//	PVedge *tmpPV;
//
//	for (i = 0; i<M; ++i)
//	{
//		std::vector<triplet> tp;
//		int tp_length;
//
//		segment[0]=0;
//		sup_end = (*Patch_list)[i].Psize;
//		segment[1]=sup_end;
//		(*Patch_list)[i].Mark=1;
//		tmpPV = (*Patch_list)[i].Member;
//		for (j = 0; j<sup_end; ++j)
//		{
//			tmpsup[j] = tmpPV->Index;
//			tmppsi[j] = (*Vertex_list)[tmpPV->Index].Phi;
//			(*Vertex_list)[tmpPV->Index].Relative = j;
//			tmpPV = tmpPV->Next;
//		}
//		Extend_A(Vertex_list, Patch_list, &tp, &tmpw, &tmpsup, 0, sup_end);
//		tp_length = tp.size();
//
//		Extend_Psi(Vertex_list, &tp, tp_length, &tmppsi, &tmpw, &tmpsup, sup_end, &segment, 1, epsilon);
//
//		for (j = 0; j<sup_end; ++j)
//		{
//			(*Vertex_list)[tmpsup[j]].Psi = tmppsi[j];
//		}
//
//		(*Patch_list)[i].Mark=0;
//
//	}
//}