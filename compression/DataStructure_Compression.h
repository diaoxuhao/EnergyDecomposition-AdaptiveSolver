#include <stdio.h>
#include <iostream>
#include <math.h>


struct Vertex;
struct Patch;
struct Vedge;
struct Pedge;
struct PVedge;
struct Psi;


struct Vertex
{
	int Index;
	int Degree;
	double Diagonal;
	double Selfloop;
	Vedge *Neighbour;
	int Pindex;
	int Relative;
	double Phi;
	double Hous;
	double Psi;
};

struct Patch
{
	int Index;
	int Mark;
	int Degree;
	int Psize;
	double Selfloop;
	PVedge *Member;
	PVedge *Lastmember;
	Pedge *Neighbour;
	Pedge *Lastneighbour;
};

struct Vedge
{
	double Diagonal;
	double Offdiag;
	int Index;
	Vedge *Next;
};

struct Pedge
{
	int Index;
	double Diagonal;
	double Offdiag;
	Pedge *Next;
};

struct PVedge
{
	int Index;
	PVedge *Next;
};

struct Psi
{
	int Index;
	double Value;
};

