#pragma once
#include <stdio.h>
#include <iostream>
#include <math.h>

struct Vertex;
struct Patch;
struct Vedge;
struct Pedge;
struct PVedge;

struct Vertex
{
	int Index;
	int Degree;
	double Selfloop;
	double Wdegree;
	double totaldiagonal;
	double totaloffdiagonal;
	Vedge *Neighbour;
	Patch *Cluster;
	int Relativecor[2];
	double Eigfun;
	double U;
};

struct Patch
{
	int Index;
	int Mark;
	int Degree;
	int Psize;
	double Delta;
	double Lambda;
	
	PVedge *Member;
	PVedge *Lastmember;
	Pedge *Neighbour;
	Pedge *Lastneighbour;
	Patch *Cluster;

};

struct Vedge
{
	double Weight;
	double diagonal;
	double offdiagonal;
	Vertex *Invertex;
	Vertex *Outvertex;
	Vedge *Next;
};

struct Pedge
{
	double Bond;
	Patch *Inpatch;
	Patch *Outpatch;
	Pedge *Next;
};

struct PVedge
{
	double Eigfun;
	Vertex *Outvertex;
	PVedge *Next;
};