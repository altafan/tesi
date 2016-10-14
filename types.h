#pragma once
#include <string>
#include <vector>

using namespace std;

#ifndef types_h
#define types_h

#define BIT_EXTERN 15
#define BIT_TAPPO 14
#define BIT_N 1
#define BIT_S 2
#define BIT_W 4
#define BIT_E 8
#define ZERO 55
#define IN 21
#define OUT 39

typedef float F;

typedef struct F4_ {
	float x,y,z,w;
} F4;

typedef struct uchar4_ {
	unsigned char x,y,z,w;
} uchar4;

typedef struct point_ {
  double x,y;
} point;

typedef struct int2_ {
	int x,y;
} int2;

typedef struct int4_ {
	int x,y,z,w;
} int4;

typedef struct ushort2_ {
	unsigned short x,y;
} ushort2;

typedef struct {
	char lev;
  	int n1;
  	int n2;
} neigh_t;

typedef struct polygon_ {
	vector<point> points;
  	vector<int> edges;    
} polygon;

typedef struct slab_ {
	point m_points;
	point M_points;
	int dx;
	int dy;
}slab;

typedef struct maps_ {
	int ncols,nrows;
	int esx,esy;
	int2 min, max;
	int slabs_nrows;
	int dx,dy;
	F* btm_map;
	uchar4* host_info;
	int tot_blocks, bound_blocks;
	int first_slab, last_slab;
	int dxs, dys;
	F4* host_grid_multi;
	unsigned char* host_grid_level_multi;
	ushort2* host_ofs_blocks;
	neigh_t* neigh;
} maps;

typedef struct global_ {
	vector<F4> punti_m;
} global;

typedef struct bc_ {
	int ifl_bc_w, numbc_w;
  	int ifl_bc_e, numbc_e;
  	int ifl_bc_n, numbc_n;
  	int ifl_bc_s, numbc_s;
} bc;

#endif
