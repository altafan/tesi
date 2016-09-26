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
	vector<int2> m_points;
	vector<int2> M_points;
	vector<int> dx;
	vector<int> dy;
}slab;

typedef struct maps_ {
	int ncols,nrows;
	int esx,esy;
	int2 min, max;
	int slabs_nrows;
	int dx,dy;
	F* btm_map;
	uchar4* host_info;
	int* host_info_x_m; 
	int* host_info_y_m;
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

#endif
