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
  float x,y;
} point;

typedef struct int_point_ {
	int x,y;
} int_point;

typedef struct polygon_ {
  vector<point> points;
  vector<int> edges;    
} polygon;

typedef struct slab_ {
	vector<int_point> m_points;
	vector<int_point> M_points;
	vector<int> dx;
	vector<int> dy;
}slab;

typedef struct maps_ {
	int ncols,nrows;
	int esx,esy;
	int minx_map,miny_map;
	int maxx_map, maxy_map;
	int dx,dy;
	F* btm_map;
	uchar4** host_info_m;
	int* host_info_x_m; 
	int* host_info_y_m;
	int tot_blocks, bound_blocks;
	int first, last;
	int dxs, dys;
	F4* host_grid_multi;
} maps;

typedef struct global_ {
	vector<F4> punti_m;
} global;

#endif
