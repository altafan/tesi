#pragma once
#include <vector>

using namespace std;

#ifndef types_h
#define types_h

typedef struct point_ {
  float x,y;
} point;

typedef struct int_point_{
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

#endif
