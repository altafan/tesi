//#include "bln_reader.h"

#include <iostream>
#include <fstream>
#include "types.h"
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

polygon pol;
slab slabs;

point pol_mp,pol_Mp;
int_point min_xy;
int_point s_max;
int_point s_min;

void readBLN(string path) {

	string pathBLN = path+".BLN";
	int npoints;

	/*strcpy(pathBLN, path);
	strcat(pathBLN, ".BLN");*/

	ifstream file(pathBLN);

	//legge poligono da file
	if(!file.fail()) {
		file >> npoints;
		pol.points.resize(npoints);
		pol.edges.resize(npoints);
		int i = 0;

		while(!file.eof()) {
			file >> pol.points[i].x;
			file >> pol.points[i].y;
			file >> pol.edges[i];
			++i;
		}
	} else
			printf("Unable to open file. %s not found.\n\n", pathBLN.c_str());

	file.close();

	pol_mp.x = pol.points[0].x;
	pol_mp.y = pol.points[0].y;
  	pol_Mp.x = pol_mp.x;
	pol_Mp.y = pol_mp.y;

	//trovo i corner del bounding box (sx-dn, dx-up)
	for(int i = 1; i < npoints; ++i) {
		if(pol.points[i].x < pol_mp.x)
			pol_mp.x = (int)pol.points[i].x;
		if(pol.points[i].y < pol_mp.y)
			pol_mp.y = (int)pol.points[i].y;
		if(pol.points[i].x > pol_Mp.x)
			pol_Mp.x = (int)pol.points[i].x;
		if(pol.points[i].y > pol_Mp.y)
			pol_Mp.y = (int)pol.points[i].y;
	}

	printf("Bounding box corners:\n");
	printf("v_min: (%f,%f) \n",pol_mp.x,pol_mp.y);
	printf("v_max: (%f,%f) \n",pol_Mp.x,pol_Mp.y);
	printf("polygon read\n");

}

void readGRD(string path, int n_slab) {

	string pathGRD = path;

	//dimensiona i vettori contenti le info delle tavolette
	slabs.m_points.resize(n_slab);
	slabs.M_points.resize(n_slab);
	slabs.dx.resize(n_slab);
	slabs.dy.resize(n_slab);

	s_max.x = 0;
	s_max.y = 0;

	s_min.x = 9999999;
	s_min.y = 9999999;

	//legge le info delle tavolette
	for(int i = 0; i < n_slab; ++i) {
		pathGRD = "";
		pathGRD = path + to_string(i+1) + ".grd";		

		ifstream file(pathGRD);

		char type[10];

		if(!file.fail()) {
			file >> type;
			file >> slabs.dx[i] >> slabs.dy[i];
			file >> slabs.m_points[i].x >> slabs.M_points[i].x;
			file >> slabs.m_points[i].y >> slabs.M_points[i].y;

			//trova coordinate della tavoletta con max(x,y)
			if(slabs.M_points[i].x > s_max.x)
				s_max.x = slabs.M_points[i].x;
			if(slabs.M_points[i].y > s_max.y)
				s_max.y = slabs.M_points[i].y;

			//trova coordinate della tavoletta con min(x,y)
			if(slabs.m_points[i].x < s_min.x)
				s_min.x = slabs.m_points[i].x;
			if(slabs.m_points[i].y < s_min.y)
				s_min.y = slabs.m_points[i].y;

		} else
			printf("Unable to open file. %s not found.\n\n", pathGRD.c_str());

		file.close();
	}
	
	//calcola il vertice in basso a sx della tavoletta più a sx interna al bounding box
	min_xy.x = ((int)(pol_mp.x - s_min.x) / slabs.dx[0]) * 
		slabs.dx[0] + s_min.x;
	min_xy.y = ((int)(pol_mp.y - s_min.y) / slabs.dy[0]) * 
		slabs.dy[0] + s_min.y;
	
	printf("s_min: %d, %d\n", s_min.x, s_min.y);
	printf("s_max: %d, %d\n", s_max.x, s_max.y);
	printf("min_xy: %d, %d\n", min_xy.x, min_xy.y);
	printf("dx dy: %d %d\n", slabs.dx[0], slabs.dy[0]);
	printf("slabs read\n");

}

int countGRD(string path) {

	//esegue script che conta il numero di tavolette in slabs/
	int status = system("./createGrid.sh");
	int n;
	if(status == 0) {
		//legge il numero di tavolette
		ifstream file("grd.info");	
		file >> n;
		file.close();
	} else
		printf("Unable to find .grd files in %s\n", path.c_str());

	printf("slab info read\n");

	return n;

}

void raster(int n_slab) {

	//calcola numero di righe e di colonne della griglia di tavolette
	int cols = (s_max.x - s_min.x) / slabs.dx[0] + 1;
	int rows = (s_max.y - s_min.y) / slabs.dy[0] + 1;

	printf("rows cols: %d %d\n", rows, cols);

	int m[rows][cols];

	//calcola gli indici della griglia di tavolette in cui cade il bounding box
	int pmi = rows - 1 - ((int)(pol_Mp.y - s_min.y)) / slabs.dy[0];
	int pmj = (min_xy.x - s_min.x) / slabs.dx[0];
	int pMi = rows-1-((min_xy.y - s_min.y)/slabs.dy[0]);
	int pMj = ((int)(pol_Mp.x - s_min.x)) / slabs.dx[0];	
	
	printf("pMi pmj: %d %d\n", pMi, pmj);
	printf("pmi pMj: %d %d\n", pmi, pMj);

	//per ogni tavoletta calcola gli indici (la posizione nella griglia) e controlla
	//se sta nell'intervallo degli indici del boundig box
    for(int k = 0; k < n_slab; ++k) {
    	int j = ceil((slabs.m_points[k].x - s_min.x) / (float)slabs.dx[k]);
    	int i = rows - 1 - ceil((slabs.m_points[k].y - s_min.y) / (float)slabs.dy[k]);

    	if(pmi <= i && i <= pMi && pmj <= j && j <= pMj) 
    		m[i][j] = k+1;
    	else 	
    		m[i][j] = 0;
    }

    printf("\nrasterizzation:\n");
    for(int i = 0; i < rows; ++i){
    	for(int j = 0; j < cols; ++j)
    		if(m[i][j])
    			printf("x ");
    		else 
    			printf("0 ");
    	printf("\n");
    }	

}

int main() {

	string BLN_path = "polygons/bln2";
	string GRD_path = "slabs/";

	readBLN(BLN_path);

	int n_slab = countGRD(GRD_path);

	readGRD(GRD_path, n_slab);
	
	raster(n_slab);

	return 0;

}
