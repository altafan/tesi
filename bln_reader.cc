#include <iostream>
#include <fstream>
#include "types.h"
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

polygon pol;
vector<slab> slabs;

point pol_mp,pol_Mp;
int2 min_xy;
int2 s_max;
int2 s_min;
int npoints;
int dbg = 1;

void read_bln(string path) {

	string pathBLN = path+".BLN";

	ifstream file(pathBLN);

	// legge poligono da file
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
	} else{
		printf("Unable to open file. %s not found.\n\n", pathBLN.c_str());
		exit(1);
	}

	file.close();

	pol_mp.x = pol.points[0].x;
	pol_mp.y = pol.points[0].y;
  	pol_Mp.x = pol_mp.x;
	pol_Mp.y = pol_mp.y;

	// trovo i corners del bounding box (sx-dn, dx-up)
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

	if(dbg) printf("Polygon read\n\n");

}

void bln_interpolation() {

	string path = "polygons/bln_raster.PTS";
	FILE * file;
	file = fopen(path.c_str(),"w");

	// per ogni coppia di vertici interpolazione del segmento che li unisce
	for(int i = 0; i < npoints; ++i) {
		double x0,x1,y0,y1;

		x0 = pol.points[i].x;
		y0 = pol.points[i].y;
		x1 = pol.points[(i+1) % npoints].x;
		y1 = pol.points[(i+1) % npoints].y;
		
		double deltax = x1 - x0;
		double deltay = y1 - y0;
		double dd = sqrt(deltax * deltax + deltay * deltay);
		double tx = x0;
		double ty = y0;
		double l = 0;

		while(l <= dd) {
			tx += deltax / dd;
			ty += deltay / dd;

			fprintf(file,"%f %f 3 0\n",tx,ty/*,pol.edges[i]*/);

			l = sqrt((tx - x0) * (tx - x0) + (ty - y0) * (ty - y0));
		}

	}

	fclose(file);

	if(dbg) printf("\nPolygon interpolation complete\n");

}

void read_grd(string path, int n_slab) {

	string pathGRD = path;

	// dimensiona i vettori contenti le info delle tavolette
	slabs.resize(n_slab);

	s_max.x = 0;
	s_max.y = 0;

	s_min.x = 9999999;
	s_min.y = 9999999;

	// legge le info delle tavolette
	for(int i = 0; i < n_slab; ++i) {
		pathGRD = "";
		pathGRD = path + to_string(i+1) + ".grd";		

		ifstream file(pathGRD);

		char type[10];

		if(!file.fail()) {
			file >> type;
			file >> slabs[i].dx >> slabs[i].dy;
			file >> slabs[i].m_points.x >> slabs[i].M_points.x;
			file >> slabs[i].m_points.y >> slabs[i].M_points.y;

			// trova coordinate della tavoletta con max(x,y)
			if(slabs[i].M_points.x > s_max.x)
				s_max.x = slabs[i].M_points.x;
			if(slabs[i].M_points.y > s_max.y)
				s_max.y = slabs[i].M_points.y;

			// trova coordinate della tavoletta con min(x,y)
			if(slabs[i].m_points.x < s_min.x)
				s_min.x = slabs[i].m_points.x;
			if(slabs[i].m_points.y < s_min.y)
				s_min.y = slabs[i].m_points.y;

		} else {
			printf("Unable to open file. %s not found.\n\n", pathGRD.c_str());
			exit(1);
		}

		file.close();
	}
	
	// calcola il vertice in basso a sx della tavoletta più a sx interna al bounding box
	min_xy.x = floor((((int)(pol_mp.x - s_min.x) / slabs[0].dx) * slabs[0].dx + s_min.x) / 10) * 10;
	min_xy.y = floor((((int)(pol_mp.y - s_min.y) / slabs[0].dy) * slabs[0].dy + s_min.y) / 10) * 10;
	
	if(dbg) {
		printf("Slabs read:\n");
		printf("s_min: %d, %d\n", s_min.x, s_min.y);
		printf("s_max: %d, %d\n", s_max.x, s_max.y);
		printf("min_xy: %d, %d\n", min_xy.x, min_xy.y);
		printf("dx dy: %d %d\n\n", slabs[0].dx, slabs[0].dy);
	}

}

int count_grd(string path) {

	// esegue script che conta il numero di tavolette in slabs/
	int status = system("./createGrid.sh");
	int n;
	if(status == 0) {
		// legge il numero di tavolette
		ifstream file("grd.info");	
		file >> n;
		file.close();
	} else {
		printf("Unable to find .grd files in %s\n", path.c_str());
		exit(1);
	}

	return n;

}

void bounding_box(int n_slab) {

	// calcola numero di righe e di colonne della griglia di tavolette
	int cols = (s_max.x - s_min.x) / slabs[0].dx + 1;
	int rows = (s_max.y - s_min.y) / slabs[0].dy + 1;

	int m[rows][cols];

	// calcola gli indici della griglia di tavolette in cui cade il bounding box
	int pmi = rows - 1 - ((int)(pol_Mp.y - s_min.y)) / slabs[0].dy;
	int pmj = (pol_mp.x - s_min.x) / slabs[0].dx;
	int pMi = rows - ((pol_mp.y - s_min.y) / slabs[0].dy) ;
	int pMj = ((int)(pol_Mp.x - s_min.x)) / slabs[0].dx;	

	// per ogni tavoletta calcola gli indici (la posizione nella griglia) e controlla
	// se sta nell'intervallo degli indici del boundig box
    for(int k = 0; k < n_slab; ++k) {
    	int j = ceil((slabs[k].m_points.x - s_min.x) / (float)slabs[k].dx);
    	int i = rows - 1 - ceil((slabs[k].m_points.y - s_min.y) / (float)slabs[k].dy);

    	if(pmi <= i && i <= pMi && pmj <= j && j <= pMj) 
    		m[i][j] = k+1;
    	else 	
    		m[i][j] = 0;
    }
    
    // salva le info della matrice di tavolette interne al bounding box:
    // nrow ncol, min_xy, dx dy
    ofstream file;
    file.open("map_info.txt");
    file << rows << "\n";
    file << (pMi - pmi + 1) * (slabs[0].dy - 1) << " " << (pMj - pmj + 1) * (slabs[0].dx - 1)<< "\n";
    file << min_xy.x << " " << min_xy.y << "\n";
    file << 1 << " " << 1 << "\n";
    file << cols * pmj + (rows - pMi) << " " << cols * pMj + (rows - pmi)<< "\n";
    file << slabs[0].dx << " " << slabs[0].dy << "\n";
    file.close();

    if(dbg) { 
		printf("Slab matrix - rows cols: %d %d\n", rows, cols);
   		printf("Bounding box min  x,y: %d %d\n", pMi, pmj);
		printf("Bounding box max x,y: %d %d\n\n", pmi, pMj);
	    for(int i = 0; i < rows; ++i){
	    	for(int j = 0; j < cols; ++j)
	    		if(m[i][j])
	    			printf("x ");
	    		else 
	    			printf("0 ");
	    	printf("\n");
	    }
	}

}

int main() {

	string BLN_path = "polygons/bln2";
	string GRD_path = "slabs/";

	read_bln(BLN_path);

	int n_slab = count_grd(GRD_path);

	read_grd(GRD_path, n_slab);
	
	bounding_box(n_slab);

	bln_interpolation();

	return 0;

}
