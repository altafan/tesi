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
int2 min_xy;
int2 s_max;
int2 s_min;
int npoints;

void readBLN(string path) {

	string pathBLN = path+".BLN";

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
	} else{
		printf("Unable to open file. %s not found.\n\n", pathBLN.c_str());
		exit(1);
	}

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

	/*printf("Bounding box corners:\n");
	printf("v_min: (%d,%d) \n",(int)pol_mp.x,(int)pol_mp.y);
	printf("v_max: (%d,%d) \n",(int)pol_Mp.x,(int)pol_Mp.y);*/
	printf("polygon read\n");

}

void blnInterpolation() {

	string path = "polygons/bln_raster.PTS";
	FILE * file;
	file = fopen(path.c_str(),"w");

	//per ogni coppia di vertici interpolazione del segmento che li unisce
	for(int i = 0; i < npoints; ++i) {
	//int i=0; {
		double x0,x1,y0,y1,m;
		if(i == npoints - 1) {
			x0 = pol.points[i].x;
			y0 = pol.points[i].y;
			x1 = pol.points[0].x;
			y1 = pol.points[0].y;
		} else {
			x0 = pol.points[i].x;
			y0 = pol.points[i].y;
			x1 = pol.points[i+1].x;
			y1 = pol.points[i+1].y;
		}
		/*m = (y1-y0) / (x1-x0);

		if(x1 > x0)
			for(int x = x0; x < x1; x += 4) {
				int y = (int)(m * (x - x0) + y0);
				//printf("%d %d\n",x,y);
				//file << x << " " << y << " 1 0\n";
				printf("%d %d\n",x,y);
			}
		else
			for(int x = x1; x < x0; x += 4) {
				int y = (int)(m * (x - x1) + y1);
				//file << x << " " << y << " 1 0\n";
				printf("%d %d\n",x,y);	
			}*/
		double deltax = x1 - x0;
		double deltay = y1 - y0;
		double dd = sqrt(deltax * deltax + deltay * deltay);
		double tx = x0;
		double ty = y0;
		double l = 0;
		while(l <= dd) {
			tx += deltax / dd;
			ty += deltay / dd;
			//file << (int)tx << " " << (int)ty << " 3 0\n";
			fprintf(file,"%f %f 3 0\n",tx,ty);
			l = sqrt((tx - x0) * (tx - x0) + (ty - y0) * (ty - y0));
		}

	}

	fclose(file);

	printf("Polygon interpolation: complete\n");

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

		} else {
			printf("Unable to open file. %s not found.\n\n", pathGRD.c_str());
			exit(1);
		}

		file.close();
	}
	
	//calcola il vertice in basso a sx della tavoletta più a sx interna al bounding box
	min_xy.x = //(int)((pol_mp.x - s_min.x) / slabs.dx[0]) * slabs.dx[0] + s_min.x;
	floor((((int)(pol_mp.x - s_min.x) / slabs.dx[0]) * slabs.dx[0] + s_min.x) / 10) * 10;
	min_xy.y = //(int)((pol_mp.y - s_min.y) / slabs.dy[0]) * slabs.dy[0] + s_min.y;
	floor((((int)(pol_mp.y - s_min.y) / slabs.dy[0]) * slabs.dy[0] + s_min.y) / 10) * 10;
	
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
	} else {
		printf("Unable to find .grd files in %s\n", path.c_str());
		exit(1);
	}

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
	int pmj = (pol_mp.x - s_min.x) / slabs.dx[0];
	int pMi = rows - ((pol_mp.y - s_min.y)/slabs.dy[0]) ;
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
    
    // salva le info della matrice di tavolette interne al bounding box:
    // nrow ncol, min_xy, dx dy
    ofstream file;
    file.open("map_info.txt");
    file << rows << "\n";
    file << (pMi - pmi + 1) * slabs.dy[0] << " " << (pMj - pmj + 1) * slabs.dx[0] << "\n";
    file << min_xy.x << " " << min_xy.y << "\n";
    file << 1 << " " << 1 << "\n";
    file << cols * pmj + (rows - pMi) << " " << cols * pMj + (rows - pmi)<< "\n";
    file << slabs.dx[0] << " " << slabs.dy[0] << "\n";
    file.close();

}

int main() {

	string BLN_path = "polygons/bln2";
	string GRD_path = "slabs/";

	readBLN(BLN_path);
	blnInterpolation();

	int n_slab = countGRD(GRD_path);

	readGRD(GRD_path, n_slab);
	
	raster(n_slab);

	return 0;

}
