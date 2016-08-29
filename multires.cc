#include <iostream>
#include <fstream>
#include "types.h"
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

maps map;
global g;

void read_map(string path) {
	
	ifstream file(path);

	if(!file.fail()) {
		file >> map.nrows;
		file >> map.ncols;
		file >> map.minx_map;
		file >> map.miny_map;
		file >> map.dx;
		file >> map.dy; 		
	} else
		printf("Unable to open file. %s not found. \n\n", path.c_str());
		
	file.close();
		
	map.maxx_map = map.minx_map + (map.ncols - 1) * map.dx;
	map.maxy_map = map.miny_map + (map.nrows - 1) * map.dy;
		
	printf("Map info read:\n");
	printf("nrows: %d \nncols: %d \nminx: %d \nminy: %d \n",map.nrows,map.ncols,map.minx_map,map.miny_map);
	printf("maxx: %d \nmaxy: %d \ndx: %d \ndy: %d \n",map.maxx_map,map.maxy_map,map.dx,map.dy);

}

void read_pts(string path) {
	
	ifstream file(path);
	
	if(!file.fail()) {
		int i = 0;
		g.punti_m.resize(3);
		
		while(!file.eof()) {
			file >> g.punti_m[i].x;
			file >> g.punti_m[i].y;
			file >> g.punti_m[i].z;
			file >> g.punti_m[i].w;
			++i;
		}
	} else 
		printf("Unable to open file. %s not found. \n\n", path.c_str());
		
	file.close();
	
	printf("\nFile .pts read:\n");
	for(int i = 0; i < g.punti_m.size(); ++i)
		printf("%d° punto: %f, %f, %f, %f\n",(i+1),g.punti_m[i].x,g.punti_m[i].y,g.punti_m[i].z,g.punti_m[i].w);
	
}

void multires(/*char* path*/) {
	
	printf("\nMultiresolution processing ------------\n");
	
    int sx = map.ncols;
    int sy = map.nrows;
    int BLOCKSIZE_X = 16;
    int BLOCKSIZE_Y = 16;
    int levels = 4; //livelli di multirisoluzione
    int lblock = (1<<(levels - 1));
    int hi = 0;  
    
    // la mappa multires e' aumentata a multipli della risoluzione minima, 
    //in modo da avere quadrati divisi esattamente 
    int esx = ((((sx - 1) / BLOCKSIZE_X + 1) - 1) / lblock + 1) * lblock * BLOCKSIZE_X;
    int esy=((((sy - 1) / BLOCKSIZE_Y + 1) - 1)/ lblock + 1) * lblock * BLOCKSIZE_Y;
    map.esx = esx;
    map.esy = esy;

    int bsx = esx / BLOCKSIZE_X; // bitmap size (ridotta di blocksize)
    int bsy = esy / BLOCKSIZE_Y;


    //aggiorno range
    //ho tolto il bordo +2
    //ma allargato anche eventualmente la parte alta
	map.minx_map += map.dx;
	map.miny_map += map.dy;
	map.maxx_map = map.minx_map + map.dx * (esx - 1);
	map.maxy_map = map.miny_map + map.dy * (esy - 1);	
	
	printf("new maxx: %d\nnew maxy: %d\n", map.maxx_map,map.maxy_map);

    // bitmask: valori -1..levels-1, 
    //    -1 = quadrato copre completamente zona esterna
    // 
    char** bitmask = (char **) malloc(levels * sizeof(char*));
    int** bitmaskC = (int **) malloc(levels * sizeof(int*));
    uchar4** bitmaskS = (uchar4 **) malloc(levels * sizeof(uchar4*));//sides

    F* Z; // supporto per downsampling

    /*host_info_m
    map.host_info_m = (uchar4 **) malloc(levels * sizeof(uchar4*));//sides
    map.host_info_x_m = (int*) malloc(levels * sizeof(int));//sides
    map.host_info_y_m = (int*) malloc(levels * sizeof(int));//sides
    */

    for(int i = 0; i < levels; i++){
		int size = bsx * bsy / (1<<i) / (1<<i);
		printf("Level --> %d, size %d\n",i,size);
		bitmask[i] = (char *) malloc(size * sizeof(char));
		bitmaskC[i] = (int *) malloc(size * sizeof(int));
		bitmaskS[i] = (uchar4 *) malloc(size * sizeof(uchar4));
		for(int j = 0; j < size; j++) {
	  		bitmask[i][j] = 0;
	  		bitmaskC[i][j] = -1;
	  		bitmaskS[i][j].x = 0;	 
	  		bitmaskS[i][j].y = 0;	 
	  		bitmaskS[i][j].z = 0;	 
	  		bitmaskS[i][j].w = 0;	 
		} 

		int nx = map.ncols;
		int ny = map.nrows;
		int scale = (1<<i);
		nx = (nx - 1) / scale + 1;
		ny = (ny - 1) / scale + 1;
		
		if (i == 0) // caso livello 0, metto direttamente dentro
	  		Z = (F*) malloc(nx * ny * sizeof(F));
    }
    
    if(hi) // speciale: se tutto ad alta risoluzione
    for(int y = 1; y < sy + 1; y++){  
		for(int x = 1; x < sx + 1; x ++)
			if(map.btm_map[y * sx + x] < 1e10) {
				int bx = (x - 1) / BLOCKSIZE_X; // shifto uno indietro (c'e' bordo!)
			    int by = (y - 1) / BLOCKSIZE_Y;
			    bitmask[0][bsx*by+bx] = 1;
			}
    }
    	
	printf("Processing PTS %d points\n",g.punti_m.size());
    for(int i = 0; i < levels; i++) {
		int sizex = bsx / (1<<i);
		int sizey = bsy / (1<<i);

		// impone richieste di risoluzione
		if(hi == 0) // se non tutto ad alta
			for(int i1 = 0; i1 < g.punti_m.size(); i1++){
				if((int) g.punti_m[i1].z - 1 == i) { // lavora al livello corrente
					int x = ((g.punti_m[i1].x - map.minx_map) / map.dx - 1) / BLOCKSIZE_X; // scalo bordo maps
					int y = ((g.punti_m[i1].y - map.miny_map) / map.dy - 1) / BLOCKSIZE_Y;
					// ora allineo rispetto patch rispetto al livello di multirisoluzione richiesto
					x = x / (1<<i) * (1<<i);
					y = y / (1<<i) * (1<<i);

					printf("Force level %d point %f %f --> %d %d (%d %d)\n",(int)g.punti_m[i1].z,g.punti_m[i1].x,g.punti_m[i1].y,x,y,bsx,bsy);
					for(int ii = 0; ii < (1<<i); ii++)
						for(int jj = 0; jj < (1<<i); jj++)
							if(x + ii >= 0 && x + ii < bsx && y + jj >= 0 && y + jj < bsy)
								if(bitmask[0][bsx*(y+jj)+x+ii]==0 || bitmask[0][bsx*(y+jj)+x+ii]>g.punti_m[i1].z) // conquista la cella solo se non c'e' una risoluzione piu' alta (in caso la lascia stare!)
									bitmask[0][bsx * (y + jj) + x + ii] = g.punti_m[i1].z;
				}
			}
    }
	
	//multires a partire dai seed point
    for (int i = 0; i < levels; i++){
		int sizex = bsx / (1 << i);
		int sizey = bsy/(1 << i);
		
		for(int x = 0; x < sizex; x++)
	  		for(int y = 0; y < sizey; y++)
	    		if(bitmask[i][sizex * y + x] > 0) // se quadrato attivo
	      			if(i < levels - 1) { // aggiunge gli altri per completare il quadrato
						if(bitmask[i][sizex*((y/2)*2+0)+((x/2)*2+0)] == 0) 
							bitmask[i][sizex*((y/2)*2+0)+((x/2)*2+0)] = bitmask[i][sizex*y+x];
						if(bitmask[i][sizex*((y/2)*2+1)+((x/2)*2+0)] == 0) 
							bitmask[i][sizex*((y/2)*2+1)+((x/2)*2+0)] = bitmask[i][sizex*y+x];
						if(bitmask[i][sizex*((y/2)*2+0)+((x/2)*2+1)] == 0) 
							bitmask[i][sizex*((y/2)*2+0)+((x/2)*2+1)] = bitmask[i][sizex*y+x];
						if(bitmask[i][sizex*((y/2)*2+1)+((x/2)*2+1)] == 0) 
							bitmask[i][sizex*((y/2)*2+1)+((x/2)*2+1)] = bitmask[i][sizex*y+x];
					}	
	  			
		// copia nel livello successivo (tengo il valore piu' alto delle 4)
		if(i < levels - 1) {
	  		for(int x = 0; x < sizex / 2; x++)
	    		for(int y = 0; y < sizey / 2; y++){
		      		int max = -1;
		      		if(max < bitmask[i][sizex * (2*y+0) + (2*x+0)]) max = bitmask[i][sizex * (2*y+0) + (2*x+0)];
		      		if(max < bitmask[i][sizex * (2*y+1) + (2*x+0)]) max = bitmask[i][sizex * (2*y+1) + (2*x+0)];
		      		if(max < bitmask[i][sizex * (2*y+0) + (2*x+1)]) max = bitmask[i][sizex * (2*y+0) + (2*x+1)];
		      		if(max < bitmask[i][sizex * (2*y+1) + (2*x+1)]) max = bitmask[i][sizex * (2*y+1) + (2*x+1)];
		      		bitmask[i + 1][(sizex / 2) * (y) + (x)] = max;
		    	}

		  	//garantisce proprieta' di adiacenza quadrati con al massimo un livello di differenza
		  	for(int x = 0; x < sizex; x++)
		    	for(int y = 0; y < sizey; y++)
		      		if(bitmask[i][sizex * y + x] == i + 1) // se quadrato attivo
			    		for(int dx = -1; dx <= 1; dx++)
			  				for(int dy = -1; dy <= 1; dy++)
			    				if(dx != 0 || dy != 0) {
						      		int nx = x + dx;
						      		int ny = y + dy;
						      		
						      		if(nx < 0 || nx >= sizex || ny < 0 || ny >= sizey) {}
			      					else 
										if(bitmask[i][sizex * ny + nx] == 0) // manca un vicino ---> deve esistere ad un livello di multires piu' alto
				  							bitmask[i + 1][(sizex / 2) * (ny / 2) + (nx / 2)] = i + 2;
										else 
				  							if(bitmask[i][sizex * ny + nx] != -1 &&
				      						   abs(bitmask[i][sizex*ny+nx] - bitmask[i][sizex*y+x]) > 1) {
													printf("error: incompatible multiresolution requests level %d in %d %d (%d %dS) (lev %d) - %d %d (lev %d)\n",
														i+1,x,y, map.minx_map+map.dx*(1<<i)*x , map.miny_map+map.dy*(1<<i)*y,
														bitmask[i][sizex*y+x],
														nx,ny,
														bitmask[i][sizex*ny+nx]);
													exit(2);
				  							}			    
		      						
		    					}
		}
	}
	
	//sostituisce gli 0 con 4 (ris. più bassa)
	for(int i = levels - 1; i >= 0; i--) {
		int sizex = bsx / (1<<i);
		int sizey = bsy / (1<<i);

		for(int x = 0; x < sizex; x++)
			for(int y = 0; y < sizey; y++) {
				// assicura ultimo livello tutto calcolato
				if(i == levels - 1 && bitmask[i][sizex * y + x] == 0)
					bitmask[i][sizex * y + x] = i + 1;
					
				if(i > 0 && (bitmask[i][sizex*y+x] >= i+1 || bitmask[i][sizex*y+x] == -1)) { // se quadrato attivo o fuori
					bitmask[i - 1][sizex * 2 * (y*2+0) + (x*2+0)] = bitmask[i][sizex * y + x];
					bitmask[i - 1][sizex * 2 * (y*2+0) + (x*2+1)] = bitmask[i][sizex * y + x];
					bitmask[i - 1][sizex * 2 * (y*2+1) + (x*2+0)] = bitmask[i][sizex * y + x];
					bitmask[i - 1][sizex * 2 * (y*2+1) + (x*2+1)] = bitmask[i][sizex * y + x];
				}
			}
    }
	
	if(0)
		//stampa bitmask
	    for(int i1 = 0; i1 < levels; i1++) {
			int sizex = bsx / (1<<i1);
			int sizey = bsy / (1<<i1);
			
			for(int y = sizey-1; y >=0; y--) {
				printf("%3d: ", y);
				for(int x = 0; x < sizex; x++)
		    		if(bitmask[i1][sizex * y + x] == -1)
		      			printf("! ");
		    		else 
						printf("%d ",bitmask[i1][sizex*y+x]);
		  		printf("\n");
			}
	    }
    
    //codifica dei blocchi
	int tot_blocks = 0;
	int bound_blocks = 0;
	int codes = 0;
	for (int i = 0; i < levels; i++){
		int dbg = 0;
		int sizex = bsx / (1<<i);
		int sizey = bsy / (1<<i);
		int ct = 0;
		for (int y = sizey-1; y >=0; y--){  
			if(dbg) printf("%d: ",y);
			for (int x = 0; x < sizex; x++){
				ct += (bitmask[i][sizex*y+x] == i+1 ? 1 : 0);
				if(bitmask[i][sizex*y+x] != i+1){
					if(dbg) printf("%d ",bitmask[i][sizex*y+x]);//<=i+1?bitmask[i][sizex*y+x]:0));
				} else {
					if(dbg) printf("%dB",bitmask[i][sizex*y+x]);
					bitmaskC[i][sizex*y+x]=codes++;
					bound_blocks++;
				}
					
			}
			if(dbg) printf("\n");
		}
		tot_blocks += ct;
		printf("Level %d: ratio %f (%d cells)\n",i,(ct+0.0)/(sizex*sizey),ct);	  
	}  
	
	printf("Celle originali %d, celle multires %d\n",map.ncols*map.nrows,tot_blocks*BLOCKSIZE_X*BLOCKSIZE_Y);
    printf("Compressione: %2.1fx\n",1.0/((0.0+tot_blocks)*BLOCKSIZE_X*BLOCKSIZE_Y/map.ncols/map.nrows));

	for(int i = 0; i < levels; i++){
		int dbg = 0;
		int sizex = bsx / (1<<i);
		int sizey = bsy / (1<<i);

		for(int y = sizey - 1; y >= 0; y--){  
			for(int x = 0; x < sizex; x++){
				if (bitmask[i][sizex*y+x]==i+1 &&
					bitmaskC[i][sizex*y+x]==-1) // se non gia' assegnato (contorni)
						bitmaskC[i][sizex*y+x] = codes++;
			
				if(dbg) {
					if (bitmaskC[i][sizex*y+x]==-1) 
						printf("... ");
					else
						printf("%3d ",bitmaskC[i][sizex*y+x]);
		    	}
		  	}
		  	if(dbg) printf("\n");
		}
	}
	
	printf("Celle %d, bound %d\n",tot_blocks,bound_blocks);
	
	map.tot_blocks = tot_blocks;
	map.bound_blocks = bound_blocks;
}


int main() {
	
	string map_file = "map_info.txt";
	string pts_file = "polygons/prova.PTS";

	read_map(map_file);
	read_pts(pts_file);
	multires();

	return 0;

}
