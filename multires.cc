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
polygon pol;
int dbg = 0;

void read_map(string path) {
	
	ifstream file(path);

	if(!file.fail()) {
		file >> map.slabs_nrows;
		file >> map.nrows;
		file >> map.ncols;
		file >> map.min.x;
		file >> map.min.y;
		file >> map.dx;
		file >> map.dy; 
		file >> map.first_slab;
		file >> map.last_slab;
		file >> map.dxs;
		file >> map.dys;		
	} else {
		printf("Unable to open file. %s not found. \n\n", path.c_str());
		exit(1);
	}
		
	file.close();
		
	map.max.x = floor((map.min.x + (map.ncols - 1) * map.dx) / 10) * 10;
	map.max.y = floor((map.min.y + (map.nrows - 1) * map.dy) / 10) * 10;
	
	if(dbg) {
		printf("Bounding box:\nrows cols:%d %d\nmin x,y: %d %d\n",map.nrows,map.ncols,map.min.x,map.min.y);
		printf("max x,y: %d %d \ndx dy: %d %d\n",map.max.x,map.max.y,map.dx,map.dy);
		printf("First and last slabs: %d %d\nSlab dx dy: %d %d\n", map.first_slab, map.last_slab, map.dxs,map.dys);
	}

}

void read_pts(string path) {
	
	ifstream file(path);
	
	if(!file.fail()) {
		int i = 0;
		g.punti_m.resize(5755);
		
		while(!file.eof()) {
			file >> g.punti_m[i].x;
			file >> g.punti_m[i].y;
			file >> g.punti_m[i].z;
			file >> g.punti_m[i].w;
			++i;
		}
	} else {
		printf("Unable to open file. %s not found. \n\n", path.c_str());
		exit(1);
	}
		
	file.close();
	
	if(dbg) printf("\nFile %s read\n\n",path.c_str());

}

void readBLN(string path) {

	string pathBLN = path+".BLN";
	int npoints;

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

	if(dbg) printf("Polygon read\n\n");

}

void multires(string path) {
	
	if(dbg) printf("Multi resolution processing:\n\n");
	
    int sx = map.ncols;
    int sy = map.nrows;
    int BLOCKSIZE_X = 16;
    int BLOCKSIZE_Y = 16;
    int levels = 4; //livelli di multirisoluzione
    int lblock = (1 << (levels - 1));
    int hi = 0;  
    
    // la mappa multires e' aumentata a multipli della risoluzione minima, 
    //in modo da avere quadrati divisi esattamente 
    int esx = ((((sx - 1) / BLOCKSIZE_X + 1) - 1) / lblock + 1) * lblock * BLOCKSIZE_X;
    int esy = ((((sy - 1) / BLOCKSIZE_Y + 1) - 1) / lblock + 1) * lblock * BLOCKSIZE_Y;
    map.esx = esx;
    map.esy = esy;

    int bsx = esx / BLOCKSIZE_X; // bitmap size (ridotta di blocksize)
    int bsy = esy / BLOCKSIZE_Y;

	map.max.x = map.min.x + map.dx * (esx);
	map.max.y = map.min.y + map.dy * (esy);	
	
	if(dbg) printf("New max x,y: %d %d\n\n", map.max.x,map.max.y);

    // bitmask: valori -1..levels-1, 
    //    -1 = quadrato copre completamente zona esterna
    // 
    char** bitmask = (char **) malloc(levels * sizeof(char*));
    int** bitmaskC = (int **) malloc(levels * sizeof(int*));
    uchar4** bitmaskS = (uchar4 **) malloc(levels * sizeof(uchar4*));//sides

    F* Z; // supporto per downsampling

    for(int i = 0; i < levels; i++){
		int size = bsx * bsy / (1<<i) / (1<<i);

		if (dbg) printf("Level --> %d, size %d\n",i,size);

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

		int scale = (1 << i);

		nx = (nx - 1) / scale + 1;
		ny = (ny - 1) / scale + 1;
		
		if (i == 0) // caso livello 0, metto direttamente dentro
	  		Z = (F*) malloc(nx * ny * sizeof(F));
    }
    
    if(hi) // speciale: se tutto ad alta risoluzione
	    for(int y = 1; y < sy + 1; y++){  
			for(int x = 1; x < sx + 1; x++)
				if(map.btm_map[y * sx + x] < 1e10) {
					int bx = (x - 1) / BLOCKSIZE_X; // shifto uno indietro (c'e' bordo!)
				    int by = (y - 1) / BLOCKSIZE_Y;

				    bitmask[0][bsx*by+bx] = 1;
				}
	    }
    	
	if(dbg) printf("\nProcessing %ld seed points:\n\n",g.punti_m.size());

    for(int i = 0; i < levels; i++) {
		int sizex = bsx / (1<<i);
		int sizey = bsy / (1<<i);

		// impone richieste di risoluzione
		if(hi == 0) // se non tutto ad alta
			for(int i1 = 0; i1 < g.punti_m.size(); i1++){
				if((int) g.punti_m[i1].z - 1 == i) { // lavora al livello corrente
					int x = ((g.punti_m[i1].x - map.min.x) / map.dx) / BLOCKSIZE_X; 
					int y = ((g.punti_m[i1].y - map.min.y) / map.dy) / BLOCKSIZE_Y;

					// ora allineo rispetto patch rispetto al livello di multirisoluzione richiesto
					x = x / (1<<i) * (1<<i);
					y = y / (1<<i) * (1<<i);

					//printf("Force level %d point %f %f --> %d %d (%d %d)\n",(int)g.punti_m[i1].z,g.punti_m[i1].x,g.punti_m[i1].y,x,y,bsx,bsy);
					for(int ii = 0; ii < (1<<i); ii++)
						for(int jj = 0; jj < (1<<i); jj++)
							if(x + ii >= 0 && x + ii < bsx && y + jj >= 0 && y + jj < bsy)
								//conquista la cella solo se non c'e' una risoluzione piu' alta (in caso la lascia stare!)
								if(bitmask[0][bsx*(y+jj)+x+ii] == 0 || bitmask[0][bsx*(y+jj)+x+ii] > g.punti_m[i1].z) 
									bitmask[0][bsx * (y + jj) + x + ii] = g.punti_m[i1].z;
				}
			}
    }
	
	//multires a partire dai seed point
    for(int i = 0; i < levels; i++) {
		int sizex = bsx / (1 << i);
		int sizey = bsy/(1 << i);
		
		for(int x = 0; x < sizex; x++)
	  		for(int y = 0; y < sizey; y++)
	    		if(bitmask[i][sizex * y + x] > 0) // se quadrato attivo
	      			if(i < levels - 1) { // aggiunge gli altri per completare il quadrato
						if(bitmask[i][sizex * ((y/2)*2+0) + ((x/2)*2+0)] == 0) 
							bitmask[i][sizex * ((y/2)*2+0) + ((x/2)*2+0)] = bitmask[i][sizex*y+x];
						if(bitmask[i][sizex * ((y/2)*2+1) + ((x/2)*2+0)] == 0) 
							bitmask[i][sizex * ((y/2)*2+1) + ((x/2)*2+0)] = bitmask[i][sizex*y+x];
						if(bitmask[i][sizex * ((y/2)*2+0) + ((x/2)*2+1)] == 0) 
							bitmask[i][sizex * ((y/2)*2+0) + ((x/2)*2+1)] = bitmask[i][sizex*y+x];
						if(bitmask[i][sizex * ((y/2)*2+1) + ((x/2)*2+1)] == 0) 
							bitmask[i][sizex * ((y/2)*2+1) + ((x/2)*2+1)] = bitmask[i][sizex*y+x];
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
		      		bitmask[i + 1][(sizex / 2) * y + x] = max;
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
													printf("Error: incompatible multiresolution requests level %d in %d %d (%d %dS) (lev %d) - %d %d (lev %d)\n",
														i+1,x,y, map.min.x+map.dx*(1<<i)*x , map.min.y+map.dy*(1<<i)*y,
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
					
				if(i > 0 && (bitmask[i][sizex*y+x] >= i + 1 || bitmask[i][sizex*y+x] == -1)) { // se quadrato attivo o fuori
					bitmask[i - 1][sizex * 2 * (y*2+0) + (x*2+0)] = bitmask[i][sizex * y + x];
					bitmask[i - 1][sizex * 2 * (y*2+0) + (x*2+1)] = bitmask[i][sizex * y + x];
					bitmask[i - 1][sizex * 2 * (y*2+1) + (x*2+0)] = bitmask[i][sizex * y + x];
					bitmask[i - 1][sizex * 2 * (y*2+1) + (x*2+1)] = bitmask[i][sizex * y + x];
				}
			}
    }
  
    //codifica dei blocchi
	int tot_blocks = 0;
	int bound_blocks = 0;
	int codes = 0;

	for (int i = 0; i < levels; i++){
		int sizex = bsx / (1 << i);
		int sizey = bsy / (1 << i);
		int ct = 0;
		for (int x = 0; x < sizex; x++){  
			for (int y = 0; y < sizey; y++){
				ct += (bitmask[i][sizex*y+x] == i + 1 ? 1 : 0);

				if(bitmask[i][sizex * y + x] == i + 1){
					bitmaskC[i][sizex * y + x] = codes++;
					//bound_blocks++;
				}
			}
		}
		tot_blocks += ct;
		if(dbg) printf("Level %d: ratio %f (%d cells)\n",i,(ct+0.0)/(sizex*sizey),ct);	  
	}  
	
	if(dbg) {
		printf("\nCelle originali %d, celle multires %d\n",map.ncols*map.nrows,tot_blocks*BLOCKSIZE_X*BLOCKSIZE_Y);
	    printf("Compressione: %2.1fx\n\n",1.0/((0.0+tot_blocks)*BLOCKSIZE_X*BLOCKSIZE_Y/map.ncols/map.nrows));

	    //stampa bitmask e bitmaskC
	    for(int i1 = 0; i1 < levels; i1++) {
			int sizex = bsx / (1 << i1);
			int sizey = bsy / (1 << i1);
			
			for(int y = sizey - 1; y >= 0; y--) {
				printf("%3d: ", y);
				for(int x = 0; x < sizex; x++)
		    		if(bitmask[i1][sizex * y + x] == -1)
		      			printf("! ");
		    		else 
		    			printf("%d ",bitmask[i1][sizex*y+x]);
				printf("\n");
			}
    	}
}
    	for(int i = 0; i < levels; i++){
			int sizex = bsx / (1<<i);
			int sizey = bsy / (1<<i);

			for(int y = sizey - 1; y >= 0; y--){ 
				printf("%d: ",y); 
				for(int x = 0; x < sizex; x++){
					if (bitmaskC[i][sizex * y + x] == -1) 
						printf("... ");
					else
						printf("%3d ",bitmaskC[i][sizex * y + x]);
			  	}
			  	printf("\n");
			}
		}

		printf("\nCelle %d, bound %d\n\n",tot_blocks,bound_blocks);
  	//}
	
	map.tot_blocks = tot_blocks;
	map.bound_blocks = bound_blocks;

	int x_blocks = 1 << (int)(ceil(log((float)map.tot_blocks) / log((float)2) / (float)2)); 
  	int y_blocks = (tot_blocks - 1) / x_blocks + 1;

  	if(dbg) {
  		printf("------------------------------------------------------\n\n");
		printf ("Blocks %d alloc: %d x %d array \n",tot_blocks,x_blocks,y_blocks);
	}

	map.host_grid_level_multi = (unsigned char*) malloc(tot_blocks * sizeof(unsigned char));
    map.host_ofs_blocks = (ushort2*) malloc(tot_blocks * sizeof(ushort2));
	map.host_grid_multi = (F4*) malloc(x_blocks * y_blocks * BLOCKSIZE_X * BLOCKSIZE_Y * sizeof(F4));
	map.host_info = (uchar4*) malloc(map.nrows * map.ncols * sizeof(uchar4));
	map.neigh = (neigh_t*) malloc(8 * tot_blocks * sizeof(neigh_t));
	int* counter = (int*) malloc(x_blocks * y_blocks * BLOCKSIZE_X * BLOCKSIZE_Y * sizeof(int));
	int* assegnato = (int*) malloc(tot_blocks * sizeof(int*));
	vector<point> r_pt_list;
	vector<int4> pt_list_info;
	vector<int4> queue;
	vector<int4> queue0;

	//inizializzo i vettori per memorizzare livello e ofs delle celle di bitmask 
	for(int i = 0; i < levels; i++) {
		int sizex = bsx / (1 << i);
		int sizey = bsy / (1 << i);
		for (int y = 0; y < sizey; y++)
			for (int x = 0; x < sizex; x ++)
		    	if (bitmaskC[i][sizex * y + x] >= 0){
		      		/// memorizza livello
		      		map.host_grid_level_multi[bitmaskC[i][sizex * y + x]] = i;
		      		//carico coordinate
		      		map.host_ofs_blocks[bitmaskC[i][sizex * y + x]].x = x * (1 << i);
		      		map.host_ofs_blocks[bitmaskC[i][sizex * y + x]].y = y * (1 << i);
		      	}
	}	

	//iniziallizzo host_info 
	for (int i = 0; i < map.ncols; i++)
    	for (int j = 0; j < map.nrows; j++) {
      		map.host_info[i + map.ncols * j].x = 0;
    	}

    //inizializzo a false un array di dimensione tot_blocks per controllare di caricare una sola volta 
    //tutti i blocchi in pt_list
    for(int i = 0; i < tot_blocks; i++) 
    	assegnato[i] = 0;

    // esploro vicini
    // formato
    // codice spostamento per trovare l'amico: 0 (-1,-1), 1 (0 -1), 2 (1 -1), 3(-1 0), 4 (1 0), 5 (-1 1), 6 (0 1), 7 (1 1)
    // livello di multires: -1, 0, 1 (delta per raggiungere il livello amico)
    // codice vicino n1 (-1 muro)
    // info aggiuntiva n2: se livello  1 -->  quadrante del quadrato amico (codifica x+2*y, 0 (0 0) 1(1 0) 2 (0 1) 3 (1 1))
    //                     se livello -1 e amico non e' corner -->  secondo quadrato adiacente

    for(int i = 0; i < levels; i++) {
    	//int dbg = 0;
		int sizex = bsx / (1 << i);
		int sizey = bsy / (1 << i);
		for(int y = 0; y < sizey; y++)
	  		for(int x = 0; x < sizex; x++)
	    		if(bitmask[i][sizex * y + x] == i + 1){
	      			int side = 0;
	      			for(int dy = -1; dy <= 1; dy++)
	      				for(int dx = -1; dx <= 1; dx++)
							if(dx != 0 || dy != 0){
								int nx = x + dx;
							  	int ny = y + dy;
							  	int idx = 8 * bitmaskC[i][sizex * y + x] + side;
							  	if (nx < 0 || nx >= sizex ||
							      	ny < 0 || ny >= sizey ||
							      	bitmask[i][sizex * ny + nx] == -1) { // quadrato fuori
							    	map.neigh[idx].lev = 0;
							    	map.neigh[idx].n1 = -1;
							    	map.neigh[idx].n2 = -1;
							    	if(dbg) 
							      		printf("%d %d %d (%d %d), %d: %d %d [%d, %d]\n",i,x,y,dx,dy,bitmaskC[i][sizex*y+x],side,0,-1,-1); // muro
							  	}
		  						else {
								    if(bitmask[i][sizex * ny + nx] == i + 2){ //multires dopo, uso secondo posto per indicare se 0 (prima meta') o 1 (seconda meta' di adiacenza)
									    int bit = -1;
									    bit = (nx % 2) + 2 * (ny % 2);
									    map.neigh[idx].lev = 1;
									    map.neigh[idx].n1 = bitmaskC[i + 1][sizex / 2 * (ny / 2) + (nx / 2)];
									    map.neigh[idx].n2 = bit;
									    if(dbg)
									    	printf("%d %d %d (%d %d), %d: %d %d [%d, %d]\n",i,x,y,dx,dy,bitmaskC[i][sizex*y+x],
										    	side,
										     	1,
										     	bitmaskC[i+1][sizex/2*(ny/2)+(nx/2)],
										     	bit);
									}
								    else {
								    	if(i == 0) { // primo livello non scendo, non ci sono ulteriori livelli
								      		map.neigh[idx].lev = 0;
								      		map.neigh[idx].n1 = bitmaskC[i][sizex * ny + nx];
								      		map.neigh[idx].n2 = -1; // non usata
											if(dbg)
			  									printf("%d %d %d (%d %d), %d: %d %d %d\n",i,x,y,dx,dy,bitmaskC[i][sizex*y+x],side,0,bitmaskC[i][sizex*(ny)+(nx)]);
		      							}
		      							else { // guardo il livello precedente
											int ex = x * 2;
											int ey = y * 2;
											if(dx < 0) ex--;// adiacente sx, scavalco quadrato
											if(dy < 0) ey--;
											if(dx > 0) ex += 2; // adiacente destro, +1 per vicino nello stesso quad (a ris i) +1 scavalco quadrato
											if(dy > 0) ey += 2; // adiacente destro, +1 per vicino nello stesso quad (a ris i) +1 scavalco quadrato
											int ex1,ey1; // secondario (usato se si va a livello -1
											ex1 = ex;
											ey1 = ey;
											if(dx == 0) ex1++; // sopra o sotto, uso anche cella adiacente a destra
											if(dy == 0) ey1++; // sx o dx, uso anche cella adiacente +1			  
											if(bitmask[i - 1][sizex * 2 * ey + ex] == i + 1) { // mio livello (garantisco che era tutto il quadrato dello stesso livello)
										  		map.neigh[idx].lev = 0;
											  	map.neigh[idx].n1 = bitmaskC[i][sizex * ny + nx];
											  	map.neigh[idx].n2 = -1; // non usata			  
											  	if(dbg)
											    	printf("%d %d %d (%d %d), %d: %d %d [%d, %d]\n",i,x,y,dx,dy,bitmaskC[i][sizex*y+x],side,0,bitmaskC[i][sizex*(ny)+(nx)],-1);
											}
											else {
												int bit = -1;
											  	if(dx * dy == 0) // quando non diagonale, metto il secondo quadrato
											    	bit = bitmaskC[i-1][sizex * 2 * ey1 + ex1];
											  		map.neigh[idx].lev = -1;
											  		map.neigh[idx].n1 = bitmaskC[i - 1][sizex * 2 * ey + ex];
											  		map.neigh[idx].n2 = bit; 
											  		if(dbg)
											    		printf("%d %d %d (%d %d), %d: %d %d [%d, %d]\n",i,x,y,dx,dy,bitmaskC[i][sizex*y+x],side,-1,
												   			bitmaskC[i-1][sizex*2*(ey)+(ex)],
												   			bit
												   		);
											}
		      							}
		      
		    						}
		  						}
		  						side++;
							}
	      
	   		 	}
    }

    //per ogni punto di ogni retta ricavo il blocco a cui appartiene e aggiungo tutto il blocco a pt_list 
    for(int i = 0; i < g.punti_m.size(); i++) {
    	//int i = 1;{
    	int ii = (int)g.punti_m[i].z - 1;
    	int sizex = bsx / (1 << ii);

    	//ricavo le coordinate di blocco
		int x = (int)((g.punti_m[i].x - map.min.x) / map.dx / BLOCKSIZE_X);
		int y = (int)((g.punti_m[i].y - map.min.y) / map.dy / BLOCKSIZE_Y);

		int codice = bitmaskC[ii][sizex * (y / (1 << ii)) + (x / (1 << ii))];

		//coordinate del blocco in multirisoluzione
		int x_multi = codice % x_blocks;
		int y_multi = codice / x_blocks;

		if(assegnato[codice] == 0) {
			for(int y1 = 0; y1 < BLOCKSIZE_Y; y1++) 
				for(int x1 = 0; x1 < BLOCKSIZE_X; x1++) {
					int idx = (y_multi*BLOCKSIZE_Y+y1) * BLOCKSIZE_X * x_blocks + (BLOCKSIZE_X*x_multi+x1);

					point real_coord;
					int4 pt_info;


					map.host_info[idx].x = 1;

					real_coord.x = map.min.x + (BLOCKSIZE_X * x * map.dx) + x1 * map.dx * (1<<ii);
					real_coord.y = map.min.y + (BLOCKSIZE_Y * y * map.dy) + y1 * map.dy * (1<<ii);
					pt_info.x = x1;
					pt_info.y = y1;
					pt_info.z = codice;
					//printf("%d %d %d %d %lf %lf %d %d\n",x,y,x_multi,y_multi,real_coord.x,real_coord.y,codice,idx);
					r_pt_list.push_back(real_coord);
					pt_list_info.push_back(pt_info);
				}
			assegnato[codice] = 1;
		}
	}

	for(int pt_it = 0; pt_it < r_pt_list.size(); pt_it++) {
		int numbc = 0;
		int i_est = 0;			
		int n = (int)pol.points.size();
		point foo = r_pt_list[pt_it];
		int4 foo_info = pt_list_info[pt_it];

		int x_multi = foo_info.z % x_blocks;
		int y_multi = foo_info.z / x_blocks;

		int idx = (y_multi*BLOCKSIZE_Y+foo_info.y) * BLOCKSIZE_X * x_blocks + (BLOCKSIZE_X*x_multi+foo_info.x);

		for(int k = 0; k < n; k++) {  
			if( (pol.points[k].x < foo.x && pol.points[(k+1)%n].x >= foo.x) ||
			    (pol.points[k].x >= foo.x && pol.points[(k+1)%n].x < foo.x) ) {
				numbc = 1;
			  	F polyXi = pol.points[k % n].x;
			  	F polyXj = pol.points[(k + 1) % n].x;
			  	F polyYi = pol.points[k % n].y;
			  	F polyYj = pol.points[(k + 1) % n].y;
			  	if(polyYi + (foo.x - polyXi) / (polyXj - polyXi) * (polyYj - polyYi) < foo.y)
			    	i_est = 1 - i_est;
			  	//printf("%f %lf --> %d\n",polyYi+(foo.x-polyXi)/(polyXj-polyXi)*(polyYj-polyYi),foo.y,i_est);
			}
		}

		if (i_est == 0 || numbc == 0) { //!condizione per celle con xp < min(xpoint) o xp > max(xpoint)
			//map.host_info[idx].x = BIT_EXTERN;
			queue.push_back(foo_info);
		}
		else {
			//map.host_info[idx].x = 0; 
			queue0.push_back(foo_info);
		}
    }

    while(queue.size() > 0) {
	    int4 foo = queue[queue.size() - 1];
	    queue.erase(queue.end() - 1);

	    int x_multi = foo.z % x_blocks;
		int y_multi = foo.z / x_blocks;
		int idx = (y_multi*BLOCKSIZE_Y+foo.y) * BLOCKSIZE_X * x_blocks + (BLOCKSIZE_X*x_multi+foo.x);
		//printf("%d %d %d %d\n", foo.z, foo.x,foo.y,map.host_info[idx].x);

		if(map.host_info[idx].x != BIT_EXTERN) {
			//printf("aaa");
		    map.host_info[idx].x = BIT_EXTERN;

		    //aggiungo vicini alla coda

		    int dx = 0;
		  	int dy = 0;
		  	int sh = 0;																																																																																																																																																																																														
		  	int ox = 0;
		  	int oy = 0;

		  	// quale rappresenta la posizione del vicino: 0 -> N, 1 -> S, 2 -> W, 3 -> E
		  	//if(foo.z == 209)
		  	for(int quale = 0; quale < 4; quale++) {
			  	switch(quale){
			  		case 0: {dy = 1; sh = 6; ox = foo.x; oy = 0;} break;
			  		case 1: {dy = -1; sh = 1; ox = foo.x; oy = BLOCKSIZE_Y - 1;} break;
			  		case 2: {dx = -1; sh = 3; ox = BLOCKSIZE_X - 1; oy = foo.y;} break;
			  		case 3: {dx = 1; sh = 4; ox = 0; oy = foo.y;} break;
			  	}

			  	if ((quale == 0 && foo.y < BLOCKSIZE_Y - 1) ||
		      		(quale == 1 && foo.y > 0) ||
		      		(quale == 2 && foo.x > 0) ||
		      		(quale == 3 && foo.x < BLOCKSIZE_X - 1)) { // sono dentro quadrato corrente -> prelevo valore
			  			int4 pt;
			  			pt.x = ox;
			  			pt.y = oy;
			  			pt.z = foo.z;
			  			queue.push_back(pt);
			  			//printf("AAA %d\n",pt.z);	    	
		  		}else {
				    neigh_t neigh = map.neigh[8 * foo.z + sh];
				    char lev = neigh.lev;
				    //printf("%d\n",neigh.n1);
				    if(neigh.n1 != -1)
					    switch(lev) {
					    	int4 pt1;
					    	int base;

					    	case 0: //blocco vicino è alla stessa risoluzione del blocco corrente      
					      		pt1.x = ox;
					      		pt1.y = oy;
					      		pt1.z = neigh.n1;
					      		queue.push_back(pt1);
					      		//printf("0 BBB %d %d %d\n",pt1.z,pt1.x,pt1.y);
					      		break;

					      	case -1: //il blocco vicino è a risoluzione piu' alta del blocco corrente
				          			 //A.F.: prendo il vicino n1 o n2 in base alla posizione della cella nel blocco corrente
					      		int4 pt2;

					      		pt1.x = 0;
					      		pt1.y = 0;
					      		pt1.z = neigh.n1;
					      		pt2.x = 0;
					      		pt2.y = 0;
					      		pt2.z = neigh.n2;

					      		//queue.push_back(pt1);
					      		//queue.push_back(pt2);
					      		//printf("-1 CCC %d %d\n",pt1.z,pt2.z);

					      		break;

					      	case 1:
					      		pt1.x = ox;
					      		pt1.y = oy;
					      		pt1.z = neigh.n1;
					      		queue.push_back(pt1);
					      		//printf("1 DDD %d %d %d\n",pt1.z,pt1.x,pt1.y);

					      		break;
					      		/*if ((quale <= 1 && foo.x < BLOCKSIZE_X / 2) || // pesco su dim corretta!
						  			(quale > 1 && foo.y < BLOCKSIZE_Y / 2)) 
									pt1.z = neigh.n1;
					      		else
									pt1.z = neigh.n2;
						      	//per costruzione idx >=0 (cioe' esiste, altrimenti trovavo muro!)
					    	  	//posizione*2, calcolo direttamente indice del quadrato corretto
							  	//A.F.: base è l'indice della prima cella (delle due visto che ho ris maggiore) che devo prendere sul blocco vicino (*2 perchè ho ris maggiore con il doppio di celle)
						      	if (quale <= 1)
						      		base = (foo.x * 2) % BLOCKSIZE_X;
						      	else
							      	base = (foo.y * 2) % BLOCKSIZE_Y;

					       		// unico caso in cui ho 2 valori in output, carico other e other_
								switch(quale){
									case 0: {ox1 = base; oy1 = 0;} break;
									case 1: {ox1 = base; oy1 = BLOCKSIZE_Y - 1;} break;
									case 2: {ox1 = BLOCKSIZE_X - 1; oy1 = base;} break;
									case 3: {ox1 = 0; oy1 = base;} break;
								}*/
					      		
					    }
				}
			}
		}
  	}

  /*while (queue0.size()>0){
    int2 foo=queue0[queue0.size()-1];
    //    printf("-%4d--0> %d %d\n",queue0.size(),foo.x,foo.y);
    queue0.erase(queue0.end()-1);
      maps.host_info[foo.y * ncol + foo.x].x=0;
    //aggiungi neigh
    int2 n;
    n.x=foo.x+1;n.y=foo.y+0;if (n.x<ncol && maps.host_info[n.y * ncol + n.x].x==1){
      maps.host_info[n.y * ncol + n.x].x=0;
      queue0.push_back(n);
    }
    n.x=foo.x-1;n.y=foo.y+0;if (n.x>=0 && maps.host_info[n.y * ncol + n.x].x==1){
      maps.host_info[n.y * ncol + n.x].x=0;
      queue0.push_back(n);
    }
    n.x=foo.x;n.y=foo.y+1;if (n.y<nrow && maps.host_info[n.y * ncol + n.x].x==1){
      maps.host_info[n.y * ncol + n.x].x=0;
      queue0.push_back(n);
    }
    n.x=foo.x;n.y=foo.y-1;if (n.y>=0 && maps.host_info[n.y * ncol + n.x].x==1){
      maps.host_info[n.y * ncol + n.x].x=0;
      queue0.push_back(n);
    }
  }*/

	printf("Multiresolution matrix size: %d\n",x_blocks * y_blocks * BLOCKSIZE_X * BLOCKSIZE_Y);

	/*int border_top = map.last_slab % map.slabs_nrows;

	for(int i = map.first_slab; i <= map.last_slab; ++i) {  
		if(i <= (i - i%map.slabs_nrows + border_top) && (i % map.slabs_nrows) != 0) {
			string file = path;
			file = file + to_string(i) + ".grd";

			ifstream f(file);

			char type[10];
			int dx, dy, max, min;
			int2 m_point, M_point;
			int row = 0;
			int col = 0;

			if(!f.fail()) {
				f >> type;
				f >> dx >> dy;
				f >> m_point.x >> M_point.x;
				f >> m_point.y >> M_point.y;
				double m1,m2;
				f >> m1 >> m2;

				while(!f.eof()) {
					double h_val;
					int2 slab;

					f >> h_val;
	
					if(col == map.dxs) {
						++row;
						col = 0;
					}
	 
					if(h_val != 1.70141e38) {
						//coordinate reali del punto che viene letto
						slab.x = m_point.x + col;
						slab.y = m_point.y + row;
						
						//coordinate di blocco
						int x = (int)(((slab.x - map.min.x) / map.dx) / BLOCKSIZE_X);
						int y = (int)(((slab.y - map.min.y) / map.dy) / BLOCKSIZE_Y);

						int found = -1;
						int codice = -1;

						for(int ii = 0; ii < levels; ii++){
						    int sizex = bsx / (1<<ii);
						    int sizey = bsy / (1<<ii);
					        if (bitmask[ii][sizex * (y / (1 << ii)) + (x / (1 << ii))] == ii + 1) {
					          	found = ii;
					          	codice = bitmaskC[ii][sizex * (y / (1 << ii)) + (x / (1 << ii))];

					          	//offset nel blocco
					          	int x1 = ((int)((slab.x - map.min.x) / map.dx) % (BLOCKSIZE_X * (1 << found))) / (1 << found);
					          	int y1 = ((int)((slab.y - map.min.y) / map.dy) % (BLOCKSIZE_Y * (1 << found))) / (1 << found);
					          	//printf("%d %d\n",x1,y1);

								//offset in host_grid_multi
								int x_multi = codice % x_blocks;
				        		int y_multi = codice / x_blocks;
				        		//printf("%d %d %d\n",x_multi, y_multi,codice);

				        		//ricavo ora l'indice corretto in host_grid_multi
				        		int idx = (y_multi*BLOCKSIZE_Y+y1) * BLOCKSIZE_X * x_blocks + (BLOCKSIZE_X*x_multi+x1);
				        		//if(codice == 0)
				        		//printf("%d %d %d %d %d %d %d %d %d\n",slab.x,slab.y,x,y, x1,y1, x_multi,y_multi,idx);

				        		map.host_grid_multi[idx].w += h_val;
				        		++counter[idx];

					        }  	
						}
					}
					++col;
				}
				f.close();
				if(dbg) printf("%s loaded.\n",file.c_str());
			} else {
				printf("Unable to open file. %s not found.\n", file.c_str());
				exit(1);
			}
		}
	}

	for(int j = 0; j < x_blocks*y_blocks*BLOCKSIZE_X*BLOCKSIZE_Y; ++j)
		if(counter[j] != 0)
			map.host_grid_multi[j].w /= counter[j];*/

}


int main() {
	
	string map_file = "map_info.txt";
	string bln_file = "polygons/bln2";
	string pts_file = "polygons/bln_raster.PTS";
	string path = "slabs/";

	read_map(map_file);
	read_pts(pts_file);
	readBLN(bln_file);
	multires(path);

	return 0;

}