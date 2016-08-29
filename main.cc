#include <iostream>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>

//#include <cuda_profiler_api.h> 

#include "stdlib.h"
#include "types.h"

#if (defined(_WIN32) || defined(__WIN32__) || defined(WIN32))
#define WIN 1
#endif

#if !defined(NOOPENGL) && !defined(WIN)
#include "3d_window.h"
#endif
#include "math.h"
#include "wavesimulator.h"
#include "timing.h"
#include "dsaa_reader.h"
#include "bln_reader.h"
#include "sez_reader.h"

#include "svnversion.h"

#include <cuda.h>
#include <cuda_runtime.h>

#ifndef max
#define max(a,b) (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

extern int BLOCKSIZE_X;
extern int BLOCKSIZE_Y;

char base_path1[256];

struct global g;
struct maps_ maps;

vector<cons> conss;  // elenco vincoli
vector<polygon> polygons;

long kernelcalls=0;

int simulate;
int stepperframe;
int kernelflops;


int argc;
char ** argv;


int hi=0; // tutto multi ad alta risoluzione



void help(){
  printf("Shallow water simulator with GPU\nUsage: %s parameters_file options\n\n",argv[0]);
  printf("Options:\n\t -debug=n\tn=0,1 flag for debugging information\n");
  printf("\t -multi=n\tn=0,1,hi flag for enabling multiresolution (hi=max multires everywhere)\n");
  printf("\t -bcc-level=n\t level at which the boundary conditions are rastered (range 1...4)\n");
  printf("\t -order=n\tn=1,2 equation order\n");
  printf("\t -opengl=n\tn=0,..n 0=step by step >0 how many iteration per shown frame\n");
  printf("\t -gpu=n\tn=0,..n selected device\n");
  printf("\t -resume=s\ts=filename suffix for resuming computation\n");
  printf("\t -por=n\tn=0,1 flag for enabling porosity matrix\n");
  printf("\nDecode utility\nUsage: %s -decode input_matrix output_matrix [-frames x y] [-all] [-res=n] [-range minx maxx miny maxy] [-point x y] [-ndigits n] [-sez x1 y1 x2 y2]\n",argv[0]);
  printf("\t -all enables the conversion of the whole domain (overrides -range)\n");
  printf("\t -res=n\t n=1..4 selects the multiresolution level to sample (1 most detailes, 4 most coarse)\n");
  printf("\t -frames x y\t selects from the filename input_matrix (must contain the -output-XXXX.YYY string) the files that match the numbers x .. y. Output maintains the frame numbering\n");
  printf("\t -ndigits n \t n specifies the number of digits of numbers in the output file\n");
  printf("\t -range minx maxx miny maxy \t specifies the rectangle rendered to file\n");
  printf("\t -point x y \t specifies the sampled point (in association with frames, output in the same file). Equivalent to -range x x y y\n");
  exit(1);
}

/*

void dbg_block(int i, int bound_blocks,int tot_blocks){
  if (i>=tot_blocks) return;
  int x_blocks=1<<(int)(ceil(log((float)maps.tot_blocks)/log((float)2)/(float)2)); 
  int y_blocks=(maps.tot_blocks-1)/x_blocks+1;
  int x=i%x_blocks;
  int y=i/x_blocks;
  printf("Block %d, xb%d %d, %d %d, ofs %d %d, lev %d\n",i,x_blocks,y_blocks,x,y,maps.host_ofs_blocks[i].x,maps.host_ofs_blocks[i].y,maps.host_grid_level_multi[i]);
  int size=BLOCKSIZE_X;
  for (int y1=size-1;y1>=0;y1--){
    for (int x1=0;x1<size;x1++) {

      int idx=(y*BLOCKSIZE_Y+y1)*x_blocks*BLOCKSIZE_X+BLOCKSIZE_X*x+x1;

      unsigned char val=maps.host_info_multi[idx].x;
      if (val==0)
	printf("  .");
      else{
	if (val==BIT_EXTERN)
	  printf("  X");
	else
	  printf("%3d",val);			
      }
    }
    // valori
    printf("\t");
    for (int x1=0;x1<size;x1++) {
      int x=i%x_blocks;
      int y=i/x_blocks;
      int idx=(y*BLOCKSIZE_Y+y1)*x_blocks*BLOCKSIZE_X+BLOCKSIZE_X*x+x1;
      F val=maps.host_grid_multi[idx].w;
      printf("%.2e ",val);			
    }
    printf("\n");
  }
  printf("\n");
}

*/
void multires(char* path) {
	printf("\nMultiresolution processing ------------\n");
    int ncols = maps.ncols;
    int nrows = maps.nrows;
    int levels = 4; //livelli di multirisoluzione
    int lblock = (1<<(levels - 1));
    // la matrice e' allargata +2
    int sx = maps.ncols - 2;
    int sy = maps.nrows - 2;
      
    // la mappa multires e' aumentata a multipli della risoluzione minima, 
    //in modo da avere quadrati divisi esattamente 
    int esx = ((((sx - 1) / BLOCKSIZE_X + 1) - 1) / lblock + 1) * lblock * BLOCKSIZE_X;
    int esy=((((sy - 1) / BLOCKSIZE_Y + 1) - 1)/ lblock + 1) * lblock * BLOCKSIZE_Y;
    maps.esx = esx;
    maps.esy = esy;

    int bsx = esx / BLOCKSIZE_X; // bitmap size (ridotta di blocksize)
    int bsy = esy / BLOCKSIZE_Y;


    //aggiorno range
    // ho tolto il bordo +2
    // ma allargato anche eventualmente la parte alta
	maps.minx_map += maps.dx;
	maps.miny_map += maps.dy;
	maps.maxx_map = maps.minx_map + maps.dx * (esx - 1);
	maps.maxy_map = maps.miny_map + maps.dy * (esy - 1);	


    printf("Orig. size %dx%d, largest block %d, new size %dx%d, blocks %dx%d, bound map %dx%d\n",sx,sy,lblock,esx,esy,bsx,bsy,
		maps.host_info_x,maps.host_info_y);

    // bitmask: valori -1..levels-1, 
    //    -1 = quadrato copre completamente zona esterna
    // 

    char** bitmask = (char **) malloc(levels * sizeof(char*));
    int** bitmaskC = (int **) malloc(levels * sizeof(int*));
    uchar4** bitmaskS = (uchar4 **) malloc(levels * sizeof(uchar4*));//sides

    F* Z; // supporto per downsampling

    //host_info_m
    maps.host_info_m = (uchar4 **) malloc(levels * sizeof(uchar4*));//sides
    maps.host_info_x_m = (int*) malloc(levels * sizeof(int));//sides
    maps.host_info_y_m = (int*) malloc(levels * sizeof(int));//sides

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

		int nx = maps.ncols;
		int ny = maps.nrows;
		int scale = (1<<i);
		nx = (nx - 3) / scale + 3;
		ny = (ny - 3) / scale + 3;

		if (i == 0) // caso livello 0, metto direttamente dentro
	  		Z = (F*) malloc(nx * ny * sizeof(F));
    }

    if(hi) // speciale: se tutto ad alta risoluzione
    for(int y = 1; y < nrows-1; y++){  
		for(int x = 1; x < ncols-1; x ++)
			if(maps.btm_map[y * ncols + x] < 1e10) {
				int bx = (x - 1) / BLOCKSIZE_X; // shifto uno indietro (c'e' bordo!)
			    int by = (y - 1) / BLOCKSIZE_Y;
			    bitmask[0][bsx*by+bx] = 1;
			}
    }

	// lancio rasterizzazione bln alla risoluzione richiesta
	// ricavo richieste blocchi bln, unisco a richieste multires gia' presenti
	// lancio suddivisione blocchi

    // per ogni livello
    //   raster bln U tappi U correzioni Z (da livello>=1 guardando livello prec)
    //   codici livello di patch non vuote
    //   copio
      
    { // lavoro sul livello di rendering bcc 
		int i = g.global_level_bcc - 1;
		int nx = maps.ncols;
		int ny = maps.nrows;
		int scale = (1<<i);
		nx = (nx - 3) / scale + 3;
		ny = (ny - 3) / scale + 3;	
		for(int x = 0; x < nx; x++)
			for(int y = 0; y < ny; y++)
		    	if(x == 0 || x == nx - 1 || y == 0 || y == ny - 1)
		    		Z[y * nx + x] = 10e10;
			else
		    	Z[y * nx + x] = 0;


		// aggiunto riferimento a bitmask con livelli di multires (qui non vogliamo correggere i tappi!)
		readBLNBCC(path, i, Z, i+1 == g.global_level_bcc, NULL, 0);

		if(0) {
			for(int j = nrows - 1; j >= 0; j--) {
				printf("%d--: ", j);
		  		for(int i = 0; i < ncols; i++)
		    		if(maps.btm_map[j * ncols + i] >= 10e10) 
		    			printf("X");
		    		else
						printf(" ");
		  		printf("\n");
			}

			for(int j = nrows - 1; j >= 0; j--){
				printf("%d--: ",j);
			  	for(int i = 0; i < ncols; i++)
			    	if(Z[j * ncols + i] >= 10e10) 
			    		printf("X");
			    	else
			      		printf(" ");
			  	printf("\n");
			}
		}

		maps.host_info_x_m[i] = maps.host_info_x;
		maps.host_info_y_m[i] = maps.host_info_y;
		maps.host_info_m[i] = (uchar4 *) malloc( maps.host_info_x * maps.host_info_y * sizeof(uchar4));//sides
		//copia valori
		//printf("copy %d x %d\n",maps.host_info_x,maps.host_info_y);
		for(int x = 0; x < maps.host_info_x; x++)
	  		for (int y = 0; y < maps.host_info_y; y++)
	    		maps.host_info_m[i][y * maps.host_info_x + x]=maps.host_info[y * maps.host_info_x + x];
    }

	/*
	//downsample: 
	// i tappi da btm sono fatti sparire e aggiunti alle richieste in g.punti_m 
	// 
	// quindi cancello tutti i tappi livelli >0 (reintrodotti tramite elenco subito dopo!)
	for (int x=1;x<nx-1;x++)
	  for (int y=1;y<ny-1;y++){
	    F temp=0;
	    int ext=0;
	    for (int dx=0;dx<scale;dx++)
	      for (int dy=0;dy<scale;dy++){
		int ix=(x-1)*scale+dx+1;
		int iy=(y-1)*scale+dy+1;
		//printf("%d %d\n",iy,ix);
		if (ix>=ncols || iy>=nrows || maps.btm_map[iy * ncols + ix]>10e10)
		  ext=1;
		else
		  temp+=maps.btm_map[iy * ncols + ix];
	      }
	    temp/=scale*scale;
	    if (ext){
	      if (i==(g.global_level_bcc-1)){ // sulla ris piu' bassa prelevo tappi e impongo richieste di risoluzione
		Z[y*nx+x]=1.70141e+038;
	      }
	      else
		Z[y*nx+x]=0; // no tappo!
	    }
	    else		
	      Z[y*nx+x]=temp;		  
	  }

	for (int i1=0;i1<g.punti_m.size();i1++)
	  if (g.punti_m[i1].w==1 && 
	      g.punti_m[i1].z==i+1
	      ) // e' tappo --> scrivo direttamente dentro Z
	    {
	      int x=((g.punti_m[i1].x-maps.minx_map)/maps.dx-1); // scalo bordo maps
	      int y=((g.punti_m[i1].y-maps.miny_map)/maps.dy-1);
	      // scalo rispetto a multirisoluzione corretta
	      x=(x)/(1<<((int)g.punti_m[i1].z-1))+1;
	      y=(y)/(1<<((int)g.punti_m[i1].z-1))+1;
	      // scrivo punto
	      printf("TAPPO manuale: %d %d livello %d\n",x,y,(int)g.punti_m[i1].z);
	      Z[y*nx+x]=1.70141e+038;
	    }


	readBLNBCC(path,i,Z,i+1==g.global_level_bcc);

	maps.host_info_x_m[i]=maps.host_info_x;
	maps.host_info_y_m[i]=maps.host_info_y;
	maps.host_info_m[i] = (uchar4 *) malloc( maps.host_info_x *maps.host_info_y * sizeof(uchar4));//sides
	//copia valori
	//printf("copy %d x %d\n",maps.host_info_x,maps.host_info_y);
	for (int x=0;x<maps.host_info_x;x++)
	  for (int y=0;y<maps.host_info_y;y++){
	    maps.host_info_m[i][y*maps.host_info_x+x]=maps.host_info[y*maps.host_info_x+x];
	  }

	if (0)
	  for(int j = 0; j < maps.host_info_y; j ++){
	  for(int ii = 0; ii < maps.host_info_x; ii ++){
	    int v=maps.host_info_m[i][ii + maps.host_info_x * j].x;
	    if (v!=0 && v!=BIT_EXTERN)
	      printf("%d",maps.host_info_m[i][ii + maps.host_info_x * j].y);
	    else{
	      if (v==0)
		printf(" ");
	      else
		printf("x");
	    }
	  }
	  printf("\n");
	}

*/

    float thr=10; // gradiente 
    for(int y = 1; y < nrows-1; y++)  
		for(int x = 1; x < ncols - 1; x++) {
	  		int high = 0;
	  		float delta;
			delta = 0;
			if (x + 1 < ncols - 1) // salto bordo esterno
				delta = fabs(maps.btm_map[y * ncols + x] - maps.btm_map[y * ncols + x + 1]);
			if (delta > thr && delta < 10e3) 
				high = 1;
			if (y + 1 < nrows - 1) // salto bordo esterno
				delta = fabs(maps.btm_map[y * ncols + x] - maps.btm_map[(y + 1) * ncols + x]);
			if (delta > thr && delta < 10e3) 
				high = 1;
		    // RV: added here to get uniform resoution in multigrid code.
		    //high=1;
			int bx = (x - 1) / BLOCKSIZE_X; // shifto uno indietro (c'e' bordo!)
			int by = (y - 1) / BLOCKSIZE_Y;
			if(high)
				bitmask[0][bsx * by + bx] = 1;
		}

    printf("Processing PTS %d points\n",g.punti_m.size());
    for(int i = 0; i < levels; i++) {
		int sizex = bsx / (1<<i);
		int sizey = bsy / (1<<i);

		// impone richieste di risoluzione
		if(hi == 0) // se non tutto ad alta
		for(int i1 = 0;i1 < g.punti_m.size(); i1++){
			//printf("%d: %f %f %d %d\n",i1,g.punti_m[i1].x,g.punti_m[i1].y,(int)g.punti_m[i1].z,(int)g.punti_m[i1].w);
		    if(//anche se tappi  
			//g.punti_m[i].w==0 &&// non e' tappo, ma solo richiesta di risoluzione
			(int) g.punti_m[i1].z - 1 == i) { // lavora al livello corrente
				int x = ((g.punti_m[i1].x - maps.minx_map) / maps.dx - 1) / BLOCKSIZE_X; // scalo bordo maps
				int y = ((g.punti_m[i1].y - maps.miny_map) / maps.dy - 1) / BLOCKSIZE_Y;
				// ora allineo rispetto patch rispetto al livello di multirisoluzione richiesto
				x = x / (1<<i) * (1<<i);
				y = y / (1<<i) * (1<<i);

				printf("Force level %d point %f %f --> %d %d (%d %d\n",(int)g.punti_m[i1].z,g.punti_m[i1].x,g.punti_m[i1].y,x,y,bsx,bsy);
				for(int ii = 0; ii < (1<<i); ii++)
			  		for(int jj = 0; jj < (1<<i); jj++)
			    		if(x + ii >= 0 && x + ii < bsx && y + jj >= 0 && y + jj < bsy)
			    			if(bitmask[0][bsx*(y+jj)+x+ii]==0 || bitmask[0][bsx*(y+jj)+x+ii]>g.punti_m[i1].z) // conquista la cella solo se non c'e' una risoluzione piu' alta (in caso la lascia stare!)
			      				bitmask[0][bsx * (y + jj) + x + ii] = g.punti_m[i1].z;
		    }
		}
    }

    for(int i = 0; i < levels; i++) {
		int sizex = bsx / (1<<i);
		int sizey = bsy / (1<<i);

		if(0)
      		for(int i1 = 0;i1 < levels; i1++) {
				int dbg = 1;
				int sizex = bsx / (1<<i1);
				int sizey = bsy / (1<<i1);
				for(int y = sizey-1; y >=0; y--) {  
	  				for(int x = 0; x < sizex; x ++) {
	    				if(bitmask[i1][sizex * y + x] == -1)
	      					if(dbg) 
	      						printf("! ");
	    				else {
	      					int ok = 1;
	      					if(bitmask[i1][sizex * y + x] == i1 + 1) {
								//verifico se c'e' almeno un valore != 0
								int size = BLOCKSIZE_X;
								for(int x1 = x * size; x1 < (x + 1) * size; x1++)
		  							for(int y1 = y * size; y1 < (y + 1) * size; y1++)
		    							if(x1+1<maps.host_info_x_m[i1]-1 && y1+1<maps.host_info_y_m[i1]-1) { // salto bordo esterno
		      								if(maps.host_info_m[i1][(y1+1) * maps.host_info_x_m[i1] + (x1+1)].x!=0 &&
			  									maps.host_info_m[i1][(y1+1) * maps.host_info_x_m[i1] + (x1+1)].x!=BIT_EXTERN) // ricordati lo shift!!
												ok=0;
		    								}
	      					}	      
	      					if(ok || bitmask[i1][sizex * y + x] != i1 + 1)
								if(dbg) 
									printf("%d ",bitmask[i1][sizex*y+x]);//<=i+1?bitmask[i][sizex*y+x]:0));
	      					else
								if(dbg) 
									printf("%dB",bitmask[i1][sizex*y+x]);			
	    				}
	  				}
	  				if(dbg) 
	  					printf("\n");
				}
      		}

		printf("lev -> %d level bcc %d\n",i,g.global_level_bcc-1);
		/// gestione bcc
		if(i == g.global_level_bcc - 1)  //importo le celle	  
	  		for(int y = 0; y < sizey; y++) {
	    		int dbg = -1;//sizey+1;
	    		for(int x = 0; x < sizex; x++) { // cella
					//cerco se c'e' una condizione nella cella
					int cond = 0;
					int ext = 1;

					//garantire interpolazione blocco vicino a risoluzione piu' bassa
					int bw = 0; //condizioni sul bordo di spessore 2 celle
					int bn = 0;
					int bs = 0;
					int be = 0;

					for(int dx = 0; dx < BLOCKSIZE_X; dx++)
		  				for(int dy = 0;dy < BLOCKSIZE_Y; dy++) {
		    				int px = x * BLOCKSIZE_X + dx + 1; // salto bordo
						    int py = y * BLOCKSIZE_Y + dy + 1;
						    if(px < maps.host_info_x_m[i] - 1 && py < maps.host_info_y_m[i] - 1) { // controllo di non uscire dalla matrice originale ridimensionata (qui lavoro con quadrati che estendono la matrice!)
						      	// host info sempre con shift di 1 
							    if(maps.host_info_m[i][py * maps.host_info_x_m[i] + px].x != BIT_EXTERN && 
								  maps.host_info_m[i][py * maps.host_info_x_m[i] + px].x != 0) {
									cond = 1;
									if (dx < 2) bw = 1;
									if (dx >= BLOCKSIZE_X - 2) be = 1;
									if (dy < 2) bs = 1;
									if (dy >= BLOCKSIZE_Y - 2) bn = 1;
			      				}
			      				if(maps.host_info_m[i][py * maps.host_info_x_m[i] + px].x != BIT_EXTERN)
									ext = 0;
			      				if(1 == 0 && y == dbg)
									printf("\n%d %d (%d %d): %d, ext %d cond %d\n",px,py,
							    		maps.host_info_x_m[i],
							    		maps.host_info_y_m[i],maps.host_info_m[i][py * maps.host_info_x_m[i] + px].x,ext,cond);
					   		}
		  				}
					if(0)
						if(dbg < 0 || dbg == y) 
	  						if(ext) 
	  							printf("x ");
	  						else
	    						printf("%d ",cond);

					int dbgc = 0;
					if(be || bw || bs || bn) {
	  					if(y > 0 && bs){
	    					if(bitmask[i][sizex * (y-1) + x] == 0 || bitmask[i][sizex * (y-1) + x] > i+1) { // forzo a stessa risoluzione (potrebbe non essere ancora stato assegnato
	      						bitmask[i][sizex * (y - 1) + x] = i + 1;
	      						if(dbgc) 
	      							printf("forzo blocco lev %d, %d %d --> %d %d\n",i,x,y,x,y-1);
	      						if(bitmask[i][sizex * (y - 1) + x] > i + 1)
									if(dbgc) 
										printf("warning: richiesta era per una risoluzione piu' bassa\n");
	    					}
		    				else 
		      					if(bitmask[i][sizex * (y-1) + x] > 0 && bitmask[i][sizex * (y-1) + x] <= i) { 
									printf("ERRORE: non posso forzare blocco (bcc vicine a bordo blocco) lev %d, %d %d <--> %d %d\n",i,x,y,x,y-1);		    
									exit(2);
		      					}
	  					}
	  					if(y < sizey - 1 && bn) {
	    					if(bitmask[i][sizex * (y+1) + x] == 0 || bitmask[i][sizex * (y+1) + x] > i+1) { // forzo a stessa risoluzione (potrebbe non essere ancora stato assegnato
	      						bitmask[i][sizex * (y + 1) + x] = i + 1;
	      						if(dbgc) 
	      							printf("forzo blocco lev %d, %d %d --> %d %d\n",i,x,y,x,y+1);
	      						if(bitmask[i][sizex * (y + 1) + x] > i + 1)
									if(dbgc)
										printf("warning: richieta era per una risoluzione piu' bassa\n");
	    					}
		    				else 
		      					if(bitmask[i][sizex * (y+1) + x] > 0 && bitmask[i][sizex * (y+1) + x] <= i) {
									printf("ERRORE: non posso forzare blocco (bcc vicine a bordo blocco) lev %d, %d %d <--> %d %d\n",i,x,y,x,y+1);		    
									exit(2);
		      					}
	  					}
	  					if(x < sizex - 1 && be) {
	    					if(bitmask[i][sizex * (y) + x+1] == 0 || bitmask[i][sizex * (y) + x+1] > i+1) { // forzo a stessa risoluzione (potrebbe non essere ancora stato assegnato
						    	bitmask[i][sizex * (y) + x + 1] = i + 1;
						      	if(dbgc) 
						      		printf("forzo blocco lev %d, %d %d --> %d %d\n",i,x,y,x+1,y);
						      	if(bitmask[i][sizex * (y) + x+1] > i + 1)
									if(dbgc) 
										printf("warning: richiesta era per una risoluzione piu' bassa\n");
	    					}
						    else
						    	if (bitmask[i][sizex * (y) + x+1] > 0 && bitmask[i][sizex * (y) + x+1] <= i) {
									printf("ERRORE: non posso forzare blocco (bcc vicine a bordo blocco) lev %d, %d %d <--> %d %d\n",i,x,y,x+1,y);		    
									exit(2);
						      	}
						}
	  					if(x > 0 && bw) {
	    					if(bitmask[i][sizex * (y) + x-1]==0 || bitmask[i][sizex* (y) + x-1] > i+1) { // forzo a stessa risoluzione (potrebbe non essere ancora stato assegnato
	    						bitmask[i][sizex* (y) + x-1] = i + 1;
	    						if(dbgc) 
	    							printf("forzo blocco lev %d, %d %d --> %d %d\n",i,x,y,x-1,y);
	    						if(bitmask[i][sizex * (y) + x-1] > i + 1)
	      							if(dbgc) 
	      								printf("warning: richiesta era per una risoluzione piu' bassa\n");
	    					}
	    					else 
	      						if (bitmask[i][sizex * (y) + x-1] > 0 && bitmask[i][sizex * (y) + x-1] <= i) {
									printf("ERRORE: non posso forzare blocco (bcc vicine a bordo blocco) lev %d, %d %d <--> %d %d\n",i,x,y,x+1,y);		    
									exit(2);			
								}
	  					}
					}

      				if(!ext && cond) { // devo forzare il livello corrente
						// se una cella contiene un valore di livello <=i e >0, errore (richiesta per bcc incompatibile)
						if(bitmask[i][sizex * y + x] > 0 && bitmask[i][sizex * y + x] <= i){
	  						printf("Error: boundary condition (level %d) over a region with level %d. Decrease to -bcc-level=%d.\n",
		 						i+1,bitmask[i][sizex*y+x],bitmask[i][sizex*y+x]);
	  						exit(2);
						}
						else
	  						bitmask[i][sizex * y + x] = i + 1;
					}
      				if(ext)
						bitmask[i][sizex * y + x] = -1;
      			}
    			if(0)
    				if(dbg < 0 || dbg == y) 
      					printf("\n");
    		}
	  

		//completa il quadrato 2x2 (espande eventuali celle vuote) per permettere il collasso al livello superiore 
		for(int x = 0; x < sizex; x++)
	  		for(int y = 0; y < sizey; y++)
	    		if(bitmask[i][sizex * y + x] > 0) { // se quadrato attivo
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
		      		if(bitmask[i][sizex * y + x] == i + 1) { // se quadrato attivo
			    		for(int dx = -1; dx <= 1; dx++)
			  				for(int dy = -1; dy <= 1; dy++)
			    				if(dx != 0 || dy != 0) {
						      		int nx = x + dx;
						      		int ny = y + dy;
						      		if(nx < 0 || nx >= sizex || ny < 0 || ny >= sizey) {
						      		}
			      					else {
										if(bitmask[i][sizex * ny + nx] == 0) // manca un vicino ---> deve esistere ad un livello di multires piu' alto
				  							bitmask[i + 1][(sizex / 2) * (ny / 2) + (nx / 2)] = i + 2;
										else {
				  							if(bitmask[i][sizex * ny + nx] != -1 &&
				      							abs(bitmask[i][sizex*ny+nx] - bitmask[i][sizex*y+x]) > 1) {
				    							printf("error: incompatible multiresolution requests level %d in %d %d (%f %f) (lev %d) - %d %d (lev %d)\n",
					   								i+1,x,y, maps.minx_map+maps.dx*(1<<i)*x , maps.miny_map+maps.dy*(1<<i)*y,
					   								bitmask[i][sizex*y+x],
					   								nx,ny,
					   								bitmask[i][sizex*ny+nx]);
				    							exit(2);
				  							}			    
										}
		      						}
		    					}
	      			}
		}
	}

    // riporta indietro i valori (per comodita' di visualizzazione) e 
    if(1)
    	for(int i = levels - 1; i >= 0; i--) {
			int sizex = bsx / (1<<i);
			int sizey = bsy / (1<<i);

			for(int x = 0; x < sizex; x++)
	  			for(int y = 0; y < sizey; y++) {
	    			// assicura ultimo livello tutto calcolato
	    			if(i == levels - 1 && bitmask[i][sizex * y + x] == 0)
	      				bitmask[i][sizex * y + x] = i + 1;

	    /// clean: se tutto il quadrato e' bit_extern --> tolgo quadrato (codice=-1)
	    // if (i+1>=g.global_level_bcc) // se sono sopra o uguale a risoluzione rasterizzazione

	    /*	    if (0)
	      { //se non cancellato
	      int dbg=0;//i==2;
	      int fuori=1;
	      // check fuori
	      //int size=(1<<i)*BLOCKSIZE_X/(1<<(g.global_level_bcc-1));
	      int size=BLOCKSIZE_X;
	      for (int x1=x*size;x1<(x+1)*size;x1++){
		for (int y1=y*size;y1<(y+1)*size;y1++)
		  if (x1+1<maps.host_info_x_m[i]-1 && y1+1<maps.host_info_y_m[i]-1){ // salto bordo esterno
		    if (dbg) printf("%d ",maps.host_info_m[i][(y1+1) * maps.host_info_x_m[i] + (x1+1)].x!=BIT_EXTERN);
		    if (maps.host_info_m[i][(y1+1) * maps.host_info_x_m[i] + (x1+1)].x!=BIT_EXTERN) // ricordati lo shift!!
		      fuori=0;
		  }
		if (dbg) printf("\n");
	      }
	      if (dbg) printf("fuori %d %d: %d\n",x,y,fuori);
	      if (0)
	      if (fuori==1)
		bitmask[i][sizex*y+x]=-1;
	      }
	    */

					if(i > 0 && (bitmask[i][sizex*y+x] >= i+1 || bitmask[i][sizex*y+x] == -1)) { // se quadrato attivo o fuori
				    	bitmask[i - 1][sizex * 2 * (y*2+0) + (x*2+0)] = bitmask[i][sizex * y + x];
				    	bitmask[i - 1][sizex * 2 * (y*2+0) + (x*2+1)] = bitmask[i][sizex * y + x];
				    	bitmask[i - 1][sizex * 2 * (y*2+1) + (x*2+0)] = bitmask[i][sizex * y + x];
				      	bitmask[i - 1][sizex * 2 * (y*2+1) + (x*2+1)] = bitmask[i][sizex * y + x];
				    }
	  			}
      	}

    for(int i1 = 0; i1 < levels; i1++) {
		int dbg = 0;
		int sizex = bsx / (1<<i1);
		int sizey = bsy / (1<<i1);
		for(int y = sizey-1; y >=0; y--) {  
			if(0)
	  			printf("%3d: ",y);
			for(int x = 0; x < sizex; x++){
	    		if(bitmask[i1][sizex * y + x] == -1)
	      			if(dbg) 
	      				printf("! ");
	    		else {
	      			int ok = 1;
			      	if(bitmask[i1][sizex * y + x] == i1 + 1 && i1 == g.global_level_bcc - 1) {
						//verifico se c'e' almeno un valore != 0
						int size = BLOCKSIZE_X;
						for(int x1 = x * size; x1 < (x + 1) * size; x1++)
				  			for(int y1 = y * size; y1 < (y + 1) * size; y1++)
				    			if(x1+1 < maps.host_info_x_m[i1]-1 && y1+1 < maps.host_info_y_m[i1]-1) // salto bordo esterno
				      				if(maps.host_info_m[i1][(y1+1) * maps.host_info_x_m[i1] + (x1+1)].x != 0 &&
					  				   maps.host_info_m[i1][(y1+1) * maps.host_info_x_m[i1] + (x1+1)].x != BIT_EXTERN) // ricordati lo shift!!
										ok = 0;
	      			}	      
	      			if(ok || bitmask[i1][sizex * y + x] != i1 + 1)
						if(dbg) 
							printf("%d ",bitmask[i1][sizex*y+x]);//<=i+1?bitmask[i][sizex*y+x]:0));
	      			else
						if(dbg) 
							printf("%dB",bitmask[i1][sizex*y+x]);
	    		}
	  		}
	  		if(dbg) 
	  			printf("\n");
		}
    }
//--------------------------------------------------------------------------------------------------------------------
      /// celle B non necessariamente usate, dato che sono solo la parte bln, 
      /// se tappi -> cambia rendering (da rifare) con nuove Z dei tappi

      // per ogni livello i=3..0
      //   raster bln U tappi U correzioni Z (da livello i<3 guardando livello i+1)

      // per ogni livello i=0..3
      //   codici livello di patch non vuote -- prima bound
      // per ogni livello i=0..3
      //   copio

      for (int i=levels-1;i>=0;i--){
	printf("compute levl %d\n",i);
	int nx=maps.ncols;
	int ny=maps.nrows;
	int scale=(1<<i);
	nx=(nx-3)/scale+3;
	ny=(ny-3)/scale+3;	

	for (int x=0;x<nx;x++){
	  Z[0*nx+x]=0;//1.70141e+038;
	  Z[(ny-1)*nx+x]=0;//1.70141e+038;
	}
	for (int y=0;y<ny;y++){
	  Z[y*nx+0]=0;//1.70141e+038;
	  Z[y*nx+nx-1]=0;//1.70141e+038;
	}

	for (int x=1;x<nx-1;x++)
	  for (int y=1;y<ny-1;y++){
	    F temp=0;
	    int ext=0;
	    for (int dx=0;dx<scale;dx++)
	      for (int dy=0;dy<scale;dy++){
		int ix=(x-1)*scale+dx+1;
		int iy=(y-1)*scale+dy+1;
		//printf("%d %d\n",iy,ix);

		if (ix>=ncols || iy>=nrows || maps.btm_map[iy * ncols + ix]>10e10)
		  ext=1;
		else
		  temp+=maps.btm_map[iy * ncols + ix];
	      }
	    temp/=scale*scale;
	    if (ext)
	      Z[y*nx+x]=1.70141e+038;
	    else		
	      Z[y*nx+x]=temp;		  
	  }

	if (0){
	for (int j=nrows-1;j>=0;j--){
	  printf("%dM: ",j);
	  for (int i=0;i<ncols;i++)
	    if (maps.btm_map[j*ncols+i]>=10e10) printf("X");
	    else
	      printf(" ");
	  printf("\n");
	}

	for (int j=nrows-1;j>=0;j--){
	  printf("%dZ: ",j);
	  for (int i=0;i<ncols;i++)
	    if (Z[j*ncols+i]>=10e10) printf("X");
	    else
	      printf(" ");
	  printf("\n");
	}
	}

	

	for (int i1=0;i1<g.punti_m.size();i1++)
	  if (g.punti_m[i1].w==1 && 
	      g.punti_m[i1].z==i+1
	      ) // e' tappo --> scrivo direttamente dentro Z
	    {
	      int x=((g.punti_m[i1].x-maps.minx_map)/maps.dx-1); // scalo bordo maps
	      int y=((g.punti_m[i1].y-maps.miny_map)/maps.dy-1);
	      // scalo rispetto a multirisoluzione corretta
	      x=(x)/(1<<((int)g.punti_m[i1].z-1))+1;
	      y=(y)/(1<<((int)g.punti_m[i1].z-1))+1;
	      // scrivo punto
	      printf("TAPPO manuale: %d %d livello %d\n",x,y,(int)g.punti_m[i1].z);
	      Z[y*nx+x]=1.70141e+038;
	    }

	printf("------------------\n");
      
	int i1=i+1; //test livello precedente
	int sizex=bsx/(1<<i1);
	int sizey=bsy/(1<<i1);
	int larg=(1<<i); // da scalare (tappato se almeno un tappo)

	if (i<levels-1){
	/// passa in rassegna divisione i-1, per ogni blocco adiacente, controlla 
	  
	  /*
		    for (int y2=ny-1;y2>=1;y2--){
		      printf(" %d: ",y2);
		      for (int x2=0;x2<nx-1;x2++){
			printf("%d ",Z[(y2+1) * nx + (x2+1)]>10e10);
		      }
		      printf("\n");
		    }
	  */
	  /// Z padding di 1 cella intorno
	  for (int y = 0; y <sizey; y ++) {  
	    for (int x = 0; x < sizex; x ++){
	      int dbg=0;//y*(1<<i1)>=157&&y*(1<<i1)<=159; 
	      //(i1==1 && x==30 && y==109) ||
	      //	(i1==2 && x==15 && y==54)
	      //	;
	      //		x*(1<<i1)*BLOCKSIZE_X <=998 && (x+1)*(1<<i1)*BLOCKSIZE_X >=988 &&
	      //	y*(1<<i1)*BLOCKSIZE_X <=3518 && (y+1)*(1<<i1)*BLOCKSIZE_X >=3510
		;//i==2;

		if (dbg){
		printf("%d %d: \n",x,y);
		for (int test1=0;test1<BLOCKSIZE_X*2;test1++){

		  for (int test2=0;test2<2*BLOCKSIZE_X*2;test2++){
		    if (test2%(BLOCKSIZE_X*2)==0) printf("  ");
		    printf("%d ",Z[(y*(BLOCKSIZE_X)*2+test1+1)* nx + (x*(BLOCKSIZE_X)*2+test2+1)]>10e10);		    
		  }
		  printf("\n");		  
		}
		}


	      if (bitmask[i1][sizex*y+x]==i1+1){ //risol precedente
		  //verifico bordi: se a fianco ho risol attuale -> controllo tappi
		int e=0,w=0,n=0,s=0;
		if (x<sizex-1 && bitmask[i1][sizex*y+x+1]==i1) //east risol piu' alta
		  e=1;
		if (x>0 && bitmask[i1][sizex*y+x-1]==i1 ) //w risol piu' alta
		  w=1;
		if (y<sizey-1 && bitmask[i1][sizex*(y+1)+x]==i1) //n risol piu' alta
		  n=1;
		if (y>0 && bitmask[i1][sizex*(y-1)+x]==i1) //s risol piu' alta
		  s=1;



		if (n+s+w+e>0){

		  int size=2*BLOCKSIZE_X; // lavoro su risoluzione maggiore!! (Z e' gia' caricata nella risoluzione piu' alta (i))

		  if (dbg) printf("bordo %d %d ris %d, %d %d %d %d\n",x,y,i1+1,n,s,w,e);
		  if (0 && dbg){
		    for (int y2=(y+1)*size;y2>=y*size;y2--){
		      printf(" %d: ",y2);
		      for (int x2=0;x2<nx-1;x2++){
			if (x2%size==0) printf("  ");
			printf("%d ",Z[(y2+1) * nx + (x2+1)]>10e10);
		      }
		      printf("\n");
		    }
		  }


		  // se ho 4 celle host_info dall'altra parte del bordo miste  --> tutto tappato su Z che sto costruendo

		  //scan di tutto il quadrato vicino ad alta (piu' affaccio del 2x2 nel quadrato a bassa)
		  

		  // se 2x2 (ad alta) del bordo del quadrato a bassa sono miste (a bassa risoluzione e' tappo) -> devo mettere tutti tappi ad alta,
		  // nel caso in cui dall'altra parte ci sia tutto libero, altrimenti non metto correttamente le condizioni di muro ad alta sull'affaccio

		  // se 2x2(ad alta) del bordo del quadrato senza tappi e dall'altra parte 2x2 misto -> tappo l'altra parte 
		  
		  // misto vs misto -> tappo tutto
		  
		  // celle in patch ad alta risoluzione (vicina)
		  for (int caso=0;caso<4;caso++){
		  int dx=0,dy=0; // cellle multires bitmask
		  if (caso==0 && s==1) dy=-1;
		  if (caso==1 && n==1) dy=+1;
		  if (caso==2 && w==1) dx=-1;
		  if (caso==3 && e==1) dx=+1;
		  if (dbg) printf("caso %d d %d %d, larg %d\n",caso,dx,dy,larg);



		  if (0){
		    for (int y2=(y+1+dy)*size-1;y2>=(y+dy)*size;y2--){
		      printf(" %d: ",y2);
		      for (int x2=(x+dx)*size;x2<(x+dx+1)*size;x2++){
			printf("%d ",Z[(y2+1) * nx + (x2+1)]>10e10);
		      }
		      printf("\n");
		    }
		  }

		  for (int x1=(x+dx)*size;x1<(x+dx+1)*size;x1+=2)
		    for (int y1=(y+dy)*size;y1<(y+dy+1)*size;y1+=2) // lavoro su livello corrente!!
		      if ( ((caso==0 && s==1 && y1>=(y+dy+1)*size-2)||
			    (caso==1 && n==1 && y1<(y+dy)*size+2)||
			    (caso==2 && w==1 && x1>=(x+dx+1)*size-2)||
			    (caso==3 && e==1 && x1<(x+dx)*size+2))){
			if (dbg) printf("Caso %d %d %d\n",caso,x1,y1);
			  // lavoro a quadrati 2x2, processo se quadrato 2x2 appartiene a bordo corretto adiacente a risol piu' alta
			      int ext=0;
			      for (int y2=y1+1;y2>=y1;y2--){
				if (dbg) printf("%d %2d ",y1,y2);
				for (int x2=x1;x2<x1+2;x2++){
				  int val=0; ///dentro
				  if (x2+1>=nx || y2+1>=ny)
				    val++;
				  else
				    if (Z[(y2+1) * nx + (x2+1)]>10e10)
				      val++;

				  if (dbg) printf("%d ",val);
				  if (val>0)
				    ext++;
				}
				if (dbg) printf("\n");
			      }
			      if (dbg) printf("\n %d\n",ext);
			      if (ext!=0 && ext!=4){ // basta mettere tutte a esterno
				if (dbg) printf("fix1 %d %d\n",x,y);								
				for (int y2=y1;y2<(y1+2);y2++)
				  for (int x2=x1;x2<(x1+2);x2++)
				    //setto Z 10e10
				    Z[(y2+1) * nx + (x2+1)]=1.70141e+038;
			      }
		      }
		  /*
		  printf("-->\n");
		  for (int y2=(y+1+dy)*size-1;y2>=(y+dy)*size;y2--){
		    printf(" %d: ",y2);
		    for (int x2=(x+dx)*size;x2<(x+dx+1)*size;x2++){
		      printf("%d ",Z[(y2+1) * nx + (x2+1)]>10e10);
		    }
		    printf("\n");
		  }
		  printf("\n");
		  printf("\n");
*/

		  }


		  if (0)
		    for (int y2=(y+1)*size-1;y2>=y*size;y2--){
		      printf(" %d: ",y2);
		      for (int x2=x*size;x2<(x+1)*size;x2++){
			printf("%d ",Z[(y2+1) * nx + (x2+1)]>10e10);
		      }
		      printf("\n");
		    }

		  //celle ad alta risoluzione nella patch a bassa (per eventuali muri nella patch a alta)
		  for (int caso=0;caso<4;caso++)
		  for (int x1=(x)*size;x1<(x+1)*size;x1+=2)
		    for (int y1=(y)*size;y1<(y+1)*size;y1+=2) // lavoro su livello corrente!!
		      if ( ((caso==0 && s==1 && y1<y*size+2)||
			    (caso==1 && n==1 && y1>=(y+1)*size-2)||
			    (caso==2 && w==1 && x1<x*size+2)||
			    (caso==3 && e==1 && x1>=(x+1)*size-2))){
			if (dbg) printf("Caso %d %d %d (patch ad alta)\n",caso,x1,y1);
			  // lavoro a quadrati 2x2, processo se quadrato 2x2 appartiene a bordo corretto adiacente a risol piu' alta
			int ext=0;

			for (int y2=y1-10;y2<y1+10;y2++){
			  if (dbg) printf("%2d %2d",x1-10,y2);
			  for (int x2=x1-10;x2<x1+10;x2++){
			    if (dbg) printf("%c ",Z[(y2+1) * nx + (x2+1)]>10e10?'X':'.');
			    if (dbg && (x2-1)%BLOCKSIZE_X==0) printf(" ");
			  }
			  if (dbg) printf("\n");
			  if (dbg && (y2-1)%BLOCKSIZE_Y==0) printf("\n");
			}

			if (dbg) printf("\n");
			for (int y2=y1+1;y2>=y1;y2--){
			  if (dbg) printf("%2d ",y2);
			  for (int x2=x1;x2<x1+2;x2++){
			    if (dbg) printf("%c ",Z[(y2+1) * nx + (x2+1)]>10e10?'X':'.');
			    if (x2+1>=nx || y2+1>=ny || Z[(y2+1) * nx + (x2+1)]>10e10)
			      ext++;
			  }
			  if (dbg) printf("\n");
			}
			if (dbg) printf("\n");
			if (ext!=0 && ext!=4){ // basta mettere tutte a esterno
			  if (dbg) printf("fix %d %d\n",x,y);
			  for (int x2=x1;x2<x1+2;x2++)
			    for (int y2=y1;y2<y1+2;y2++)
			      //setto Z 10e10
			      Z[(y2+1) * nx + (x2+1)]=1.70141e+038;
			}
		      }
		  
		  if (0){
		      printf("------>\n");
		    for (int y2=(y+1)*size-1;y2>=y*size;y2--){
		      printf(" %d: ",y2);
		      for (int x2=x*size;x2<(x+1)*size;x2++){
			printf("%d ",Z[(y2+1) * nx + (x2+1)]>10e10);
		      }
		      printf("\n");
		    }
		  }

		  if (dbg){
		    printf("%d %d: \n",x,y);
		    for (int test1=0;test1<BLOCKSIZE_X*2;test1++){
		      
		      for (int test2=0;test2<2*BLOCKSIZE_X*2;test2++){
			if (test2%(BLOCKSIZE_X*2)==0) printf("  ");
			printf("%d ",Z[(y*(BLOCKSIZE_X)*2+test1+1)* nx + (x*(BLOCKSIZE_X)*2+test2+1)]>10e10);		    
		      }
		      printf("\n");		  
		    }

		  }
		  
		 
		}
	      }



	    }
	  }

	  
	  /*
	    printf("dopo------------------\n");
	  for (int y2=ny-1;y2>=1;y2--){
	    printf(" %d: ",y2);
	    for (int x2=0;x2<nx-1;x2++){
	      printf("%d ",Z[(y2+1) * nx + (x2+1)]>10e10);
	    }
	    printf("\n");
	  }
	  */



	}

	if (0)	
	for (int j=159;j>=157;j--){
	  printf("%d: ",j);
	  for (int i=0;i<nx;i++)
	    if (Z[j*nx+i]>=10e10) printf("X");
	    else
	      printf(" ");
	  printf("\n");
	}

      
	readBLNBCC(path,i,Z,i+1==g.global_level_bcc,bitmask[i],sizex);


	/*	for (int j=47*BLOCKSIZE_X;j>=46*BLOCKSIZE_X;j--){
	  printf("%d: ",j);
	  for (int i=11*BLOCKSIZE_X;i<12*BLOCKSIZE_X;i++)
	    printf("%3d ",maps.host_info[j*ncols+i].x);
	  printf("\n");
	}
	for (int j=47*BLOCKSIZE_X;j>=46*BLOCKSIZE_X;j--){
	  printf("%d: ",j);
	  for (int i=11*BLOCKSIZE_X;i<12*BLOCKSIZE_X;i++)
	    if (Z[j*ncols+i]>=10e10) printf("X");
            else
              printf(" ");
	  printf("\n");
	}
	for (int j=47*BLOCKSIZE_X;j>=46*BLOCKSIZE_X;j--){
	  printf("%d: ",j);
	  for (int i=11*BLOCKSIZE_X;i<12*BLOCKSIZE_X;i++)
	    if (maps.btm_map[j*ncols+i]>=10e10) printf("X");
            else
              printf(" ");
	  printf("\n");
	}
	*/


	if (0){
	for (int j=nrows-1;j>=0;j--){
	  printf("%dM: ",j);
	  for (int i=0;i<ncols;i++)
	    if (maps.btm_map[j*ncols+i]>=10e10) printf("X");
	    else
	      printf(" ");
	  printf("\n");
	}

	for (int j=nrows-1;j>=0;j--){
	  printf("%dZ: ",j);
	  for (int i=0;i<ncols;i++)
	    if (Z[j*ncols+i]>=10e10) printf("X");
	    else
	      printf(" ");
	  printf("\n");
	}
	}


	maps.host_info_x_m[i]=maps.host_info_x;
	maps.host_info_y_m[i]=maps.host_info_y;
	if (i!=g.global_level_bcc-1) // gia' allocato
	  maps.host_info_m[i] = (uchar4 *) malloc( maps.host_info_x *maps.host_info_y * sizeof(uchar4));//sides
	//copia valori
	printf("copy %d x %d\n",maps.host_info_x,maps.host_info_y);
	for (int x=0;x<maps.host_info_x_m[i];x++)
	  for (int y=0;y<maps.host_info_y_m[i];y++){
	    maps.host_info_m[i][y*maps.host_info_x_m[i]+x]=maps.host_info[y*maps.host_info_x_m[i]+x];
	  }	

	int dbgx=15;
	int dbgy=54;
	if (0 && i==1)
		{
		printf("AFTER %d %d: \n",dbgx,dbgy);
		for (int test1=0;test1<BLOCKSIZE_X*2;test1++){

                  for (int test2=0;test2<2*BLOCKSIZE_X*2;test2++){
                    if (test2%(BLOCKSIZE_X*2)==0) printf("  ");
                    printf("%d ",Z[(dbgy*(BLOCKSIZE_X)*2+test1+1)* nx + (dbgx*(BLOCKSIZE_X)*2+test2+1)]>10e10);
                  }
		  printf("\n");
                }


		for (int test1=0;test1<BLOCKSIZE_X*2;test1++){

		  for (int test2=0;test2<2*BLOCKSIZE_X*2;test2++){
		    if (test2%(BLOCKSIZE_X*2)==0) printf("  ");

		    unsigned char val=maps.host_info[(dbgy*(BLOCKSIZE_X)*2+test1+1)* nx + (dbgx*(BLOCKSIZE_X)*2+test2+1)].x;
		    if (val==0)
		      printf("  .");
		    else{
		      if (val==BIT_EXTERN)
			printf("  X");
		      else
			printf("%3d",val);
		    }

		  }
		  printf("\n");		  
		}
		}

      }


      /// ora codifico blocchi

if(0)
  for (int i=0;i<levels;i++)
	if (i==2)
	{
	int dbg=1;
	int sizex=bsx/(1<<i);
	int sizey=bsy/(1<<i);
	for (int y = 0; y <sizey*BLOCKSIZE_X; y ++)  {
	  printf("%d: ",y);
	  for (int x = 0; x < sizex*BLOCKSIZE_X; x ++){
	    if (x<maps.host_info_y_m[i])
	      printf("%d ",
		     maps.host_info_m[i][y*maps.host_info_x_m[i]+x].x==15?
		     9:maps.host_info_m[i][y*maps.host_info_x_m[i]+x].x);
	  }
	  printf("\n");
	}
	}

      // per ogni blocco, esclude se tutte celle sono fuori
      if (1)
      for (int i=0;i<levels;i++){
	int dbg=0;
	int sizex=bsx/(1<<i);
	int sizey=bsy/(1<<i);
	for (int y = 0; y <sizey; y ++)  
	  for (int x = 0; x < sizex; x ++)
	    if (bitmask[i][sizex*y+x]==i+1){
	      int dbg=0;//y>=sizey-2;
	      if (dbg)printf("liv %d, xy %d %d\n",i,x,y);
	      int fuori=1;
	      int px,py;
	      for (int y1=BLOCKSIZE_X-1;y1>=0;y1--){
		py=((y*BLOCKSIZE_Y+y1)+1 );
		if (dbg) printf("%d: ",py);
		for (int x1=0;x1<BLOCKSIZE_X;x1++){
		  px=((x*BLOCKSIZE_X+x1)+1 );

		  if (dbg){
		    if (py >= maps.host_info_y_m[i]-1 ||
			px >= maps.host_info_x_m[i]-1)
			printf("F ");
		    else
		      printf("%d ",
			     maps.host_info_m[i][py*maps.host_info_x_m[i]+px].x==15?9:
			     maps.host_info_m[i][py*maps.host_info_x_m[i]+px].x);
		  }

		  if (py < maps.host_info_y_m[i]-1 &&
		      px < maps.host_info_x_m[i]-1 && 
		      maps.host_info_m[i][py*maps.host_info_x_m[i]+px].x!=BIT_EXTERN)
		    fuori=0;
		}
		if (dbg) printf("\n");
	      }
	      if (fuori)
		bitmask[i][sizex*y+x]=-1; // cancello	      
	    }
      }

      printf("go numerazione\n");





      //numerazione blocchi

      // stats: quante celle per livello
      // primo pass: assegno codici blocchi con contorno
      int tot_blocks=0;
      int bound_blocks=0;
      int codes=0;

      for (int i=0;i<levels;i++){
	int dbg=1;
	int sizex=bsx/(1<<i);
	int sizey=bsy/(1<<i);
	int ct=0;
	for (int y = sizey-1; y >=0; y --){  
	  if (dbg) printf("%d: ",y);
	  for (int x = 0; x < sizex; x ++){
	    ct+=(bitmask[i][sizex*y+x]==i+1?1:0);	    
	    if (bitmask[i][sizex*y+x]==-1) {
	      if (dbg) printf("! ");
	    }
	    else{
	      int ok=1;
	      if (bitmask[i][sizex*y+x]==i+1){
		//verifico se c'e' almeno un valore != 0
		int size=BLOCKSIZE_X;
		for (int x1=x*size;x1<(x+1)*size;x1++)
		  for (int y1=y*size;y1<(y+1)*size;y1++)
		    if (x1+1<maps.host_info_x_m[i]-1 && y1+1<maps.host_info_y_m[i]-1){ // salto bordo esterno
		      if (maps.host_info_m[i][(y1+1) * maps.host_info_x_m[i] + (x1+1)].x!=0 &&
			  maps.host_info_m[i][(y1+1) * maps.host_info_x_m[i] + (x1+1)].x!=BIT_EXTERN) // ricordati lo shift!!
			ok=0;
		    }
	      }
	      
	      if (ok || bitmask[i][sizex*y+x]!=i+1){
		if (dbg) printf("%d ",bitmask[i][sizex*y+x]);//<=i+1?bitmask[i][sizex*y+x]:0));
	      }
	      else{
		if (dbg) printf("%dB",bitmask[i][sizex*y+x]);
		bitmaskC[i][sizex*y+x]=codes++;
		bound_blocks++;
	      }
	    }
	  }
	  if (dbg) printf("\n");
	}
	tot_blocks+=ct;
	printf("Level %d: ratio %f (%d cells)\n",i,(ct+0.0)/(sizex*sizey),ct);
      }
      printf("Celle originali %d, celle multires %d\n",ncols*nrows,tot_blocks*BLOCKSIZE_X*BLOCKSIZE_Y);
      printf("Compressione: %2.1fx\n",1.0/((0.0+tot_blocks)*BLOCKSIZE_X*BLOCKSIZE_Y/ncols/nrows));

      // assegno codici di blocchi: 2nd pass
      // prima codici contorno e poi resto
      for (int i=0;i<levels;i++){
	int dbg=1;
	int sizex=bsx/(1<<i);
	int sizey=bsy/(1<<i);

	for (int y = sizey-1; y >=0; y --){  
	  for (int x = 0; x < sizex; x ++){
	    if (bitmask[i][sizex*y+x]==i+1 &&
		bitmaskC[i][sizex*y+x]==-1){ // se non gia' assegnato (contorni)
	      bitmaskC[i][sizex*y+x]=codes++;
	    }
	    if (dbg){
	      if (bitmaskC[i][sizex*y+x]==-1) 
		printf("... ");
	      else
		printf("%3d ",bitmaskC[i][sizex*y+x]);//==i+1?i+1:0));
	    }
	  }
	  if (dbg) printf("\n");
	}
      }
      printf("Celle %d, bound %d\n",tot_blocks,bound_blocks);
      /// setup array celle e host_info_multi
      maps.tot_blocks=tot_blocks;

      /// host_info --> sono i primi bound_blocks blocchi (cosi' passo direttamente a kernel senza ulteriore numerazione)
int x_blocks=1<<(int)(ceil(log((float)maps.tot_blocks)/log((float)2)/(float)2)); 
      int y_blocks=(tot_blocks-1)/x_blocks+1;
      printf ("blocks %d alloc: %d x %d array \n",tot_blocks,x_blocks,y_blocks);

      maps.host_info_multi=(uchar4*)malloc(x_blocks*y_blocks*BLOCKSIZE_X*BLOCKSIZE_Y*sizeof(uchar4));

      // resetto (soprattutto la parte di non bound -> deve essere a 0 per compatibilita')
      for (int b=tot_blocks;b<x_blocks*y_blocks;b++){
	for (int y1=0;y1<BLOCKSIZE_Y;y1++){
	  for (int x1=0;x1<BLOCKSIZE_X;x1++) {
	    int x=b%x_blocks;
	    int y=b/x_blocks;
	    int idx=(y*BLOCKSIZE_Y+y1)*x_blocks*BLOCKSIZE_X+BLOCKSIZE_X*x+x1;
	    maps.host_info_multi[idx].x=BIT_EXTERN;
	    maps.host_info_multi[idx].y=0;
	    maps.host_info_multi[idx].z=0;
	    maps.host_info_multi[idx].w=0;
	  }
	}
      }
      for (int b=0;b<tot_blocks;b++){
	for (int y1=0;y1<BLOCKSIZE_Y;y1++){
	  for (int x1=0;x1<BLOCKSIZE_X;x1++) {
	    int x=b%x_blocks;
	    int y=b/x_blocks;
	    int idx=(y*BLOCKSIZE_Y+y1)*x_blocks*BLOCKSIZE_X+BLOCKSIZE_X*x+x1;
	    maps.host_info_multi[idx].x=0;
	    maps.host_info_multi[idx].y=0;
	    maps.host_info_multi[idx].z=0;
	    maps.host_info_multi[idx].w=0;
	  }
	}
      }

      maps.host_grid_multi=(F4*)malloc(x_blocks*y_blocks*BLOCKSIZE_X*BLOCKSIZE_Y*sizeof(F4));
      maps.host_man_map_multi=(F*)malloc(x_blocks*y_blocks*BLOCKSIZE_X*BLOCKSIZE_Y*sizeof(F));
      maps.host_grid_level_multi=(unsigned char*)malloc(tot_blocks*sizeof(unsigned char));
      maps.host_ofs_blocks=(ushort2*)malloc(tot_blocks*sizeof(ushort2));

      if (g.global_por){
	maps.host_hlx_map_multi=(F*)malloc(x_blocks*y_blocks*BLOCKSIZE_X*BLOCKSIZE_Y*sizeof(F));
	maps.host_hly_map_multi=(F*)malloc(x_blocks*y_blocks*BLOCKSIZE_X*BLOCKSIZE_Y*sizeof(F));
	maps.host_por_map_multi=(F*)malloc(x_blocks*y_blocks*BLOCKSIZE_X*BLOCKSIZE_Y*sizeof(F));
      }
      //host_grid_flxL_multi=(F4*)malloc(tot_blocks*BLOCKSIZE_X*BLOCKSIZE_Y*sizeof(F4));
      //host_grid_flxR_multi=(F4*)malloc(tot_blocks*BLOCKSIZE_X*BLOCKSIZE_Y*sizeof(F4));
      //host_grid_flyL_multi=(F4*)malloc(tot_blocks*BLOCKSIZE_X*BLOCKSIZE_Y*sizeof(F4));
      //host_grid_flyR_multi=(F4*)malloc(tot_blocks*BLOCKSIZE_X*BLOCKSIZE_Y*sizeof(F4));

      for (int i=0;i<levels;i++) {
	int sizex=bsx/(1<<i);
	int sizey=bsy/(1<<i);
	for (int y = 0; y < sizey; y ++)
	  for (int x = 0; x < sizex; x ++)
	    if (bitmaskC[i][sizex*y+x]>=0){
	      /// memorizza livello
	      maps.host_grid_level_multi[bitmaskC[i][sizex*y+x]]=i;
	      //carico coordinate
	      maps.host_ofs_blocks[bitmaskC[i][sizex*y+x]].x=x*(1<<i);
	      maps.host_ofs_blocks[bitmaskC[i][sizex*y+x]].y=y*(1<<i);

	      int dbg=0;//bitmaskC[i][sizex*y+x]==3 ;//i==0;
	      if (dbg) printf("%d xy %d %d bl:%d\n",i,x,y,bitmaskC[i][sizex*y+x]);
	      int idx=bitmaskC[i][sizex*y+x]*BLOCKSIZE_X; // base del quadrato
	      for (int x1=0;x1<BLOCKSIZE_X;x1++)
		for (int y1=0;y1<BLOCKSIZE_X;y1++){
		  int dbg=0;//bitmaskC[i][sizex*y+x]==3;
		  F h=0;
		  F z=0;
		  F vvx=0;
		  F vvy=0;
		  F man=0;
		  F hlx=0;
		  F hly=0;
		  F por=0;
		  int cth=0;
		  int cthh=0;
		  int tappo=0;
		  int bit=0;
		  if (maps.host_info_m[i][(y*BLOCKSIZE_Y+y1 +1)*maps.host_info_x_m[i]+(x*BLOCKSIZE_X+x1 +1)].x==BIT_EXTERN)
		    tappo=1;
		  if (tappo==0)
		  for (int x2=0;x2<(1<<i);x2++)  // calcola altezza media nei livelli di multires
		    for (int y2=0;y2<(1<<i);y2++){
		      int px,py;
		      px=((x*BLOCKSIZE_X+x1)*(1<<i)+x2 +1 );
		      py=((y*BLOCKSIZE_Y+y1)*(1<<i)+y2 +1 );
		      if (py >= nrows-1 ||
			  px >= ncols-1){
			if (dbg) printf("xy1: %d %d --> pxy %d %d TAPPO\n",x1,y1,px,py);
			tappo=1;
		      }
		      else {
			if (dbg) printf("xy1: %d %d --> pxy %d %d:",x1,y1,px,py);
			if (dbg)  printf(" (%dx%d)\n",ncols,nrows);
			if (maps.inh_map[py*ncols+px]-maps.btm_map[py*ncols+px]>=g.global_YEPS){
			  h+=maps.inh_map[py*ncols+px]; // +1 perche' salto bordo padding
			  vvx+=maps.vvx_map[py*ncols+px]; // +1 perche' salto bordo padding
			  vvy+=maps.vvy_map[py*ncols+px]; // +1 perche' salto bordo padding
			  cthh++;
			}
			bit=maps.host_info_m[0][py*maps.host_info_x_m[0]+px].x;
			z+=maps.btm_map[py*ncols+px]; // +1 perche' salto bordo padding
			man+=maps.man_map[py*ncols+px]; // +1 perche' salto bordo padding
			if (g.global_por){
			  por+=maps.por_map[py*ncols+px]; // +1 perche' salto bordo padding
			  hlx+=maps.hlx_map[py*ncols+px]; // +1 perche' salto bordo padding
			  hly+=maps.hly_map[py*ncols+px]; // +1 perche' salto bordo padding
			}
			cth++;
		      }
		    }		    
		  if (tappo){
		    h=0;
		    z=1.70141e+038;
		  }
		  else{
		    if (cthh>0){
		      h/=cthh;
		      vvx/=cthh;
		      vvy/=cthh;
		    }
		    z/=cth;
		    man/=cth;
		    if (g.global_por){
		      por/=cth;
		    hlx/=cth;
		    hly/=cth;
		    }
		  }
		  if (h<z) h=z;
		  if (dbg) printf("%g %g, bit %d\n",h,z,bit);
		  int b=bitmaskC[i][sizex*y+x];
		  int x=b%x_blocks;
		  int y=b/x_blocks;

		  int idx1=(y*BLOCKSIZE_Y+y1)*x_blocks*BLOCKSIZE_X+BLOCKSIZE_X*x+x1;
		  if (dbg) printf("idx %d,  ",idx1);
		  if (dbg)  printf("Block %d, xb%d %d, %d %d\n",i,x_blocks,y_blocks,x,y);

		  maps.host_grid_multi[idx1].x=h;
		  maps.host_grid_multi[idx1].y=vvx;
		  maps.host_grid_multi[idx1].z=vvy;
		  maps.host_grid_multi[idx1].w=z;
		  maps.host_man_map_multi[idx1]=man;
		  if (g.global_por){
		  maps.host_por_map_multi[idx1]=por;
		  maps.host_hlx_map_multi[idx1]=hlx;
		  maps.host_hly_map_multi[idx1]=hly;
		  }
		}
	    }
      }

      /*
	output debug: posizione x y res di ogni cella
       */
      char path1[256];
      strcpy(path1,g.base_path);strcat(path1,".PMAP");   
      printf("Write %s, with points coordinates\n",path1);
      FILE* fo=fopen(path1,"w+");
      for (int i=0;i<tot_blocks;i++){
	int lev=maps.host_grid_level_multi[i];
	//carico coordinate
	int x=maps.host_ofs_blocks[i].x; // gia' in coordinate ris. massima
	int y=maps.host_ofs_blocks[i].y;
	//printf("bl %d: %d %d\n",i,x,y);
	for (int px=0;px<BLOCKSIZE_X;px++) 
	  for (int py=0;py<BLOCKSIZE_Y;py++){
	    double ofsx=maps.dx/2.0*(1<<lev); // baricentro cella (anche con multires) //PROVA!!
	    double ofsy=maps.dy/2.0*(1<<lev);
	    fprintf(fo,"%f %f %d\n",
		   maps.minx_map+maps.dx*BLOCKSIZE_X*x+maps.dx*px*(1<<lev)+ofsx,
		   maps.miny_map+maps.dy*BLOCKSIZE_Y*y+maps.dy*py*(1<<lev)+ofsy,
		   lev);
	       }
      }
      fclose(fo);



      //copio grid da multires
      for (int i=0;i<levels;i++) {
	int sizex=bsx/(1<<i);
	int sizey=bsy/(1<<i);
	for (int y = 0; y < sizey; y ++)
	  for (int x = 0; x < sizex; x ++)
	    if (bitmaskC[i][sizex*y+x]>=0 &&
		bitmaskC[i][sizex*y+x]<bound_blocks){
	      int idx=bitmaskC[i][sizex*y+x]*BLOCKSIZE_X; // base del quadrato
	      int b=bitmaskC[i][sizex*y+x];
	      int xb=b%x_blocks;
	      int yb=b/x_blocks;
	      
	      for (int x1=0;x1<BLOCKSIZE_X;x1++)
		for (int y1=0;y1<BLOCKSIZE_X;y1++){
		  int idx1=(yb*BLOCKSIZE_Y+y1)*x_blocks*BLOCKSIZE_X+BLOCKSIZE_X*xb+x1;
		  if ((x*BLOCKSIZE_X+x1 +1)<maps.host_info_x_m[i]-1 &&
		      (y*BLOCKSIZE_Y+y1 +1)<maps.host_info_y_m[i]-1)		      
		    maps.host_info_multi[idx1]=
		      maps.host_info_m[i][(y*BLOCKSIZE_Y+y1 +1)*maps.host_info_x_m[i]+(x*BLOCKSIZE_X+x1 +1)]; // +1 perche' salto bordo padding;
		  else
		    maps.host_info_multi[idx1].x=BIT_EXTERN;
		}      
	    }      
      }      

      // esploro vicini
      // formato
      // codice spostamento per trovare l'amico: 0 (-1,-1), 1 (0 -1), 2 (1 -1), 3(-1 0), 4 (1 0), 5 (-1 1), 6 (0 1), 7 (1 1)
      // livello di multires: -1, 0, 1 (delta per raggiungere il livello amico)
      // codice vicino n1 (-1 muro)
      // info aggiuntiva n2: se livello  1 -->  quadrante del quadrato amico (codifica x+2*y, 0 (0 0) 1(1 0) 2 (0 1) 3 (1 1))
      //                     se livello -1 e amico non e' corner -->  secondo quadrato adiacente

      maps.neigh=(neigh_t*)malloc(8*tot_blocks*sizeof(neigh_t));

      for (int i=0;i<levels;i++) {
	int dbg=0;
	int sizex=bsx/(1<<i);
	int sizey=bsy/(1<<i);
	for (int y = 0; y < sizey; y ++)
	  for (int x = 0; x < sizex; x ++)
	    if (bitmask[i][sizex*y+x]==i+1){
	      int side=0;
	      for (int dy=-1;dy<=1;dy++)
	      for (int dx=-1;dx<=1;dx++)
		if (dx!=0 || dy!=0){
		  int nx=x+dx;
		  int ny=y+dy;
		  int idx=8*bitmaskC[i][sizex*y+x]+side;
		  if (nx<0 || nx>=sizex ||
		      ny<0 || ny>=sizey ||
		      bitmask[i][sizex*ny+nx]==-1){ // quadrato fuori
		    maps.neigh[idx].lev=0;
		    maps.neigh[idx].n1=-1;
		    maps.neigh[idx].n2=-1;
		    if (dbg) 
		      printf("%d %d %d (%d %d), %d: %d %d [%d, %d]\n",i,x,y,dx,dy,bitmaskC[i][sizex*y+x],side,0,-1,-1); // muro
		  }
		  else{
		    if (bitmask[i][sizex*ny+nx]==i+2){ //multires dopo, uso secondo posto per indicare se 0 (prima meta') o 1 (seconda meta' di adiacenza)
		      int bit=-1;
		      bit=(nx%2)+2*(ny%2);
		      maps.neigh[idx].lev=1;
		      maps.neigh[idx].n1=bitmaskC[i+1][sizex/2*(ny/2)+(nx/2)];
		      maps.neigh[idx].n2=bit;
		      if (dbg)
		      printf("%d %d %d (%d %d), %d: %d %d [%d, %d]\n",i,x,y,dx,dy,bitmaskC[i][sizex*y+x],
			     side,
			     1,
			     bitmaskC[i+1][sizex/2*(ny/2)+(nx/2)],
			     bit);
		    }
		    else{
		      if (i==0){ // primo livello non scendo, non ci sono ulteriori livelli
		      maps.neigh[idx].lev=0;
		      maps.neigh[idx].n1=bitmaskC[i][sizex*(ny)+(nx)];
		      maps.neigh[idx].n2=-1; // non usata
			if (dbg)
			  printf("%d %d %d (%d %d), %d: %d %d %d\n",i,x,y,dx,dy,bitmaskC[i][sizex*y+x],side,0,bitmaskC[i][sizex*(ny)+(nx)]);
		      }
		      else{ // guardo il livello precedente
			int ex=x*2;
			int ey=y*2;
			if (dx<0) ex--;// adiacente sx, scavalco quadrato
			if (dy<0) ey--;
			if (dx>0) ex+=2; // adiacente destro, +1 per vicino nello stesso quad (a ris i) +1 scavalco quadrato
			if (dy>0) ey+=2; // adiacente destro, +1 per vicino nello stesso quad (a ris i) +1 scavalco quadrato
			int ex1,ey1; // secondario (usato se si va a livello -1
			ex1=ex;
			ey1=ey;
			if (dx==0) ex1++; // sopra o sotto, uso anche cella adiacente a destra
			if (dy==0) ey1++; // sx o dx, uso anche cella adiacente +1			  
			if (bitmask[i-1][sizex*2*(ey)+(ex)]==i+1){ // mio livello (garantisco che era tutto il quadrato dello stesso livello)
			  maps.neigh[idx].lev=0;
			  maps.neigh[idx].n1=bitmaskC[i][sizex*(ny)+(nx)];
			  maps.neigh[idx].n2=-1; // non usata			  
			  if (dbg)
			    printf("%d %d %d (%d %d), %d: %d %d [%d, %d]\n",i,x,y,dx,dy,bitmaskC[i][sizex*y+x],side,0,bitmaskC[i][sizex*(ny)+(nx)],-1);
			}
			else{
			  int bit=-1;
			  if (dx*dy==0) // quando non diagonale, metto il secondo quadrato
			    bit=bitmaskC[i-1][sizex*2*(ey1)+(ex1)];
			  maps.neigh[idx].lev=-1;
			  maps.neigh[idx].n1=bitmaskC[i-1][sizex*2*(ey)+(ex)];
			  maps.neigh[idx].n2=bit; 
			  if (dbg)
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



      //completo informazione host_info per blocchi adiacenti mancanti
      //se manca blocco -> aggiungo muro sulla cella adiacente
      //versione incrementale, eventuale riassegnamento codici e .y z w
      //in questo modo evito di chiedere info blocco vicino se c'e' muro
      int size=BLOCKSIZE_X;
      for (int i=0;i<tot_blocks;i++){
	int dbg=i==23;
	for (int y1=0;y1<size;y1++)
	  for (int x1=0;x1<size;x1++) {

	    int x=i%x_blocks;
	    int y=i/x_blocks;
	    int idx=(y*BLOCKSIZE_Y+y1)*x_blocks*BLOCKSIZE_X+BLOCKSIZE_X*x+x1;

	    unsigned char val=maps.host_info_multi[idx].x;
      // codice spostamento per trovare l'amico: 0 (-1,-1), 1 (0 -1), 2 (1 -1), 3(-1 0), 4 (1 0), 5 (-1 1), 6 (0 1), 7 (1 1)

	    if (maps.neigh[8*i+ 3 ].n1==-1) // fisso west
	      if (x1==0 && val!=BIT_EXTERN){ // check west
		if (val==0){
		  maps.host_info_multi[idx].x=BIT_W;
		  maps.host_info_multi[idx].y=1; // muro
		  maps.host_info_multi[idx].z=0;
		}
		if (val>>4 !=0){ //ci sono gia' due condizioni!! controlla che ci sia W, altrimenti err
		  if ((val>>4)!=BIT_W &&  //non in seconda
		      (val & 15)!=BIT_W){ //non in prima

		    int coef=1<<maps.host_grid_level_multi[i];	  
		    int ofsx=maps.host_ofs_blocks[i].x;
		    int ofsy=maps.host_ofs_blocks[i].y;
		    printf("ERROR W: 3 conditions around block %d, position%d %d\n",i,
			   BLOCKSIZE_X*ofsx+coef*x1,
			   BLOCKSIZE_X*ofsy+coef*y1);
		    for (int i1=0;i1<8;i1++)
		      printf("%d %d: lev %d, %d %d\n",i,i1,maps.neigh[8*i+i1].lev,maps.neigh[8*i+i1].n1,maps.neigh[8*i+i1].n2);
		    dbg_block(i,bound_blocks,tot_blocks);
		    exit(2);
		  }
		} 
		if(val!=0 && val>>4 ==0){ // se c'e' gia' W, lascio (bound c hanno precedenza)
		  if (val!=BIT_W){ //devo aggiungere west
		    if (val==BIT_E){ // shift e somma
		      maps.host_info_multi[idx].x=(val<<4)+BIT_W;
		      maps.host_info_multi[idx].y=(maps.host_info_multi[idx].y<<4)+1;
		      maps.host_info_multi[idx].z=0;
		      maps.host_info_multi[idx].w=maps.host_info_multi[idx].z;
		    }
		    else{ // aggiungo
		      maps.host_info_multi[idx].x+=BIT_W<<4;
		      maps.host_info_multi[idx].y+=1<<16;
		      maps.host_info_multi[idx].w=0;  
		    }
		  }
		}
	      }
	    val=maps.host_info_multi[idx].x;
	    if (maps.neigh[8*i+ 4 ].n1==-1) // fisso EAST
	      if (x1==BLOCKSIZE_X-1 && val!=BIT_EXTERN){ // check east
		if (dbg) printf("E %d %d %d\n",x1,y1,val);
		if (val==0){
		  if (dbg) printf("E1\n");
		  maps.host_info_multi[idx].x=BIT_E;
		  maps.host_info_multi[idx].y=1; // muro
		  maps.host_info_multi[idx].z=0;
		} 
		if (val>>4 !=0){ //ci sono gia' due condizioni!! controlla che la seconda sia E, altrimenti err
		  if (val>>4!=BIT_E){ //E e' sempre nella seconda pos
		    int coef=1<<maps.host_grid_level_multi[i];	  
		    int ofsx=maps.host_ofs_blocks[i].x;
		    int ofsy=maps.host_ofs_blocks[i].y;
		    printf("ERROR E: 3 conditions around block %d, position%d %d\n",i,
			   BLOCKSIZE_X*ofsx+coef*x1,
			   BLOCKSIZE_X*ofsy+coef*y1);
		    for (int i1=0;i1<8;i1++)
		      printf("%d %d: lev %d, %d %d\n",i,i1,maps.neigh[8*i+i1].lev,maps.neigh[8*i+i1].n1,maps.neigh[8*i+i1].n2);
		    dbg_block(i,bound_blocks,tot_blocks);
		    exit(2);
		  }
		}
		if(val!=0 && val>>4 ==0){ // se c'e' gia' E, lascio (bound c hanno precedenza)
		  if (dbg) printf("E2\n");
		  if (val!=BIT_E){ //devo aggiungere east
		    // aggiungo
		      maps.host_info_multi[idx].x+=BIT_E<<4;
		      maps.host_info_multi[idx].y+=1<<16;
		      maps.host_info_multi[idx].w=0;  
		  }
		}
	      }

	    val=maps.host_info_multi[idx].x;
	    if (maps.neigh[8*i+ 6 ].n1==-1) // fisso N
	      if (y1==BLOCKSIZE_Y-1 && val!=BIT_EXTERN){ // check north
		if (dbg) printf("N %d %d %d\n",x1,y1,val);
		if (val==0){
		  if (dbg) printf("N0\n");
		  maps.host_info_multi[idx].x=BIT_N;
		  maps.host_info_multi[idx].y=1; // muro
		  maps.host_info_multi[idx].z=0;
		}
		if (val>>4 !=0){ //ci sono gia' due condizioni!! controlla che ci sia N, altrimenti err
		  if (dbg) printf("N1\n");
		  if ((val>>4)!=BIT_N &&  //non in seconda
		      (val & 15)!=BIT_N){ //non in prima
		    int coef=1<<maps.host_grid_level_multi[i];	  
		    int ofsx=maps.host_ofs_blocks[i].x;
		    int ofsy=maps.host_ofs_blocks[i].y;
		    printf("ERROR N: 3 conditions around block %d (%d %d), position %d %d\n",i,val & 15, val>>4,
			   BLOCKSIZE_X*ofsx+coef*x1,
			   BLOCKSIZE_X*ofsy+coef*y1);
		    for (int i1=0;i1<8;i1++)
		      printf("%d %d: lev %d, %d %d\n",i,i1,maps.neigh[8*i+i1].lev,maps.neigh[8*i+i1].n1,maps.neigh[8*i+i1].n2);
		    dbg_block(i,bound_blocks,tot_blocks);
		    exit(2);
		  }
		}
		if(val!=0 && val>>4 ==0){ // se c'e' gia' N, lascio (bound c hanno precedenza)
		  if (dbg) printf("N2\n");
		  if (val!=BIT_N){ //devo aggiungere north (sempre in prima pos
		    // shift e somma
		    maps.host_info_multi[idx].x=(val<<4)+BIT_N;
		    maps.host_info_multi[idx].y=(maps.host_info_multi[idx].y<<4)+1;
		    maps.host_info_multi[idx].z=0;
		    maps.host_info_multi[idx].w=maps.host_info_multi[idx].z;
		  }
		}
	      }


	    val=maps.host_info_multi[idx].x;
	    if (maps.neigh[8*i+ 1 ].n1==-1) // fisso S
	      if (y1==0 && val!=BIT_EXTERN){ // check south
		if (val==0){
		  maps.host_info_multi[idx].x=BIT_S;
		  maps.host_info_multi[idx].y=1; // muro
		  maps.host_info_multi[idx].z=0;
		}
		if (val>>4 !=0){ //ci sono gia' due condizioni!! controlla che ci sia S, altrimenti err
		  if ((val>>4)!=BIT_S &&  //non in seconda
		      (val & 15)!=BIT_S){ //non in prima
		    int coef=1<<maps.host_grid_level_multi[i];	  
		    int ofsx=maps.host_ofs_blocks[i].x;
		    int ofsy=maps.host_ofs_blocks[i].y;
		    printf("ERROR S: 3 conditions around block %d (%d %d), position %d %d\n",i,val & 15, val>>4,
			   BLOCKSIZE_X*ofsx+coef*x1,
			   BLOCKSIZE_X*ofsy+coef*y1);
		    for (int i1=0;i1<8;i1++)
		      printf("%d %d: lev %d, %d %d\n",i,i1,maps.neigh[8*i+i1].lev,maps.neigh[8*i+i1].n1,maps.neigh[8*i+i1].n2);
		    dbg_block(i,bound_blocks,tot_blocks);
		    exit(2);
		  }
		}
		if(val!=0 && val>>4 ==0){ // se c'e' gia' W, lascio (bound c hanno precedenza)
		  if (val!=BIT_S){ //devo aggiungere west
		    if (val!=BIT_N){ // shift e somma
		      maps.host_info_multi[idx].x=(val<<4)+BIT_S;
		      maps.host_info_multi[idx].y=(maps.host_info_multi[idx].y<<4)+1;
		      maps.host_info_multi[idx].z=0;
		      maps.host_info_multi[idx].w=maps.host_info_multi[idx].z;
		    }
		    else{ // aggiungo
		      maps.host_info_multi[idx].x+=BIT_W<<4;
		      maps.host_info_multi[idx].y+=1<<16;
		      maps.host_info_multi[idx].w=0;  
		    }
		  }
		}
	      }


	  }
      }
      
      /*
	dbg_block(96,bound_blocks,tot_blocks);
	dbg_block(272,bound_blocks,tot_blocks);

	for (int j=47*BLOCKSIZE_X;j>=46*BLOCKSIZE_X;j--){
	  printf("%d: ",j);
	  for (int i=11*BLOCKSIZE_X;i<12*BLOCKSIZE_X;i++)
	    printf("%3d ",maps.host_info[j*ncols+i].x);
	  printf("\n");
	}
      */

      /*
      dbg_block(0,bound_blocks,tot_blocks);
      int x_blocks1=1<<(int)(ceil(log((float)maps.tot_blocks)/log((float)2)/(float)2));
      
      int idx=(208/BLOCKSIZE_X)+x_blocks1*(122/BLOCKSIZE_Y);

      dbg_block(idx,bound_blocks,tot_blocks);

      int coef1=1<<maps.host_grid_level_multi[idx];	  
      int ofsx1=maps.host_ofs_blocks[idx].x;
      int ofsy1=maps.host_ofs_blocks[idx].y;
      printf(" block %d, position%d %d\n",idx,
	     BLOCKSIZE_X*ofsx1+coef1*208%BLOCKSIZE_X,
	     BLOCKSIZE_X*ofsy1+coef1*122%BLOCKSIZE_Y);
      
      /// debug
      if (1){
	int dbg_i=idx;
	for (int i=0;i<8;i++){
	  printf("%d %d: lev %d, %d %d\n",dbg_i,i,maps.neigh[8*dbg_i+i].lev,maps.neigh[8*dbg_i+i].n1,maps.neigh[8*dbg_i+i].n2);
	}
      }
      dbg_block(1226,bound_blocks,tot_blocks);

      if (1){
	int dbg_i=149;
	for (int i=0;i<8;i++){
	  printf("%d %d: lev %d, %d %d\n",dbg_i,i,maps.neigh[8*dbg_i+i].lev,maps.neigh[8*dbg_i+i].n1,maps.neigh[8*dbg_i+i].n2);
	}
      }
      */


      maps.bound_blocks=bound_blocks;

}

int main(int ac, char ** av){

  // valori di default
  g.global_selecteddev=0; // prima scheda video trovata
  g.global_multi=0;  // default singolo livello
  g.global_level_bcc=1; // a che livello le boundary condition sono rasterizzate (>=1)
  g.global_time=0;  // non caricato
  g.global_MAXDEP=0;
  g.global_vel_eps=0;
  g.global_printstat=100;
  g.global_pend_far_field=0;
  g.global_ordine=2;
  g.global_debug=0;
  g.global_por=0; // porosita' (di default disattivata)
  
  argc = ac;
  argv = av;

  printf("%s: SVN Version %d\n",argv[0],SVNVERSION);

    if (argc<2 || 
	(argv[1][0]=='-' &&
	 strcmp(argv[1],"-decode")!=0)
	)
      help();



    std::string input_base;
    std::string input;
    char dummy[256];    
    char path[256];


    if (strcmp(argv[1],"-decode")==0){
      /// conversione
      printf("Decode utility\n");
      if (argc<4)
	help();
      char inpath[256];
      strcpy(inpath,argv[2]);
      char outpath[256];
      char outpaths[256]; // per sezioni
      char outpathsx[256]; // per sezioni
      char outpathsy[256]; // per sezioni
      strcpy(outpath,argv[3]);
      printf("decoding file %s --to--> %s\n",inpath,outpath);
      char test[256];
      char temp[256];
      int all=0;
      int res=0;
      F minx=0,maxx=0,miny=0,maxy=0;
      int okrange=0;
      int from=-1;
      int to=-1;
      int ndigits=5;
      int sezioni=0;
      for (int i=3;i<ac;i++){
	if (strcmp(argv[i],"-all")==0){
	  all=1;
	  okrange=1;
	  printf("Converting the whole map\n");
	}

	strcpy(test,"-res=");
	strcpy(temp,argv[i]);
	temp[strlen(test)]=0;
	if (strcmp(test,temp)==0){
	  res=atoi(argv[i]+strlen(test));
	  printf("Sample at resolution: %d\n",res);
	}

	if (strcmp(argv[i],"-frames")==0){
	  if (argc>i+2){
	    from=atoi(argv[i+1]);
	    to=atoi(argv[i+2]);
	    printf("Convert files from  %04d to %04d\n",from,to);
	  }
	}

	if (strcmp(argv[i],"-ndigits")==0){
	  if (argc>i+1){
	    ndigits=atoi(argv[i+1]);
	    printf("Using %d digits\n",ndigits);
	  }
	}

	if (strcmp(argv[i],"-range")==0){
	  if (argc>i+4){
	  okrange=1;
#ifndef doubleprecision
     sscanf(argv[i+1],"%f",&minx);
	 sscanf(argv[i+2],"%f",&maxx);
	 sscanf(argv[i+3],"%f",&miny);
	 sscanf(argv[i+4],"%f",&maxy);
#else
     sscanf(argv[i+1],"%lf",&minx);
	 sscanf(argv[i+2],"%lf",&maxx);
	 sscanf(argv[i+3],"%lf",&miny);
	 sscanf(argv[i+4],"%lf",&maxy);
#endif
	  printf("Using bounding box: %f %f %f %f\n",minx,maxx,miny,maxy);
	  }
	}

	if (strcmp(argv[i],"-sez")==0){
	  if (argc>i+5){
	  okrange=1;
#ifndef doubleprecision
     sscanf(argv[i+1],"%f",&minx);
	 sscanf(argv[i+2],"%f",&miny);
	 sscanf(argv[i+3],"%f",&maxx);
	 sscanf(argv[i+4],"%f",&maxy);
#else
     sscanf(argv[i+1],"%lf",&minx);
	 sscanf(argv[i+2],"%lf",&miny);
	 sscanf(argv[i+3],"%lf",&maxx);
	 sscanf(argv[i+4],"%lf",&maxy);
#endif
	  printf("Using section: %f %f %f %f\n",minx,maxx,miny,maxy);
	  sezioni=1;
	  }
	}

	if (strcmp(argv[i],"-point")==0){
	  if (argc>i+2){
	  okrange=1;
#ifndef doubleprecision
	  sscanf(argv[i+1],"%f",&minx);
	  sscanf(argv[i+2],"%f",&miny);
#else
	  sscanf(argv[i+1],"%lf",&minx);
	  sscanf(argv[i+2],"%lf",&miny);
#endif
	  maxx=minx;
	  maxy=miny;
	  res=1;
	  printf("Using point: %f %f\n",minx,miny);
	  }
	}
      }
      if (okrange==0){
	printf("Please specify the range\n");
	exit(1);
      }
	
      if (from!=-1){ // multipli
	// controlla che input abbia XXXX.formato
	int pos=strlen(inpath)-8;
	if ( inpath[pos]<'0' || inpath[pos] >'9' ||
	     inpath[pos+1]<'0' || inpath[pos+1] >'9' ||
	     inpath[pos+2]<'0' || inpath[pos+2] >'9' ||
	     inpath[pos+3]<'0' || inpath[pos+3] >'9' ||
	     inpath[pos+4]!='.'){
	  printf("error: file non valid for range of frames\n");
	  exit(1);
	}

	if (sezioni==0 && (all==1 || minx!=maxx || miny!=maxy))
	  for (int i=from;i<=to;i++) {
	    char num[5];
	    sprintf(num,"%04d",i);
	    inpath[pos]=num[0];
	    inpath[pos+1]=num[1];
	    inpath[pos+2]=num[2];
	    inpath[pos+3]=num[3];
	    printf("convert %s\n",inpath);
	    char outpath1[256];
	    strcpy(outpath1,outpath);
	    strcat(outpath1,"-");
	    strcat(outpath1,num);
	    strcat(outpath1,inpath+pos+4);	  
	    printf("--> %s\n",outpath1);
	    decode(all, res, minx,  maxx, miny, maxy, inpath, outpath1,"w+",ndigits);      
	  } else{ if (sezioni==0){/// caso punto -> 
	    printf("sample point\n");
	    FILE* fo=fopen(outpath,"w+"); // reset file
	    fclose(fo);
	    for (int i=from;i<=to;i++) {
	      char num[5];
	      sprintf(num,"%04d",i);
	      inpath[pos]=num[0];
	      inpath[pos+1]=num[1];
	      inpath[pos+2]=num[2];
	      inpath[pos+3]=num[3];
	      printf("convert %s\n",inpath);
	      decode(all, res, minx,  maxx, miny, maxy, inpath, outpath,"a+",ndigits);
	    }	  
	  }
	  else{
	    //sezioni su range
	    strcpy(outpaths,outpath);
	    strcat(outpaths,".SEZ");	  
	    FILE* fo=fopen(outpaths,"w+");
	    fclose(fo);
	    strcpy(outpathsx,outpath);
	    strcat(outpathsx,".SEZX");	  
	    fo=fopen(outpathsx,"w+");
	    fclose(fo);
	    strcpy(outpathsy,outpath);
	    strcat(outpathsy,".SEZY");	  
	    fo=fopen(outpathsy,"w+");
	    fclose(fo);
	    char inpathc[256];
	    strcpy(inpathc,inpath);
	    int l=strlen(inpath);
	    if (inpath[l-4]!='.'){
	      printf("expected *.VVX as input path (read %s)\n",inpath);
	      exit(2);
	    }
	    for (int i=from;i<=to;i++) {
	      char num[5];
	      sprintf(num,"%04d",i);
	      inpathc[pos+0]=num[0];
	      inpathc[pos+1]=num[1];
	      inpathc[pos+2]=num[2];
	      inpathc[pos+3]=num[3];

	    inpathc[l-3]='D';
	    inpathc[l-2]='E';
	    inpathc[l-1]='P';
	      decode(all, 1, minx, maxx, miny, maxy, inpathc, outpath,"w+",ndigits,3,"");

	    inpathc[l-3]='V';
	    inpathc[l-2]='V';
	    inpathc[l-1]='X';

	      decode(all, 1, minx, maxx, miny, maxy, inpathc, outpath,"w+",ndigits,1,outpathsx);
	      inpathc[l-1]='Y';	    
	      decode(all, 1, minx, maxx, miny, maxy, inpathc, outpath,"w+",ndigits,2,outpathsy);
	    }
	    std::ifstream ifsx(outpathsx);
	    std::ifstream ifsy(outpathsy);
	    fo=fopen(outpaths,"w+");
	    F x,y,dummy;
	    while (!ifsy.eof()){
	      ifsy >> x >> dummy;
	      if (!ifsy.eof()){
		ifsx >> dummy >> y;
		fprintf(fo,"%f %f\n",x,y);
	      }
	    }
	    fclose(fo);	    	    
	  }
	}

      }
      else{
	if (sezioni==0)
	  decode(all, res, minx,  maxx, miny, maxy, inpath, outpath,"w+",ndigits);      
	else{
	  strcpy(outpaths,outpath);
	  strcat(outpaths,".SEZ");	  
	  FILE* fo=fopen(outpaths,"w+");
	  fclose(fo);
	  strcpy(outpathsx,outpath);
	  strcat(outpathsx,".SEZX");	  
	  fo=fopen(outpathsx,"w+");
	  fclose(fo);
	  char inpathc[256];
	  strcpy(inpathc,inpath);
	  int l=strlen(inpath);
	  if (inpath[l-4]!='.'){
	    printf("expected *.VVX as input path (read %s)\n",inpath);
	    exit(2);
	  }
	  inpathc[l-3]='D';
	  inpathc[l-2]='E';
	  inpathc[l-1]='P';
	  decode(all, 1, minx,  maxx, miny, maxy, inpathc, outpath,"w+",ndigits,3,"");
	  inpathc[l-3]='V';
	  inpathc[l-2]='V';
	  inpathc[l-1]='X';
	  decode(all, 1, minx,  maxx, miny, maxy, inpathc, outpath,"w+",ndigits,1,outpathsx);
	  strcpy(outpathsy,outpath);
	  strcat(outpathsy,".SEZY");	  
	  fo=fopen(outpathsy,"w+");
	  fclose(fo);
	  inpathc[l-1]='Y';	    
	  decode(all, 1, minx, maxx, miny, maxy, inpathc, outpath,"w+",ndigits,2,outpathsy);
	  // manca merge
	}
      }
	exit(0);
    }


    strcpy(base_path1,argv[1]);
    int found=0;
    for (int i=strlen(base_path1)-1;i>=0;i--)
      if (base_path1[i]=='/' ||
	  base_path1[i]=='\\'){
	base_path1[i+1]=0;
	found=1;
	i=0; // break
      }
    if (found==0)
      base_path1[0]=0; // kill string

    printf("Base path: %s\n",base_path1);
    strcpy(g.base_path,base_path1);
    std::ifstream ifs(argv[1]);             /// open input file
    ifs >> input_base;                      /// read name for map files
    strcat(g.base_path,input_base.c_str());   /// full path
    printf("Base path for files: %s.*\n",g.base_path);
    ifs.getline(dummy, 256);//printf("comment: %s\n",dummy);
    

    while (!ifs.eof()){
    ifs >> input_base;
    if (!ifs.eof()){
    printf("KEY %s: ",input_base.c_str());
    if (strcmp("START",input_base.c_str())==0){
      ifs >> g.global_time_start;
      printf("%f\n",g.global_time_start);
    }
    if (strcmp("END",input_base.c_str())==0){
      ifs >> g.global_time_end;
      printf("%f\n",g.global_time_end);
    }
    if (strcmp("CR",input_base.c_str())==0){
      ifs >> g.global_CR;
      printf("%f\n",g.global_CR);
    }
    if (strcmp("LIMITER",input_base.c_str())==0){
      ifs >> g.global_limiter;
      printf("%d\n",g.global_limiter);
    }
    if (strcmp("YEPS",input_base.c_str())==0){
      ifs >> g.global_YEPS;
      printf("%f\n",g.global_YEPS);
    }
    if (strcmp("MUSCL",input_base.c_str())==0){
      ifs >> g.global_MUSCL;
      printf("%d\n",g.global_MUSCL);
    }
    if (strcmp("IBINARY",input_base.c_str())==0){
      ifs >> g.global_ibinary;
      printf("%d\n",g.global_ibinary);
    }
    if (strcmp("AL",input_base.c_str())==0){
      ifs >> g.global_al;
      printf("%d\n",g.global_al);
    }
    if (strcmp("SR",input_base.c_str())==0){
      ifs >> g.global_sr;
      printf("%d\n",g.global_sr);
    }
    if (strcmp("DTSOGLIA",input_base.c_str())==0){
      ifs >> g.global_dtsoglia;
      printf("%f\n",g.global_dtsoglia);
    }
    if (strcmp("EXPON",input_base.c_str())==0){
      ifs >> g.global_expon;
      printf("%f\n",g.global_expon);
    }
    if (strcmp("AW",input_base.c_str())==0){
      ifs >> g.global_AWSDGM;
      printf("%f\n",g.global_AWSDGM);
    }
    if (strcmp("BW",input_base.c_str())==0){
      ifs >> g.global_BWSDGM;
      printf("%f\n",g.global_BWSDGM);
    }
    if (strcmp("METODO",input_base.c_str())==0){
      ifs >> g.global_metodo;
      printf("%d\n",g.global_metodo);
    }
    if (strcmp("NITER",input_base.c_str())==0){
      ifs >> g.global_niter;
      printf("%d\n",g.global_niter);
    }
    if (strcmp("DTOUTPUT",input_base.c_str())==0){
      ifs >> g.global_dtoutput;
      printf("%f\n",g.global_dtoutput);
    }
    if (strcmp("MAXDEP",input_base.c_str())==0){
      ifs >> g.global_MAXDEP;
      printf("%f\n",g.global_MAXDEP);
    }
    if (strcmp("VELEPS",input_base.c_str())==0){
      ifs >> g.global_vel_eps;
      printf("%f\n",g.global_vel_eps);
    }
    if (strcmp("PRINTSTAT",input_base.c_str())==0){
      ifs >> g.global_printstat;
      printf("%d\n",g.global_printstat);
    }
    if (strcmp("PENDFARFIELD",input_base.c_str())==0){
      ifs >> g.global_pend_far_field;
      printf("%f\n",g.global_pend_far_field);
    }




    ifs.getline(dummy, 256);//printf("comment: %s\n",dummy);
    }
    }


    ifs.close();


    //file BRE
    strcpy(path,g.base_path);strcat(path,".BRE");   
    std::ifstream ifsb(path);             /// open input file
    if (!ifsb){
      g.br_time_start=10e10; // non comincia mai
      g.br_time_end=10e10; // non comincia mai
    }
    else
    while (!ifsb.eof()){
      ifsb >> input_base;
      if (!ifsb.eof()){
	printf("KEY %s: ",input_base.c_str());
	if (strcmp("TIMESTART",input_base.c_str())==0){
	  ifsb >> g.br_time_start;
	  printf("%f\n",g.br_time_start);
	}
	if (strcmp("TIMEEND",input_base.c_str())==0){
	  ifsb >> g.br_time_end;
	  printf("%f\n",g.br_time_end);
	}
	if (strcmp("TIMEENDL",input_base.c_str())==0){
	  ifsb >> g.br_time_end_l;
	  printf("%f\n",g.br_time_end_l);
	}
	if (strcmp("BSTART",input_base.c_str())==0){
	  ifsb >> g.B_iniz;
	  printf("%f\n",g.B_iniz);
	}
	if (strcmp("BEND",input_base.c_str())==0){
	  ifsb >> g.B_fin;
	  printf("%f\n",g.B_fin);
	}
	if (strcmp("BETA",input_base.c_str())==0){
	  ifsb >> g.beta;
	  printf("%f\n",g.beta);
	}
	if (strcmp("ZSTART",input_base.c_str())==0){
	  ifsb >> g.z_iniz;
	  g.z_iniz = g.z_iniz-g.offset_btm;				//con offset_btm abbasso quota terreno zstart e zend 
	  printf("%f\n",g.z_iniz);
	}
	if (strcmp("ZEND",input_base.c_str())==0){
	  ifsb >> g.z_fin;
	  g.z_fin = g.z_fin-g.offset_btm;
	  printf("%f\n",g.z_fin);
	}
	if (strcmp("MAXWIDTH",input_base.c_str())==0){
	  ifsb >> g.width_max;
	  printf("%f\n",g.width_max);
	}
	if (strcmp("XBR",input_base.c_str())==0){
	  ifsb >> g.x_br;
	  printf("%f\n",g.x_br);
	}
	if (strcmp("YBR",input_base.c_str())==0){
	  ifsb >> g.y_br;
	  printf("%f\n",g.y_br);
	}
	if (strcmp("DIRXBR",input_base.c_str())==0){
	  ifsb >> g.dir_br_x;
	  printf("%f\n",g.dir_br_x);
	}
	if (strcmp("DIRYBR",input_base.c_str())==0){
	  ifsb >> g.dir_br_y;
	  printf("%f\n",g.dir_br_y);
	}
	ifsb.getline(dummy, 256);//printf("comment: %s\n",dummy);
      }
    }
    ifsb.close();

    /*
    float br_time_start=3600; //tempo di inizio evoluzione della breccia
    float br_time_end=3600*3; //tempo di fine della evoluzione della breccia

    float B_iniz=5; // larghezza breccia al tempo time_iniz
    float B_fin=50; // larghezza della base della breccia al tempo time_fin
    float beta=0.1745; //angolo di inclinazione delle sponde della breccia in radianti
    float z_iniz=29.5; //quota minima della breccia al tempo time_iniz
    float z_fin=22; //quota minima della breccia al tempo time_br

    float width_max=12; //massima larghezza della breccia (direzione normale a quella di evoluzione)
    //coordinate del centro della breccia
    float x_br=500;
    float y_br=500;

//versore che individua la direzione di propagazione della breccia (45 gradi)
    float dir_br_x=0.7071;
    float dir_br_y=0.7071;
*/


    char inhpath[256];
    char vhxpath[256];
    char vhypath[256];
    char btmpath[256];
    char inhpath_[256]; // 
    char vhxpath_[256];
    char vhypath_[256];
    char btmpath_[256];
    int correctvx=0; // se bisogna moltiplicare per h -->1
    strcpy(inhpath,g.base_path);strcat(inhpath,".INH");
    strcpy(vhxpath,g.base_path);strcat(vhxpath,".VHX");
    strcpy(vhypath,g.base_path);strcat(vhypath,".VHY");
    strcpy(btmpath,g.base_path);strcat(btmpath,".BTM");
    strcpy(inhpath_,"");
    strcpy(vhxpath_,"");
    strcpy(vhypath_,"");
    strcpy(btmpath_,"");

    /// lettura commandline
    simulate=1; // default: no gpu
    /// process input parameters
    char temp[256];
    char test[256];
    for (int i=1;i<argc;i++){

      //      printf("%d: %s\n",i,argv[i]);

      strcpy(test,"-h");
      strcpy(temp,argv[i]);
      temp[strlen(test)]=0;
      if (strcmp(test,temp)==0){
	help();
      }
      strcpy(test,"--help");
      strcpy(temp,argv[i]);
      temp[strlen(test)]=0;
      if (strcmp(test,temp)==0){
	help();
      }

      strcpy(test,"-debug=");
      strcpy(temp,argv[i]);
      temp[strlen(test)]=0;
      if (strcmp(test,temp)==0){
	g.global_debug=atoi(argv[i]+strlen(test));
	printf("Debug: %d\n",g.global_debug);
      }

      strcpy(test,"-por=");
      strcpy(temp,argv[i]);
      temp[strlen(test)]=0;
      if (strcmp(test,temp)==0){
	g.global_por=atoi(argv[i]+strlen(test));
	printf("Porosity: %d\n",g.global_por);
      }

      strcpy(test,"-multi=");
      strcpy(temp,argv[i]);
      temp[strlen(test)]=0;
      if (strcmp(test,temp)==0){
	if (argv[i][7]=='h' &&
	    argv[i][8]=='i') {
	  g.global_multi=1; // codice speciale per tutto a massima res
	  hi=1;
	}
	else
	g.global_multi=atoi(argv[i]+strlen(test));
	printf("Multiresolution: %d (10=hi)\n",g.global_multi);
      }

      strcpy(test,"-bcc-level=");
      strcpy(temp,argv[i]);
      temp[strlen(test)]=0;
      if (strcmp(test,temp)==0){
	if (g.global_multi!=0){
	  g.global_level_bcc=atoi(argv[i]+strlen(test));
	  if (g.global_level_bcc<1){
	    printf("Error: multiresolution range 1..4\n");
	    exit(2);
	  }
	    
	  printf("Boundary condition rasterization at level: %d\n",g.global_level_bcc);
	}
	else
	  printf("Warning: Ignoring bcc-level (no multiresolution enabled)\n");
      }



      strcpy(test,"-order=");
      strcpy(temp,argv[i]);
      temp[strlen(test)]=0;
      if (strcmp(test,temp)==0){
	g.global_ordine=atoi(argv[i]+strlen(test));
	printf("Order: %d\n",g.global_ordine);
      }

      strcpy(test,"-gpu=");
      strcpy(temp,argv[i]);
      temp[strlen(test)]=0;
      if (strcmp(test,temp)==0){
	g.global_selecteddev=atoi(argv[i]+strlen(test));
	printf("GPU selected: %d\n",g.global_selecteddev);
      }

      strcpy(test,"-opengl=");
      strcpy(temp,argv[i]);
      temp[strlen(test)]=0;
      if (strcmp(test,temp)==0){
	simulate = 0;
	stepperframe = atoi(argv[i]+strlen(test));
#if defined(NOOPENGL) || defined(WIN)
	simulate=1;
	printf("Sorry: no opengl support. Recompile by enabling the option (not supported in Windows)\n");
	exit(1);
#endif
	printf("Step: %d\n",stepperframe);
      }

      strcpy(test,"-resume=");
      strcpy(temp,argv[i]);
      temp[strlen(test)]=0;
      if (strcmp(test,temp)==0){
	correctvx=1;
	char name[256];
	//if (g.global_multi==0)
	  {
	  strcpy(name,argv[i]+strlen(test));

	  strcpy(inhpath_,g.base_path);
	  strcat(inhpath_,"-output-");
	  strcat(inhpath_,name);
	  strcat(inhpath_,".WSE");

	  strcpy(vhxpath_,g.base_path);
	  strcat(vhxpath_,"-output-");
	  strcat(vhxpath_,name);
	  strcat(vhxpath_,".VVX");

	  strcpy(vhypath_,g.base_path);
	  strcat(vhypath_,"-output-");
	  strcat(vhypath_,name);
	  strcat(vhypath_,".VVY");	
        
	  g.global_time_start=g.global_dtoutput*(1+atoi(name)); // corregge tempo inizio

	  printf("time start =%f, breach time start=%f \n", g.global_time_start, g.br_time_start);
 
	  if (g.global_time_start >=g.br_time_start/3600 ) {
	    strcpy(btmpath_,g.base_path);
	    strcat(btmpath_,"-output-");
	    strcat(btmpath_,name);
	    strcat(btmpath_,".BTM");

	    printf("btmpath %s \n", btmpath_);
	  }	
	}
      }
    }

#ifdef doubleprecision
    printf("Using DOUBLE precision\n");
#else
    printf("Using SINGLE precision\n");
#endif
    // lettura mappe
    if (inhpath_[0]!=0 && g.global_multi==0) // se c'e' resume senza multires -> carico nuova mappa, altrimenti con resume multires, intanto carico dati originali, e poi sovrascrivo valori corretti
      strcpy(inhpath,inhpath_);
    if (vhxpath_[0]!=0 && g.global_multi==0) // se c'e' resume senza multires -> carico nuova mappa, altrimenti con resume multires, intanto carico dati originali, e poi sovrascrivo valori corretti
      strcpy(vhxpath,vhxpath_);
    if (vhypath_[0]!=0 && g.global_multi==0) // se c'e' resume senza multires -> carico nuova mappa, altrimenti con resume multires, intanto carico dati originali, e poi sovrascrivo valori corretti
      strcpy(vhypath,vhypath_);
    if (btmpath_[0]!=0 && g.global_multi==0) // se c'e' resume senza multires -> carico nuova mappa, altrimenti con resume multires, intanto carico dati originali, e poi sovrascrivo valori corretti
      strcpy(btmpath,btmpath_);


    printf("read inh: ");
    readDSAA(inhpath,0);//nrows, ncols, inh_map, minx_map, maxx_map, miny_map, maxy_map,minv_map, maxv_map);

    maps.dx=(maps.maxx_map-maps.minx_map)/(maps.ncols-1);
    maps.dy=(maps.maxy_map-maps.miny_map)/(maps.nrows-1);
    printf("delta %f %f\n",maps.dx,maps.dy);
    //    printf("delta %f %f %d \n",maps.maxx_map,maps.minx_map,maps.ncols);

    strcpy(path,g.base_path);
    readSEZ(path);//nrows, ncols, minx_map, maxx_map, miny_map, maxy_map); 

    printf("read vhx: ");
    readDSAA(vhxpath,1);
    printf("read vhy: ");
    readDSAA(vhypath,2);
        
    if (g.global_por){
      printf("read por: ");
      strcpy(path,g.base_path);strcat(path,".POR");
      readDSAA(path, 3);

      for(int y = 1; y < maps.nrows-1; y ++)
	for(int x = 1; x < maps.ncols-1; x ++)
	  if(maps.por_map[y * maps.ncols + x]==0){
	    printf("ERROR: zero porosity non allowed, please verify cell %d %d\n",x-1,y-1);
	    //exit(2);
	  }
      printf("read hlx: ");
      strcpy(path,g.base_path);strcat(path,".HLX");
      readDSAA(path,  4);
      printf("read hly: ");
      strcpy(path,g.base_path);strcat(path,".HLY");
      readDSAA(path, 5);
    }

    printf("read btm: ");
    //strcpy(path,g.base_path);strcat(path,".BTM");
    g.offset_btm=10e10;   
	readDSAA(btmpath, 6);      //AF leggo BTM e calcolo offset_btm
	printf("Valor minimo batimetria (offset_btm): %3.2f \n", g.offset_btm);

	//------abbasso BTM e INH dell'offset_btm----------
	for(int y = 1; y < maps.nrows-1; y++){
       for(int x = 1; x < maps.ncols-1; x++){
				maps.inh_map[x + y * maps.ncols] =  maps.inh_map[x + y * maps.ncols] - g.offset_btm;
				maps.btm_map[x + y * maps.ncols] = maps.btm_map[x + y * maps.ncols] - g.offset_btm;
		   }
	}	 
	 //----------------
	 
    printf("read man: ");
    strcpy(path,g.base_path);strcat(path,".MAN");
    readDSAA(path, 7);
    strcpy(path,g.base_path);




    //necessario anche per multirisoluzione (uso subito per sistemare inh)
    readBLNBCC(path,maps.btm_map); //raster livello 0



    /*
// reset output file list
    char path1[256];
    strcpy(path1,g.base_path);
    strcat(path1,"-output.txt");
    printf("Erase %s\n",path1);
    FILE* f1=fopen(path1,"w");
    fclose(f1);
    printf("erase ok\n");
    */
    
    //legge punti per multirisoluzione
    //impone richieste specifiche
    strcpy(path,g.base_path);strcat(path,".PTS");
    printf("read pts %s\n",path);


    std::ifstream ifpts(path);             /// open input file
    if (!ifpts.fail()){
      int ctp=0;
      while (!ifpts.eof()){
      F4 temp;
      ifpts >> temp.x; //x
      if (!ifpts.eof()){
      ifpts >> temp.y; //y
      ifpts >> temp.z; // livello (1..4) int
      ifpts >> temp.w; // bool tappo
      printf("Read multires point: %f %f %d %d\n",temp.x,temp.y,(int)temp.z,(int)temp.w);
      ctp++;
      g.punti_m.push_back(temp);
      }

    }
      printf("read %d multires points\n",ctp);
   }
      ////////////////////////////////////////////////

    for(int y = 0; y < maps.nrows; y ++){
      for(int x = 0; x < maps.ncols; x ++)
	if (maps.man_map[y * maps.ncols + x]>10e10){
	  maps.man_map[y * maps.ncols + x]=0;
	}
    }


    if (0)
    for(int y = 0; y < maps.nrows; y ++){
      for(int x = 0; x < maps.ncols; x ++){
	printf("%3.2f ",maps.btm_map[y * maps.ncols + x]);
      }
      printf("\n");
    }
    if(0)
    for(int y = 0; y < maps.nrows; y ++){ 
      for(int x = 0; x < maps.ncols; x ++){
	maps.inh_map[y * maps.ncols + x]=40+10*(F)x/maps.ncols;
	//maps.btm_map[y * maps.ncols + x]=10*(F)x/maps.ncols;
      }
    }

    if(0)
    for(int y = 0; y < maps.nrows; y ++){ 
      for(int x = 0; x < maps.ncols; x ++)
	if (maps.btm_map[y * maps.ncols + x]<1000){
	  maps.btm_map[y * maps.ncols + x]=20;
	//maps.btm_map[y * maps.ncols + x]=10*(F)x/maps.ncols;
      }
    }
    
    printf("init gpu\n");
    initGPU(); // block setup


    int x_blocks=-1;
    int y_blocks=-1;
    if (g.global_multi){
      multires(g.base_path);
      x_blocks=1<<(int)(ceil(log((float)maps.tot_blocks)/log((float)2)/(float)2)); 
      y_blocks=(maps.tot_blocks-1)/x_blocks+1;

      int sizet=4*x_blocks*y_blocks*BLOCKSIZE_X*BLOCKSIZE_Y;
      printf("alloc size %d\n",sizet*sizeof(F4));
      if (sizet<maps.nrows * maps.ncols)
	sizet=maps.nrows * maps.ncols;
      maps.debug = (F4 *) malloc(sizet * sizeof(F4)); // usato per varie cose... per ora alloco largo
      maps.temp = (F *) malloc(sizet * sizeof(F));

      maps.watersurfacevertices = (F4 *) malloc(sizet * sizeof(F4)); // usato per varie cose... per ora alloco largo
      maps.colors = (rgb *) malloc(sizet * sizeof(rgb));

      maps.landscape = (F4 *) malloc(sizet * sizeof(F4));  // per visualizzazione
      maps.landscapecolors = (rgb *) malloc(sizet * sizeof(rgb));  // per visualizzazione
      
      if (inhpath_[0]!=0){ // resume multires, sovrascrivo valori corretti
	readDSAA(inhpath_,8);
	for(int y = 0; y < y_blocks*BLOCKSIZE_Y; y ++)  // sembra che usino coordinata x come cartesiana e -y (non le coordiate matrice riga, colonna)
	  for(int x = 0; x < x_blocks*BLOCKSIZE_X; x ++){
	    maps.host_grid_multi[ y*x_blocks*BLOCKSIZE_X+x].x = maps.temp[ y*x_blocks*BLOCKSIZE_X+x];
	  }
      }
      if (vhxpath_[0]!=0){ // resume multires, sovrascrivo valori corretti
	readDSAA(vhxpath_,8);
	for(int y = 0; y < y_blocks*BLOCKSIZE_Y; y ++)  // sembra che usino coordinata x come cartesiana e -y (non le coordiate matrice riga, colonna)
	  for(int x = 0; x < x_blocks*BLOCKSIZE_X; x ++){
	    maps.host_grid_multi[ y*x_blocks*BLOCKSIZE_X+x].y = maps.temp[ y*x_blocks*BLOCKSIZE_X+x];
	  }
      }
      if (vhypath_[0]!=0){ // resume multires, sovrascrivo valori corretti
	readDSAA(vhypath_,8);
	for(int y = 0; y < y_blocks*BLOCKSIZE_Y; y ++)  // sembra che usino coordinata x come cartesiana e -y (non le coordiate matrice riga, colonna)
	  for(int x = 0; x < x_blocks*BLOCKSIZE_X; x ++){
	    maps.host_grid_multi[ y*x_blocks*BLOCKSIZE_X+x].z = maps.temp[ y*x_blocks*BLOCKSIZE_X+x];
	  }
      }

      if (btmpath_[0]!=0){ // resume multires, sovrascrivo valori corretti
	readDSAA(btmpath_,8);
	for(int y = 0; y < y_blocks*BLOCKSIZE_Y; y ++)  // sembra che usino coordinata x come cartesiana e -y (non le coordiate matrice riga, colonna)
	  for(int x = 0; x < x_blocks*BLOCKSIZE_X; x ++){
	    maps.host_grid_multi[ y*x_blocks*BLOCKSIZE_X+x].w = maps.temp[ y*x_blocks*BLOCKSIZE_X+x];
	  }
      }

      if (0)
      for (int i=0;i<maps.tot_blocks;i++)
	printf("%d: %d %d\n",i,maps.host_ofs_blocks[i].x,maps.host_ofs_blocks[i].y);

      for(int y = 0; y < y_blocks*BLOCKSIZE_Y; y ++)  // sembra che usino coordinata x come cartesiana e -y (non le coordiate matrice riga, colonna)
	for(int x = 0; x < x_blocks*BLOCKSIZE_X; x ++){
	  F4 v;
	  int width=x_blocks*BLOCKSIZE_X;
	  F dx1=maps.dx;
	  F dy1=maps.dy;
	  F ofsx=maps.host_ofs_blocks[y/BLOCKSIZE_Y*x_blocks+x/BLOCKSIZE_X].x;
	  F ofsy=maps.host_ofs_blocks[y/BLOCKSIZE_Y*x_blocks+x/BLOCKSIZE_X].y;
	  F dx,dy;
	  int x1=x%BLOCKSIZE_X;
	  int y1=y%BLOCKSIZE_X;
	  //correzione dx per multirisoluzione
	  dx=dx1*(1<<maps.host_grid_level_multi[y/BLOCKSIZE_Y*x_blocks+x/BLOCKSIZE_X]);
	  dy=dy1*(1<<maps.host_grid_level_multi[y/BLOCKSIZE_Y*x_blocks+x/BLOCKSIZE_X]);
	  //	  printf("%d %d -> idx %d: %f %f\n",x,y,y/BLOCKSIZE_Y*x_blocks+x/BLOCKSIZE_X,ofsx,ofsy);
	  v.y=maps.host_grid_multi[ y*x_blocks*BLOCKSIZE_X+x].w;
	  v.x = (dx1*ofsx+dx *x1 / F(BLOCKSIZE_X)) * 10 - 5;
	  v.z = (dx1*ofsy+dy *y1 / F(BLOCKSIZE_Y)) * 10 - 5;
	  maps.landscape[4*(y * width + x)] = v;
	  v.x = (dx1*ofsx+dx *(1+x1) / F(BLOCKSIZE_X)) * 10 - 5;
	  v.z = (dx1*ofsy+dy *y1 / F(BLOCKSIZE_Y)) * 10 - 5;
	  maps.landscape[4*(y * width + x)+1] = v;
	  v.x = (dx1*ofsx+dx *(1+x1) / F(BLOCKSIZE_X)) * 10 - 5;
	  v.z = (dx1*ofsy+dy *(1+y1) / F(BLOCKSIZE_Y)) * 10 - 5;
	  maps.landscape[4*(y * width + x)+2] = v;
	  v.x = (dx1*ofsx+dx *x1 / F(BLOCKSIZE_X)) * 10 - 5;
	  v.z = (dx1*ofsy+dy *(1+y1) / F(BLOCKSIZE_Y)) * 10 - 5;
	  maps.landscape[4*(y * width + x)+3] = v;
	}

    }// end multires preproc
    else{
      int nrows=maps.nrows;
      int ncols=maps.ncols;
      maps.watersurfacevertices = (F4 *) malloc(nrows * ncols * sizeof(F4));
      maps.debug = (F4 *) malloc(nrows * ncols * sizeof(F4));
      maps.temp = (F *) malloc(nrows * ncols * sizeof(F4));
      maps.colors = (rgb *) malloc(nrows * ncols * sizeof(rgb));
      maps.landscape = (F4 *) malloc(nrows * ncols * sizeof(F4));
      maps.landscapecolors = (rgb *) malloc(nrows * ncols * sizeof(rgb));  // per visualizzazione

      for(int y = 0; y < nrows; y ++)  // sembra che usino coordinata x come cartesiana e -y (non le coordiate matrice riga, colonna)
	for(int x = 0; x < ncols; x ++){
	  F deltax=(maps.maxx_map-maps.minx_map);
	  F deltay=(maps.maxy_map-maps.miny_map);
	  F ratio=deltax/deltay;
	  maps.landscape[y * ncols + x].x = x / F(ncols - 1) * 10 - 5;
	  maps.landscape[y * ncols + x].z = (1-y / F(nrows - 1)) * 10/ratio - 5/ratio;
	  if (maps.btm_map[y * ncols + x]>10e9 || maps.host_info[y * ncols + x].x==BIT_EXTERN)
	    maps.landscape[y * ncols + x].y=0;	
	  else
	    maps.landscape[y * ncols + x].y=(maps.btm_map[y * ncols + x]-maps.min_btm)/maps.scale;	
	}

      for(int y = 0; y < nrows; y ++)  // sembra che usino coordinata x come cartesiana e -y (non le coordiate matrice riga, colonna)
	for(int x = 0; x < ncols; x ++){
	  rgb c;
	  if (maps.btm_map[y * ncols + x]>10E9 || maps.host_info[y * ncols + x].x==BIT_EXTERN){
	    c.y=0;
	    c.x = 0;
	  }
	  else{
	    c.x = 128;
	    F temp=(maps.btm_map[y * ncols + x]-maps.min_btm)/(maps.max_btm-maps.min_btm);
	    if (temp<0) temp=0;
	    if (temp>1) temp=1;
	    c.y = (unsigned char)(255*temp);
	  }
	  c.z = 0;
	  maps.landscapecolors[y * ncols + x] = c;
	}
    }

    /// fix inh_map e min max btm
    maps.min_btm=10000;
    maps.max_btm=-10000;
    maps.max_inh=-10000;
    maps.min_inh=10000;
    maps.max_dep=0;
    maps.max_ini_uh=0;
    maps.max_uh=0;

    int ncols;
    int nrows;

    if (g.global_multi==0){
    ncols = maps.ncols;
    nrows= maps.nrows;
    }
    else{
      ncols=x_blocks*BLOCKSIZE_X;
      nrows=y_blocks*BLOCKSIZE_Y;
    }

    if (g.global_multi==0) {
    for(int y = 0; y < nrows; y ++){ 
      for(int x = 0; x < ncols; x ++)
	if(maps.host_info[y * ncols + x].x!=BIT_EXTERN){ // non e' cella exclusa
	  if (maps.inh_map[y * ncols + x]<maps.btm_map[y * ncols + x] || maps.inh_map[y * ncols + x] > 1.E10){
	    //	    printf("warning h negativo: riga %d col %d: %f %f\n",y,x,inh_map[y * ncols + x],btm_map[y * ncols + x]);
	    maps.inh_map[y * ncols + x]=maps.btm_map[y * ncols + x];
	  } 
	  F val=maps.btm_map[y * ncols + x];
	  if (val<10E9 && val>=0){
	    if (val>maps.max_btm)
	      maps.max_btm=val;
	    if (val<maps.min_btm)
	      maps.min_btm=val;
	  }
	  if (maps.min_inh>maps.inh_map[y * ncols + x])
	    maps.min_inh=maps.inh_map[y * ncols + x];
	  if (maps.max_inh<maps.inh_map[y * ncols + x] && maps.inh_map[y * ncols + x]<10000)
	    maps.max_inh=maps.inh_map[y * ncols + x]; 	  
	}
    }

    for(int y = 0; y < nrows; y ++){ 
      for(int x = 0; x < ncols; x ++)
	if(maps.host_info[y * ncols + x].x!=BIT_EXTERN){ // non e' cella exclusa
	  if (maps.max_dep<maps.inh_map[y * ncols + x]-maps.btm_map[y * ncols + x] && maps.btm_map[y * ncols + x]<10000)
	    maps.max_dep=maps.inh_map[y * ncols + x]-maps.btm_map[y * ncols + x];
	}
    }
    if (maps.max_dep<g.global_MAXDEP){
      maps.max_dep=g.global_MAXDEP;
    }

      for(int y = 0; y < nrows; y ++){ 
        for(int x = 0; x < ncols; x ++){
          if (maps.inh_map[y * ncols + x]-maps.btm_map[y * ncols + x]<g.global_YEPS){
	    maps.vvx_map[y * ncols + x]=0;
	    maps.vvy_map[y * ncols + x]=0;
	  }	
      } 
    }

    if (correctvx){ // moltiplica per h
      for(int y = 0; y < nrows; y ++){ 
        for(int x = 0; x < ncols; x ++){
	
          if (maps.inh_map[y * ncols + x]-maps.btm_map[y * ncols + x]>g.global_YEPS){
	    maps.vvx_map[y * ncols + x]*=( maps.inh_map[y * ncols + x]-maps.btm_map[y * ncols + x] );
	    maps.vvy_map[y * ncols + x]*=( maps.inh_map[y * ncols + x]-maps.btm_map[y * ncols + x] );
	  }else{
            maps.vvx_map[y * ncols + x]=0.;
	    maps.vvy_map[y * ncols + x]=0.;
        }	
      } 
    }
    }

    for(int y = 0; y < nrows; y ++){ 
      for(int x = 0; x < ncols; x ++)
	if(maps.host_info[y * ncols + x].x!=BIT_EXTERN && 
	   maps.inh_map[y * ncols + x]<10e10
	   ){ // non e' cella exclusa
	  F h=maps.inh_map[y * ncols + x]-maps.btm_map[y * ncols + x];
	  if (h>g.global_YEPS){
	    F v=fabs(maps.vvx_map[y * ncols + x]);
	    if (maps.max_ini_uh<h*v)
	      maps.max_ini_uh=h*v;
	    v=fabs(maps.vvy_map[y * ncols + x]);
	    if (maps.max_ini_uh<h*v)
	      maps.max_ini_uh=h*v;
	  }
	}
    }



    }
    else { // multires
    for(int y = 0; y < nrows; y ++){ 
      for(int x = 0; x < ncols; x ++)
	if(maps.host_info_multi[y * ncols + x].x!=BIT_EXTERN){ // non e' cella exclusa
	  if (maps.host_grid_multi[y * ncols + x].x<maps.host_grid_multi[y * ncols + x].w || 
	      maps.host_grid_multi[y * ncols + x].x > 1.E10){
	    maps.host_grid_multi[y * ncols + x].x=maps.host_grid_multi[y * ncols + x].w;
	  } 
	  F val=maps.host_grid_multi[y * ncols + x].w;
	  if (val<10E9 && val>=0){
	    if (val>maps.max_btm)
	      maps.max_btm=val;
	    if (val<maps.min_btm)
	      maps.min_btm=val;
	  }
	  if (maps.min_inh>maps.host_grid_multi[y * ncols + x].x)
	    maps.min_inh=maps.host_grid_multi[y * ncols + x].x;
	  if (maps.max_inh<maps.host_grid_multi[y * ncols + x].x && 
	      maps.host_grid_multi[y * ncols + x].x<10000)
	    maps.max_inh=maps.host_grid_multi[y * ncols + x].x; 	  
	}
    }

    for(int y = 0; y < nrows; y ++){ 
      for(int x = 0; x < ncols; x ++)
	if(maps.host_info_multi[y * ncols + x].x!=BIT_EXTERN){ // non e' cella exclusa
	  if (maps.max_dep<maps.host_grid_multi[y * ncols + x].x-maps.host_grid_multi[y * ncols + x].w && 
	      maps.host_grid_multi[y * ncols + x].w<10000)
	    maps.max_dep=max(maps.host_grid_multi[y * ncols + x].x, maps.host_grid_multi[y * ncols + x].w);
	}
    }
    if (maps.max_dep<g.global_MAXDEP){
      maps.max_dep=g.global_MAXDEP;
    }

      for(int y = 0; y < nrows; y ++){ 
        for(int x = 0; x < ncols; x ++){
          if (maps.host_grid_multi[y * ncols + x].x-maps.host_grid_multi[y * ncols + x].w<g.global_YEPS){
	    maps.host_grid_multi[y * ncols + x].y=0;
	    maps.host_grid_multi[y * ncols + x].z=0;
	  }	
      } 
    }

    if (correctvx){ // moltiplica per h
      for(int y = 0; y < nrows; y ++){ 
        for(int x = 0; x < ncols; x ++){
	
          if (maps.host_grid_multi[y * ncols + x].x-maps.host_grid_multi[y * ncols + x].w>g.global_YEPS){
	    maps.host_grid_multi[y * ncols + x].y*=( maps.host_grid_multi[y * ncols + x].x-maps.host_grid_multi[y * ncols + x].w );
	    maps.host_grid_multi[y * ncols + x].z*=( maps.host_grid_multi[y * ncols + x].x-maps.host_grid_multi[y * ncols + x].w );
	  }else{
            maps.host_grid_multi[y * ncols + x].y=0.;
	    maps.host_grid_multi[y * ncols + x].z=0.;
        }	
      } 
    }
    }

    for(int y = 0; y < nrows; y ++){ 
      for(int x = 0; x < ncols; x ++)
	if(maps.host_info_multi[y * ncols + x].x!=BIT_EXTERN && 
	   maps.host_grid_multi[y * ncols + x].x<10e10
	   ){ // non e' cella exclusa
	  F h=maps.host_grid_multi[y * ncols + x].x-maps.host_grid_multi[y * ncols + x].w;
	  if (h>g.global_YEPS){
	    F v=fabs(maps.host_grid_multi[y * ncols + x].y);
	    if (maps.max_ini_uh<h*v)
	      maps.max_ini_uh=h*v;
	    v=fabs(maps.host_grid_multi[y * ncols + x].z);
	    if (maps.max_ini_uh<h*v)
	      maps.max_ini_uh=h*v;
	  }
	}
    }


    }

    /////////////////






    maps.max_uh=maps.max_dep*2*pow(maps.max_dep*9.81,0.5);
    if (maps.max_uh<maps.max_ini_uh)
      maps.max_uh=maps.max_ini_uh;


    F maxx=max(maps.maxv_map,maps.max_btm);
    maxx=max(maps.max_inh,maps.max_btm);
    
    maps.scale=0.5*(maps.maxx_map-maps.minx_map)/(maxx-maps.min_btm);
    maps.scale=(maxx-maps.min_btm)/8.0;
    printf("opengl z scale %f \n",maps.scale);
    printf("gridsize %d x %d\n", ncols,nrows);



    initWaterSurface();
    printf("init ok\n");

    /*
    /// ora posso cancellare valori alti da landscape, dato che e' tutto copiato su gpu
    for(int y = 0; y < nrows; y ++)  // sembra che usino coordinata x come cartesiana e -y (non le coordiate matrice riga, colonna)
      for(int x = 0; x < ncols; x ++)
	if (maps.landscape[y * ncols + x].y>10E9)
	  maps.landscape[y * ncols + x].y=0;	
	else
	  maps.landscape[y * ncols + x].y=(maps.landscape[y * ncols + x].y-maps.min_btm)/maps.scale;	
    */
      
    /////////////
    // errori macchina
    printf("ranges:\n");
    printf("inh: %f %f\n",maps.min_inh,maps.max_inh);
    printf("btm: %f %f\n",maps.min_btm,maps.max_btm);
    printf("maxdep: %f\n",maps.max_dep);
    printf("max ini u: %f\n",maps.max_ini_uh);
    printf("max u: %f\n",maps.max_uh);   
    if (maps.max_uh<maps.max_ini_uh)
      maps.max_uh=maps.max_ini_uh;

#ifndef doubleprecision
     maps.err_h=maps.max_dep/pow(2.0f,23);
     maps.err_uh=maps.max_uh/pow(2.0f,23);
#else
     maps.err_h=maps.max_dep/pow(2.0f,52);
     maps.err_uh=maps.max_uh/pow(2.0f,52);
#endif

    printf("Float precision h: %g\n",maps.err_h);
    printf("Float precision uh: %g\n",maps.err_uh);

#if !defined(NOOPENGL) && !defined(WIN)
    if (simulate == 0)
      createWindow(argc, argv, 800, 600, stepperframe, kernelflops);
    else
#endif
    {
      initTimer();
      //cudaProfilerStart();
      kernelcalls+=computeNext(0);
      //cudaProfilerEnd();
      cudaDeviceReset();
      F runtime = timeSinceInit() / 1000.0f;
      long threads = nrows* ncols;
      long flops = kernelcalls * threads * kernelflops;

      std::cout << "Grid: width " << ncols << "x height" << nrows << std::endl;
      std::cout << "Launched the main kernel " << kernelcalls << " times." << std::endl;
      printf("Execution time: %2.4fs\n", runtime);
      if(flops > 0)
        {
	  std::cout << "Total computed flops: " << flops << std::endl;
	  printf("Speed: %6.2f GFlop/s\n", flops / runtime / 1000000000);
        }
    }
    return 0;
}

/*// solo switch tra tipi di dato
int main(int argc, char ** argv){
    /// lettura commandline
    /// process input parameters
    char temp[256];
    char test[256];
    int prec=0;
    for (int i=1;i<argc;i++){
      strcpy(test,"-precision=");
      strcpy(temp,argv[i]);
      temp[strlen(test)]=0;
      if (strcmp(test,temp)==0){
	prec=atoi(argv[i]+strlen(test));
	printf("Select precision: %d\n",prec);
      }
    }
    if (prec==0)
      return main1<float,float4>(argc,argv);
    else
      return main1<double,double4>(argc,argv);
}
*/
