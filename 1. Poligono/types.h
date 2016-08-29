#pragma once
#include "vector_types.h"
#include <string>
#include <vector>
using namespace std;

#ifndef types_h
#define types_h

//#define WINDBG 1

typedef uchar4 rgb;

//#define doubleprecision
//#if (__CUDA_ARCH__ < 200 && defined(doubleprecision))
//#error no support for type double in CUDA <2.0 
//#endif

#ifndef doubleprecision
typedef float F; 
typedef float4 F4; 
#else
typedef double F; 
typedef double4 F4; 
#endif

typedef int reflection;
// trucco per marcare extern: tutti i 4 bit a 1
#define BIT_EXTERN 15
#define BIT_TAPPO 14
#define BIT_N 1
#define BIT_S 2
#define BIT_W 4
#define BIT_E 8

typedef struct point_ {
  float x,y;
} point;

typedef struct polygon_ {
  int fill; // 1 = orario, !=1 antiorario
  vector<point> points;
  vector<int> edges;    // edges[i] refers to points[i]...points[i+1%n]
} polygon;


typedef struct cons_edge_ {
  float time,h,uh,vh;
} cons_edge;

typedef struct cons_ {
  int edge; // in riferimenti interni (parto da 0)
  vector<cons_edge> cons_edges;
} cons;

typedef struct {
  char lev;
  int n1;
  int n2;
} neigh_t;

struct maps_{
  //maps
  int ncols,nrows;  // dimensione matrice

  F scalex,scaley;
  F scale;
  F4 *landscape;
  rgb *landscapecolors;
  F4 *watersurfacevertices;
  F4 *debug;
  F *temp;
  rgb *colors;
  F* inh_map;     // mappa INH
  F* vvx_map;     // mappa velocita iniziale x
  F* vvy_map;     // mappa velocita iniziale y
  F* btm_map;     // mappa quote terreno
  F* man_map;     // mappa manning
  F* por_map;     // mappa porosita
  F* hlx_map;     // mappa hlc
  F* hly_map;     // mappa hlc
  F minx_map, maxx_map, miny_map, maxy_map, minv_map, maxv_map;
  F min_btm;
  F max_btm;
  F min_inh;
  F max_inh;
  F min_dep;
  F max_dep;
  F max_ini_uh;
  F max_uh;
  uchar4* host_info;     // mappa contorno (ricavata da BLN con scanline) (x=tipo di cella, y=4bit+4bit per orientazione bordi, z w sono i riferimenti a condizioni variabili nel tempo)
  int host_info_x; // dimensione utile matrice
  int host_info_y;


  F err_h;
  F err_uh;

  F dx,dy;

 /* vector <F4> sezioni;
  vector < vector <int2> > sezioni_port;
  vector < vector <int2> > sezioni_port_N;
  vector < vector <int2> > sezioni_port_E;
*/
  ///multires maps
  int bound_blocks;
  int tot_blocks;

  int esx,esy; // extended map size (max res)
  /// var appoggio per ricavare i blocchi di bcc
  uchar4** host_info_m;     // mappa contorno (ricavata da BLN con scanline) (x=tipo di cella, y=4bit+4bit per orientazione bordi, z w sono i riferimenti a condizioni variabili nel tempo)
  int* host_info_x_m; // dimensione utile matrice
  int* host_info_y_m;

  uchar4* host_info_multi;
  F4* host_grid_multi;
  unsigned char* host_grid_level_multi;
  ushort2* host_ofs_blocks; // riferimento per posizione blocco in multires rispetto  a griglia max ris
  F* host_man_map_multi;
  F* host_por_map_multi;
  F* host_hlx_map_multi;
  F* host_hly_map_multi;
  neigh_t* neigh;

};

struct global{
/////////// global parameters
//multirisoluzione

  int global_selecteddev;  // default scheda GPU da usare

  int global_multi;  // default singolo livello
  int global_level_bcc; // a che livello le boundary condition sono rasterizzate (>=1)

  F global_time;  // non caricato
  F global_time_start;
  F global_time_end;
  F global_CR;
  int global_limiter;
  F global_YEPS;
  int global_MUSCL;
  int global_ibinary;
  int global_al;
  int global_sr;
  F global_dtsoglia;
  F global_expon;
  F global_AWSDGM;
  F global_BWSDGM;
  int global_metodo;
  int global_niter;
  F global_dtoutput;
  F global_MAXDEP;
  F global_vel_eps;
  int global_printstat;
  F global_pend_far_field;
  int global_ordine;
  int global_debug;
  int global_por; // porosita' (di default disattivata)

  char base_path[256];

  // per display
  F display_dt;
  F display_vol;
  F display_global_time;

  //forzature multirisoluzione
  vector<F4> punti_m; //formato x,y, z=multires level [1--4], w=tappo (0 NO, 1 SI)

  //BRE

    F br_time_start; //tempo di inizio evoluzione della breccia
    F br_time_end; //tempo di fine della evoluzione della breccia
    F br_time_end_l; //tempo di fine della evoluzione della breccia per allargamento

    F B_iniz; // larghezza breccia al tempo time_iniz
    F B_fin; // larghezza della base della breccia al tempo time_fin
    F beta; //angolo di inclinazione delle sponde della breccia in radianti
    F z_iniz; //quota minima della breccia al tempo time_iniz
    F z_fin; //quota minima della breccia al tempo time_br

    F width_max; //massima larghezza della breccia (direzione normale a quella di evoluzione)
    //coordinate del centro della breccia
    F x_br;
    F y_br;

//versore che individua la direzione di propagazione della breccia (45 gradi)
    F dir_br_x;
    F dir_br_y;


} ;



#endif
