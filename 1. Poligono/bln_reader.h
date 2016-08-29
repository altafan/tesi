#pragma once
//#include "types.h"

#ifndef bln_reader
#define bln_reader

typedef struct bc_{
    int ifl_bc_w; // tipo di condizione
    int ifl_bc_e;
    int ifl_bc_n;
    int ifl_bc_s;
    int numbc_w;  // riferimento a numero di segmento polilinea caricata
    int numbc_e;
    int numbc_n;
    int numbc_s;
  } bc;

void readBLNBCC(char* path);
//void readBLNBCC(char* path, int level, F* Z,int );

#endif
