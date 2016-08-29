#include "bln_reader.h"

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "types.h"

using namespace std;


extern vector<cons> conss;  // elenco vincoli
vector<cons> conss_final;  // elenco vincoli interessanti
extern vector<polygon> polygons;

extern struct global g;
extern struct maps_ maps;


int assign_bc(F xp, F yp){

  // restituisce -1 se non c'era condizione (es bordo dovuto a btm)
  // in questo caso si mette muro!

  int numbc;  
  F seglung, alpha, xp_p, yp_p, xp_pp, yp_pp, dist, dist_min;
  
  dist_min=1.E16;
  numbc=-9999;

  for (int i=0;i<1;i++)//(int)polygons.size();i++)
    for (int k=0;k<(int)polygons[i].points.size();k++){
      //printf("Pol %d Ed %d: %f %f %d\n",i,k,polygons[i].points[k].x,polygons[i].points[k].y,polygons[i].edges[k]);
      int n=(int)polygons[i].points.size();
      F deltax=polygons[i].points[(k+1)%n].x-polygons[i].points[k].x;
      F deltay=polygons[i].points[(k+1)%n].y-polygons[i].points[k].y;
      seglung=pow(deltax*deltax+deltay*deltay,0.5f);
      //rotate
      alpha=atan2(deltay,deltax);
      xp_p = xp-polygons[i].points[k].x;
      yp_p = yp-polygons[i].points[k].y;
      xp_pp = xp_p*cos(alpha)+yp_p*sin(alpha);
      yp_pp = -xp_p*sin(alpha)+yp_p*cos(alpha);          
      dist=fabs(yp_pp);
      //printf("%d: %f<%f %f d %f, k %d\n",k,xp_pp,seglung, yp_pp,dist,k);
      if ( xp_pp>0 && xp_pp<seglung &&
	   dist<dist_min ){
	//printf("*\n");
	numbc=k;
	dist_min=dist;
      }  
    }
  
  if (numbc==-9999){
    //printf("ERROR: point with coordinate =( %f, %f )  has no boundary condition assigned\n",xp,yp);
    //exit(1);
    return -1;
  }      
  return numbc;
}

vector<bc> bc_list; // per ora non metto il riferimento a cella

void readBLNBCC(char* path, F* Z){//int nrow, int ncol, float mx, float Mx, float my, float My, uchar4*& maps.host_info, float* Z){
  readBLNBCC(path,0,Z,1); // multi e' a 0, quindi rasterizzazione normale, forza tutte bln
}

void readBLNBCC(char* path,int multi_level, F* Z, int force_BLN){//int nrow, int ncol, float mx, float Mx, float my, float My, uchar4*& maps.host_info, float* Z){

  /// parse input polygons and load maps.host_info
  /// ncol,nrow are input!
  /// Z mappa di battimetria (btm), per sovrascrivere punti alti
  /// restituisce il numero di constraints da copiare (
  F mx=maps.minx_map;
  F my=maps.miny_map;
  F Mx=maps.maxx_map;
  F My=maps.maxy_map;

  //  F* Z=maps.btm_map;

  int nrow=maps.nrows;
  int ncol=maps.ncols;

  if (ncol==0){
    printf("ncol=0!\n");
    return;
  }

  if (maps.host_info==NULL){
    //printf("alloc %dx%d\n",nrow,ncol);
    maps.host_info = (uchar4*)malloc(nrow* ncol * sizeof(uchar4));
  }

  char pathBLN[256];
  char pathBCC[256];

  strcpy(pathBLN,path);
  strcat(pathBLN,".BLN");
  std::ifstream ifs(pathBLN);

  int must_read_BCC=0;
  int nocond=0;
  if (ifs.fail()){
    printf("No boundary conditions! %s not found\n",pathBLN);
    polygons.resize(1);
    nocond=1;
    if (force_BLN && nocond && g.global_multi!=0) // senza bln, anche sulla risoluzione richiesta, non forzo
      force_BLN=0;
  }
  else{
  // load polygons
    polygons.clear();
  while (!ifs.eof()){
    int npoints;
    ifs >> npoints;
    if (!ifs.eof()){ // read polygon
      polygons.resize(polygons.size()+1);
      
      // parametro per sbiancare??
      //ifs >> polygons[polygons.size()-1].fill;
      polygons[polygons.size()-1].points.resize(npoints);
      polygons[polygons.size()-1].edges.resize(npoints);
      for (int i=0;i<npoints;i++){
	ifs >> polygons[polygons.size()-1].points[i].x;
	ifs >> polygons[polygons.size()-1].points[i].y;
	ifs >> polygons[polygons.size()-1].edges[i];
	if (polygons[polygons.size()-1].edges[i]==21)
	  polygons[polygons.size()-1].edges[i]=6;
	if (polygons[polygons.size()-1].edges[i]==2||
	    polygons[polygons.size()-1].edges[i]==6||
	    polygons[polygons.size()-1].edges[i]==3||
	    polygons[polygons.size()-1].edges[i]==5
	    )  // non e' muro
	  must_read_BCC=1;
      }	
    }
  }
  }
  ifs.close();
  
  printf("Polygons read\n");
  for (int i=0;i<(int)polygons.size();i++)
    for (int j=0;j<(int)polygons[i].points.size();j++)
      printf("Pol %d Ed %d: %f %f %d\n",i,j,polygons[i].points[j].x,polygons[i].points[j].y,polygons[i].edges[j]);
  

  int ncons=0;

  if (must_read_BCC){
    strcpy(pathBCC,path);
    strcat(pathBCC,".BCC");
    std::ifstream ifs(pathBCC);
    ifs >> ncons;
    conss.clear();
    conss.resize(ncons);
    for (int i=0;i<ncons;i++){
      ifs >> conss[i].edge;
      conss[i].edge--;  // internal reference is shifted to 0!
      int npoints;
      ifs >> npoints;
      conss[i].cons_edges.resize(npoints);
      for (int j=0;j<npoints;j++){
	conss[i].cons_edges[j].time=0;
	if (polygons[0].edges[conss[i].edge]==2 ||
            polygons[0].edges[conss[i].edge]==3 ||
            polygons[0].edges[conss[i].edge]==6
            ) /// funziona solo se c'e' un poligono	  
	  ifs >> conss[i].cons_edges[j].time;
	ifs >> conss[i].cons_edges[j].h;
	ifs >> conss[i].cons_edges[j].uh;
	ifs >> conss[i].cons_edges[j].vh;
      }	    
    }
    ifs.close();

    printf("Vincoli:\n");
    for (int i=0;i<(int)conss.size();i++){
      printf("edge %d ref to segm %d\n",i,conss[i].edge);
      for (int j=0;j<(int)conss[i].cons_edges.size();j++){
	printf("%f %f %f %f\n",
	       conss[i].cons_edges[j].time,
	       conss[i].cons_edges[j].h,
	       conss[i].cons_edges[j].uh,
	       conss[i].cons_edges[j].vh);
      }
    }
  }

  F dx=maps.dx; //1000
  F dy=maps.dy; //1000

  maps.host_info_x=ncol;
  maps.host_info_y=nrow;
  int scale=1<<(multi_level); //perchè shift?

  // con multi tolgo il padding di 1

  if (g.global_multi!=0){
    printf("Multiresolution scale for bounds: %d\n",scale);
    maps.host_info_x=(ncol-3)/scale+3; // sottraggo padding di 1 intorno a mappa, scalo (devo sottrarre 1 per divisione corretta) e riaggiungo 3
    maps.host_info_y=(nrow-3)/scale+3;
    ncol=maps.host_info_x;
    nrow=maps.host_info_y;
    printf("use size %d %d\n",maps.host_info_x,maps.host_info_y);
    //regolo in passo di rasterizzazione
    mx=mx+dx-(dx*scale);
    my=my+dy-(dy*scale); // scalo correttamente il padding
    dx*=(F)scale;
    dy*=(F)scale;
    if (force_BLN==0) // se non richiedo di rasterizzare bln, come se non ci fossero condizioni
      nocond=1;
  }
    printf("force bln %d, no cond %d\n",force_BLN,nocond);

  F dets[100];
  /// polygon rasterization 
  
  memset(maps.host_info,0,sizeof(uchar4)*ncol*nrow); //
  // Come modifico?
  vector<int2> pt_list;
  int n=(int)polygons[0].points.size();
  for (int k=0;k<n;k++){  
    int xmp=(polygons[0].points[k].x-mx)/dx;
    int xMp=(polygons[0].points[(k+1)%n].x-mx)/dx;
    if (xmp>xMp){
      xMp=(polygons[0].points[k].x-mx)/dx;
      xmp=(polygons[0].points[(k+1)%n].x-mx)/dx;
    }
    int ymp=(polygons[0].points[k].y-my)/dy;
    int yMp=(polygons[0].points[(k+1)%n].y-my)/dy;
    if (ymp>yMp){
      yMp=(polygons[0].points[k].y-my)/dy;
      ymp=(polygons[0].points[(k+1)%n].y-my)/dy;
    }/*
    xmp-=1;
    ymp-=1;
    xMp+=1;
    yMp+=1;
    if (xmp<=0) xmp=1;
    if (ymp<=0) ymp=1;
    if (xMp>=ncol) xMp=ncol-1;
    if (yMp>=nrow) yMp=nrow-1;*/
    for (int x=xmp;x<=xMp;x++)
      for (int y=ymp;y<=yMp;y++){
	int2 temp;
	temp.x=x;
	temp.y=y;
	if (maps.host_info[x + ncol * y].x==0){
	  maps.host_info[x + ncol * y].x=1;
	  pt_list.push_back(temp);
	}
      }      
  }
  
  //printf("buffer size %d\n",pt_list.size());

  for (int x=0;x<ncol;x++)
    for (int y=0;y<nrow;y++)
      maps.host_info[x + ncol * y].x=1; // temp per dire "non processato" --> 2

  vector<int2> queue; //out
  vector<int2> queue0; //in
  queue.reserve(3*pt_list.size()); //spazio enough/too much/insuff?
  queue0.reserve(3*pt_list.size());
  int2 foo; // rappresenta un punto. Bisogna cambiare il tipo per la tavoletta?
  printf("aa\n");

  // casi speciali
  if (nocond) {
    if (g.global_multi==0 || force_BLN)
      { // se c'e' multirisoluzione, il muro diventa implicito

      for (int i=1;i<ncol-1;i++){
	maps.host_info[1 * ncol + i].x=0; // ---> 1
	foo.x=i;
	foo.y=1;
	queue0.push_back(foo);
	maps.host_info[(nrow-2) * ncol + i].x=0;
	foo.x=i;
	foo.y=nrow-2;
	queue0.push_back(foo);
      }
      for (int j=1;j<nrow-1;j++){
	maps.host_info[j * ncol + 1].x=0;
	foo.x=1;
	foo.y=j;
	queue0.push_back(foo);
	maps.host_info[j * ncol + (ncol-2)].x=0;
	foo.x=ncol-2;
	foo.y=j;
	queue0.push_back(foo);
      }

      for (int i=0;i<ncol;i++){
	maps.host_info[0 * ncol + i].x=BIT_EXTERN; //---> 0
	foo.x=i;
	foo.y=0;
	queue.push_back(foo);
	maps.host_info[(nrow-1) * ncol + i].x=BIT_EXTERN;
	foo.x=i;
	foo.y=nrow-1;
	queue.push_back(foo);
      }
      for (int j=0;j<nrow;j++){
	maps.host_info[j * ncol + 0].x=BIT_EXTERN;
	foo.x=0;
	foo.y=j;
	queue.push_back(foo);
	maps.host_info[j * ncol + (ncol-1)].x=BIT_EXTERN;
	foo.x=ncol-1;
	foo.y=j;
	queue.push_back(foo);
      }
    }
    else{
      printf("multires senza condizioni\n");
    }
  }
  else {
    for (int pt_it=0;pt_it<pt_list.size();pt_it++) {
      int i=pt_list[pt_it].x;
      int j=pt_list[pt_it].y;
      F xp=mx+(i)*dx; // scalo per livello
      F yp=my+(j)*dy;
      //printf("%f %f: \n",xp,yp);
      int numbc=0;
      int i_est=0;			
      int n=(int)polygons[0].points.size();
      for (int k=0;k<n;k++){  
	//	printf("test %f %f: %f %f, %f %f\n",xp,yp,polygons[0].points[k].x,polygons[0].points[k].y,polygons[0].points[(k+1)%n].x,polygons[0].points[(k+1)%n].y);

	if ( (polygons[0].points[k].x < xp && polygons[0].points[(k+1)%n].x >= xp) ||
	     (polygons[0].points[k].x >= xp && polygons[0].points[(k+1)%n].x < xp) ){
	  numbc=1;
	  F polyXi=polygons[0].points[(k)%n].x;
	  F polyXj=polygons[0].points[(k+1)%n].x;
	  F polyYi=polygons[0].points[(k)%n].y;
	  F polyYj=polygons[0].points[(k+1)%n].y;
	  if (polyYi+(xp-polyXi)/(polyXj-polyXi)*(polyYj-polyYi)<yp)
	    i_est=1-i_est;
	  //printf("%d %d %f %f --> %d\n",i,j,polyYi+(xp-polyXi)/(polyXj-polyXi)*(polyYj-polyYi),yp,i_est);
	}
      }
      if (i_est==0 || numbc==0){ //!condizione per celle con xp < min(xpoint) o xp > max(xpoint)
	maps.host_info[j * ncol + i].x=BIT_EXTERN; //!cella tappo
	foo.x=i;
	foo.y=j;
	queue.push_back(foo);
      }
      else{
	maps.host_info[j * ncol + i].x=0; //!cella tappo
	foo.x=i;
	foo.y=j;
	queue0.push_back(foo);
      }
    }

    // aggiungo a mano il bordo esterno
    if (g.global_multi==0 || force_BLN){ // se c'e' multirisoluzione, il muro diventa implicito
      for (int i=0;i<ncol;i++){
	maps.host_info[0 * ncol + i].x=BIT_EXTERN;
	foo.x=i;
	foo.y=0;
	queue.push_back(foo);
	maps.host_info[(nrow-1) * ncol + i].x=BIT_EXTERN;
	foo.x=i;
	foo.y=nrow-1;
	queue.push_back(foo);
      }
      for (int j=0;j<nrow;j++){
	maps.host_info[j * ncol + 0].x=BIT_EXTERN;
	foo.x=0;
	foo.y=j;
	queue.push_back(foo);
	maps.host_info[j * ncol + (ncol-1)].x=BIT_EXTERN;
	foo.x=ncol-1;
	foo.y=j;
	queue.push_back(foo);
      }        
    }
  }

  // ora completa con simil floodfill
  while (queue.size()>0){
    int2 foo=queue[queue.size()-1];
    //    printf("--%3d-1> %d %d: %d\n",queue.size(),foo.x,foo.y,maps.host_info[foo.y * ncol + foo.x].x);
    queue.erase(queue.end()-1);
    maps.host_info[foo.y * ncol + foo.x].x=BIT_EXTERN;
    //aggiungi neigh
    int2 n;
    n.x=foo.x+1;n.y=foo.y+0;if (n.x<ncol && maps.host_info[n.y * ncol + n.x].x==1){
      maps.host_info[n.y * ncol + n.x].x=BIT_EXTERN;
      queue.push_back(n);
    }
    n.x=foo.x-1;n.y=foo.y+0;if (n.x>=0 && maps.host_info[n.y * ncol + n.x].x==1){
      maps.host_info[n.y * ncol + n.x].x=BIT_EXTERN;
      queue.push_back(n);
    }
    n.x=foo.x;n.y=foo.y+1;if (n.y<nrow && maps.host_info[n.y * ncol + n.x].x==1){
      maps.host_info[n.y * ncol + n.x].x=BIT_EXTERN;
      queue.push_back(n);
    }
    n.x=foo.x;n.y=foo.y-1;if (n.y>=0 && maps.host_info[n.y * ncol + n.x].x==1){
      maps.host_info[n.y * ncol + n.x].x=BIT_EXTERN;
      queue.push_back(n);
    }
  }

  while (queue0.size()>0){
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
  }
  //---------------------------------------------------------------------------------------///
  if (0) 
  for (int x=0;x<nrow;x++){
    for (int y=0;y<ncol;y++)
      printf("%2d ",maps.host_info[x + ncol * y]);
    printf("\n");
  }


  // riporta tappi da z, considerando anche resampling (basta un tappo per forzare a tappo tutta la patch!)
  for (int i=1;i<ncol-1;i++)
    for (int j=1;j<nrow-1;j++){
      if (Z[j * ncol + i]>10E10)
	maps.host_info[j * ncol + i].x=BIT_TAPPO;
    }

  if(0){
    for (int y=0;y<30;y++){
      for (int x=0;x<ncol;x++)
	if (maps.host_info[x + ncol * y].x==BIT_EXTERN)
	  printf("x");
	else{
	  if (maps.host_info[x + ncol * y].x==BIT_TAPPO)
	    printf("T");
	  else
	    printf(".");
	}
      printf("\n");
    }
  }
  //  printf("fix\n");
  // fix 
  int modified=1;
  //  if (!nocond)
    while (modified){ // aggiunto punto fisso, altrimenti non si garantisce che la propagazione abbia finito
      modified=0;
      for (int j=0;j<nrow;j++)
	for (int i=0;i<ncol;i++){        
	  int num_bou=0; //contatore del numero di cc per ciascuna cella
	if (maps.host_info[i + ncol * j].x<BIT_TAPPO){// non e' esterna
	  if ((i-1)>=0 && maps.host_info[(i-1) + ncol * j].x>=BIT_TAPPO) num_bou=num_bou+1;
	  if ((i+1)<ncol && maps.host_info[(i+1) + ncol * j].x>=BIT_TAPPO) num_bou=num_bou+1;
	  if ((j+1)<nrow && maps.host_info[(i) + ncol *( j+1)].x>=BIT_TAPPO) num_bou=num_bou+1;
	  if ((j-1)>=0 && maps.host_info[(i) + ncol * (j-1)].x>=BIT_TAPPO) num_bou=num_bou+1;           
	  if (num_bou>2){
	    maps.host_info[i + ncol * j].x=BIT_EXTERN;
	    //Z[i + ncol * j]=1.7014E+38;
	    //printf("DBG: tappo %d %d, lev %d\n",i,j,multi_level);
	    // aggiorno modifica nella mappa btm globale (ricavo la giusta risoluzione
	    int idx=((j-1)*scale+1)* maps.ncols +((i-1)*scale+1);
	      maps.btm_map[idx]=1.7014E+38;
	    modified=1;
	  }
	  if (num_bou>0 && num_bou<=2){ //cella di contorno (1 o 2 vicini)
	    maps.host_info[i + ncol * j].x=2;
	  }	
	  if (num_bou==0 && (multi_level>1 || Z[j * ncol + i]<10E10)
	      ){ //cella interna
	    maps.host_info[i + ncol * j].x=0;
	    //ncell_int=ncell_int+1
	    //         i_int(ncell_int)=i
	    //     j_int(ncell_int)=j     
	  }           
	}
      }
  }
  if (1==0)
    for (int j=0;j<nrow;j+=1){
      for (int i=0;i<ncol;i++){
	//if (maps.host_info[i + ncol * j].x==BIT_EXTERN)
	  printf("%2d.",maps.host_info[i + ncol * j].x);
	  //  else
	  //printf(". ");
	//	printf("%3.2f ",Z[i + ncol * j]);
      }
      printf("\n");
    }
  //flag per le condizioni al contorno sulle celle: 
  // se iflag_bc_w=0 no bc west
  // se iflag_bc_w=1 condizione di muro west
  // se iflag_bc_w=2 condizione di inflow west
  // se iflag_bc_w=3 condizione di outflow west
  // se iflag_bc_w=4 condizione di far field west

  int ncell_bc=0; //numero di celle con condizioni al contorno

  /*
  struct bc{
    ifl_bc_w=0; // tipo di condizione
    ifl_bc_e=0;
    ifl_bc_n=0;
    ifl_bc_s=0;
    numbc_w=0;  // riferimento a numero di segmento polilinea caricata
    numbc_e=0;
    numbc_n=0;
    numbc_s=0;
  };
  */
  /*
  dbgx=1371;
  dbgy=903;
  for (int x=dbgx-2;x<=dbgx+2;x++){
    for (int y=dbgy-2;y<=dbgy+2;y++)
      printf("%d %d %2d ",x,y,maps.host_info[x + ncol * y]);
    printf("\n");
  }
  */

  if (1==0)
    for (int y=0;y<nrow;y++){
      for (int x=0;x<ncol;x++)
	printf("%2d ",maps.host_info[x + ncol * y].x);
      printf("\n");
    }

  int dbgx=ncol-2;
  int dbgy=nrow-2;


  // assumo che i vincoli speciali siano solo sul primo poligono 
  // if (!nocond)
  for (int j=0;j<nrow;j++)
    for (int i=0;i<ncol;i++){        
      if (maps.host_info[j * ncol + i].x==2){ // cella bordo
	int if_bc=0; 
	bc cell;

	maps.host_info[j * ncol + i].x=0; // resetto, ora contano i bit settati (che saltano 1, riservato per cella fuori)

	cell.ifl_bc_w=0;cell.numbc_w=0;
	cell.ifl_bc_e=0;cell.numbc_e=0;
	cell.ifl_bc_n=0;cell.numbc_n=0;
	cell.ifl_bc_s=0;cell.numbc_s=0;
	

	//	printf("%d %d %d\n",i,j,maps.host_info[i + ncol * j]);

	F xp,yp;
	
	if ((i-1)>=0 && (maps.host_info[j * ncol +(i-1)].x==BIT_EXTERN || maps.host_info[j * ncol +(i-1)].x==BIT_TAPPO)){
	  if_bc=if_bc+1; 
	  //posizione dell'intercella
	  xp=mx+(i-0.5)*dx;
	  yp=my+(j)*dy;
	  int kk=assign_bc(mx+(i-0.5)*dx, my+(j)*dy); //!kk è il numero della condizione al contorno per la cella i,j
	  if (kk<0 || maps.host_info[j * ncol +(i-1)].x==BIT_TAPPO){
	    cell.ifl_bc_w=1;
	    cell.numbc_w=0;
	  }
	  else{
	    cell.ifl_bc_w=polygons[0].edges[kk]; //tipo di condizine al contorno per la cella i,j sul lato west (1 muro, 2 inflow, 3 outflow, 4 far field)
	    // se condizione 1 -> non c'e' segmento con dati, se >1 cerco il codice del segmento bcc
	    int found=-1;
	    if (polygons[0].edges[kk]>1){
	      for (int j=0;j<ncons;j++)
		if (conss[j].edge==kk)
		  found=j;
	      if (found==-1 && polygons[0].edges[kk]!=1 &&polygons[0].edges[kk]!=4)
		printf("ERROR: reference segment not found\n");
	    }
	    else 
	      found=0;
	    cell.numbc_w=found; //condizione variabile nel tempo
	  }
	  //Z[(i-1) + ncol * j]=Z[i + ncol * j]; //devo aggiornare la quota per non avere tappi nelle celle fantasma
	}               
	if ((i+1)<ncol && (maps.host_info[(i+1) + ncol * j].x==BIT_EXTERN || maps.host_info[j * ncol +(i+1)].x==BIT_TAPPO)){
	  if_bc=if_bc+1;
	  //posizione dell'intercella
	  xp=mx+(i+0.5)*dx;
	  yp=my+(j)*dy;
	  int kk=assign_bc(mx+(i+0.5)*dx, my+(j)*dy); //!kk è il numero della condizione al contorno per la cella i,j
	  //	  printf("east %d %d: %d\n",i,j,kk);
	  if (kk<0 || maps.host_info[j * ncol +(i+1)].x==BIT_TAPPO){
	    cell.ifl_bc_e=1;
	    cell.numbc_e=0;
	  }
	  else{
	    cell.ifl_bc_e=polygons[0].edges[kk]; //tipo di condizine al contorno per la cella i,j sul lato west (1 muro, 2 inflow, 3 outflow, 4 far field)
	    //printf("east --> %d\n",polygons[0].edges[kk]);
	    // se condizione 1 -> non c'e' segmento con dati, se >1 cerco il codice del segmento bcc
	    int found=-1;
	    if (polygons[0].edges[kk]>1){
	      for (int j=0;j<ncons;j++)
		if (conss[j].edge==kk)
		  found=j;
	      if (found==-1 && polygons[0].edges[kk]!=1 && polygons[0].edges[kk]!=4)
		printf("ERROR: reference segment not found\n");
	    }
	    else 
	      found=0;
	    cell.numbc_e=found; //condizione variabile nel tempo
	  }
	  //Z[(i+1) + ncol * j]=Z[i + ncol * j]; //devo aggiornare la quota per non avere tappi nelle celle fantasma
	}               
	if ((j+1)<nrow && (maps.host_info[i + ncol *( j + 1)].x==BIT_EXTERN || maps.host_info[i + ncol *( j + 1)].x==BIT_TAPPO)){
	  if_bc=if_bc+1;
	  //posizione dell'intercella
	  xp=mx+(i)*dx;
	  yp=my+(j+0.5)*dy;
	  int kk=assign_bc(mx+(i)*dx, my+(j+0.5)*dy); //!kk è il numero della condizione al contorno per la cella i,j
	  if (kk<0 || maps.host_info[i + ncol *( j + 1)].x==BIT_TAPPO){
	    cell.ifl_bc_n=1;
	    cell.numbc_n=0;
	  }
	  else{
	    cell.ifl_bc_n=polygons[0].edges[kk]; //tipo di condizine al contorno per la cella i,j sul lato west (1 muro, 2 inflow, 3 outflow, 4 far field)
	    // se condizione 1 -> non c'e' segmento con dati, se >1 cerco il codice del segmento bcc
	    int found=-1;
	    if (polygons[0].edges[kk]>1){
	      for (int j=0;j<ncons;j++)
		if (conss[j].edge==kk)
		  found=j;
	      if (found==-1 && polygons[0].edges[kk]!=1 &&polygons[0].edges[kk]!=4)
		printf("ERROR: reference segment not found\n");
	    }
	    else 
	      found=0;
	    cell.numbc_n=found; //condizione variabile nel tempo
	  }
	  //Z[i + ncol *(j + 1)]=Z[i + ncol * j]; //devo aggiornare la quota per non avere tappi nelle celle fantasma
	}

	if ((j-1)>=0 && (maps.host_info[i + ncol * (j - 1)].x==BIT_EXTERN || maps.host_info[i + ncol * (j - 1)].x==BIT_TAPPO)){
	  if_bc=if_bc+1;
	  xp=mx+(i)*dx;
	  yp=my+(j-0.5)*dy;
	  int kk=assign_bc(mx+(i)*dx, my+(j-0.5)*dy); //!kk è il numero della condizione al contorno per la cella i,j
	  if (kk<0 || maps.host_info[i + ncol * (j - 1)].x==BIT_TAPPO){
	    cell.ifl_bc_s=1;	
	    cell.numbc_s=0;
	  }
	  else{
	    cell.ifl_bc_s=polygons[0].edges[kk]; //tipo di condizine al contorno per la cella i,j sul lato west (1 muro, 2 inflow, 3 outflow, 4 far field)
	    // se condizione 1 -> non c'e' segmento con dati, se >1 cerco il codice del segmento bcc
	    int found=-1;
	    if (polygons[0].edges[kk]>1){
	      for (int j=0;j<ncons;j++)
		if (conss[j].edge==kk)
		  found=j;
	      if (found==-1 && polygons[0].edges[kk]!=1 &&polygons[0].edges[kk]!=4)
		printf("ERROR: reference segment not found\n");
	    }
	    else 
	      found=0;
	    cell.numbc_s=found; //condizione variabile nel tempo
	  }
	  //Z[i + ncol *( j - 1)]=Z[i + ncol * j]; //devo aggiornare la quota per non avere tappi nelle celle fantasma
	}

	if (if_bc!=0){
	  // cerca se gia' presente
	  int found=-1;
	  for (int i1=0;i1<(int)bc_list.size();i1++){
	    bc cellc=bc_list[i1];
	    if (cell.ifl_bc_w==cellc.ifl_bc_w &&
		cell.ifl_bc_e==cellc.ifl_bc_e &&
		cell.ifl_bc_s==cellc.ifl_bc_s &&
		cell.ifl_bc_n==cellc.ifl_bc_n &&
		cell.numbc_n ==cellc.numbc_n &&
		cell.numbc_s ==cellc.numbc_s &&
		cell.numbc_e ==cellc.numbc_e &&
		cell.numbc_w ==cellc.numbc_w)
	      found=i1;
	  }

	  if (found==-1){
	    bc_list.push_back(cell);
	    ncell_bc=ncell_bc+1;
	    found=(int)bc_list.size()-1;
	  }

	  //	  maps.host_info[i * ncol + j].x=found; // riferimento ai dati 

	  // parse and store .y .z .w
	  // in .x store flags 4 bordi
	  // in .y ci sono 4 bit per ciascun <= 2 bordo

	  maps.host_info[i + ncol * j].x=0;
	  if (bc_list[found].ifl_bc_w!=0){
	    maps.host_info[i + ncol * j].x|=BIT_W;
	  }
	  if (bc_list[found].ifl_bc_e!=0){
	    maps.host_info[i + ncol * j].x|=BIT_E;
	  }
	  if (bc_list[found].ifl_bc_s!=0){
	    maps.host_info[i + ncol * j].x|=BIT_S;
	  }
	  if (bc_list[found].ifl_bc_n!=0){
	    maps.host_info[i + ncol * j].x|=BIT_N;
	  }
	  
	  if (1==0)
	  printf("%d, %d %d, %d %d, %d %d, %d %d\n",maps.host_info[i + ncol * j].x,
		 bc_list[found].ifl_bc_n,bc_list[found].numbc_n,
		 bc_list[found].ifl_bc_s,bc_list[found].numbc_s,
		 bc_list[found].ifl_bc_e,bc_list[found].numbc_e,
		 bc_list[found].ifl_bc_w,bc_list[found].numbc_w);

	  switch(maps.host_info[i + ncol * j].x){
	  case BIT_N: // solo nord
	    maps.host_info[i + ncol * j].y=bc_list[found].ifl_bc_n;
	    maps.host_info[i + ncol * j].z=bc_list[found].numbc_n;
	    break;
	  case BIT_S: // solo sud
	    maps.host_info[i + ncol * j].y=bc_list[found].ifl_bc_s;
	    maps.host_info[i + ncol * j].z=bc_list[found].numbc_s;
	    break;
	  case BIT_W: // solo w
	    maps.host_info[i + ncol * j].y=bc_list[found].ifl_bc_w;
	    maps.host_info[i + ncol * j].z=bc_list[found].numbc_w;
	    break;
	  case BIT_E: // solo e
	    maps.host_info[i + ncol * j].y=bc_list[found].ifl_bc_e;
	    maps.host_info[i + ncol * j].z=bc_list[found].numbc_e;
	    break;
	  case BIT_N+BIT_S: //due bordi
	    maps.host_info[i + ncol * j].x=BIT_N+16*BIT_S;
	    maps.host_info[i + ncol * j].y=bc_list[found].ifl_bc_n+16*bc_list[found].ifl_bc_s;
	    maps.host_info[i + ncol * j].z=bc_list[found].numbc_n;
	    maps.host_info[i + ncol * j].w=bc_list[found].numbc_s;
	    break;
	  case BIT_N+BIT_E: //due bordi
	    maps.host_info[i + ncol * j].x=BIT_N+16*BIT_E;
	    maps.host_info[i + ncol * j].y=bc_list[found].ifl_bc_n+16*bc_list[found].ifl_bc_e;
	    maps.host_info[i + ncol * j].z=bc_list[found].numbc_n;
	    maps.host_info[i + ncol * j].w=bc_list[found].numbc_e;
	    break;
	  case BIT_N+BIT_W: //due bordi
	    maps.host_info[i + ncol * j].x=BIT_N+16*BIT_W;
	    maps.host_info[i + ncol * j].y=bc_list[found].ifl_bc_n+16*bc_list[found].ifl_bc_w;
	    maps.host_info[i + ncol * j].z=bc_list[found].numbc_n;
	    maps.host_info[i + ncol * j].w=bc_list[found].numbc_w;
	    break;
	  case BIT_S+BIT_W: //due bordi
	    maps.host_info[i + ncol * j].x=BIT_S+16*BIT_W;
	    maps.host_info[i + ncol * j].y=bc_list[found].ifl_bc_s+16*bc_list[found].ifl_bc_w;
	    maps.host_info[i + ncol * j].z=bc_list[found].numbc_s;
	    maps.host_info[i + ncol * j].w=bc_list[found].numbc_w;
	    break;
	  case BIT_S+BIT_E: //due bordi
	    maps.host_info[i + ncol * j].x=BIT_S+16*BIT_E;
	    maps.host_info[i + ncol * j].y=bc_list[found].ifl_bc_s+16*bc_list[found].ifl_bc_e;
	    maps.host_info[i + ncol * j].z=bc_list[found].numbc_s;
	    maps.host_info[i + ncol * j].w=bc_list[found].numbc_e;
	    break;
	  case BIT_W+BIT_E: //due bordi
	    maps.host_info[i + ncol * j].x=BIT_W+16*BIT_E;
	    maps.host_info[i + ncol * j].y=bc_list[found].ifl_bc_w+16*bc_list[found].ifl_bc_e;
	    maps.host_info[i + ncol * j].z=bc_list[found].numbc_w;
	    maps.host_info[i + ncol * j].w=bc_list[found].numbc_e;
	    break;
	    
	      
	  }

	}
      }
    }
  
  for (int j=0;j<nrow;j++)
    for (int i=0;i<ncol;i++){
      if (maps.host_info[i + ncol * j].x==BIT_TAPPO)
        maps.host_info[i + ncol * j].x=BIT_EXTERN;
     //if (maps.host_info[i + ncol * j].x==BIT_EXTERN)
	//Z[i + ncol * j]=1.7014E+38;
 	}
  
  /*  for (int j=0;j<nrow;j++)
    for (int i=0;i<ncol;i++)
	maps.host_info[j * ncol + i].x=j%128;
*/
  
  int ct=0;


  char path1[256];
  char pathcode[2];

  if (g.global_debug==1 && multi_level==0)
  for (int side=0;side<4;side++){
    strcpy(path1,path);
    strcat(path1,"-debug.");
    if (side==0)
      pathcode[0]='N';
    if (side==1)
      pathcode[0]='S';
    if (side==2)
      pathcode[0]='E';
    if (side==3)
      pathcode[0]='W';
    pathcode[1]=0;
    strcat(path1,pathcode);

    FILE* fo=fopen(path1,"w+");
    fprintf(fo,"DSAA\n");
    fprintf(fo,"%d %d\n",ncol-2,nrow-2);
    F deltax=dx;
    F deltay=dy;
    fprintf(fo,"%f %f\n",maps.minx_map+deltax, maps.maxx_map-deltax);
    fprintf(fo,"%f %f\n",maps.miny_map+deltay, maps.maxy_map-deltay);
    F minv=10E10;
    F maxv=0;
    for(int j = 1; j < nrow-1; j ++)  
      for(int i = 1; i < ncol-1; i ++){
	int scrivo=0;
	if (maps.host_info[i + ncol * j].x==BIT_EXTERN)
	  scrivo=maps.host_info[i + ncol * j].x;
	else{
	  int val=0;
	  for (int k=0;k<2;k++){
	    char info,vincolo;
	    if (k==0){
	      info=maps.host_info[i + ncol * j].x & 15; 
	      vincolo=maps.host_info[i + ncol * j].y & 15; 
	    }else{
	      info=maps.host_info[i + ncol * j].x >> 4;	      
	      vincolo=maps.host_info[i + ncol * j].y >> 4;	      
	    }

	    if ((side==0 && info==BIT_N)||
		(side==1 && info==BIT_S)||
		(side==2 && info==BIT_E)||
		(side==3 && info==BIT_W))
	      val=vincolo;
	  }
	  scrivo=val;
	}
	if (minv>scrivo) minv=scrivo;
	if (maxv<scrivo) maxv=scrivo;
      }

    fprintf(fo,"%f %f\n",minv, maxv);
    for(int j = 1; j < nrow-1; j ++){
      for(int i = 1; i < ncol-1; i ++){
	int scrivo=0;
	if (maps.host_info[i + ncol * j].x==BIT_EXTERN)
	  scrivo=maps.host_info[i + ncol * j].x;
	else{
	  int val=0;
	  for (int k=0;k<2;k++){
	    char info,vincolo;
	    if (k==0){
	      info=maps.host_info[i + ncol * j].x & 15; 
	      vincolo=maps.host_info[i + ncol * j].y & 15; 
	    }else{
	      info=maps.host_info[i + ncol * j].x >> 4;	      
	      vincolo=maps.host_info[i + ncol * j].y >> 4;	      
	    }

	    if ((side==0 && info==BIT_N)||
		(side==1 && info==BIT_S)||
		(side==2 && info==BIT_E)||
		(side==3 && info==BIT_W))
	      val=vincolo;
	  }
	  scrivo=val;
	}
	fprintf(fo,"%.1f ",(F)scrivo);
      }
      fprintf(fo,"\n");
    }
    fclose(fo);
  }

  printf("map %d %d\n",ncol,nrow);
  if (0)// && g.global_multi && multi_level>0)
    for(int j = nrow-2; j>=1; j --){
      printf("%3d: ",j);
      for(int i = 1; i < ncol-1; i ++){
	int v=maps.host_info[i + ncol * j].x;
	if (v!=0 && v!=BIT_EXTERN)
	  printf("C");
	else{
	  if (v==0)
	    printf(" ");
	  else
	    printf("!");
	}
      }
      printf("\n");
    }
    
  
  
  if (1==0)
  for (int i=0;i<bc_list.size();i++){
    bc cell=bc_list[i];
    printf("%d %d %d %d %d %d %d %d\n",
	   cell.ifl_bc_w,cell.numbc_w,
	   cell.ifl_bc_e,cell.numbc_e,
	   cell.ifl_bc_n,cell.numbc_n,
	   cell.ifl_bc_s,cell.numbc_s);
  }
  printf("done\n");
}
