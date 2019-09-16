#include "mem_seg_levelset120214sf.h"
#include <mpi.h>

DIR* dir;
struct dirent* dp;
int count_file = 0;
char ddirname[MAX_CHAR_NUM] =    "/home/niwayama/memsegmentation/140212sparsefield/140212sfm/img/";
char diroutname[MAX_CHAR_NUM] =  "/home/niwayama/memsegmentation/140212sparsefield/140212sfm/dst/";
char FILENAME[TOTFILENAME][MAX_CHAR_NUM];

double dd_out[totalpixels];
double d_phi[totalpixels];
double dd_g[totalpixels];
double d_label[totalpixels];
double d_mask[totalpixels];
//Global Vvariables for Eichonal
//solveMatrixEq

int NumDetectedNuc;
int CircleMap[WBAND+1][300*(WBAND+1)][3];
int NCircleMap[WBAND+1];
int i_ResetInval;
int i_NucX[MAXCELLNUM];
int i_NucY[MAXCELLNUM];
int i_NucZ[MAXCELLNUM];

int GRIDSIZE_Z;
int GRIDSIZE_Y;
int GRIDSIZE_X;
long    dims[5];
long   mdims[3];
#define USAGE(argv)\
{\
  fprintf(stderr,"Usage %s imagefile \n",argv[0]);\
  exit(-1);\
}
int Indix(int ii, int kk, int jj){
  int idx=(int)(ii*DIMXY+jj*DIMY+kk);
  return idx;
}
/*
void InitializeCircleMap(){  
  int i, x, y, z;
  for(i=0;i<=WBAND;i++){
    NCircleMap[i]=0;
  }
  for(x=-WBAND;x<=WBAND;x++){
    for(y=-WBAND;y<=WBAND;y++){
      for(z=-WBAND;z<=WBAND;z++){
	double dist2 = double(x*x+y*y+z*z);
	if(dist2<0.0){
	  dist2=0.0;
	}
	int d =int(sqrt(dist2));
	if(d<=WBAND){
	  CircleMap[d][NCircleMap[d]][0]=z;
	  CircleMap[d][NCircleMap[d]][1]=y;
	  CircleMap[d][NCircleMap[d]][2]=x;
	  NCircleMap[d]++;
	  //printf("%dNC\n",NCircleMap[d]);
	  if(NCircleMap[d]>300*(WBAND+1)-2){
	    // printf("%dNC out\n",NCircleMap[d]);
	    exit(-1);
	  }
	}
      }
    }
  }
  for(i=0;i<=WBAND;i++){
    // printf("%dNC\n",NCircleMap[i]);
  }
}

void ReInitialization(){
  for(int h=0; h<GRIDSIZE_Z;h++){
    for(int i=0; i<GRIDSIZE_Y;i++){
      for(int j=0; j<GRIDSIZE_X;j++){
	if(d_phi[Indix(h,i,j)]>0.0){
	  d_phi[Indix(h,i,j)]=1.0;
	}
	else{
	  d_phi[Indix(h,i,j)]=-0.000001;
	}
      }
    }
  }
}
*/

void Calcg(){
  for(int h=0; h<GRIDSIZE_Z;h++){
    for(int i=0; i<GRIDSIZE_Y;i++){
      for(int j=0; j<GRIDSIZE_X;j++){
	dd_g[Indix(h,i,j)] =1.0;
	if(h>=GRIDSIZE_Z-2||i>=GRIDSIZE_Y-2||j>=GRIDSIZE_X-2||h==0||i==0||j==0){
	  dd_g[Indix(h,i,j)]=1.0;
	}
	else{
	  double d_Uppo, d_Upmo, d_Umpo, d_Ummo;
	  double d_Upop, d_Upom, d_Umop, d_Umom;
	  double d_Uopp, d_Uopm, d_Uomp, d_Uomm;
	  if(h>0&&h<GRIDSIZE_Z-1&&i>0&&i<GRIDSIZE_Y-1){
	    d_Uppo=(dd_out[Indix(h,i,j)]+dd_out[Indix(h,i+1,j)]+dd_out[Indix(h+1,i,j)]+dd_out[Indix(h+1,i+1,j)])/4.0;
	    d_Upmo=(dd_out[Indix(h,i,j)]+dd_out[Indix(h,i-1,j)]+dd_out[Indix(h+1,i,j)]+dd_out[Indix(h+1,i-1,j)])/4.0;
	    d_Umpo=(dd_out[Indix(h,i,j)]+dd_out[Indix(h,i+1,j)]+dd_out[Indix(h-1,i,j)]+dd_out[Indix(h-1,i+1,j)])/4.0;
	    d_Ummo=(dd_out[Indix(h,i,j)]+dd_out[Indix(h,i-1,j)]+dd_out[Indix(h-1,i,j)]+dd_out[Indix(h-1,i-1,j)])/4.0;
	  }
	  else{
	    d_Uppo=(dd_out[Indix(h,i,j)]);
	    d_Upmo=(dd_out[Indix(h,i,j)]);
	    d_Umpo=(dd_out[Indix(h,i,j)]);
	    d_Ummo=(dd_out[Indix(h,i,j)]);
	  }
	  if(h>0&&h<GRIDSIZE_Z-1&&j>0&&j<GRIDSIZE_X-1){
	    d_Upop=(dd_out[Indix(h,i,j)]+dd_out[Indix(h,i,j+1)]+dd_out[Indix(h+1,i,j)]+dd_out[Indix(h+1,i,j+1)])/4.0;
	    d_Upom=(dd_out[Indix(h,i,j)]+dd_out[Indix(h,i,j-1)]+dd_out[Indix(h+1,i,j)]+dd_out[Indix(h+1,i,j-1)])/4.0;
	    d_Umop=(dd_out[Indix(h,i,j)]+dd_out[Indix(h,i,j+1)]+dd_out[Indix(h-1,i,j)]+dd_out[Indix(h-1,i,j+1)])/4.0;
	    d_Umom=(dd_out[Indix(h,i,j)]+dd_out[Indix(h,i,j-1)]+dd_out[Indix(h-1,i,j)]+dd_out[Indix(h-1,i,j-1)])/4.0;
	  }
	  else{
	    d_Upop=(dd_out[Indix(h,i,j)]);
	    d_Upom=(dd_out[Indix(h,i,j)]);
	    d_Umop=(dd_out[Indix(h,i,j)]);
	    d_Umom=(dd_out[Indix(h,i,j)]);
	  }
	  if(i>0&&i<GRIDSIZE_Y-1&&j>0&&j<GRIDSIZE_X-1){
	    d_Uopp=(dd_out[Indix(h,i,j)]+dd_out[Indix(h,i,j+1)]+dd_out[Indix(h,i+1,j)]+dd_out[Indix(h,i+1,j+1)])/4.0;
	    d_Uopm=(dd_out[Indix(h,i,j)]+dd_out[Indix(h,i,j-1)]+dd_out[Indix(h,i+1,j)]+dd_out[Indix(h,i+1,j-1)])/4.0;
	    d_Uomp=(dd_out[Indix(h,i,j)]+dd_out[Indix(h,i,j+1)]+dd_out[Indix(h,i-1,j)]+dd_out[Indix(h,i-1,j+1)])/4.0;
	    d_Uomm=(dd_out[Indix(h,i,j)]+dd_out[Indix(h,i,j-1)]+dd_out[Indix(h,i-1,j)]+dd_out[Indix(h,i-1,j-1)])/4.0;
	  }
	  else{
	    d_Uopp=(dd_out[Indix(h,i,j)]);
	    d_Uopm=(dd_out[Indix(h,i,j)]);
	    d_Uomp=(dd_out[Indix(h,i,j)]);
	    d_Uomm=(dd_out[Indix(h,i,j)]);
	  }
	  double d_delpooz,d_delpooy,d_delpoox,d_delmooz,d_delmooy,d_delmoox;
	  double d_delopoz,d_delopoy,d_delopox,d_delomoz,d_delomoy,d_delomox;
	  double d_deloopz,d_deloopy,d_deloopx,d_deloomz,d_deloomy,d_deloomx;	
	  d_delpooz=dd_out[Indix(h+1,i,j)]-dd_out[Indix(h,i,j)];
	  d_delmooz=dd_out[Indix(h,i,j)]-dd_out[Indix(h-1,i,j)];
	  d_delpooy=d_Uppo-d_Upmo;
	  d_delmooy=d_Umpo-d_Ummo;
	  d_delpoox=d_Upop-d_Upom;
	  d_delmoox=d_Umop-d_Umom;
	  d_delopoz=d_Uppo-d_Umpo;
	  d_delomoz=d_Upmo-d_Ummo;
	  d_delopoy=dd_out[Indix(h,i+1,j)]-dd_out[Indix(h,i,j)];
	  d_delomoy=dd_out[Indix(h,i,j)]-dd_out[Indix(h,i-1,j)];
	  d_delopox=d_Uopp-d_Uopm;
	  d_delomox=d_Uomp-d_Uomm;
	  d_deloopz=d_Upop-d_Umop;
	  d_deloomz=d_Upom-d_Umom;
	  d_deloopy=d_Uopp-d_Uomp;
	  d_deloomy=d_Uopm-d_Uomm;
	  d_deloopx=dd_out[Indix(h,i,j+1)]-dd_out[Indix(h,i,j)];
	  d_deloomx=dd_out[Indix(h,i,j)]-dd_out[Indix(h,i,j-1)];
	  double d_Qbarkarisum=0.0;
	  double d_Qbarkari=d_delpooz*d_delpooz+d_delpooy*d_delpooy+d_delpoox*d_delpoox;
	  if(d_Qbarkari<0.0){
	    d_Qbarkari=0.0;
	  }
	  d_Qbarkarisum+=sqrt(d_Qbarkari);
	  
	  d_Qbarkari=d_delmooz*d_delmooz+d_delmooy*d_delmooy+d_delmoox*d_delmoox;
	  if(d_Qbarkari<0.0){
	    d_Qbarkari=0.0;
	  }
	  d_Qbarkarisum+=sqrt(d_Qbarkari);
	  
	  d_Qbarkari=d_delopoz*d_delopoz+d_delopoy*d_delopoy+d_delopox*d_delopox;
	  if(d_Qbarkari<0.0){
	    d_Qbarkari=0.0;
	  }
	  d_Qbarkarisum+=sqrt(d_Qbarkari);
	  
	  d_Qbarkari=d_delomoz*d_delomoz+d_delomoy*d_delomoy+d_delomox*d_delomox;
	  if(d_Qbarkari<0.0){
	    d_Qbarkari=0.0;
	  }
	  d_Qbarkarisum+=sqrt(d_Qbarkari);
	  
	  d_Qbarkari=d_deloopz*d_deloopz+d_deloopy*d_deloopy+d_deloopx*d_deloopx;
	  if(d_Qbarkari<0.0){
	    d_Qbarkari=0.0;
	  }
	  d_Qbarkarisum+=sqrt(d_Qbarkari);
	  
	  d_Qbarkari=d_deloomz*d_deloomz+d_deloomy*d_deloomy+d_deloomx*d_deloomx;
	  if(d_Qbarkari<0.0){
	    d_Qbarkari=0.0;
	  }
	  d_Qbarkarisum+=sqrt(d_Qbarkari);
	  d_Qbarkarisum/=6.0;
	  dd_g[Indix(h,i,j)]=1.0/(1.0+LARGE_K*(d_Qbarkarisum*d_Qbarkarisum));
	}
	if(h==0||h==GRIDSIZE_Z-1||i==0||i==GRIDSIZE_Y-1||j==0||j==GRIDSIZE_X-1){
	  
	}
      }
    }
  }
}


void CalcRight_SF(LL *Lz){
  int x,y,z,idx;
  ll_init(Lz); //begining of list;
  while(Lz->curr != NULL){ 
    x = Lz->curr->x;
    y = Lz->curr->y; z = Lz->curr->z; idx = Lz->curr->idx;
    int h=z;
    int i=y;
    int j=x;
    Lz->curr->d_Right=d_phi[Indix(h,i,j)];
    if(h<1||h>GRIDSIZE_Z-3||i<1||i>GRIDSIZE_Y-3||j<1||j>GRIDSIZE_X-3){
    }
    else{
      double d_Grax=(dd_g[Indix(h,i,j+1)]-dd_g[Indix(h,i,j-1)])/2.0;
      double d_Gray=(dd_g[Indix(h,i+1,j)]-dd_g[Indix(h,i-1,j)])/2.0;
      double d_Graz=(dd_g[Indix(h+1,i,j)]-dd_g[Indix(h-1,i,j)])/2.0;
      Lz->curr->d_Right+=-wa*DT*max(-d_Grax,0.0)*(d_phi[Indix(h,i,j)]-d_phi[Indix(h,i,j-1)]);
      Lz->curr->d_Right+=-wa*DT*min(-d_Grax,0.0)*(d_phi[Indix(h,i,j+1)]-d_phi[Indix(h,i,j)]);	
      Lz->curr->d_Right+=-wa*DT*max(-d_Gray,0.0)*(d_phi[Indix(h,i,j)]-d_phi[Indix(h,i-1,j)]);	
      Lz->curr->d_Right+=-wa*DT*min(-d_Gray,0.0)*(d_phi[Indix(h,i+1,j)]-d_phi[Indix(h,i,j)]);
      Lz->curr->d_Right+=-wa*DT*max(-d_Graz,0.0)*(d_phi[Indix(h,i,j)]-d_phi[Indix(h-1,i,j)]);
      Lz->curr->d_Right+=-wa*DT*min(-d_Graz,0.0)*(d_phi[Indix(h+1,i,j)]-d_phi[Indix(h,i,j)]);
    }
    ll_step(Lz);
  }
}

void CalcQs_SF(LL *Lz){
  int x,y,z,idx;
  ll_init(Lz); //begining of list;
  while(Lz->curr != NULL){
    x = Lz->curr->x;
    y = Lz->curr->y;
    z = Lz->curr->z;
    idx = Lz->curr->idx;
    int h=z;
    int i=y;
    int j=x;
    double d_Uppo, d_Upmo, d_Umpo, d_Ummo;
    double d_Upop, d_Upom, d_Umop, d_Umom;
    double d_Uopp, d_Uopm, d_Uomp, d_Uomm;
    if(h>0&&h<GRIDSIZE_Z-1&&i>0&&i<GRIDSIZE_Y-1){
      d_Uppo=(d_phi[Indix(h,i,j)]+d_phi[Indix(h,i+1,j)]+d_phi[Indix(h+1,i,j)]+d_phi[Indix(h+1,i+1,j)])/4.0;
      d_Upmo=(d_phi[Indix(h,i,j)]+d_phi[Indix(h,i-1,j)]+d_phi[Indix(h+1,i,j)]+d_phi[Indix(h+1,i-1,j)])/4.0;
      d_Umpo=(d_phi[Indix(h,i,j)]+d_phi[Indix(h,i+1,j)]+d_phi[Indix(h-1,i,j)]+d_phi[Indix(h-1,i+1,j)])/4.0;
      d_Ummo=(d_phi[Indix(h,i,j)]+d_phi[Indix(h,i-1,j)]+d_phi[Indix(h-1,i,j)]+d_phi[Indix(h-1,i-1,j)])/4.0;
    }
    else{
      d_Uppo=(d_phi[Indix(h,i,j)]);
      d_Upmo=(d_phi[Indix(h,i,j)]);
      d_Umpo=(d_phi[Indix(h,i,j)]);
      d_Ummo=(d_phi[Indix(h,i,j)]);
    }
    if(h>0&&h<GRIDSIZE_Z-1&&j>0&&j<GRIDSIZE_X-1){
      d_Upop=(d_phi[Indix(h,i,j)]+d_phi[Indix(h,i,j+1)]+d_phi[Indix(h+1,i,j)]+d_phi[Indix(h+1,i,j+1)])/4.0;
      d_Upom=(d_phi[Indix(h,i,j)]+d_phi[Indix(h,i,j-1)]+d_phi[Indix(h+1,i,j)]+d_phi[Indix(h+1,i,j-1)])/4.0;
      d_Umop=(d_phi[Indix(h,i,j)]+d_phi[Indix(h,i,j+1)]+d_phi[Indix(h-1,i,j)]+d_phi[Indix(h-1,i,j+1)])/4.0;
      d_Umom=(d_phi[Indix(h,i,j)]+d_phi[Indix(h,i,j-1)]+d_phi[Indix(h-1,i,j)]+d_phi[Indix(h-1,i,j-1)])/4.0;
    }
    else{
      d_Upop=(d_phi[Indix(h,i,j)]);
      d_Upom=(d_phi[Indix(h,i,j)]);
      d_Umop=(d_phi[Indix(h,i,j)]);
      d_Umom=(d_phi[Indix(h,i,j)]);
    }
    if(i>0&&i<GRIDSIZE_Y-1&&j>0&&j<GRIDSIZE_X-1){
      d_Uopp=(d_phi[Indix(h,i,j)]+d_phi[Indix(h,i,j+1)]+d_phi[Indix(h,i+1,j)]+d_phi[Indix(h,i+1,j+1)])/4.0;
      d_Uopm=(d_phi[Indix(h,i,j)]+d_phi[Indix(h,i,j-1)]+d_phi[Indix(h,i+1,j)]+d_phi[Indix(h,i+1,j-1)])/4.0;
      d_Uomp=(d_phi[Indix(h,i,j)]+d_phi[Indix(h,i,j+1)]+d_phi[Indix(h,i-1,j)]+d_phi[Indix(h,i-1,j+1)])/4.0;
      d_Uomm=(d_phi[Indix(h,i,j)]+d_phi[Indix(h,i,j-1)]+d_phi[Indix(h,i-1,j)]+d_phi[Indix(h,i-1,j-1)])/4.0;
    }
    else{
      d_Uopp=(d_phi[Indix(h,i,j)]);
      d_Uopm=(d_phi[Indix(h,i,j)]);
      d_Uomp=(d_phi[Indix(h,i,j)]);
      d_Uomm=(d_phi[Indix(h,i,j)]);
    }
    double d_delpooz,d_delpooy,d_delpoox,d_delmooz,d_delmooy,d_delmoox;
    double d_delopoz,d_delopoy,d_delopox,d_delomoz,d_delomoy,d_delomox;
    double d_deloopz,d_deloopy,d_deloopx,d_deloomz,d_deloomy,d_deloomx;	
    d_delpooz=d_phi[Indix(h+1,i,j)]-d_phi[Indix(h,i,j)];
    d_delmooz=d_phi[Indix(h,i,j)]-d_phi[Indix(h-1,i,j)];
    d_delpooy=d_Uppo-d_Upmo;
    d_delmooy=d_Umpo-d_Ummo;
    d_delpoox=d_Upop-d_Upom;
    d_delmoox=d_Umop-d_Umom;
    d_delopoz=d_Uppo-d_Umpo;
    d_delomoz=d_Upmo-d_Ummo;
    d_delopoy=d_phi[Indix(h,i+1,j)]-d_phi[Indix(h,i,j)];
    d_delomoy=d_phi[Indix(h,i,j)]-d_phi[Indix(h,i-1,j)];
    d_delopox=d_Uopp-d_Uopm;
    d_delomox=d_Uomp-d_Uomm;
    d_deloopz=d_Upop-d_Umop;
    d_deloomz=d_Upom-d_Umom;
    d_deloopy=d_Uopp-d_Uomp;
    d_deloomy=d_Uopm-d_Uomm;
    d_deloopx=d_phi[Indix(h,i,j+1)]-d_phi[Indix(h,i,j)];
    d_deloomx=d_phi[Indix(h,i,j)]-d_phi[Indix(h,i,j-1)];
    double d_Qpoo,d_Qmoo,d_Qopo,d_Qomo,d_Qoop,d_Qoom;
    d_Qpoo=sqrt(EPSILON*EPSILON+d_delpooz*d_delpooz+d_delpooy*d_delpooy+d_delpoox*d_delpoox);
    d_Qmoo=sqrt(EPSILON*EPSILON+d_delmooz*d_delmooz+d_delmooy*d_delmooy+d_delmoox*d_delmoox);
    d_Qopo=sqrt(EPSILON*EPSILON+d_delopoz*d_delopoz+d_delopoy*d_delopoy+d_delopox*d_delopox);
    d_Qomo=sqrt(EPSILON*EPSILON+d_delomoz*d_delomoz+d_delomoy*d_delomoy+d_delomox*d_delomox);
    d_Qoop=sqrt(EPSILON*EPSILON+d_deloopz*d_deloopz+d_deloopy*d_deloopy+d_deloopx*d_deloopx);
    d_Qoom=sqrt(EPSILON*EPSILON+d_deloomz*d_deloomz+d_deloomy*d_deloomy+d_deloomx*d_deloomx);
    double d_Qbarkarisum=0.0;
    double d_Qbarkari=d_delpooz*d_delpooz+d_delpooy*d_delpooy+d_delpoox*d_delpoox;
    if(d_Qbarkari<0.0){
      d_Qbarkari=0.0;
    }
    d_Qbarkarisum+=sqrt(d_Qbarkari);
    
    d_Qbarkari=d_delmooz*d_delmooz+d_delmooy*d_delmooy+d_delmoox*d_delmoox;
    if(d_Qbarkari<0.0){
      d_Qbarkari=0.0;
    }
    d_Qbarkarisum+=sqrt(d_Qbarkari);
    
    d_Qbarkari=d_delopoz*d_delopoz+d_delopoy*d_delopoy+d_delopox*d_delopox;
    if(d_Qbarkari<0.0){
      d_Qbarkari=0.0;
    }
    d_Qbarkarisum+=sqrt(d_Qbarkari);
    
    d_Qbarkari=d_delomoz*d_delomoz+d_delomoy*d_delomoy+d_delomox*d_delomox;
    if(d_Qbarkari<0.0){
      d_Qbarkari=0.0;
    }
    d_Qbarkarisum+=sqrt(d_Qbarkari);
    
    d_Qbarkari=d_deloopz*d_deloopz+d_deloopy*d_deloopy+d_deloopx*d_deloopx;
    if(d_Qbarkari<0.0){
      d_Qbarkari=0.0;
    }
    d_Qbarkarisum+=sqrt(d_Qbarkari);
    
    d_Qbarkari=d_deloomz*d_deloomz+d_deloomy*d_deloomy+d_deloomx*d_deloomx;
    if(d_Qbarkari<0.0){
      d_Qbarkari=0.0;
    }
    d_Qbarkarisum+=sqrt(d_Qbarkari);
    d_Qbarkarisum/=6.0;
    double d_Qbar;
    d_Qbar=sqrt(EPSILON*EPSILON+d_Qbarkarisum*d_Qbarkarisum);
    Lz->curr->d_czp =    -wd*dd_g[Indix(h,i,j)]*DT*(d_Qbar/d_Qpoo);
    Lz->curr->d_czm =    -wd*dd_g[Indix(h,i,j)]*DT*(d_Qbar/d_Qmoo);
    Lz->curr->d_cyp =    -wd*dd_g[Indix(h,i,j)]*DT*(d_Qbar/d_Qopo);
    Lz->curr->d_cym =    -wd*dd_g[Indix(h,i,j)]*DT*(d_Qbar/d_Qomo);
    Lz->curr->d_cxp =    -wd*dd_g[Indix(h,i,j)]*DT*(d_Qbar/d_Qoop);
    Lz->curr->d_cxm =    -wd*dd_g[Indix(h,i,j)]*DT*(d_Qbar/d_Qoom);
    Lz->curr->d_cce = 1.0+wd*dd_g[Indix(h,i,j)]*DT*(d_Qbar/d_Qpoo+d_Qbar/d_Qmoo+d_Qbar/d_Qopo+d_Qbar/d_Qomo+d_Qbar/d_Qoop+d_Qbar/d_Qoom);  
    ll_step(Lz);
  }
}


	

void SolveUsingSOR_SF(int my_rank,LL *Lz){
  
  int enu;
  int i,j,k;
  double delp,pdel;
  double omg,t1;
  int flagg = 1;
  int i1,j1,k1;
  int ie,je,ke;
  int ii,jj,kk;
  double pdelmax;
  // I have to solve from 0toGRIDSIZE-1
  // in this program I set 0 and GRIDSIZE-1 Dirichlet boundary
  i1 = GRIDSIZE_Z-1-1;
  j1 = GRIDSIZE_Y-1-1;
  k1 = GRIDSIZE_X-1-1;
  for(enu=1; enu<=nend;enu++){
    delp = 0.0;
    pdelmax =-1.0;    
    if(flagg==1){
      int x,y,z,idx;
      ll_init(Lz); //begining of list;
      while(Lz->curr != NULL){
	x = Lz->curr->x;
	y = Lz->curr->y;
	z = Lz->curr->z;
	idx = Lz->curr->idx;
	int i=z;
	int j=y;
	int k=x;
	if(i<1||i>GRIDSIZE_Z-3||j<1||j>GRIDSIZE_Y-3||k<1||k>GRIDSIZE_X-3)continue;
	omg = om/Lz->curr->d_cce;
	t1= (Lz->curr->d_cce*d_phi[Indix(i,j,k)]
	     +Lz->curr->d_czp*d_phi[Indix(i+1,j,k)]
	     +Lz->curr->d_czm*d_phi[Indix(i-1,j,k)]
	     +Lz->curr->d_cyp*d_phi[Indix(i,j+1,k)]
	     +Lz->curr->d_cym*d_phi[Indix(i,j-1,k)]
	     +Lz->curr->d_cxp*d_phi[Indix(i,j,k+1)]
	     +Lz->curr->d_cxm*d_phi[Indix(i,j,k-1)]-Lz->curr->d_Right);
	pdel = -omg*t1;
	if(abs(t1)>pdelmax){
	  pdelmax=t1;
	}
	delp = delp + pdel*pdel;
	d_phi[Indix(i,j,k)] =d_phi[Indix(i,j,k)] +pdel;
	ll_step(Lz);
      }
      delp = sqrt(delp/(double)((GRIDSIZE_X-2)*(GRIDSIZE_Y-2)*(GRIDSIZE_Z-2)));
      //rintf("res %f\n",delp);
      if(pdelmax<convp){
	//printf("fin\n");
	flagg=0;
      }
  
    }
  }
  if(flagg==1){
    printf("notfin %d\n",my_rank);
  }
}

int main(int argc, char* argv[]){
  int my_rank;
  int n_proc;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
  long    dimx, dimy, dimz;
  //InitializeCircleMap();
  //initialize_list();
  //dims[2] = 1;
  //dims[1] = 1;
  //switch(numdims){
  //case 3: dimz = (int)mdims[2]; dims[2] = dimz;
  //case 2: dimx = (int)mdims[1]; dims[1] = dimx;
  //case 1: dimy = (int)mdims[0]; dims[0] = dimy;
  //}
  dims[3] = dims[0]*dims[1];
  dims[4] = dims[0]*dims[1]*dims[2];
  int saki;
  int n1,n2,n3,n4;
  int img_num=0;
  SetFilename();
  int reset = 1;
  double Fs = 0.0;
  FILE *fp;
  
  char infile[50];
  FILE *fq;
  //todo set of mask
  MPI_Barrier(MPI_COMM_WORLD);
  for(int kuri=0;kuri<1;kuri++){
    printf("unko\n");
    fp=fopen("/g/hiiragi/Users/Ritsuya/corrected190114/nuclei/t580-nuclei","r");
    
    double usum, vsum, dthresh, lambda, rad;
    int    iter,countdown,display;
    //ScanData(fp, my_rank);//readthenucleiposition
    ReadTIFFile2forSF(my_rank); 
    dims[2] = 1;
    dims[1] = 1;
    switch(numdims){
    case 3: dimz = (int)mdims[2]; dims[2] = dimz;
    case 2: dimx = (int)mdims[1]; dims[1] = dimx;
    case 1: dimy = (int)mdims[0]; dims[0] = dimy;
    }
    dims[3] = dims[0]*dims[1];
    dims[4] = dims[0]*dims[1]*dims[2];
    for(int i=79;i<119;i++){
      for(int j=43;j<109;j++){
	d_mask[Indix(0,j,i)]=1.0;
      }
    }
    lambda=0.1;
    iter=1000;
    rad=15;
    display=1;
    dthresh =500;
    printf("%ld %ld %ld %ld %ld\n",dims[0],dims[1],dims[2],dims[3],dims[4]);
    printf("dth%f lam %f rad%f\n",dthresh,lambda,rad);
    printf("%d %d %d\n",GRIDSIZE_Y,GRIDSIZE_X,GRIDSIZE_Z);
    int yy,xx;
    for (yy=0;yy<(int)100; yy++){
      for (xx=0;xx<(int)100;xx++){
	if(xx>50&&xx<60&&yy>50&&yy<60){
	  printf("img_info %d %d %f\n",xx,yy,dd_out[Indix(0,yy,xx)]);
	}
      }
    }
    LL *Lz, *Ln1, *Ln2, *Lp1, *Lp2;
    LL *Lin2out, *Lout2in;
    Lz  = ll_create();
    Ln1 = ll_create();
    Ln2 = ll_create();
    Lp1 = ll_create();
    Lp2 = ll_create();
    Lin2out = ll_create();
    Lout2in = ll_create();
    
    //Calcg();
    ls_mask2phi3c(d_mask,d_phi,d_label,dims,Lz,Ln1,Ln2,Lp1,Lp2);
   
    
    lrbac_chanvese(dd_out,d_phi,d_label,dims,Lz,Ln1,Lp1,Ln2,Lp2,Lin2out,Lout2in,iter,rad,lambda,display);
    
    WriteFrontPixelsforSF(my_rank,img_num,Lz);
    ll_destroy(Lz);
    ll_destroy(Ln1);
    ll_destroy(Ln2);
    ll_destroy(Lp1);
    ll_destroy(Lp2);
    ll_destroy(Lin2out);
    ll_destroy(Lout2in);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}

