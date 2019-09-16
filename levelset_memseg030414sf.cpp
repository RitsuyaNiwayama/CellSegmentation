#include "mem_seg_levelset030414sf.h"
#include <mpi.h>

DIR* dir;
struct dirent* dp;
DIR* dir_g;
struct dirent* dp_g;
int count_file = 0;
int count_file_g = 0;
//char ddirname[MAX_CHAR_NUM] =    "/scratch/niwayama/dst15/";
//char dnucdirname[MAX_CHAR_NUM] = "/scratch/niwayama/nuclei/";

//char diroutname[MAX_CHAR_NUM] =  "/home/niwayama/memsegmentation/140212sparsefield/140219sfm/dst5/";
//char diroutname[MAX_CHAR_NUM] =  "/g/hiiragi/Users/Ritsuya/membraneMPI/segment1_180314/";
char diroutname[MAX_CHAR_NUM];
char dirlogname[MAX_CHAR_NUM];
char headername[MAX_CHAR_NUM] = "/g/hiiragi/Users/Ritsuya/140403sfm140711modsacked/mem_seg_levelset030414sf.h";
//char ddirname[MAX_CHAR_NUM] =    "/g/hiiragi/Users/Ritsuya/260514track/c4dst4mod/";
//char dnucdirname[MAX_CHAR_NUM] = "/g/hiiragi/Users/Ritsuya/260514track/nuclei/";
//char globalname[MAX_CHAR_NUM] = "/g/hiiragi/Users/Ritsuya/260514track/c5dst4alpha/";
char ddirname[MAX_CHAR_NUM] =    "/g/hiiragi/Users/Ritsuya/outv9/mdst5/";
char dnucdirname[MAX_CHAR_NUM] = "/g/hiiragi/Users/Ritsuya/outv9/ndst17/";
char globalname[MAX_CHAR_NUM] = "/g/hiiragi/Users/Ritsuya/outv9/mdst7/";
//char ddirname[MAX_CHAR_NUM] =    "/g/hiiragi/Users/Ritsuya/spinningdisk/processing/2016-02-22/numi20/mdst5/";
//char dnucdirname[MAX_CHAR_NUM] = "/g/hiiragi/Users/Ritsuya/spinningdisk/processing/2016-02-22/numi20/mdst8/";
//char globalname[MAX_CHAR_NUM] =  "/g/hiiragi/Users/Ritsuya/spinningdisk/processing/2016-02-22/numi20/mdst7/";
char FILENAME[TOTFILENAME][MAX_CHAR_NUM];
char FILENAME_G[TOTFILENAME][MAX_CHAR_NUM];

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
int i_NucT[MAXCELLNUM];
int i_NucX[MAXCELLNUM];
int i_NucY[MAXCELLNUM];
int i_NucZ[MAXCELLNUM];
int i_NucC[MAXCELLNUM];
int GRIDSIZE_Z;
int GRIDSIZE_Y;
int GRIDSIZE_X;
long    dims[5];
long   mdims[3];
TOOL_t Tools;

double new_count;
double old_count;
FILE *FileLog;
//140711added
time_t g_time_t1;
time_t g_time_t2;

#define USAGE(argv)\
{\
  fprintf(stderr,"Usage %s imagefile \n",argv[0]);\
  exit(-1);\
}
int Indix(int ii, int kk, int jj){
  int idx=(int)(ii*DIMXY+jj*DIMY+kk);
  return idx;
}


double ReInitialization(LL *Lz,double setd){
  int x,y,z,idx;
  int ret_valu=0;
  ll_init(Lz); //begining of list;
  double existn=0.0;
  while(Lz->curr != NULL){ 
    x = Lz->curr->x;
    y = Lz->curr->y; z = Lz->curr->z; idx = Lz->curr->idx;
    int h=z;
    int i=y;
    int j=x;
    d_phi[Indix(h,i,j)]=setd;
    existn=1.0;
    ll_step(Lz);
  }
  return existn;
}

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
double Check_Variation(LL *Lz){
  int x,y,z,idx;
  int ret_valu=0;
  ll_init(Lz); //begining of list;
  double max=-100000000000.0;
  while(Lz->curr != NULL){ 
    x = Lz->curr->x;
    y = Lz->curr->y; z = Lz->curr->z; idx = Lz->curr->idx;
    int h=z;
    int i=y;
    int j=x;
    double sa =fabs(Lz->curr->d_oldphi-d_phi[Indix(h,i,j)]);
    if(sa>max){
      max=sa;
    }
    ll_step(Lz);
  }
  return max;
}


int TestifbeingOut(LL *Lz){
  int x,y,z,idx;
  int flagging=1;
  ll_init(Lz); //begining of list;
  while(Lz->curr != NULL){ 
    x = Lz->curr->x;
    y = Lz->curr->y; z = Lz->curr->z; idx = Lz->curr->idx;
    int h=z;
    int i=y;
    int j=x;
    Lz->curr->d_Right=d_phi[Indix(h,i,j)];
    Lz->curr->d_oldphi=d_phi[Indix(h,i,j)];
    if(h<2||h>GRIDSIZE_Z-3||i<2||i>GRIDSIZE_Y-3||j<2||j>GRIDSIZE_X-3){
      flagging=0;
    }
    ll_step(Lz);
  }
  return flagging;
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
    Lz->curr->d_oldphi=d_phi[Indix(h,i,j)];
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
    double d_Dpoo,d_Dmoo,d_Dopo,d_Domo,d_Doop,d_Doom;//zyx
    d_Dpoo =min(d_phi[Indix(h+1,i,j)]-d_phi[Indix(h,i,j)],0);
    d_Dpoo*=min(d_phi[Indix(h+1,i,j)]-d_phi[Indix(h,i,j)],0);
    d_Dmoo =min(d_phi[Indix(h-1,i,j)]-d_phi[Indix(h,i,j)],0);
    d_Dmoo*=min(d_phi[Indix(h-1,i,j)]-d_phi[Indix(h,i,j)],0);
    d_Dopo =min(d_phi[Indix(h,i+1,j)]-d_phi[Indix(h,i,j)],0);
    d_Dopo*=min(d_phi[Indix(h,i+1,j)]-d_phi[Indix(h,i,j)],0);
    d_Domo =min(d_phi[Indix(h,i-1,j)]-d_phi[Indix(h,i,j)],0);
    d_Domo*=min(d_phi[Indix(h,i-1,j)]-d_phi[Indix(h,i,j)],0);
    d_Doop =min(d_phi[Indix(h,i,j+1)]-d_phi[Indix(h,i,j)],0);
    d_Doop*=min(d_phi[Indix(h,i,j+1)]-d_phi[Indix(h,i,j)],0);
    d_Doom =min(d_phi[Indix(h,i,j-1)]-d_phi[Indix(h,i,j)],0);
    d_Doom*=min(d_phi[Indix(h,i,j-1)]-d_phi[Indix(h,i,j)],0);
    double d_Mioo = max(d_Dpoo,d_Dmoo);
    double d_Moio = max(d_Dopo,d_Domo);
    double d_Mooi = max(d_Doop,d_Doom);
    double ssq = d_Mioo+d_Moio+d_Mooi;
    if(ssq<0.0){
      ssq = 0.0;
    }
    else{
      ssq = sqrt(ssq);
    }
    Lz->curr->d_Right+=DELTA*DT*dd_g[Indix(h,i,j)]*ssq;
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
      g_time_t2=time(NULL);
      int dif_t=(int)(g_time_t2-g_time_t1);
      int d_TIMELIMITITR=TIMELIMITITR;
      if(dif_t>TIMELIMITITR)goto CHINKO;
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
	if(i<1||i>GRIDSIZE_Z-3||j<1||j>GRIDSIZE_Y-3||k<1||k>GRIDSIZE_X-3)goto CHINKO;
	ll_step(Lz);
      }
      ll_init(Lz); //begining of list;
      while(Lz->curr != NULL){
	x = Lz->curr->x;
	y = Lz->curr->y;
	z = Lz->curr->z;
	idx = Lz->curr->idx;
	int i=z;
	int j=y;
	int k=x;
	//if(i<1||i>GRIDSIZE_Z-3||j<1||j>GRIDSIZE_Y-3||k<1||k>GRIDSIZE_X-3)goto chinko;
	//if(d_mask[Indix(i,j,k)]!=0.0)continue;
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
	//	printf("fin\n");
	flagg=0;
      }
      
    }
  }
  if(flagg==1){
    //  printf("notfin %d\n",my_rank);
  }
 CHINKO:int gggg=1;

}

int main(int argc, char* argv[]){
  int my_rank;
  int n_proc;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_proc); 
  strcpy(diroutname,argv[1]);
  strcpy(dirlogname,argv[2]);
   
  long    dimx, dimy, dimz;
  char filelogname[MAX_CHAR_NUM];
  filelogname[0] = '\0';
  strcpy(filelogname,dirlogname);
  strcat(filelogname,"log.text");
  FileLog=fopen(filelogname,"w");
  MPI_Barrier(MPI_COMM_WORLD);
  fprintf(FileLog,"first np %d myrank %d\n", n_proc,my_rank);
  if(my_rank==0){
   fprintf(FileLog,"outname %s\n",diroutname);
   fprintf(FileLog,"logname %s\n",dirlogname);
  }
  Writeparameters(my_rank,FileLog);
  
  dims[3] = dims[0]*dims[1];
  dims[4] = dims[0]*dims[1]*dims[2];
  int saki;
  int n1,n2,n3,n4;
  int img_num=0;
  SetFilename();
  SetFilename_G();
  int reset = 1;
  double Fs = 0.0;
  FILE *fp;
  char infile[50];
  FILE *fq;
  MPI_Barrier(MPI_COMM_WORLD);
  MakeTools (&Tools, 30,z_factor);
  MPI_Barrier(MPI_COMM_WORLD);
  ScanDataALL(my_rank);
  MPI_Barrier(MPI_COMM_WORLD);
 
  int success_of_ini_nuc_pos=1;
  for(int kuri=0;kuri<KURILIM;kuri++){
    
    int yobi= (kuri+my_rank+1)%n_proc;
    int current_nuc =kuri*n_proc+yobi;
    fprintf(FileLog,"np %d myrank %d yobi %d cur %d\n", n_proc,my_rank,yobi,current_nuc);
		printf("argforTIFFREAD  %d \n", i_NucT[current_nuc]+OFFTIME);
    if(current_nuc<NumDetectedNuc){
      ReadTIFFile2forSF(i_NucT[current_nuc]+OFFTIME);
		  ReadTIFFilemaskforSF(i_NucT[current_nuc]+OFFTIME);
      g_time_t1=time(NULL);
      for(int opi=0;opi<3;opi++){
					success_of_ini_nuc_pos=OptimizeNucleiPosition(current_nuc,&Tools);
      }
			img_num=i_NucT[current_nuc]+INITIALTIME-1;
			InitializeFrontPositionCellbyCellMPIforSF(2,current_nuc);
      
      LL *Lz, *Ln1, *Ln2, *Lp1, *Lp2;
      LL *Lin2out, *Lout2in;
      Lz  = ll_create();
      Ln1 = ll_create();
      Ln2 = ll_create();
      Lp1 = ll_create();
      Lp2 = ll_create();
      Lin2out = ll_create();
      Lout2in = ll_create();
  
      Calcg();
      //d_mask ok
      ls_mask2phi3c(d_mask,d_phi,d_label,dims,Lz,Ln1,Ln2,Lp1,Lp2);
      int flag_ii=1;
      old_count=0.0;
      int step_num=0;
      //
      if(success_of_ini_nuc_pos!=1){
	flag_ii=0;
      }
      for(int ii =1;ii< ILIM;ii++){
	double remaining=1.0;
	if(flag_ii==1){
	  CalcRight_SF(Lz);
	  CalcQs_SF(Lz);
	  SolveUsingSOR_SF(my_rank,Lz);
	  
	  double maxsa= Check_Variation(Lz);
	  if(maxsa>0.5){
	    fprintf(FileLog,"specify lower dt:fin %d cell %d time %d\n",current_nuc,i_NucC[current_nuc],i_NucT[current_nuc]);
	    break;
	  }
	  ls_iteration_RN(d_phi, d_label,dims, Lz, Ln1,Lp1,Ln2,Lp2, Lin2out,Lout2in);
	}
	if(ii%OUTPUTINVAL==2){
		// printf("%d \n",img_num);
		WriteFrontPixelsforSF(i_NucC[current_nuc],img_num,Lz,step_num);
	  step_num++;
	  double valu = (new_count-old_count)/new_count;
	  if(valu<fin_valu){
	    flag_ii=0;
	  }
	  if(new_count>maxpixels){
	    flag_ii=0;
	  }
	    
	  old_count=new_count;
	}
	
      }
      ll_destroy(Lz);
      ll_destroy(Ln1);
      ll_destroy(Ln2);
      ll_destroy(Lp1);
      ll_destroy(Lp2);
      ll_destroy(Lin2out);
      ll_destroy(Lout2in);
	
    }
    
  }
  fclose(FileLog);
  CleanUpTools (&Tools,30);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}

