#include "mem_seg_levelset030414sf.h"
bool startsWith(const char *pre, const char *str)
{
    size_t lenpre = strlen(pre),
           lenstr = strlen(str);
    return lenstr < lenpre ? false : strncmp(pre, str, lenpre) == 0;
}
void Writeparameters(int myrank,FILE* filp){
  if(myrank==0){
    /*
    int curint;
    double curdouble;
    curint= TOTFILENAME;
    fprintf(filp,"TOTFILENAME %d\n",curint);
    curint= MAX_CHAR_NUM;
    fprintf(filp,"MAX_CHAR_NUM %d\n",curint);
    curint=OTHER_FILE_NUM;
    fprintf(filp,"OTHER_FILE_NUM %d\n",curint);
    curint=WBAND;
    fprintf(filp,"WBAND %d\n",curint);
    curint=WRESET;
    fprintf(filp,"WRESET %d\n",curint);
    curint=numdims;
    fprintf(filp,"numdims %d\n",curint);
    curint=totalpixels;
    fprintf(filp,"totalpixels %d\n",curint);
  
#define DIMX dims[1]
#define DIMY dims[0]
#define DIMZ dims[2]
#define DIMXY dims[3]
#define NUMEL dims[4]

#define OFFX dims[0]
#define OFFY 1
#define OFFZ dims[3]
  

//mainly for solving equation
    curint= ILIM;
    fprintf(filp,"ILIM %d\n",curint);
    
    curint=OUTPUTINVAL;
    fprintf(filp,"OUTPUTINVAL %d\n",curint);
    curdouble=EPSILON;
    fprintf(filp,"EPSILON %f.10\n",curdouble);
    curdouble= DT;
    fprintf(filp,"DT %f\n",curdouble);
    curint=KURILIM;
    fprintf(filp,"KURILIM %d\n",curint);
    curdouble=EPS_ITEL;
    fprintf(filp,"EPS_ITEL %f\n",curdouble);
    curdouble=wa;
    fprintf(filp,"wa %f\n",curdouble);
    curdouble=wd;
    fprintf(filp,"wd %f\n",curdouble);
    curdouble=DELTA;
    fprintf(filp,"DELTA %f\n",curdouble);
    curdouble=LARGE_K;
    fprintf(filp,"LARGE_K %f\n",curdouble);
    curint=MAXCELLNUM;
    fprintf(filp,"MAXCELLNUM %d\n",curint);
    
//SOR
    curdouble=om;
    fprintf(filp,"om %f\n",curdouble);
    curint=nend;
    fprintf(filp,"nend %d\n",curint);
    
    curdouble=convp;
    fprintf(filp,"convp %f\n",curdouble);
//not used
//for InnerCell
    curdouble=DWBAND;
    fprintf(filp,"DWBAND %f\n",curdouble);
    curint=LIM_RAD;
    fprintf(filp,"LIM_RAD %d\n",curint);
    curdouble=z_factor;
    fprintf(filp,"z_factor %f\n",curdouble);
    curdouble=fin_valu;
    fprintf(filp,"fin_valu %f\n",curdouble);
    curint=INITIALTIME;
    fprintf(filp,"INITIALTIME %d\n",curint);
    curint= ENDTIME;
    fprintf(filp,"ENDTIME %d\n",curint);
    */
  }
  if(myrank==0){ 
    FILE * fp=fopen(headername,"r");
    char row[256];
    char gig[256]="#define";
    while(fgets(row,sizeof(row),fp)!=NULL){
      if(startsWith(gig,row)){
	fprintf(filp,"%s\n",row);
      }
    }
  } 
}
void ScanDataALL(int  my_rank){
 int i;
  double minekus = 100.0;
  double maxzee = -100.0;
  int aka, ao, mido, kuro, kii;
  int tate_count, VY_NUM, VT_NUM;
  int num,ok,pred,suc1,suc2,size;
  int gx,gy;
  double gz;
  char name[100];
  int weight;
  double x,y, z1,z2;
  char s2[100];
  char s[100];
  FILE *fp;
  NumDetectedNuc=0;
  for(int current_t=INITIALTIME;current_t<=ENDTIME;current_t++){
    char nucfile[MAX_CHAR_NUM];
    nucfile[0] = '\0';
    strcpy(nucfile,dnucdirname);
    char localfilename[MAX_CHAR_NUM];
    localfilename[0] = '\0';
    sprintf(localfilename,"t%03d-nuclei",current_t);
    strcat(nucfile,localfilename);
    fp=fopen(nucfile,"r");
    if(my_rank == 0){
      while(fscanf(fp,"%d, %d, %d, %d, %d, %d, %d, %lf, %d,%[^,], %d, %d, %d, %d,%[^,],",&num,&ok,&pred,&suc1,&suc2,&gx,&gy,&gz,&size,name,&weight,&ao,&mido,&kuro,s2)!=EOF){
	if(strlen(name)>2){
	  i_NucT[NumDetectedNuc]=current_t-INITIALTIME;
	  i_NucC[NumDetectedNuc]=num;
	  double d_gx=gx*MAG_XY;
	  double d_gy=gy*MAG_XY;
	  i_NucX[NumDetectedNuc]=(int)d_gx;
	  i_NucY[NumDetectedNuc]=(int)d_gy;
	  double dz =(gz+0.5)*MAG_Z;
	  int iz = (int)dz;
	  i_NucZ[NumDetectedNuc]=iz;
	  //fprintf(FileLog,"%d %d and%s and %s a\n",current_t,num,name,s2);	
	  NumDetectedNuc++;
	}
      }
    }
    fclose(fp);
  }
  if(my_rank==0){
    fprintf(FileLog,"NumDetectedNuc and xyz%d %d %d %d %d\n",NumDetectedNuc,i_NucT[NumDetectedNuc-1],i_NucX[NumDetectedNuc-1],i_NucY[NumDetectedNuc-1],i_NucZ[NumDetectedNuc-1]);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&NumDetectedNuc,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(i_NucC,MAXCELLNUM,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(i_NucT,MAXCELLNUM,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(i_NucX,MAXCELLNUM,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(i_NucY,MAXCELLNUM,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(i_NucZ,MAXCELLNUM,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

}
void ScanData(FILE *fp, int  my_rank){
  int i;
  double minekus = 100.0;
  double maxzee = -100.0;
  int aka, ao, mido, kuro, kii;
  int tate_count, VY_NUM, VT_NUM;
  int num,ok,pred,suc1,suc2,size;
  int gx,gy;
  double gz;
  char name[100];
  int weight;
  double x,y, z1,z2;
  char s2[100];
  char s[100];
  NumDetectedNuc=0;
  if(my_rank == 0){
    while(fscanf(fp,"%d, %d, %d, %d, %d, %d, %d, %lf, %d,%[^,], %d, %d, %d, %d,%[^,],",&num,&ok,&pred,&suc1,&suc2,&gx,&gy,&gz,&size,name,&weight,&ao,&mido,&kuro,s2)!=EOF){
      if(strlen(name)>1){
	//double d_gx=(gx*5.0/4.0-30.0)*1.1/2.0;
	//double d_gy=(gy*5.0/4.0-54.0)*1.1/2.0;
	double d_gx=gx*MAG_XY;
	double d_gy=gy*MAG_XY;
	i_NucX[NumDetectedNuc]=(int)d_gx;
	i_NucY[NumDetectedNuc]=(int)d_gy;
	double dz =(gz+0.5)*MAG_Z;
	int iz = (int)dz;
	//iz=iz-20;
	i_NucZ[NumDetectedNuc]=iz;
	printf("%d and%s and %s a\n",num,name,s2);	
	NumDetectedNuc++;
      }
    }
  }
  fclose(fp);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(i_NucX,MAXCELLNUM,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(i_NucY,MAXCELLNUM,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(i_NucZ,MAXCELLNUM,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

}

void WriteFrontPixelsforSF(int my_rank,int img_num,LL *Lz,int step_num){
  char infile[MAX_CHAR_NUM];
  infile[0] = '\0';
  strcpy(infile,diroutname);
  char infile2[MAX_CHAR_NUM];
  sprintf(infile2,"out_t%04d_cell%03dFR_Step%04d.txt",img_num+1,my_rank,step_num);
  strcat(infile, infile2);
  if(step_num%5==0){
    FILE* fq;
    fq=fopen(infile,"w");
    ll_init(Lz);
    int count = 0;
    if(Lz!=NULL){
      while(Lz->curr != NULL){
	fprintf(fq,"%ld\t%ld\t%ld\n", Lz->curr->z,Lz->curr->y,Lz->curr->x);
	ll_step(Lz);
	// printf("bullsshit\n");
	count++;
      }
    }
    fclose(fq);
  }
  ll_init(Lz);
  int countn = 0;
  if(Lz!=NULL){
    while(Lz->curr != NULL){
      ll_step(Lz);
      // printf("bullsshit\n");
      countn++;
    }
  }
  new_count = countn;
  
}
