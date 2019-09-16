#include "mem_seg_levelset030414sf.h"
int OptimizeNucleiPosition(int my_rank,TOOL_t * const tools){
  int az = i_NucZ[my_rank];
  int ay = i_NucY[my_rank];
  int ax = i_NucX[my_rank];
  //  int directions =8;
  // from 140711
  int directions=18;
  int max_go =5;
  int quant_rad = 5;
  //def and initialization
  int Xm[directions][max_go];
  int Ym[directions][max_go];
  int Zm[directions][max_go];
  double Intm[directions][max_go];
  int Index_MinDirec=0;
  int CheckD[directions];
  int MaXDevNum[directions];
  double IntMax[directions];
  double IntMin[directions];
  int IntMinIndex[directions];
  double Intorig=Inten(tools,ax,ay,az,quant_rad);
  //printf("oring %d %d %d %f\n",ax,ay,az,Intorig);
  for(int i=0;i<directions;i++){
    for(int j=0;j<max_go;j++){
      Xm[i][j]=0;
      Ym[i][j]=0;
      Zm[i][j]=0;
      Intm[i][j]=0.0;
    }
    CheckD[i]=1;
    MaXDevNum[i]=0;
    IntMax[i]= -100000000000.0;
    IntMin[i]=  100000000000.0;
    IntMinIndex[i]=0;
  }
  //direction 1
  for(int i=1;i<max_go;i++){
    Xm[0][i]=ax+i;
    Ym[0][i]=ay;
    Zm[0][i]=az;
    MaXDevNum[0]+=1;
  }
  //direction 2
  for(int i=1;i<max_go;i++){
    Xm[1][i]=ax-i;
    Ym[1][i]=ay;
    Zm[1][i]=az;
    MaXDevNum[1]+=1;
  }
  //direction 3
  for(int i=1;i<max_go;i++){
    Xm[2][i]=ax;
    Ym[2][i]=ay+i;
    Zm[2][i]=az;
    MaXDevNum[2]+=1;
  }
  //direction 4
  for(int i=1;i<max_go;i++){
    Xm[3][i]=ax;
    Ym[3][i]=ay-i;
    Zm[3][i]=az;
    MaXDevNum[3]+=1;
  }
  
  //direction 5
  for(int i=1;i<max_go;i++){
    int dx=(int)(i*1.0/sqrt(2.0));
    int dy=(int)(i*1.0/sqrt(2.0));
    
    Xm[4][i]=ax+dx;
    Ym[4][i]=ay+dy;
    Zm[4][i]=az;
    MaXDevNum[4]+=1;
  }
  //direction 6
  for(int i=1;i<max_go;i++){
    int dx=-(int)(i*1.0/sqrt(2.0));
    int dy=(int)(i*1.0/sqrt(2.0));
    
    Xm[5][i]=ax+dx;
    Ym[5][i]=ay+dy;
    Zm[5][i]=az;
    MaXDevNum[5]+=1;
  }
  //direction 7
  for(int i=1;i<max_go;i++){
    int dx=(int)(i*1.0/sqrt(2.0));
    int dy=-(int)(i*1.0/sqrt(2.0));
    
    Xm[6][i]=ax+dx;
    Ym[6][i]=ay+dy;
    Zm[6][i]=az;
    MaXDevNum[6]+=1;
  }
  //direction 8
  for(int i=1;i<max_go;i++){
    int dx=-(int)(i*1.0/sqrt(2.0));
    int dy=-(int)(i*1.0/sqrt(2.0));
    
    Xm[7][i]=ax+dx;
    Ym[7][i]=ay+dy;
    Zm[7][i]=az;
    MaXDevNum[7]+=1;
  }
  //added 140711
  //direction 9
  for(int i=1;i<max_go;i++){
    Xm[8][i]=ax;
    Ym[8][i]=ay;
    Zm[8][i]=az+i;
    MaXDevNum[8]+=1;
  }
  //direction 10
  for(int i=1;i<max_go;i++){
    Xm[9][i]=ax;
    Ym[9][i]=ay;
    Zm[9][i]=az-i;
    MaXDevNum[9]+=1;
  }
  //direction 11
  for(int i=1;i<max_go;i++){
    int dx=(int)(i*1.0/sqrt(2.0));
    int dy=(int)(i*1.0/sqrt(2.0));
    
    Xm[10][i]=ax;
    Ym[10][i]=ay+dy;
    Zm[10][i]=az+dx;
    MaXDevNum[10]+=1;
  }
  //direction 12
  for(int i=1;i<max_go;i++){
    int dx=-(int)(i*1.0/sqrt(2.0));
    int dy=(int)(i*1.0/sqrt(2.0));
    
    Xm[11][i]=ax;
    Ym[11][i]=ay+dy;
    Zm[11][i]=az+dx;
    MaXDevNum[11]+=1;
  }
  //direction 13
  for(int i=1;i<max_go;i++){
    int dx=(int)(i*1.0/sqrt(2.0));
    int dy=-(int)(i*1.0/sqrt(2.0));
    
    Xm[12][i]=ax;
    Ym[12][i]=ay+dy;
    Zm[12][i]=az+dx;
    MaXDevNum[12]+=1;
  }
  //direction 14
  for(int i=1;i<max_go;i++){
    int dx=-(int)(i*1.0/sqrt(2.0));
    int dy=-(int)(i*1.0/sqrt(2.0));
    
    Xm[13][i]=ax;
    Ym[13][i]=ay+dy;
    Zm[13][i]=az+dx;
    MaXDevNum[13]+=1;
  }
  //direction 15
  for(int i=1;i<max_go;i++){
    int dx=(int)(i*1.0/sqrt(2.0));
    int dy=(int)(i*1.0/sqrt(2.0));
    
    Xm[14][i]=ax+dy;
    Ym[14][i]=ay;
    Zm[14][i]=az+dx;
    MaXDevNum[14]+=1;
  }
  //direction 16
  for(int i=1;i<max_go;i++){
    int dx=-(int)(i*1.0/sqrt(2.0));
    int dy=(int)(i*1.0/sqrt(2.0));
    
    Xm[15][i]=ax+dy;
    Ym[15][i]=ay;
    Zm[15][i]=az+dx;
    MaXDevNum[15]+=1;
  }
  //direction 17
  for(int i=1;i<max_go;i++){
    int dx=(int)(i*1.0/sqrt(2.0));
    int dy=-(int)(i*1.0/sqrt(2.0));
    
    Xm[16][i]=ax+dy;
    Ym[16][i]=ay;
    Zm[16][i]=az+dx;
    MaXDevNum[16]+=1;
  }
  //direction 18
  for(int i=1;i<max_go;i++){
    int dx=-(int)(i*1.0/sqrt(2.0));
    int dy=-(int)(i*1.0/sqrt(2.0));
    
    Xm[17][i]=ax+dy;
    Ym[17][i]=ay;
    Zm[17][i]=az+dx;
    MaXDevNum[17]+=1;
  }
  for(int i=0;i<directions;i++){
    for(int j=1;j<max_go;j++){
      int x = Xm[i][j];
      int y = Ym[i][j];
      int z = Zm[i][j];
      Intm[i][j]=Inten(tools,x,y,z,quant_rad);
      //printf("moved %d %d %d %f\n",x,y,z,Intm[i][j]);
    }
  }
  //search max and min 
  for(int i=0;i<directions;i++){
    for(int j=1;j<max_go;j++){
      if(Intm[i][j]>IntMax[i]){
	IntMax[i]=Intm[i][j];
      }
      if(Intm[i][j]<IntMin[i]){
	IntMin[i]=Intm[i][j];
	IntMinIndex[i]=j;
      }
    }
    //printf("min max %f %f\n",IntMax[i],IntMin[i]);
  }
  //check if max is lower than a threshold
  int check_count=0;
  for(int i=0;i<directions;i++){
    if(IntMax[i]>=1.3*Intorig){
      CheckD[i]=0;
      check_count++;
    }
  }
  //find theminimumof the minimum
  double iniint=10000000000000000.0;
  for(int i=0;i<directions;i++){
    double curint=Intm[i][IntMinIndex[i]];
    if(curint<iniint&&CheckD[i]!=0){
      iniint=curint;
      Index_MinDirec=i;
    }
  }
  //check if all direction CheckD is not zero
  if(check_count==directions){
    //    printf("unko\n");
    return 1;
  }
  //move 
  int ii =Index_MinDirec;
  int x = Xm[ii][IntMinIndex[ii]];
  int y = Ym[ii][IntMinIndex[ii]];
  int z = Zm[ii][IntMinIndex[ii]];
  i_NucZ[my_rank] = z;
  i_NucY[my_rank] = y;
  i_NucX[my_rank] = x;
  double intee=Inten(tools,x,y,z,quant_rad);
  if(intee>240){
    return 1;
  }
  return 1;
  //printf("%d %d %d new_nuclei_position \n",x,y,z);
}

double Inten(const TOOL_t * const tools, int xx, int yy, int zz,int rad){
  const int plane_number = DIMZ;
  const int my_image_width = DIMX;
  const int my_image_height = DIMY;
  int j;
  int dist;
  int x, y, z;
  int base, x_start, x_end;
  double signal;
  double tot=0.0;
  int nuc_size = 2*rad;
  double count_intens=0.0;
  int *const *const x_range = tools->spheres[nuc_size];
  const int *const radii = tools->s_radii[nuc_size];
  const int slices = tools->s_layers[nuc_size];
  x = xx;
  y = yy;
  z = zz;
  int base_line = 0;
  /* base_line is set at least once */
  if(slices <= 0){
    printf("sliceout\n");
    exit(0);
  }
  //see z planes
  for (j = 0 - slices; j <= slices; j++)
  {
    int k;
    if (z + j < 0 || z + j >= plane_number) continue;
    const int r = radii[slices + j];
    base_line = y;
    //ydirection movement
    for (k = 0 - r; k <= r; k++)
    {
      int m;
      if (base_line + k < 0)
        continue;
      if (base_line + k >= my_image_height)
        break;
      dist = x_range[slices + j][r + k];
      x_start = x - dist > 0 ? x - dist : 0;
      x_end = x + dist < my_image_width ? x + dist : my_image_width;
      //x direction movement
      for (m = x_start; m <= x_end; m++)
      {
        tot +=  dd_out[Indix(z+j,y+k,m)];
        count_intens=count_intens+1.0;
      }
    }
  }
  tot/=count_intens;
  return tot;
}



void InitializeFrontPositionCellbyCellMPIforSF(int ca,int my_rank){
  //sphere
  
  int n=0;
  switch (ca) {
  case 1:

    break;
  case 2:
    int az = i_NucZ[my_rank];
    int ay = i_NucY[my_rank];
    int ax = i_NucX[my_rank];
    //printf("%d %d %d\n",az,ay,ax);
    int flagg=0;
    for(int ii=0;ii<GRIDSIZE_Z;ii++){
      for(int jj=0;jj<GRIDSIZE_Y;jj++){
	for(int kk=0;kk<GRIDSIZE_X;kk++){
	  int rk = my_rank;
	  int z = i_NucZ[rk];
	  int y = i_NucY[rk];
	  int x = i_NucX[rk];
	  float dx = (float)(x-kk);
	  float dy = (float)(y-jj);
	  float dz = (float)(z-ii);
	  dz/=z_factor;
	  float dist2=dx*dx+dy*dy+dz*dz;
	  if(dist2<0.0){
	    dist2=0.0;
	  }
	  if(Indix(ii,jj,kk)>=dims[4]){
	    printf("out range\n");
	    exit(0);
	  }
	  float dist=sqrt(dist2);
	   d_mask[Indix(ii,jj,kk)]=1;
	  if(dist<LIM_RAD){
	    d_mask[Indix(ii,jj,kk)]=0;
	    n++;
	    flagg=1;
	  }
	  else{ 
	    //phi[ii][jj][kk]=DWBAND;
	  }
	}   
      } 
    }
    //printf("flagg %d %d\n",flagg,n);
    //    SetSpeedFunction(1,2);
    break;
  }
}

