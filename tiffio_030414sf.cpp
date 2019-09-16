#include "mem_seg_levelset030414sf.h"
#include "tiffio.h"
#include "tiff.h"
#include "tiffvers.h"
void SetFilename(){
  int i,j;
  dir = opendir(ddirname);
  if(dir==NULL){
    printf("unko\n");
    exit(0);
  }
  for(i = 0; NULL != (dp = readdir(dir)); i++){
    strcpy(FILENAME[i],dp->d_name);
    count_file++;
  }
  char temp[MAX_CHAR_NUM];
  for(i =0; i < count_file-1; i++){
    for(j = i+1; j < count_file; j++){
      if ((strcmp(FILENAME[j],FILENAME[i])) < 0){
	strcpy(temp,FILENAME[i]);
	strcpy(FILENAME[i],FILENAME[j]);
	strcpy(FILENAME[j],temp);
      }
    }
  }
}
void SetFilename_G(){
  int i,j;
  dir_g = opendir(globalname);
  if(dir_g==NULL){
    printf("unko\n");
    exit(0);
  }
  for(i = 0; NULL != (dp_g = readdir(dir_g)); i++){
    strcpy(FILENAME_G[i],dp_g->d_name);
    count_file_g++;
  }
  char temp[MAX_CHAR_NUM];
  for(i =0; i < count_file_g-1; i++){
    for(j = i+1; j < count_file_g; j++){
      if ((strcmp(FILENAME_G[j],FILENAME_G[i])) < 0){
	strcpy(temp,FILENAME_G[i]);
	strcpy(FILENAME_G[i],FILENAME_G[j]);
	strcpy(FILENAME_G[j],temp);
      }
    }
  }
}
/*

*/
void ReadTIFFile2forSF(int kurikaeshi){
  int r=0;
  int rows=0;
  int columns;
  uint16 BitsPerSample;
  uint16 SamplesPerPixels;
  uint16 PhotoMetric;
  uint32 TifSize;
  uint8 *in_image;
  TIFF *in_filep;
  int check_io;
  char infile[MAX_CHAR_NUM];
  char outfile[MAX_CHAR_NUM];
  int n1, n2, n3, n4;
  int m1,m2,m3,m4;
  int  x, y,z;
  int slice_tot;
  int slice;
  int i,j;
  infile[0] = '\0'; 
  strcpy(infile, ddirname);
  strcat(infile, FILENAME[kurikaeshi+2]);
  //printf("in_name %s\n", infile);
  in_filep = TIFFOpen(infile, "r");
  
  if (in_filep == NULL) {
    cout << "could not open input file" << endl;
  }
  slice_tot = TIFFNumberOfDirectories(in_filep);
  TIFFGetField(in_filep, TIFFTAG_IMAGELENGTH, &rows);
  TIFFGetField(in_filep, TIFFTAG_IMAGEWIDTH, &columns);
  TIFFGetField(in_filep, TIFFTAG_BITSPERSAMPLE, &BitsPerSample);
  TIFFGetField(in_filep, TIFFTAG_SAMPLESPERPIXEL, &SamplesPerPixels);
  TIFFGetField(in_filep, TIFFTAG_PHOTOMETRIC, &PhotoMetric);
  //printf("rows = %d\n", rows);
  // printf("columns = %d %d\n", columns,slice_tot);
  // printf("BitsPerSample = %d\n", BitsPerSample);
  //printf("SamplesPerPixel = %d\n", SamplesPerPixels);
  //printf("PhotoMetric = %d\n", PhotoMetric);
  mdims[1]=(long)columns;
  mdims[0]=(long)rows;
  mdims[2]=(long)slice_tot;
  GRIDSIZE_X=(int)columns;
  GRIDSIZE_Y=(int)rows;
  GRIDSIZE_Z=(int)slice_tot;
  //
  int dimz, dimx,dimy;
  dims[2] = 1;
  dims[1] = 1;
  switch(numdims){
  case 3: dimz = (int)mdims[2]; dims[2] = dimz;
  case 2: dimx = (int)mdims[1]; dims[1] = dimx;
  case 1: dimy = (int)mdims[0]; dims[0] = dimy;
  }
  dims[3] = dims[0]*dims[1];
  dims[4] = dims[0]*dims[1]*dims[2];
  //fprintf(FileLog,"dims %ld %ld %ld %ld %ld\n",dims[0],dims[1],dims[2],dims[3],dims[4]);
  TifSize = rows* columns * sizeof( uint8 ); 
  in_image = ( uint8 *) malloc(TifSize);
  if (in_image == NULL) {
    cout << "could not alloctate memory" << endl;
  }
  
  for(slice = 0; slice < slice_tot; slice++){
    
    if ( ! TIFFSetDirectory( in_filep, slice ) ){
      printf("ERROR: Can not set TIFF directory [TIFFSetDirectory] \n");
      exit( EXIT_FAILURE );
    }
    
    for (r = 0; r <(int)rows; r++) {
      check_io = TIFFReadScanline(in_filep, &in_image[r*columns], r, 1);
      if (check_io != 1){
	cout<< "not read image from files" << endl;
      }
    }
    for (y=0;y<(int)rows; y++){
      for (x=0;x<(int)columns;x++){
	//if(Indix(slice,y,x)==32000-1){
	// exit(0);
	//}
	dd_out[Indix(slice,y,x)]=(double)in_image[y*columns+x]/256.0;
	
      }
    }
  }
 
  _TIFFfree(in_image);
  TIFFClose(in_filep);
  
}
void ReadTIFFilemaskforSF(int kurikaeshi){
  int r=0;
  int rows=0;
  int columns;
  uint16 BitsPerSample;
  uint16 SamplesPerPixels;
  uint16 PhotoMetric;
  uint32 TifSize;
  uint8 *in_image;
  TIFF *in_filep;
  int check_io;
  char infile[MAX_CHAR_NUM];
  char outfile[MAX_CHAR_NUM];
  int n1, n2, n3, n4;
  int m1,m2,m3,m4;
  int  x, y,z;
  int slice_tot;
  int slice;
  int i,j;
  infile[0] = '\0'; 
  strcpy(infile, globalname);
  strcat(infile, FILENAME_G[kurikaeshi+2]);
  //printf("g_name %s\n", infile);
  in_filep = TIFFOpen(infile, "r");
  
  if (in_filep == NULL) {
    cout << "could not open input file" << endl;
  }
  slice_tot = TIFFNumberOfDirectories(in_filep);
  TIFFGetField(in_filep, TIFFTAG_IMAGELENGTH, &rows);
  TIFFGetField(in_filep, TIFFTAG_IMAGEWIDTH, &columns);
  TIFFGetField(in_filep, TIFFTAG_BITSPERSAMPLE, &BitsPerSample);
  TIFFGetField(in_filep, TIFFTAG_SAMPLESPERPIXEL, &SamplesPerPixels);
  TIFFGetField(in_filep, TIFFTAG_PHOTOMETRIC, &PhotoMetric);
  //printf("rows = %d\n", rows);
  // printf("columns = %d %d\n", columns,slice_tot);
  // printf("BitsPerSample = %d\n", BitsPerSample);
  //printf("SamplesPerPixel = %d\n", SamplesPerPixels);
  //printf("PhotoMetric = %d\n", PhotoMetric);
  //
  int dimz, dimx,dimy;
 
  TifSize = rows* columns * sizeof( uint8 ); 
  in_image = ( uint8 *) malloc(TifSize);
  if (in_image == NULL) {
    cout << "could not alloctate memory" << endl;
  }
  
  for(slice = 0; slice < slice_tot; slice++){
    
    if ( ! TIFFSetDirectory( in_filep, slice ) ){
      printf("ERROR: Can not set TIFF directory [TIFFSetDirectory] \n");
      exit( EXIT_FAILURE );
    }
    
    for (r = 0; r <(int)rows; r++) {
      check_io = TIFFReadScanline(in_filep, &in_image[r*columns], r, 1);
      if (check_io != 1){
	cout<< "not read image from files" << endl;
      }
    }
    for (y=0;y<(int)rows; y++){
      for (x=0;x<(int)columns;x++){
	if(in_image[y*columns+x]>100){
	  //outside
	  //140403ver
	  //or 140406ver which means outline
	  // dd_out[Indix(slice,y,x)]=1.0;
	  //d_mask[Indix(slice,y,x)]=1.0;
	  //if(x>10&&x<columns-10&&y>10&&y<rows-10){
	  //dd_out[Indix(slice,y,x)]=1.0;
	    //}
	  //140528
	  dd_out[Indix(slice,y,x)]=1.0;
	}
	else{
	  //globalmembrane
       
	 
	}
      }
    }
  }
 
  _TIFFfree(in_image);
  TIFFClose(in_filep);
  
}
/*
void WriteTIFFile(int kurikaeshi){
  int    i,j;
  char   tiftitle[256] = { '\0' };
  float  resolution;
  uint16 samplesperpixel, bitspersample;
  uint16 *data;
  uint16 i_u16, j_u16;
  uint16 page_max;
  uint32 width, height;
  uint32 rowsperstrip;
  int32  tifsize; // tiffio.h
  TIFF   *tif;
  int n1, n2, n3, n4;
  int r=0;
  int rows=0;
  int columns;
  uint16 *in_image;
  TIFF *in_filep;
  int check_io;
  char infile[MAX_CHAR_NUM];
  char infile2[MAX_CHAR_NUM];
  int  x, y,z;
  width           = GRIDSIZE_X;
  height          = GRIDSIZE_Y;
  samplesperpixel = 1;
  bitspersample   = 16;
  resolution      = 72;
  page_max        = GRIDSIZE_Z;
  n1=kurikaeshi/1000;
  n2=(kurikaeshi-n1*1000)/100;
  n3=(kurikaeshi-n1*1000-n2*100)/10;
  n4=kurikaeshi-n1*1000-n2*100-n3*10;
  strcpy(infile, diroutname);
 
  sprintf(infile2,"output%d%d%d%d.tif",n1,n2,n3,n4);
  strcat(infile,infile2);
  tif = TIFFOpen(infile, "w");
  if ( tif == NULL ){
    printf("ERROR: Can not open tiff file \n");
    return;
  }
  tifsize = width * height * sizeof( uint16 );
  data    = ( uint16 * )malloc( tifsize );
  if ( data == NULL ){
    printf("ERROR: Can not get memory [data] \n");
    exit( EXIT_FAILURE );
  }
  double phi_max=-100000000000000000.0;
  double phi_min= 100000000000000000.0;
  for (i_u16 = 0; i_u16 < page_max; i_u16++){
    for(y=0;y<(int)height;y++){
      for(x=0;x<(int)width;x++){
	if(phi_max<phi[(int)i_u16][y][x]){
	  phi_max=phi[(int)i_u16][y][x];
	}
	if(phi_min>phi[(int)i_u16][y][x]){
	  phi_min=phi[(int)i_u16][y][x];
	}
      }
    }
  }
  for (i_u16 = 0; i_u16 < page_max; i_u16++){
    for(y=0;y<(int)height;y++){
      for(x=0;x<(int)width;x++){
	phi[(int)i_u16][y][x]=30000.0*(phi[(int)i_u16][y][x]-phi_min)/(phi_max-phi_min)+1.0;
      }
    }
  }
  for (i_u16 = 0; i_u16 < page_max; i_u16++){
    TIFFSetField(tif, TIFFTAG_IMAGELENGTH, height);
    TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, width);
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, bitspersample);
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, samplesperpixel);
    TIFFSetField(tif, TIFFTAG_XRESOLUTION, resolution);
    TIFFSetField(tif, TIFFTAG_YRESOLUTION, resolution);
    TIFFSetField(tif, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
    TIFFSetField(tif, TIFFTAG_PAGENUMBER, i_u16, page_max);
    for (r = 0; r <(int)rows; r++) {
      check_io = TIFFReadScanline(in_filep, &in_image[r*columns], r, 1);
      if (check_io != 1){
	cout<< "not read image from files" << endl;
      }
    }
    for (y=0;y<(int)height;y++){
      for (x=0;x<(int)width;x++){
	data[y*width+x]=phi[(int)i_u16][y][x];
	if(data[y*width+x]<0){
	  data[y*width+x]=0;
	}
	if(data[y*width+x]>65535){
	  data[y*width+x]=65535;
	}
      }
    }
    TIFFWriteEncodedStrip(tif, 0, data, tifsize);
    TIFFWriteDirectory( tif );
  }
}
*/
