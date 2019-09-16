#include <stdio.h>
#include <sys/stat.h>
#include <dirent.h>

int main(){
  FILE * fp;
  printf("mkdir%d\n", mkdir("test_dir_01", 0777));
  fp = fopen("test_dir_01/test_file_01.txt", "w");
  fprintf(fp,"unko");
  fclose(fp);
}
