#include <iostream>
#include "tiffio.h"
#include "tiff.h"
#include "tiffvers.h"
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <cstdio>
#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <time.h>
bool startsWith(const char *pre, const char *str)
{
    size_t lenpre = strlen(pre),
           lenstr = strlen(str);
    return lenstr < lenpre ? false : strncmp(pre, str, lenpre) == 0;
}


int main(void){
  char unko[100]="abd";
  char shikko[100]="abcde";
  
  if(startsWith(unko,shikko)){
    printf("true\n");
  }
  else{
    printf("false\n");
  }
}
