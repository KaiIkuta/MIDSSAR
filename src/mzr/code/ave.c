#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fitsio.h"

#define LINEPIX 4848
#define LINEPIY 1280

void printerror( int status);

int main(int argc,char *argv[])
{
  fitsfile *fp;
  int i,j,k,status,nfound,anynull,bitpix=-32;
  long fpixel,naxis=2,naxes[2]={LINEPIX, LINEPIY};
  float nullval,exptime;
  static float z[LINEPIY][LINEPIX];
  static double zs[LINEPIY][LINEPIX];
  char jststr[16],date[16],gain[16],fene[16],object[32];
  
  if(argc<4){
    printf("Usage: %s input1 inpout2 ... output \n",argv[0]);
    exit(1);
  }
  for(i=0;i<LINEPIY;i++){
    for(j=0;j<LINEPIX;j++) zs[i][j]=0./0.;
  }
  for(k=1;k<argc-1;k++){
    fprintf(stderr,"Read: %s ... ",argv[k]);
    status=0;
    if(fits_open_file(&fp, argv[k], READONLY, &status)) printerror(status);
    if(fits_read_keys_lng(fp, "NAXIS", 1, 2, naxes, &nfound, &status)) printerror(status);
    fprintf(stderr,"%ld x %ld format ",naxes[0],naxes[1]);
    if((naxes[0]>LINEPIX)||(naxes[1]>LINEPIY)){
      fprintf(stderr,"\nIllegal image format\n");
      exit(0);
    }
    if(fits_read_key(fp, TSTRING, "DATE-OBS", date, NULL, &status)) printerror(status);
    if(fits_read_key(fp, TSTRING, "JST-END", jststr, NULL, &status)) printerror(status);
    if(fits_read_key(fp, TSTRING, "OBJECT", object, NULL, &status)) printerror(status);
    if(fits_read_key(fp, TFLOAT, "EXPTIME", &exptime, NULL, &status)) printerror(status);
    if(fits_read_key(fp, TSTRING, "GAIN", gain, NULL, &status)) printerror(status);
    if(fits_read_key(fp, TSTRING, "FENE-LMP", fene, NULL, &status)) printerror(status);
    fpixel = 1;
    nullval = 0;
    for(i=0;i<naxes[1];i++){
      if (fits_read_img(fp, TFLOAT, fpixel, naxes[0], &nullval, z[i], &anynull, &status)) printerror(status);
      for(j=0;j<naxes[0];j++){
	if(isnan(z[i][j])) continue;
	if(isnan(zs[i][j])) zs[i][j]=0.;
	zs[i][j]+=z[i][j];
      }
      fpixel+=naxes[0];
    }
    if(fits_close_file(fp, &status)) printerror(status);
    fprintf(stderr,"OK\n");
  }
  for(i=0;i<naxes[1];i++){
    for(j=0;j<naxes[0];j++) z[i][j]=zs[i][j]/(argc-2);
  }
  
  remove(argv[argc-1]);
  fprintf(stderr,"Write %s ... ",argv[argc-1]);
  status=0;
  if (fits_create_file(&fp, argv[argc-1], &status)) printerror(status);
  if (fits_create_img(fp, bitpix, naxis, naxes, &status)) printerror(status);
  if(fits_write_key(fp, TSTRING, "DATE-OBS", date, NULL, &status)) printerror(status);
  if(fits_write_key(fp, TSTRING, "JST-END", jststr, NULL, &status)) printerror(status);
  if(fits_write_key(fp, TSTRING, "OBJECT", object, NULL, &status)) printerror(status);
  if(fits_write_key(fp, TFLOAT, "EXPTIME", &exptime, NULL, &status)) printerror(status);
  if(fits_write_key(fp, TSTRING, "GAIN", gain, NULL, &status)) printerror(status);
  if(fits_write_key(fp, TSTRING, "FENE-LMP", fene, NULL, &status)) printerror(status);
  fpixel=1;
  fprintf(stderr,"%ld x %ld format ",naxes[0],naxes[1]);
  for(i=0;i<naxes[1];i++){
    if(fits_write_img(fp, TFLOAT, fpixel, naxes[0], z[i], &status)) printerror(status);
    fpixel+=naxes[0];
  }
  if (fits_close_file(fp, &status)) printerror(status);
  fprintf(stderr,"OK\n");
  return 0;
}

void printerror( int status)
{
    /*****************************************************/
    /* Print out cfitsio error messages and exit program */
    /*****************************************************/

    char status_str[FLEN_STATUS], errmsg[FLEN_ERRMSG];
  
    if (status)
      fprintf(stderr, "\n*** Error occurred during program execution ***\n");

    fits_get_errstatus(status, status_str);   /* get the error description */
    fprintf(stderr, "\nstatus = %d: %s\n", status, status_str);

    /* get first message; null if stack is empty */
    if ( fits_read_errmsg(errmsg) ) 
    {
         fprintf(stderr, "\nError message stack:\n");
         fprintf(stderr, " %s\n", errmsg);

         while ( fits_read_errmsg(errmsg) )  /* get remaining messages */
             fprintf(stderr, " %s\n", errmsg);
    }

    exit( status );       /* terminate the program, returning error status */
}
