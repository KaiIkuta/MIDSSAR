#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fitsio.h"

#define LINEPIX 2160
#define LINEPIY 1280
#define CX 2000
#define CY 200
#define THRE 0.5
#define MAXODR 9
#define BUFFER 0

void printerror( int status);

void sweep(double p[MAXODR][MAXODR],double q[MAXODR][MAXODR],int n)
{
  int i,j,k,im;
  double w,at,r[MAXODR][MAXODR];

  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      r[i][j]=p[i][j];
      q[i][j]=0.;
    }
    q[i][i]=1.;
  }
  for(k=0;k<n;k++){
    w=0.;
    for(i=k;i<n;i++){
      if(fabs(w)<fabs(r[i][k])){
        w=r[i][k];
        im=i;
      }
    }
    if(fabs(w)<1e-10){
      fprintf(stderr,"No inverse matrix ...\n");
      exit(0);
    }
    for(j=0;j<n;j++){
      at=r[im][j];r[im][j]=r[k][j];r[k][j]=at;
      at=q[im][j];q[im][j]=q[k][j];q[k][j]=at;
    }
    for(j=0;j<n;j++){
      r[k][j]/=w;
      q[k][j]/=w;
    }    
    for(i=0;i<n;i++){
      if(i==k) continue;
      w=r[i][k];
      for(j=0;j<n;j++){
        r[i][j]-=r[k][j]*w;
        q[i][j]-=q[k][j]*w;
      }
    }
    //for(i=0;i<n;i++){
    //for(j=0;j<n;j++) fprintf(stderr," %g",r[i][j]);
    //for(j=0;j<n;j++) fprintf(stderr," %g",q[i][j]);
    //fprintf(stderr,"\n");
    //}getchar();
  }
}

double lsfit(double y[], int ns, int ne, int od, double f[])
{
  int i,j,n;
  double sxn[MAXODR*2],xn[MAXODR*2],dy,dy2s;
  double p[MAXODR][MAXODR],q[MAXODR][MAXODR],sxny[MAXODR],a[MAXODR];

  xn[0]=1;
  sxn[0]=1;
  for(i=1;i<MAXODR*2;i++) sxn[i]=0;
  for(i=0;i<MAXODR;i++) sxny[i]=0;
  for(n=ns;n<ne;n++){
    //fprintf(stderr,"%d %g %g",n,(n-(ns+ne-1)/2),y[n]);getchar();
    for(i=1;i<MAXODR*2;i++) {
      xn[i]=xn[i-1]*(n-(ns+ne-1)/2);
      sxn[i]+=xn[i]/(ne-ns);}
    for(i=0;i<MAXODR;i++) sxny[i]+=xn[i]*y[n]/(ne-ns);
  }
  for(i=0;i<MAXODR;i++){
    for(j=0;j<MAXODR;j++) p[i][j]=sxn[i+j];
  }
  sweep(p,q,od+1);
  for(i=0;i<=od;i++){
    a[i]=0.;
    for(j=0;j<=od;j++){
      a[i]+=q[i][j]*sxny[j];
    }
  }
  dy2s=0.;
  for(n=0;n<LINEPIX;n++){
    f[n]=a[od];
    for(i=od-1;i>=0;i--) f[n]+=f[n]*(n-(ns+ne-1)/2)+a[i];
    if(ns<=n && n<ne){
      dy=y[n]-f[n];
      dy2s+=dy*dy;
    }
  }
  dy2s=sqrt(dy2s/(ne-ns));
  return(dy2s);
}

int main(int argc,char *argv[])
{
  fitsfile *fp;
  int i,j,status,nfound,anynull,bitpix=-32;
  long fpixel,naxis=2,naxes[2]={LINEPIX, LINEPIY};
  float nullval,exptime;
  static float z[LINEPIY][LINEPIX],zxy[LINEPIY][LINEPIX];
  double zx[LINEPIX],zy[LINEPIY],zf[LINEPIX];
  double cf[4][4],cfs;
  int i1,i2,j1,j2,i1t,i2t,j1t,j2t;
  char fname[64],fname0[54],jststr[16],date[16],gain[16],fene[16],object[32];
  
  if(argc!=2){
    printf("Usage: %s input\n",argv[0]);
    exit(1);
  }
  sprintf(fname0,"%s",argv[1]);
  sprintf(fname,"%s.fits",fname0);
  fprintf(stderr,"Read: %s ... ",fname);
  status=0;
  if(fits_open_file(&fp, fname, READONLY, &status)) printerror(status);
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
    fpixel+=naxes[0];
  }
  if(fits_close_file(fp, &status)) printerror(status);
  fprintf(stderr,"OK\n");
  for(i=0;i<naxes[1];i++){
    zy[i]=0.;
    for(j=0;j<naxes[0];j++) zy[i]+=z[i][j];
    zy[i]/=naxes[0];
  }
  for(i=0;i<naxes[1];i++) {
    if(zy[i]>THRE) break;
  }
  i1=i;
  for(i=naxes[1]-3;i>=0;i--) {
    if(zy[i]>THRE) break;
  }
  i2=i;
  if(i1>=i2){
    fprintf(stderr,"No bright area found ... Y[%d:%d]\n",i1,i2);
    exit(1);
  }
  for(j=0;j<naxes[0];j++){
    zx[j]=0.;
    for(i=i1;i<=i2;i++) zx[j]+=z[i][j];
    zx[j]/=i2-i1+1;
  }
  for(j=0;j<naxes[0];j++) {
    if(zx[j]>THRE) break;
  }
  j1=j;
  for(j=naxes[0]-1;j>=0;j--) {
    if(zx[j]>THRE) break;
  }
  j2=j;
  if(j1>=j2){
    fprintf(stderr,"No bright area found ... X[%d:%d]\n",j1,j2);
    exit(1);
  }
  i1t=(i1+i2-CY)/2;
  i2t=(i1+i2+CY)/2;
  j1t=(j1+j2-CX)/2;
  j2t=(j1+j2+CX)/2;
  if(j1t<0) j1t=j2t%4;
  if(j2t>naxes[0]) j2t=naxes[0]-4+j1t%4;
  fprintf(stderr,"Bright : [%d:%d,%d:%d]\n",j1+1,j2+1,i1+1,i2+1);
  fprintf(stderr,"Center : [%d:%d,%d:%d]\n",j1t+1,j2t+1,i1t+1,i2t+1);
  for(i=0;i<4;i++){
    for(j=0;j<4;j++) cf[i][j]=0.;
  }
  cfs=0.;
  for(i=i1t;i<i2t;i++){
    for(j=j1t;j<j2t;j++){
      cf[i%4][j%4]+=z[i][j];
      cfs+=z[i][j];
    }
  }
  for(i=0;i<4;i++){
    for(j=0;j<4;j++) cf[i][j]/=cfs/16;
  }
  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      fprintf(stderr," %lf",cf[i][j]);
    }
    fprintf(stderr,"\n");
  }
  for(i=0;i<naxes[1];i++){
    zy[i]=0.;
    for(j=0;j<naxes[0];j++){
      z[i][j]/=cf[i%4][j%4];
      zy[i]+=z[i][j];
    }
    zy[i]/=naxes[0];
  }
  for(j=0;j<naxes[0];j++){
    zx[j]=0.;
    for(i=i1t;i<=i2t;i++) zx[j]+=z[i][j];
    zx[j]/=i2-i1+1;
  }
  for(i=0;i<naxes[1];i++){
    if(i1-BUFFER<i && i<i2+BUFFER){
      for(j=0;j<naxes[0];j++){
	zxy[i][j]=zy[i]*zx[j];
	z[i][j]/=zxy[i][j];
      }
    }
    else{
      for(j=0;j<naxes[0];j++) z[i][j]=1.;
    }
  }
  for(i=i1-BUFFER;i<=i2+BUFFER;i++){
    for(j=0;j<naxes[0];j++) zx[j]=z[i][j];
    lsfit(zx,0,naxes[0],4,zf);
    for(j=0;j<naxes[0];j++){
      zxy[i][j]*=zf[j];
      z[i][j]/=zf[j];
    }
  }
  for(j=0;j<naxes[0];j++){
    for(i=0;i<naxes[1];i++) zy[i]=z[i][j];
    lsfit(zy,i1t,i2t,1,zf);
    for(i=i1-BUFFER;i<=i2+BUFFER;i++){
      zxy[i][j]*=zf[i];
      z[i][j]/=zf[i];
    }
  }
  for(i=0;i<naxes[1];i++){
    for(j=0;j<naxes[0];j++) z[i][j]*=cf[i%4][j%4];
  }
  sprintf(fname,"%s_xy.fits",fname0);
  remove(fname);
  fprintf(stderr,"Write %s ... ",fname);
  status=0;
  if (fits_create_file(&fp, fname, &status)) printerror(status);
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
  sprintf(fname,"%s_res.fits",fname0);
  remove(fname);
  fprintf(stderr,"Write %s ... ",fname);
  status=0;
  if (fits_create_file(&fp, fname, &status)) printerror(status);
  if (fits_create_img(fp, bitpix, naxis, naxes, &status)) printerror(status);
  fpixel=1;
  fprintf(stderr,"%ld x %ld format ",naxes[0],naxes[1]);
  for(i=0;i<naxes[1];i++){
    if(fits_write_img(fp, TFLOAT, fpixel, naxes[0], zxy[i], &status)) printerror(status);
    fpixel+=naxes[0];
  }
  if (fits_close_file(fp, &status)) printerror(status);
  fprintf(stderr,"OK\n");
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
