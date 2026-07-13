#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fitsio.h"

#define LINEPIX 2160
#define LINEPIY 1280
#define NCUBE 1000
#define SATURATE 64000
#define CX 2000
#define CY 200
#define THRE 100
#define CONVF 0.083

void printerror( int status);

int main(int argc,char *argv[])
{
  fitsfile *fp[NCUBE],*fpd;
  int i,j,k,imin,jmin,i1,i2,j1,j2,i1t,i2t,j1t,j2t,ns,status,nfound,bitpix=-32;
  long fpixel,naxis=3,naxes[3]={LINEPIX, LINEPIY, NCUBE};
  float sig=3.,dk[LINEPIX],av0,sg0,exptime,tmpval;
  unsigned short zi[LINEPIX],nsmin,argc_dark=0,argc_flat=0;
  unsigned long nsat;
  double zs,z2s,ks,k2s,kzs,zx[LINEPIX],zy[LINEPIY];
  static float a[LINEPIY][LINEPIX],b[LINEPIY][LINEPIX],sg[LINEPIY][LINEPIX],z[NCUBE][LINEPIX];
  FILE *fpt;
  int lp,lpmax,nsig=1;
  char fname[64],fname0[54],buf[256],jstend[16],date[16],gain[16],fene[16],object[32];
 
  if(argc!=3&&argc!=4&&argc!=5&&argc!=6){
    printf("Usage: %s input output [sigma] [darkimage] [flatimage]\n",argv[0]);
    exit(1);
  }
  sscanf(argv[1],"%s",fname0);
  if(argc>=4){
    if(strcmp(argv[3]+strlen(argv[3])-4,"fits")==0){
      argc_dark=3;
      if(argc>=5){
	if(strcmp(argv[4]+strlen(argv[4])-4,"fits")==0){
	  argc_flat=4;
	  if(argc==6) sscanf(argv[5],"%f",&sig);
	}
	else{
	  sscanf(argv[4],"%f",&sig);
	  if(argc==6) argc_flat=5;
	}
      }
    }
    else{
      sscanf(argv[3],"%f",&sig);
      if(argc>=5) argc_dark=4;
      if(argc==6) argc_flat=5;
    }
  }
  fprintf(stderr,"Clipping threshold: %g sigma\n",fabs(sig));
  if(sig<0.) fprintf(stderr,"No fit mode\n");
  if(argc_dark) fprintf(stderr,"Dark image: %s\n",argv[argc_dark]);
  if(argc_flat) fprintf(stderr,"Flat image: %s\n",argv[argc_flat]);
  
  k=0;
  status=0;
  if(fits_open_file(&fp[k], fname0, READONLY, &status)==0){
    fprintf(stderr,"\rRead: %s ... ",fname0);
    if(fits_read_keys_lng(fp[k], "NAXIS", 1, 2, naxes, &nfound, &status)) printerror(status);
    //fprintf(stderr,"%ld x %ld format ",naxes[0],naxes[1]);
    if(naxes[0]>LINEPIX||naxes[1]>LINEPIY){
      fprintf(stderr,"\nIllegal image format @%s (%ld x %ld)\n",fname,naxes[0],naxes[1]);
      exit(0);
    }
    k++;
  }
  if(k==0){
    while(1){
      sprintf(fname,"%s_%03d.fits",fname0,k+1);
      status=0;
      if(fits_open_file(&fp[k], fname, READONLY, &status)) break;
      fprintf(stderr,"\rRead: %s ... ",fname);
      if(fits_read_keys_lng(fp[k], "NAXIS", 1, 2, naxes, &nfound, &status)) printerror(status);
      //fprintf(stderr,"%ld x %ld format ",naxes[0],naxes[1]);
      if(naxes[0]>LINEPIX||naxes[1]>LINEPIY){
	fprintf(stderr,"\nIllegal image format @%s (%ld x %ld)\n",fname,naxes[0],naxes[1]);
	exit(0);
      }
      k++;
    }
  }  
  if(k==0){
    fprintf(stderr,"Read: %s ... \n",argv[1]);
    if((fpt=fopen(argv[1],"r"))==NULL){
      fprintf(stderr,"Can not open file: %s\n",argv[1]);
      exit(0);
    }
    k=0;
    while(!feof(fpt)){
      if(fgets(buf,sizeof(buf),fpt)==NULL) break;
      sscanf(buf,"%s",fname);
      fprintf(stderr,"\rRead: %s ... ",fname);
      status=0;
      if(fits_open_file(&fp[k], fname, READONLY, &status)) break;
      if(fits_read_keys_lng(fp[k], "NAXIS", 1, 2, naxes, &nfound, &status)) printerror(status);
      //fprintf(stderr,"%ld x %ld format ",naxes[0],naxes[1]);
      if(naxes[0]>LINEPIX||naxes[1]>LINEPIY){
	fprintf(stderr,"\nIllegal image format @%s (%ld x %ld)\n",fname,naxes[0],naxes[1]);
	exit(0);
      }
      k++;
    }
    fclose(fpt);
    fprintf(stderr,"OK\n");
    nsig=0;
  }
  if(k==0) printerror(status);
  else status=0;
  naxes[2]=k;
  if(fits_read_key(fp[0], TSTRING, "DATE-OBS", date, NULL, &status)) printerror(status);
  if(fits_read_key(fp[0], TSTRING, "JST-END", jstend, NULL, &status)) printerror(status);
  if(fits_read_key(fp[0], TSTRING, "OBJECT",object, NULL, &status)) printerror(status);
  if(fits_read_key(fp[0], TFLOAT, "EXPTIME", &exptime, NULL, &status)) printerror(status);
  if(fits_read_key(fp[0], TSTRING, "GAIN", gain, NULL, &status)) printerror(status);
  if(gain[0]!='x'){
    sscanf(gain,"%f",&tmpval);
    sprintf(gain,"x%d",(int)(16/tmpval*CONVF+0.5));
  }
  if(fits_read_key(fp[0], TSTRING, "FENE-LMP", fene, NULL, &status)) printerror(status);
  for(k=0;k<naxes[2];k++){
    lpmax=naxes[2]-1;
    if(lpmax>3) lpmax=3;
    for(i=0;i<naxes[1];i++){
      for(j=0;j<naxes[0];j++){
	a[i][j]=0.;
	b[i][j]=0.;
	sg[i][j]=1e10;
      }
    }
    nsmin=naxes[2];
    nsat=0;
    for(i=0;i<naxes[1];i++){
      fpixel=i*naxes[0]+1;
      for(k=0;k<naxes[2];k++){
	if(fits_read_img(fp[k], TUSHORT, fpixel, naxes[0], NULL, zi, NULL, &status)){
	  status=0;
	  if(fits_read_img(fp[k], TFLOAT, fpixel, naxes[0], NULL, z[k], NULL, &status)) printerror(status);	  
	}
	else{
	  for(j=0;j<naxes[0];j++) z[k][j]=zi[j];
	}
      }
      for(j=0;j<naxes[0];j++){
	for(lp=lpmax;lp>=0;lp--){
	  zs=0.;
	  z2s=0.;
	  ks=0.;
	  k2s=0.;
	  kzs=0.;
	  ns=0;
	  for(k=0;k<naxes[2];k++){
	    //if(j+1==1084 && i+1==605) fprintf(stderr,"%d/%d: %d %d %g %g\n",lp,lpmax,k,z[k][j],sg[i][j],a[i][j]*(k-naxes[2]/2)+b[i][j]-z[k][j]);
	    if(lp==0&&z[k][j]>SATURATE) nsat++;
	    if(fabs(a[i][j]*(k-naxes[2]/2)+b[i][j]-z[k][j])<fabs(sig)*sg[i][j]+1e-5){
	      zs+=z[k][j];
	      z2s+=(double)z[k][j]*z[k][j];
	      ks+=(k-naxes[2]/2);
	      k2s+=(k-naxes[2]/2)*(k-naxes[2]/2);
	      kzs+=(double)(k-naxes[2]/2)*z[k][j];
	      ns++;
	    }
	    //if(j+1==955 && i+1==1279) fprintf(stderr,"%ld %lf %lf %d %lf %d\n",k-naxes[2]/2,a[i][j],b[i][j],z[k][j],sg[i][j],ns);
	  }
	  zs/=ns;
	  z2s/=ns;
	  ks/=ns;
	  k2s/=ns;
	  kzs/=ns;
	  if(naxes[2]!=1&&sig>0.){
	    a[i][j]=(kzs-ks*zs)/(k2s-ks*ks);
	    b[i][j]=(k2s*zs-kzs*ks)/(k2s-ks*ks);
	    sg[i][j]=sqrt(k2s*a[i][j]*a[i][j]+(double)b[i][j]*b[i][j]+z2s-2*kzs*a[i][j]-2*zs*b[i][j]+2*ks*a[i][j]*b[i][j]+1e-6);
	    //if(j+1==539 && i+1==717) fprintf(stderr,"%d %g %g %g %g %g %g %g %g\n",ns,zs,z2s,ks,k2s,kzs,a[i][j],b[i][j],sg[i][j]);
	  }
	  else{
	    a[i][j]=0.;
	    b[i][j]=zs;
	    sg[i][j]=sqrt(z2s-zs*zs);
	  }
	  if(isnan(a[i][j])||isnan(b[i][j])||isnan(sg[i][j])) fprintf(stderr,"NaN: (%d,%d) %g %g %g\n",j+1,i+1,a[i][j],b[i][j],sg[i][j]);
	  //if(j+1==955 && i+1==1279){
	  //fprintf(stderr,"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",a[i][j],b[i][j],sg[i][j],k2s,z2s,kzs,zs,ks,a[i][j]*a[i][j]*k2s+(double)b[i][j]*b[i][j]+z2s-2*a[i][j]*kzs-2*b[i][j]*zs+2*a[i][j]*b[i][j]*ks);
	  //fprintf(stderr,"%lf %lf %lf %lf %lf %lf\n",a[i][j]*a[i][j]*k2s,(double)b[i][j]*b[i][j],z2s,-2*a[i][j]*kzs,-2*b[i][j]*zs,2*a[i][j]*b[i][j]*ks);
	  //}
	  //if(lp==0 && isnan(sg[i][j])) fprintf(stderr,"%d %d %g %g\n",j,i,a[i][j],b[i][j]);
	}
	//if(i==naxes[1]/2 && j==naxes[0]/2) fprintf(stderr,"y=%g*x%+g  %g %g %g %g %g %g\n",a[i][j],b[i][j],zs,z2s,ks,k2s,kzs,sg[i][j]);
	nsmin=ns;
	jmin=j;
	imin=i;
	//if(i==naxes[1]/2 && j==naxes[0]/2) fprintf(stderr,"\n %d %g %g %g %g %g\n",ns,zs,z2s,a[i][j],b[i][j],sg[i][j]);
      }
    }
  }
  for(k=0;k<naxes[2];k++){
    if(fits_close_file(fp[k], &status)) printerror(status);
  }
  
  if(lpmax) fprintf(stderr,"Maximum rejected pixel: [%d,%d] %ld/%ld rejected\n",jmin+1,imin+1,naxes[2]-nsmin,naxes[2]);
  fprintf(stderr,"Saturated (>%d ADU) pixels: %g pix/frame\n",SATURATE,(double)nsat/naxes[2]);
  //zs=0.;
  //for(i=0;i<naxes[1];i++){
  //for(j=0;j<naxes[0];j++) zs+=b[i][j];
  //}
  //fprintf(stderr,"average=%g\n",zs/naxes[1]/naxes[0]);

  i1t=0;i2t=naxes[1];
  j1t=0;j2t=naxes[0];
  if(argc_dark){
    fprintf(stderr,"Read %s  ... ",argv[argc_dark]);
    status=0;
    if(fits_open_file(&fpd, argv[argc_dark], READONLY, &status)) printerror(status);
    if(fits_read_keys_lng(fpd, "NAXIS", 1, 2, naxes, &nfound, &status)) printerror(status);
    fprintf(stderr,"%ld x %ld format ",naxes[0],naxes[1]);
    if(naxes[0]>LINEPIX||naxes[1]>LINEPIY){
      fprintf(stderr,"\nIllegal image format\n");
      exit(0);
    }
    fpixel=1;
    for(i=0;i<naxes[1];i++){
      if(fits_read_img(fpd, TFLOAT, fpixel, naxes[0], NULL, dk, NULL, &status)) printerror(status);
      fpixel+=naxes[0];
      for(j=0;j<naxes[0];j++) b[i][j]-=dk[j];
    }
    if(fits_close_file(fpd, &status)) printerror(status);
    fprintf(stderr,"Subtraction OK\n");
    if(!argc_flat){
      for(i=0;i<naxes[1];i++){
	zy[i]=0.;
	for(j=0;j<naxes[0];j++) zy[i]+=b[i][j];
	zy[i]/=naxes[0];
      }
      for(i=0;i<naxes[1];i++) {
	//fprintf(stderr,"%d: %g",i,zy[i]);getchar();
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
	for(i=i1;i<=i2;i++) zx[j]+=b[i][j];
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
      j1t=0;
      j2t=naxes[0];
      fprintf(stderr,"Bright : [%d:%d,%d:%d]\n",j1+1,j2+1,i1+1,i2+1);
      fprintf(stderr,"Center : [%d:%d,%d:%d]\n",j1t+1,j2t+1,i1t+1,i2t+1);
    }
  }
  av0=0.;sg0=20000.;
  for(lp=5;lp>=0;lp--){
    zs=0.;
    z2s=0.;
    ns=0;
    for(i=i1t;i<i2t;i++){
      for(j=j1t;j<j2t;j++){
	if(sg[i][j]<av0+fabs(sig)*sg0+1e-5){
	  zs+=sg[i][j];
	  z2s+=(double)sg[i][j]*sg[i][j];
	  ns++;
	}
      }
    }
    zs/=ns;
    sg0=sqrt(z2s/ns-zs*zs);
    av0=zs;
    //fprintf(stderr,"%d: Sigma  : %g ± %g (Clipped: %d pix)\n",lp,av0,sg0,(i2t-i1t)*(j2t-j1t)-ns);
  }
  fprintf(stderr,"Sigma  : %g ± %g (Clipped: %d pix)\n",av0,sg0,(i2t-i1t)*(j2t-j1t)-ns);
  fpt=fopen("sigclip.log","a");
  fprintf(fpt,"%s %s %g %lf %lf %d ",fname0,argv[2],(double)nsat/naxes[2],av0,sg0,(i2t-i1t)*(j2t-j1t)-ns);
  if(nsig){
    for(i=0;i<naxes[1];i++){
      for(j=0;j<naxes[0];j++){
	if(i1t<=i&&i<i2t&&j1t<=j&&j<j2t){
	  if(isnan(sg[i][j])) sg[i][j]=1000;
	  else sg[i][j]=(sg[i][j]-av0)/(sg0+1e-5);
	}
	else sg[i][j]=0.;
      }
    }
  }
  av0=0.;sg0=20000.;
  for(lp=5;lp>=0;lp--){
    zs=0.;
    z2s=0.;
    ns=0;
    for(i=i1t;i<i2t;i++){
      for(j=j1t;j<j2t;j++){
	if(b[i][j]<av0+fabs(sig)*sg0+1e-5){
	  zs+=b[i][j];
	  z2s+=(double)b[i][j]*b[i][j];
	  ns++;
	}
      }
    }
    zs/=ns;
    sg0=sqrt(z2s/ns-zs*zs);
    av0=zs;
  }
  fprintf(stderr,"Average: %g ± %g (Clipped: %d pix)\n",av0,sg0,(i2t-i1t)*(j2t-j1t)-ns);
  fprintf(stderr,"Append: sigclip.log ... ");
  fprintf(fpt,"%lf %lf %d\n",av0,sg0,(i2t-i1t)*(j2t-j1t)-ns);
  fclose(fpt);
  fprintf(stderr,"OK\n");
  if(argc_dark){
    if(argc_flat){
      fprintf(stderr,"Read %s  ... ",argv[argc_flat]);
      status=0;
      if(fits_open_file(&fpd, argv[argc_flat], READONLY, &status)) printerror(status);
      if(fits_read_keys_lng(fpd, "NAXIS", 1, 2, naxes, &nfound, &status)) printerror(status);
      fprintf(stderr,"%ld x %ld format ",naxes[0],naxes[1]);
      if(naxes[0]>LINEPIX||naxes[1]>LINEPIY){
	fprintf(stderr,"\nIllegal image format\n");
	exit(0);
      }
      fpixel=1;
      for(i=0;i<naxes[1];i++){
	if(fits_read_img(fpd, TFLOAT, fpixel, naxes[0], NULL, dk, NULL, &status)) printerror(status);
	fpixel+=naxes[0];
	for(j=0;j<naxes[0];j++) b[i][j]/=dk[j];
      }
      if(fits_close_file(fpd, &status)) printerror(status);
      fprintf(stderr,"Flat fielding OK\n"); 
      for(i=naxes[1]-2;i<naxes[1];i++){
	for(j=0;j<naxes[0];j++) b[i][j]=0.;
      }
    }
    else{
      for(i=0;i<naxes[1];i++){
	for(j=0;j<naxes[0];j++) b[i][j]/=av0;
      }
      fprintf(stderr,"Normalization OK\n");
    }
  }
  else{
    for(i=naxes[1]-2;i<naxes[1];i++){
      for(j=0;j<naxes[0];j++) b[i][j]=0.;
    }
  }

  naxis=2;
  sscanf(argv[2],"%s",fname);
  if(strcmp(fname+strlen(fname)-4,"fits")!=0){
    strcat(fname,".fits");
  }
  remove(fname);
  fprintf(stderr,"Write %s ... ",fname);
  status=0;
  if(fits_create_file(&fpd, fname, &status)) printerror(status);
  if(fits_create_img(fpd, bitpix, naxis, naxes, &status)) printerror(status);
  if(fits_write_key(fpd, TSTRING, "DATE-OBS", date, NULL, &status)) printerror(status);
  if(fits_write_key(fpd, TSTRING, "JST-END", jstend, NULL, &status)) printerror(status);
  if(fits_write_key(fpd, TSTRING, "OBJECT",object, NULL, &status)) printerror(status);
  if(fits_write_key(fpd, TFLOAT, "EXPTIME", &exptime, NULL, &status)) printerror(status);
  if(fits_write_key(fpd, TSTRING, "GAIN", gain, NULL, &status)) printerror(status);
  if(fits_write_key(fpd, TSTRING, "FENE-LMP", fene, NULL, &status)) printerror(status);
  fpixel=1;
  fprintf(stderr,"%ld x %ld format ",naxes[0],naxes[1]);
  for(i=0;i<naxes[1];i++){
    if(fits_write_img(fpd, TFLOAT, fpixel, naxes[0], b[i], &status)) printerror(status);
    fpixel+=naxes[0];
  }
  if(fits_close_file(fpd, &status)) printerror(status);
  fprintf(stderr,"OK\n");
  remove("sigclip.fits");
  fprintf(stderr,"Write %s ... ","sigclip.fits");
  status=0;
  if(fits_create_file(&fpd, "sigclip.fits", &status)) printerror(status);
  if(fits_create_img(fpd, bitpix, naxis, naxes, &status)) printerror(status);
  fpixel=1;
  fprintf(stderr,"%ld x %ld format ",naxes[0],naxes[1]);
  for(i=0;i<naxes[1];i++){
    if(fits_write_img(fpd, TFLOAT, fpixel, naxes[0], sg[i], &status)) printerror(status);
    fpixel+=naxes[0];
  }
  if(fits_close_file(fpd, &status)) printerror(status);
  fprintf(stderr,"OK\n");
  return 0;
}

void printerror( int status)
{
    /*****************************************************/
    /* Print out cfitsio error messages and exit program */
    /*****************************************************/

    char status_str[FLEN_STATUS], errmsg[FLEN_ERRMSG];
  
    if(status)
      fprintf(stderr, "\n*** Error occurred during program execution ***\n");

    fits_get_errstatus(status, status_str);   /* get the error description */
    fprintf(stderr, "\nstatus = %d: %s\n", status, status_str);

    /* get first message; null ifstack is empty */
    if( fits_read_errmsg(errmsg) ) 
    {
         fprintf(stderr, "\nError message stack:\n");
         fprintf(stderr, " %s\n", errmsg);

         while ( fits_read_errmsg(errmsg) )  /* get remaining messages */
             fprintf(stderr, " %s\n", errmsg);
    }

    exit( status );       /* terminate the program, returning error status */
}
