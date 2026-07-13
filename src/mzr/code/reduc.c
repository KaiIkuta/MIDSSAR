#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fitsio.h"

#define LINEPIX 2160
#define LINEPIY 1280
#define THRE 15000
#define HW 10
#define CONVF 0.083
#define RN 3.2
#define NAG 3

void printerror( int status);

void sweep(double p[NAG+1][NAG+1],double q[NAG+1][NAG+1])
{
  int i,j,k,im;
  double w,at,r[NAG+1][NAG+1];

  for(i=0;i<=NAG;i++){
    for(j=0;j<=NAG;j++){
      r[i][j]=p[i][j];
      q[i][j]=0.;
    }
    q[i][i]=1.;
  }
  for(k=0;k<=NAG;k++){
    w=0.;
    for(i=k;i<=NAG;i++){
      if(fabs(w)<fabs(r[i][k])){
        w=r[i][k];
        im=i;
      }
    }
    if(fabs(w)<1e-10){
      fprintf(stderr,"No inverse matrix ...\n");
      exit(0);
    }
    for(j=0;j<=NAG;j++){
      at=r[im][j];r[im][j]=r[k][j];r[k][j]=at;
      at=q[im][j];q[im][j]=q[k][j];q[k][j]=at;
    }
    for(j=0;j<=NAG;j++){
      r[k][j]/=w;
      q[k][j]/=w;
    }    
    for(i=0;i<=NAG;i++){
      if(i==k) continue;
      w=r[i][k];
      for(j=0;j<=NAG;j++){
        r[i][j]-=r[k][j]*w;
        q[i][j]-=q[k][j]*w;
      }
    }
    //for(i=0;i<=NAG;i++){
    //for(j=0;j<=NAG;j++) fprintf(stderr," %g",r[i][j]);
    //for(j=0;j<=NAG;j++) fprintf(stderr," %g",q[i][j]);
    //fprintf(stderr,"\n");
    //}getchar();
  }
}

void lsfit(float z[LINEPIY][LINEPIX], float za[NAG][LINEPIY][LINEPIX], double a[], int icen)
{
  int i,j,ii,jj,n;
  double p[NAG+1][NAG+1],q[NAG+1][NAG+1],r[NAG+1];

  for(ii=0;ii<NAG+1;ii++){
    r[ii]=0;
    for(jj=0;jj<NAG+1;jj++) p[ii][jj]=0.;
  }
  n=0;
  for(i=0;i<LINEPIY;i++){
    if(abs(i-icen)<150) continue;
    for(j=0;j<LINEPIX;j++){
      if(isnan(z[i][j])==0){
	for(ii=0;ii<NAG;ii++){
	  for(jj=ii;jj<NAG;jj++) p[ii][jj]+=za[ii][i][j]*za[jj][i][j];
	  p[ii][NAG]+=za[ii][i][j];
	  r[ii]+=za[ii][i][j]*z[i][j];
	}
	r[NAG]+=z[i][j];
	n++;
      }
    }
  }
  for(ii=0;ii<=NAG;ii++){
    r[ii]/=n;
    for(jj=ii;jj<=NAG;jj++) p[ii][jj]/=n;
  }
  for(ii=0;ii<=NAG;ii++){
    for(jj=0;jj<ii;jj++) p[ii][jj]=p[jj][ii];
  }
  p[NAG][NAG]=1.;
  sweep(p,q);
  for(ii=0;ii<=NAG;ii++){
    a[ii]=0.;
    for(jj=0;jj<=NAG;jj++) a[ii]+=q[ii][jj]*r[jj];
  }
}

float med(float z[LINEPIY][LINEPIX],int ic,int jc){
  float m[HW+1],md;
  int j,k,kk,n;

  for(k=0;k<=HW;k++) m[k]=-1e10;
  n=-1;
  for(j=jc-HW;j<=jc+HW;j++){
    if(j<HW||j>LINEPIX-HW-1) continue;
    if(isnan(z[ic][j])) continue;
    n++;
    for(k=0;k<=HW;k++){
      if(m[k]<z[ic][j]){
	for(kk=HW;kk>k;kk--) m[kk]=m[kk-1];
	m[k]=z[ic][j];
	break;
      }
    }
  }
  if(n%2) md=(m[n/2+1]+m[n/2])/2;
  else md=m[n/2];
  //for(k=0;k<=HW;k++) fprintf(stderr,"%d: %g\n",k,m[k]);
  //fprintf(stderr,"n=%d md=%g",n,md);getchar();
  return(md);
}

int medcut(float z[LINEPIY][LINEPIX],int ic,int jc){
  int i,j,k,n;

  k=0;
  if(z[ic][jc]>med(z,ic,jc)+THRE/3){
    z[ic][jc]=0./0.;
    for(k=1;k<=LINEPIX;k++){
      n=0;
      for(i=ic-k;i<=ic+k;i++){
	if(i<=0 || i>=LINEPIY-1) continue;
	for(j=jc-k;j<=jc+k;j++){
	  if(j<=0 || j>=LINEPIX-1) continue;
	  if(i==ic-k||i==ic+k||j==jc-k||j==jc+k){
	    if(z[i][j]>med(z,i,j)+THRE/3&&(isnan(z[i-1][j])||isnan(z[i+1][j])||isnan(z[i][j-1])||isnan(z[i][j+1]))){z[i][j]=0./0.;n++;}
	  }
	}
      }
      for(i=ic+k;i>=ic-k;i--){
	if(i<=0 || i>=LINEPIY-1) continue;
	for(j=jc+k;j>=jc-k;j--){
	  if(j<=0 || j>=LINEPIX-1) continue;
	  if(i==ic-k||i==ic+k||j==jc-k||j==jc+k){
	    if(z[i][j]>med(z,i,j)+THRE/3&&(isnan(z[i-1][j])||isnan(z[i+1][j])||isnan(z[i][j-1])||isnan(z[i][j+1]))){z[i][j]=0./0.;n++;}
	  }
	}
      }
      if(n==0) break;
    }
    for(i=ic-k;i<=ic+k;i++){
      if(i<=0 || i>=LINEPIY-1) continue;
      for(j=jc-k;j<=jc+k;j++){
	if(j<=0 || j>=LINEPIX-1) continue;
	if(isnan(z[i][j])==0&&(isnan(z[i-1][j])||isnan(z[i+1][j])||isnan(z[i][j-1])||isnan(z[i][j+1]))) z[i][j]=1./0.;
      }
    }
    for(i=ic-k;i<=ic+k;i++){
      if(i<=0 || i>=LINEPIY-1) continue;
      for(j=jc-k;j<=jc+k;j++){
	if(j<=0 || j>=LINEPIX-1) continue;
	if(isinf(z[i][j])) z[i][j]=0./0.;
      }
    }
  }
  return(k);
}

int main(int argc,char *argv[])
{
  fitsfile *fp;
  int i,j,ii,jj,k,m,n,r,rmax,cam,icen,gn,status,nfound,bitpix;
  long fpixel,naxis=2,naxes[2]={LINEPIX,LINEPIY};
  unsigned short zi[LINEPIX];
  static float zo[LINEPIY][LINEPIX],zd[LINEPIY][LINEPIX],za[NAG][LINEPIY][LINEPIX],zf[LINEPIY][LINEPIX];
  double bg[2][LINEPIX];
  FILE *fpt;
  char fname[64],fname0[54],fname1[54],buf[256],dir[32],jstend[16],date[16],gain[16],cn[2]={'C','H'},fene[16],object[32];
  double a[NAG+1],gf,dmy;
  float exptime,tmpval;
 
  if(argc!=3&&argc!=4&&argc!=5){
    printf("Usage: %s input output [darkimage [flatimage]]\n",argv[0]);
    exit(1);
  }
  if(argc>3){
    if((fpt = fopen("camnum","r")) == NULL){
      fprintf(stderr,"camnum file is not found...\n");
      exit(0);
    }
    fscanf(fpt,"%d",&cam);
    if(cam!=1) cam=0;
    fclose(fpt);
    sscanf(argv[0],"%s",dir);
    dir[strlen(dir)-5]='\0';
    sprintf(fname,"%scal_%c.dat",dir,cn[cam]);
    fprintf(stderr,"Read: %s ... ",fname);
    if((fpt=fopen(fname,"r"))==NULL){
      fprintf(stderr,"Can not open file: %s\n",fname);
      exit(0);
    }
    if(fgets(buf,sizeof(buf),fpt)==NULL){
      fprintf(stderr,"No data !\n");
      exit(1);
    }
    if(sscanf(buf,"%lf %lf %d",&dmy,&dmy,&icen)!=3){
      fprintf(stderr,"Broken file ?\n");
      exit(1);
    }
    fclose(fpt);
    fprintf(stderr,"OK\n");
    fprintf(stderr,"icen=%d\n",icen);
    fprintf(stderr,"Read %s ... ",argv[3]);
    status=0;
    if(fits_open_file(&fp, argv[3], READONLY, &status)) printerror(status);
    if(fits_read_keys_lng(fp, "NAXIS", 1, 2, naxes, &nfound, &status)) printerror(status);
    fprintf(stderr,"%ld x %ld format ",naxes[0],naxes[1]);
    if(naxes[0]!=LINEPIX||naxes[1]!=LINEPIY){
      fprintf(stderr,"\nIllegal image format\n");
      exit(0);
    }
    fpixel=1;
    for(i=0;i<naxes[1];i++){
      if(fits_read_img(fp, TFLOAT, fpixel, naxes[0], NULL, zd[i], NULL, &status)) printerror(status);
      fpixel+=naxes[0];
    }
    if(fits_close_file(fp, &status)) printerror(status);
    fprintf(stderr,"OK\n");
    for(k=0;k<NAG;k++){
      sprintf(fname,"%sampg%c%d.fits",dir,cn[cam],k+1);
      fprintf(stderr,"Read %s ... ",fname);
      status=0;
      if(fits_open_file(&fp, fname, READONLY, &status)) printerror(status);
      if(fits_read_keys_lng(fp, "NAXIS", 1, 2, naxes, &nfound, &status)) printerror(status);
      fprintf(stderr,"%ld x %ld format ",naxes[0],naxes[1]);
      if(naxes[0]!=LINEPIX||naxes[1]!=LINEPIY){
	fprintf(stderr,"\nIllegal image format\n");
	exit(0);
      }
      fpixel=1;
      for(i=0;i<naxes[1];i++){
	if(fits_read_img(fp, TFLOAT, fpixel, naxes[0], NULL, za[k][i], NULL, &status)) printerror(status);
	fpixel+=naxes[0];
      }
      if(fits_close_file(fp, &status)) printerror(status);
      fprintf(stderr,"OK\n");
    }
    if(argc>4){
      fprintf(stderr,"Read %s  ... ",argv[4]);
      status=0;
      if(fits_open_file(&fp, argv[4], READONLY, &status)) printerror(status);
      if(fits_read_keys_lng(fp, "NAXIS", 1, 2, naxes, &nfound, &status)) printerror(status);
      fprintf(stderr,"%ld x %ld format ",naxes[0],naxes[1]);
      if(naxes[0]!=LINEPIX||naxes[1]!=LINEPIY){
	fprintf(stderr,"\nIllegal image format\n");
	exit(0);
      }
      fpixel=1;
      for(i=0;i<naxes[1];i++){
	if(fits_read_img(fp, TFLOAT, fpixel, naxes[0], NULL, zf[i], NULL, &status)) printerror(status);
	fpixel+=naxes[0];
      }
      if(fits_close_file(fp, &status)) printerror(status);
      fprintf(stderr,"OK\n");
    }
    else{
      for(i=0;i<naxes[1];i++){
	for(j=0;j<naxes[0];j++) zf[i][j]=1.;
      }
    }
  }
  else{
    for(i=0;i<naxes[1];i++){
      for(j=0;j<naxes[0];j++) zd[i][j]=0.;
    }
  }
  sscanf(argv[1],"%s",fname0);
  sscanf(argv[2],"%s",fname1);
  
  k=0;
  while(1){
    if(k==0){
      status=0;
      fprintf(stderr,"Read: %s ... \n",fname0);
      if(fits_open_file(&fp, fname0, READONLY, &status)){
	fprintf(stderr,"Read: %s ... \n",fname0);
	if((fpt=fopen(fname0,"r"))==NULL){
	  fprintf(stderr,"Can not open file: %s\n",fname0);
	  exit(0);
	}
	if(fgets(buf,sizeof(buf),fpt)==NULL||feof(fpt)){
	  fprintf(stderr,"Empty file.\n");
	  fclose(fpt);
	  exit(0);
	}
	sscanf(buf,"%s",fname);
	fprintf(stderr,"\rRead: %s ... ",fname);
	status=0;
	if(fits_open_file(&fp, fname, READONLY, &status)) printerror(status);
	k++;
      }
    }
    else{
      if(fgets(buf,sizeof(buf),fpt)==NULL||feof(fpt)){
	fclose(fpt);
	exit(0);
      }
      sscanf(buf,"%s",fname);
      fprintf(stderr,"\rRead: %s ... ",fname);
      status=0;
      if(fits_open_file(&fp, fname, READONLY, &status)) printerror(status);
      k++;
    }
    if(fits_read_keys_lng(fp, "NAXIS", 1, 2, naxes, &nfound, &status)) printerror(status);
    //fprintf(stderr,"%ld x %ld format ",naxes[0],naxes[1]);
    if(naxes[0]!=LINEPIX||naxes[1]!=LINEPIY){
      fprintf(stderr,"\nIllegal image format @%s (%ld x %ld)\n",fname,naxes[0],naxes[1]);
      exit(0);
    }
    if(fits_read_key(fp, TINT, "BITPIX", &bitpix, NULL, &status)) printerror(status);
    if(bitpix!=16){
      fprintf(stderr,"Datatype is not TUSHORT. Maybe reduced data ...\n");
      sprintf(buf,"cp %s %s",fname0,fname1);
      fprintf(stderr,"\n%s\n",buf);
      system(buf);
      status=0;
      if(fits_close_file(fp, &status)) printerror(status);
      fprintf(stderr,"OK\n");
      exit(0);
    }
    fpixel=1;
    for(i=0;i<naxes[1];i++){
      if(fits_read_key(fp, TSTRING, "DATE-OBS", date, NULL, &status)) printerror(status);
      if(fits_read_key(fp, TSTRING, "JST-END", jstend, NULL, &status)) printerror(status);
      if(fits_read_key(fp, TSTRING, "OBJECT", object, NULL, &status)) printerror(status);
      if(fits_read_key(fp, TFLOAT, "EXPTIME", &exptime, NULL, &status)) printerror(status);
      if(fits_read_key(fp, TSTRING, "GAIN", gain, NULL, &status)) printerror(status);
      if(gain[0]!='x'){
	sscanf(gain,"%f",&tmpval);
	sprintf(gain,"x%d",(int)(16/tmpval*CONVF+0.5));
      }
      if(fits_read_key(fp, TSTRING, "FENE-LMP", fene, NULL, &status)) printerror(status);
      if(fits_read_img(fp, TUSHORT, fpixel, naxes[0], NULL, zi, NULL, &status)) printerror(status);
      for(j=0;j<naxes[0];j++) zo[i][j]=zi[j]-zd[i][j];
      fpixel+=naxes[0];
    }
    if(fits_close_file(fp, &status)) printerror(status);
    if(strncmp(gain,"x16",3)==0) gn=16;
    else if(strncmp(gain,"x8",2)==0) gn=8;
    else if(strncmp(gain,"x4",2)==0) gn=4;
    else if(strncmp(gain,"x2",2)==0) gn=2;
    else gn=1;
    gf=0.77+0.23*(16./gn)*(16./gn);
    fprintf(stderr,"gain=x%d (%g) OK\n",gn,gf);

    if(argc>3){
      rmax=0;n=0;
      for(i=1;i<naxes[1]-1;i++){
	for(j=1;j<naxes[0]-1;j++){
	  if(zo[i][j]>(zo[i-1][j]+zo[i+1][j]+zo[i][j-1]+zo[i][j+1])/4+THRE){
	    //fprintf(stderr,"%d %d: %g %g ",j+1,i+1,zo[i][j],(zo[i-1][j]+zo[i+1][j]+zo[i][j-1]+zo[i][j+1])/4);
	    r=medcut(zo,i,j);
	    //fprintf(stderr,"r=%d\n",r);
	    if(rmax<r) rmax=r;
	    n++;
	  }
	}
      }
      fprintf(stderr,"N_CR=%d  R_CRmax=%d\n",n,rmax);
      if(argc>4){
	for(j=0;j<naxes[0];j++){
	  for(m=0;m<2;m++){
	    i=icen+400*(m-0.5);
	    bg[m][j]=0.;n=0;
	    for(jj=j-10;jj<=j+10;jj++){
	      if(jj<0||jj>=naxes[0]) continue;
	      for(ii=i-50;ii<=i+50;ii++){
		if(isnan(zo[ii][jj])) continue;
		bg[m][j]+=(zo[ii][jj]+RN*RN/CONVF*gf)/zf[ii][jj];
		n++;
	      }
	    }
	    bg[m][j]/=n;
	  }
	}
      }
      lsfit(zo,za,a,icen);
      fprintf(stderr,"Ampglow correction factor: ");
      for(m=0;m<NAG;m++) fprintf(stderr,"%+g*ampg%c%d.fits",a[m],cn[cam],m+1);
      fprintf(stderr,"%+g\n",a[NAG]);
      for(i=0;i<naxes[1];i++){
	for(j=0;j<naxes[0];j++){
	  for(m=0;m<NAG;m++) zo[i][j]-=a[m]*za[m][i][j];
	  zo[i][j]-=a[NAG];
	}
      }
      if(argc>4){
	for(i=0;i<naxes[1];i++){
	  if(i>=naxes[1]-2){
	    for(j=0;j<naxes[0];j++) zo[i][j]=bg[i-(naxes[1]-2)][j];
	  }
	  else{
	    for(j=0;j<naxes[0];j++) zo[i][j]/=zf[i][j];
	  }
	}
      }
    }
    if(k==0) sprintf(fname,"%s",fname1);
    else sprintf(fname,"%s_%03d.fits",fname1,k);
    remove(fname);
    fprintf(stderr,"Write %s ... ",fname);
    status=0;bitpix=-32;
    if(fits_create_file(&fp, fname, &status)) printerror(status);
    if(fits_create_img(fp, bitpix, naxis, naxes, &status)) printerror(status);
    if(fits_write_key(fp, TSTRING, "DATE-OBS", date, NULL, &status)) printerror(status);
    if(fits_write_key(fp, TSTRING, "JST-END", jstend, NULL, &status)) printerror(status);
    if(fits_write_key(fp, TSTRING, "OBJECT", object, NULL, &status)) printerror(status);
    if(fits_write_key(fp, TFLOAT, "EXPTIME", &exptime, NULL, &status)) printerror(status);
    if(fits_write_key(fp, TSTRING, "GAIN", gain, NULL, &status)) printerror(status);
    if(fits_write_key(fp, TSTRING, "FENE-LMP", fene, NULL, &status)) printerror(status);
    fpixel=1;
    fprintf(stderr,"%ld x %ld format ",naxes[0],naxes[1]);
    for(i=0;i<naxes[1];i++){
      if(fits_write_img(fp, TFLOAT, fpixel, naxes[0], zo[i], &status)) printerror(status);
      fpixel+=naxes[0];
    }
    if(fits_close_file(fp, &status)) printerror(status);
    fprintf(stderr,"OK\n");
    if(k==0) break;
  }
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
