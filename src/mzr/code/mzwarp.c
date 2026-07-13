#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fitsio.h"
#include "warp.c"

#define NCAL 80
#define MAXODR 9
#define YOD 4
#define XOD 7
#define THRE 100

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

double lsft(double x[], double y[], int ns, int ne, int od, double clip, double a[])
{
  int i,j,n,nd,lp=3;
  double sxn[MAXODR*2],xn[MAXODR*2],dy2s,sig=1e10;
  double p[MAXODR][MAXODR],q[MAXODR][MAXODR],sxny[MAXODR],f[LINEPIX2];

  xn[0]=1;
  sxn[0]=1;
  for(n=ns;n<ne;n++) f[n]=0.;
  if(clip>9) lp=1;
  for(;lp>0;lp--){
    for(i=1;i<MAXODR*2;i++) sxn[i]=0;
    for(i=0;i<MAXODR;i++) sxny[i]=0;
    nd=0.;
    for(n=ns;n<ne;n++){
      //fprintf(stderr,"%d %g %g %g",n,x[n],y[n],f[n]);getchar();
      if(isnan(x[n])||isnan(y[n])) continue;
      if(fabs(f[n])<clip*sig){
	for(i=1;i<MAXODR*2;i++){
	  xn[i]=xn[i-1]*x[n];
	  sxn[i]+=xn[i];
	}
	for(i=0;i<MAXODR;i++) sxny[i]+=xn[i]*y[n];
	nd++;
      }
    }
    for(i=1;i<MAXODR*2;i++) sxn[i]/=nd;
    for(i=0;i<MAXODR;i++){
      sxny[i]/=nd;
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
    nd=0;
    for(n=ns;n<ne;n++){
      f[n]=a[od];
      for(i=od-1;i>=0;i--) f[n]=f[n]*x[n]+a[i];
      f[n]=y[n]-f[n];
      if(fabs(f[n])<clip*sig){
	dy2s+=f[n]*f[n];
	nd++;
      }
    }
    sig=sqrt(dy2s/nd);
  }
  return(sig);
}

double lsfit(double y[], int ns, int ne, int od, double clip, double a[])
{
  int n;
  double x[LINEPIX2];

  for(n=ns;n<ne;n++) x[n]=n-(ns+ne-1.)/2;
  return(lsft(x,y,ns,ne,od,clip,a));
}

void wrp(float zi[LINEPIY2][LINEPIX2],float zx[LINEPIY2][LINEPIX2],float zy[LINEPIY2][LINEPIX2])
{
  int i,j;
  static float zo[LINEPIY2][LINEPIX2];

  warp(zi,zo,zx,zy);
  for(i=0;i<LINEPIY2;i++){
    for(j=0;j<LINEPIX2;j++) zi[i][j]=zo[i][j];
  }
}

int main(int argc,char *argv[])
{
  fitsfile *fp;
  int i,j,status,nfound,anynull,bitpix=-32;
  long fpixel,naxis=2,naxes[3]={LINEPIX2, LINEPIY2, 2};
  float nullval,tmpval;
  static float zi[LINEPIY][LINEPIX2],zo[LINEPIY2][LINEPIX2],zx[LINEPIY2][LINEPIX2],zy[LINEPIY2][LINEPIX2];
  int cam,nf,ii,it,icen,jj,k,nc,nd,n,nn,of,m,mmax,ofmax,sn,cl;
  double yc[NFIBER][LINEPIX2],zz[5],a[MAXODR],af[NFIBER][MAXODR],sig,sigs,aa[MAXODR][MAXODR],at[NFIBER],x,dy,dx;
  double kn[9]={0.2,0.5,0.8,1.,1.,1.,0.8,0.5,0.2},zt[LINEPIX2],ztb[LINEPIX2],cx[NCAL][NFIBER+1],cxt[NCAL],cw[NCAL],wt,wc,dw,dwt,zs;
  double xc[NFIBER][LINEPIX2],day[4],aw[4],ax[4],ay[4],th,yav[LINEPIY2],at2s,at2s0;
  FILE *fpt;
  char buf[256],fname[64],lname[8],dir[32];
  
  if(argc!=2&&argc!=3){
    printf("Usage: %s cont cal\n",argv[0]);
    printf("  or : %s sky\n",argv[0]);
    exit(0);
  }
  if((fpt = fopen("camnum","r")) == NULL){
    fprintf(stderr,"camnum file is not found...\n");
    exit(0);
  }
  fscanf(fpt,"%d",&cam);
  if(cam!=1) cam=0;
  fclose(fpt);
  sscanf(argv[0],"%s",dir);
  dir[strlen(dir)-6]='\0';
  if(cam){
    sprintf(fname,"%scal_H.dat",dir);
    aw[0]=665.206311;day[0]=-8.691166;
    aw[1]=665.206336;day[1]=9.427719;
    aw[2]=644.471912;day[2]=-7.485728;
    aw[3]=644.463510;day[3]=9.549897;
  }
  else{
    sprintf(fname,"%scal_C.dat",dir);
    aw[0]=385.989442;day[0]=-7.921143;
    aw[1]=385.991504;day[1]=8.662075;
    aw[2]=404.579419;day[2]=-7.926201;
    aw[3]=404.580738;day[3]=8.523750;
  }
  fprintf(stderr,"Read: %s ... ",fname);
  if((fpt=fopen(fname,"r"))==NULL){
    fprintf(stderr,"Can not open file: %s\n",fname);
    exit(0);
  }
  if(fgets(buf,sizeof(buf),fpt)==NULL){
    fprintf(stderr,"No data !\n");
    exit(1);
  }
  if(sscanf(buf,"%lf %lf %d",&wc,&dw,&icen)!=3){
    fprintf(stderr,"Broken file ?\n");
    exit(1);
  }
  nc=0;
  while(!feof(fpt)){
    if(fgets(buf,sizeof(buf),fpt)==NULL) break;
    if(buf[0]=='#') continue;
    n=sscanf(buf,"%lf %lf %s",&cx[nc][NFIBER],&cw[nc],lname);
    //fprintf(stderr,"%d: %s",n,buf);
    if(n!=3) continue;
    if(argc==2){if(strcmp(lname,"Sky")!=0) continue;}
    else{
      if(cam){if(strcmp(lname,"Ne")!=0) continue;}
      else if(strcmp(lname,"Fe")!=0) continue;
    }
    cx[nc][NFIBER]-=1.;
    //fprintf(stderr,"%d: %lf %lf\n",nc,cx[nc][NFIBER],cw[nc]);
    nc++;
    if(nc==NCAL){
      fprintf(stderr,"NCAL is too small.\n");
      exit(0);
    }
  }
  fclose(fpt);
  fprintf(stderr,"OK\n");
  fprintf(stderr,"Read: %s ... ",argv[1]);
  status=0;
  if(fits_open_file(&fp, argv[1], READONLY, &status)) printerror(status);
  if(fits_read_keys_lng(fp, "NAXIS", 1, 2, naxes, &nfound, &status)) printerror(status);
  fprintf(stderr,"%ld x %ld format ",naxes[0],naxes[1]);
  if(((naxes[0]!=LINEPIX)||(naxes[1]!=LINEPIY))&&((naxes[0]!=LINEPIX2)||(naxes[1]!=LINEPIY2))){
    fprintf(stderr,"\nIllegal image format\n");
    exit(0);
  }
  fpixel = 1;
  nullval = 0;
  for(i=0;i<naxes[1];i++){
    if (fits_read_img(fp, TFLOAT, fpixel, naxes[0], &nullval, zi[i], &anynull, &status)) printerror(status);
    fpixel+=naxes[0];
  }
  if(fits_close_file(fp, &status)) printerror(status);
  fprintf(stderr,"OK\n");
  th=0.;
  for(i=0;i<naxes[1];i++){
    ztb[i]=0.;
    for(j=0;j<naxes[0];j++) ztb[i]+=zi[i][j];
    ztb[i]/=naxes[0];
    th+=ztb[i];
  }
  th/=naxes[1];
  for(i=0;i<naxes[1];i++){
    if(ztb[i]>th) break;
  }
  it=i;
  for(i=naxes[1]-1;i>=0;i--){
    if(ztb[i]>th) break;
  }
  icen=(i+it)/2;
  fprintf(stderr,"Write: tmp.dat\n");
  fpt=fopen("tmp.dat","w");
  fprintf(fpt,"%g %g %d\n",wc,dw,icen);
  fclose(fpt);
  sprintf(buf,"tail -n+2 %s >> tmp.dat && mv tmp.dat %s",fname,fname);
  fprintf(stderr,"%s\n",buf);
  system(buf);
  resize(zi,naxes,icen);
  for(j=0;j<naxes[0];j++){
    zs=0.;
    for(i=0;i<naxes[1];i++) zs+=zi[i][j];
    zs/=naxes[1];
    for(i=0;i<naxes[1];i++) zi[i][j]*=1000./zs;
  }
  for(i=0;i<naxes[1];i++){
    for(j=0;j<naxes[0];j++) zx[i][j]=0.;
  }
  j=naxes[0]/2;
  nf=0;
  sigs=0.;
  for(i=4;i<naxes[1]-4;i++){
    if(zi[i][j]-(zi[i-4][j]+zi[i+4][j])/2>THRE){
      if(zi[i-2][j]<zi[i-1][j]&&zi[i-1][j]<zi[i][j]&&zi[i+2][j]<zi[i+1][j]&&zi[i+1][j]<zi[i][j]){
	//fprintf(stderr,"%d %d: %g %g %g %g %g\n",nf,i,zi[i-2][j],zi[i-1][j],zi[i][j],zi[i+1][j],zi[i+2][j]);
	if(nf==NFIBER){
	  fprintf(stderr,"Extra fiber? @(%d,%d)\n",j+1,i+1);
	  exit(1);
	}
	for(ii=0;ii<5;ii++) zz[ii]=zi[i-2+ii][j];
	lsfit(zz,0,5,2,10.,a);
	yc[nf][j]=i-a[1]/2/a[2];
	//if(nf) fprintf(stderr,"%d: %d %f %lf %lf\n",nf,i,zi[i][j],yc[nf][j],yc[nf][j]-yc[nf-1][j]);
	//else fprintf(stderr,"%d: %d %f %lf\n",nf,i,zi[i][j],yc[nf][j]);
	for(jj=j-1;jj>=0;jj--){
	  it=floor(yc[nf][jj+1]+0.5);
	  if(isnan(zi[it][jj])){
	    for(;jj>=0;jj--) yc[nf][jj]=0./0.;
	    continue;
	  }
	  for(ii=0;ii<5;ii++) zz[ii]=zi[it-2+ii][jj];
	  lsfit(zz,0,5,2,10.,a);
	  yc[nf][jj]=it-a[1]/2/a[2];
	  if(isnan(yc[nf][jj])){
	    yc[nf][jj]=yc[nf][jj+1];
	    continue;
	  }
	  if(a[2]>0){
	    if(yc[nf][jj]>yc[nf][j]){
	      lsfit(zz,0,4,2,10.,a);
	      yc[nf][jj]=it-a[1]/2/a[2]-0.5;
	      if(a[2]>0){
		lsfit(zz,0,3,2,10.,a);
		yc[nf][jj]=it-a[1]/2/a[2]-1;
		if(a[2]>0){
		  //fprintf(stderr,"a[2]>0 error @(%d,%d)\n",jj+1,it+1);
		  yc[nf][jj]=yc[nf][jj+1];
		}
	      }
	    }
	    else{
	      lsfit(zz,1,5,2,10.,a);
	      yc[nf][jj]=it-a[1]/2/a[2]+0.5;
	      if(a[2]>0){
		lsfit(zz,2,5,2,10.,a);
		yc[nf][jj]=it-a[1]/2/a[2]+1;
		if(a[2]>0){
		  //fprintf(stderr,"a[2]>0 error @(%d,%d)\n",jj+1,it+1);
		  yc[nf][jj]=yc[nf][jj+1];
		}
	      }
	    }	    
	  }
	  if(fabs(yc[nf][jj]-yc[nf][jj+1])>0.5){
	    //fprintf(stderr,"Broken spectrum @(%d %lf)-(%d %lf)\n",jj+1,yc[nf][jj]+1,jj+2,yc[nf][jj+1]+1);
	    if(yc[nf][jj]>yc[nf][jj+1]) yc[nf][jj]=yc[nf][jj+1]+0.5;
	    else yc[nf][jj]=yc[nf][jj+1]-0.5;
	  }
	}
	for(jj=j+1;jj<naxes[0];jj++){
	  it=floor(yc[nf][jj-1]+0.5);
	  if(isnan(zi[it][jj])){
	    for(;jj<naxes[0];jj++) yc[nf][jj]=0./0.;
	    continue;
	  }
	  for(ii=0;ii<5;ii++) zz[ii]=zi[it-2+ii][jj];
	  lsfit(zz,0,5,2,10.,a);
	  yc[nf][jj]=it-a[1]/2/a[2];
	  if(isnan(yc[nf][jj])){
	    yc[nf][jj]=yc[nf][jj-1];
	    continue;
	  }
	  if(a[2]>0){
	    if(yc[nf][jj]>yc[nf][j]){
	      lsfit(zz,0,4,2,10.,a);
	      yc[nf][jj]=it-a[1]/2/a[2]-0.5;
	      if(a[2]>0){
		lsfit(zz,0,3,2,10.,a);
		yc[nf][jj]=it-a[1]/2/a[2]-1;
		if(a[2]>0){
		  //fprintf(stderr,"a[2]>0 error @(%d,%d)\n",jj+1,it+1);
		  yc[nf][jj]=yc[nf][jj+1];
		}
	      }
	    }	    
	    else{
	      lsfit(zz,1,5,2,10.,a);
	      yc[nf][jj]=it-a[1]/2/a[2]+0.5;
	      if(a[2]>0){
		lsfit(zz,2,5,2,10.,a);
		yc[nf][jj]=it-a[1]/2/a[2]+1;
		if(a[2]>0){
		  //fprintf(stderr,"a[2]>0 error @(%d,%d)\n",jj+1,it+1);
		  yc[nf][jj]=yc[nf][jj+1];
		}
	      }
	    }	    
	  }
	  if(fabs(yc[nf][jj]-yc[nf][jj-1])>0.5){
	    //fprintf(stderr,"Broken spectrum @(%d %lf)-(%d %lf)\n",jj+1,yc[nf][jj]+1,jj,yc[nf][jj-1]+1);
	    if(yc[nf][jj]>yc[nf][jj-1]) yc[nf][jj]=yc[nf][jj-1]+0.5;
	    else yc[nf][jj]=yc[nf][jj-1]-0.5;
	  }
	}
	//fprintf(stderr,"%02d: %lf %lf\n",nf,yc[nf][0]-yc[nf][j],yc[nf][naxes[0]-1]-yc[nf][j]);
	sig=lsfit(yc[nf],0,naxes[0],YOD,3.,af[nf]);
	sigs+=sig;
	//fprintf(stderr,"%d: %g %g %g %g",nf,yc[nf][naxes[0]/2],af[nf][0],af[nf][0]-yc[nf][naxes[0]/2],af[nf][0]-(LINE0+PSPEC*nf));getchar();
	af[nf][0]-=LINE0+PSPEC*nf;
	//printf("%02d",nf);
	//for(k=YOD;k>=0;k--) printf(" %9.6lf",af[nf][k]*pow(10.,k*3));
	//printf(" %lf\n",sig);
	nf++;
      }
    }
  }
  if(nf!=NFIBER){
    fprintf(stderr,"Number of fibers (%d) is not %d ...\n",nf,NFIBER);
    exit(1);
  }
  sigs/=NFIBER;
  fprintf(stderr,"Y sigma = %g\n",sigs);
  for(k=1;k<=YOD;k++){
    for(nf=0;nf<NFIBER;nf++) at[nf]=af[nf][k];
    lsfit(at,0,NFIBER,2,3.,aa[k]);
  }
  for(nf=0;nf<NFIBER;nf++){
    x=nf-(NFIBER-1)/2;
    for(k=1;k<=YOD;k++) af[nf][k]=(aa[k][2]*x+aa[k][1])*x+aa[k][0];
    for(j=0;j<naxes[0];j++){
      x=j-(naxes[0]-1)/2;
      yc[nf][j]=af[nf][YOD];
      for(k=YOD-1;k>=0;k--) yc[nf][j]=yc[nf][j]*x+af[nf][k];
      //fprintf(stderr,"%4d %g",j,yc[nf][j]);getchar();
    }
    //printf("%02d",nf);
    //for(k=YOD;k>=0;k--) printf(" %9.6lf",af[nf][k]*pow(10.,k*3));
    //printf(" %lf\n",sig);
  }
  at2s0=1e10;
  while(1){
    for(j=0;j<naxes[0];j++){
      //dy=(yc[0][j]-(yc[1][j]-yc[0][j])*2)/(LINE0+yc[0][j]-(yc[1][j]-yc[0][j]+PSPEC)*2);
      //for(i=0;i<LINE0+yc[0][j]-(yc[1][j]-yc[0][j]+PSPEC)*2;i++) zy[i][j]+=-dy*i;//edge 0
      //for(i=0;i<LINE0+yc[0][j]-(yc[1][j]-yc[0][j]+PSPEC)*2;i++) zy[i][j]+=-yc[0][j]+(yc[1][j]-yc[0][j])*2;//const
      i=0;
      for(nf=0;nf<NFIBER-1;nf++){
	dy=(yc[nf+1][j]-yc[nf][j])/(yc[nf+1][j]-yc[nf][j]+PSPEC);
	for(;i<LINE0+PSPEC*(nf+1)+yc[nf+1][j];i++) zy[i][j]=-yc[nf][j]-dy*(i-(LINE0+PSPEC*nf+yc[nf][j]));
      }
      //for(;i<LINE0+PSPEC*(NFIBER-1)+yc[NFIBER-1][j]+(yc[NFIBER-1][j]-yc[NFIBER-2][j]+PSPEC)*2;i++) zy[i][j]=-yc[NFIBER-1][j]-dy*(i-(LINE0+PSPEC*(NFIBER-1)+yc[NFIBER-1][j]));
      for(;i<naxes[1];i++) zy[i][j]=-yc[NFIBER-1][j]-dy*(i-(LINE0+PSPEC*(NFIBER-1)+yc[NFIBER-1][j]));
      //dy=(yc[NFIBER-1][j]+(yc[NFIBER-1][j]-yc[NFIBER-2][j])*2)/(naxes[1]-1-(LINE0+PSPEC*(NFIBER-1)+yc[NFIBER-1][j]+(yc[NFIBER-1][j]-yc[NFIBER-2][j]+PSPEC)*2));
      //for(;i<naxes[1];i++) zy[i][j]=-dy*(naxes[1]-1-i);//edge 0
      //for(;i<naxes[1];i++) zy[i][j]=-yc[NFIBER-1][j]-(yc[NFIBER-1][j]-yc[NFIBER-2][j])*2;//const
      //for(i=0;i<naxes[1];i++){fprintf(stderr,"%d: %g",i,zy[i][j]);getchar();}
    }
    warp(zi,zo,zx,zy);
    for(i=0;i<naxes[1];i++){
      yav[i]=0.;n=0;
      for(j=1;j<naxes[0]-1;j++){
	if(isnan(zo[i][j-1])==0&&isnan(zo[i][j])==0&&isnan(zo[i][j+1])==0){
	  yav[i]+=zo[i][j];
	  n++;
	}
      }
      yav[i]/=n;
    }
    at2s=0.;
    for(nf=0;nf<NFIBER;nf++){
      i=LINE0+PSPEC*nf;
      for(ii=0;ii<5;ii++) zz[ii]=yav[i-2+ii];
      lsfit(zz,0,5,2,10.,a);
      at[nf]=a[1]/2/a[2];
      at2s+=at[nf]*at[nf];
      //fprintf(stderr,"%d: %g %g\n",nf,at[nf],yc[nf][naxes[0]/2]);
    }
    at2s=sqrt(at2s/NFIBER);
    fprintf(stderr,"Yave sigma = %g\r",at2s);
    if(at2s0>at2s+0.001){
      at2s0=at2s;
      for(nf=0;nf<NFIBER;nf++){
	for(j=0;j<naxes[0];j++) yc[nf][j]-=at[nf];
      }
    }
    else{
      fprintf(stderr,"\n");
      break;
    }
  }
  if(argc==3){
    cl=2;sn=1;
    for(i=0;i<naxes[1];i++){
      for(j=0;j<naxes[0];j++) zo[i][j]=zi[i][j];
    }
  }
  else{
    cl=1;sn=-1;
    for(i=0;i<naxes[1];i++){
      for(j=0;j<naxes[0];j++) zo[i][j]=0.;
    }
  }
  fprintf(stderr,"Read: %s ... ",argv[cl]);
  status=0;
  if(fits_open_file(&fp, argv[cl], READONLY, &status)) printerror(status);
  if(fits_read_keys_lng(fp, "NAXIS", 1, 2, naxes, &nfound, &status)) printerror(status);
  fprintf(stderr,"%ld x %ld format ",naxes[0],naxes[1]);
  if(((naxes[0]!=LINEPIX)||(naxes[1]!=LINEPIY))&&((naxes[0]!=LINEPIX2)||(naxes[1]!=LINEPIY2))){
    fprintf(stderr,"\nIllegal image format\n");
    exit(0);
  }
  fpixel = 1;
  nullval = 0;
  for(i=0;i<naxes[1];i++){
    if (fits_read_img(fp, TFLOAT, fpixel, naxes[0], &nullval, zi[i], &anynull, &status)) printerror(status);
    fpixel+=naxes[0];
  }
  if(fits_close_file(fp, &status)) printerror(status);
  fprintf(stderr,"OK\n");
  resize(zi,naxes,icen);
  for(i=0;i<naxes[1];i++){
    for(j=0;j<naxes[0];j++) zo[i][j]+=zi[i][j];
  }
  wrp(zi,zx,zy);
  for(nf=0;nf<NFIBER;nf++){
    i=LINE0+PSPEC*nf-PSPEC/2;
    for(j=0;j<naxes[0];j++){
      zt[j]=0.;
      for(ii=0;ii<PSPEC;ii++) zt[j]+=zi[i+ii][j]*kn[ii];
    }
    if(argc==3){
      for(j=0;j<naxes[0];j++) ztb[j]=1.;
    }
    else{
      for(j=0;j<naxes[0];j++){
	if(j<50 || j>naxes[0]-1-50){ztb[j]=0./0.;continue;}
	ztb[j]=0.;
	for(jj=j-50;jj<=j+50;jj++) ztb[j]+=zt[jj];
	ztb[j]/=100*1000;
      }
    }
    th=THRE;nd=0;
    while(nd<nc||NCAL<=nd){
      nd=0;
      for(j=4;j<naxes[0]-4;j++){
	//if(j==2843) fprintf(stderr,"%g %g %g %g %g %g %g %g %g\n",zt[j-4],zt[j-3],zt[j-2],zt[j-1],zt[j],zt[j+1],zt[j+2],zt[j+3],zt[j+4]);
	if((zt[j]-(zt[j-4]+zt[j+4])/2)*sn>th*ztb[j]){
	  if(zt[j-2]*sn<zt[j-1]*sn&&zt[j-1]*sn<zt[j]*sn&&zt[j+2]*sn<zt[j+1]*sn&&zt[j+1]*sn<zt[j]*sn) nd++;
	}
      }
      if(nd>=NCAL) th*=1.1;
      else if(nd<nc) th*=0.9;
      //fprintf(stderr,"%d: %d %g\n",nf,nd,th);
    }//getchar();
    //fprintf(stderr,"%d: %d %g\n",nf,nd,th);
    nd=0;
    for(j=4;j<naxes[0]-4;j++){
      if((zt[j]-(zt[j-4]+zt[j+4])/2)*sn>th*ztb[j]){
	if(zt[j-2]*sn<zt[j-1]*sn&&zt[j-1]*sn<zt[j]*sn&&zt[j+2]*sn<zt[j+1]*sn&&zt[j+1]*sn<zt[j]*sn){
	  for(jj=0;jj<5;jj++) zz[jj]=zt[j-2+jj]*sn;
	  lsfit(zz,0,5,2,10.,a);
	  cxt[nd]=j-a[1]/2/a[2];
	  //fprintf(stderr,"%d %d: %g %g %g %g %g %g %g\n",nd,j+1,zz[0],zz[1],zz[2],zz[3],zz[4],zt[j]-(zt[j-4]+zt[j+4])/2,cxt[nd]+1);
	  //if(nf==NFIBER/2) fprintf(stderr,"%g\n",cxt[nd]+1);
	  //if(nf==NFIBER/2) fprintf(stderr,"%lf %lf\n",cxt[nd]+1,(cxt[nd]-LINEPIX2/2)*dw+wc);
	  nd++;
	}
      }
    }
    mmax=0;
    for(of=-100;of<=100;of++){
      m=0;
      for(n=0;n<nc;n++){
	for(nn=0;nn<nd;nn++){
	  if(fabs(cxt[nn]-cx[n][NFIBER]-of)<2.){
	    m++;
	    break;
	  }
	}
      }
      if(mmax<m){
	mmax=m;
	ofmax=of;
	if(mmax==nc) break;
      }
    }
    if(nf==NFIBER/2) fprintf(stderr,"Number of detected lines: %d/%d\n",mmax,nc);
    for(n=0;n<nc;n++){
      for(nn=0;nn<nd;nn++){
	if(fabs(cxt[nn]-cx[n][NFIBER]-ofmax)<2.){
	  cx[n][nf]=cxt[nn];
	  break;
	}
      }
      if(nn==nd) cx[n][nf]=0./0.;
    }
  }
  //for(nf=0;nf<NFIBER;nf++){
  //fprintf(stderr,"%2d:",nf+1);
  //for(n=0;n<nc;n++){
  //fprintf(stderr," %7.2lf",cx[n][nf]);
  //}
  //fprintf(stderr,"\n");
  //}
  sigs=0.;
  for(nf=0;nf<NFIBER;nf++){
    for(n=0;n<nc;n++) cxt[n]=cx[n][nf]-naxes[0]/2;
    sig=lsft(cxt,cw,0,nc,XOD,10.,af[nf]);
    sigs+=sig/af[nf][1];
    //printf("%02d",nf);
    //for(k=XOD;k>=0;k--) printf(" %9.6lf",af[nf][k]*pow(10.,k*3));
    //printf(" %lf\n",sig);
  }
  sigs/=NFIBER;
  fprintf(stderr,"X sigma = %g\n",sigs);
  fprintf(stderr,"Write: mzwarp.dat ... ");
  fpt=fopen("mzwarp.dat","w");    
  nf=NFIBER/2;
  for(n=0;n<nc;n++){
    wt=af[nf][XOD];
    dwt=af[nf][XOD]*XOD;
    for(k=XOD-1;k>=0;k--){
      wt=wt*(cx[n][nf]-naxes[0]/2)+af[nf][k];
      if(k>0) dwt=dwt*(cx[n][nf]-naxes[0]/2)+af[nf][k]*k;
    }
    fprintf(fpt,"%lf %lf %lf %+lf %lf\n",cx[n][nf]+1,cw[n],wt,(wt-cw[n])/dwt,dwt);
  }
  fclose(fpt);
  fprintf(stderr,"OK\n");
  sprintf(buf,"gnuplot < %smzwarp.plt\n",dir);
  fprintf(stderr,"%s",buf);
  system(buf);
  if(cam)sprintf(buf,"convert -density 65 -rotate 90 mzwarp.ps /mnt/QL/inst2.png 2>/dev/null\n");
  else sprintf(buf,"convert -density 65 -rotate 90 mzwarp.ps /ssd/QL/inst.png 2>/dev/null\n");
  //fprintf(stderr,"%s",buf);
  system(buf);  
  
  //nf=60;
  //for(j=0;j<naxes[0];j++){
  //wt=af[nf][XOD];
  //for(k=XOD-1;k>=0;k--) wt=wt*(j-naxes[0]/2)+af[nf][k];
  //zt[j]=wt;
  //}
  for(k=1;k<=XOD;k++){
    for(nf=0;nf<NFIBER;nf++) at[nf]=af[nf][k];
    lsfit(at,0,NFIBER,1,2.,aa[k]);
    //lsfit(at,0,NFIBER,0,2.,aa[k]);aa[k][1]=0.;
  }
  for(nf=0;nf<NFIBER;nf++){
    x=nf-(NFIBER-1)/2;
    for(k=1;k<=XOD;k++) af[nf][k]=aa[k][1]*x+aa[k][0];
    //fprintf(stderr,"%2d:",nf+1);
    //for(n=0;n<nc;n++){
    //wt=af[nf][XOD];
    //for(k=XOD-1;k>=0;k--) wt=wt*(cx[n][nf]-naxes[0]/2)+af[nf][k];
    //fprintf(stderr," %7.3lf",wt);
    //fprintf(stderr," %7.2lf",(wt-wc)/dw+LINEPIX2/2);
    //}
    //fprintf(stderr,"\n");
    //printf("%02d",nf);
    //for(k=XOD;k>=0;k--) printf(" %9.6lf",af[nf][k]*pow(10.,k*3));
    //printf(" %lf\n",sig);
    for(j=0;j<naxes[0];j++){
      wt=af[nf][XOD];
      for(k=XOD-1;k>=0;k--) wt=wt*(j-naxes[0]/2)+af[nf][k];
      //if(nf==60) printf("%d %g %g\n",j,zt[j],wt);
      xc[nf][j]=j-((wt-wc)/dw+LINEPIX2/2);
    }
  }
  if(sn>0){
    for(j=0;j<naxes[0];j++){
      for(i=0;i<LINE0+yc[0][j];i++) zx[i][j]=-xc[0][j];
      for(nf=0;nf<NFIBER-1;nf++){
	dx=(xc[nf+1][j]-xc[nf][j])/(yc[nf+1][j]-yc[nf][j]+PSPEC);
	for(;i<LINE0+PSPEC*(nf+1)+yc[nf+1][j];i++) zx[i][j]=-xc[nf][j]-dx*(i-(LINE0+PSPEC*nf+yc[nf][j]));
      }
      for(;i<naxes[1];i++) zx[i][j]=-xc[NFIBER-1][j];
    }
  }
  else{
    at2s0=1e10;
    while(1){
      for(j=0;j<naxes[0];j++){
	for(i=0;i<LINE0+yc[0][j];i++) zx[i][j]=-xc[0][j];
	for(nf=0;nf<NFIBER-1;nf++){
	  dx=(xc[nf+1][j]-xc[nf][j])/(yc[nf+1][j]-yc[nf][j]+PSPEC);
	  for(;i<LINE0+PSPEC*(nf+1)+yc[nf+1][j];i++) zx[i][j]=-xc[nf][j]-dx*(i-(LINE0+PSPEC*nf+yc[nf][j]));
	}
	for(;i<naxes[1];i++) zx[i][j]=-xc[NFIBER-1][j];
      }
      warp(zo,zi,zx,zy);
      for(nf=0;nf<NFIBER;nf++) at[nf]=0.;
      nd=0;
      for(n=0;n<nc;n++){
	for(nf=0;nf<NFIBER;nf++){
	  if(isnan(cx[n][nf])) break;
	}
	if(nf<NFIBER) continue;
	for(nf=0;nf<NFIBER;nf++){
	  i=LINE0+PSPEC*nf-PSPEC/2;
	  j=floor((cw[n]-wc)/dw+0.5)+LINEPIX2/2;
	  for(jj=0;jj<5;jj++){
	    zz[jj]=0.;
	    for(ii=0;ii<PSPEC;ii++) zz[jj]+=zi[i+ii][j-2+jj]*kn[ii]*sn;
	  }
	  lsfit(zz,0,5,2,10.,a);
	  at[nf]+=a[1]/2/a[2];
	}
	nd++;
      }
      at2s=0.;
      for(nf=0;nf<NFIBER;nf++){
	at[nf]/=nd;
	at2s+=at[nf]*at[nf];
	//fprintf(stderr,"%d: %g %g\n",nf,at[nf],xc[nf][naxes[0]/2]);
      }
      at2s=sqrt(at2s/NFIBER);
      fprintf(stderr,"Xave sigma = %g\r",at2s);
      if(at2s0>at2s+0.001){
	at2s0=at2s;	
	for(nf=0;nf<NFIBER;nf++){
	  for(j=0;j<naxes[0];j++) xc[nf][j]-=at[nf];
	}
      }
      else{
	fprintf(stderr,"\n");
	break;
      }
    }
  }

  for(i=0;i<naxes[1];i++){
    for(j=naxes[0]/2-1;j>=0;j--){
      if(isnan(zi[i][j])){
        zx[i][j]=zx[i][j+1];zy[i][j]=zy[i][j+1];
      }
    }
    for(j=naxes[0]/2+1;j<naxes[0];j++){
      if(isnan(zi[i][j])){
        zx[i][j]=zx[i][j-1];zy[i][j]=zy[i][j-1];
      }
    }
  }
  //double bx[4]={854.63,854.62,4219.47,4219.67},by[4]={116.81,599.47,116.47,597.69};
  //double bx[4]={689.13,691.18,4051.42,4053.04},by[4]={88.34,630.27,88.77,629.92};
  for(n=0;n<2;n++){
    nf=(NFIBER-1)*n;
    x=(aw[n]-wc)/dw+LINEPIX2/2;
    j=LINEPIX2/2;
    while(j-xc[nf][j]>x) j--;
    dx=(x-(j-xc[nf][j]))/(1-xc[nf][j+1]+xc[nf][j]);
    ax[n]=j+dx;
    ay[n]=LINE0+PSPEC*nf+yc[nf][j]*(1-dx)+yc[nf][j+1]*dx+day[n];
    fprintf(stderr,"Anchor%d: %lf %lf\n",n,ax[n]+1,ay[n]+1);
    //fprintf(stderr,"Anchor%d: %lf %lf\n",n,aw[n]-(ax[n]+1-bx[n])*dw,by[n]-(ay[n]+1-day[n]));
  }
  for(n=0;n<2;n++){
    nf=(NFIBER-1)*n;
    x=(aw[n+2]-wc)/dw+LINEPIX2/2;
    j=LINEPIX2/2;
    while(j-xc[nf][j]<x) j++;
    dx=((j-xc[nf][j])-x)/(1-xc[nf][j]+xc[nf][j-1]);
    ax[n+2]=j-dx;
    ay[n+2]=LINE0+PSPEC*nf+yc[nf][j]*(1-dx)+yc[nf][j-1]*dx+day[n+2];
    fprintf(stderr,"Anchor%d: %lf %lf\n",n+2,ax[n+2]+1,ay[n+2]+1);
    //fprintf(stderr,"Anchor%d: %lf %lf\n",n+2,aw[n+2]-(ax[n+2]+1-bx[n+2])*dw,by[n+2]-(ay[n+2]+1-day[n+2]));
  }
  wrp(zo,zx,zy);
  remove("mzwarp.fits");
  fprintf(stderr,"Write mzwarp.fits ... ");
  status=0;
  if (fits_create_file(&fp, "mzwarp.fits", &status)) printerror(status);
  if (fits_create_img(fp, bitpix, naxis, naxes, &status)) printerror(status);
  if (fits_write_key(fp, TSTRING, "CTYPE1", "WAVELENGTH", "", &status)) printerror(status);
  if (fits_write_key(fp, TSTRING, "CUNIT1", "nm", "", &status)) printerror(status);
  tmpval=naxes[0]/2+1.;
  if (fits_write_key(fp, TFLOAT, "CRPIX1", &tmpval, "", &status)) printerror(status);
  tmpval=wc;
  if (fits_write_key(fp, TFLOAT, "CRVAL1", &tmpval, "", &status)) printerror(status);
  tmpval=dw;
  if (fits_write_key(fp, TFLOAT, "CDELT1", &tmpval, "", &status)) printerror(status);
  if (fits_write_key(fp, TFLOAT, "CD1_1", &tmpval, "", &status)) printerror(status);
  tmpval=1.;
  if (fits_write_key(fp, TFLOAT, "CD2_2", &tmpval, "", &status)) printerror(status);
  fpixel=1;
  fprintf(stderr,"%ld x %ld format ",naxes[0],naxes[1]);
  for(i=0;i<naxes[1];i++){
    if(fits_write_img(fp, TFLOAT, fpixel, naxes[0], zo[i], &status)) printerror(status);
    fpixel+=naxes[0];
  }
  if (fits_close_file(fp, &status)) printerror(status);
  fprintf(stderr,"OK\n");
  naxis=3;
  remove("mzwarp_dxdy.fits");
  fprintf(stderr,"Write mzwarp_dxdy.fits ... ");
  status=0;
  if (fits_create_file(&fp, "mzwarp_dxdy.fits", &status)) printerror(status);
  if (fits_create_img(fp, bitpix, naxis, naxes, &status)) printerror(status);
  if (fits_write_key(fp, TSTRING, "CTYPE1", "WAVELENGTH", "", &status)) printerror(status);
  if (fits_write_key(fp, TSTRING, "CUNIT1", "nm", "", &status)) printerror(status);
  tmpval=naxes[0]/2+1.;
  if (fits_write_key(fp, TFLOAT, "CRPIX1", &tmpval, "", &status)) printerror(status);
  tmpval=wc;
  if (fits_write_key(fp, TFLOAT, "CRVAL1", &tmpval, "", &status)) printerror(status);
  tmpval=dw;
  if (fits_write_key(fp, TFLOAT, "CDELT1", &tmpval, "", &status)) printerror(status);
  if (fits_write_key(fp, TFLOAT, "CD1_1", &tmpval, "", &status)) printerror(status);
  tmpval=1.;
  if (fits_write_key(fp, TFLOAT, "CD2_2", &tmpval, "", &status)) printerror(status);
  tmpval=1.;
  if (fits_write_key(fp, TFLOAT, "CD3_3", &tmpval, "", &status)) printerror(status);
  tmpval=ax[0]+1;
  if (fits_write_key(fp, TFLOAT, "MZ_ANCX1", &tmpval, "", &status)) printerror(status);
  tmpval=ay[0]+1;
  if (fits_write_key(fp, TFLOAT, "MZ_ANCY1", &tmpval, "", &status)) printerror(status);
  tmpval=ax[1]+1;
  if (fits_write_key(fp, TFLOAT, "MZ_ANCX2", &tmpval, "", &status)) printerror(status);
  tmpval=ay[1]+1;
  if (fits_write_key(fp, TFLOAT, "MZ_ANCY2", &tmpval, "", &status)) printerror(status);
  tmpval=ax[2]+1;
  if (fits_write_key(fp, TFLOAT, "MZ_ANCX3", &tmpval, "", &status)) printerror(status);
  tmpval=ay[2]+1;
  if (fits_write_key(fp, TFLOAT, "MZ_ANCY3", &tmpval, "", &status)) printerror(status);
  tmpval=ax[3]+1;
  if (fits_write_key(fp, TFLOAT, "MZ_ANCX4", &tmpval, "", &status)) printerror(status);
  tmpval=ay[3]+1;
  if (fits_write_key(fp, TFLOAT, "MZ_ANCY4", &tmpval, "", &status)) printerror(status);
  fpixel=1;
  fprintf(stderr,"%ld x %ld x %ld format ",naxes[0],naxes[1],naxes[2]);
  for(i=0;i<naxes[1];i++){
    if(fits_write_img(fp, TFLOAT, fpixel, naxes[0], zx[i], &status)) printerror(status);
    fpixel+=naxes[0];
  }
  for(i=0;i<naxes[1];i++){
    if(fits_write_img(fp, TFLOAT, fpixel, naxes[0], zy[i], &status)) printerror(status);
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
