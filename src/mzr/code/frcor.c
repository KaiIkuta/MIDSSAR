#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define LINEPIX2 4848
#define PER 14
#define TK 1.113
#define AOI 37.4
#define MAXODR 4
#define THRE 0.1
#define BUFF 28

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
      for(i=0;i<n;i++){
	for(j=0;j<n;j++) fprintf(stderr," %lf",p[i][j]);
	fprintf(stderr,"\n");
      }
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

int main(int argc,char *argv[])
{
  FILE *fp;
  char fname[64],buf[512],calfile[64];
  int n,ndat,j,jmax,k,rv[384],i;
  double wl[2][LINEPIX2],fl[2][LINEPIX2],fs[LINEPIX2],fw[LINEPIX2],ff[LINEPIX2],no[LINEPIX2],flc[2][LINEPIX2],pk[384],fws,wc,dw,wt,per,d,ph;
  double f[MAXODR][LINEPIX2],p[MAXODR][MAXODR],q[MAXODR][MAXODR],r[MAXODR],a[MAXODR];

  if(argc!=2&&argc!=3){
    fprintf(stderr,"Usage : %s finename",argv[0]);
    fprintf(stderr,"             or phase plot mode\n");
    fprintf(stderr,"             %s finename -",argv[0]);
    exit(0);
  }
  sscanf(argv[0],"%s",calfile);
  sprintf(calfile+strlen(calfile)-5,"cal_H.dat");
  sscanf(argv[1],"%s",fname);
  fprintf(stderr,"Read: %s ... ",calfile);
  if((fp = fopen(calfile,"r"))== NULL){
    fprintf(stderr,"Cannot open.\n");
    exit(1);
  }
  if(fscanf(fp,"%lf %lf",&wc,&dw)!=2){
    fprintf(stderr,"Broken file ?\n");
    exit(1);
  }
  fclose(fp);
  fprintf(stderr,"wc=%g dw=%g OK\n",wc,dw);
  fprintf(stderr,"Read: %s ... ",fname);
  if((fp = fopen(fname,"r"))== NULL){
    fprintf(stderr,"Cannot open.\n");
    exit(1);
  }
  n=0;
  while(!feof(fp)){
    if(fgets(buf,sizeof(buf),fp)==NULL) break;
    if(buf[0]=='#') continue;
    if(sscanf(buf,"%lf %lf %lf %lf %lf",&wl[0][n],&wl[1][n],&fl[0][n],&fl[1][n],&no[n])!=5){
      fprintf(stderr,"Truncated data: %s",buf);
      exit(1);
    }
    n++;
  }
  ndat=n;
  fclose(fp);
  fprintf(stderr,"ndat=%d OK\n",ndat);
  for(k=0;k<2;k++){
    for(n=0;n<ndat;n++) fs[n]=0.;
    for(n=PER;n<ndat-PER;n++){
      for(j=-PER;j<=PER;j++) fs[n]+=fl[k][n+j];
      fs[n]/=2*PER+1;
    }
    for(n=0;n<ndat;n++){
      fw[n]=fl[k][n]/fs[n]-1;
      //printf("%lf %lf %lf\n",wl[0][n],fl[k][n],fw[n]);
      if(fabs(fw[n])>THRE){
	for(j=0;j<=BUFF;j++) fw[n-j]=0./0.;
      }
    }
    for(n=ndat-1;n>=0;n--){
      if(isnan(fw[n])){
	for(j=1;j<=BUFF;j++) fw[n+j]=0./0.;
      }
    }
    d=2*TK*1e6*cos(AOI/180*M_PI)*1.456;//Fused Silica
    if(argc==3 && k==0){
      j=-2;
      for(n=1;n<ndat;n++){
	if(fw[n-1]*fw[n]<0){
	  j++;
	  if(j<0) continue;
	  if(j==0) rv[j]=n;
	  else{
	    if(n-rv[j-1]>0){//PER-3){
	      rv[j]=n;
	      //fprintf(stderr,"%d %d %d\n",j,n,rv[j]);
	    }
	    else j--;
	  }
	}
      }
      jmax=j-1;
      fprintf(stderr,"jmax=%d\n",jmax);
      for(j=0;j<jmax;j++){
	//if(j==101) fprintf(stderr,"%d %d %d\n",j,rv[j],rv[j+1]);
	pk[j]=0.;fws=0.;
	for(n=rv[j];n<rv[j+1];n++){
	  pk[j]+=n*fw[n];
	  fws+=fw[n];
	  //if(j==101) fprintf(stderr,"%d %lf %lf %lf\n",n,fw[n],pk[j],fws);
	}
	pk[j]/=fws;
      }
      //fprintf(stderr,"pk[0]=%lf\n",pk[0]);
      fprintf(stderr,"Write: frcor.dat ... ");
      fp = fopen("frcor.dat","w");
      for(j=10;j<jmax-10;j++){
	wt=wc+(pk[j]-ndat/2)*dw;
	per=wt*wt/d;
	//fprintf(fp,"%d %lf %lf %lf\n",j,pk[j],pk[j+1]-pk[j-1],(pk[j+10]-pk[j-10])/10);
	fprintf(fp,"%lf %lf %lf %lf\n",wt,(pk[j-1]-pk[j+1])*dw,(pk[j-10]-pk[j+10])/10*dw,per);
      }
      fclose(fp);
      fprintf(stderr,"OK\n");
    }
    for(i=0;i<MAXODR;i++){
      r[i]=0.;
      for(j=0;j<MAXODR;j++) p[i][j]=0.;
    }
    for(n=0;n<ndat;n++){
      ph=2*M_PI*d*(1/wc-1/wl[0][n]);
      f[0][n]=sin(ph);
      f[1][n]=cos(ph);
      f[2][n]=(wl[0][n]-wc)/10*f[0][n];
      f[3][n]=(wl[0][n]-wc)/10*f[1][n];
    }
    j=0;
    for(n=0;n<ndat;n++){
      if(isnan(fw[n])) continue;
      p[0][0]+=f[0][n]*f[0][n];
      p[0][1]+=f[0][n]*f[1][n];
      p[0][2]+=f[0][n]*f[2][n];
      p[0][3]+=f[0][n]*f[3][n];
      p[1][1]+=f[1][n]*f[1][n];
      p[1][2]+=f[1][n]*f[2][n];
      p[1][3]+=f[1][n]*f[3][n];
      p[2][2]+=f[2][n]*f[2][n];
      p[2][3]+=f[2][n]*f[3][n];
      p[3][3]+=f[3][n]*f[3][n];
      r[0]+=fw[n]*f[0][n];
      r[1]+=fw[n]*f[1][n];
      r[2]+=fw[n]*f[2][n];
      r[3]+=fw[n]*f[3][n];
      j++;
    }
    p[0][0]/=j;
    p[0][1]/=j;
    p[0][2]/=j;
    p[0][3]/=j;
    p[1][0]=p[0][1];
    p[1][1]/=j;
    p[1][2]/=j;
    p[1][3]/=j;
    p[2][0]=p[0][2];
    p[2][1]=p[1][2];
    p[2][2]/=j;
    p[2][3]/=j;
    p[3][0]=p[0][3];
    p[3][1]=p[1][3];
    p[3][2]=p[2][3];
    p[3][3]/=j;
    r[0]/=j;
    r[1]/=j;
    r[2]/=j;
    r[3]/=j;
    sweep(p,q,MAXODR);
    for(i=0;i<MAXODR;i++){
      a[i]=0.;
      for(j=0;j<MAXODR;j++){
	a[i]+=q[i][j]*r[j];
      }
    }
    for(n=0;n<ndat;n++){
      ff[n]=a[0]*f[0][n]+a[1]*f[1][n]+a[2]*f[2][n]+a[3]*f[3][n];
      flc[k][n]=fl[k][n]-ff[n]*fs[n];
      //printf("%lf %lf %lf %lf %lf %lf %lf\n",wl[0][n],fl[k][n],fw[n],fs[n]*(1+ff[n]),ff[n],flc[k][n],fw[n]-ff[n]);    
    }
  }
  for(n=0;n<ndat;n++) printf("%lf %lf %lf %lf %lf\n",wl[0][n],wl[1][n],flc[0][n],flc[1][n],no[n]);
  return 0;
}
