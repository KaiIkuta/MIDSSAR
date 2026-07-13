#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <dirent.h>
#include <sys/stat.h>
#include <time.h>
#include "fitsio.h"
#include "warp.c"
#define JMAX 40
#define IMAX 40
#define BMPZ1 -1//0.9
#define BMPZ2 10//1.1
#define SKY 0.66666
#define MAXODR 5
#define SATU 60000
#define CONVF 0.083
#define MW 40
#define NCAL 80
#define HWZI 300

char dir[32]=".";

void printerror( int status);

float ave(float z[LINEPIY2][LINEPIX2],int jc,int ic){
  int i,j,is,ie,n,lp;
  double zs,z2s,av=0.,sg=1e99;

  if(ic){
    is=1;
    ie=LINE0-PSPEC/2;
  }
  else{
    is=LINE0+PSPEC*NFIBER-PSPEC/2;
    ie=LINEPIY2-1;
  }
  for(lp=0;lp<4;lp++){
    n=0;zs=0.;z2s=0.;
    for(j=jc-100;j<=jc+100;j++){
      for(i=is;i<ie;i++){
	if(fabs(z[i][j]-av)<sg*3){
	  zs+=z[i][j];
	  z2s+=z[i][j]*z[i][j];
	  n++;
	}
      }
    }
    av=zs/n;
    sg=sqrt(z2s/n-av*av);
    //if(jc==2424) fprintf(stderr,"%d: %g %g\n",lp,av,sg);
  }
  return (av);
}

float med(float sp[LINEPIX2],int jc){
  float m[MW+1],md;
  int j,k,kk,n;

  for(k=0;k<=MW;k++) m[k]=-1e10;
  n=-1;
  for(j=jc-MW;j<=jc+MW;j++){
    if(j<MW||j>LINEPIX2-MW-1) continue;
    if(isnan(sp[j])) continue;
    n++;
    for(k=0;k<=MW;k++){
      if(m[k]<sp[j]){
        for(kk=MW;kk>k;kk--) m[kk]=m[kk-1];
        m[k]=sp[j];
        break;
      }
    }
  }
  if(n%2) md=(m[n/2+1]+m[n/2])/2;
  else md=m[n/2];
  //if(j==2424){
  //for(k=0;k<=MW;k++) fprintf(stderr,"%d: %g %g\n",k,m[k],sp[jc+k]);
  //fprintf(stderr,"n=%d md=%g",n,md);getchar();
  //}
  return(md);
}

float med2(double f[]){
  float m[NFIBER/2+1],md;
  int j,k,kk;

  for(k=0;k<=NFIBER/2;k++) m[k]=-1e10;
  for(j=0;j<NFIBER;j++){
    for(k=0;k<=NFIBER/2;k++){
      if(m[k]<f[j]){
        for(kk=NFIBER/2;kk>k;kk--) m[kk]=m[kk-1];
        m[k]=f[j];
        break;
      }
    }
  }
  //for(k=0;k<=NFIBER/2;k++) fprintf(stderr,"%d: %g\n",k,m[k]);
  //for(j=0;j<NFIBER;j++) fprintf(stderr,"%d: %g\n",j,f[j]);
  //md=m[NFIBER/2];
  md=m[NFIBER/3];
  return(md);
}

int sweep(double p[MAXODR][MAXODR],double q[MAXODR][MAXODR],int n)
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
      return 1;
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
  return 0;
}

double lsfit(double x[LINEPIX2], double y[LINEPIX2], int od, double clip, double a[])
{
  int i,j,n,nd,lp=3;
  double sxn[MAXODR*2],xn[MAXODR*2],dy2s,sig=1e10;
  double p[MAXODR][MAXODR],q[MAXODR][MAXODR],sxny[MAXODR],f[LINEPIX2];

  xn[0]=1;
  sxn[0]=1;
  for(n=0;n<LINEPIX2;n++) f[n]=0.;
  if(clip>9) lp=1;
  for(;lp>0;lp--){
    for(i=1;i<MAXODR*2;i++) sxn[i]=0;
    for(i=0;i<MAXODR;i++) sxny[i]=0;
    nd=0.;
    for(n=0;n<LINEPIX2;n++){
      if(isnan(x[n])||isnan(y[n])) continue;
      //fprintf(stderr,"%d %g %g %g",n,x[n],y[n],f[n]);getchar();
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
    for(n=0;n<LINEPIX2;n++){
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

double lsft(float z[LINEPIX2], double y[2][LINEPIX2], double clip, float f[LINEPIX2])
{
  int i,j,n,nd,lp=3,od;
  double x[LINEPIX2],x1,x2,x3,x4,dy2s,sig=1e10;
  double sx1,sx2,sx3,sx4,sy[2],sy2[2],sxy[2],sx2y[2],syy,sz[MAXODR];
  double p[MAXODR][MAXODR],q[MAXODR][MAXODR],a[MAXODR];

  od=MAXODR-2;
  for(n=0;n<LINEPIX2;n++){
    f[n]=0.;
    x[n]=(double)(n-LINEPIX/2)/(LINEPIX/2);
    if(fabs(y[1][n])>1e-5) od=MAXODR-1;
  }
  if(clip>9) lp=1;
  for(;lp>0;lp--){
    sx1=0.;sx2=0.;sx3=0.;sx4=0.;
    sy[0]=0.;sy2[0]=0.;sxy[0]=0.;sx2y[0]=0.;
    sy[1]=0.;sy2[1]=0.;sxy[1]=0.;sx2y[1]=0.;
    syy=0.;sz[0]=0.;sz[1]=0.;sz[2]=0.;sz[3]=0.;sz[4]=0.;
    nd=0.;
    for(n=0;n<LINEPIX2;n++){
      //fprintf(stderr,"%d %g %g %g",n,x[n],y[n],f[n]);getchar();
      if(isnan(z[n])||isnan(y[0][n])||isnan(y[1][n])) continue;
      if(fabs(z[n]-f[n])<clip*sig){
	x1=x[n];
	x2=x1*x[n];
	x3=x2*x[n];
	x4=x3*x[n];
	sx1+=x1;
	sx2+=x2;
	sx3+=x3;
	sx4+=x4;
	sy[0]+=y[0][n];
	sy2[0]+=y[0][n]*y[0][n];
	sxy[0]+=x1*y[0][n];
	sx2y[0]+=x2*y[0][n];
	sy[1]+=y[1][n];
	sy2[1]+=y[1][n]*y[1][n];
	sxy[1]+=x1*y[1][n];
	sx2y[1]+=x2*y[1][n];
	syy+=y[0][n]*y[1][n];
	sz[0]+=z[n];
	sz[1]+=x1*z[n];
	sz[2]+=x2*z[n];
	sz[3]+=y[0][n]*z[n];
	sz[4]+=y[1][n]*z[n];
	nd++;
      }
    }
    sx1/=nd;sx2/=nd;sx3/=nd;sx4/=nd;
    sy[0]/=nd;sy2[0]/=nd;sxy[0]/=nd;sx2y[0]/=nd;
    sy[1]/=nd;sy2[1]/=nd;sxy[1]/=nd;sx2y[1]/=nd;
    syy/=nd;sz[0]/=nd;sz[1]/=nd;sz[2]/=nd;sz[3]/=nd;sz[4]/=nd;
    p[0][0]=1.;
    p[0][1]=sx1;
    p[0][2]=sx2;
    p[0][3]=sy[0];
    p[0][4]=sy[1];
    p[1][0]=sx1;
    p[1][1]=sx2;
    p[1][2]=sx3;
    p[1][3]=sxy[0];
    p[1][4]=sxy[1];
    p[2][0]=sx2;
    p[2][1]=sx3;
    p[2][2]=sx4;
    p[2][3]=sx2y[0];
    p[2][4]=sx2y[1];
    p[3][0]=sy[0];
    p[3][1]=sxy[0];
    p[3][2]=sx2y[0];
    p[3][3]=sy2[0];
    p[3][4]=syy;
    p[4][0]=sy[1];
    p[4][1]=sxy[1];
    p[4][2]=sx2y[1];
    p[4][3]=syy;
    p[4][4]=sy2[1];
    if(sweep(p,q,od+1)) return(-1.);
    a[4]=0.;
    for(i=0;i<=od;i++){
      a[i]=0.;
      for(j=0;j<=od;j++){
	a[i]+=q[i][j]*sz[j];
      }
    }
    dy2s=0.;
    nd=0;
    for(n=0;n<LINEPIX2;n++){
      f[n]=a[0]+(a[1]+a[2]*x[n])*x[n]+a[3]*y[0][n]+a[4]*y[1][n];
      if(fabs(z[n]-f[n])<clip*sig){
	dy2s+=(z[n]-f[n])*(z[n]-f[n]);
	nd++;
      }
    }
    sig=sqrt(dy2s/nd);
  }
  for(n=0;n<LINEPIX2;n++) z[n]-=f[n];
  //fprintf(stderr," %d",LINEPIX2-nd);
  return(sig);
}

int time_sort(const struct dirent **a, const struct dirent **b) {
  struct stat fs[2];
  char fname[2][512];
  sprintf(fname[0],"%s/%s",dir,(*a)->d_name);
  sprintf(fname[1],"%s/%s",dir,(*b)->d_name);
  if(stat(fname[0],&fs[0])==-1||stat(fname[1],&fs[1])==-1) return 0; 
  if(fs[0].st_mtime < fs[1].st_mtime) return 1;
  if(fs[0].st_mtime > fs[1].st_mtime) return -1;
  return 0;
}

int fits_filter(const struct dirent *entry) {
  const char *d_name = entry->d_name;
  size_t len = strlen(d_name);

  struct stat st;
  if (stat(d_name, &st) == 0 && S_ISDIR(st.st_mode)) return 0;
  if (len >= 5 && strcmp(d_name + len - 5, ".fits") == 0){
    //if(strcmp(d_name,"reduc.fits")==0||strcmp(d_name,"double.fits")==0||strcmp(d_name,"double_im.fits")==0||strcmp(d_name,"double_rf.fits")==0||strcmp(d_name,"sigclip.fits")==0||strcmp(d_name,"flat.fits")==0||strcmp(d_name,"dark.fits")==0) return 0;
    if(strncmp(d_name,"mzr",3)==0) return 1;
  }
  return 0;
}

double airindex(double w)
{
  double   u,n3,t=288,ps=1012.5,pw=0,ds,dw;

  u=1000/w;
  ds=ps/t*(1+ps*(57.90e-8-9.3250e-4/t+0.25844/t/t));
  dw=pw/t*(1+pw*(1+3.7e-4*pw)*(-2.37321e-3+2.23366/t-710.792/t/t+7.75141e4/t/t/t));
  n3=(2371.34+683939.7/(130-u*u)+4547.3/(38.9-u*u))*ds+(6487.31+58.058*u*u-0.71150*u*u*u*u+0.08851*u*u*u*u*u*u)*dw;
  return (1+n3*1e-8);
}

double zv(float z[LINEPIY2][LINEPIX2], double x, double y){
  int i,j,n;
  double dx,dy,zt;

  j=(int)x;i=(int)y;dx=x-j;dy=y-i;
  n=0;
  if(isnan(z[i][j])) n++;
  if(isnan(z[i][j+1])) n++;
  if(isnan(z[i+1][j])) n++;
  if(isnan(z[i+1][j+1])) n++;
  if(n==0) zt=z[i][j]*(1-dx)*(1-dy)+z[i][j+1]*dx*(1-dy)+z[i+1][j]*(1-dx)*dy+z[i+1][j+1]*dx*dy;
  else if(n==1){
    if(isnan(z[i][j])) zt=z[i][j+1]*(1-dy)+z[i+1][j]*(1-dx)+z[i+1][j+1]*(dx+dy-1);
    else if(isnan(z[i][j+1])) zt=z[i][j]*(1-dy)+z[i+1][j]*(dy-dx)+z[i+1][j+1]*dx;
    else if(isnan(z[i+1][j])) zt=z[i][j]*(1-dx)+z[i][j+1]*(dx-dy)+z[i+1][j+1]*dy;
    else if(isnan(z[i+1][j+1])) zt=z[i][j]*(1-dx-dy)+z[i][j+1]*dx+z[i+1][j]*dy;
  }
  else{
    if(x<0.5){
      if(y<0.5){
	if(isnan(z[i][j])==0&&isnan(z[i][j+1])==0) zt=z[i][j]*(1-dx)+z[i][j+1]*dx;
	else if(isnan(z[i][j])==0&&isnan(z[i+1][j])==0) zt=z[i][j]*(1-dy)+z[i+1][j]*dy;
	else zt=z[i][j];
      }
      else{
	if(isnan(z[i+1][j])==0&&isnan(z[i+1][j+1])==0) zt=z[i+1][j]*(1-dx)+z[i+1][j+1]*dx;
	else if(isnan(z[i][j])==0&&isnan(z[i+1][j])==0) zt=z[i][j]*(1-dy)+z[i+1][j]*dy;
	else zt=z[i+1][j];
      }
    }
    else{
      if(y<0.5){
	if(isnan(z[i][j])==0&&isnan(z[i][j+1])==0) zt=z[i][j]*(1-dx)+z[i][j+1]*dx;
	else if(isnan(z[i][j+1])==0&&isnan(z[i+1][j+1])==0) zt=z[i][j+1]*(1-dy)+z[i+1][j+1]*dy;
	else zt=z[i][j+1];
      }
      else{
	if(isnan(z[i+1][j])==0&&isnan(z[i+1][j+1])==0) zt=z[i+1][j]*(1-dx)+z[i+1][j+1]*dx;
	else if(isnan(z[i][j+1])==0&&isnan(z[i+1][j+1])==0) zt=z[i][j+1]*(1-dy)+z[i+1][j+1]*dy;
	else zt=z[i+1][j+1];
      }
    }
  }
  return(zt);
}

double plfit(double x[], double y[], double z[], int ndata, double a[])
{
  double sx,sy,sz,sx2,sy2,sz2,sxy,sxz,syz;
  double tx2,ty2,txy,txz,tyz;
  int    n;

  sx=0.;sy=0.;sz=0.;sx2=0.;sy2=0.;sz2=0.;sxy=0.;sxz=0.;syz=0.;
  for(n=0;n<ndata;n++){
    sx+=x[n];
    sy+=y[n];
    sz+=z[n];
    sx2+=x[n]*x[n];
    sy2+=y[n]*y[n];
    sz2+=z[n]*z[n];
    sxy+=x[n]*y[n];
    sxz+=x[n]*z[n];
    syz+=y[n]*z[n];
  }
  sx/=ndata;
  sy/=ndata;
  sz/=ndata;
  sx2/=ndata;
  sy2/=ndata;
  sz2/=ndata;
  sxy/=ndata;
  sxz/=ndata;
  syz/=ndata;
  tx2=sx2-sx*sx;
  ty2=sy2-sy*sy;
  txy=sxy-sx*sy;
  txz=sxz-sx*sz;
  tyz=syz-sy*sz;
  a[0]=(txz*ty2-txy*tyz)/(tx2*ty2-txy*txy);
  a[1]=(tyz*tx2-txy*tyz)/(tx2*ty2-txy*txy);
  a[2]=sz-(a[0]*sx+a[1]*sy);
  return (sqrt(a[0]*a[0]*sx2+a[1]*a[1]*sy2+a[2]*a[2]+sz2+2*a[0]*a[1]*sxy+2*a[0]*a[2]*sx+2*a[1]*a[2]*sy-2*a[0]*sxz-2*a[1]*syz-2*a[2]*sz));
}

void putbytes(FILE *f, int n, unsigned long x)
{
    while (--n >= 0) {
        fputc(x & 255, f);  x >>= 8;
    }
}

void gr_BMP(char *fname, long gr_screen[], int ni, int nj)
{
  int i,j;
  FILE *fp;
    
  fprintf(stderr,"Write: %s ... ",fname);
  fp = fopen(fname, "wb");
  fputs("BM", fp);
  putbytes(fp, 4, nj * ni * 4 + 54);
  putbytes(fp, 4, 0);
  putbytes(fp, 4, 54);
  putbytes(fp, 4, 40);
  putbytes(fp, 4, nj);
  putbytes(fp, 4, ni);
  putbytes(fp, 2, 1);
  putbytes(fp, 2, 32);
  putbytes(fp, 4, 0);
  putbytes(fp, 4, nj * ni * 4);
  putbytes(fp, 4, 3780);
  putbytes(fp, 4, 3780);
  putbytes(fp, 4, 0);
  putbytes(fp, 4, 0);
  for (i = 0; i < ni; i++){
    for (j = 0; j < nj; j++) putbytes(fp, 4, gr_screen[i*nj+j]);
  }
  fclose(fp);
  fprintf(stderr,"OK\n");
}

/*
void loadcmp(long c[]){
  FILE *fp;
  char fname[32];
  int n;

  //sprintf(fname,"viridis.cmp");
  sprintf(fname,"turbo.cmp");
  fprintf(stderr,"Read: %s ... ",fname);
  if((fp=fopen(fname,"r"))==NULL){
    fprintf(stderr,"Can not open file: %s\n",fname);
    exit(1);
  }
  for(n=0;n<256;n++){
    if(fscanf(fp,"%lx",&c[n])!=1){
      fprintf(stderr,"File Truncated ... (%d)\n",n);;
      exit(1);
    }
  }
  fclose(fp);
  fprintf(stderr,"OK\n");
}
*/

void loadcmp(long c[]){
  int n;
  double t,r,g,b;
  
  for (n=0;n<256;n++) {
    t=n/255.;
    if(t<1./3) r=t*3;else r=1.;
    if(t<1./3) g=0.;else if(t<2./3) g=t*3-1.;else g=1.;
    if(t<2./3) b=0.;else b=t*3-2.;
    c[n]=(int)(r*256-0.01)*65536+(int)(g*256-0.01)*256+(int)(b*256-0.01);
  }
  c[255]=255;
}

double tweight(char date[3][16], char jstend[3][16], float exptime){
  struct tm tm[3]={0};
  time_t t[3];
  double dt[3];

  tm[0].tm_isdst=-1;
  tm[1].tm_isdst=-1;
  tm[2].tm_isdst=-1;
  sscanf(date[0],"%d-%d-%d",&tm[0].tm_year,&tm[0].tm_mon,&tm[0].tm_mday);
  sscanf(date[1],"%d-%d-%d",&tm[1].tm_year,&tm[1].tm_mon,&tm[1].tm_mday);
  sscanf(date[2],"%d-%d-%d",&tm[2].tm_year,&tm[2].tm_mon,&tm[2].tm_mday);
  sscanf(jstend[0],"%d:%d:%d",&tm[0].tm_hour,&tm[0].tm_min,&tm[0].tm_sec);
  sscanf(jstend[1],"%d:%d:%d",&tm[1].tm_hour,&tm[1].tm_min,&tm[1].tm_sec);
  sscanf(jstend[2],"%d:%d:%d",&tm[2].tm_hour,&tm[2].tm_min,&tm[2].tm_sec);
  //fprintf(stderr,"%d-%d-%d %d:%d:%d\n",tm[0].tm_year,tm[0].tm_mon,tm[0].tm_mday,tm[0].tm_hour,tm[0].tm_min,tm[0].tm_sec);
  //fprintf(stderr,"%d-%d-%d %d:%d:%d\n",tm[1].tm_year,tm[1].tm_mon,tm[1].tm_mday,tm[1].tm_hour,tm[1].tm_min,tm[1].tm_sec);
  //fprintf(stderr,"%d-%d-%d %d:%d:%d\n",tm[2].tm_year,tm[2].tm_mon,tm[2].tm_mday,tm[2].tm_hour,tm[2].tm_min,tm[2].tm_sec);
  tm[0].tm_year-=1900;
  tm[1].tm_year-=1900;
  tm[2].tm_year-=1900;
  tm[0].tm_mon-=1;
  tm[1].tm_mon-=1;
  tm[2].tm_mon-=1;
  t[0]=mktime(&tm[0]);
  t[1]=mktime(&tm[1]);
  t[2]=mktime(&tm[2]);
  //fprintf(stderr,"%ld %ld %ld\n",t[0],t[2],t[1]);
  dt[0]=difftime(t[2],t[1]);
  dt[1]=difftime(t[2],t[0]);//+exptime/2;
  dt[2]=difftime(t[0],t[1]);//-exptime/2;
  //fprintf(stderr,"%lf %lf %lf\n",dt[0],dt[2],dt[1]);
  return (dt[1]/dt[0]);
}

int main(int argc,char *argv[])
{
  fitsfile *fp[3];
  int i,j,ii,jj,status,nfound,anynull,bitpix=-32;
  long fpixel,naxis=2,naxes[3]={LINEPIX, LINEPIY, 2};
  float nullval,tmpval,exptime;
  static float zi[LINEPIY][LINEPIX2],zo[LINEPIY2][LINEPIX2],zn[LINEPIY2][LINEPIX2],zx[LINEPIY2][LINEPIX2],zy[LINEPIY2][LINEPIX2],zp[LINEPIY2][LINEPIX2],bn[40][40],rfs[NFIBER][LINEPIX2/48],rf0[NFIBER][LINEPIX2],rf[NFIBER][LINEPIX2],sp[NFIBER][LINEPIX2],bg[NFIBER][LINEPIX2],bg0[2][LINEPIX2];
  double kn[PSPEC]={0.2,0.5,0.8,1.,1.,1.,0.8,0.5,0.2},wc,dw,ax[4][3],ay[4][3],xav[3][LINEPIX2],yav[LINEPIY2],fav[2][NFIBER],fav0[NFIBER],kns,bav,bavav,kn1[PSPEC],cw[NCAL],dmy,ww[LINEPIX2],aa[MAXODR],ew,ewn;
  double xct,yct,zt,r,r2,r2s,zs,dxzs,dyzs,drzs,dx2zs,dy2zs,dr2zs,zmax,dx,dy,x[4],y[4],z[4],a[2][3],favav[2],favs,wt,xc[NFIBER],yc[NFIBER],favmin,favmax,bmpz1,bmpz2,ws,wk,sk,sg,sgmax,x1,x2,sig,noi,ymax,ymin;//,bkg;
  int cam,im,jm,lp,n,nf,icen,kn2[4][4]={{0,1,1,0},{1,1,1,1},{1,1,1,1},{0,1,1,0}},dae=0,new=1,jc[NFIBER],ic[NFIBER],ix,jx,nsat,k,m,gn,npix[4],tpix[4]={30000,40000,50000,60000},sft,na1,na2,nc,n1,n2;
  FILE *fpt;
  char buf[1024],fname[5][64],fnamet[96],fnamet2[96],date[3][16],jstend[3][16],gain[16],cn[2]={'C','H'},sdir[32],fene[16],lname[8],object[32];
  long gr_screen[LINEPIY/4*LINEPIX/4],c[256];
  struct stat fs[6];
  struct dirent **namelist;

  if(argc!=2&&argc!=3&&argc!=4&&argc!=5){
    printf("Usage: %s daemon\n",argv[0]);
    printf("  or : %s input output [mzwarp_dxdy.fits [anc.fits or anclist]]\n",argv[0]);
    exit(1);
  }
  kns=0.;
  for(ii=0;ii<PSPEC;ii++){
    //kn[ii]=1.;
    kns+=kn[ii];
  }
  if((fpt = fopen("camnum","r")) == NULL){
    fprintf(stderr,"camnum file is not found...\n");
    exit(0);
  }
  fscanf(fpt,"%d",&cam);
  if(cam!=1) cam=0;
  fclose(fpt);
  for(i=0;i<LINEPIY2;i++){
    for(j=0;j<LINEPIX2;j++) zo[i][j]=0./0.;
  }
  sscanf(argv[0],"%s",sdir);
  sdir[strlen(sdir)-6]='\0';
  sprintf(fnamet,"%smzr_bundle.pos",sdir);
  fprintf(stderr,"Read: %s ... ",fnamet);
  if((fpt=fopen(fnamet,"r"))==NULL){
    fprintf(stderr,"Can not open file: %s\n",fnamet);
    exit(1);
  }
  for(nf=0;nf<NFIBER;nf++){ 
    if(fgets(buf,sizeof(buf),fpt)==NULL){
      fprintf(stderr,"Broken file ?\n");
      exit(1);
    }
    if(sscanf(buf,"%d %d %lf %lf",&jc[nf],&ic[nf],&xc[nf],&yc[nf])!=4){
      fprintf(stderr,"Broken file ?\n");
      exit(1);
    }
    xc[nf]/=110;yc[nf]/=110;
  }
  fclose(fpt);
  fprintf(stderr,"OK\n");
  sprintf(fnamet,"%scal_%c.dat",sdir,cn[cam]);
  fprintf(stderr,"Read: %s ... ",fnamet);
  if((fpt=fopen(fnamet,"r"))==NULL){
    fprintf(stderr,"Can not open file: %s\n",fnamet);
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
    n=sscanf(buf,"%lf %lf %s",&dmy,&cw[nc],lname);
    //fprintf(stderr,"%d: %s",n,buf);
    if(n!=3) continue;
    if(strcmp(lname,"Sky")!=0) continue;
    nc++;
    if(nc==NCAL){
      fprintf(stderr,"NCAL is too small.\n");
      exit(0);
    }
  }
  fclose(fpt);
  fprintf(stderr,"OK\n");
  if(argc==2){
    if(strcmp(argv[1],"daemon")==0) dae=1;
    else{
      printf("Usage: %s daemon\n",argv[0]);
      exit(1);
    }
    system("touch daemon");
    stat("double", &fs[2]);
    dir[0]='\0';
    sprintf(fname[0],"reduc.fits");
    sprintf(fname[1],"double.fits");
    sprintf(fname[2],"mzwarp_dxdy.fits");
    sprintf(fname[3],"anc.fits");
    argc=5;
    //system("sxiv spec.png &");
  }
  else{
    for(n=1;n<argc;n++) sscanf(argv[n],"%s",fname[n-1]);
  }
  fprintf(stderr,"Read: flat.dat ... ");
  if((fpt=fopen("flat.dat","r"))==NULL){
    fprintf(stderr,"Can not open.\n");
    fprintf(stderr,"Relative flux correction will not be applied.\n");
    for(nf=0;nf<NFIBER;nf++) fav0[nf]=1.;
  }
  else{
    for(nf=0;nf<NFIBER;nf++){
      if(fscanf(fpt,"%lf",&fav0[nf])!=1){
	fprintf(stderr,"Broken file ? (%d)\n",nf);
	exit(1);
      }
    }
    fclose(fpt);
    fprintf(stderr,"OK\n");
  }
  fprintf(stderr,"Read: flats.fits ... ");
  status=0;
  if(fits_open_file(&fp[1], "flats.fits", READONLY, &status)){
    fprintf(stderr,"Can not open.\n");
    fprintf(stderr,"Relative spectum correction will not be applied.\n");
    for(i=0;i<NFIBER;i++){
      for(j=0;j<LINEPIX2;j++) rf0[i][j]=1.;
    }
  }
  else{
    if(fits_read_keys_lng(fp[1], "NAXIS", 1, 2, naxes, &nfound, &status)) printerror(status);
    fprintf(stderr,"%ld x %ld format ",naxes[0],naxes[1]);
    if((naxes[0]!=LINEPIX2)||(naxes[1]!=NFIBER)){
      fprintf(stderr,"\nIllegal image format\n");
      exit(0);
    }
    fpixel = 1;
    nullval = 0;
    for(i=0;i<naxes[1];i++){
      if (fits_read_img(fp[1], TFLOAT, fpixel, naxes[0], &nullval, rf0[i], &anynull, &status)) printerror(status);
      for(j=0;j<naxes[0];j++){
	if(isinf(rf0[i][j])) rf0[i][j]=1.;
      }
      fpixel+=naxes[0];
    }
    if(fits_close_file(fp[1], &status)) printerror(status);
    fprintf(stderr,"OK\n");
  }
  loadcmp(c);
  if(dae){
    if(stat(fname[3], &fs[0])) {
      fprintf(stderr,"Cannot read node information of %s\n", fname[3]);
      exit(1);
    }
  }
  k=0;
  while(1){
    if(!dae){
      if(k==0){
	fprintf(stderr,"Read: %s ... ",fname[0]);
	status=0;
	if(fits_open_file(&fp[0], fname[0], READONLY, &status)){
	  sprintf(fnamet,"%s_%03d.fits",fname[0],k+1);
	  status=0;
	  if(fits_open_file(&fp[0], fnamet, READONLY, &status)){
	    fprintf(stderr,"Can not open file: %s , %s_001.fits\n",fname[0],fname[0]);
	    exit(0);
	  }
	  k++;
	}
      }
      else{
	sprintf(fnamet,"%s_%03d.fits",fname[0],k+1);
	status=0;
	if(fits_open_file(&fp[0], fnamet, READONLY, &status)){
	  fprintf(stderr,"Loop ends. Number of files: %d\n",k);
	  exit(0);
	}
	k++;
      }
      if(fits_read_key(fp[0], TSTRING, "DATE-OBS", date[0], NULL, &status)) printerror(status);
      if(fits_read_key(fp[0], TSTRING, "JST-END", jstend[0], NULL, &status)) printerror(status);
      if(fits_read_key(fp[0], TSTRING, "OBJECT", object, NULL, &status)) printerror(status);
      if(fits_read_key(fp[0], TFLOAT, "EXPTIME", &exptime, NULL, &status)) printerror(status);
    }
    if(new&&argc>=4){
      fprintf(stderr,"Read: %s ... ",fname[2]);
      status=0;
      if(fits_open_file(&fp[1], fname[2], READONLY, &status)) printerror(status);
      if(fits_read_keys_lng(fp[1], "NAXIS", 1, 2, naxes, &nfound, &status)) printerror(status);
      if(fits_read_key(fp[1], TFLOAT, "MZ_ANCX1", &tmpval, NULL, &status)) printerror(status);
      ax[0][0]=tmpval-1;ax[0][1]=tmpval-1;
      if(fits_read_key(fp[1], TFLOAT, "MZ_ANCY1", &tmpval, NULL, &status)) printerror(status);
      ay[0][0]=tmpval-1;ay[0][1]=tmpval-1;
      if(fits_read_key(fp[1], TFLOAT, "MZ_ANCX2", &tmpval, NULL, &status)) printerror(status);
      ax[1][0]=tmpval-1;ax[1][1]=tmpval-1;
      if(fits_read_key(fp[1], TFLOAT, "MZ_ANCY2", &tmpval, NULL, &status)) printerror(status);
      ay[1][0]=tmpval-1;ay[1][1]=tmpval-1;
      if(fits_read_key(fp[1], TFLOAT, "MZ_ANCX3", &tmpval, NULL, &status)) printerror(status);
      ax[2][0]=tmpval-1;ax[2][1]=tmpval-1;
      if(fits_read_key(fp[1], TFLOAT, "MZ_ANCY3", &tmpval, NULL, &status)) printerror(status);
      ay[2][0]=tmpval-1;ay[2][1]=tmpval-1;
      if(fits_read_key(fp[1], TFLOAT, "MZ_ANCX4", &tmpval, NULL, &status)) printerror(status);
      ax[3][0]=tmpval-1;ax[3][1]=tmpval-1;
      if(fits_read_key(fp[1], TFLOAT, "MZ_ANCY4", &tmpval, NULL, &status)) printerror(status);
      ay[3][0]=tmpval-1;ay[3][1]=tmpval-1;
      fprintf(stderr,"%ld x %ld x %ld format ",naxes[0],naxes[1],naxes[2]);
      if((naxes[0]!=LINEPIX2)||(naxes[1]!=LINEPIY2)||(naxes[2]!=2)){
	fprintf(stderr,"\nIllegal image format\n");
	exit(0);
      }
      fpixel = 1;
      nullval = 0;
      for(i=0;i<naxes[1];i++){
	if (fits_read_img(fp[1], TFLOAT, fpixel, naxes[0], &nullval, zx[i], &anynull, &status)) printerror(status);
	fpixel+=naxes[0];
      }
      for(i=0;i<naxes[1];i++){
	if (fits_read_img(fp[1], TFLOAT, fpixel, naxes[0], &nullval, zy[i], &anynull, &status)) printerror(status);
	fpixel+=naxes[0];
      }
      if(fits_close_file(fp[1], &status)) printerror(status);
      fprintf(stderr,"OK\n");
      if(argc==5){
	fprintf(stderr,"*** Anchor correction ***\n");
	status=0;
	na1=0;na2=0;
	if(fits_open_file(&fp[1], fname[3], READONLY, &status)==0){
	  fprintf(stderr,"Read: %s ... ",fname[3]);
	  if(fits_read_key(fp[1], TSTRING, "DATE-OBS", date[1], NULL, &status)) printerror(status);
	  if(fits_read_key(fp[1], TSTRING, "JST-END", jstend[1], NULL, &status)) printerror(status);
	}
	else{
	  sprintf(fnamet,"%s_%03d.fits",fname[3],na1);
	  //fprintf(stderr,"Read: %s ...\n",fnamet);
	  status=0;
	  if(fits_open_file(&fp[1], fnamet, READONLY, &status)){
	    fprintf(stderr,"Can not open file: %s\n",fnamet);
	    exit(0);
	  }
	  if(fits_read_key(fp[1], TSTRING, "DATE-OBS", date[1], NULL, &status)) printerror(status);
	  if(fits_read_key(fp[1], TSTRING, "JST-END", jstend[1], NULL, &status)) printerror(status);
	  na2=1;
	  while(1){
	    sprintf(fnamet,"%s_%03d.fits",fname[3],na2);
	    //fprintf(stderr,"Read: %s ...\n",fnamet);
	    status=0;
	    if(fits_open_file(&fp[2], fnamet, READONLY, &status)){
	      if(na2==1){
		fprintf(stderr,"Can not open file: %s\n",fnamet);
		exit(0);
	      }
	      else{
		na2=na1;
		status=0;
		break;
	      }
	    }
	    if(fits_read_key(fp[2], TSTRING, "DATE-OBS", date[2], NULL, &status)) printerror(status);
	    if(fits_read_key(fp[2], TSTRING, "JST-END", jstend[2], NULL, &status)) printerror(status);
	    wt=tweight(date,jstend,exptime);
	    //fprintf(stderr,"wt=%g %s %s %s",wt,jstend[0],jstend[1],jstend[2]);getchar();
	    if(wt>=0.){
	      if(wt>1.){
		if(fits_close_file(fp[2], &status)) printerror(status);
		na2=na1;
	      }
	      break;
	    }
	    if(fits_close_file(fp[1], &status)) printerror(status);
	    fp[1]=fp[2];
	    strcpy(date[1],date[2]);
	    strcpy(jstend[1],jstend[2]);
	    na1++;na2++;
	  }
	}
	fprintf(stderr,"OK\n");
	if(na1!=na2) fprintf(stderr,"Anchor: %g*%s_%03d.fits%+g*%s_%03d.fits\n",wt,fname[3],na1,1-wt,fname[3],na2);
	for(m=1;m<=na2-na1+1;m++){
	  if(fits_read_keys_lng(fp[m], "NAXIS", 1, 2, naxes, &nfound, &status)) printerror(status);
	  //fprintf(stderr,"%ld x %ld format ",naxes[0],naxes[1]);
	  if((naxes[0]!=LINEPIX)||(naxes[1]!=LINEPIY)){
	    fprintf(stderr,"\nIllegal image format\n");
	    exit(0);
	  }
	  fpixel = 1;
	  nullval = 0;
	  for(i=0;i<naxes[1];i++){
	    if (fits_read_img(fp[m], TFLOAT, fpixel, naxes[0], &nullval, zi[i], &anynull, &status)) printerror(status);
	    fpixel+=naxes[0];
	  }
	  if(fits_close_file(fp[m], &status)) printerror(status);
	  //fprintf(stderr,"OK\n");
	  resize(zi,naxes,icen);
	  for(n=0;n<4;n++){
	    //fprintf(stderr,"%d %d %d\n",(int)ax[n][0],(int)ay[n][0],HW);
	    zmax=0.;
	    for(i=(int)ay[n][0]-HW*20;i<=(int)ay[n][0]+HW*20;i++){
	      if(i<0) i=0;
	      if(i>=naxes[1]) break;
	      for(j=(int)ax[n][0]-HW*20;j<=(int)ax[n][0]+HW*20;j++){
		if(j<0) j=0;
		if(j>=naxes[0]) break;
		//fprintf(stderr," %6.1lf",zi[i][j]);
		//fprintf(stderr," (%d,%d)",j+1,i+1);
		if(zmax<zi[i][j]){
		  jm=j;
		  im=i;
		  zmax=zi[i][j];
		}
	      }
	      //fprintf(stderr,"\n");
	    }
	    fprintf(stderr,"Peak [%d,%d]: %g ",jm+1,im+1,zmax);
	    if(zmax<100){
	      fprintf(stderr,"No spot found.\n");
	      break;
	    }
	    xct=jm;yct=im;
	    for(lp=0;lp<5;lp++){
	      zs=0.;dxzs=0.;dyzs=0.;drzs=0.;dx2zs=0.;dy2zs=0.;dr2zs=0.;
	      for(i=-HW;i<=HW;i++){
		for(j=-HW;j<=HW;j++){
		  zt=zv(zi,xct+j,yct+i);
		  if(i==0&&j==0) zmax=zt;
		  //fprintf(stderr," %g",zt);
		  r=sqrt((double)(i*i+j*j));
		  zs+=zt;
		  dxzs+=j*zt;
		  dyzs+=i*zt;
		  if(r>0.1) drzs+=zt/r;//weight
		  dx2zs+=j*j*zt;
		  dy2zs+=i*i*zt;
		  dr2zs+=r*zt;//r*r*(zt/r)
		}
		//fprintf(stderr,"\n");
	      }
	      dx2zs/=zs;
	      dy2zs/=zs;
	      dr2zs/=drzs;
	      dxzs/=zs;
	      dyzs/=zs;
	      drzs=(zs-zmax)/drzs;
	      dx2zs=sqrt(dx2zs);//sqrt(dx2zs-dxzs*dxzs);
	      dy2zs=sqrt(dy2zs);//sqrt(dy2zs-dyzs*dyzs);
	      dr2zs=sqrt(dr2zs);
	      xct+=dxzs;
	      yct+=dyzs;
	      //fprintf(stderr,"%g %g %g %g",xct,yct,dx2zs,dy2zs);getchar();
	    }
	    ax[n][m]=xct;ay[n][m]=yct;
	    fprintf(stderr,"dr/dx/dy=%.3lf/%.3lf/%.3lf (%.3lf %.3lf) @%s\n",dr2zs,dx2zs,dy2zs,ax[n][m],ay[n][m],jstend[m]);
	  }
	  if(n<4) break;
	}
	if(na1!=na2){
	  for(n=0;n<4;n++){
	    ax[n][1]=ax[n][1]*wt+ax[n][2]*(1-wt);
	    ay[n][1]=ay[n][1]*wt+ay[n][2]*(1-wt);
	  }
	}
	dy=0.;
	for(n=0;n<4;n++) dy+=ay[n][1]-ay[n][0];
	fprintf(stderr,"icen=%d => ",icen);
	icen+=(int)(dy/8);
	fprintf(stderr,"%d\n",icen);
	for(n=0;n<4;n++) ay[n][1]-=(int)(dy/8)*2;
	for(n=0;n<4;n++){
	  fprintf(stderr,"(%.3lf %.3lf) => (%.3lf %.3lf)\n",ax[n][1],ay[n][1],ax[n][0],ay[n][0]);
	  x[n]=ax[n][0];
	  y[n]=ay[n][0];
	  z[n]=ax[n][1]-ax[n][0];
	}
	plfit(x,y,z,4,a[0]);
	//for(n=0;n<4;n++) fprintf(stderr,"%dx: %g => %g\n",n,z[n],a[0][0]*x[n]+a[0][1]*y[n]+a[0][2]);
	for(n=0;n<4;n++) z[n]=ay[n][1]-ay[n][0];
	plfit(x,y,z,4,a[1]);
	//for(n=0;n<4;n++) fprintf(stderr,"%dy: %g => %g\n",n,z[n],a[1][0]*x[n]+a[1][1]*y[n]+a[1][2]);
	for(i=0;i<naxes[1];i++){
	  for(j=0;j<naxes[0];j++) zi[i][j]=zx[i][j];
	}
	for(i=0;i<naxes[1];i++){
	  for(j=0;j<naxes[0];j++){
	    xct=j-(a[0][0]*j+a[0][1]*i+a[0][2]);
	    yct=i-(a[1][0]*j+a[1][1]*i+a[1][2]);
	    if(xct<0) xct=0.;
	    else if(xct>naxes[0]-1) xct=naxes[0]-1;
	    if(yct<0) yct=0.;
	    else if(yct>naxes[1]-1) yct=naxes[1]-1;
	    jj=(int)xct;dx=xct-jj;
	    ii=(int)yct;dy=yct-ii;
	    zx[i][j]=zi[ii][jj]*(1-dx)*(1-dy)+zi[ii][jj+1]*dx*(1-dy)+zi[ii+1][jj]*(1-dx)*dy+zi[ii+1][jj+1]*dx*dy-(a[0][0]*j+a[0][1]*i+a[0][2]);
	  }
	}
	for(i=0;i<naxes[1];i++){
	  for(j=0;j<naxes[0];j++) zi[i][j]=zy[i][j];
	}
	for(i=0;i<naxes[1];i++){
	  for(j=0;j<naxes[0];j++){
	    xct=j-(a[0][0]*j+a[0][1]*i+a[0][2]);
	    yct=i-(a[1][0]*j+a[1][1]*i+a[1][2]);
	    if(xct<0) xct=0.;
	    else if(xct>naxes[0]-1) xct=naxes[0]-1;
	    if(yct<0) yct=0.;
	    else if(yct>naxes[1]-1) yct=naxes[1]-1;
	    jj=(int)xct;dx=xct-jj;
	    ii=(int)yct;dy=yct-ii;
	    zy[i][j]=zi[ii][jj]*(1-dx)*(1-dy)+zi[ii][jj+1]*dx*(1-dy)+zi[ii+1][jj]*(1-dx)*dy+zi[ii+1][jj+1]*dx*dy-(a[1][0]*j+a[1][1]*i+a[1][2]);
	  }
	}
      }
    }
    if(dae){
      while(1){
	if(stat("daemon", &fs[5])) exit(0);
	sprintf(fnamet,"%sBH67.path",sdir);
	if(stat(fnamet,&fs[4])&&strcmp(dir,".")!=0){
	  fprintf(stderr,"BH67.path does not exist. Use current directory.\n");
	  strcpy(dir,".");
	}
	else{
	  if(fs[2].st_mtime < fs[4].st_mtime || dir[0]=='\0'){
	    fs[2]=fs[4];
	    fprintf(stderr,"Read: %s ... ",fnamet);
	    if((fpt = fopen(fnamet,"r")) == NULL){
	      fprintf(stderr,"Cannot open. Use current directory.\n");
	      strcpy(dir,".");
	    }
	    else{
	      if(fscanf(fpt,"%s",dir)!=1){
		fprintf(stderr,"Empty file. Use current directory.\n");
		strcpy(dir,".");
	      }
	      else fprintf(stderr,"%s ... OK\n",dir);
	      fclose(fpt);
	    }
	  }
	}
	usleep(200000);
	if(stat(fname[3], &fs[1])) {
	  fprintf(stderr,"Cannot read node information of %s\n", fname[3]);
	  exit(1);
	}
	if(fs[0].st_mtime < fs[1].st_mtime){
	  fs[0]=fs[1];
	  new=1;
	  break;
	}
	else new=0;
	n=scandir(dir, &namelist, fits_filter, time_sort);
	if(n<0) {
	  fprintf(stderr,"Scandir Error (%s)\n",dir);
	  exit(1);
	}
	if(n>0){
	  sprintf(buf,"%s/%s",dir,namelist[0]->d_name);
	  stat(buf,&fs[3]);
	  if(fs[2].st_mtime < fs[3].st_mtime){
	    fs[2]=fs[3];
	    sprintf(buf,"%sreduc %s/%s %s bias.fits flat.fits",sdir,dir,namelist[0]->d_name,fname[0]);
	    fprintf(stderr,"%s\n",buf);
	    if(system(buf)!=0) continue;
	    break;
	  }
	}
	usleep(200000);
      }
      if(new) continue;
      status=0;
      if(fits_open_file(&fp[0], fname[0], READONLY, &status)) printerror(status);
      if(fits_read_key(fp[0], TSTRING, "DATE-OBS", date[0], NULL, &status)) printerror(status);
      if(fits_read_key(fp[0], TSTRING, "JST-END", jstend[0], NULL, &status)) printerror(status);
      if(fits_read_key(fp[0], TSTRING, "OBJECT", object, NULL, &status)) printerror(status);
      if(fits_read_key(fp[0], TFLOAT, "EXPTIME", &exptime, NULL, &status)) printerror(status);
    }
    jm=strlen(object);
    for(j=0;j<jm;j++){
      //fprintf(stderr,"%d: %s\n",j,object);
      if(object[j]=='"'){
	jm--;
	for(jj=j;jj<=jm;jj++) object[jj]=object[jj+1];
      }
    }
    if(k==0) fprintf(stderr,"Read: %s ... ",fname[0]);
    else fprintf(stderr,"%s_%03d.fits ... ",fname[0],k);
    if(fits_read_keys_lng(fp[0], "NAXIS", 1, 2, naxes, &nfound, &status)) printerror(status);
    fprintf(stderr,"%ld x %ld format ",naxes[0],naxes[1]);
    if((naxes[0]!=LINEPIX)||(naxes[1]!=LINEPIY)){
      fprintf(stderr,"\nIllegal image format\n");
      exit(0);
    }
    if(fits_read_key(fp[0], TSTRING, "GAIN", gain, NULL, &status)) printerror(status);
    if(fits_read_key(fp[0], TSTRING, "FENE-LMP", fene, NULL, &status)) printerror(status);
    fpixel = 1;
    nullval = 0;
    zs=0.;sg=0.;n=0;
    for(i=0;i<naxes[1];i++){
      if(fits_read_img(fp[0], TFLOAT, fpixel, naxes[0], &nullval, zi[i], &anynull, &status)) printerror(status);
      fpixel+=naxes[0];
      for(j=0;j<naxes[0];j++){
	if(isnan(zi[i][j])==0){
	  zs+=zi[i][j];
	  sg+=zi[i][j]*zi[i][j];
	  n++;
	}
      }
    }
    npix[0]=0.;npix[1]=0.;npix[2]=0.;npix[3]=0.;
    for(i=icen-200;i<icen+200;i++){
      for(j=0;j<naxes[0];j++){
	if(zi[i][j]>(float)tpix[0]) npix[0]++;
	if(zi[i][j]>(float)tpix[1]) npix[1]++;
	if(zi[i][j]>(float)tpix[2]) npix[2]++;
	if(zi[i][j]>(float)tpix[3]) npix[3]++;
      }
    }
    if(fits_close_file(fp[0], &status)) printerror(status);
    if(strncmp(gain,"x16",3)==0) gn=16;
    else if(strncmp(gain,"x8",2)==0) gn=8;
    else if(strncmp(gain,"x4",2)==0) gn=4;
    else if(strncmp(gain,"x2",2)==0) gn=2;
    else gn=1;
    fprintf(stderr,"gain=x%d OK\n",gn);
    zs/=n;
    sg=sqrt(sg/n-zs*zs);
    fprintf(stderr,"av:%g sg:%g\n",zs,sg);
    if(naxes[0]==LINEPIX2&&naxes[1]==LINEPIY2){jx=8;jm=264;ix=2;im=40;}
    else {jx=4;jm=0;ix=4;im=0;}
    nsat=0;
    for(i=0;i<LINEPIY/4;i++){
      for(j=0;j<LINEPIX/4;j++){
	zt=0.;
	for(ii=i*ix+im;ii<(i+1)*ix+im;ii++){
	  for(jj=j*jx+jm;jj<(j+1)*jx+jm;jj++){
	    if(zi[ii][jj]>SATU&&abs(ii-icen)<LINEPIY2/4){zt=0./0.;nsat++;}
	    else zt+=zi[ii][jj];
	  }
	}
	zt/=16.;
	if(isnan(zt)) n=255;
	else{
	  n=(zt-zs+sg)/(sg*4)*255;
	  if(n<0) n=0;
	  else if(n>254) n=254;
	}
	gr_screen[i*LINEPIX/4+j]=c[n];
      }
    }
    gr_BMP("double.bmp",gr_screen,LINEPIY/4,LINEPIX/4);
    resize(zi,naxes,icen);
    if(argc>=4){
      warp(zi,zo,zx,zy);
      for(j=0;j<naxes[0];j++){bg0[0][j]=0.;bg0[1][j]=0.;}
      jj=0;
      for(j=0;j<LINEPIX;j++){
	jm=j*2-LINEPIX+naxes[0]/2;
	x1=(zx[icen][jm-2]-2+zx[icen][jm])/2+jm;
	x2=(zx[icen][jm+2]+2+zx[icen][jm])/2+jm;
	dy=(zy[icen+1][jm]-zy[icen-1][jm])/2;
	r=(x2-x1)*(1+dy)*2;
	while(jj-0.5<x2){
	  if(x1<jj+0.5){
	    if(jj-0.5<x1){
	      bg0[0][jj]+=zi[LINEPIY-2][j]/r*(jj+0.5-x1);
	      bg0[1][jj]+=zi[LINEPIY-1][j]/r*(jj+0.5-x1);
	      //if(jj>2424){fprintf(stderr,"L: %d %d: %g %g %g %g %g %g",jj,j,x1,x2,bg0[0][jj],zi[LINEPIY-2][j],bg0[1][jj],zi[LINEPIY-1][j]);getchar();}
	    }
	    else if(x2<jj+0.5){
	      bg0[0][jj]+=zi[LINEPIY-2][j]/r*(x2-(jj-0.5));
	      bg0[1][jj]+=zi[LINEPIY-1][j]/r*(x2-(jj-0.5));
	      //if(jj>2424){fprintf(stderr,"R: %d %d: %g %g %g %g %g %g",jj,j,x1,x2,bg0[0][jj],zi[LINEPIY-2][j],bg0[21][jj],zi[LINEPIY-1][j]);getchar();}
	      break;
	    }
	    else{
	      bg0[0][jj]+=zi[LINEPIY-2][j]/r;
	      bg0[1][jj]+=zi[LINEPIY-1][j]/r;
	      //if(jj>2424){fprintf(stderr,"C: %d %d: %g %g %g %g %g %g",jj,j,x1,x2,bg0[0][jj],zi[LINEPIY-2][j],bg0[1][jj],zi[LINEPIY-1][j]);getchar();}
	    }
	  }
	  jj++;
	}
      }
      for(i=0;i<naxes[1];i++){
	ws=((double)naxes[1]/2+400-i)/800;
	for(j=0;j<naxes[0];j++){
	  zn[i][j]=bg0[0][j]*ws+bg0[1][j]*(1-ws);
	  //if(i==naxes[1]/2&&j==naxes[0]/2){fprintf(stderr,"%g %g %g %g",ws,bg0[0][j],bg0[1][j],zn[i][j]);getchar();}
	}
      }
      for(nf=0;nf<NFIBER;nf++){
	i=LINE0+PSPEC*nf-PSPEC/2;
	for(ii=0;ii<PSPEC;ii++){
	  for(j=0;j<naxes[0];j++){
	    zo[i+ii][j]/=rf0[nf][j];
	    zn[i+ii][j]/=rf0[nf][j]*rf0[nf][j];
	    //if(j==naxes[0]/2){fprintf(stderr,"%d %d %d %g %g",j,nf,i+ii,zo[i+ii][j],rf0[nf][j]);getchar();}
	  }
	}
      }
      for(i=0;i<naxes[1];i++){
	for(j=0;j<naxes[0];j++) zp[i][j]=zo[i][j];
      }
      if(strcmp(fene,"ON")){
	for(j=0;j<naxes[0];j++){
	  ws=ave(zo,j,0);
	  wt=ave(zo,j,1);
	  //if(j==naxes[0]/2) fprintf(stderr,"%g %g\n",ws,wt);
	  for(i=0;i<naxes[1];i++){
	    zi[i][j]=(wt-ws)/(naxes[1]-LINE0+PSPEC)*(i-(LINE0-PSPEC)/2)+ws;
	  }
	}
	for(i=0;i<naxes[1];i++){
	  for(j=0;j<naxes[0];j++){
	    zo[i][j]-=zi[i][j];
	    //if(j==700){fprintf(stderr,"%d %d %g %g",i,j,zo[i][j],zi[i][j]);getchar();}
	  }
	}
      }
      for(lp=0;lp<4;lp++){
	if(lp==3){
	  for(i=0;i<naxes[1];i++){
	    for(j=0;j<naxes[0];j++){
	      //jm=j-zx[i][j];im=i-zy[i][j];
	      //dx=(zx[im][jm+1]-zx[im][jm-1])/2;
	      //dy=(zy[im+1][jm]-zy[im-1][jm])/2;
	      //r=(1+dx)*(1+dy)*4;
	      //if(i==naxes[1]/2&&j==naxes[0]/2){fprintf(stderr,"%d %d %g %g %g %d %d %g %g %g %g %g",j,jm,zx[im][jm+1],zx[im][jm-1],dx,i,im,zy[im+1][jm],zy[im-1][jm],dy,r,zn[i][j]);getchar();}
	      //zn[i][j]/=r; // pixel 細分化でノイズが増加しないことを無視する場合はここはコメントアウト
	      //if(j==naxes[0]/2){fprintf(stderr,"%d %d %g",i,j,zn[i][j]);getchar();}
	      if(zo[i][j]>0.) zn[i][j]+=zo[i][j];
	      zn[i][j]/=CONVF*16/gn;
	    }
	  }
	}
	for(i=0;i<naxes[1];i++){
	  yav[i]=0.;n=0;
	  for(j=1;j<naxes[0]-1;j++){
	    if(isnan(zo[i][j-1])==0&&isnan(zo[i][j])==0&&isnan(zo[i][j+1])==0){
	      yav[i]+=zo[i][j];
	      n++;
	    }
	  }
	  yav[i]/=n;
	  //fprintf(stderr,"%d: %g\n",i,yav[i]);
	}//getchar();
	favav[0]=0.;favav[1]=0.;bavav=0.;
	favmin=1e99;favmax=-1e99;
	for(nf=0;nf<NFIBER;nf++){
	  i=LINE0+PSPEC*nf-PSPEC/2;
	  fav[0][nf]=0.;bav=0.;
	  for(ii=0;ii<PSPEC;ii++){
	    fav[0][nf]+=yav[i+ii]*kn[ii];
	    bav+=yav[i+ii]*(1-kn[ii]);
	  }
	  fav[1][nf]=fav[0][nf]/fav0[nf];
	  bav/=fav0[nf];
	  if(favmin>fav[1][nf]) favmin=fav[1][nf];
	  if(favmax<fav[1][nf]) favmax=fav[1][nf];
	  favav[0]+=fav[0][nf];
	  favav[1]+=fav[1][nf];
	  bavav+=bav;
	  //fprintf(stderr,"%d: %g %g %g\n",nf,fav[0][nf],fav[1][nf],fav0[nf]);
	}//getchar();
	favav[0]/=NFIBER;
	favav[1]/=NFIBER;
	bavav*=kns/(PSPEC-kns)/NFIBER;
	fprintf(stderr,"Fiber flux (min/ave/max/bave): %g/%g/%g/%g\n",favmin,favav[1],favmax,bavav);
	for(i=0;i<IMAX;i++){
	  for(j=0;j<JMAX;j++) bn[i][j]=0.;
	}
	for(nf=0;nf<NFIBER;nf++){
	  for(ii=0;ii<4;ii++){
	    for(jj=0;jj<4;jj++) bn[ic[nf]+ii][jc[nf]+jj]+=fav[1][nf]*kn2[ii][jj];
	  }
	}
	if(lp==3){
	  for(ii=0;ii<PSPEC;ii++) kn1[ii]=0.;
	  for(nf=0;nf<NFIBER;nf++){
	    i=LINE0+PSPEC*nf-PSPEC/2;
	    for(ii=0;ii<PSPEC;ii++) kn1[ii]+=yav[i+ii];
	  }
	  if(favav[1]>100 && (kn1[PSPEC/2-1]<kn1[PSPEC/2-2]||kn1[PSPEC/2+1]<kn1[PSPEC/2+2])){
	    fprintf(stderr,"Position Shift is detected !\n");
	    for(ii=0;ii<PSPEC;ii++) fprintf(stderr," %g",kn1[ii]);
	    fprintf(stderr,"\n");
	    sft=1;
	  }
	  else sft=0;
	  remove("double_im.fits");
	  fprintf(stderr,"Write double_im.fits ... ");
	  status=0;
	  naxes[0]=JMAX;naxes[1]=IMAX;
	  if (fits_create_file(&fp[1], "double_im.fits", &status)) printerror(status);
	  if (fits_create_img(fp[1], bitpix, naxis, naxes, &status)) printerror(status);
	  fpixel=1;
	  fprintf(stderr,"%ld x %ld format ",naxes[0],naxes[1]);
	  for(i=0;i<naxes[1];i++){
	    if(fits_write_img(fp[1], TFLOAT, fpixel, naxes[0], bn[i], &status)) printerror(status);
	    fpixel+=naxes[0];
	  }
	  if (fits_close_file(fp[1], &status)) printerror(status);
	  fprintf(stderr,"OK\n");
	  naxes[0]=LINEPIX2;naxes[1]=LINEPIY2;
	}
	if(favav[1]<10){
	  fprintf(stderr,"Dark image. Equal weight will be applied.\n");
	  for(nf=0;nf<NFIBER;nf++){
	    fav[0][nf]=1.;
	    fav[1][nf]=1.;
	  }
	  sk=0.;
	}
	else{
	  for(nf=0;nf<NFIBER;nf++){
	    fav[0][nf]/=favav[0];
	    fav[1][nf]/=favav[1];
	    //fprintf(stderr,"%d: %g %g %g\n",nf,fav[0][nf],fav0[nf],fav[1][nf]);
	  }
	  if(lp<3){
	    n=0;sk=0.;
	    for(nf=0;nf<NFIBER;nf++){
	      if(fav[1][nf]<SKY){
		sk+=fav[1][nf];
		n++;
	      }
	    }
	    if(n>0) sk/=n;
	    fprintf(stderr,"Sky level: %g\n",sk);
	    for(nf=0;nf<NFIBER;nf++){
	      if(fav[1][nf]<SKY) fav[0][nf]=0.;
	      else fav[0][nf]*=1-sk/fav[1][nf];
	    }
	  }
	  else sk=med2(fav[0]);
	}
	for(j=0;j<naxes[0];j++){
	  xav[0][j]=0.;xav[1][j]=0.;xav[2][j]=0.;
	  ws=0.,wt=0.;
	  for(nf=0;nf<NFIBER;nf++){
	    sp[nf][j]=0.;wk=0.;
	    if(lp==3) bg[nf][j]=0.;
	    i=LINE0+PSPEC*nf-PSPEC/2;
	    for(ii=0;ii<PSPEC;ii++){
	      if(isnan(zo[i+ii][j])==0){
		sp[nf][j]+=zo[i+ii][j]*kn[ii];
		wk+=kn[ii];
		//if(j==3057){fprintf(stderr,"%d %d %g %g %g %g",nf,j,sp[nf][j],zo[i+ii][j],kn[ii],wk);getchar();}
		if(lp==3){
		  bg[nf][j]+=zn[i+ii][j]*kn[ii]*kn[ii];
		  //if(j==naxes[0]/2){fprintf(stderr,"%d %d %g %g %g",nf,i+ii,zo[i+ii][j],zn[i+ii][j],zn[i+ii][j]/zo[i+ii][j]);getchar();}
		}
	      }
	    }
	    if(wk<1.) sp[nf][j]=0./0.;
	    else sp[nf][j]*=kns/wk;
	    if(lp==3){
	      if(wk<1.) bg[nf][j]=0./0.;
	      else bg[nf][j]*=kns/wk*kns/wk;
	      if(fav[0][nf]>1e-5){
		xav[0][j]+=sp[nf][j]*fav[0][nf];
		wt+=fav[0][nf];
		if(isnan(sp[nf][j])==0) xav[2][j]+=bg[nf][j]*fav[0][nf]*fav[0][nf];
		//if(j==naxes[0]/2){fprintf(stderr,"%d %g %g %g %g %g",nf,sp[nf][j],bg[nf][j],fav[0][nf],xav[0][j],xav[2][j]);getchar();}
	      }
	      if(fav[0][nf]-sk<1e-5){
		  xav[1][j]+=sp[nf][j];
		  ws+=1.;
		  //if(j==naxes[0]/2) fprintf(stderr,"%d: %g %g %g\n",nf,ws,fav[0][nf],sk);
	      }
	    }
	    else{
	      if(fav[0][nf]>1e-5){
		if(isnan(sp[nf][j])==0){
		  xav[0][j]+=sp[nf][j]*fav[0][nf];
		  wt+=fav[0][nf];
		  //fprintf(stderr,"%d %d %g %g %g %g",nf,j,xav[0][j],sp[nf][j],rf0[nf][j],fav[0][nf]);getchar();
		}
	      }
	      else{
		if(isnan(sp[nf][j])==0){
		  xav[1][j]+=sp[nf][j];
		  ws+=1.;
		  //if(j==288 || j==289) fprintf(stderr,"%d %d: %g %g %g\n",j,nf,xav[1][j],sp[nf][j],ws);
		  //fprintf(stderr,"%d %d %g %g %g %g",nf,j,xav[1][j],sp[nf][j],rf0[nf][j],fav[0][nf]);getchar();
		}
	      }
	    }
	  }
	  xav[0][j]/=wt;
	  if(ws>0.5) xav[1][j]/=ws;
	  if(lp==3){
	    xav[2][j]=sqrt(xav[2][j])/wt;
	    //fprintf(stderr,"%d %g %g %g",j,xav[0][j],xav[2][j],wt);getchar();
	    //if(j==naxes[0]/2){fprintf(stderr,"%d: %g %g %g %g",j,xav[0][j],wt,xav[2][j],xav[0][j]/xav[2][j]);getchar();}
	  }
	}
	if(lp==3) break;
	//j=289;fprintf(stderr,"%d: %g\n",j,xav[1][j]);
	sgmax=0.;
	for(nf=0;nf<NFIBER;nf++){
	  i=LINE0+PSPEC*nf-PSPEC/2;
	  //fprintf(stderr,"lsft %02d %g %g\n",nf,xav[0][4800],xav[1][4800]);
	  for(ii=0;ii<PSPEC;ii++){
	    sg=lsft(zo[i+ii],xav,6,zi[i+ii]);
	    if(sg<0) break;
	    if(sgmax<sg) sgmax=sg;
	  }
	}
	if(sg<0) break;
	fprintf(stderr,"sgmax=%g\n",sgmax);
	for(i=0;i<naxes[1];i++){
	  for(j=0;j<naxes[0];j++){
	    if(zo[i][j]>sgmax*10||zo[i][j]<-sgmax*5) zo[i][j]=0./0.;
	    //if(j==2184&&i==233){fprintf(stderr,"%d %d: %g %g %g\n",j,i,zo[i][j],sgmax,yav[i]);}
	  }
	}
	for(i=0;i<naxes[1];i++){
	  for(j=0;j<naxes[0];j++){
	    if(isnan(zo[i][j])==0&&(
				    isnan(zo[i-1][j])||isnan(zo[i+1][j])||isnan(zo[i][j-1])||isnan(zo[i][j+1])
				    ||isnan(zo[i-2][j])||isnan(zo[i+2][j])||isnan(zo[i][j-2])||isnan(zo[i][j+2])
				    ||isnan(zo[i-1][j-1])||isnan(zo[i+1][j-1])||isnan(zo[i-1][j+1])||isnan(zo[i+1][j+1])
				    )) zo[i][j]=0.;
	  }
	}
	for(i=0;i<naxes[1];i++){
	  for(j=0;j<naxes[0];j++){
	    if(isnan(zo[i][j])) zo[i][j]=0.;
	  }
	}
	for(nf=0;nf<NFIBER;nf++){
	  i=LINE0+PSPEC*nf-PSPEC/2;
	  for(ii=0;ii<PSPEC;ii++){
	    for(j=0;j<naxes[0];j++){
	      zo[i+ii][j]+=zi[i+ii][j];
	    }
	  }
	}
      }
      if(sg<0){
	if(!dae) break;
	continue;
      }
      for(nf=0;nf<NFIBER;nf++){
	for(j=0;j<naxes[0]/48;j++) rfs[nf][j]=0.;
	ws=0.;n=0;
	for(j=0;j<naxes[0];j++){
	  if(isnan(sp[nf][j])==0){
	    rf[nf][j]=sp[nf][j]/xav[0][j];
	    ws+=rf[nf][j];
	    n++;
	  }
	}
	ws/=n;
	for(j=0;j<naxes[0];j++){
	  if(isnan(sp[nf][j])) rf[nf][j]=1.;
	  else rf[nf][j]/=ws;
	  if(isnan(rf[nf][j])==0) rfs[nf][j/48]+=rf[nf][j]/48;
	  //if(j==2450) fprintf(stderr,"%2d: %lf %lf %lf %lf\n",nf,ws,rf[nf][j],xav[0][j],fav[0][nf]);
	}
      }
      fprintf(stderr,"Write: double_y.dat ... ");
      fpt=fopen("double_y.dat","w");
      for(nf=0;nf<NFIBER;nf++) fprintf(fpt,"%g\n",fav[0][nf]);
      fclose(fpt);
      fprintf(stderr,"OK\n");
      for(nf=0;nf<NFIBER;nf++){
	for(j=naxes[0]/2-1;j>=0;j--){
	  if(isnan(rf[nf][j])) rf[nf][j]=rf[nf][j+1];
	}
	for(j=naxes[0]/2;j<naxes[0];j++){
	  if(isnan(rf[nf][j])) rf[nf][j]=rf[nf][j-1];
	}
	for(j=0;j<naxes[0];j++) sp[nf][j]=rf[nf][j]*rf0[nf][j];
	for(j=0;j<naxes[0];j++){
	  rf[nf][j]=med(sp[nf],j);
	  //if(j==naxes[0]/2) fprintf(stderr,"%d: %g\n",nf,rf[nf][j]);
	}
      }
      for(i=0;i<IMAX;i++){
	for(j=0;j<JMAX;j++) bn[i][j]=0.;
      }
      favs=0.;xct=0.;yct=0.;
      for(nf=0;nf<NFIBER;nf++){
	for(ii=0;ii<4;ii++){
	  for(jj=0;jj<4;jj++) bn[ic[nf]+ii][jc[nf]+jj]+=fav[1][nf]*kn2[ii][jj];
	}
	xct+=fav[1][nf]*xc[nf];
	yct+=fav[1][nf]*yc[nf];
	favs+=fav[1][nf];
      }
      xct/=favs;
      yct/=favs;
      zs=0.;zt=0.;r2s=0.;favmin=999.;favmax=-999.;
      for(nf=0;nf<NFIBER;nf++){
	r2=(xc[nf]-xct)*(xc[nf]-xct)+(yc[nf]-yct)*(yc[nf]-yct);
	r=exp(-r2/2.25);
	r2s+=fav[1][nf]*r2;
	//fprintf(stderr,"%g %g %g\n",r2,r2s,fav[1][nf]);
	zs+=fav[1][nf]*r;
	zt+=r;
	if(favmin>fav[1][nf]) favmin=fav[1][nf];
	if(favmax<fav[1][nf]) favmax=fav[1][nf];
      }
      r2s=sqrt(r2s/favs)*2;
      sig=0.;noi=0.;
      for(j=naxes[0]/4;j<naxes[0]/4*3;j++){
	sig+=xav[0][j];
	noi+=xav[2][j];
      }
      sig/=naxes[0]/2;
      noi/=naxes[0]/2;
      //sig=(zs-(favs-zs)*zt/(61-zt))*favav[1];
      //bkg=(favs-zs)/(61-zt)*favav[1];
      //fprintf(stderr,"Center / Sigma : %g %g / %g\n",xct,yct,r2s);
      //fprintf(stderr,"Fave / Signal / Backgr: %g / %g / %g\n",favav[0],sig,bkg);
      fprintf(stderr,"Center : %g %g\n",xct,yct);
      fprintf(stderr,"Sigma :  %g\n",r2s);
      fprintf(stderr,"Fave  : %g\n",favav[0]);
      fprintf(stderr,"S/N    : %g\n",sig/noi);
      fprintf(stderr,"Write: double_x.dat ... ");
      fpt=fopen("double_x.dat","w");
      //fprintf(fpt,"# %g %g %g %g %g %g\n",xct,yct,r2s,favav[0],sig,bkg);
      fprintf(fpt,"# %g %g %g %g %g\n",xct,yct,r2s,favav[0],sig/noi);
      r=CONVF*16/gn*1.98645e-9/exptime/1.04367e5/fabs(dw*10)*1e16;
      n=1;ew=0.;ewn=0.;
      for(j=0;j<naxes[0];j++){
	if(isnan(xav[0][j])) xav[0][j]=0.;
	if(isnan(xav[2][j])) xav[2][j]=0.;
	ww[j]=wc+(j-naxes[0]/2)*dw;
	zs=0.;
	for(i=LINE0-PSPEC/2;i<LINE0+PSPEC*NFIBER-PSPEC/2;i++) zs+=zo[i][j];
	zs*=r/ww[j];
	if(isnan(zs)) zs=0.;
	fprintf(fpt,"%g\t%g\t%g\t%g\t%g\t%g\n",ww[j],ww[j]*airindex(ww[j]),xav[0][j],xav[0][j]-xav[1][j],xav[2][j],zs);
	while(cw[n]/dw<ww[j]/dw){
	  if(n==nc-1) break;
	  n++;
	}
	if(cam){
	  if(fabs((ww[j]-656.28)/dw)<HWZI){
	    if(fabs((ww[j]-656.28)/dw)<HWZI/10) ewn+=xav[0][j];
	    ew+=xav[0][j];
	    ww[j]=0./0.;
	  }
	  else if(fabs((ww[j]-cw[n])/dw)<10 || fabs((ww[j]-cw[n-1])/dw)<10) ww[j]=0./0.;
	  else if(ww[j]<642 || ww[j]>667) ww[j]=0./0.;
	  else ww[j]-=wc;
	}
	else{
	  if(fabs((ww[j]-393.38)/dw)<HWZI || fabs((ww[j]-396.85)/dw)<HWZI){
	    if(fabs((ww[j]-393.38)/dw)<HWZI/10 || fabs((ww[j]-396.85)/dw)<HWZI/10) ewn+=xav[0][j];
	    ew+=xav[0][j];
	    ww[j]=0./0.;
	  }
	  else if(fabs((ww[j]-cw[n])/dw)<10 || fabs((ww[j]-cw[n-1])/dw)<10) ww[j]=0./0.;
	  else if(ww[j]<385 || ww[j]>405) ww[j]=0./0.;
	  else ww[j]-=wc;
	}
	//if(isnan(ww[j])==0){fprintf(stderr,"%d: %g",j,ww[j]);getchar();}
      }
      lsfit(ww,xav[0],1,3,aa);
      if(cam){
	ew=(ew/((aa[0]+aa[1]*(656.28-wc))*(2*HWZI-1))-1)*(2*HWZI-1)*fabs(dw);
	ewn=(ewn/((aa[0]+aa[1]*(656.28-wc))*(2*HWZI/10-1))-1)*(2*HWZI/10-1)*fabs(dw);
	n1=(656.28-wc)/dw+naxes[0]/2;
	n2=n1;
      }
      else{
	ew=(ew/((aa[0]+aa[1]*((393.38+396.85)/2-wc))*(2*HWZI-1)*2)-1)*(2*HWZI-1)*2*fabs(dw);
	ewn=(ewn/((aa[0]+aa[1]*((393.38+396.85)/2-wc))*(2*HWZI/10-1)*2)-1)*(2*HWZI/10-1)*2*fabs(dw);
	n1=(393.38-wc)/dw+naxes[0]/2;
	n2=(396.85-wc)/dw+naxes[0]/2;
      }
      fprintf(stderr,"F_cont=%g%+g*(w-%g), EW_line=%g / %g ",aa[0],aa[1],wc,ew,ewn);
      ymin=1e99;ymax=-1e99;
      for(j=naxes[0]/4;j<naxes[0]/4*3;j++){
	if(ymin>xav[0][j]) ymin=xav[0][j];
	if(ymax<xav[0][j]) ymax=xav[0][j];
      }
      wt=(ymin+ymax)/2;
      ymin+=ymin-wt;
      ymax+=ymax-wt;
      fclose(fpt);
      fprintf(stderr,"OK\n");
      if(k==0){
	strcpy(fnamet,fname[1]);
	fnamet[strlen(fnamet)-5]='\0';
      }
      else sprintf(fnamet,"%s_%03d",fname[1],k);
      if(strcmp(fnamet,"double")!=0){
	sprintf(buf,"cp double_x.dat %s_x.dat",fnamet);
	fprintf(stderr,"%s\n",buf);
	system(buf);
      }
      if(BMPZ1<0) bmpz1=favmin; else bmpz1=BMPZ1;
      if(BMPZ2>9.9) bmpz2=favmax; else bmpz2=BMPZ2;
      if(fabs(bmpz2-bmpz1)<0.02){bmpz1-=0.01;bmpz2+=0.01;}
      for(i=0;i<IMAX;i++){
	for(j=0;j<JMAX;j++){
	  n=(bn[i][j]-bmpz1)/(bmpz2-bmpz1)*256;
	  if(n<0) n=0;
	  else if(n>255) n=255;
	  gr_screen[i*JMAX+j]=c[n];
	}
      }
      gr_BMP("double_im.bmp",gr_screen,IMAX,JMAX);
      remove("double_rf.fits");
      fprintf(stderr,"Write double_rf.fits ... ");
      status=0;
      naxes[0]=LINEPIX2;naxes[1]=NFIBER;
      if (fits_create_file(&fp[1], "double_rf.fits", &status)) printerror(status);
      if (fits_create_img(fp[1], bitpix, naxis, naxes, &status)) printerror(status);
      fpixel=1;
      fprintf(stderr,"%ld x %ld format ",naxes[0],naxes[1]);
      for(i=0;i<naxes[1];i++){
	if(fits_write_img(fp[1], TFLOAT, fpixel, naxes[0], rf[i], &status)) printerror(status);
	fpixel+=naxes[0];
      }
      if (fits_close_file(fp[1], &status)) printerror(status);
      fprintf(stderr,"OK\n");
      for(i=0;i<NFIBER;i++){
	for(j=0;j<LINEPIX2/48;j++){
	  n=(rfs[i][j]-bmpz1)/(bmpz2-bmpz1)*256;
	  if(n<0) n=0;
	  else if(n>255) n=255;
	  gr_screen[i*LINEPIX2/48+j]=c[n];
	}
      }
      gr_BMP("double_rf.bmp",gr_screen,NFIBER,LINEPIX2/48);
      remove("double_rs.fits");
      fprintf(stderr,"Write double_rs.fits ... ");
      status=0;
      naxes[0]=LINEPIX2;naxes[1]=LINEPIY2;
      if (fits_create_file(&fp[1], "double_rs.fits", &status)) printerror(status);
      if (fits_create_img(fp[1], bitpix, naxis, naxes, &status)) printerror(status);
      fpixel=1;
      fprintf(stderr,"%ld x %ld format ",naxes[0],naxes[1]);
      for(i=0;i<naxes[1];i++){
	for(j=0;j<naxes[0];j++) zp[i][j]-=zo[i][j];
	if(fits_write_img(fp[1], TFLOAT, fpixel, naxes[0], zp[i], &status)) printerror(status);
	fpixel+=naxes[0];
      }
      if (fits_close_file(fp[1], &status)) printerror(status);
      fprintf(stderr,"OK\n");
    }
    else{
      for(i=0;i<naxes[1];i++){
	for(j=0;j<naxes[0];j++) zo[i][j]=zi[i][j];
      }
      if(k==0){
	strcpy(fnamet,fname[1]);
	fnamet[strlen(fnamet)-5]='\0';
      }
      else sprintf(fnamet,"%s_%03d",fname[1],k);
    }
    strcat(fnamet,".fits");
    remove(fnamet);
    fprintf(stderr,"Write %s ... ",fnamet);
    status=0;
    naxes[0]=LINEPIX2;naxes[1]=LINEPIY2;
    naxis=3;
    if (fits_create_file(&fp[1], fnamet, &status)) printerror(status);
    if (fits_create_img(fp[1], bitpix, naxis, naxes, &status)) printerror(status);
    if(argc>=4){
      if (fits_write_key(fp[1], TSTRING, "CTYPE1", "WAVELENGTH", "", &status)) printerror(status);
      if (fits_write_key(fp[1], TSTRING, "CTYPE2", "FIBERNUM", "", &status)) printerror(status);
      if (fits_write_key(fp[1], TSTRING, "CUNIT1", "nm", "", &status)) printerror(status);
      tmpval=naxes[0]/2+1.;
      if (fits_write_key(fp[1], TFLOAT, "CRPIX1", &tmpval, "", &status)) printerror(status);
      tmpval=wc;
      if (fits_write_key(fp[1], TFLOAT, "CRVAL1", &tmpval, "", &status)) printerror(status);
      tmpval=dw;
      if (fits_write_key(fp[1], TFLOAT, "CDELT1", &tmpval, "", &status)) printerror(status);
      if (fits_write_key(fp[1], TFLOAT, "CD1_1", &tmpval, "", &status)) printerror(status);
      tmpval=naxes[1]/2+1.;
      if (fits_write_key(fp[1], TFLOAT, "CRPIX2", &tmpval, "", &status)) printerror(status);
      tmpval=((float)NFIBER)/2+1;
      if (fits_write_key(fp[1], TFLOAT, "CRVAL2", &tmpval, "", &status)) printerror(status);
      tmpval=1./9;
      if (fits_write_key(fp[1], TFLOAT, "CDELT2", &tmpval, "", &status)) printerror(status);
      if (fits_write_key(fp[1], TFLOAT, "CD2_2", &tmpval, "", &status)) printerror(status);
    }
    fpixel=1;
    fprintf(stderr,"%ld x %ld x %ld format ",naxes[0],naxes[1],naxes[2]);
    for(i=0;i<naxes[1];i++){
      if(fits_write_img(fp[1], TFLOAT, fpixel, naxes[0], zo[i], &status)) printerror(status);
      fpixel+=naxes[0];
    }
    for(i=0;i<naxes[1];i++){
      if(fits_write_img(fp[1], TFLOAT, fpixel, naxes[0], zn[i], &status)) printerror(status);
      fpixel+=naxes[0];
    }
    if (fits_close_file(fp[1], &status)) printerror(status);
    fprintf(stderr,"OK\n");
    if(argc>=4){
      if(dae) strcpy(fnamet,namelist[0]->d_name);
      jm=strlen(fnamet);
      for(j=0;j<jm;j++){
	//fprintf(stderr,"%d: %s\n",j,fnamet);
	if(fnamet[j]=='_'||fnamet[j]=='^'){
	  jm+=2;
	  for(jj=jm;jj>j+1;jj--) fnamet[jj]=fnamet[jj-2];
	  fnamet[j]='\\';fnamet[j+1]='\\';
	j+=2;
	}
      }
      if(sft&&strcmp(object,"INSTFLAT")!=0) sprintf(fnamet+strlen(fnamet),"! need ANCHOR");
      sprintf(buf,"sed -e 's/FILENAME/%s/g' %sdouble.plt.template | sed -e 's/YMIN/%lf/g' | sed -e 's/YMAX/%lf/g' | sed -e 's/XCT/%+5.2lf/g' | sed -e 's/YCT/%+5.2lf/g' | sed -e 's/R2S/%4.2lf/g' | sed -e 's/FAV/%g/g' | sed -e 's/SN/%g/g' | sed -e 's/T0/%5d/g' | sed -e 's/N0/%d/g' | sed -e 's/T1/%5d/g' | sed -e 's/N1/%d/g' | sed -e 's/T2/%5d/g' | sed -e 's/N2/%d/g' | sed -e 's/T3/%5d/g' | sed -e 's/N3/%d/g' | sed -e 's/A0/%g/g' | sed -e 's/A1/%+g/g' | sed -e 's/WC/%g/g' | sed -e 's/S1/%d/g' | sed -e 's/E1/%d/g' | sed -e 's/S2/%d/g' | sed -e 's/E2/%d/g' | sed -e 's/S3/%d/g' | sed -e 's/E3/%d/g' | sed -e 's/S4/%d/g' | sed -e 's/E4/%d/g' > double.plt",fnamet,sdir,ymin,ymax,xct,yct,r2s,favav[0],sig/noi,tpix[0],npix[0],tpix[1],npix[1],tpix[2],npix[2],tpix[3],npix[3],aa[0],aa[1],wc,n1-HWZI,n1+HWZI,n2-HWZI,n2+HWZI,n1-HWZI/10,n1+HWZI/10,n2-HWZI/10,n2+HWZI/10);
      fprintf(stderr,"%s\n",buf);
      system(buf);
      sprintf(buf,"gnuplot < double.plt");
      fprintf(stderr,"%s\n",buf);
      system(buf);
      if(strcmp(object,"")&&strcmp(object,"NONE")&&strcmp(object,"FOCUSING")&&strcmp(object,"ANCHOR")&&strcmp(object,"DARK")&&strcmp(object,"BIAS")&&strcmp(object,"INSTFLAT")&&strcmp(object,"SKYFLAT")){
	if(cam) sprintf(fnamet,"%s/%s_H.log",dir,object);
	else sprintf(fnamet,"%s/%s_C.log",dir,object);
	fprintf(stderr,"Append: %s ... ",fnamet);
	fpt=fopen(fnamet,"a");
	fprintf(fpt,"%s %s %lf %lf %lf %lf %lf %d %d %d %d %lf %lf %g %lf %lf\n",date[0],jstend[0],xct,yct,r2s,favav[0],sig/noi,npix[0],npix[1],npix[2],npix[3],aa[0],aa[1],wc,ew,ewn);
	fclose(fpt);
	fprintf(stderr,"OK\n");
	jm=strlen(fnamet);
	for(j=0;j<jm;j++){
	  //fprintf(stderr,"%d: %s\n",j,fnamet);
	  if(fnamet[j]=='/'){
	    jm++;
	    for(jj=jm;jj>j;jj--) fnamet[jj]=fnamet[jj-1];
	    fnamet[j]='\\';
	    j++;
	  }
	}
	if(cam) sprintf(fnamet2,"%s_H.log",object);
	else sprintf(fnamet2,"%s_C.log",object);
	jm=strlen(fnamet2);
	for(j=0;j<jm;j++){
	  //fprintf(stderr,"%d: %s\n",j,fnamet2);
	  if(fnamet2[j]=='_'||fnamet2[j]=='^'){
	    jm+=2;
	    for(jj=jm;jj>j+1;jj--) fnamet2[jj]=fnamet2[jj-2];
	    fnamet2[j]='\\';fnamet2[j+1]='\\';
	    j+=2;
	  }
	}
	sprintf(buf,"sed -e 's/FNAME/%s/g' %sewlog.plt.template | sed -e 's/FILENAME/%s/g' > ewlog.plt",fnamet2,sdir,fnamet);
	fprintf(stderr,"%s\n",buf);
	system(buf);
	sprintf(buf,"gnuplot < ewlog.plt");
	fprintf(stderr,"%s\n",buf);
	system(buf);
	if(dae){
	  if(cam) sprintf(buf,"cp ewlog.png %s/%s_ewH.png",dir,object);
	  else sprintf(buf,"cp ewlog.png %s/%s_ewC.png",dir,object);
	  fprintf(stderr,"%s\n",buf);
	  system(buf);
	  if(cam) sprintf(buf,"cp ewlog.png /mnt/QL/inst2.png");
	  else sprintf(buf,"cp ewlog.png inst.png");
	  fprintf(stderr,"%s\n",buf);
	  system(buf);
	}
      }
      if(dae){
	sprintf(buf,"cp spec.png %s/%s",dir,namelist[0]->d_name);
	sprintf(buf+strlen(buf)-4,"png");
	fprintf(stderr,"%s\n",buf);
	system(buf);
	if(cam){
	  //sprintf(buf,"scp spec.png 192.168.1.160:/ssd/QL/spec2.png");
	  sprintf(buf,"cp spec.png /mnt/QL/spec2.png");
	  fprintf(stderr,"%s\n",buf);
	  system(buf);
	}
      }
    }
    if(!dae&&k==0) break;
  }
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
