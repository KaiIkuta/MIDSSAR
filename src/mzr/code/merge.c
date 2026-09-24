#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define NMAX 512
#define LINEPIX2 4848

int main(int argc,char *argv[])
{
  FILE *fp;
  char fname[NMAX][64],fnamet[64],buf[1024],dir[32];
  int n,ndat,j,jj,jm,ns,n1,n2;
  float x[NMAX],y[NMAX],s[NMAX],f[NMAX],sn[NMAX],w[2][LINEPIX2];
  static float fl[NMAX][LINEPIX2],flc[NMAX][LINEPIX2],no[NMAX][LINEPIX2],ef[NMAX][LINEPIX2];
  double xs,ys,ss,fs,fs2,wt,ws,fav,nav,sx,sy,sx2,sy2,sxy,a,b;
  static double fls[LINEPIX2],nos[LINEPIX2],flcs[LINEPIX2];
  
  if(argc!=2&&argc!=4){
    fprintf(stderr,"Usage : %s prefix [from# to#]\n",argv[0]);
    exit(0);
  }
  if(argc==4){
    sscanf(argv[2],"%d",&n1);
    sscanf(argv[3],"%d",&n2);
  }
  else{
    n1=1;n2=10000;
  }
  sscanf(argv[0],"%s",dir);
  dir[strlen(dir)-5]='\0';
  xs=0.;ys=0.;ss=0.;fs=0.;fs2=0.;
  for(j=0;j<LINEPIX2;j++){fls[j]=0.;flcs[j]=0.;nos[j]=0.;}
  n=0;ws=0.;ns=0;
  while(1){
    sprintf(fname[n],"%s_%03d_x.dat",argv[1],n+ns+n1);
    if(n+ns+n1>n2){
      fprintf(stderr,"Loop ends. Number of files: %d\n",n);
      break;
   }
    if((fp = fopen(fname[n],"r"))== NULL){
      if(ns>=20){
	fprintf(stderr,"Loop ends. Number of files: %d\n",n);
	break;
      }
      ns++;
      continue;
    }
    fprintf(stderr,"Read: %s ... ",fname[n]);
    fgets(buf,sizeof(buf),fp);
    if(sscanf(buf+1,"%f %f %f %f %f",&x[n],&y[n],&s[n],&f[n],&sn[n])!=5){
      fprintf(stderr,"Broken file (%s)\n",buf);
      fclose(fp);
      exit(1);
    }
    xs+=x[n]*f[n];
    ys+=y[n]*f[n];
    ss+=s[n]*f[n];
    fs+=f[n];
    fs2+=f[n]*f[n];
    for(j=0;j<LINEPIX2;j++){
      if(fgets(buf,sizeof(buf),fp)==NULL){
	fprintf(stderr,"Truncated file (%s @line %d)\n",buf,j+1);
	fclose(fp);
	exit(1);
      }
      if(sscanf(buf,"%f %f %f %f %f %f",&w[0][j],&w[1][j],&fl[n][j],&flc[n][j],&no[n][j],&ef[n][j])!=6){
	fprintf(stderr,"Broken file (%s @line %d)\n",buf,j+1);
	fclose(fp);
	exit(1);
      }
    }
    fav=0.;nav=0.;
    for(j=LINEPIX2/4;j<LINEPIX2/4*3;j++){
      fav+=fl[n][j];
      nav+=no[n][j];
    }
    fav/=LINEPIX2/2;
    nav/=LINEPIX2/2;
    wt=1.;//fav/nav/nav;
    ws+=wt;
    for(j=0;j<LINEPIX2;j++){
      fls[j]+=fl[n][j]*wt;
      flcs[j]+=flc[n][j]*wt;
      nos[j]+=no[n][j]*no[n][j]*wt*wt;
      //if(j==LINEPIX2/4){fprintf(stderr,"%g %g %g %g %g\n",fls[j],fl[n][j],f[n],nos[j],no[n][j]);}
    }
    fclose(fp);
    fprintf(stderr,"%f %f %f %f OK\n",x[n],y[n],s[n],f[n]);
    n++;
  }
  ndat=n;
  fprintf(stderr,"ndat=%d\n",ndat);
  if(ndat<2){fprintf(stderr,"More than 2 files are required.\n");exit(0);}
  xs/=fs;
  ys/=fs;
  ss/=fs;
  fs2/=fs;
  sprintf(fnamet,"%s_x.dat",argv[1]);
  fprintf(stderr,"Write: %s ... ",fnamet);
  fp=fopen(fnamet,"w");
  fprintf(fp,"# %lf %lf %lf %lf\n",xs,ys,ss,fs2);
  for(j=0;j<LINEPIX2;j++){
    //if(j==LINEPIX2/4) fprintf(stderr,"%g %g %g\n",fls[j],fs,nos[j]);
    fls[j]/=ws;
    flcs[j]/=ws;
    nos[j]=sqrt(nos[j])/ws;
    //if(j==LINEPIX2/4) fprintf(stderr,"%g %g %g\n",fls[j],fs,nos[j]);
    fprintf(fp,"%lf %lf %lf %lf %lf\n",w[0][j],w[1][j],fls[j],flcs[j],nos[j]);
  }
  fclose(fp);
  fprintf(stderr,"OK\n");
  for(n=0;n<ndat;n++){
    sx=0.;sy=0.;sx2=0.;sxy=0.;sy2=0.;
    for(j=LINEPIX2/4;j<LINEPIX2/4*3;j++){
      sx+=(j-LINEPIX2/2);
      sy+=fl[n][j]/fls[j];
      sx2+=(j-LINEPIX2/2)*(j-LINEPIX2/2);
      sxy+=(j-LINEPIX2/2)*fl[n][j]/fls[j];
      sy2+=fl[n][j]/fls[j]*fl[n][j]/fls[j];
    }
    sx/=LINEPIX2/2;
    sy/=LINEPIX2/2;
    sx2/=LINEPIX2/2;
    sxy/=LINEPIX2/2;
    sy2/=LINEPIX2/2;
    a=(sxy-sx*sy)/(sx2-sx*sx);
    b=(sx2*sy-sxy*sx)/(sx2-sx*sx);
    ws=sqrt(a*a*sx2+b*b+sy2+2*a*b*sx-2*a*sxy-2*b*sy);
    //printf("%s %g %g %g\n",fname[n],a,b,ws);
    fprintf(stderr,"Write: %s ... ",fname[n]);
    fp=fopen(fname[n],"w");
    fprintf(fp,"# %g %g %g %g %g %g\n",x[n],y[n],s[n],f[n],sn[n],sy/ws);
    for(j=0;j<LINEPIX2;j++){
      if(fls[j]>1e-5) fprintf(fp,"%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",w[0][j],w[1][j],fl[n][j],flc[n][j],no[n][j],ef[n][j],fl[n][j]/fls[j],no[n][j]/fls[j],fl[n][j]/fls[j]/sy,no[n][j]/fls[j]/sy);
      else fprintf(fp,"%g\t%g\t%g\t%g\t%g\t%g\t0\t0\t0\t0\n",w[0][j],w[1][j],fl[n][j],flc[n][j],no[n][j],ef[n][j]);
    }
    fclose(fp);
    fprintf(stderr,"OK\n");
    strcpy(fnamet,fname[n]);
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
    //fprintf(stderr,"%s\n",fnamet);
    sprintf(buf,"sed -e 's/FILENAME/%s/g' %smerge.plt.template | sed -e 's/FNAME/%s/g' > merge.plt",fnamet,dir,fname[n]);
    fprintf(stderr,"%s\n",buf);
    system(buf);
    sprintf(buf,"gnuplot < merge.plt");
    fprintf(stderr,"%s\n",buf);
    system(buf);
  }
  return 0;
}
