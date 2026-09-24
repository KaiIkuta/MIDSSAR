#define LINEPIX 2160
#define LINEPIY 1280
#define LINEPIX2 4848
#define LINEPIY2 720
#define NFIBER 61
#define LINE0 90
#define PSPEC 9
#define NX 1
#define NY 1
#define HW 4

void resize(float zi[LINEPIY][LINEPIX2],long naxes[2],int icen)
{
  static double zo[LINEPIY2][LINEPIX2],zw[LINEPIY2][LINEPIX2];
  int i,j,i2,j2,ii,jj,kn[3][3]={{4,2,4},{2,1,2},{4,2,4}};
  
  if(naxes[0]==LINEPIX&&naxes[1]==LINEPIY){
    for(i=0;i<LINEPIY2;i++){
      for(j=0;j<LINEPIX2;j++){zo[i][j]=0.;zw[i][j]=0.;}
    }
    for(i=icen-LINEPIY2/4;i<icen+LINEPIY2/4;i++){
      i2=(i-icen)*2+LINEPIY2/2;
      for(j=0;j<naxes[0];j++){
	j2=j*2-naxes[0]+LINEPIX2/2;
	for(ii=i2-1;ii<=i2+1;ii++){
	  if(ii<0) continue;
	  for(jj=j2-1;jj<=j2+1;jj++){
	    if(isnan(zi[i][j])==0){
	      zo[ii][jj]+=zi[i][j]/kn[ii-(i2-1)][jj-(j2-1)];
	      zw[ii][jj]+=1./kn[ii-(i2-1)][jj-(j2-1)];
	    }
	    //if(ii==LINEPIY2/2&&jj==LINEPIX2/2) fprintf(stderr,"%d %d %d %d: %g %d %g %g\n",jj,ii,j,i,zi[i][j],kn[ii-(i2-1)][jj-(j2-1)],zo[ii][jj],zw[ii][jj]);
	  }
	}
	//if(i==naxes[1]/2 && j==naxes[0]/2) fprintf(stderr,"%d %d %g %g %g %g %g %g %g %g %g\n",j2,i2,zo[i2-1][j2-1],zo[i2-1][j2],zo[i2-1][j2+1],zo[i2][j2-1],zo[i2][j2],zo[i2][j2+1],zo[i2+1][j2-1],zo[i2+1][j2],zo[i2+1][j2+1]);
      }
    }
    naxes[0]=LINEPIX2;
    naxes[1]=LINEPIY2;
    for(i=0;i<naxes[1];i++){
      for(j=0;j<naxes[0];j++){
	if(zw[i][j]>0.4) zi[i][j]=zo[i][j]/zw[i][j]/4;
	else zi[i][j]=0./0.;
      }
    }
  }
}

double tr1(double a,double b,double c,double d,double e,double f,int md)
{
  double s;
  if(md==1) s=(a/(a-c)*(d-b)+a/(a-e)*(b-f))*a/2;
  else if(md==2) s=(b/(b-d)*(a-c)+b/(b-f)*(e-a))*b/2;
  return (fabs(s));
}

double tr2(double a,double b,double c,double d)
{
  a=fabs(a);b=fabs(b);c=fabs(c);d=fabs(d);
  return((a*d-b*c)*(a*d-b*c)/(a+c)/(b+d)/2);
}

double rc1(double a,double b,double c,double d,double e,double f)
{
  return(fabs(a*d-b*c)/fabs(a-c)*fabs(a)/2+fabs(b*e-a*f)/fabs(b-f)*fabs(b)/2);
}

double rc2(double a,double b,double c,double d,double e,double f,double g,double h,int md)
{
  double s;
  if(md==1) s=(a+g)*(b-h)/2+a*a*(d-b)/(a-c)/2+g*g*(h-f)/(g-e)/2;
  else if(md==2) s=(b+d)*(a-c)/2+b*b*(g-a)/(b-h)/2+d*d*(c-e)/(d-f)/2;
  else if(md==3) s=(c+e)*(f-d)/2+c*c*(b-d)/(a-c)/2+e*e*(f-h)/(g-e)/2;
  else if(md==4) s=(f+h)*(e-g)/2+h*h*(a-g)/(b-h)/2+f*f*(e-c)/(d-f)/2;
  if(s<0.){fprintf(stderr,"rc2: %g %d\n(%g,%g)\n(%g,%g)\n(%g,%g)\n(%g,%g)\n%g %g %g\n",s,md,a,b,c,d,e,f,g,h,(c+e)*(f-d)/2,c*c*(b-d)/(a-c)/2,e*e*(f-h)/(g-e)/2);getchar();}
  return (s);
}

void warp(float zi[LINEPIY2][LINEPIX2],float zo[LINEPIY2][LINEPIX2],float zx[LINEPIY2][LINEPIX2],float zy[LINEPIY2][LINEPIX2])
{
  int i,j,ni[3],nj[3],i1,i2,i3,i4,j1,j2,j3,j4,iq,jq,is,js,ic,jc,ii1,ii2,ii3,ii4,jj1,jj2,jj3,jj4,nan;
  float fx1,fx2,fy1,fy2;
  double gx[3][3],gy[3][3],x1,x2,x3,x4,y1,y2,y3,y4,ff1,ff2,ff3,ff4,ds;

  for(i=0;i<LINEPIY2;i++){
    for(j=0;j<LINEPIX2;j++) zo[i][j]=0./0.;
  }
  for(i=0;i<LINEPIY2;i++){
    ni[1]=i;
    if(i==0) ni[0]=i; else ni[0]=i-1;
    if(i==LINEPIY2-1) ni[2]=i; else ni[2]=i+1;
    for(j=0;j<LINEPIX2;j++){
      if(isnan(zi[i][j])||isnan(zx[i][j])||isnan(zy[i][j])) continue;
      nj[1]=j;
      if(j==0) nj[0]=j; else nj[0]=j-1;
      if(j==LINEPIX2-1) nj[2]=j; else nj[2]=j+1;
      nan=0;
      for(iq=0;iq<3;iq++){
	for(jq=0;jq<3;jq++){
	  if(isnan(zy[ni[iq]][nj[jq]])) {gy[iq][jq]=zy[i][j]+ni[iq];nan++;} else gy[iq][jq]=zy[ni[iq]][nj[jq]]+ni[iq];
	  if(isnan(zx[ni[iq]][nj[jq]])) {gx[iq][jq]=zx[i][j]+nj[jq];nan++;} else gx[iq][jq]=zx[ni[iq]][nj[jq]]+nj[jq];
	  if(i==0&&iq==0)  gy[iq][jq]-=1.;
	  else if(i==LINEPIY2-1&&iq==2)  gy[iq][jq]+=1.;
	  if(j==0&&jq==0)  gx[iq][jq]-=1.;
	  else if(j==LINEPIX2-1&&jq==2)  gx[iq][jq]+=1.;
	}
      }
      if(nan) continue;
      for(iq=0;iq<3;iq++){
	for(jq=0;jq<3;jq++){
	  gx[iq][jq]=(gx[iq][jq]+gx[1][1])/2;
	  gy[iq][jq]=(gy[iq][jq]+gy[1][1])/2;
	}
      }
      for(iq=0;iq<2;iq++){
	for(jq=0;jq<2;jq++){
	  if((gx[iq][jq]>gx[iq+1][jq+1]&&gy[iq][jq]>gy[iq+1][jq+1])||(gx[iq+1][jq]>gx[iq][jq+1]&&gy[iq][jq+1]>gy[iq+1][jq])){
	    fprintf(stderr,"Crushed grid: %d %d %d %d %d %d\n",nj[0],nj[1],nj[2],ni[0],ni[1],ni[2]);
	    fprintf(stderr,"(%g,%g)(%g,%g)(%g,%g)(%g,%g)\n",gx[iq][jq],gy[iq][jq],gx[iq][jq+1],gy[iq][jq+1],gx[iq+1][jq],gy[iq+1][jq],gx[iq+1][jq+1],gy[iq+1][jq+1]);getchar();
	  }
	  for(is=0;is<NY;is++){
	    fy1=(float)is/NY;
	    fy2=(float)(is+1)/NY;
	    for(js=0;js<NX;js++){
	      fx1=(float)js/NX;
	      fx2=(float)(js+1)/NX;
	      y1=gy[iq][jq]*(1-fy2)*(1-fx2)+gy[iq][jq+1]*(1-fy2)*fx2+gy[iq+1][jq]*fy2*(1-fx2)+gy[iq+1][jq+1]*fy2*fx2+0.5;
	      x1=gx[iq][jq]*(1-fy2)*(1-fx2)+gx[iq][jq+1]*(1-fy2)*fx2+gx[iq+1][jq]*fy2*(1-fx2)+gx[iq+1][jq+1]*fy2*fx2+0.5;
	      y2=gy[iq][jq]*(1-fy2)*(1-fx1)+gy[iq][jq+1]*(1-fy2)*fx1+gy[iq+1][jq]*fy2*(1-fx1)+gy[iq+1][jq+1]*fy2*fx1+0.5;
	      x2=gx[iq][jq]*(1-fy2)*(1-fx1)+gx[iq][jq+1]*(1-fy2)*fx1+gx[iq+1][jq]*fy2*(1-fx1)+gx[iq+1][jq+1]*fy2*fx1+0.5;
	      y3=gy[iq][jq]*(1-fy1)*(1-fx1)+gy[iq][jq+1]*(1-fy1)*fx1+gy[iq+1][jq]*fy1*(1-fx1)+gy[iq+1][jq+1]*fy1*fx1+0.5;
	      x3=gx[iq][jq]*(1-fy1)*(1-fx1)+gx[iq][jq+1]*(1-fy1)*fx1+gx[iq+1][jq]*fy1*(1-fx1)+gx[iq+1][jq+1]*fy1*fx1+0.5;
	      y4=gy[iq][jq]*(1-fy1)*(1-fx2)+gy[iq][jq+1]*(1-fy1)*fx2+gy[iq+1][jq]*fy1*(1-fx2)+gy[iq+1][jq+1]*fy1*fx2+0.5;
	      x4=gx[iq][jq]*(1-fy1)*(1-fx2)+gx[iq][jq+1]*(1-fy1)*fx2+gx[iq+1][jq]*fy1*(1-fx2)+gx[iq+1][jq+1]*fy1*fx2+0.5;
	      i1=y1;i2=y2;i3=y3;i4=y4;
	      if(y1<0) i1--;
	      if(y2<0) i2--;
	      if(y3<0) i3--;
	      if(y4<0) i4--;
	      j1=x1;j2=x2;j3=x3;j4=x4;
	      if(x1<0) j1--;
	      if(x2<0) j2--;
	      if(x3<0) j3--;
	      if(x4<0) j4--;
	      if(i1-i3>1||i2-i4>1){
		fprintf(stderr,"NY is too small: (%d,%d)(%d,%d)(%d,%d)(%d,%d)\n",j1,i1,j2,i2,j3,i3,j4,i4);
		fprintf(stderr,"(%d %d) %g %g %g\n",j-1,i+1,zi[i+1][j-1],zx[i+1][j-1]-zx[i][j],zy[i+1][j-1]-zy[i][j]);
		fprintf(stderr,"(%d %d) %g %g %g\n",j  ,i+1,zi[i+1][j  ],zx[i+1][j  ]-zx[i][j],zy[i+1][j  ]-zy[i][j]);
		fprintf(stderr,"(%d %d) %g %g %g\n",j+1,i+1,zi[i+1][j+1],zx[i+1][j+1]-zx[i][j],zy[i+1][j+1]-zy[i][j]);
		fprintf(stderr,"(%d %d) %g %g %g\n",j-1,i  ,zi[i  ][j-1],zx[i  ][j-1]-zx[i][j],zy[i  ][j-1]-zy[i][j]);
		fprintf(stderr,"(%d %d) %g %g %g\n",j  ,i  ,zi[i  ][j  ],zx[i  ][j  ],zy[i][j]);
		fprintf(stderr,"(%d %d) %g %g %g\n",j+1,i  ,zi[i  ][j+1],zx[i  ][j+1]-zx[i][j],zy[i  ][j+1]-zy[i][j]);
		fprintf(stderr,"(%d %d) %g %g %g\n",j-1,i-1,zi[i-1][j-1],zx[i-1][j-1]-zx[i][j],zy[i-1][j-1]-zy[i][j]);
		fprintf(stderr,"(%d %d) %g %g %g\n",j  ,i-1,zi[i-1][j  ],zx[i-1][j  ]-zx[i][j],zy[i-1][j  ]-zy[i][j]);
		fprintf(stderr,"(%d %d) %g %g %g\n",j+1,i-1,zi[i-1][j+1],zx[i-1][j+1]-zx[i][j],zy[i-1][j+1]-zy[i][j]);
		exit(0);
	      }
	      if(j1-j3>1||j4-j2>1){
		fprintf(stderr,"NX is too small: (%d,%d)(%d,%d)(%d,%d)(%d,%d)\n",j1,i1,j2,i2,j3,i3,j4,i4);
		fprintf(stderr,"(%d %d) %g %g %g\n",j-1,i+1,zi[i+1][j-1],zx[i+1][j-1]-zx[i][j],zy[i+1][j-1]-zy[i][j]);
		fprintf(stderr,"(%d %d) %g %g %g\n",j  ,i+1,zi[i+1][j  ],zx[i+1][j  ]-zx[i][j],zy[i+1][j  ]-zy[i][j]);
		fprintf(stderr,"(%d %d) %g %g %g\n",j+1,i+1,zi[i+1][j+1],zx[i+1][j+1]-zx[i][j],zy[i+1][j+1]-zy[i][j]);
		fprintf(stderr,"(%d %d) %g %g %g\n",j-1,i  ,zi[i  ][j-1],zx[i  ][j-1]-zx[i][j],zy[i  ][j-1]-zy[i][j]);
		fprintf(stderr,"(%d %d) %g %g %g\n",j  ,i  ,zi[i  ][j  ],zx[i  ][j  ],zy[i][j]);
		fprintf(stderr,"(%d %d) %g %g %g\n",j+1,i  ,zi[i  ][j+1],zx[i  ][j+1]-zx[i][j],zy[i  ][j+1]-zy[i][j]);
		fprintf(stderr,"(%d %d) %g %g %g\n",j-1,i-1,zi[i-1][j-1],zx[i-1][j-1]-zx[i][j],zy[i-1][j-1]-zy[i][j]);
		fprintf(stderr,"(%d %d) %g %g %g\n",j  ,i-1,zi[i-1][j  ],zx[i-1][j  ]-zx[i][j],zy[i-1][j  ]-zy[i][j]);
		fprintf(stderr,"(%d %d) %g %g %g\n",j+1,i-1,zi[i-1][j+1],zx[i-1][j+1]-zx[i][j],zy[i-1][j+1]-zy[i][j]);
		exit(0);
	      }
	      if(x1+x2+x3+x4<-2.) jc=(x1+x2+x3+x4)/4-0.5;else jc=(x1+x2+x3+x4)/4+0.5;
	      if(y1+y2+y3+y4<-2.) ic=(y1+y2+y3+y4)/4-0.5;else ic=(y1+y2+y3+y4)/4+0.5;
	      x1-=jc;x2-=jc;x3-=jc;x4-=jc;y1-=ic;y2-=ic;y3-=ic;y4-=ic;
	      ds=(x1*y2-x2*y1+x2*y3-x3*y2+x3*y4-x4*y3+x4*y1-x1*y4)/2;
	      ii1=i1;ii2=i2;ii3=i3;ii4=i4;jj1=j1;jj2=j2;jj3=j3;jj4=j4;ff1=0.;ff2=0.;ff3=0.;ff4=0.;
	      //1234
	      if(i1==i2&&i2==i3&&i3==i4&&j1==j2&&j2==j3&&j3==j4){ff1=ds;}
	      //123-4
	      else if(i1==i2&&i2==i3&&i3==i4&&j1==j2&&j2==j3){ff4=tr1(x4,y4,x1,y1,x3,y3,1);ff1=ds-ff4;}
	      else if(i1==i2&&i2==i3&&j1==j2&&j2==j3&&j3==j4){ff4=tr1(x4,y4,x3,y3,x1,y1,2);ff1=ds-ff4;}
	      else if(i1==i2&&i2==i3&&j1==j2&&j2==j3){
		if(x3*y4-x4*y3<0.){jj1++;ff1=tr2(x4,y4,x1,y1)-tr2(x3,y3,x4,y4);ff4=tr1(x4,y4,x3,y3,x1,y1,2);ff2=ds-ff1-ff4;}
		else if(x4*y1-x1*y4<0.){ii3--;ff3=tr2(x3,y3,x4,y4)-tr2(x4,y4,x1,y1);ff4=tr1(x4,y4,x1,y1,x3,y3,1);ff2=ds-ff3-ff4;}
		else {jj1++;ii3--;ff3=tr2(x3,y3,x4,y4);ff1=tr2(x4,y4,x1,y1);ff4=rc1(x4,y4,x3,y3,x1,y1);ff2=ds-ff1-ff3-ff4;}
	      }
	      //124-3
	      else if(i1==i2&&i2==i3&&i3==i4&&j1==j2&&j2==j4){ff3=tr1(x3,y3,x2,y2,x4,y4,1);ff1=ds-ff3;}
	      else if(i1==i2&&i2==i4&&j1==j2&&j2==j3&&j3==j4){ff3=tr1(x3,y3,x4,y4,x2,y2,2);ff1=ds-ff3;}
	      else if(i1==i2&&i2==i4&&j1==j2&&j2==j4){
		if(x2*y3-x3*y2<0.){ii4--;ff4=tr2(x3,y3,x4,y4)-tr2(x2,y2,x3,y3);ff3=tr1(x3,y3,x2,y2,x4,y4,1);ff1=ds-ff3-ff4;}
		else if(x3*y4-x4*y3<0.){jj2--;ff2=tr2(x2,y2,x3,y3)-tr2(x3,y3,x4,y4);ff3=tr1(x3,y3,x4,y4,x2,y2,2);ff1=ds-ff2-ff3;}
		else {jj2--;ii4--;ff2=tr2(x2,y2,x3,y3);ff4=tr2(x3,y3,x4,y4);ff3=rc1(x3,y3,x4,y4,x2,y2);ff1=ds-ff2-ff3-ff4;}
	      }
	      //134-2
	      else if(i1==i2&&i2==i3&&i3==i4&&j1==j3&&j3==j4){ff2=tr1(x2,y2,x3,y3,x1,y1,1);ff1=ds-ff2;}
	      else if(i1==i3&&i3==i4&&j1==j2&&j2==j3&&j3==j4){ff2=tr1(x2,y2,x1,y1,x3,y3,2);ff1=ds-ff2;}
	      else if(i1==i3&&i3==i4&&j1==j3&&j3==j4){
		if(x1*y2-x2*y1<0.){jj3--;ff3=tr2(x2,y2,x3,y3)-tr2(x1,y1,x2,y2);ff2=tr1(x2,y2,x1,y1,x3,y3,2);ff4=ds-ff2-ff3;}
		else if(x2*y3-x3*y2<0.){ii1++;ff1=tr2(x1,y1,x2,y2)-tr2(x2,y2,x3,y3);ff2=tr1(x2,y2,x3,y3,x1,y1,1);ff4=ds-ff1-ff2;}
		else {ii1++;jj3--;ff1=tr2(x1,y1,x2,y2);ff3=tr2(x2,y2,x3,y3);ff2=rc1(x2,y2,x1,y1,x3,y3);ff4=ds-ff1-ff2-ff3;}
	      }
	      //234-1
	      else if(i1==i2&&i2==i3&&i3==i4&&j2==j3&&j3==j4){ff1=tr1(x1,y1,x4,y4,x2,y2,1);ff2=ds-ff1;}
	      else if(i2==i3&&i3==i4&&j1==j2&&j2==j3&&j3==j4){ff1=tr1(x1,y1,x2,y2,x4,y4,2);ff2=ds-ff1;}
	      else if(i2==i3&&i3==i4&&j2==j3&&j3==j4){
		if(x4*y1-x1*y4<0.){ii2++;ff2=tr2(x1,y1,x2,y2)-tr2(x4,y4,x1,y1);ff1=tr1(x1,y1,x4,y4,x2,y2,1);ff3=ds-ff1-ff2;}
		else if(x1*y2-x2*y1<0.){jj4++;ff4=tr2(x4,y4,x1,y1)-tr2(x1,y1,x2,y2);ff1=tr1(x1,y1,x2,y2,x4,y4,2);ff3=ds-ff1-ff4;}
		else {ii2++;jj4++;ff4=tr2(x4,y4,x1,y1);ff2=tr2(x1,y1,x2,y2);ff1=rc1(x1,y1,x2,y2,x4,y4);ff3=ds-ff1-ff2-ff4;}
	      }
	      //14-23
	      else if(i1==i2&&i2==i3&&i3==i4&&j1==j4&&j2==j3){ff1=rc2(x1,y1,x2,y2,x3,y3,x4,y4,1);ff3=rc2(x1,y1,x2,y2,x3,y3,x4,y4,3);}
	      else if(i1==i4&&i2==i3&&i2<i4&&j1==j2&&j2==j3&&j3==j4){ff1=rc2(x4,y4,x1,y1,x2,y2,x3,y3,2);ff3=rc2(x4,y4,x1,y1,x2,y2,x3,y3,4);}
	      else if(i1==i4&&i2==i3&&i2<i4&&j1==j4&&j2==j3){
		if(x1*y2-x2*y1<0.){ii4--;ff3=rc2(x1,y1,x2,y2,x3,y3,x4,y4,3);ff1=rc2(x4,y4,x1,y1,x2,y2,x3,y3,2);ff4=tr2(x3,y3,x4,y4)-tr2(x1,y1,x2,y2);}
		else if(x3*y4-x4*y3<0.){ii2++;ff1=rc2(x1,y1,x2,y2,x3,y3,x4,y4,1);ff3=rc2(x4,y4,x1,y1,x2,y2,x3,y3,4);ff2=tr2(x1,y1,x2,y2)-tr2(x3,y3,x4,y4);}
		else {ii2++;ii4--;ff2=tr2(x1,y1,x2,y2);ff4=tr2(x3,y3,x4,y4);ff1=rc2(x1,y1,x2,y2,x3,y3,x4,y4,1)-ff4;ff3=rc2(x1,y1,x2,y2,x3,y3,x4,y4,3)-ff2;}
	      }
	      else if(i1==i4&&i2==i3&&i1<i3&&j1==j2&&j2==j3&&j3==j4){ff2=rc2(x2,y2,x3,y3,x4,y4,x1,y1,2);ff4=rc2(x2,y2,x3,y3,x4,y4,x1,y1,4);}
	      else if(i1==i4&&i2==i3&&i1<i3&&j1==j4&&j2==j3){
		if(x1*y2-x2*y1<0.){ii3--;ff4=rc2(x1,y1,x2,y2,x3,y3,x4,y4,1);ff2=rc2(x2,y2,x3,y3,x4,y4,x1,y1,2);ff3=tr2(x3,y3,x4,y4)-tr2(x1,y1,x2,y2);}
		else if(x3*y4-x4*y3<0.){ii1++;ff2=rc2(x1,y1,x2,y2,x3,y3,x4,y4,3);ff4=rc2(x2,y2,x3,y3,x4,y4,x1,y1,4);ff1=tr2(x1,y1,x2,y2)-tr2(x3,y3,x4,y4);}
		else {ii1++;ii3--;ff1=tr2(x1,y1,x2,y2);ff3=tr2(x3,y3,x4,y4);ff4=rc2(x1,y1,x2,y2,x3,y3,x4,y4,1)-ff1;ff2=rc2(x1,y1,x2,y2,x3,y3,x4,y4,3)-ff3;}
	      }
	      //12-34
	      else if(i1==i2&&i3==i4&&j1==j2&&j2==j3&&j3==j4){ff2=rc2(x1,y1,x2,y2,x3,y3,x4,y4,2);ff4=rc2(x1,y1,x2,y2,x3,y3,x4,y4,4);}
	      else if(i1==i2&&i2==i3&&i3==i4&&j1==j2&&j3==j4&&j4<j2){ff1=rc2(x2,y2,x3,y3,x4,y4,x1,y1,1);ff3=rc2(x2,y2,x3,y3,x4,y4,x1,y1,3);}
	      else if(i1==i2&&i3==i4&&j1==j2&&j3==j4&&j4<j2){
		if(x2*y3-x3*y2<0.){jj4++;ff1=rc2(x1,y1,x2,y2,x3,y3,x4,y4,2);ff3=rc2(x2,y2,x3,y3,x4,y4,x1,y1,3);ff4=tr2(x4,y4,x1,y1)-tr2(x2,y2,x3,y3);}
		else if(x4*y1-x1*y4<0.){jj2--;ff3=rc2(x1,y1,x2,y2,x3,y3,x4,y4,4);ff1=rc2(x2,y2,x3,y3,x4,y4,x1,y1,1);ff2=tr2(x2,y2,x3,y3)-tr2(x4,y4,x1,y1);}
		else {jj2--;jj4++;ff2=tr2(x2,y2,x3,y3);ff4=tr2(x4,y4,x1,y1);ff1=rc2(x1,y1,x2,y2,x3,y3,x4,y4,2)-ff2;ff3=rc2(x1,y1,x2,y2,x3,y3,x4,y4,4)-ff4;}
	      }
	      else if(i1==i2&&i2==i3&&i3==i4&&j1==j2&&j3==j4&&j1<j3){ff4=rc2(x4,y4,x1,y1,x2,y2,x3,y3,1);ff2=rc2(x4,y4,x1,y1,x2,y2,x3,y3,3);}
	      else if(i1==i2&&i3==i4&&j1==j2&&j3==j4&&j1<j3){
		if(x2*y3-x3*y2<0.){jj1++;ff4=rc2(x1,y1,x2,y2,x3,y3,x4,y4,4);ff2=rc2(x4,y4,x1,y1,x2,y2,x3,y3,3);ff1=tr2(x4,y4,x1,y1)-tr2(x2,y2,x3,y3);}
		else if(x4*y1-x1*y4<0.){jj3--;ff2=rc2(x1,y1,x2,y2,x3,y3,x4,y4,2);ff4=rc2(x4,y4,x1,y1,x2,y2,x3,y3,1);ff3=tr2(x2,y2,x3,y3)-tr2(x4,y4,x1,y1);}
		else {jj1++;jj3--;ff3=tr2(x2,y2,x3,y3);ff1=tr2(x4,y4,x1,y1);ff2=rc2(x1,y1,x2,y2,x3,y3,x4,y4,2)-ff1;ff4=rc2(x1,y1,x2,y2,x3,y3,x4,y4,4)-ff3;}
	      }
	      //12-3-4
	      else if(i1==i2&&i3==i4&&j1==j2&&j2==j4){
		if(x2*y3-x3*y2>0.) {jj2--;ff2=tr2(x2,y2,x3,y3);ff3=rc1(x3,y3,x4,y4,x2,y2);ff4=rc1(x4,y4,x3,y3,x1,y1);ff1=rc2(x1,y1,x2,y2,x3,y3,x4,y4,2)-ff2;}
		else {ff1=rc2(x1,y1,x2,y2,x3,y3,x4,y4,2);ff3=tr1(x3,y3,x2,y2,x4,y4,1);ff4=rc2(x1,y1,x2,y2,x3,y3,x4,y4,4)-ff3;}
	      }
	      else if(i1==i2&&i2==i3&&j1==j2&&j2==j4){
		if(x3*y4-x4*y3>0.) {jj2--;ii3--;ff3=tr2(x3,y3,x4,y4);ff2=rc1(x3,y3,x2,y2,x4,y4);ff4=rc1(x4,y4,x3,y3,x1,y1);ff1=ds-ff2-ff3-ff4;}
		else {ff4=tr1(x4,y4,x3,y3,x1,y1,2);ff3=tr1(x3,y3,x2,y2,x4,y4,1);ff1=ds-ff3-ff4;}
	      }
	      else if(i1==i2&&i2==i3&&j1==j2&&j3==j4&&j4<j2){
		if(x4*y1-x1*y4>0.) {jj2--;ii3--;jj4++;ff4=tr2(x4,y4,x1,y1);ff2=rc1(x3,y3,x2,y2,x4,y4);ff3=rc1(x4,y4,x1,y1,x3,y3);ff1=rc2(x2,y2,x3,y3,x4,y4,x1,y1,1)-ff4;}
		else {ff4=tr1(x4,y4,x3,y3,x1,y1,2);ff3=rc1(x3,y3,x2,y2,x4,y4)-tr2(x4,y4,x1,y1);ff1=rc2(x2,y2,x3,y3,x4,y4,x1,y1,1);}
	      }
	      else if(i1==i2&&i3==i4&&j1==j2&&j2==j3){
		if(x4*y1-x1*y4>0.) {jj1++;ff1=tr2(x4,y4,x1,y1);ff4=rc1(x4,y4,x3,y3,x1,y1);ff3=rc1(x3,y3,x4,y4,x2,y2);ff2=rc2(x1,y1,x2,y2,x3,y3,x4,y4,2)-ff1;}
		else {ff2=rc2(x1,y1,x2,y2,x3,y3,x4,y4,2);ff4=tr1(x4,y4,x1,y1,x3,y3,1);ff3=rc1(x3,y3,x4,y4,x2,y2)-tr2(x4,y4,x1,y1);}
	      }
	      else if(i1==i2&&i2==i4&&j1==j2&&j2==j3){
		if(x3*y4-x4*y3>0.) {jj1++;ii4--;ff4=tr2(x3,y3,x4,y4);ff1=rc1(x4,y4,x1,y1,x3,y3);ff3=rc1(x3,y3,x4,y4,x2,y2);ff2=ds-ff1-ff3-ff4;}
		else {ff3=tr1(x3,y3,x4,y4,x2,y2,2);ff4=tr1(x4,y4,x1,y1,x3,y3,1);ff2=ds-ff3-ff4;}
	      }
	      else if(i1==i2&&i2==i4&&j1==j2&&j3==j4&&j3>j1){
		if(x2*y3-x3*y2>0.) {jj1++;ii4--;jj3--;ff3=tr2(x2,y2,x3,y3);ff1=rc1(x4,y4,x1,y1,x3,y3);ff4=rc1(x3,y3,x2,y2,x4,y4);ff2=rc2(x4,y4,x1,y1,x2,y2,x3,y3,3)-ff3;}
		else {ff3=tr1(x3,y3,x4,y4,x2,y2,2);ff4=rc1(x4,y4,x1,y1,x3,y3)-tr2(x2,y2,x3,y3);ff2=rc2(x4,y4,x1,y1,x2,y2,x3,y3,3);}
	      }
	      //23-1-4
	      else if(i1==i2&&i2==i3&&j1==j4&&j2==j3){
		if(x3*y4-x4*y3>0.) {ii3--;ff3=tr2(x3,y3,x4,y4);ff4=rc1(x4,y4,x3,y3,x1,y1);ff1=rc1(x1,y1,x2,y2,x4,y4);ff2=rc2(x1,y1,x2,y2,x3,y3,x4,y4,3)-ff3;}
		else {ff2=rc2(x1,y1,x2,y2,x3,y3,x4,y4,3);ff4=tr1(x4,y4,x3,y3,x1,y1,2);ff1=rc1(x1,y1,x2,y2,x4,y4)-tr2(x3,y3,x4,y4);}
	      }
	      else if(i1==i2&&i2==i3&&j2==j3&&j3==j4){
		if(x4*y1-x1*y4>0.) {jj4++;ii3--;ff4=tr2(x4,y4,x1,y1);ff3=rc1(x4,y4,x1,y1,x3,y3);ff1=rc1(x1,y1,x2,y2,x4,y4);ff2=ds-ff1-ff3-ff4;}
		else {ff1=tr1(x1,y1,x4,y4,x2,y2,1);ff4=tr1(x4,y4,x3,y3,x1,y1,2);ff2=ds-ff1-ff4;}
	      }
	      else if(i1==i4&&i2==i3&&j2==j3&&j3==j4&&i1<i3){
		if(x1*y2-x2*y1>0.) {jj4++;ii3--;ii1++;ff1=tr2(x1,y1,x2,y2);ff3=rc1(x4,y4,x1,y1,x3,y3);ff4=rc1(x1,y1,x4,y4,x2,y2);ff2=rc2(x2,y2,x3,y3,x4,y4,x1,y1,2)-ff1;}
		else {ff1=tr1(x1,y1,x4,y4,x2,y2,1);ff4=rc1(x4,y4,x1,y1,x3,y3)-tr2(x1,y1,x2,y2);ff2=rc2(x2,y2,x3,y3,x4,y4,x1,y1,2);}
	      }
	      else if(i2==i3&&i3==i4&&j1==j4&&j2==j3){
		if(x1*y2-x2*y1>0.) {ii2++;ff2=tr2(x1,y1,x2,y2);ff1=rc1(x1,y1,x2,y2,x4,y4);ff4=rc1(x4,y4,x3,y3,x1,y1);ff3=rc2(x1,y1,x2,y2,x3,y3,x4,y4,3)-ff2;}
		else {ff3=rc2(x1,y1,x2,y2,x3,y3,x4,y4,3);ff1=tr1(x1,y1,x2,y2,x4,y4,2);ff4=rc1(x4,y4,x3,y3,x1,y1)-tr2(x1,y1,x2,y2);}
	      }
	      else if(i2==i3&&i3==i4&&j1==j2&&j2==j3){
		if(x4*y1-x1*y4>0.) {jj1++;ii2++;ff1=tr2(x4,y4,x1,y1);ff2=rc1(x1,y1,x4,y4,x2,y2);ff4=rc1(x4,y4,x3,y3,x1,y1);ff3=ds-ff1-ff2-ff4;}
		else {ff4=tr1(x4,y4,x1,y1,x3,y3,1);ff1=tr1(x1,y1,x2,y2,x4,y4,2);ff3=ds-ff1-ff4;}
	      }
	      else if(i1==i4&&i2==i3&&j1==j2&&j2==j3&&i4>i2){
		if(x3*y4-x4*y3>0.) {jj1++;ii2++;ii4--;ff4=tr2(x3,y3,x4,y4);ff2=rc1(x1,y1,x4,y4,x2,y2);ff1=rc1(x4,y4,x1,y1,x3,y3);ff3=rc2(x4,y4,x1,y1,x2,y2,x3,y3,4)-ff4;}
		else {ff4=tr1(x4,y4,x1,y1,x3,y3,1);ff1=rc1(x1,y1,x4,y4,x2,y2)-tr2(x3,y3,x4,y4);ff3=rc2(x4,y4,x1,y1,x2,y2,x3,y3,4);}
	      }
	      //34-1-2
	      else if(i1==i2&&i3==i4&&j2==j3&&j3==j4){
		if(x4*y1-x1*y4>0.) {jj4++;ff4=tr2(x4,y4,x1,y1);ff1=rc1(x1,y1,x2,y2,x4,y4);ff2=rc1(x2,y2,x1,y1,x3,y3);ff3=rc2(x1,y1,x2,y2,x3,y3,x4,y4,4)-ff4;}
		else {ff3=rc2(x1,y1,x2,y2,x3,y3,x4,y4,4);ff1=tr1(x1,y1,x4,y4,x2,y2,1);ff2=rc1(x2,y2,x1,y1,x3,y3)-tr2(x4,y4,x1,y1);}
	      }
	      else if(i1==i3&&i3==i4&&j2==j3&&j3==j4){
		if(x1*y2-x2*y1>0.) {jj4++;ii1++;ff1=tr2(x1,y1,x2,y2);ff4=rc1(x1,y1,x4,y4,x2,y2);ff2=rc1(x2,y2,x1,y1,x3,y3);ff3=ds-ff1-ff2-ff4;}
		else {ff2=tr1(x2,y2,x1,y1,x3,y3,2);ff1=tr1(x1,y1,x4,y4,x2,y2,1);ff3=ds-ff1-ff2;}
	      }
	      else if(i1==i3&&i3==i4&&j1==j2&&j3==j4&&j2>j4){
		if(x2*y3-x3*y2>0.) {jj4++;ii1++;jj2--;ff2=tr2(x2,y2,x3,y3);ff4=rc1(x1,y1,x4,y4,x2,y2);ff1=rc1(x2,y2,x3,y3,x1,y1);ff3=rc2(x2,y2,x3,y3,x4,y4,x1,y1,3)-ff2;}
		else {ff2=tr1(x2,y2,x1,y1,x3,y3,2);ff1=rc1(x1,y1,x4,y4,x2,y2)-tr2(x2,y2,x3,y3);ff3=rc2(x2,y2,x3,y3,x4,y4,x1,y1,3);}
	      }
	      else if(i1==i2&&i3==i4&&j1==j3&&j3==j4){
		if(x2*y3-x3*y2>0.) {jj3--;ff3=tr2(x2,y2,x3,y3);ff2=rc1(x2,y2,x1,y1,x3,y3);ff1=rc1(x1,y1,x2,y2,x4,y4);ff4=rc2(x1,y1,x2,y2,x3,y3,x4,y4,4)-ff3;}
		else {ff4=rc2(x1,y1,x2,y2,x3,y3,x4,y4,4);ff2=tr1(x2,y2,x3,y3,x1,y1,1);ff1=rc1(x1,y1,x2,y2,x4,y4)-tr2(x2,y2,x3,y3);}
	      }
	      else if(i2==i3&&i3==i4&&j1==j3&&j3==j4){
		if(x1*y2-x2*y1>0.) {jj3--;ii2++;ff2=tr2(x1,y1,x2,y2);ff3=rc1(x2,y2,x3,y3,x1,y1);ff1=rc1(x1,y1,x2,y2,x4,y4);ff4=ds-ff1-ff2-ff3;}
		else {ff1=tr1(x1,y1,x2,y2,x4,y4,2);ff2=tr1(x2,y2,x3,y3,x1,y1,1);ff4=ds-ff1-ff2;}
	      }
	      else if(i2==i3&&i3==i4&&j1==j2&&j3==j4&&j1<j3){
		if(x4*y1-x1*y4>0.) {jj3--;ii2++;jj1++;ff1=tr2(x4,y4,x1,y1);ff3=rc1(x2,y2,x3,y3,x1,y1);ff2=rc1(x1,y1,x4,y4,x2,y2);ff4=rc2(x4,y4,x1,y1,x2,y2,x3,y3,1)-ff1;}
		else {ff1=tr1(x1,y1,x2,y2,x4,y4,2);ff2=rc1(x2,y2,x3,y3,x1,y1)-tr2(x4,y4,x1,y1);ff4=rc2(x4,y4,x1,y1,x2,y2,x3,y3,1);}
	      }
	      //14-2-3
	      else if(i1==i3&&i3==i4&&j1==j4&&j2==j3){
		if(x1*y2-x2*y1>0.) {ii1++;ff1=tr2(x1,y1,x2,y2);ff2=rc1(x2,y2,x1,y1,x3,y3);ff3=rc1(x3,y3,x4,y4,x2,y2);ff4=rc2(x1,y1,x2,y2,x3,y3,x4,y4,1)-ff1;}
		else {ff4=rc2(x1,y1,x2,y2,x3,y3,x4,y4,1);ff2=tr1(x2,y2,x1,y1,x3,y3,2);ff3=rc1(x3,y3,x4,y4,x2,y2)-tr2(x1,y1,x2,y2);}
	      }
	      else if(i1==i3&&i3==i4&&j1==j2&&j2==j4){
		if(x2*y3-x3*y2>0.) {jj2--;ii1++;ff2=tr2(x2,y2,x3,y3);ff1=rc1(x2,y2,x3,y3,x1,y1);ff3=rc1(x3,y3,x4,y4,x2,y2);ff4=ds-ff1-ff2-ff3;}
		else {ff3=tr1(x3,y3,x2,y2,x4,y4,1);ff2=tr1(x2,y2,x1,y1,x3,y3,2);ff4=ds-ff2-ff3;}
	      }
	      else if(i1==i4&&i2==i3&&j1==j2&&j2==j4&&i3>i1){
		if(x3*y4-x4*y3>0.) {jj2--;ii1++;i3--;ff3=tr2(x3,y3,x4,y4);ff1=rc1(x2,y2,x3,y3,x1,y1);ff2=rc1(x3,y3,x2,y2,x4,y4);ff4=rc2(x2,y2,x3,y3,x4,y4,x1,y1,4)-ff3;}
		else {ff3=tr1(x3,y3,x2,y2,x4,y4,1);ff2=rc1(x2,y2,x3,y3,x1,y1)-tr2(x3,y3,x4,y4);ff4=rc2(x2,y2,x3,y3,x4,y4,x1,y1,4);}
	      }
	      else if(i1==i2&&i2==i4&&j1==j4&&j2==j3){
		if(x3*y4-x4*y3>0.) {ii4--;ff4=tr2(x3,y3,x4,y4);ff3=rc1(x3,y3,x4,y4,x2,y2);ff2=rc1(x2,y2,x1,y1,x3,y3);ff1=rc2(x1,y1,x2,y2,x3,y3,x4,y4,1)-ff4;}
		else {ff1=rc2(x1,y1,x2,y2,x3,y3,x4,y4,1);ff3=tr1(x3,y3,x4,y4,x2,y2,2);ff2=rc1(x2,y2,x1,y1,x3,y3)-tr2(x3,y3,x4,y4);}
	      }
	      else if(i1==i2&&i2==i4&&j1==j3&&j3==j4){
		if(x2*y3-x3*y2>0.) {jj3--;ii4--;ff3=tr2(x2,y2,x3,y3);ff4=rc1(x3,y3,x2,y2,x4,y4);ff2=rc1(x2,y2,x1,y1,x3,y3);ff1=ds-ff2-ff3-ff4;}
		else {ff2=tr1(x2,y2,x3,y3,x1,y1,1);ff3=tr1(x3,y3,x4,y4,x2,y2,2);ff1=ds-ff2-ff3;}
	      }
	      else if(i1==i4&&i2==i3&&j1==j3&&j3==j4&&i2<i4){
		if(x1*y2-x2*y1>0.) {jj3--;ii4--;ii2++;ff2=tr2(x1,y1,x2,y2);ff4=rc1(x3,y3,x2,y2,x4,y4);ff3=rc1(x2,y2,x3,y3,x1,y1);ff1=rc2(x4,y4,x1,y1,x2,y2,x3,y3,2)-ff2;}
		else {ff2=tr1(x2,y2,x3,y3,x1,y1,1);ff3=rc1(x3,y3,x2,y2,x4,y4)-tr2(x1,y1,x2,y2);ff1=rc2(x4,y4,x1,y1,x2,y2,x3,y3,2);}
	      }
	      //13-2-4
	      else if(i1==i2&& i2==i3&&j1==j3&&j3==j4){ff2=tr1(x2,y2,x3,y3,x1,y1,1);ff4=tr1(x4,y4,x3,y3,x1,y1,2);ff1=ds-ff2-ff4;}
	      else if(i1==i3&& i3==i4&&j1==j2&&j2==j3){ff2=tr1(x2,y2,x1,y1,x3,y3,2);ff4=tr1(x4,y4,x1,y1,x3,y3,1);ff3=ds-ff2-ff4;}
	      //24-1-3
	      else if(i1==i2&& i2==i4&&j2==j3&&j3==j4){ff1=tr1(x1,y1,x4,y4,x2,y2,1);ff3=tr1(x3,y3,x4,y4,x2,y2,2);ff2=ds-ff1-ff3;}
	      else if(i2==i3&& i3==i4&&j1==j2&&j2==j4){ff1=tr1(x1,y1,x2,y2,x4,y4,2);ff3=tr1(x3,y3,x2,y2,x4,y4,1);ff4=ds-ff1-ff3;}
	      //1-2-3-4
	      else if(i1==i2&&i3==i4&&j1==j4&&j2==j3){
		ff1=rc1(x1,y1,x2,y2,x4,y4);ff2=rc1(x2,y2,x1,y1,x3,y3);ff3=rc1(x3,y3,x4,y4,x2,y2);ff4=rc1(x4,y4,x3,y3,x1,y1);
	      }
	      else{
		fprintf(stderr,"Undefined configuration: %d %d: (%d,%d)(%d,%d)(%d,%d)(%d,%d)(%g,%g)(%g,%g)(%g,%g)(%g,%g)\n",j,i,j1,i1,j2,i2,j3,i3,j4,i4,x1,y1,x2,y2,x3,y3,x4,y4);
		exit(0);
	      }
	      if(isinf(ff1)||isinf(ff2)||isinf(ff3)||isinf(ff4)){
		fprintf(stderr,"Inf factor: %d %d: (%d,%d)(%d,%d)(%d,%d)(%d,%d)(%g,%g)(%g,%g)(%g,%g)(%g,%g) %g %g %g %g\n",j,i,j1,i1,j2,i2,j3,i3,j4,i4,x1,y1,x2,y2,x3,y3,x4,y4,ff1,ff2,ff3,ff4);getchar();
	      }
	      if(fabs(ff1+ff2+ff3+ff4-ds)>1e-6||ff1<0.||ff2<0.||ff3<0.||ff4<0.){
		fprintf(stderr,"Something Wrong...: %d %d: (%d,%d)(%d,%d)(%d,%d)(%d,%d)(%g,%g)(%g,%g)(%g,%g)(%g,%g) %g %g %g %g %g %g\n",j,i,j1,i1,j2,i2,j3,i3,j4,i4,x1,y1,x2,y2,x3,y3,x4,y4,ff1,ff2,ff3,ff4,ff1+ff2+ff3+ff4,ds);getchar();
	      }
	      if(ii1<0) ii1=0; else if(ii1>LINEPIY2-1) ii1=LINEPIY2-1;
	      if(ii2<0) ii2=0; else if(ii2>LINEPIY2-1) ii2=LINEPIY2-1;
	      if(ii3<0) ii3=0; else if(ii3>LINEPIY2-1) ii3=LINEPIY2-1;
	      if(ii4<0) ii4=0; else if(ii4>LINEPIY2-1) ii4=LINEPIY2-1;
	      if(jj1<0) jj1=0; else if(jj1>LINEPIX2-1) jj1=LINEPIX2-1;
	      if(jj2<0) jj2=0; else if(jj2>LINEPIX2-1) jj2=LINEPIX2-1;
	      if(jj3<0) jj3=0; else if(jj3>LINEPIX2-1) jj3=LINEPIX2-1;
	      if(jj4<0) jj4=0; else if(jj4>LINEPIX2-1) jj4=LINEPIX2-1;
	      if(isnan(zo[ii1][jj1])) zo[ii1][jj1]=0.;
	      if(isnan(zo[ii2][jj2])) zo[ii2][jj2]=0.;
	      if(isnan(zo[ii3][jj3])) zo[ii3][jj3]=0.;
	      if(isnan(zo[ii4][jj4])) zo[ii4][jj4]=0.;
	      zo[ii1][jj1]+=zi[i][j]/NX/NY/4*ff1/ds;
	      zo[ii2][jj2]+=zi[i][j]/NX/NY/4*ff2/ds;
	      zo[ii3][jj3]+=zi[i][j]/NX/NY/4*ff3/ds;
	      zo[ii4][jj4]+=zi[i][j]/NX/NY/4*ff4/ds;
	    }
	  }
	}
      }
    }
  }
  /*
  int nzi,nzo;
  double szi,szo;
  szi=0.;nzi=0;
  szo=0.;nzo=0;
  for(i=0;i<LINEPIY2;i++){
    for(j=0;j<LINEPIX2;j++){
      if(isnan(zi[i][j])==0&&isnan(zx[i][j])==0&&isnan(zy[i][j])==0){
	szi+=zi[i][j];
	nzi++;
      }
      if(isnan(zo[i][j])==0){
	szo+=zo[i][j];
	nzo++;
	if(isnan(szo)){fprintf(stderr,"zo[%d][%d]=%g szo=%g nzo=%d\n",i,j,zo[i][j],szo,nzo);getchar();}
      }
      //zi[i][j]=zo[i][j];
    }
  }
  nzi=LINEPIX2*LINEPIY2-nzi;nzo=LINEPIX2*LINEPIY2-nzo;
  fprintf(stderr,"Sum : %lf => %lf\n",szi,szo);
  fprintf(stderr,"NAN : %d => %d\n",nzi,nzo);
  */
}
