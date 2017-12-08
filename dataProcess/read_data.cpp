static char help[] = "Solves a linear system in parallel with KSP.\n\
Input parameters include:\n\
  -random_exact_sol : use a random exact solution vector\n\
  -view_exact_sol   : write exact solution vector to stdout\n\
  -m <mesh_x>       : number of mesh points in x-direction\n\
  -n <mesh_n>       : number of mesh points in y-direction\n\n";

#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <petscksp.h>
using namespace std;
// ndat+1=(npsi-1)*(n2th-1)
// n2th=2*(nthe-1)=2*(101-1)=200

double psia,psmin,psmax,qmin,qmax,q0;
double alphar=0;
double prho=1.;
double arho=1.;


double rhom(double psival){
    double rsq,rho_re;
    rsq=(psival-psmin)/(psia-psmin);
    rho_re=(1.00000-alphar*rsq);
    return rho_re;
}

double rhomp(double psival){
    double rsq,rho_re;
    rsq=(psival-psmin)/(psia-psmin);
    rho_re=-arho*alphar*prho/(psia-psmin);
    return rho_re;
}

void interp1d3l(double x1,double x2,double x3,double x4,double y1,double y2,double y3,double y4,double y,double &ans){
  double d1 = (y1-y2)*(y1-y3)*(y1-y4);
  double d2 = (y2-y1)*(y2-y3)*(y2-y4);
  double d3 = (y3-y1)*(y3-y2)*(y3-y4);
  double d4 = (y4-y1)*(y4-y2)*(y4-y3);
  ans = x1*(y-y2)*(y-y3)*(y-y4)/d1 + x2*(y-y1)*(y-y3)*(y-y4)/d2 + x3*(y-y1)*(y-y2)*(y-y4)/d3 + x4*(y-y1)*(y-y2)*(y-y3)/d4;
  return ;
}


int main(int argc,char **args){
  FILE           *fp1,*fp2,*fp3,*fp4,*fp5,*fp6,*fp7,*fp8,*fp9,*fp10,*fp;
  PetscErrorCode ierr;
  PetscInitialize(&argc,&args,(char*)0,help);

  int npsi=111,nthe = 101;
  int n2th=2*(nthe-1),nphi = 8,ndat = (npsi-1)*(n2th-1),ndat12=(npsi-1)*(nthe+4),ndat34=(npsi-1)*(nthe+4);
  double pi = acos(-1);
  double xplmin,xplmax,zplmax,aguess,xzero,xmag,xmaj,xzmax,xatpi,xofset,aratio,bzero,curtotal,curnorm;
  double thst[n2th+5+1];
  double psst[n2th+5+1][npsi+1],xxst[n2th+5+1][npsi+1],zzst[n2th+5+1][npsi+1],tpst[n2th+5+1][npsi+1],tst[n2th+5+1][npsi+1],rst[n2th+5+1][npsi+1];
  double bxst[n2th+5+1][npsi+1],bxdxst[n2th+5+1][npsi+1],bxdzst[n2th+5+1][npsi+1],bzst[n2th+5+1][npsi+1],bzdxst[n2th+5+1][npsi+1],bzdzst[n2th+5+1][npsi+1];
  double byst[n2th+5+1][npsi+1],bydxst[n2th+5+1][npsi+1],cyst[n2th+5+1][npsi+1],uyst[n2th+5+1][npsi+1],uydxst[n2th+5+1][npsi+1],uydzst[n2th+5+1][npsi+1];
  double cxst[n2th+5+1][npsi+1],czst[n2th+5+1][npsi+1];
  int ipi;
  double ps_NOVA[ndat+1],xx_NOVA[ndat+1],zz_NOVA[ndat+1],bx_NOVA[ndat+1],bz_NOVA[ndat+1],bxdx_NOVA[ndat+1],bzdx_NOVA[ndat+1],bxdz_NOVA[ndat+1],bzdz_NOVA[ndat+1],th_NOVA[ndat+1];
  double pt_NOVA[ndat+1],rh_NOVA[ndat+1];
  double by_NOVA[ndat+1],bydx_NOVA[ndat+1],bydz_NOVA[ndat+1],pdx_NOVA[ndat+1],pdz_NOVA[ndat+1],cx_NOVA[ndat+1],cz_NOVA[ndat+1];
  double cy_NOVA[ndat+1],uy_NOVA[ndat+1];
  double psival_NOVA[npsi+1],q_NOVA[npsi+1],qp_NOVA[npsi+1],p_NOVA[npsi+1],pp_NOVA[npsi+1],g_NOVA[npsi+1],gp_NOVA[npsi+1],f_NOVA[npsi+1];
  double fp_NOVA[npsi+1],fb_NOVA[npsi+1],fbp_NOVA[npsi+1],omrot_NOVA[npsi+1],omprot_NOVA[npsi+1];
  double theta_NOVA[ndat+1],r_NOVA[ndat+1];
  for(int jt = 1;jt <= n2th+5;jt++){
    thst[jt]=pi*(jt-3)/(nthe-1);
  }
  ipi=nthe+2;

  ifstream filereader;
  filereader.open("psi_xz.dat",ios::in);
  string linestring;
  getline(filereader,linestring);
  string xplmin_str,xplmax_str,zplmax_str,aguess_str,xzero_str,xmag_str,xmaj_str,xzmax_str,xatpi_str,xofset_str,aratio_str,bzero_str,curtotal_str,curnorm_str;
  filereader >> xplmin_str >> xplmax_str >> zplmax_str >> aguess_str >> xzero_str >> xmag_str >> xmaj_str;
  filereader >> xzmax_str >> xatpi_str >> xofset_str >> aratio_str >> bzero_str >> curtotal_str >> curnorm_str;
  const char *xplmin_cstr  = xplmin_str.c_str();  xplmin = atof(xplmin_cstr);
  const char *xplmax_cstr = xplmax_str.c_str();  xplmax = atof(xplmax_cstr);
  const char *zplmax_cstr  = zplmax_str.c_str();  zplmax = atof(zplmax_cstr);
  const char *aguess_cstr = aguess_str.c_str();  aguess = atof(aguess_cstr);
  const char *xzero_cstr = xzero_str.c_str();  xzero = atof(xzero_cstr);
  const char *xmag_cstr  = xmag_str.c_str();  xmag = atof(xmag_cstr);
  const char *xmaj_cstr = xmaj_str.c_str();  xmaj = atof(xmaj_cstr);
  const char *xzmax_cstr  = xzmax_str.c_str();  xzmax = atof(xzmax_cstr);
  const char *xatpi_cstr = xatpi_str.c_str();  xatpi = atof(xatpi_cstr);
  const char *xofset_cstr  = xofset_str.c_str();  xofset = atof(xofset_cstr);
  const char *aratio_cstr = aratio_str.c_str();  aratio = atof(aratio_cstr);
  const char *bzero_cstr  = bzero_str.c_str();  bzero = atof(bzero_cstr);
  const char *curtotal_cstr = curtotal_str.c_str();  curtotal = atof(curtotal_cstr);
  const char *curnorm_cstr  = curnorm_str.c_str();  curnorm = atof(curnorm_cstr);
//=================== 归一化处理 ==================
  double aa = aguess;
  double b0 = 1;
  xzero=xzero/aa; //环坐标的轴
  double xmg=xmag/aa; //环磁场的轴
  double xmin=xplmin/aa;
  double xmax=xplmax/aa;
  double zmax=zplmax/aa;
  double zmin=-zmax/aa;
  double zmg=0.0;
  double cIp=curnorm*xzero*b0;
  int js,jt;
  double epsilon=1./aratio;
  string js_str,jt_str,psst_str,xxst_str,zzst_str,bxst_str,bxdxst_str,bxdzst_str,bzst_str,bzdxst_str,bzdzst_str;
  for(int j = 2;j <= npsi;j++){
    for(int i = 3;i <= nthe+1;i++){
      filereader >> js_str >> jt_str >> psst_str >> xxst_str >> zzst_str >> bxst_str;
      filereader >> bxdxst_str >> bxdzst_str >> bzst_str >> bzdxst_str >> bzdzst_str;
      const char *js_cstr = js_str.c_str();  js = atof(js_cstr);
      const char *jt_cstr = jt_str.c_str();  jt = atof(jt_cstr);
      const char *psst_cstr = psst_str.c_str();  psst[i][j] = atof(psst_cstr);
      const char *xxst_cstr = xxst_str.c_str();  xxst[i][j] = atof(xxst_cstr);
      const char *zzst_cstr = zzst_str.c_str();  zzst[i][j] = atof(zzst_cstr);
      const char *bxst_cstr = bxst_str.c_str();  bxst[i][j] = atof(bxst_cstr);
      const char *bxdxst_cstr = bxdxst_str.c_str();  bxdxst[i][j] = atof(bxdxst_cstr);
      const char *bxdzst_cstr = bxdzst_str.c_str();  bxdzst[i][j] = atof(bxdzst_cstr);
      const char *bzst_cstr = bzst_str.c_str();  bzst[i][j] = atof(bzst_cstr);
      const char *bzdxst_cstr = bzdxst_str.c_str();  bzdxst[i][j] = atof(bzdxst_cstr);
      const char *bzdzst_cstr = bzdzst_str.c_str();  bzdzst[i][j] = atof(bzdzst_cstr);
      xxst[i][j]=xxst[i][j]/aa;
      zzst[i][j]=zzst[i][j]/aa;
      psst[i][j]=psst[i][j]/(b0*aa*aa);
      bxst[i][j]=bxst[i][j]/b0;
      bzst[i][j]=bzst[i][j]/b0;
      bxdxst[i][j]=bxdxst[i][j]/(b0/aa);
      bxdzst[i][j]=bxdzst[i][j]/(b0/aa);
      bzdxst[i][j]=bzdxst[i][j]/(b0/aa);
      bzdzst[i][j]=bzdzst[i][j]/(b0/aa);
      // ======================= 划分的转换 =========================
      // ndat+1=(npsi-1)*(n2th-1)
      // n2th=2*(nthe-1)=2*(101-1)=200
      int jd=(j-2)*(n2th-1)+i-2;
      tst[i][j]=atan2(zzst[i][j],xxst[i][j]-xmg);
      rst[i][j]=sqrt(zzst[i][j]*zzst[i][j]+(xxst[i][j]-xmg)*(xxst[i][j]-xmg));

      th_NOVA[jd]=thst[i];
      xx_NOVA[jd]=xxst[i][j];
      zz_NOVA[jd]=zzst[i][j];
      ps_NOVA[jd]=psst[i][j];
      bx_NOVA[jd]=bxst[i][j];
      bz_NOVA[jd]=bzst[i][j];
      bxdx_NOVA[jd]=bxdxst[i][j];
      bxdz_NOVA[jd]=bxdzst[i][j];
      bzdx_NOVA[jd]=bzdxst[i][j];
      bzdz_NOVA[jd]=bzdzst[i][j];

      if(i >= 3){
          int im=2*nthe+2-(i-2);
          xxst[im][j]=xxst[i][j];
          zzst[im][j]=-zzst[i][j];
          psst[im][j]=psst[i][j];
          bxst[im][j]=-bxst[i][j];
          bzst[im][j]=bzst[i][j];
          bxdxst[im][j]=-bxdxst[i][j];
          bxdzst[im][j]=bxdzst[i][j];
          bzdxst[im][j]=bzdxst[i][j];
          bzdzst[im][j]=-bzdzst[i][j];
          tst[im][j]=2*pi-tst[i][j];
          rst[im][j]=rst[i][j];

          int jdm=(j-2)*(n2th-1)+im-3;
          th_NOVA[jdm]=thst[im];
          xx_NOVA[jdm]=xxst[im][j];
          zz_NOVA[jdm]=zzst[im][j];
          ps_NOVA[jdm]=psst[im][j];
          bx_NOVA[jdm]=bxst[im][j];
          bz_NOVA[jdm]=bzst[im][j];
          bxdx_NOVA[jdm]=bxdxst[im][j];
          bxdz_NOVA[jdm]=bxdzst[im][j];
          bzdx_NOVA[jdm]=bzdxst[im][j];
          bzdz_NOVA[jdm]=bzdzst[im][j];
      }
    }
//更改下标，便于差分操作。
    xxst[1][j]=xxst[1+n2th][j];
    zzst[1][j]=zzst[1+n2th][j];
    psst[1][j]=psst[1+n2th][j];
    bxst[1][j]=bxst[1+n2th][j];
    bzst[1][j]=bzst[1+n2th][j];
    bxdxst[1][j]=bxdxst[1+n2th][j];
    bzdxst[1][j]=bzdxst[1+n2th][j];
    bxdzst[1][j]=bxdzst[1+n2th][j];
    bzdzst[1][j]=bzdzst[1+n2th][j];
    tst[1][j]=tst[1+n2th][j]-2*pi;
    rst[1][j]=rst[1+n2th][j];

    xxst[2][j]=xxst[2+n2th][j];
    zzst[2][j]=zzst[2+n2th][j];
    psst[2][j]=psst[2+n2th][j];
    bxst[2][j]=bxst[2+n2th][j];
    bzst[2][j]=bzst[2+n2th][j];
    bxdxst[2][j]=bxdxst[2+n2th][j];
    bzdxst[2][j]=bzdxst[2+n2th][j];
    bxdzst[2][j]=bxdzst[2+n2th][j];
    bzdzst[2][j]=bzdzst[2+n2th][j];

    tst[2][j]=tst[2+n2th][j]-2*pi;
    rst[2][j]=rst[2+n2th][j];

    xxst[3+n2th][j]=xxst[3][j];
    zzst[3+n2th][j]=zzst[3][j];
    psst[3+n2th][j]=psst[3][j];
    bxst[3+n2th][j]=bxst[3][j];
    bzst[3+n2th][j]=bzst[3][j];
    bxdxst[3+n2th][j]=bxdxst[3][j];
    bzdxst[3+n2th][j]=bzdxst[3][j];
    bxdzst[3+n2th][j]=bxdzst[3][j];
    bzdzst[3+n2th][j]=bzdzst[3][j];
    tst[3+n2th][j]=tst[3][j]+2*pi;
    rst[3+n2th][j]=rst[3][j];

    xxst[4+n2th][j]=xxst[4][j];
    zzst[4+n2th][j]=zzst[4][j];
    psst[4+n2th][j]=psst[4][j];
    bxst[4+n2th][j]=bxst[4][j];
    bzst[4+n2th][j]=bzst[4][j];
    bxdxst[4+n2th][j]=bxdxst[4][j];
    bzdxst[4+n2th][j]=bzdxst[4][j];
    bxdzst[4+n2th][j]=bxdzst[4][j];
    bzdzst[4+n2th][j]=bzdzst[4][j];
    tst[4+n2th][j]=tst[4][j]+2*pi;
    rst[4+n2th][j]=rst[4][j];

    xxst[5+n2th][j]=xxst[5][j];
    zzst[5+n2th][j]=zzst[5][j];
    psst[5+n2th][j]=psst[5][j];
    bxst[5+n2th][j]=bxst[5][j];
    bzst[5+n2th][j]=bzst[5][j];
    bxdxst[5+n2th][j]=bxdxst[5][j];
    bzdxst[5+n2th][j]=bzdxst[5][j];
    bxdzst[5+n2th][j]=bxdzst[5][j];
    bzdzst[5+n2th][j]=bzdzst[5][j];
    tst[5+n2th][j]=tst[5][j]+2*pi;
    rst[5+n2th][j]=rst[5][j];

  //在ipi处补充赋值
    zzst[ipi][j]=0;
    interp1d3l(xxst[ipi-2][j],xxst[ipi-1][j],xxst[ipi+1][j],xxst[ipi+2][j],zzst[ipi-2][j],zzst[ipi-1][j],zzst[ipi+1][j],zzst[ipi+2][j],zzst[ipi][j],xxst[ipi][j]);
    interp1d3l(psst[ipi-2][j],psst[ipi-1][j],psst[ipi+1][j],psst[ipi+2][j],zzst[ipi-2][j],zzst[ipi-1][j],zzst[ipi+1][j],zzst[ipi+2][j],zzst[ipi][j],psst[ipi][j]);
    interp1d3l(bxst[ipi-2][j],bxst[ipi-1][j],bxst[ipi+1][j],bxst[ipi+2][j],zzst[ipi-2][j],zzst[ipi-1][j],zzst[ipi+1][j],zzst[ipi+2][j],zzst[ipi][j],bxst[ipi][j]);
    interp1d3l(bzst[ipi-2][j],bzst[ipi-1][j],bzst[ipi+1][j],bzst[ipi+2][j],zzst[ipi-2][j],zzst[ipi-1][j],zzst[ipi+1][j],zzst[ipi+2][j],zzst[ipi][j],bzst[ipi][j]);
    interp1d3l(bxdxst[ipi-2][j],bxdxst[ipi-1][j],bxdxst[ipi+1][j],bxdxst[ipi+2][j],zzst[ipi-2][j],zzst[ipi-1][j],zzst[ipi+1][j],zzst[ipi+2][j],zzst[ipi][j],bxdxst[ipi][j]);
    interp1d3l(bxdzst[ipi-2][j],bxdzst[ipi-1][j],bxdzst[ipi+1][j],bxdzst[ipi+2][j],zzst[ipi-2][j],zzst[ipi-1][j],zzst[ipi+1][j],zzst[ipi+2][j],zzst[ipi][j],bxdzst[ipi][j]);
    interp1d3l(bzdxst[ipi-2][j],bzdxst[ipi-1][j],bzdxst[ipi+1][j],bzdxst[ipi+2][j],zzst[ipi-2][j],zzst[ipi-1][j],zzst[ipi+1][j],zzst[ipi+2][j],zzst[ipi][j],bzdxst[ipi][j]);
    interp1d3l(bzdzst[ipi-2][j],bzdzst[ipi-1][j],bzdzst[ipi+1][j],bzdzst[ipi+2][j],zzst[ipi-2][j],zzst[ipi-1][j],zzst[ipi+1][j],zzst[ipi+2][j],zzst[ipi][j],bzdzst[ipi][j]);
    tst[ipi][j]=pi;
    rst[ipi][j]=fabs(xxst[ipi][j]-xmg);
  }
  filereader.close();
  filereader.open("q_p_g.dat",ios::in);
  int jj;
  for(int j = 1;j <= npsi;j++){
    string jj_str,psival_NOVA_str,q_NOVA_str,qp_NOVA_str,p_NOVA_str,pp_NOVA_str,g_NOVA_str;
    string gp_NOVA_str,f_NOVA_str,fp_NOVA_str,fb_NOVA_str,fbp_NOVA_str,omrot_NOVA_str,omprot_NOVA_str;
    filereader >> jj_str >> psival_NOVA_str >> q_NOVA_str >> qp_NOVA_str >> p_NOVA_str >> pp_NOVA_str >> g_NOVA_str;
    filereader >> gp_NOVA_str >> f_NOVA_str >> fp_NOVA_str >> fb_NOVA_str >> fbp_NOVA_str >> omrot_NOVA_str >> omprot_NOVA_str;
    const char *jj_cstr = jj_str.c_str();  jj = atoi(jj_cstr);
    const char *psival_NOVA_cstr = psival_NOVA_str.c_str();  psival_NOVA[j] = atof(psival_NOVA_cstr);
    const char *q_NOVA_cstr = q_NOVA_str.c_str();  q_NOVA[j] = atof(q_NOVA_cstr);
    const char *qp_NOVA_cstr = qp_NOVA_str.c_str();  qp_NOVA[j] = atof(qp_NOVA_cstr);
    const char *p_NOVA_cstr = p_NOVA_str.c_str();  p_NOVA[j] = atof(p_NOVA_cstr);
    const char *pp_NOVA_cstr = pp_NOVA_str.c_str();  pp_NOVA[j] = atof(pp_NOVA_cstr);
    const char *g_NOVA_cstr = g_NOVA_str.c_str();  g_NOVA[j] = atof(g_NOVA_cstr);
    const char *gp_NOVA_cstr = gp_NOVA_str.c_str();  gp_NOVA[j] = atof(gp_NOVA_cstr);
    const char *f_NOVA_cstr = f_NOVA_str.c_str();  f_NOVA[j] = atof(f_NOVA_cstr);
    const char *fp_NOVA_cstr = fp_NOVA_str.c_str();  fp_NOVA[j] = atof(fp_NOVA_cstr);
    const char *fb_NOVA_cstr = fb_NOVA_str.c_str();  fb_NOVA[j] = atof(fb_NOVA_cstr);
    const char *fbp_NOVA_cstr = fbp_NOVA_str.c_str();  fbp_NOVA[j] = atof(fbp_NOVA_cstr);
    const char *omrot_NOVA_cstr = omrot_NOVA_str.c_str();  omrot_NOVA[j] = atof(omrot_NOVA_cstr);
    const char *omprot_NOVA_cstr = omprot_NOVA_str.c_str();  omprot_NOVA[j] = atof(omprot_NOVA_cstr);
    psival_NOVA[j]=psival_NOVA[j]/(b0*aa*aa);
    p_NOVA[j]=p_NOVA[j]/(b0*b0);
    g_NOVA[j]=g_NOVA[j]/b0;
    qp_NOVA[j]=qp_NOVA[j]*(b0*aa*aa);
    pp_NOVA[j]=pp_NOVA[j]/(b0*b0)*(b0*aa*aa);
    gp_NOVA[j]=gp_NOVA[j]/b0*(b0*aa*aa);
    omrot_NOVA[j]=omrot_NOVA[j]/(b0*aa);
    omprot_NOVA[j]=omprot_NOVA[j]/(b0*aa)*(b0*aa*aa);
  }
  filereader.close();
  for(int i = 1;i <= n2th+5;i++){
    xxst[i][1]=xmg;
    zzst[i][1]=0;
    psst[i][1]=psival_NOVA[1];
    bxst[i][1]=0;
    bzst[i][1]=0;
    tst[i][1]=tst[i][2];
    rst[i][1]=0;
  }

  psia =psival_NOVA[npsi];
  psmin=*min_element(ps_NOVA+1,ps_NOVA+ndat+1);
  psmax=*max_element(ps_NOVA+1,ps_NOVA+ndat+1);
  qmin=*min_element(q_NOVA+1,q_NOVA+ndat+1);
  qmax=*max_element(q_NOVA+1,q_NOVA+ndat+1);
  q0=q_NOVA[1];

  interp1d3l(bxdxst[ipi][3],bxdxst[ipi][2],bxdxst[3][2],bxdxst[3][3],xxst[ipi][3],xxst[ipi][2],xxst[3][2],xxst[3][3],xmg,bxdxst[3][1]);
  interp1d3l(bxdzst[ipi][3],bxdzst[ipi][2],bxdzst[3][2],bxdzst[3][3],xxst[ipi][3],xxst[ipi][2],xxst[3][2],xxst[3][3],xmg,bxdzst[3][1]);
  interp1d3l(bzdxst[ipi][3],bzdxst[ipi][2],bzdxst[3][2],bzdxst[3][3],xxst[ipi][3],xxst[ipi][2],xxst[3][2],xxst[3][3],xmg,bzdxst[3][1]);
  interp1d3l(bzdzst[ipi][3],bzdzst[ipi][2],bzdzst[3][2],bzdzst[3][3],xxst[ipi][3],xxst[ipi][2],xxst[3][2],xxst[3][3],xmg,bzdzst[3][1]);
  for(int i = 1;i <= n2th+5;i++){
    bxdxst[i][1]=bxdxst[3][1];
    bxdzst[i][1]=bxdzst[3][1];
    bzdxst[i][1]=bzdxst[3][1];
    bzdzst[i][1]=bzdzst[3][1];
  }

  char fname[] = "data.dat";
  PetscFOpen(PETSC_COMM_SELF,fname,"w",&fp);
  if (!fp) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot open file");
  int jd;
  double temp1,temp2,temp3;
  for(int j = 1;j <= npsi;j++){
    for(int i = 1;i <= n2th+5;i++){
      byst[i][j]=xzero*g_NOVA[j]/xxst[i][j];
      // bydxst[i][j]=-xzero*(gp_NOVA[j]*bzst[i][j]+g_NOVA[j]/xxst[i][j]*xxst[i][j]);
      // bydzst[i][j]=xzero*gp_NOVA[j]*bxst[i][j];
      // pdxst[i][j]=-pp_NOVA[j]*bzst[i][j]*xxst[i][j];
      // pdzst[i][j]=pp_NOVA[j]*bxst[i][j]*xxst[i][j];

      cxst[i][j]=-xzero*gp_NOVA[j]*bxst[i][j];
      czst[i][j]=-xzero*gp_NOVA[j]*bzst[i][j];
      cyst[i][j]=bxdzst[i][j]-bzdxst[i][j];

      uyst[i][j]=omrot_NOVA[j]*xxst[i][j];
      // uydxst[i][j]=-omprot_NOVA[j]*bzst[i][j]*xxst[i][j]*xxst[i][j]+omrot_NOVA[j];
      // uydzst[i][j]= omprot_NOVA[j]*bxst[i][j]*xxst[i][j]*xxst[i][j];

      if(j >=2 && i >= 3 && i <= n2th+2 && i != ipi){
        if(i <= ipi) jd=(j-2)*(n2th-1)+i-2;
        if(i >= ipi) jd=(j-2)*(n2th-1)+i-3;

        by_NOVA[jd]=byst[i][j];
        // bydx_NOVA[jd]=bydxst[i][j];
        // bydz_NOVA[jd]=bydzst[i][j];
        // pdx_NOVA[jd]=pdxst[i][j];
        // pdz_NOVA[jd]=pdzst[i][j];
        cy_NOVA[jd]=cyst[i][j];
        cx_NOVA[jd]=cxst[i][j];
        cz_NOVA[jd]=czst[i][j];

        uy_NOVA[jd]=uyst[i][j];
        // uydx_NOVA[jd]=uydxst[i][j];
        // uydz_NOVA[jd]=uydzst[i][j];

        rh_NOVA[jd]=rhom(ps_NOVA[jd]);
        // rhdx_NOVA[jd]=-rhomp(ps_NOVA[jd])*bzst[i][j]*xxst[i][j];
        // rhdz_NOVA[jd]=rhomp(ps_NOVA[jd])*bxst[i][j]*xxst[i][j];

        pt_NOVA[jd]=p_NOVA[j];
        // ptdx_NOVA[jd]=pdx_NOVA[jd];
        // ptdz_NOVA[jd]=pdz_NOVA[jd];

        r_NOVA[jd] = rst[i][j];
        //坐标变换
        temp1 = bx_NOVA[jd]*cos(th_NOVA[jd])-bz_NOVA[jd]/r_NOVA[jd]*sin(th_NOVA[jd]);
        temp2 = bx_NOVA[jd]*sin(th_NOVA[jd])+bz_NOVA[jd]/r_NOVA[jd]*cos(th_NOVA[jd]);
        temp3 = by_NOVA[jd];
        PetscFPrintf(PETSC_COMM_SELF,fp,"%10e ",(double)temp1);
        PetscFPrintf(PETSC_COMM_SELF,fp,"%10e ",(double)temp2);
        PetscFPrintf(PETSC_COMM_SELF,fp,"%10e ",(double)temp3);
        PetscFPrintf(PETSC_COMM_SELF,fp,"%10e ",0);
        PetscFPrintf(PETSC_COMM_SELF,fp,"%10e ",0);
        PetscFPrintf(PETSC_COMM_SELF,fp,"%10e ",(double)uy_NOVA[jd]);
        PetscFPrintf(PETSC_COMM_SELF,fp,"%10e ",(double)rh_NOVA[jd]);
        PetscFPrintf(PETSC_COMM_SELF,fp,"%10e ",(double)pt_NOVA[jd]);
        PetscFPrintf(PETSC_COMM_SELF,fp,"%10e ",(double)r_NOVA[jd]);
        PetscFPrintf(PETSC_COMM_SELF,fp,"%10e ",(double)th_NOVA[jd]);
        PetscFPrintf(PETSC_COMM_SELF,fp,"%10e ",(double)xx_NOVA[jd]);
        PetscFPrintf(PETSC_COMM_SELF,fp,"%10e ",(double)zz_NOVA[jd]);
        PetscFPrintf(PETSC_COMM_SELF,fp,"\n");
      }
    }
  }
  PetscFClose(PETSC_COMM_SELF,fp);

  

  // char fname1[] = "rho_nova.dat";
  // char fname2[] = "P_nova.dat";
  // char fname3[] = "Br_nova.dat";
  // char fname4[] = "Bt_nova.dat";
  // char fname5[] = "Bp_nova.dat";
  // char fname6[] = "vr_nova.dat";
  // char fname7[] = "vt_nova.dat";
  // char fname8[] = "vp_nova.dat";
  // char fname9[] = "r_nova.dat";
  // char fname10[] = "theta_nova.dat";
  
  // PetscFOpen(PETSC_COMM_SELF,fname1,"w",&fp1);
  // PetscFOpen(PETSC_COMM_SELF,fname2,"w",&fp2);
  // PetscFOpen(PETSC_COMM_SELF,fname3,"w",&fp3);
  // PetscFOpen(PETSC_COMM_SELF,fname4,"w",&fp4);
  // PetscFOpen(PETSC_COMM_SELF,fname5,"w",&fp5);
  // PetscFOpen(PETSC_COMM_SELF,fname6,"w",&fp6);
  // PetscFOpen(PETSC_COMM_SELF,fname7,"w",&fp7);
  // PetscFOpen(PETSC_COMM_SELF,fname8,"w",&fp8);
  // PetscFOpen(PETSC_COMM_SELF,fname9,"w",&fp9);
  // PetscFOpen(PETSC_COMM_SELF,fname10,"w",&fp10);
  
  // if (!fp1) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot open file");
  // if (!fp2) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot open file");
  // if (!fp3) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot open file");
  // if (!fp4) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot open file");
  // if (!fp5) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot open file");
  // if (!fp6) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot open file");
  // if (!fp7) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot open file");
  // if (!fp8) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot open file");
  // if (!fp9) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot open file");
  // if (!fp10) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_USER,"Cannot open file");
  
  // for(int j = 1;j <= npsi-1;j++){
  //   for(int i = 1;i <= n2th-1;i++){
  //     for(int k = 1;k <= nphi;k++){
  //       PetscFPrintf(PETSC_COMM_SELF,fp1,"%10e ",(double)rho[i][j][k]);
  //       PetscFPrintf(PETSC_COMM_SELF,fp2,"%10e ",(double)P[i][j][k]);
  //       PetscFPrintf(PETSC_COMM_SELF,fp3,"%10e ",(double)Br[i][j][k]);
  //       PetscFPrintf(PETSC_COMM_SELF,fp4,"%10e ",(double)Bt[i][j][k]);
  //       PetscFPrintf(PETSC_COMM_SELF,fp5,"%10e ",(double)Bp[i][j][k]);
  //       PetscFPrintf(PETSC_COMM_SELF,fp6,"%10e ",(double)vr[i][j][k]);
  //       PetscFPrintf(PETSC_COMM_SELF,fp7,"%10e ",(double)vt[i][j][k]);
  //       PetscFPrintf(PETSC_COMM_SELF,fp8,"%10e ",(double)vp[i][j][k]);
  //       PetscFPrintf(PETSC_COMM_SELF,fp9,"%10e ",(double)r[i][j][k]);
  //       PetscFPrintf(PETSC_COMM_SELF,fp10,"%10e ",(double)theta[i][j][k]);
  //     }
  //   }
  //   PetscFPrintf(PETSC_COMM_SELF,fp1,"\n");
  //   PetscFPrintf(PETSC_COMM_SELF,fp2,"\n");
  //   PetscFPrintf(PETSC_COMM_SELF,fp3,"\n");
  //   PetscFPrintf(PETSC_COMM_SELF,fp4,"\n");
  //   PetscFPrintf(PETSC_COMM_SELF,fp5,"\n");
  //   PetscFPrintf(PETSC_COMM_SELF,fp6,"\n");
  //   PetscFPrintf(PETSC_COMM_SELF,fp7,"\n");
  //   PetscFPrintf(PETSC_COMM_SELF,fp8,"\n");
  //   PetscFPrintf(PETSC_COMM_SELF,fp9,"\n");
  //   PetscFPrintf(PETSC_COMM_SELF,fp10,"\n");
  // }
  // PetscFClose(PETSC_COMM_SELF,fp1);
  // PetscFClose(PETSC_COMM_SELF,fp2);
  // PetscFClose(PETSC_COMM_SELF,fp3);
  // PetscFClose(PETSC_COMM_SELF,fp4);
  // PetscFClose(PETSC_COMM_SELF,fp5);
  // PetscFClose(PETSC_COMM_SELF,fp6);
  // PetscFClose(PETSC_COMM_SELF,fp7);
  // PetscFClose(PETSC_COMM_SELF,fp8);
  // PetscFClose(PETSC_COMM_SELF,fp9);
  // PetscFClose(PETSC_COMM_SELF,fp10);
  // if(nrank.eq.0) then
  // do j=1,npsi
  // do i=3,nthe+1
  // fffst(i,j,1)=cyst[i][j]*bzst[i][j]-czst[i][j]*byst[i][j]-pdxst[i][j]-uyst[i][j]*uydxst[i][j]+uyst[i][j]**2/xxst[i][j]
  // fffst(i,j,2)=czst[i][j]*bxst[i][j]-cxst[i][j]*bzst[i][j]
  // fffst(i,j,3)=cxst[i][j]*byst[i][j]-cyst[i][j]*bxst[i][j]-pdzst[i][j]-uyst[i][j]*uydzst[i][j]
  // enddo
  // enddo
  //
  // open(unit=101,file='fst_NOVA.dat',status='unknown',form='formatted')
  // write(101,300)(((fffst(i,j,m),m=1,3),i=3,nthe+1),j=2,npsi)
  // 300  format(3(1x,e12.5))
  //
  // endif



  // do j=2,npsi
  // do i=1,n2th+5
  // tpst[i][j]=atan2(bxst[i][j],-bzst[i][j])
  // if(i .gt. ipi) tpst[i][j]=tpst[i][j]+2*pi
  // if(i .eq. ipi) tpst[i][j]=pi
  // bpst[i][j]=sqrt(bxst[i][j]**2+bzst[i][j]**2)
  // wstx2r[i][j]=-bzst[i][j]/bpst[i][j]
  // wstz2r[i][j]=bxst[i][j]/bpst[i][j]
  // wstx2p[i][j]=-bxst[i][j]/bpst[i][j]
  // wstz2p[i][j]=-bzst[i][j]/bpst[i][j]
  // enddo
  // enddo
  ierr = PetscFinalize( );
  return 0;
}