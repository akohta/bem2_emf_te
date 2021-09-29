/*
 * d2te_field.c
 *
 *  Created on: Sep 13, 2018
 *      Author: ohta
 */
#include "bem2_emf_te.h"

int Ez_s(double complex *Ez,double *rt,int type,DOMD *md)
{
  double complex CC[7],k,A=0.0,B=0.0;
  double F=0.0;
  int s,id,i;

  id=domain_id(rt,md);
  k=md->wd.k0*md->n[id];
  for(s=1;s<=md->bd.sb[id].Ne;s++){
    if     (type==0) coefficient_GL(CC,rt,md->bd.sb[id].eid[s],k,md);
    else if(type==1) coefficient_GK(CC,rt,md->bd.sb[id].eid[s],k,md);
    else if(type==2) coefficient_HP(CC,rt,md->bd.sb[id].eid[s],k,md);
    else if(type==3) coefficient_DE(CC,rt,md->bd.sb[id].eid[s],k,md);
    else             coefficient_NV(CC,rt,md->bd.sb[id].eid[s],k,md);

    for(i=0;i<3;i++){
      A+= md->bd.sb[id].dudn[s][i]*CC[i];
      B+=-md->bd.sb[id].u[s][i]*CC[3+i];
    }
    F+=creal(CC[6]);
  }
  if(id==0) *Ez=(A+B)/(1.0+F);
  else *Ez=(A+B)/F;

  return id;
}

int Ez_t(double complex *Ez,double *rt,int type,DOMD *md)
{
  int id;

  id=Ez_s(Ez,rt,type,md);
  if(id==0) *Ez+=infd_Ez(rt,&(md->wd));
  return id;
}

int Ez_i(double complex *Ez,double *rt,int type,DOMD *md)
{
  *Ez=infd_Ez(rt,&(md->wd));
  return domain_id(rt,md);
}

int Hr_s(double complex *Hr,double *rt,int type,DOMD *md)
{
  double complex CC[7],k,Ax=0.0,Bx=0.0,Ay=0.0,By=0.0;
  double F=0.0;
  int s,id,i;

  id=domain_id(rt,md);
  k=md->wd.k0*md->n[id];
  for(s=1;s<=md->bd.sb[id].Ne;s++){
    if     (type==0) coefficient_GL(CC,rt,md->bd.sb[id].eid[s],k,md);
    else if(type==1) coefficient_GK(CC,rt,md->bd.sb[id].eid[s],k,md);
    else if(type==2) coefficient_HP(CC,rt,md->bd.sb[id].eid[s],k,md);
    else if(type==3) coefficient_DE(CC,rt,md->bd.sb[id].eid[s],k,md);
    else             coefficient_NV(CC,rt,md->bd.sb[id].eid[s],k,md);

    for(i=0;i<3;i++){
      Ax+= md->bd.sb[id].dvdn[s][i]*CC[i];
      Bx+=-md->bd.sb[id].v[s][i]*CC[3+i];
      Ay+= md->bd.sb[id].dwdn[s][i]*CC[i];
      By+=-md->bd.sb[id].w[s][i]*CC[3+i];
    }
    F+=creal(CC[6]);
  }
  if(id==0) F+=1.0;
  Hr[0]=(Ax+Bx)/F;
  Hr[1]=(Ay+By)/F;

  return id;
}

int Hr_t(double complex *Hr,double *rt,int type,DOMD *md)
{
  double complex hi[2];
  int id;

  id=Hr_s(Hr,rt,type,md);
  if(id==0){
    infd_Hr(hi,rt,&(md->wd));
    Hr[0]+=hi[0];
    Hr[1]+=hi[1];
  }
  return id;
}

int Hr_i(double complex *Hr,double *rt,int type,DOMD *md)
{
  infd_Hr(Hr,rt,&(md->wd));
  return domain_id(rt,md);
}

int EH_s(double complex *Ez,double complex *Hr,double *rt,int type,DOMD *md)
{
  double complex CC[7],k;
  double complex A=0.0,Ax=0.0,Ay=0.0,B=0.0,Bx=0.0,By=0.0;
  double F=0.0;
  int s,id,i;

  id=domain_id(rt,md);
  k=md->wd.k0*md->n[id];
  for(s=1;s<=md->bd.sb[id].Ne;s++){
    if     (type==0) coefficient_GL(CC,rt,md->bd.sb[id].eid[s],k,md);
    else if(type==1) coefficient_GK(CC,rt,md->bd.sb[id].eid[s],k,md);
    else if(type==2) coefficient_HP(CC,rt,md->bd.sb[id].eid[s],k,md);
    else if(type==3) coefficient_DE(CC,rt,md->bd.sb[id].eid[s],k,md);
    else             coefficient_NV(CC,rt,md->bd.sb[id].eid[s],k,md);

    for(i=0;i<3;i++){
      A += md->bd.sb[id].dudn[s][i]*CC[i];      B +=-md->bd.sb[id].u[s][i]*CC[3+i];
      Ax+= md->bd.sb[id].dvdn[s][i]*CC[i];      Bx+=-md->bd.sb[id].v[s][i]*CC[3+i];
      Ay+= md->bd.sb[id].dwdn[s][i]*CC[i];      By+=-md->bd.sb[id].w[s][i]*CC[3+i];
    }
    F+=creal(CC[6]);
  }
  if(id==0) F+=1.0;
  *Ez=(A+B)/F;
  Hr[0]=(Ax+Bx)/F;
  Hr[1]=(Ay+By)/F;

  return id;
}

int EH_t(double complex *Ez,double complex *Hr,double *rt,int type,DOMD *md)
{
  double complex te,th[2];
  int id;

  id=EH_s(Ez,Hr,rt,type,md);
  if(id==0){
    infd_EH(&te,th,rt,&(md->wd));
    *Ez+=te;
    Hr[0]+=th[0];
    Hr[1]+=th[1];
  }
  return id;
}

int EH_i(double complex *Ez,double complex *Hr,double *rt,int type,DOMD *md)
{
  infd_EH(Ez,Hr,rt,&(md->wd));
  return domain_id(rt,md);
}

int EH_s_dbieq(double complex *Ez,double complex *Hr,double *rt,DOMD *md)
{
  double complex A=0.0,B=0.0,Ax=0.0,Bx=0.0,Ay=0.0,By=0.0;
  double complex CC[7],Cx[7],Cy[7],k,i_k2;
  double F=0.0,dFx=0.0,dFy=0.0;
  int s,i,id;

  id=domain_id(rt,md);
  k=md->wd.k0*md->n[id];
  i_k2=md->wd.lambda0/(I*2.0*M_PI);

  for(s=1;s<=md->bd.sb[id].Ne;s++){
    d_coef_grad_GL(CC,Cx,Cy,rt,md->bd.sb[id].eid[s],k,md);
    
    for(i=0;i<3;i++){
      A += md->bd.sb[id].dudn[s][i]*CC[i];        B +=-md->bd.sb[id].u[s][i]*CC[3+i];
      Ax+= md->bd.sb[id].dudn[s][i]*Cx[i];        Bx+=-md->bd.sb[id].u[s][i]*Cx[3+i];
      Ay+= md->bd.sb[id].dudn[s][i]*Cy[i];        By+=-md->bd.sb[id].u[s][i]*Cy[3+i];
    }
    F  +=creal(CC[6]);
    dFx+=creal(Cx[6]);
    dFy+=creal(Cy[6]);
  }
  if(id==0)  *Ez=(A+B)/(1.0+F);
  else *Ez=(A+B)/F;

  if( fabs(dFx)<CDFV && fabs(dFy)<CDFV){
    if(id==0) F+=1.0;
    Hr[0]= i_k2*(Ay+By-*Ez*dFy)/F;
    Hr[1]=-i_k2*(Ax+Bx-*Ez*dFx)/F;
    return id;
  }
  else {
    *Ez=0.0;
    Hr[0]=0.0;
    Hr[1]=0.0;
    return -1;
  }
}

int EH_t_dbieq(double complex *Ez,double complex *Hr,double *rt,DOMD *md)
{
  double complex te,th[2];
  int id;

  id=EH_s_dbieq(Ez,Hr,rt,md);
  if(id==0){
    infd_EH(&te,th,rt,&(md->wd));
    *Ez+=te;
    Hr[0]+=th[0];
    Hr[1]+=th[1];
  }
  return id;
}

int EH_i_dbieq(double complex *Ez,double complex *Hr,double *rt,DOMD *md)
{
  infd_EH(Ez,Hr,rt,&(md->wd));
  return domain_id(rt,md);
}

double complex Ez_bv(int did,double eta_t,int t,DOMD *md)
{
  double complex CC[7],k;
  double complex A=0.0,B=0.0;
  double rt[2],rs[2][3],F=0.0;
  int i,s,tid,atid,sid;

  if(-1.0 < eta_t && eta_t < 1.0){
    k=md->wd.k0*md->n[did];
    tid=md->bd.sb[did].eid[t];
    atid=abs(tid);

    for(i=0;i<3;i++){
      rs[0][i]=md->bd.x[atid][i];
      rs[1][i]=md->bd.y[atid][i];
    }
    if(tid>0) rs_eta(rt,rs, eta_t);
    else     rs_eta(rt,rs,-eta_t);

    for(s=1;s<=md->bd.sb[did].Ne;s++){
      sid=md->bd.sb[did].eid[s];
      if(s!=t) coefficient_NV(CC,rt,sid,k,md);
      else coefficient_bd_eta(CC,eta_t,sid,k,md);

      for(i=0;i<3;i++){
        A+= md->bd.sb[did].dudn[s][i]*CC[i];
        B+=-md->bd.sb[did].u[s][i]*CC[3+i];
      }
      F+=creal(CC[6]);
    }
    if(did==0) return (A+B)/(1.0+F);
    else return (A+B)/F;
  }
  else {
    printf("d2te_field.c, Ez_bv() Domain error. eta_t=%15.14e. Exit...\n",eta_t);
    exit(1);
  }
}

void EH_bv(double complex *Ez,double complex *Hr,int did,double eta_t,int t,DOMD *md)
{
  double complex CC[7],k;
  double complex A=0.0,B=0.0,Ax=0.0,Bx=0.0,Ay=0.0,By=0.0;
  double rt[2],rs[2][3],F=0.0;
  int i,s,tid,atid,sid;

  if(-1.0 < eta_t && eta_t < 1.0){
    k=md->wd.k0*md->n[did];
    tid=md->bd.sb[did].eid[t];
    atid=abs(tid);

    for(i=0;i<3;i++){
      rs[0][i]=md->bd.x[atid][i];
      rs[1][i]=md->bd.y[atid][i];
    }
    if(tid>0) rs_eta(rt,rs, eta_t);
    else     rs_eta(rt,rs,-eta_t);

    for(s=1;s<=md->bd.sb[did].Ne;s++){
      sid=md->bd.sb[did].eid[s];
      if(s!=t) coefficient_NV(CC,rt,sid,k,md);
      else coefficient_bd_eta(CC,eta_t,sid,k,md);

      for(i=0;i<3;i++){
        A += md->bd.sb[did].dudn[s][i]*CC[i];  B +=-md->bd.sb[did].u[s][i]*CC[3+i];
        Ax+= md->bd.sb[did].dvdn[s][i]*CC[i];  Bx+=-md->bd.sb[did].v[s][i]*CC[3+i];
        Ay+= md->bd.sb[did].dwdn[s][i]*CC[i];  By+=-md->bd.sb[did].w[s][i]*CC[3+i];
      }
      F+=creal(CC[6]);
    }
    if(did==0){
      *Ez=(A+B)/(1.0+F);
      Hr[0]=(Ax+Bx)/(1.0+F);
      Hr[1]=(Ay+By)/(1.0+F);
    }
    else{
      *Ez=(A+B)/F;
      Hr[0]=(Ax+Bx)/F;
      Hr[1]=(Ay+By)/F;
    }
  }
  else {
    printf("d2te_field.c, EH_bv() Domain error. eta_t=%15.14e. Exit...\n",eta_t);
    exit(1);
  }
}

double complex dEzdt_bv_node_ndmtd(int did,int t,int tn,DOMD *md)
{
  double complex FD_1(int order,double h,double x0,int did,int s,DOMD *md);
  
  double dx,dy,i_aJ;
  int tid,atid;

  tid=md->bd.sb[did].eid[t];
  atid=abs(tid);

  if(tid>0){
    dx=md->bd.dx[tid][tn];
    dy=md->bd.dy[tid][tn];
  }
  else {
    if(tn==0){
      dx=-md->bd.dx[atid][1];
      dy=-md->bd.dy[atid][1];
    }
    else if(tn==1){
      dx=-md->bd.dx[atid][0];
      dy=-md->bd.dy[atid][0];
    }
    else {
      dx=-md->bd.dx[atid][2];
      dy=-md->bd.dy[atid][2];
    }
  }
  i_aJ=1.0/sqrt(dx*dx+dy*dy);

  return i_aJ*FD_1(FDO,FDH,md->bd.et[tn],did,t,md);
}

double complex dEzdt_bv_node_dbieq(int did,int t,int tn,DOMD *md)
{
  double complex Cv[7],RH,k;
  double rt[2],drt[2],aJ,vt[2],F;
  int s,sn,tid,atid,ct,ctn;

  k=md->wd.k0*md->n[did];
  tid=md->bd.sb[did].eid[t];
  atid=abs(tid);

  if(tid>0){
    rt[0]=md->bd.x[atid][tn];
    rt[1]=md->bd.y[atid][tn];
    drt[0]=md->bd.dx[atid][tn];
    drt[1]=md->bd.dy[atid][tn];
    ct=1;
  }
  else {
    if(tn==0) ctn=1;
    else if(tn==1) ctn=0;
    else ctn=tn;
    rt[0]=md->bd.x[atid][ctn];
    rt[1]=md->bd.y[atid][ctn];
    drt[0]=md->bd.dx[atid][ctn];
    drt[1]=md->bd.dy[atid][ctn];
    ct=-1;
  }
  aJ=sqrt(drt[0]*drt[0]+drt[1]*drt[1]);
  vt[0]=(double)ct*drt[0]/aJ;
  vt[1]=(double)ct*drt[1]/aJ;

  F=0.0;
  RH=0.0;
  for(s=1;s<=md->bd.sb[did].Ne;s++){
    if(s!=t) d_coef_NV(Cv,vt,rt,md->bd.sb[did].eid[s],k,md);
    else     d_coef_bd_node_t(Cv,md->bd.sb[did].eid[t],tn,k,md);

    for(sn=0;sn<3;sn++)     RH+= md->bd.sb[did].dudn[s][sn]*Cv[sn]-md->bd.sb[did].u[s][sn]*Cv[sn+3];
    F+=creal(Cv[6]);
  }
  if(did==0)   return RH/(1.0+F);
  else return RH/F;
}

///////////////////////////////////////////////////////////////////////
double complex FD_1(int order,double h,double x0,int did,int s,DOMD *md)
{
  double complex E;
  
  if(order==2){
    E=0.5*Ez_bv(did,x0+h,s,md)-0.5*Ez_bv(did,x0-h,s,md);
    return E/h;
  }
  if(order==4){
    E= 2.0/ 3.0*(Ez_bv(did,x0+1.0*h,s,md)-Ez_bv(did,x0-1.0*h,s,md))
              -1.0/12.0*(Ez_bv(did,x0+2.0*h,s,md)-Ez_bv(did,x0-2.0*h,s,md));
    return E/h;
  }
  if(order==6){
    E= 3.0/ 4.0*(Ez_bv(did,x0+1.0*h,s,md)-Ez_bv(did,x0-1.0*h,s,md))
              -3.0/20.0*(Ez_bv(did,x0+2.0*h,s,md)-Ez_bv(did,x0-2.0*h,s,md))
            +1.0/60.0*(Ez_bv(did,x0+3.0*h,s,md)-Ez_bv(did,x0-3.0*h,s,md));
    return E/h;
  }
  else {
    printf("d2te_field.c, complex FD_1(), order number error! order=%d. exit\n",order);
    exit(1);
  }
}
