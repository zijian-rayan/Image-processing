#ifndef _IMAFLOT_H
#define _IMAFLOT_H

#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;
#include <stdlib.h>
#include <math.h>
//#include "constantes.h"
//#include "def.h"
//#include "fonctions.h"
//#include "pixel.h"
//#include "imasites.h"
#include "imadata.h"
//#include "statclass.h"


template <class T> class imaflot : public imadata<T> {
 public:
/* constructeur par defaut */
	imaflot (int nl=1, int nc=1, int nd=1) : imadata<T>(nl,nc,nd+2) {}
/* constructeur a partir d'une image de donnees */
	imaflot (imadata<T> &imadon) : imadata<T>(imadon.nblig,imadon.nbcol,imadon.nbcanaux+2) {}

/* constructeur Horn and Schunck */
	imaflot (imadata<T> &imadon1, imadata<T> &imadon2, float alpha) : 
			 imadata<T>(mini(imadon1.nlig(),imadon2.nlig()),mini(imadon1.ncol(),imadon2.ncol()),mini(imadon1.ncanaux(),imadon2.ncanaux())+2) {
		cout<<" flot optique de Horn and Schunck\n";
		int i,j,k,ii,ii0,ii2,jj,jj0,jj2,n,m,p,it,idec0,jdec0,idec1,jdec1;
		const int dim=ncanaux()-2, k_u=dim, k_v=dim+1, itmax=1000, nechelles=3;
		const float u_init=0.f, v_init=0.f, dsommax=1.e-7f;
		const bool intraechelle=1;
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				for (k=0; k<dim; k++) (*this)(i,j,k)=imadon2(i,j,k);
				(*this)(i,j,k_u)=u_init; (*this)(i,j,k_v)=v_init;
			}
		bool fini;
		float u_moy, v_moy, u_it, v_it, xx, xx0;
		double alpha2=alpha*alpha, grad2, dxx, dsom;
		float *dI_x=NULL, *dI_y=NULL, *dI_t=NULL;
		dI_x=new float[dim]; dI_y=new float[dim]; dI_t=new float[dim];
		imadata<float> T_imadon1bis[nechelles], T_imadon2bis[nechelles], T_ima_u[nechelles], T_ima_v[nechelles]; 
		T_imadon1bis[0]=imadon1; T_imadon2bis[0]=imadon2; 
		T_ima_u[0]=imadata<float>(nblig,nbcol); T_ima_u[0].mise_a_zero(); T_ima_u[0]=T_ima_u[0]+u_init; 
		T_ima_v[0]=imadata<float>(nblig,nbcol); T_ima_v[0].mise_a_zero(); T_ima_v[0]=T_ima_v[0]+v_init; 
		int nblig2=nblig, nbcol2=nbcol;
		for (m=1; m<nechelles; m++) {
			nblig2=nblig2/2; nbcol2=nbcol2/2; 
			T_imadon1bis[m]=imadata<float>(nblig2,nbcol2,dim); T_imadon2bis[m]=imadata<float>(nblig2,nbcol2,dim);
			for (i=0; i<nblig2; i++) {
				ii=i*2;
				for (j=0; j<nbcol2; j++) {
					jj=j*2;
					for (k=0; k<dim; k++) {
						T_imadon1bis[m](i,j,k)=(T_imadon1bis[m-1](ii,jj,k)+T_imadon1bis[m-1](mini(ii+1,2*nblig2),mini(jj+1,2*nbcol2),k)+
											   T_imadon1bis[m-1](mini(ii+1,2*nblig2),jj,k)+T_imadon1bis[m-1](ii,mini(jj+1,2*nbcol2),k))/4;
						T_imadon2bis[m](i,j,k)=(T_imadon2bis[m-1](ii,jj,k)+T_imadon2bis[m-1](mini(ii+1,2*nblig2),mini(jj+1,2*nbcol2),k)+
											   T_imadon2bis[m-1](mini(ii+1,2*nblig2),jj,k)+T_imadon2bis[m-1](ii,mini(jj+1,2*nbcol2),k))/4;
					}
				}
			}
			T_imadon1bis[m].statbasic(1); T_imadon2bis[m].statbasic(1); 
			T_ima_u[m]=imadata<float>(nblig2,nbcol2); T_ima_u[m].mise_a_zero(); T_ima_u[m]=T_ima_u[m]+u_init;
			T_ima_v[m]=imadata<float>(nblig2,nbcol2); T_ima_v[m].mise_a_zero(); T_ima_v[m]=T_ima_v[m]+v_init;
			cout<<" echelle "<<m<<" : \n"; T_ima_u[m].statbasic(1);  T_ima_v[m].statbasic(1);
		}
		char nomimadec[15]="./imadec_0.dat", *extension=".dat", no_ech[2]="0";
		const int deb_no=strlen(nomimadec)-strlen(no_ech)-strlen(extension);
		for (m=nechelles-1; m>=0; m--) {
			it=0; cout<<" echelle "<<m<<" : \n"; fini=0;
			imadata<float> ima_dec(nblig2,nbcol2,2); ima_dec.mise_a_zero();
			if (m<nechelles-1) {
				int m2, nbli2=nblig2/2, nbco2=nbcol2/2;
				k=2;
				for (m2=m+1; m2<nechelles; m2++) {
					for (i=0; i<nbli2; i++) {
						ii2=(i+1)*k;
						for (j=0; j<nbco2; j++) {
							u_moy=(float)mini(maxi(around(T_ima_u[m2](i,j)*k),-k),k); v_moy=(float)mini(maxi(around(T_ima_v[m2](i,j)*k),-k),k);
//							u_moy=(float)around(T_ima_u[m2](i,j)*k); v_moy=(float)around(T_ima_v[m2](i,j)*k);
							T_ima_u[m2](i,j)=u_moy/k; T_ima_v[m2](i,j)=v_moy/k;
							jj2=(j+1)*k;
							for (ii=i*k; ii<ii2; ii++)
								for (jj=j*k; jj<jj2; jj++) {
									ima_dec(ii,jj,0)=ima_dec(ii,jj,0)+u_moy;
									ima_dec(ii,jj,1)=ima_dec(ii,jj,1)+v_moy;
								}
						}
					}
					k=k*2; nbli2=nbli2/2; nbco2=nbco2/2;
				}
				*no_ech=0X30+m; strncpy_s(nomimadec+deb_no,_countof(nomimadec)-deb_no,no_ech,1); ima_dec.sauve_Ima(nomimadec);
/*				for (i=0; i<nblig2; i++) {
					for (j=0; j<nbcol2; j++) cout<<(float)((int)(ima_dec(i,j,0)*10))/10.f<<" ";cout<<"\n";} cout<<"\n";
				for (i=0; i<nblig2; i++) {
					for (j=0; j<nbcol2; j++) cout<<(float)((int)(ima_dec(i,j,1)*10))/10.f<<" ";cout<<"\n";}
				char aa; cin>>aa;*/
			}
			imadata<float> imadI_x(nblig2,nbcol2,dim), imadI_y(nblig2,nbcol2,dim), imadI_t(nblig2,nbcol2,dim); 
			imadI_x.mise_a_zero(); imadI_y.mise_a_zero(); imadI_t.mise_a_zero(); 
			for (i=0; i<nblig2; i++)
				for (j=0; j<nbcol2; j++) {
					ii0=maxi(0,i-1); ii2=mini(nblig2-1,i+1);
					jj0=maxi(0,j-1); jj2=mini(nbcol2-1,j+1);
					for (k=0; k<dim; k++) dI_x[k]=dI_y[k]=dI_t[k]=0; 
					if (j<nbcol2-1) {
						n=0;                         // calcul gradient spatial selon Horn & Schunck
						for (ii=i; ii<=ii2; ii++) {
							n+=2;
							for (k=0; k<dim; k++) {
								idec1=mini(maxi(around(ima_dec(ii,j+1,1)),ii-nblig2+1),ii); 
								idec0=mini(maxi(around(ima_dec(ii,j,1)),ii-nblig2+1),ii);
								jdec1=mini(maxi(around(ima_dec(ii,j+1,0)),j+1-nbcol2+1),(j+1)); 
								jdec0=mini(maxi(around(ima_dec(ii,j,0)),j-nbcol2+1),j);
								dI_x[k]+=(T_imadon1bis[m](ii-idec1,j+1-jdec1,k)-T_imadon1bis[m](ii-idec0,j-jdec0,k));
								dI_x[k]+=(T_imadon2bis[m](ii,j+1,k)-T_imadon2bis[m](ii,j,k));
							}
						}
						if (n>0) for (k=0; k<dim; k++) imadI_x(i,j,k)=dI_x[k]/n;
/*					for (k=0; k<dim; k++) { // version modifiee calcul du gradient sur 1 pixel de l'image 2 ?????????????????
								dI_x[k]=(T_imadon2bis[m](i,j+1,k)-T_imadon2bis[m](i,j,k));
						}
						for (k=0; k<dim; k++) { // version modifiee calcul du gradient sur 2 pixels : 1 de l'image 1 et son correspondant de l'image 2
								idec1=mini(maxi(around(ima_dec(i,j+1,1)),i-nblig2+1),i); idec0=mini(maxi(around(ima_dec(i,j,1)),i-nblig2+1),i);
								jdec1=mini(maxi(around(ima_dec(i,j+1,0)),j+1-nbcol2+1),(j+1)); jdec0=mini(maxi(around(ima_dec(i,j,0)),j-nbcol2+1),j);
								xx=(T_imadon1bis[m](i-idec1,j+1-jdec1,k)-T_imadon1bis[m](i-idec0,j-jdec0,k));
								if (fabs(xx)>fabs(dI_x[k]) dI_x[k]=xx;
						}*/
					}
					if (i<nblig2-1) {
						n=0;                          // calcul gradient spatial selon Horn & Schunck
						for (jj=j; jj<=jj2; jj++) {
							n+=2;
							for (k=0; k<dim; k++) {
								idec1=mini(maxi(around(ima_dec(i+1,jj,1)),i+1-nblig2+1),i+1); 
								idec0=mini(maxi(around(ima_dec(i,jj,1)),i-nblig2+1),i);
								jdec1=mini(maxi(around(ima_dec(i+1,jj,0)),jj-nbcol2+1),jj); 
								jdec0=mini(maxi(around(ima_dec(i,jj,0)),jj-nbcol2+1),jj);
								dI_y[k]+=(T_imadon1bis[m](i+1-idec1,jj-jdec1,k)-T_imadon1bis[m](i-idec0,jj-jdec0,k));
								dI_y[k]+=(T_imadon2bis[m](i+1,jj,k)-T_imadon2bis[m](i,jj,k));
							}
						}
						if (n>0)for (k=0; k<dim; k++) imadI_y(i,j,k)=dI_y[k]/n;
/*					for (k=0; k<dim; k++) { // version modifiee calcul du gradient sur 1 pixel de l'image 2 ?????????????????
							dI_y[k]=(T_imadon2bis[m](i+1,j,k)-T_imadon2bis[m](i,j,k));
						}
						for (k=0; k<dim; k++) { // version modifiee calcul du gradient sur 2 pixels : 1 de l'image 1 et son correspondant de l'image 2
								idec1=mini(maxi(around(ima_dec(i+1,j,1)),i+1-nblig2+1),i+1); idec0=mini(maxi(around(ima_dec(i,j,1)),i-nblig2+1),i);
								jdec1=mini(maxi(around(ima_dec(i+1,j,0)),j-nbcol2+1),j); jdec0=mini(maxi(around(ima_dec(i,j,0)),j-nbcol2+1),j);
								xx=(T_imadon1bis[m](i+1-idec1,j-jdec1,k)-T_imadon1bis[m](i-idec0,j-jdec0,k));
								if (fabs(xx)>fabs(dI_y[k]) dI_y[k]=xx;
						}*/
					}
					n=0;                             // calcul gradient spatial selon Horn & Schunck
					for (ii=i; ii<=ii2; ii++)  
						for (jj=j; jj<=jj2; jj++) {
							n+=1;
							idec0=mini(maxi(around(ima_dec(ii,jj,1)),ii-nblig2+1),ii); 
							jdec0=mini(maxi(around(ima_dec(ii,jj,0)),jj-nbcol2+1),jj);
							for (k=0; k<dim; k++)
								dI_t[k]+=(T_imadon2bis[m](ii,jj,k)-T_imadon1bis[m](ii-idec0,jj-jdec0,k));
						}
					if (n>0) for (k=0; k<dim; k++) imadI_t(i,j,k)=dI_t[k]/n;
					idec0=mini(maxi(around(ima_dec(i,j,1)),i-nblig2+1),i); // version modifiee calcul du gradient temporel
					jdec0=mini(maxi(around(ima_dec(i,j,0)),j-nbcol2+1),j);
					for (k=0; k<dim; k++) {
						xx0=(T_imadon2bis[m](i,j,k)-T_imadon1bis[m](i-idec0,j-jdec0,k));
						if (fabs(imadI_x(i,j,k))>dsommax) {
							ii=i; jj=jj2;
							idec0=mini(maxi(around(ima_dec(ii,jj,1)),ii-nblig2+1),ii); 
							jdec0=mini(maxi(around(ima_dec(ii,jj,0)),jj-nbcol2+1),jj);
							xx=(T_imadon2bis[m](ii,jj,k)-T_imadon1bis[m](ii-idec0,jj-jdec0,k));
							fabs(xx)>fabs(xx0)?dI_t[k]+=xx:dI_t[k]+=xx0;
						}
						if (fabs(imadI_y(i,j,k))>dsommax) {
							ii=ii2; jj=j;
							idec0=mini(maxi(around(ima_dec(ii,jj,1)),ii-nblig2+1),ii); 
							jdec0=mini(maxi(around(ima_dec(ii,jj,0)),jj-nbcol2+1),jj);
							xx=(T_imadon2bis[m](ii,jj,k)-T_imadon1bis[m](ii-idec0,jj-jdec0,k));
							fabs(xx)>fabs(xx0)?dI_t[k]+=xx:dI_t[k]+=xx0;
						}
/*						for (ii=i; ii<=ii2; ii++)        
							for (jj=j; jj<=jj2; jj++)
//								if (ii!=i || jj!=j) {
								if ((ii!=jj) && (ii==i || jj=j)) {
									idec0=mini(maxi(around(ima_dec(ii,jj,1)),ii-nblig2+1),ii); 
									jdec0=mini(maxi(around(ima_dec(ii,jj,0)),jj-nbcol2+1),jj);
									xx=(T_imadon2bis[m](ii,jj,k)-T_imadon1bis[m](ii-idec0,jj-jdec0,k));
									fabs(xx)>fabs(xx0)?dI_t[k]+=xx:dI_t[k]+=xx0;
								}*/
					}
				}

			while (!fini) {
				if (it%100==0) cout<<" iteration "<<++it; dsom=0;
				for (i=0; i<nblig2; i++)
					for (j=0; j<nbcol2; j++) {
						ii0=maxi(0,i-1); ii2=mini(nblig2-1,i+1);
						jj0=maxi(0,j-1); jj2=mini(nbcol2-1,j+1);
						u_moy=v_moy=0; n=0;
						for (ii=ii0; ii<=ii2; ii++)
							for (jj=jj0; jj<=jj2; jj++) {
								p=3-(abs(ii-i)+abs(jj-j)); // coef 2 devant les voisins 4-connexes et 1 devant les voisins 8-connexes
								u_moy+=(intraechelle*ima_dec(ii,jj,0)+T_ima_u[m](ii,jj))*p;
								v_moy+=(intraechelle*ima_dec(ii,jj,1)+T_ima_v[m](ii,jj))*p;
//								u_moy+=T_ima_u[m](ii,jj)*p;
//								v_moy+=T_ima_v[m](ii,jj)*p;
								n+=p;
							}
						u_moy-=3*(intraechelle*ima_dec(i,j,0)+T_ima_u[m](i,j));
						v_moy-=3*(intraechelle*ima_dec(i,j,1)+T_ima_v[m](i,j));
//						u_moy-=3*(T_ima_u[m](i,j));
//						v_moy-=3*(T_ima_v[m](i,j));
						n-=3;
						if (n>0) {u_moy/=n; v_moy/=n;}

						for (k=0; k<dim; k++) {dI_x[k]=imadI_x(i,j,k); dI_y[k]=imadI_y(i,j,k); dI_t[k]=imadI_t(i,j,k);}
						grad2=0; dxx=0;
						for (k=0; k<dim; k++) {
							grad2+=dI_x[k]*dI_x[k]+dI_y[k]*dI_y[k];
							dxx=(dI_x[k]*(u_moy-intraechelle*ima_dec(i,j,0))+dI_y[k]*(v_moy-intraechelle*ima_dec(i,j,1))+dI_t[k]);
						}
						dxx/=(alpha2*dim+grad2);
						u_it=v_it=0.;
						for (k=0; k<dim; k++) {
							u_it+=(u_moy-intraechelle*ima_dec(i,j,0))-dI_x[k]*(float)dxx;
							v_it+=(v_moy-intraechelle*ima_dec(i,j,1))-dI_y[k]*(float)dxx;
						}
						u_it/=dim; dsom+=pow(T_ima_u[m](i,j)-u_it,2); T_ima_u[m](i,j)=u_it;
						v_it/=dim; dsom+=pow(T_ima_v[m](i,j)-v_it,2); T_ima_v[m](i,j)=v_it;
					}
				if (it%100==0) cout<<" norme L2 des differences d'images a it et it-1 = "<<dsom<<" \n";
				if (it>=itmax || dsom<dsommax) fini=1;
				if (it%100==0) {
					T_ima_u[m].sauve_Ima ("D:/RECHERCHE/MoezAmmar/ImagesSivic/imaflotHS_u.tmp");
					T_ima_v[m].sauve_Ima ("D:/RECHERCHE/MoezAmmar/ImagesSivic/imaflotHS_v.tmp");
				}
			}
			cout<<" m = "<<m<<" T_ima_u[m] & T_ima_v[m] =\n";
			for (i=0; i<nblig2; i++) {
				for (j=0; j<nbcol2; j++) cout<<(float)((int)(T_ima_u[m](i,j)*10))/10.f<<" ";cout<<"\n";} cout<<"\n";
			for (i=0; i<nblig2; i++) {
				for (j=0; j<nbcol2; j++) cout<<(float)((int)(T_ima_v[m](i,j)*10))/10.f<<" ";cout<<"\n";} cout<<"\n";
			nblig2=nblig2*2; nbcol2=nbcol2*2; 
		}
		if (dI_x!=NULL) {delete[] dI_x; dI_x=NULL;}
		if (dI_y!=NULL) {delete[] dI_y; dI_y=NULL;}
		if (dI_t!=NULL) {delete[] dI_t; dI_t=NULL;}
		k=1; nblig2=nblig; nbcol2=nbcol;
		for (m=0; m<nechelles; m++) {	 
			for (i=0; i<nblig2; i++) {
				ii2=(i+1)*k;
				for (j=0; j<nbcol2; j++) {
					u_moy=T_ima_u[m](i,j)*k; v_moy=T_ima_v[m](i,j)*k;
					jj2=(j+1)*k;
					for (ii=i*k; ii<ii2; ii++)
						for (jj=j*k; jj<jj2; jj++) {
							(*this)(ii,jj,k_u)=(*this)(ii,jj,k_u)+u_moy;
							(*this)(ii,jj,k_v)=(*this)(ii,jj,k_v)+v_moy;
						}
				}
			}
			k=k*2; nblig2=nblig2/2; nbcol2=nbcol2/2;
		}
		cout<<" au final flot obtenu composantes u et v :\n";
		for (i=0; i<nblig; i++) {
			for (j=0; j<nbcol; j++) cout<<(float)((int)((*this)(i,j,k_u)*10))/10.f<<" ";cout<<"\n";} cout<<"\n";
		for (i=0; i<nblig; i++) {
			for (j=0; j<nbcol; j++) cout<<(float)((int)((*this)(i,j,k_v)*10))/10.f<<" ";cout<<"\n";}
		char aa; cin>>aa;
	}

/* constructeur Weickert and Schnorr */
	imaflot (imadata<T> &imadon1, imadata<T> &imadon2, float alpha, float lambda) : 
	    imadata<T>(mini(imadon1.nlig(),imadon2.nlig()),mini(imadon1.ncol(),imadon2.ncol()),
	               mini(imadon1.ncanaux(),imadon2.ncanaux())+2) {
		cout<<" flot optique de Weickert and Schnorr sans regularisation temporelle\n";
		int i,j,k,ii,ii0,ii2,jj,jj0,jj2,n/*, m*/,it=0;
		const int dim=ncanaux()-2, k_u=dim, k_v=dim+1, itmax=1000;
		const float dsommax=(nblig*nbcol)*1.e-5f, eps=1.e-6f, eps2=eps*eps;
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				for (k=0; k<dim; k++) (*this)(i,j,k)=imadon2(i,j,k);
				(*this)(i,j,k_u)=0; (*this)(i,j,k_v)=0;
			}
		bool fini=0;
		float u_moy, v_moy, u_it, v_it;
		double alpha2=alpha*alpha, lambda2=lambda*lambda, grad2, dpsi_i, dpsi_j, dsom;
		imadata<float> dI_x(nblig,nbcol,dim), dI_y(nblig,nbcol,dim), dI_t(nblig,nbcol,dim);
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				ii0=maxi(0,i-1); ii2=mini(nblig-1,i+1);
				jj0=maxi(0,j-1); jj2=mini(nbcol-1,j+1);
				for (k=0; k<dim; k++) dI_x(i,j,k)=0; n=0;
				if (j<nbcol-1) {
					for (ii=i; ii<=ii2; ii++) {
						n+=2;
						for (k=0; k<dim; k++) {
							dI_x(i,j,k)+=(imadon1(ii,j+1,k)-imadon1(ii,j,k));
							dI_x(i,j,k)+=(imadon2(ii,j+1,k)-imadon2(ii,j,k));
						}
					}
					if (n>0)
						for (k=0; k<dim; k++) dI_x(i,j,k)/=n;
				}
				for (k=0; k<dim; k++) dI_y(i,j,k)=0; n=0;
				if (i<nblig-1) {
					for (jj=j; jj<=jj2; jj++) {
						n+=2;
						for (k=0; k<dim; k++) {
							dI_y(i,j,k)+=(imadon1(i+1,jj,k)-imadon1(i,jj,k));
							dI_y(i,j,k)+=(imadon2(i+1,jj,k)-imadon2(i,jj,k));
						}
					}
					if (n>0)
						for (k=0; k<dim; k++) dI_y(i,j,k)/=n;
				}
				for (k=0; k<dim; k++) dI_t(i,j,k)=0; n=0;
				for (ii=i; ii<=ii2; ii++)
					for (jj=j; jj<=jj2; jj++) {
						n+=1;
						for (k=0; k<dim; k++)
							dI_t(i,j,k)+=(imadon2(ii,jj,k)-imadon1(ii,jj,k));
					}
				if (n>0)
					for (k=0; k<dim; k++) dI_t(i,j,k)/=n;
			}
//		dI_x.sauve_Ima ("D:/RECHERCHE/MoezAmmar/ImagesSivic/imadI_x.tmp",0);
//		dI_y.sauve_Ima ("D:/RECHERCHE/MoezAmmar/ImagesSivic/imadI_y.tmp",0);
//		dI_t.sauve_Ima ("D:/RECHERCHE/MoezAmmar/ImagesSivic/imadI_t.tmp",0);
		float *Sa_u=NULL, *Sa_v=NULL;
		Sa_u=new float[dim]; Sa_v=new float[dim];
		while (!fini) {
			cout<<" iteration "<<++it; dsom=0;
			for (i=0; i<nblig; i++)
				for (j=0; j<nbcol; j++) {
					ii0=maxi(0,i-1); ii2=mini(nblig-1,i+1);
					jj0=maxi(0,j-1); jj2=mini(nbcol-1,j+1);
/*					u_moy=v_moy=0; n=0;
					for (ii=ii0; ii<=ii2; ii++)
						for (jj=jj0; jj<=jj2; jj++) {
							m=3-(abs(ii-i)+abs(jj-j)); //  coef 2 devant les voisins 4-connexes et 1 devant les voisins 8-connexes
							u_moy+=(*this)(i,j,k_u)*m;
							v_moy+=(*this)(i,j,k_v)*m;
							n+=m;
						}
					u_moy-=3*(*this)(i,j,k_u);
					v_moy-=3*(*this)(i,j,k_v);
					if (n>0) {
						u_moy/=n;
						v_moy/=n;
					}*/
					u_moy=(*this)(i,j,k_u);
					v_moy=(*this)(i,j,k_v);
					u_it=v_it=0.;
					for (k=0; k<dim; k++) {
						Sa_u[k]=0; Sa_v[k]=0; 
					}
					for (k=0; k<dim; k++) {
						grad2=dI_x(i,j,k)*dI_x(i,j,k)+dI_y(i,j,k)*dI_y(i,j,k);
						dpsi_i=eps+(1-eps2)/2./pow(1+grad2/lambda2,0.5);
//						dpsi_i=1;
						n=0;
						for (ii=ii0; ii<=ii2; ii++) {
							n+=1;
							grad2=dI_x(ii,j,k)*dI_x(ii,j,k)+dI_y(ii,j,k)*dI_y(ii,j,k);
							dpsi_j=eps+(1-eps2)/2./pow(1+grad2/lambda2,0.5);
//							dpsi_j=1;
							Sa_u[k]+=0.5*(dpsi_i+dpsi_j)*((*this)(ii,j,k_u)-(*this)(i,j,k_u));
							Sa_v[k]+=0.5*(dpsi_i+dpsi_j)*((*this)(ii,j,k_v)-(*this)(i,j,k_v));
						}
						n--;
						for (jj=jj0; jj<=jj2; jj++) {
							n+=1;
							grad2=dI_x(i,jj,k)*dI_x(i,jj,k)+dI_y(i,jj,k)*dI_y(i,jj,k);
							dpsi_j=eps+(1-eps2)/2./pow(1+grad2/lambda2,0.5);
//							dpsi_j=1;
							Sa_u[k]+=0.5*(dpsi_i+dpsi_j)*((*this)(i,jj,k_u)-(*this)(i,j,k_u));
							Sa_v[k]+=0.5*(dpsi_i+dpsi_j)*((*this)(i,jj,k_v)-(*this)(i,j,k_v));
						}
						n--;
						if (n>0) {
							Sa_u[k]/=n;
							Sa_v[k]/=n;
						}
						u_it+=(u_moy+Sa_u[k]-dI_x(i,j,k)*(dI_y(i,j,k)*v_moy+dI_t(i,j,k))/alpha2)/(1+pow(dI_x(i,j,k),2)/alpha2);
						v_it+=(v_moy+Sa_v[k]-dI_y(i,j,k)*(dI_x(i,j,k)*u_moy+dI_t(i,j,k))/alpha2)/(1+pow(dI_y(i,j,k),2)/alpha2);
					}
					u_it/=dim; dsom+=pow((*this)(i,j,k_u)-u_it,2); (*this)(i,j,k_u)=u_it;
					v_it/=dim; dsom+=pow((*this)(i,j,k_v)-v_it,2); (*this)(i,j,k_v)=v_it;
				}
			cout<<" norme L2 des differences d'images a it et it-1 = "<<dsom<<" \n";
			if (it>=itmax || dsom<dsommax) fini=1;
			if (it%100==0) {
				sauve_Ima ("D:/RECHERCHE/MoezAmmar/ImagesSivic/imaflotWS_u.tmp",k_u);
				sauve_Ima ("D:/RECHERCHE/MoezAmmar/ImagesSivic/imaflotWS_v.tmp",k_v);
			}
		}
		if (Sa_u!=NULL) {delete[] Sa_u; Sa_u=NULL;}
		if (Sa_v!=NULL) {delete[] Sa_v; Sa_v=NULL;}
	}

/* constructeur Weickert and Schnorr */
	imaflot (imaflot<T> &imaflo1, imadata<T> &imadon2, float alpha, float lambda) : 
		imadata<T>(mini(imaflo1.nlig(),imadon2.nlig()),mini(imaflo1.ncol(),imadon2.ncol()),
				   mini(imaflo1.ncanaux(),imadon2.ncanaux()+2)) {
		cout<<" flot optique de Weickert and Schnorr avec regularisation temporelle\n";
		int i,j,k,ii,ii0,ii2,jj,jj0,jj2,n/*, m*/,it=0;
		const int dim=ncanaux()-2, k_u=dim, k_v=dim+1, itmax=10000;
		const float dsommax=(nblig*nbcol)*1.e-5f, eps=1.e-6f, eps2=eps*eps;
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				for (k=0; k<dim; k++) (*this)(i,j,k)=imadon2(i,j,k);
				for (k=k_u; k<=k_v; k+=k_v-k_u) (*this)(i,j,k)=imaflo1(i,j,k);
				}
		bool fini=0;
		float u_moy, v_moy, u_it, v_it;
		double alpha2=alpha*alpha, lambda2=lambda*lambda, grad2, dpsi_i, dpsi_j, dsom;
		imadata<float> dI_x(nblig,nbcol,dim), dI_y(nblig,nbcol,dim), dI_t(nblig,nbcol,dim);
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				ii0=maxi(0,i-1); ii2=mini(nblig-1,i+1);
				jj0=maxi(0,j-1); jj2=mini(nbcol-1,j+1);
				for (k=0; k<dim; k++) dI_x(i,j,k)=0; n=0;
				if (j<nbcol-1) {
					for (ii=i; ii<=ii2; ii++) {
						n+=2;
						for (k=0; k<dim; k++) {
							dI_x(i,j,k)+=(imaflo1(ii,j+1,k)-imaflo1(ii,j,k));
							dI_x(i,j,k)+=(imadon2(ii,j+1,k)-imadon2(ii,j,k));
						}
					}
					if (n>0)
						for (k=0; k<dim; k++) dI_x(i,j,k)/=n;
				}
				for (k=0; k<dim; k++) dI_y(i,j,k)=0; n=0;
				if (i<nblig-1) {
					for (jj=j; jj<=jj2; jj++) {
						n+=2;
						for (k=0; k<dim; k++) {
							dI_y(i,j,k)+=(imaflo1(i+1,jj,k)-imaflo1(i,jj,k));
							dI_y(i,j,k)+=(imadon2(i+1,jj,k)-imadon2(i,jj,k));
						}
					}
					if (n>0)
						for (k=0; k<dim; k++) dI_y(i,j,k)/=n;
				}
				for (k=0; k<dim; k++) dI_t(i,j,k)=0; n=0;
				for (ii=i; ii<=ii2; ii++)
					for (jj=j; jj<=jj2; jj++) {
						n+=1;
						for (k=0; k<dim; k++)
							dI_t(i,j,k)+=(imadon2(ii,jj,k)-imaflo1(ii,jj,k));
					}
				if (n>0)
					for (k=0; k<dim; k++) dI_t(i,j,k)/=n;
			}
//		dI_x.sauve_Ima ("D:/RECHERCHE/MoezAmmar/ImagesSivic/imadI_x.tmp",0);
//		dI_y.sauve_Ima ("D:/RECHERCHE/MoezAmmar/ImagesSivic/imadI_y.tmp",0);
//		dI_t.sauve_Ima ("D:/RECHERCHE/MoezAmmar/ImagesSivic/imadI_t.tmp",0);
		float *Sa_u=NULL, *Sa_v=NULL;
		Sa_u=new float[dim]; Sa_v=new float[dim];
		while (!fini) {
			cout<<" iteration "<<++it; dsom=0;
			for (i=0; i<nblig; i++)
				for (j=0; j<nbcol; j++) {
					ii0=maxi(0,i-1); ii2=mini(nblig-1,i+1);
					jj0=maxi(0,j-1); jj2=mini(nbcol-1,j+1);
/*					u_moy=v_moy=0; n=0;
					for (ii=ii0; ii<=ii2; ii++)
						for (jj=jj0; jj<=jj2; jj++) {
							m=3-(abs(ii-i)+abs(jj-j)); //  coef 2 devant les voisins 4-connexes et 1 devant les voisins 8-connexes
							u_moy+=(*this)(i,j,k_u)*m;
							v_moy+=(*this)(i,j,k_v)*m;
							n+=m;
						}
					u_moy-=3*(*this)(i,j,k_u);
					v_moy-=3*(*this)(i,j,k_v);
					if (n>0) {
						u_moy/=n;
						v_moy/=n;
					}*/
					u_moy=(*this)(i,j,k_u);
					v_moy=(*this)(i,j,k_v);
					u_it=v_it=0.;
					for (k=0; k<dim; k++) {
						Sa_u[k]=0; Sa_v[k]=0; 
					}
					for (k=0; k<dim; k++) {
						grad2=dI_x(i,j,k)*dI_x(i,j,k)+dI_y(i,j,k)*dI_y(i,j,k);
						dpsi_i=eps+(1-eps2)/2./pow(1+grad2/lambda2,0.5);
//						dpsi_i=1;
						n=0;
						for (ii=ii0; ii<=ii2; ii++) {
							n++;
							grad2=dI_x(ii,j,k)*dI_x(ii,j,k)+dI_y(ii,j,k)*dI_y(ii,j,k);
							dpsi_j=eps+(1-eps2)/2./pow(1+grad2/lambda2,0.5);
//							dpsi_j=1;
							Sa_u[k]+=0.5*(dpsi_i+dpsi_j)*((*this)(ii,j,k_u)-(*this)(i,j,k_u));
							Sa_v[k]+=0.5*(dpsi_i+dpsi_j)*((*this)(ii,j,k_v)-(*this)(i,j,k_v));
						}
						n--;
						for (jj=jj0; jj<=jj2; jj++) {
							n++;
							grad2=dI_x(i,jj,k)*dI_x(i,jj,k)+dI_y(i,jj,k)*dI_y(i,jj,k);
							dpsi_j=eps+(1-eps2)/2./pow(1+grad2/lambda2,0.5);
//							dpsi_j=1;
							Sa_u[k]+=0.5*(dpsi_i+dpsi_j)*((*this)(i,jj,k_u)-(*this)(i,j,k_u));
							Sa_v[k]+=0.5*(dpsi_i+dpsi_j)*((*this)(i,jj,k_v)-(*this)(i,j,k_v));
						}
						n--;
						dpsi_j=dpsi_i=1;
						Sa_u[k]+=0.5*(dpsi_i+dpsi_j)*(imaflo1(i,j,k_u)-(*this)(i,j,k_u));
						Sa_v[k]+=0.5*(dpsi_i+dpsi_j)*(imaflo1(i,j,k_v)-(*this)(i,j,k_v));
						n++;
						if (n>0) {
							Sa_u[k]/=n;
							Sa_v[k]/=n;
						}
						u_it+=(u_moy+Sa_u[k]-dI_x(i,j,k)*(dI_y(i,j,k)*v_moy+dI_t(i,j,k))/alpha2)/(1+pow(dI_x(i,j,k),2)/alpha2);
						v_it+=(v_moy+Sa_v[k]-dI_y(i,j,k)*(dI_x(i,j,k)*u_moy+dI_t(i,j,k))/alpha2)/(1+pow(dI_y(i,j,k),2)/alpha2);
					}
					u_it/=dim; dsom+=pow((*this)(i,j,k_u)-u_it,2); (*this)(i,j,k_u)=u_it;
					v_it/=dim; dsom+=pow((*this)(i,j,k_v)-v_it,2); (*this)(i,j,k_v)=v_it;
				}
			cout<<" norme L2 des differences d'images a it et it-1 = "<<dsom<<" \n";
			if (it>=itmax || dsom<dsommax) fini=1;
			if (it%100==0) {
				sauve_Ima ("D:/RECHERCHE/MoezAmmar/ImagesSivic/imaflotWS_u.tmp",k_u);
				sauve_Ima ("D:/RECHERCHE/MoezAmmar/ImagesSivic/imaflotWS_v.tmp",k_v);
			}
		}
		if (Sa_u!=NULL) {delete[] Sa_u; Sa_u=NULL;}
		if (Sa_v!=NULL) {delete[] Sa_v; Sa_v=NULL;}
		cout<<" fin du calcul flot optique\n";
	}

};

/* constructeur regularisation  */
/*	imaflot (imadata<T> &imadon1, imadata<T> &imadon2, float alpha, float lambda) : 
	    imadata<T>(mini(imadon1.nlig(),imadon2.nlig()),mini(imadon1.ncol(),imadon2.ncol()),
	               mini(imadon1.ncanaux(),imadon2.ncanaux())+2) {
		cout<<" flot optique de Weickert and Schnorr sans regularisation temporelle\n";
		int i,j,k,ii,ii0,ii2,jj,jj0,jj2,n/*, m*//*,it=0;
		const int dim=ncanaux()-2, k_u=dim, k_v=dim+1, itmax=1000;
		const float dsommax=(nblig*nbcol)*1.e-5f, eps=1.e-6f, eps2=eps*eps;
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				for (k=0; k<dim; k++) (*this)(i,j,k)=imadon2(i,j,k);
				(*this)(i,j,k_u)=0; (*this)(i,j,k_v)=0;
			}
		bool fini=0;
		float u_moy, v_moy, u_it, v_it;
		double alpha2=alpha*alpha, lambda2=lambda*lambda, grad2, dpsi_i, dpsi_j, dsom;
		imadata<float> dI_x(nblig,nbcol,dim), dI_y(nblig,nbcol,dim), dI_t(nblig,nbcol,dim);
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				ii0=maxi(0,i-1); ii2=mini(nblig-1,i+1);
				jj0=maxi(0,j-1); jj2=mini(nbcol-1,j+1);
				for (k=0; k<dim; k++) dI_x(i,j,k)=0; n=0;
				if (j<nbcol-1) {
					for (ii=i; ii<=ii2; ii++) {
						n+=2;
						for (k=0; k<dim; k++) {
							dI_x(i,j,k)+=(imadon1(ii,j+1,k)-imadon1(ii,j,k));
							dI_x(i,j,k)+=(imadon2(ii,j+1,k)-imadon2(ii,j,k));
						}
					}
					if (n>0)
						for (k=0; k<dim; k++) dI_x(i,j,k)/=n;
				}
				for (k=0; k<dim; k++) dI_y(i,j,k)=0; n=0;
				if (i<nblig-1) {
					for (jj=j; jj<=jj2; jj++) {
						n+=2;
						for (k=0; k<dim; k++) {
							dI_y(i,j,k)+=(imadon1(i+1,jj,k)-imadon1(i,jj,k));
							dI_y(i,j,k)+=(imadon2(i+1,jj,k)-imadon2(i,jj,k));
						}
					}
					if (n>0)
						for (k=0; k<dim; k++) dI_y(i,j,k)/=n;
				}
				for (k=0; k<dim; k++) dI_t(i,j,k)=0; n=0;
				for (ii=i; ii<=ii2; ii++)
					for (jj=j; jj<=jj2; jj++) {
						n+=1;
						for (k=0; k<dim; k++)
							dI_t(i,j,k)+=(imadon2(ii,jj,k)-imadon1(ii,jj,k));
					}
				if (n>0)
					for (k=0; k<dim; k++) dI_t(i,j,k)/=n;
			}
//		dI_x.sauve_Ima ("D:/RECHERCHE/MoezAmmar/ImagesSivic/imadI_x.tmp",0);
//		dI_y.sauve_Ima ("D:/RECHERCHE/MoezAmmar/ImagesSivic/imadI_y.tmp",0);
//		dI_t.sauve_Ima ("D:/RECHERCHE/MoezAmmar/ImagesSivic/imadI_t.tmp",0);
		float *Sa_u=NULL, *Sa_v=NULL;
		Sa_u=new float[dim]; Sa_v=new float[dim];
		while (!fini) {
			cout<<" iteration "<<++it; dsom=0;
			for (i=0; i<nblig; i++)
				for (j=0; j<nbcol; j++) {
					ii0=maxi(0,i-1); ii2=mini(nblig-1,i+1);
					jj0=maxi(0,j-1); jj2=mini(nbcol-1,j+1);
/*					u_moy=v_moy=0; n=0;
					for (ii=ii0; ii<=ii2; ii++)
						for (jj=jj0; jj<=jj2; jj++) {
							m=3-(abs(ii-i)+abs(jj-j)); //  coef 2 devant les voisins 4-connexes et 1 devant les voisins 8-connexes
							u_moy+=(*this)(i,j,k_u)*m;
							v_moy+=(*this)(i,j,k_v)*m;
							n+=m;
						}
					u_moy-=3*(*this)(i,j,k_u);
					v_moy-=3*(*this)(i,j,k_v);
					if (n>0) {
						u_moy/=n;
						v_moy/=n;
					}*//*
					u_moy=(*this)(i,j,k_u);
					v_moy=(*this)(i,j,k_v);
					u_it=v_it=0.;
					for (k=0; k<dim; k++) {
						Sa_u[k]=0; Sa_v[k]=0; 
					}
					for (k=0; k<dim; k++) {
						grad2=dI_x(i,j,k)*dI_x(i,j,k)+dI_y(i,j,k)*dI_y(i,j,k);
						dpsi_i=eps+(1-eps2)/2./pow(1+grad2/lambda2,0.5);
//						dpsi_i=1;
						n=0;
						for (ii=ii0; ii<=ii2; ii++) {
							n+=1;
							grad2=dI_x(ii,j,k)*dI_x(ii,j,k)+dI_y(ii,j,k)*dI_y(ii,j,k);
							dpsi_j=eps+(1-eps2)/2./pow(1+grad2/lambda2,0.5);
//							dpsi_j=1;
							Sa_u[k]+=0.5*(dpsi_i+dpsi_j)*((*this)(ii,j,k_u)-(*this)(i,j,k_u));
							Sa_v[k]+=0.5*(dpsi_i+dpsi_j)*((*this)(ii,j,k_v)-(*this)(i,j,k_v));
						}
						n--;
						for (jj=jj0; jj<=jj2; jj++) {
							n+=1;
							grad2=dI_x(i,jj,k)*dI_x(i,jj,k)+dI_y(i,jj,k)*dI_y(i,jj,k);
							dpsi_j=eps+(1-eps2)/2./pow(1+grad2/lambda2,0.5);
//							dpsi_j=1;
							Sa_u[k]+=0.5*(dpsi_i+dpsi_j)*((*this)(i,jj,k_u)-(*this)(i,j,k_u));
							Sa_v[k]+=0.5*(dpsi_i+dpsi_j)*((*this)(i,jj,k_v)-(*this)(i,j,k_v));
						}
						n--;
						if (n>0) {
							Sa_u[k]/=n;
							Sa_v[k]/=n;
						}
						u_it+=(u_moy+Sa_u[k]-dI_x(i,j,k)*(dI_y(i,j,k)*v_moy+dI_t(i,j,k))/alpha2)/(1+pow(dI_x(i,j,k),2)/alpha2);
						v_it+=(v_moy+Sa_v[k]-dI_y(i,j,k)*(dI_x(i,j,k)*u_moy+dI_t(i,j,k))/alpha2)/(1+pow(dI_y(i,j,k),2)/alpha2);
					}
					u_it/=dim; dsom+=pow((*this)(i,j,k_u)-u_it,2); (*this)(i,j,k_u)=u_it;
					v_it/=dim; dsom+=pow((*this)(i,j,k_v)-v_it,2); (*this)(i,j,k_v)=v_it;
				}
			cout<<" norme L2 des differences d'images a it et it-1 = "<<dsom<<" \n";
			if (it>=itmax || dsom<dsommax) fini=1;
			if (it%100==0) {
				sauve_Ima ("D:/RECHERCHE/MoezAmmar/ImagesSivic/imaflotWS_u.tmp",k_u);
				sauve_Ima ("D:/RECHERCHE/MoezAmmar/ImagesSivic/imaflotWS_v.tmp",k_v);
			}
		}
		if (Sa_u!=NULL) delete[] Sa_u;
		if (Sa_v!=NULL) delete[] Sa_v;
	}*/

#endif