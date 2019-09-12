#include "imacontours.h"

imacontours::imacontours (imadata<float> &ima, char *masq, unsigned short int lgmax, bool isauv)
						  : imasites(ima.nlig(),ima.ncol()) { //cout<<" entree dans imacontours 1\n";
	Tima=NULL;
	nblayer=ima.ncanaux();
	valnul=0;
	Tima=new BYTE[nbpix*nblayer];
	int i,j,k;
	long int l;
	for (l=0; l<nbpix; l++) {
		for (j=0; j<nblayer; j++) Tima[l*nblayer+j]=valnul;
		Tsites[l]=&(Tima[l*nblayer]);
	}
	imadata<float> ima_grd(nblig,nbcol,nblayer), ima_dir(nblig,nbcol,nblayer);
	ima_grd=ima.gradient(ima_dir,masq);
	if (isauv) {
//		string nomfich="imagrd_"; nomfich+=masq; nomfich+=".dat"; cout<<" image de la norme du gradient dans "<<nomfich<<"\n";
//		ima_grd.imaunsignedchar(1).sauve_ImaBSQ(nomfich);
		char nomfich[80]; strcpy_s(nomfich,"imagrd_"); strcat_s(nomfich,masq); strcat_s(nomfich,".pgm");
		cout<<" image de la norme du gradient dans "<<nomfich<<"\n";
		ima_grd.sauve_ImaPGM(nomfich);
	}
	imadata<BYTE> imaRes(nblig,nbcol,nblayer);
	for (k=0; k<nblayer; k++) {
		float sh=(float)ima_grd.moyI(k)+(float)sqrt(ima_grd.varI(k));
		float sb=(float)ima_grd.moyI(k);
//		ima_grd.seuil_ima (sb,imaRes,k);
		ima_grd.seuil_hysteresis (sh,sb,imaRes,k);
	}
	for (k=0; k<nblayer; k++)
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) (*this)(i,j,k)=imaRes(i,j,k);
//	imagrd_maxloc(ima_grd,ima_dir);
//	if (lgmax>0) prolongecontours(ima_grd,lgmax);
}

imacontours::imacontours (imadata<float> &ima, char *masqL, char *masqG, bool isauv)
						 : imasites(ima.nlig(),ima.ncol()) { //cout<<" entree dans imacontours 2\n";
	Tima=NULL;
	nblayer=ima.ncanaux();
	valnul=0;
	Tima=new BYTE[nbpix*nblayer];
	int i,j,k,i2,j2/*,i0,j0,ii,jj,dir*/;
	long int l;
/*	char nomfich[80];*/
	for (l=0; l<nbpix; l++) {
		for (j=0; j<nblayer; j++) Tima[l*nblayer+j]=valnul;
		Tsites[l]=&(Tima[l*nblayer]);
	}
/*	double sum, coef;*/
	imadata<float> ima_gdVH(nblig,nbcol,2), ima_grd(nblig,nbcol,nblayer), ima_dir(nblig,nbcol,nblayer),
				   ima_Lpc(nblig,nbcol,nblayer);
	imadata<BYTE> imaRes(nblig,nbcol,nblayer);
// calcul de l'image de la norme du gradient
	ima_grd=ima.gradient(ima_dir,masqG);
	for (k=0; k<nblayer; k++) {
		float sh=(float)ima_grd.moyI(k)+(float)sqrt(ima_grd.varI(k));
		float sb=(float)ima_grd.moyI(k);
//		ima_grd.seuil_ima (sb,imaRes,k);
		ima_grd.seuil_hysteresis (sh,sb,imaRes,k);
	}
// calcul de l'image du Laplacien
	ima_Lpc=ima.laplacien(masqL);
/*	double *noyauL=NULL;
	int taille_masqL=3;
	const int dimL=taille_masqL, dimL2=dimL/2;
	noyauL=new double[dimL*dimL];
	for (i=0; i<dimL*dimL; i++) noyauL[i]=-1;
	noyauL[dimL2*dimL+dimL2]=8;
	if (stricmp(masqL,"4connex")==0) {
		noyauL[0]=noyauL[dimL-1]=noyauL[(dimL-1)*dimL]=noyauL[dimL*dimL-1]=0;
		noyauL[dimL2*dimL+dimL2]=4;
	}
	for (k=0; k<nblayer; k++) {
		for (i=0; i<nblig; i++) {
			i0=maxi(0,i-dimL2); i2=mini(i+dimL2,nblig-1);
			ii=mini(i-i0,i2-i); i0=i-ii; i2=i+ii;
			for (j=0; j<nbcol; j++) {
				j0=maxi(0,j-dimL2); j2=mini(j+dimL2,nbcol-1);
				jj=mini(j-j0,j2-j); j0=j-jj; j2=j+jj;
				sum=0;
				for (ii=i0; ii<=i2; ii++)
					for (jj=j0; jj<=j2; jj++) {
						coef=noyauL[(ii-i+dimL2)*dimL+jj-j+dimL2];
						sum+=ima(ii,jj,k)*coef;
					}
				ima_Lpc(i,j,k)=sum;
			}
		}
	}
	if (isauv) {
//		string nomfich="imaLpc_"; nomfich+=masqL; nomfich+=".dat";
//		ima_Lpc.imaunsignedchar(1).sauve_ImaBSQ(nomfich);
		strcpy(nomfich,"imaLpc_"); strcat(nomfich,masqL); strcat(nomfich,".pgm");
		ima_Lpc.sauve_ImaPGM(nomfich);
	}
	if (noyauL!=NULL) delete[] noyauL;
	ima_Lpc.statbasic(1);*/
	double epsilon=0.01;
	for (k=0; k<nblayer; k++) {
		for (i=0; i<nblig; i++) {
			i2=mini(i+1,nblig-1);
			for (j=0; j<nbcol; j++) {
				if (imaRes(i,j,k)==1) {
					j2=mini(j+1,nbcol-1);
					if (ima_Lpc(i,j,k)>=epsilon) {
						if (ima_Lpc(i2,j,k)<epsilon || ima_Lpc(i,j2,k)<epsilon ||
							ima_Lpc(i2,j2,k)<epsilon) (*this)(i,j,k)=1;
					} else {
						if (ima_Lpc(i2,j,k)>=epsilon || ima_Lpc(i,j2,k)>=epsilon ||
							ima_Lpc(i2,j2,k)>=epsilon) (*this)(i,j,k)=1;
					}
				}
			}
		}
	}
}

imacontours::imacontours (imadata<float> &ima, float sigma2, float coef, bool isauv) 
						 : imasites(ima.nlig(),ima.ncol()) { //cout<<" entree dans imacontours 3\n";
	Tima=NULL;
	nblayer=ima.ncanaux();
	valnul=0;
	Tima=new BYTE[nbpix*nblayer];
	int i,j;
	long int l;
	for (l=0; l<nbpix; l++) {
		for (j=0; j<nblayer; j++) Tima[l*nblayer+j]=valnul;
		Tsites[l]=&(Tima[l*nblayer]);
	}
	int dim2=(int)pow(2.*sigma2*log(60.),0.5), dim=2*dim2+1; cout<<"dim2 "<<dim2<<"\n";
	int *g1=new int[dim2+1];
	int *g2=new int[dim2+1];
	int k,i0,i2,j0,j2,ii,jj;
	double fact, efact, sum;
	for (i=0; i<dim2+1; i++) {
		fact=-pow((double)i,2)/sigma2; efact=pow(2./coef,-0.5)*exp(fact/2);
		g1[i]=around((1+fact)*efact); g2[i]=around(efact);
		}
	cout<<"g1 : "; for (i=0; i<dim2+1; i++) cout<<" "<<g1[i]; cout<<"\n";
	cout<<"g2 : "; for (i=0; i<dim2+1; i++) cout<<" "<<g2[i]; cout<<"\n";
	imadata<float> imaLpg(ima), imag1g2(ima), imag2g1(ima);
	for (k=0; k<nblayer; k++) {
		for (i=0; i<nblig; i++) {
			i0=maxi(0,i-dim2);
			i2=mini(i+dim2,nblig-1);
			for (j=0; j<nbcol; j++) {
				sum=0.;
				for (ii=i0; ii<=i2; ii++) sum+=imag1g2(ii,j,k)*g1[abs(ii-i)];
				imaLpg(i,j,k)=(float)sum;
			}
		}
		for (i=0; i<nblig; i++) 
			for (j=0; j<nbcol; j++) imag1g2(i,j,k)=imaLpg(i,j,k);
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				j0=maxi(0,j-dim2);
				j2=mini(j+dim2,nbcol-1);
				sum=0.;
				for (jj=j0; jj<=j2; jj++) sum+=imag1g2(i,jj,k)*g2[abs(jj-j)];
				imaLpg(i,j,k)=(float)sum;
			}
		for (i=0; i<nblig; i++) 
			for (j=0; j<nbcol; j++) imag1g2(i,j,k)=imaLpg(i,j,k);
		for (i=0; i<nblig; i++) {
			i0=maxi(0,i-dim2);
			i2=mini(i+dim2,nblig-1);
			for (j=0; j<nbcol; j++) {
				sum=0.;
				for (ii=i0; ii<=i2; ii++) sum+=imag2g1(ii,j,k)*g2[abs(ii-i)];
				imaLpg(i,j,k)=(float)sum;
			}
		}
		for (i=0; i<nblig; i++) 
			for (j=0; j<nbcol; j++) imag2g1(i,j,k)=imaLpg(i,j,k);
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				j0=maxi(0,j-dim2);
				j2=mini(j+dim2,nbcol-1);
				sum=0.;
				for (jj=j0; jj<=j2; jj++) sum+=imag2g1(i,jj,k)*g1[abs(jj-j)];
				imaLpg(i,j,k)=(float)sum;
			}
		for (i=0; i<nblig; i++) 
			for (j=0; j<nbcol; j++) imag2g1(i,j,k)=imaLpg(i,j,k);
	}
	imaLpg=imag1g2+imag2g1;
	if (isauv) imaLpg.sauve_ImaBSQ("ima_LaplacienHuertas.dat");
	imaLpg.statbasic(1);
	for (k=0; k<nblayer; k++) {
		coef=(float)pow(imaLpg.varI(k),0.5);
		for (i=0; i<nblig; i++) {
			i2=mini(i+1,nblig-1);
			for (j=0; j<nbcol; j++) {
				(*this)(i,j,k)=0;
				j2=mini(j+1,nbcol-1);
				if (imaLpg(i,j,k)>=coef) {
					if (imaLpg(i2,j,k)<coef || imaLpg(i,j2,k)<coef ||
						imaLpg(i2,j2,k)<coef) (*this)(i,j,k)=1;
				} else {
					if (imaLpg(i2,j,k)>=coef || imaLpg(i,j2,k)>=coef ||
						imaLpg(i2,j2,k)>=coef) (*this)(i,j,k)=1;
				}
			}
		}
	}
	if (g1!=NULL) delete[] g1;
	if (g2!=NULL) delete[] g2;
}


/*imacontours::imacontours (imadata<float> &ima, float gamma, bool isauv, char *detecteur) 
						 : imasites(ima.nlig(),ima.ncol()) {cout<<" entree dans detecteur de Harris-Laplace de la classe 'imacontours'\n";
	Tima=NULL; isauv=1;
	nblayer=ima.ncanaux();
	valnul=0;
	Tima=new BYTE[nbpix*nblayer]; //cout<<nblig<<" "<<nbcol<<" "<<nbpix<<" "<<nblayer<<"\n";
	long int ll;
	int i,j,k,l,j0,j2,jj,i0,i2,ii;
	for (ll=0; ll<nbpix; ll++) {
		for (j=0; j<nblayer; j++) Tima[ll*nblayer+j]=valnul;
		Tsites[ll]=&(Tima[ll*nblayer]);
	}
	const int nsigma=3, nsigmapas=5; 
	const float kapa=0.1f, sigmamin=1.f/nsigma, sigmamax=6.f/nsigma, sigmapas=(sigmamax-sigmamin)/(nsigmapas-1), nsig_seuil=0.5f;
	double *G,sigma,sigma2,coefG,xx,s,yy,xxLx2,xxLy2,xxLxLy;
	float seuil;
	int ncoef;
	bool OK,OK2;
	imadata<float> imatemp(ima),imaL(nblig,nbcol,nblayer),imaLx(nblig,nbcol,nblayer),imaLy(nblig,nbcol,nblayer);
  imadata<double> imaLx2(nblig,nbcol,nblayer),imaLy2(nblig,nbcol,nblayer),imaLxLy(nblig,nbcol,nblayer),
				   imaLx2temp(nblig,nbcol,nblayer),imaLy2temp(nblig,nbcol,nblayer),imaLxLytemp(nblig,nbcol,nblayer),
           imaM(nblig,nbcol,nblayer),imaLapl(nblig,nbcol,nblayer*nsigmapas);
	imadata<bool> ima_scale(nblig,nbcol,nblayer*nsigmapas);
	for (sigma=sigmamin; sigma<=sigmamax; sigma+=sigmapas) {//cout<<sigmamin<<" "<<sigmamax<<" "<<sigmapas<<" "<<sigma<<"\n";
		ncoef=around(2.f*(float)nsigma*sigma+1.f); // calcul du noyau gaussien à l'échelle considérée
		G=new double[ncoef];
		sigma2=sigma*sigma;
		coefG=1./2./PI/sigma2;
		for (k=0; k<ncoef; k++) {xx=k-(ncoef-1)/2.; G[k]=coefG*exp(-xx*xx/2/sigma2);}
		for (i=0; i<nblig; i++)                    // calcul de l'image filtrée par le noyau gaussien
			for (j=0; j<nbcol; j++) { 
				j0=maxi(j-(ncoef-1)/2,0); j2=mini(j+(ncoef-1)/2,nbcol-1);
				s=0;
				for (jj=j0; jj<=j2; jj++) s+=G[jj-j+(ncoef-1)/2];
				if (s>0) {
					for (l=0; l<nblayer; l++) {
						xx=0.;
						for (jj=j0; jj<=j2; jj++) xx+=G[jj-j+(ncoef-1)/2]*imatemp(i,jj,l);
						imaL(i,j,l)=(float)(xx/s);
					}
				}
			}
		imatemp=imaL;
		for (j=0; j<nbcol; j++)
			for (i=0; i<nblig; i++) {
				i0=maxi(i-(ncoef-1)/2,0); i2=mini(i+(ncoef-1)/2,nblig-1); 
				s=0;
				for (ii=i0; ii<=i2; ii++) s+=G[ii-i+(ncoef-1)/2];
				if (s>0) {
					for (l=0; l<nblayer; l++) {
						xx=0.;
						for (ii=i0; ii<=i2; ii++) xx+=G[ii-i+(ncoef-1)/2]*imatemp(ii,j,l);
						imaL(i,j,l)=(float)(xx/s);
					}
				}
			}
		for (l=0; l<nblayer; l++) {                // calcul des dérivées partielles
			for (i=0; i<nblig; i++) {
				imaLx(i,0,l)=0.;
				for (j=1; j<nbcol; j++)
					imaLx(i,j,l)=imaL(i,j,l)-imaL(i,j-1,l);
			}
			for (j=0; j<nbcol; j++) {
				imaLy(0,j,l)=0.;
				for (i=1; i<nblig; i++)
					imaLy(i,j,l)=imaL(i,j,l)-imaL(i-1,j,l);
			}
		}
		if (G!=NULL) delete[] G;
		sigma2=sigma*sigma/(gamma*gamma);          // calcul du nouveau noyau gaussien (relation entre l'échelle et le paramètre d'intégration : s=t*gamma^2)
		coefG=1./2./PI/sigma2;
		ncoef=around(2.f*(float)nsigma*pow(sigma2,0.5)+1.f); 
		G=new double[ncoef];
		for (k=0; k<ncoef; k++) {
			xx=k-(ncoef-1)/2.;
			G[k]=coefG*exp(-xx*xx/2/sigma2);
		}
		for (i=0; i<nblig; i++)                    // calcul des carrés des dérivées
			for (j=0; j<nbcol; j++) 
				for (l=0; l<nblayer; l++) {
					xx=imaLx(i,j,l); imaLx2(i,j,l)=imaLx2temp(i,j,l)=xx*xx;
					yy=imaLy(i,j,l); imaLy2(i,j,l)=imaLy2temp(i,j,l)=yy*yy;
					imaLxLy(i,j,l)=imaLxLytemp(i,j,l)=xx*yy;
				}
		for (i=0; i<nblig; i++)                    // filtrage gaussien des dérivées au carré
			for (j=0; j<nbcol; j++) { 
				j0=maxi(j-(ncoef-1)/2,0); j2=mini(j+(ncoef-1)/2,nbcol-1);
				s=0;
				for (jj=j0; jj<=j2; jj++) s+=G[jj-j+(ncoef-1)/2];
				if (s>0) {
					for (l=0; l<nblayer; l++) {
						xxLx2=xxLy2=xxLxLy=0.;
						for (jj=j0; jj<=j2; jj++) {
							xxLx2+=G[jj-j+(ncoef-1)/2]*imaLx2(i,jj,l);
							xxLy2+=G[jj-j+(ncoef-1)/2]*imaLy2(i,jj,l);
							xxLxLy+=G[jj-j+(ncoef-1)/2]*imaLxLy(i,jj,l);
						}
						imaLx2temp(i,j,l)=xxLx2/s; imaLy2temp(i,j,l)=xxLy2/s; imaLxLytemp(i,j,l)=xxLxLy/s;
					}
				}
			}
		for (j=0; j<nbcol; j++)
			for (i=0; i<nblig; i++) {
				i0=maxi(i-(ncoef-1)/2,0); i2=mini(i+(ncoef-1)/2,nblig-1); 
				s=0;
				for (ii=i0; ii<=i2; ii++) s+=G[ii-i+(ncoef-1)/2];
				if (s>0) {
					for (l=0; l<nblayer; l++) {
						xxLx2=xxLy2=xxLxLy=0.;
						for (ii=i0; ii<=i2; ii++) {
							xxLx2+=G[ii-i+(ncoef-1)/2]*imaLx2temp(ii,j,l);
							xxLy2+=G[ii-i+(ncoef-1)/2]*imaLy2temp(ii,j,l);
							xxLxLy+=G[ii-i+(ncoef-1)/2]*imaLxLytemp(ii,j,l);
						}
						imaLx2(i,j,l)=xxLx2/s; imaLy2(i,j,l)=xxLy2/s; imaLxLy(i,j,l)=xxLxLy/s;
					}
				}
			}
		for (l=0; l<nblayer; l++) {                // calcul de l'image de Mc (déterminant - carré de la trace de la matrice de Harris)
			for (i=0; i<nblig; i++)
				for (j=0; j<nbcol; j++) {
					xxLx2=imaLx2(i,j,l); xxLy2=imaLy2(i,j,l); xxLxLy=imaLxLy(i,j,l);
//					imaM(i,j,l)=xxLx2*xxLy2-2*xxLxLy-kapa*(xxLx2+xxLy2)*(xxLx2+xxLy2);
					imaM(i,j,l)=xxLx2*xxLy2-xxLxLy*xxLxLy-kapa*(xxLx2+xxLy2)*(xxLx2+xxLy2);
				}
			imaM.statbasic(1);
			cout<<" imaM.moyI(l) = "<<imaM.moyI(l)<<" pow(imaM.varI(l),0.5) = "<<pow(imaM.varI(l),0.5)<<"\n";
			seuil=(float)(imaM.moyI(l)+nsig_seuil*pow(imaM.varI(l),0.5)); cout<<" seuil = "<<seuil<<"\n";
			for (i=0; i<nblig; i++) {                // elimination des non maxima locaux spatialement et calcul de l'image du laplacien 
				i0=maxi(i-1,0); i2=mini(i+1,nblig-1);
				for (j=0; j<nbcol; j++) {
					j0=maxi(j-1,0); j2=mini(j+1,nbcol-1);
					xx=imaM(i,j,l); 
					if (xx<seuil) OK=0;
					else {
						OK=1;
						for (ii=i0; ii<=i2; ii++)
							for (jj=j0; jj<=j2; jj++) if (OK && (ii!=i || jj!=j) && imaM(ii,jj,l)>xx) OK=0;
						if (OK) Tima[(i*nbcol+j)*nblayer+l]=(BYTE)OK;
					}
					ima_scale(i,j,(int)((sigma-sigmamin)/sigmapas*nblayer+l))=OK;
					xx=0; if (i>0) xx=imaLx(i,j,l)-imaLx(i-1,j,l);
					yy=0; if (j>0) yy=imaLy(i,j,l)-imaLy(i,j-1,l);
					imaLapl(i,j,(int)((sigma-sigmamin)/sigmapas*nblayer+l))=sigma2*(xx+yy);
				}
			}
		}
		if (isauv) imaM.sauve_ImaPGM("ima_M.pgm");
		if (G!=NULL) delete[] G;
		for (l=0; l<nblayer; l++) {
			for (i=0; i<nblig; i++) {
				for (j=0; j<nbcol; j++) {
					OK2=OK=(Tima[(i*nbcol+j)*nblayer+l]!=0);
					if (OK) {
						OK2=0;
						for (k=0; k<nsigmapas; k++) {
							if (!OK2 && ima_scale(i,j,k*nblayer+l)) {
								xx=imaLapl(i,j,k*nblayer+l);
								OK2=1; 
								if (k>0 && imaLapl(i,j,(k-1)*nblayer+l)>xx) OK2=0;
								if (k<nsigmapas-1 && imaLapl(i,j,(k+1)*nblayer+l)>xx) OK2=0;
								if (!OK2) {
									OK2=1;
									if (k>0 && imaLapl(i,j,(k-1)*nblayer+l)<xx) OK2=0;
									if (k<nsigmapas-1 && imaLapl(i,j,(k+1)*nblayer+l)<xx) OK2=0;
								}
							}
						}
					}
					Tima[(i*nbcol+j)*nblayer+l]=OK2;
				}
			}
		}
	}
	cout<<" sortie imacontours 4\n";
}*/

imacontours::imacontours (imadata<float> &ima, float alpha, char *filtre, bool isauv) : imasites(ima.nlig(),ima.ncol()) {
	Tima=NULL;
	nblayer=ima.ncanaux();
	valnul=0;
	Tima=new BYTE[nbpix*nblayer];
	int i,j,k;
	long int l;
	for (l=0; l<nbpix; l++) {
		for (j=0; j<nblayer; j++) Tima[l*nblayer+j]=valnul;
		Tsites[l]=&(Tima[l*nblayer]);
		}
	imadata<float> ima_grd(nblig,nbcol,nblayer),ima_dir(nblig,nbcol,nblayer);
  imadata<double> ima_gdVH(nblig,nbcol,2);
	imadata<BYTE> imaRes(nblig,nbcol,nblayer);
	double *vect_B1=new double[maxi(nblig,nbcol)];
	double *vect_B2=new double[maxi(nblig,nbcol)];
	double e_alpha=exp(-alpha), /*coef_s=1-e_alpha, */coef_c=-pow(1-e_alpha,2)/e_alpha, 
		   coef_b=pow(1-e_alpha,2)/(1+2*alpha*e_alpha-pow(e_alpha,2));
	double cflA=coef_b, cflAm1=coef_b*e_alpha*(alpha-1), cflBm1=2*e_alpha, cflBm2=-pow(e_alpha,2),
		   cflAp1=coef_b*e_alpha*(alpha+1), cflAp2=-coef_b*pow(e_alpha,2), cflBp1=cflBm1, 
		   cflBp2=cflBm2, cfdAm1=coef_c*e_alpha, cfdBm1=cflBm1, cfdBm2=cflBm2, cfdAp1=-cfdAm1,
		   cfdBp1=cflBm1, cfdBp2=cflBm2;
	for (k=0; k<nblayer; k++) {
		if (strcmp(filtre,"Deriche")==0) {
// gradient horizontal = lissage vertical puis derivation horizontale
			for (j=0; j<nbcol; j++) { // lissage vertical
				vect_B1[0]=0;
				for (i=1; i<nblig; i++) {
					vect_B1[i]=cflA*ima(i,j,k)+cflAm1*ima(i-1,j,k)+cflBm1*vect_B1[i-1];
					if (i>1) vect_B1[i]+=cflBm2*vect_B1[i-2];
					}
				vect_B2[nblig-1]=0;
				for (i=nblig-2; i>=0; i--) {
					vect_B2[i]=cflAp1*ima(i+1,j,k)+cflBp1*vect_B2[i+1];
					if (i<nblig-2) vect_B2[i]+=cflAp2*ima(i+2,j,k)+cflBp2*vect_B2[i+2];
					}
				for (i=0; i<nblig; i++) ima_gdVH(i,j,0)=vect_B1[i]+vect_B2[i];
				}
			for (i=0; i<nblig; i++) { // derivation horizontale
				vect_B1[0]=0;
				for (j=1; j<nbcol; j++) {
					vect_B1[j]=cfdAm1*ima_gdVH(i,j-1,0)+cfdBm1*vect_B1[j-1];
					if (j>1) vect_B1[j]+=cfdBm2*vect_B1[j-2];
					}
				vect_B2[nbcol-1]=0;
				for (j=nbcol-2; j>=0; j--) {
					vect_B2[j]=cfdAp1*ima_gdVH(i,j+1,0)+cfdBp1*vect_B2[j+1];
					if (j<nbcol-2) vect_B2[j]+=cfdBp2*vect_B2[j+2];
					}
				for (j=0; j<nbcol; j++) ima_gdVH(i,j,0)=vect_B1[j]+vect_B2[j];
				}
// gradient vertical = lissage horizontal puis derivation verticale
			for (i=0; i<nblig; i++) { // lissage horizontal
				vect_B1[0]=0;
				for (j=1; j<nbcol; j++) {
					vect_B1[j]=cflA*ima(i,j,k)+cflAm1*ima(i,j-1,k)+cflBm1*vect_B1[j-1];
					if (j>1) vect_B1[j]+=cflBm2*vect_B1[j-2];
					}
				vect_B2[nbcol-1]=0;
				for (j=nbcol-2; j>=0; j--) {
					vect_B2[j]=cflAp1*ima(i,j+1,k)+cflBp1*vect_B2[j+1];
					if (j<nbcol-2) vect_B2[j]+=cflAp2*ima(i,j+2,k)+cflBp2*vect_B2[j+2];
					}
				for (j=0; j<nbcol; j++) ima_gdVH(i,j,1)=vect_B1[j]+vect_B2[j];
				}
			for (j=0; j<nbcol; j++) { // derivation verticale
				vect_B1[0]=0;
				for (i=1; i<nblig; i++) {
					vect_B1[i]=cfdAm1*ima_gdVH(i-1,j,1)+cfdBm1*vect_B1[i-1];
					if (i>1) vect_B1[i]+=cfdBm2*vect_B1[i-2];
					}
				vect_B2[nblig-1]=0;
				for (i=nblig-2; i>=0; i--) {
					vect_B2[i]=cfdAp1*ima_gdVH(i+1,j,1)+cfdBp1*vect_B2[i+1];
					if (i<nblig-2) vect_B2[i]+=cfdBp2*vect_B2[i+2];
					}
				for (i=0; i<nblig; i++) ima_gdVH(i,j,1)=vect_B1[i]+vect_B2[i];
				}
			}
		else if (strcmp(filtre,"Shen")==0 || strcmp(filtre,"Shen-Castan")==0) {
// gradient horizontal
			for (i=0; i<nblig; i++) { // derivation horizontale
				vect_B1[0]=0;
				for (j=1; j<nbcol; j++) {
					vect_B1[j]=e_alpha*(vect_B1[j-1]-ima(i,j,k))+ima(i,j,k);
					}
				vect_B2[nbcol-1]=0;
				for (j=nbcol-2; j>=0; j--) {
					vect_B2[j]=e_alpha*(vect_B2[j+1]-ima(i,j,k))+ima(i,j,k);
					}
				for (j=0; j<nbcol; j++) ima_gdVH(i,j,0)=vect_B1[j]-vect_B2[j];
				}
			for (j=0; j<nbcol; j++) { // lissage vertical
				vect_B1[0]=0;
				for (i=1; i<nblig; i++) {
					vect_B1[i]=e_alpha*(vect_B1[i-1]-ima_gdVH(i,j,0))+ima_gdVH(i,j,0);
					}
				for (i=0; i<nblig; i++) ima_gdVH(i,j,0)=vect_B1[i];
				vect_B2[nblig-1]=0;
				for (i=nblig-2; i>=0; i--) {
					vect_B2[i]=e_alpha*(vect_B2[i+1]-ima_gdVH(i,j,0))+ima_gdVH(i,j,0);
					}
				for (i=0; i<nblig; i++) ima_gdVH(i,j,0)=vect_B2[i];
				}
// gradient vertical
			for (j=0; j<nbcol; j++) { // derivation verticale
				vect_B1[0]=0;
				for (i=1; i<nblig; i++) {
					vect_B1[i]=e_alpha*(vect_B1[i-1]-ima(i,j,k))+ima(i,j,k);
					}
				vect_B2[nblig-1]=0;
				for (i=nblig-2; i>=0; i--) {
					vect_B2[i]=e_alpha*(vect_B2[i+1]-ima(i,j,k))+ima(i,j,k);
					}
				for (i=0; i<nblig; i++) ima_gdVH(i,j,1)=vect_B1[i]-vect_B2[i];
				}
			for (i=0; i<nblig; i++) { // lissage horizontal
				vect_B1[0]=0;
				for (j=1; j<nbcol; j++) {
					vect_B1[j]=e_alpha*(vect_B1[j-1]-ima_gdVH(i,j,1))+ima_gdVH(i,j,1);
					}
				for (j=0; j<nbcol; j++) ima_gdVH(i,j,1)=vect_B1[j];
				vect_B2[nbcol-1]=0;
				for (j=nbcol-2; j>=0; j--) {
					vect_B2[j]=e_alpha*(vect_B2[j+1]-ima_gdVH(i,j,1))+ima_gdVH(i,j,1);
					}
				for (j=0; j<nbcol; j++) ima_gdVH(i,j,1)=vect_B2[j];
				}
			}
//			ima_gdVH.sauve_ImaBSQ("imagrdVH_Shen.dat");
// norme et direction du gradient
//		cout<<" calcul norme et orientation gradient\n";
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				ima_grd(i,j,k)=(float)(pow(pow(ima_gdVH(i,j,0),2.)+pow(ima_gdVH(i,j,1),2.),0.5));
				if (ima_gdVH(i,j,0)==0.) ima_dir(i,j,k)=(float)(PI/2.);
				else ima_dir(i,j,k)=(float)atan(ima_gdVH(i,j,1)/ima_gdVH(i,j,0));
			}
		}
	if (vect_B1!=NULL) delete[] vect_B1;
	if (vect_B2!=NULL) delete[] vect_B2;
	if (isauv) {
/*		string nomfich="imagrd_"; nomfich+=filtre; nomfich+=".dat";
		cout<<" image de la norme du gradient dans "<<nomfich<<"\n";
		ima_grd.imaunsignedchar(1).sauve_ImaBSQ(nomfich);*/
		char nomfich[80]; strcpy_s(nomfich,"imagrd_"); strcat_s(nomfich,filtre); strcat_s(nomfich,".pgm");
		cout<<" image de la norme du gradient dans "<<nomfich<<"\n";
		ima_grd.imaunsignedchar().sauve_ImaPGM(nomfich);
		strcpy_s(nomfich,"imadir_"); strcat_s(nomfich,filtre); strcat_s(nomfich,".pgm");
		cout<<" image de la norme du gradient dans "<<nomfich<<"\n";
		((ima_dir*(180.f/(float)PI))+90.f).imaunsignedchar().sauve_ImaPGM(nomfich);
		}
	ima_grd.statbasic(1);
// seuillage a hysteresis
	for (k=0; k<nblayer; k++) {
		float sh=(float)ima_grd.moyI(k)+(float)sqrt(ima_grd.varI(k));
		float sb=(float)ima_grd.moyI(k);
//		ima_grd.seuil_ima (sb,imaRes,k);
		ima_grd.seuil_hysteresis (sh,sb,imaRes,k);
		}
	for (k=0; k<nblayer; k++)
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) (*this)(i,j,k)=imaRes(i,j,k);
// affinement des contours en ne gardant que max locaux dans la direction du gradient
//	imagrd_maxloc (ima_grd,ima_dir);
	}

//void imacontours::imagrd_maxloc (imadata<float> & ima_grd, imadata<float> & ima_dir,
//							 	 imadata<BYTE> & imaRes) {
void imacontours::imagrd_maxloc (imadata<float> & ima_grd, imadata<float> & ima_dir) {
	int i,j,k;
	float orient, norm, coef, x1, x2, x;
	for (k=0; k<nblayer; k++)
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
//				bool icontour=!!imaRes(i,j,k);
//				bool icontour=!!(*this)(i,j,k);
				bool icontour=((*this)(i,j,k)>0?1:0);
				if (icontour) {
					orient=ima_dir(i,j,k);
					norm=ima_grd(i,j,k);
// calcul grossier des positions des voisins dans la direction du gradient
/*					if (orient<=-3*PI/8. || orient>3*PI/8.) {
						if (i>0 && ima_grd(i-1,j,k)>norm) icontour=0;
						if (i<nblig-1 && ima_grd(i+1,j,k)>norm) icontour=0;
						}
					else {
						if (orient<=-PI/8.) {
							if (i>0 && j>0 && ima_grd(i-1,j-1,k)>norm) icontour=0;
							if (i<nblig-1 && j<nbcol-1 && ima_grd(i+1,j+1,k)>norm) icontour=0;
							}
						else {
							if (orient<=PI/8.) {
								if (j>0 && ima_grd(i,j-1,k)>norm) icontour=0;
								if (j<nbcol-1 && ima_grd(i,j+1,k)>norm) icontour=0;
								}
							else {
								if (i>0 && j<nbcol-1 && ima_grd(i-1,j+1,k)>norm) icontour=0;
								if (i<nblig-1 && j>0 && ima_grd(i+1,j-1,k)>norm) icontour=0;
								}
							}
						}*/
// calcul par interpolation linéaire des valeurs des voisins dans la direction du gradient
					x1=-1; x2=-1;
					if (orient>=0) {
						if (orient<=PI/4.) {
							coef=tan(orient); 
							if (coef<0 || coef>1) coef=mini(1.f,maxi(0.f,coef));
							if (i>0 && j<nbcol-1) x1=(1-coef)*ima_grd(i,j+1,k)+coef*ima_grd(i-1,j+1,k);
							if (i<nblig-1 && j>0) x2=(1-coef)*ima_grd(i,j-1,k)+coef*ima_grd(i+1,j-1,k);
							}
						else {
							coef=(float)tan(PI/2.-orient);
							if (coef<0 || coef>1) coef=mini(1.f,maxi(0.f,coef));
							if (i>0 && j<nbcol-1) x1=(1-coef)*ima_grd(i-1,j,k)+coef*ima_grd(i-1,j+1,k);
							if (i<nblig-1 && j>0) x2=(1-coef)*ima_grd(i+1,j,k)+coef*ima_grd(i+1,j-1,k);
							}
						}
					else {
						if (orient>=-PI/4.) {
							coef=tan(-orient);
							if (coef<0 || coef>1) coef=mini(1.f,maxi(0.f,coef));
							if (i<nblig-1 && j<nbcol-1) x1=(1-coef)*ima_grd(i,j+1,k)+coef*ima_grd(i+1,j+1,k);
							if (i>0 && j>0) x2=(1-coef)*ima_grd(i,j-1,k)+coef*ima_grd(i-1,j-1,k);
							}
						else {
							coef=(float)tan(PI/2.+orient);
							if (coef<0 || coef>1) coef=mini(1.f,maxi(0.f,coef));
							if (i<nblig-1 && j<nbcol-1) x1=(1-coef)*ima_grd(i+1,j,k)+coef*ima_grd(i+1,j+1,k);
							if (i>0 && j>0) x2=(1-coef)*ima_grd(i-1,j,k)+coef*ima_grd(i-1,j-1,k);
							}
						}
					x=maxi(x1,x2); x2=mini(x1,x2); x1=x;
					if (x1<0 || x1>norm) icontour=0;
					else {
						if (x1=norm) {
							if (x2<0 || (x2>=0 && x2>=norm)) icontour=0;
							}
						else {
							if (x2>norm) icontour=0;
							}
						}
// affectation à l'image résultat
					(*this)(i,j,k)=icontour;
					}
				}
	}

void imacontours::prolongecontours(imadata<float> &ima_grd, unsigned short int lgmax, float seuil_grd) {cout<<" entree dans prolongecontours\n";
	int i,j,k,l,m,n,p,ii,jj,kk,max,maxt,i0,j0,n0,pf,mf; cout<<" lgmax = "<<lgmax<<", seuil_grd = "<<seuil_grd<<"\n";
	unsigned int prof0;
	double cout0, coutf;
	imacontours imabis(*this);
// construction de la "look-up table"
	int T[256][7];
	for (n=0; n<256; n++) {T[n][0]=0; for (j=1; j<7; j++) T[n][j]=-1;}
	for (i=0; i<8; i++) {
		n=(int)pow(2.f,i);
		T[n][0]=1;
		switch (i) { 
			case 0 : T[n][1]=-1;T[n][2]=-1;T[n][3]= 0;T[n][4]=-1;T[n][5]=+1;T[n][6]=-1;break;
			case 1 : T[n][1]= 0;T[n][2]=-1;T[n][3]=+1;T[n][4]=-1;T[n][5]=+1;T[n][6]= 0;break;
			case 2 : T[n][1]=+1;T[n][2]=-1;T[n][3]=+1;T[n][4]= 0;T[n][5]=+1;T[n][6]=+1;break;
			case 3 : T[n][1]=+1;T[n][2]= 0;T[n][3]=+1;T[n][4]=+1;T[n][5]= 0;T[n][6]=+1;break;
			case 4 : T[n][1]=+1;T[n][2]=+1;T[n][3]= 0;T[n][4]=+1;T[n][5]=-1;T[n][6]=+1;break;
			case 5 : T[n][1]= 0;T[n][2]=+1;T[n][3]=-1;T[n][4]=+1;T[n][5]=-1;T[n][6]= 0;break;
			case 6 : T[n][1]=-1;T[n][2]=+1;T[n][3]=-1;T[n][4]= 0;T[n][5]=-1;T[n][6]=-1;break;
			case 7 : T[n][1]=-1;T[n][2]= 0;T[n][3]=-1;T[n][4]=-1;T[n][5]= 0;T[n][6]=-1;break;
		}
	}
// algorithme de fermeture
	const int max_noeuds=(int)pow(3.f,2); cout<<" fermeture sur un distance max de lgmax = "<<lgmax<<"\n";
	for (k=0; k<nblayer; k++) {
		float sb=seuil_grd, bonus;
		if (sb<0) sb=maxi((float)ima_grd.moyI(k)/*-sqrt((float)ima_grd.varI(k))*/,0.f);
		bonus=2*sb; cout<<" sb = "<<sb<<", bonus = "<<bonus<<"\n";
		for (i=1; i<nblig-1; i++)
			for (j=1; j<nbcol-1; j++) {
				if ((*this)(i,j,k)==1) {
					n=(*this)(i,j+1,k)+2*(*this)(i-1,j+1,k)+4*(*this)(i-1,j,k)+8*(*this)(i-1,j-1,k);
					n+=16*(*this)(i,j-1,k)+32*(*this)(i+1,j-1,k)+64*(*this)(i+1,j,k)+128*(*this)(i+1,j+1,k);
					if (T[n][0]==1) { cout<<" pixel ("<<i<<","<<j<<") configuration du voisinage : n = "<<n<<"\n"; //{char aa; cin>>aa;}
						noeud **arbre=new noeud*[lgmax+1];
						for (p=0; p<=lgmax; p++) arbre[p]=new noeud[mini((int)pow(3.f,p),max_noeuds)]; // on ne va garder que les noeuds les plus intéressants
						for (p=0; p<=lgmax; p++) 
							for (m=0; m<mini((int)pow(3.f,p),max_noeuds); m++) {arbre[p][m].stop=1; arbre[p][m].x=arbre[p][m].y=arbre[p][m].p=-1; arbre[p][m].cout=-1.;}
						arbre[0][0].x=i; arbre[0][0].y=j; arbre[0][0].code=n; arbre[0][0].prof=0; arbre[0][0].cout=0.; arbre[0][0].stop=0; // racine
						p=0; bool arbrefini=0;
						while (p<lgmax && !arbrefini) {
//							max=mini((int)pow(3.f,p),max_noeuds); cout<<p<<" max = "<<max<<"\n";
							max=mini((int)pow(3.f,p),max_noeuds); maxt=3*max; mf=-1; //cout<<p<<" max = "<<max<<" 3*max = "<<maxt<<"\n";
							noeud *noeudstemp=new noeud[maxt];
							for (m=0; m<max; m++) {
								if (arbre[p][m].stop==0) {
									i0=arbre[p][m].x; j0=arbre[p][m].y; n0=arbre[p][m].code;
									prof0=arbre[p][m].prof; pf=p+1; cout0=arbre[p][m].cout; //cout<<" *** examen noeud terminal numero "<<m<<" en ("<<i0<<","<<j0<<")\n";
									for (l=0; l<3; l++) {
										ii=i0+T[n0][2*l+1]; jj=j0+T[n0][2*l+2];
//										mf=3*m+l;
//										mf++;	arbre[pf][mf].x=ii;	arbre[pf][mf].y=jj;	arbre[pf][mf].prof=prof0+1; arbre[pf][mf].stop=0;
										mf++;	noeudstemp[mf].x=ii;	noeudstemp[mf].y=jj;	noeudstemp[mf].p=m; noeudstemp[mf].prof=prof0+1; noeudstemp[mf].stop=0;
//										cout<<" examen candidat suivant numero "<<mf<<" en ("<<ii<<","<<jj<<")\n";
//										if (ii>=0 && ii<nblig-1 && jj>=0 && jj<nbcol-1) {
										if (ii>0 && ii<nblig-1 && jj>0 && jj<nbcol-1) { // bords de l'image ????????
											if (T[n0][2*l+1]==+1) {
//												arbre[pf][mf].code=(BYTE)pow(2.f,2+T[n0][2*l+2]);
												noeudstemp[mf].code=(BYTE)pow(2.f,2+T[n0][2*l+2]);
//												if (T[n0][2*l+2]==+1) arbre[p+1][m*3+l].code=pow(2,3);
//												if (T[n0][2*l+2]== 0) arbre[p+1][m*3+l].code=pow(2,2);
//												if (T[n0][2*l+2]==-1) arbre[p+1][m*3+l].code=pow(2,1);
											}
											if (T[n0][2*l+1]==0) {
//												arbre[pf][mf].code=(BYTE)pow(2.f,2+2*T[n0][2*l+2]);
												noeudstemp[mf].code=(BYTE)pow(2.f,2+2*T[n0][2*l+2]);
//												if (T[n0][2*l+2]==+1) arbre[p+1][m*3+l].code=pow(2,4);
//												if (T[n0][2*l+2]==-1) arbre[p+1][m*3+l].code=pow(2,0);
											}
											if (T[n0][2*l+1]==-1) {
//												arbre[pf][mf].code=(BYTE)pow(2.f,6-T[n0][2*l+2]);
												noeudstemp[mf].code=(BYTE)pow(2.f,6-T[n0][2*l+2]);
//												if (T[n0][2*l+2]==+1) arbre[p+1][m*3+l].code=pow(2,5);
//												if (T[n0][2*l+2]== 0) arbre[p+1][m*3+l].code=pow(2,6);
//												if (T[n0][2*l+2]==-1) arbre[p+1][m*3+l].code=pow(2,7);
											}
											coutf=(cout0*prof0+ima_grd(ii,jj,k));
//											cout<<" cout = "<<cout0<<"*"<<prof0<<"+"<<ima_grd(ii,jj,k)<<" a diviser par "<<noeudstemp[mf].prof<<" = "<<coutf/noeudstemp[mf].prof<<"\n";
											noeudstemp[mf].cout=coutf/noeudstemp[mf].prof;
//											arbre[pf][mf].cout=coutf/arbre[pf][mf].prof;
//											if ((*this)(ii,jj,k)==1) arbre[pf][mf].stop=1;
											if ((*this)(ii,jj,k)==1) {noeudstemp[mf].stop=1; noeudstemp[mf].cout+=bonus; 
												cout<<" !!!!!!!!!!!!!! pixel ("<<ii<<","<<jj<<") rejoint -> bonus !!!!!!!!!!!!!!! -> "<<noeudstemp[mf].cout<<"\n";}
										}
//										else arbre[pf][mf].stop=1;
										else noeudstemp[mf].stop=1;
									}
								} else {cout<<" &&& examen noeud terminal numero "<<m<<" en ("<<arbre[p][m].x<<","<<arbre[p][m].y<<") stoppe\n";}
							}
//							for (l=0; l<max; l++) if (!arbre[pf][l].stop) cout<<" || c = "<<arbre[pf][l].cout; cout<<"\n";
//							for (l=0; l<=mf; l++) cout<<" || c = "<<noeudstemp[l].cout; cout<<"\n";
//							kk=0; for (l=0; l<max; l++) kk+=1-arbre[pf][l].stop;
//							bool icoupe=(3*kk>=max_noeuds?1:0); cout<<" niveau "<<pf<<", "<<kk<<" branches restent : coupe ? "<<(int)icoupe<<"\n";
							bool icoupe=(mf>=max_noeuds?1:0); //cout<<" niveau "<<pf<<", "<<l<<" branches : coupe ? "<<(int)icoupe<<"\n";
							if (!icoupe) {for (l=0; l<=mf; l++) arbre[pf][l]=noeudstemp[l];}
							else {
//								while (icoupe) {
//									double *tab_cout=new double[max];
//									kk=0;
//									for (l=0; l<max; l++) if (!arbre[pf][l].stop) tab_cout[kk++]=arbre[pf][l].cout; cout<<kk<<" branches -> ";
//									tri_rapide(tab_cout,kk,0); // tri par ordre croissant
//									double cout_min=tab_cout[2*kk/3]; cout<<" on coupe les "<<2*kk/3<<" branches de cout inferieur a "<<cout_min;
//									for (l=0; l<max; l++) if (!arbre[pf][l].stop && arbre[pf][l].cout<=cout_min) arbre[pf][l].stop=1;
//									kk=0; for (l=0; l<max; l++) kk+=1-arbre[pf][l].stop; cout<<" -> "<<kk<<" branches restantes\n";
//									if (tab_cout!=NULL) delete[] tab_cout;
//									if (3*kk<max_noeuds) icoupe=0;
//								}
								double *tab_cout=new double[mf+1];
								for (kk=0; kk<=mf; kk++) tab_cout[kk]=noeudstemp[kk].cout; //cout<<kk<<" branches -> ";
								tri_rapide(tab_cout,mf+1,0); // tri par ordre croissant
								int argseuil=mf-max_noeuds;
								double cout_min=tab_cout[argseuil]; //cout<<" on coupe les "<<argseuil<<" branches de cout inferieur a "<<cout_min;
								l=0;
								for (kk=0; kk<=mf; kk++) 
									if (noeudstemp[kk].cout>cout_min) {arbre[pf][l]=noeudstemp[kk]; l++;}
								if (tab_cout!=NULL) delete[] tab_cout;
							}
							if (noeudstemp!=NULL) delete[] noeudstemp; 
							kk=0; for (l=0; l<mini(maxt,max_noeuds); l++) kk+=1-arbre[pf][l].stop;
							p++; if (kk==0) arbrefini=1;
						}
						coutf=0; pf=0; mf=0;
						for (p=0; p<=lgmax; p++) 
							for (m=0; m<mini((int)pow(3.f,p),max_noeuds); m++) //{
//								cout<<" noeud "<<p<<" "<<m<<" ("<<arbre[p][m].x<<","<<arbre[p][m].y<<") "<<arbre[p][m].cout<<"\n";
								if (arbre[p][m].cout>coutf) {pf=p; mf=m; coutf=arbre[p][m].cout; }
//							}
						cout<<" depuis ("<<i<<","<<j<<"), meilleur cout "<<coutf<<" en "<<pf<<" "<<mf<<" >? sb "<<sb<<"\n";
						sb=0;
						if (coutf>sb) {cout<<" chemin en sens inverse : ";
							m=mf;
							for (p=pf; p>0; p--) {
								ii=arbre[p][m].x;	jj=arbre[p][m].y;	imabis(ii,jj,k)=1; n0=arbre[p][m].code; 
								cout<<" prof = "<<p<<" : ("<<ii<<","<<jj<<") -> ";
								switch (n0) {
									case 1: i0=ii; j0=jj+1; break; case 2: i0=ii-1; j0=jj+1; break;
									case 4: i0=ii-1; j0=jj; break; case 8: i0=ii-1; j0=jj-1; break;
									case 16: i0=ii; j0=jj-1; break; case 32: i0=ii+1; j0=jj-1; break;
									case 64: i0=ii+1; j0=jj; break; case 128: i0=ii+1; j0=jj+1; break;
									default: cout<<" Pb dans decodage de direction d'arrivee : "<<n0<<" ??? \n"; char aa; cin>>aa; break;
								}
								bool trouve=0; kk=0;
								while (!trouve && kk<mini((int)pow(3.f,p-1),max_noeuds)) {
									if (arbre[p-1][kk].x==i0 && arbre[p-1][kk].y==j0) {m=kk; trouve=1;}
									kk++;} cout<<" | ";
								if (!trouve) {cout<<" Pb : predecesseur non trouve ????\n"; char aa; cin>>aa; m=arbre[p][m].p;}
//								else {
//									if (m!=arbre[p][mf].p) {cout<<" conflit dans predecesseurs : "<<m<<" != "<<arbre[p][mf].p<<" ????\n"; char aa; cin>>aa; }
//									else mf=m;}
							}
						}
						for (p=0; p<=lgmax; p++) if (arbre[p]!=NULL) delete[] arbre[p];
						if (arbre!=NULL) delete[] arbre;
					}
				}
		}
	}
	(*this)=imabis;
}

void imacontours::edgedrawing(imadata<float> &ima_don, unsigned short int lgmax, unsigned short int kfact,
															float th, bool i_ima_grd) {cout<<" entree dans edgedrawing\n";
	int i,j,n,i0,j0,i2,j2,ii,jj,idir,k; cout<<" lgmax = "<<lgmax<<", detail ratio = "<<kfact<<"\n";
	const double cos45=1./pow(2.,0.5);
	const bool ifiltre2=1;
	float x0,x1,x2,xmax;
	imadata<float> ima_dir, ima_grd=ima_don.gradient(ima_dir,"Sobel"); ima_grd.statbasic(1);
	if (i_ima_grd) {ima_grd=ima_don; ima_dir=ima_dir.filtremedian(7,7); }
	imabin ima_anchor(nblig,nbcol), imab(ima_grd,(float)ima_grd.moyI());
	if (th>=0) imab=imabin(ima_grd,th);
	else {float s=ima_grd.seuil_otsu (imab); cout<<" seuil otsu norme du gradient = "<<s<<"\n"; }
	//imab=imab.fermeture (eltstruct(lgmax/2,lgmax/2));
	imab.imaunsignedchar(1).sauve_ImaPGM("ima_grd_seuillee.pgm"); if (!ifiltre2) imab.mise_a_un();
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) {
			if ((i%kfact==0 || j%kfact==0) && (*this)(i,j)>0 && imab(i,j)) ima_anchor(i,j)=1; //else ima_anchor(i,j)=0;
		}
	ima_anchor.imaunsignedchar(1).sauve_ImaPGM("ima_anchor.pgm"); 
	imadata<BYTE> ima_dir_filtree(nblig,nbcol);
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++)
			if (imab(i,j)) {if (fabs(cos(ima_dir(i,j)))>cos45) ima_dir_filtree(i,j)=2; else ima_dir_filtree(i,j)=1;}
	if (ifiltre2) {
		imadata<BYTE> ima_dir_filtree2(ima_dir_filtree);
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) 
				if (ima_dir_filtree2(i,j)>0) {
					i0=maxi(i-1,0); i2=mini(i+1,nblig-1); j0=maxi(j-1,0); j2=mini(j+1,nbcol-1); x0=x1=x2=0.f;
					for (ii=i0; ii<=i2; ii++) for (jj=j0; jj<=j2; jj++) {if (ima_dir_filtree2(ii,jj)==1) x1+=1; if (ima_dir_filtree2(ii,jj)==2) x2+=1; }
					if (x1>x2) {ima_dir_filtree(i,j)=1;} else {if (x1<x2) ima_dir_filtree(i,j)=2;}
				}
	}
	ima_dir_filtree.sauve_ImaPGM("ima_dir_filtree.pgm");
// algorithme pour relier les anchors
	imabin imacontbis(ima_anchor);
	unsigned int tab_coord[50][2], isuiv;
	bool fini, trouve;
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++)
			if (ima_anchor(i,j)) {cout<<" anchor ("<<i<<","<<j<<") ";
				idir=ima_dir_filtree(i,j); // 4 directions possibles: 1=Ouest, 2=Nord, 3=Est, 4=Sud
				fini=trouve=0; n=0; i0=i; j0=j; 
				do {
					x0=x1=x2=0.f;
					switch (idir) {
						case 1: if (j0>0) {x1=ima_grd(i0,j0-1); if (i0>0) x0=ima_grd(i0-1,j0-1); if (i0<nblig-1) x2=ima_grd(i0+1,j0-1); } break;
						case 2: if (i0>0) {x1=ima_grd(i0-1,j0); if (j0>0) x0=ima_grd(i0-1,j0-1); if (j0<nbcol-1) x2=ima_grd(i0-1,j0+1); } break;
						case 3: if (j0<nbcol-1) {x1=ima_grd(i0,j0+1); if (i0>0) x0=ima_grd(i0-1,j0+1); if (i0<nblig-1) x2=ima_grd(i0+1,j0+1); } break;
						case 4: if (i0<nblig-1) {x1=ima_grd(i0+1,j0); if (j0>0) x0=ima_grd(i0+1,j0-1); if (j0<nbcol-1) x2=ima_grd(i0+1,j0+1); } break;
					}
					xmax=maxi(x0,maxi(x1,x2)); isuiv=99;
					if (isuiv>2 && xmax==x1) isuiv=1; if (isuiv>2 && xmax==x0) isuiv=0; if (isuiv>2 && xmax==x2) isuiv=2;
					if (xmax>0 && isuiv<=2) {
						switch (idir) {
							case 1: j0=j0-1; i0=i0-1+isuiv; break; case 2: i0=i0-1; j0=j0-1+isuiv; break;
							case 3: j0=j0+1; i0=i0-1+isuiv; break; case 4: i0=i0+1; j0=j0-1+isuiv; break;
						}
						if (i0<0 || j0<0 || i0>=nblig || j0>=nbcol) cout<<" i0 = "<<i0<<" j0 = "<<j0<<"\n";
						tab_coord[n][0]=i0; tab_coord[n][1]=j0; n++;
						if (ima_anchor(i0,j0)) trouve=1;
						if (trouve) {
							for (k=0; k<n; k++) imacontbis(tab_coord[k][0],tab_coord[k][1])=1; fini=1; cout<<" prolongation de "<<n<<" pixels\n";
						}
					} else fini=1;
					if (n>lgmax) fini=1;
					if (fini && idir>0 && idir<3) {idir+=2; i0=i; j0=j; n=0; fini=trouve=0;}
				} while (!fini);
			}
	(*this)=imacontours(imacontbis);
	imadata<BYTE> imacol(nblig,nbcol,3); imacol.copiecanal(0,imacontbis.imaunsignedchar()); 
	imacol.copiecanal(1,ima_anchor.imaunsignedchar(1)); imacol.sauve_ImaPGM("ima_edgedraw.ppm");
}

void imacontours::sauve_Ima(char *nomfich, int icanal) const {
	ofstream sortie;
	sortie.open(nomfich,ios::out|ios::binary);
	if (!sortie) {
		cout<<" ouverture de "<<nomfich<<" impossible\n";
		exit (-1);}
	BYTE* adval=(BYTE*)Tsites[0];
	BYTE val=0;
	int itype, sizeval;
	itype=0;
	sizeval=sizeof(BYTE);
//	cout<<" image "<<nblig<<" lig.& "<<nbcol<<" col., de type ";
//	cout<<" bytes => itype = "<<itype<<"\n";
	sortie.write((char *)&nblig,sizeof(int));
	sortie.write((char *)&nbcol,sizeof(int));
	sortie.write((char *)&itype,sizeof(int));
	for (long int i=0; i<nbpix; i++) {
		adval=(BYTE*)Tsites[i]+icanal;
		val=(*adval)*255;
		sortie.write((char *)&val,sizeval);
	}
	sortie.close();
}

keypoints::keypoints (imadata<float> &ima, float gamma, float threshold, bool isauv, bool iaffich) {
	cout<<" entree dans detecteur de Harris de la classe 'keypoints'\n";
//	float sigma=1.6; // fct de gamma ????
	float sigma=gamma;
	n_key_points_max=0; n_key_points=0; lg_descr=0;
	const float kappa=0.04f/*0.15f*/; 
	const unsigned int i_calc_deriv=1; // 0 pour différences finies, 1 pour filtrage linéaire passe haut, 2 pour filtrage optimal...
	int nblig=ima.nlig(), nbcol=ima.ncol(), nb_ech=ima.ncanaux(); // sous-entend que ima contient les differents filtrages gaussiens de l'image initiale 
	imadata<float> imasauv(ima), imaIxIy(nblig,nbcol,3*nb_ech), imaMc(nblig,nbcol,nb_ech); 
	imaIxIy.mise_a_zero(); imaMc.mise_a_zero(); 
	float Lx,Ly,xx;
	int i,j,k,n,i0,i2,ii,j0,j2,jj;

/************************** debut a programmer ****************************/




/***************************  fin a programmer ****************************/

	imadata<BYTE> ima_keypoints_Harris(nblig,nbcol); ima_keypoints_Harris.mise_a_zero();
	bool extremI, fini;
	ima_keypoints_Harris.imaunsignedchar(0).sauve_ImaPGM("ima_keypoints_Harris.pgm");
	cout<<" sauvegarde des "<<n_key_points<<" points-cles\n";
	n_key_points_max=n_key_points;
	lg_descr=0;
	Tdescr_key_points=new float*[n_key_points_max];
	for (n=0; n<n_key_points_max; n++) Tdescr_key_points[n]=new float[lg_descr];
	ima=imasauv;
}