#include "imaregions.h"

/*const unsigned int nmaxarreteparsommet=6000;
struct sommet {
	liste_pixels Lpix;
	double val;
	unsigned int n_arretes;
	unsigned int T_arretes[nmaxarreteparsommet];
	bool valid;
};

struct arrete {
	unsigned long int pred, succ;
	double val, length;
	bool valid;
};

struct region {
	liste_pixels Lpix;
	double val;
//	unsigned long int npixB;
	unsigned int n_arretes;
	unsigned long int *T_arretes;
	bool valid;
};*/

imaregions::imaregions (imabin &imab, bool ifond0) : imasites(imab.nlig(),imab.ncol()) {
	nblayer=1;
	Tima=new long int[nbpix*nblayer];
	valnul=/*-1*/0;
	int i,j;
	for (long int l=0; l<nbpix; l++) {
		for (j=0; j<nblayer; j++) Tima[l*nblayer+j]=valnul;
		Tsites[l]=&(Tima[l*nblayer]);
	}
	int ncc=0;
	imadata<int> imacc=imab.negatif().composantes_connexes (ncc);
	if (ncc>1) {
		cout<<" fond compose de "<<ncc<<" composantes connexes !\n";
		if (ifond0) {
			cout<<"confirmez-vous de tout mettre sous le label 0 ?"; bool ifondX; cin>>ifondX;
			ifond0=!ifondX;
		}
	}
	nbregions=1;
	if (ifond0) {
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) if (!imab(i,j)) (*this)(i,j)=0;
		nbregions=0;
	} else {
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++)
				if (!imab(i,j)) (*this)(i,j)=imacc(i,j)+nbregions;
		nbregions+=ncc;
	}
	imacc=imab.composantes_connexes (ncc);
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++)
			if (imab(i,j)) (*this)(i,j)=imacc(i,j)+nbregions;
	nbregions+=ncc;
	cout<<nbregions<<" regions construites \n";
	s_homog=-1;
}

imaregions::imaregions (imalabels &imalab, int k) : imasites(imalab.nlig(),imalab.ncol()) {
	nblayer=1;
	Tima=new long int[nbpix*nblayer];
	valnul=/*-1*/0;
	int i,j,l=0,n0;
	for (long int ii=0; ii<nbpix; ii++) {
		for (j=0; j<nblayer; j++) Tima[ii*nblayer+j]=valnul;
		Tsites[ii]=&(Tima[ii*nblayer]);
	}
	s_homog=-1;
	int nclas=0,ncc=0,n;
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++)
			if (imalab(i,j,k)>nclas) nclas=imalab(i,j,k);
	imabin imal(nblig,nbcol);
	nbregions=0; n0=1-valnul;
	for (l=1; l<=nclas; l++) {
		n=0;
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				if (imalab(i,j,k)==l) {imal(i,j)=1; n++;}
				else imal(i,j)=0;
			}
		if (n>0) {
			cout<<" classe "<<l<<" : "<<n<<"pixels\n";
			imadata<int> imacc=imal.composantes_connexes (ncc);
			for (i=0; i<nblig; i++)
				for (j=0; j<nbcol; j++)
					if (imalab(i,j,k)==l) (*this)(i,j)=imacc(i,j)-1+n0;
			nbregions+=ncc; n0+=ncc;
			cout<<nbregions<<" regions construites a partir des "<<l<<" premieres classes\n";
		}
	}
	cout<<nbregions<<" regions construites a partir des "<<nclas<<" classes\n";
}

imaregions::imaregions (imaregions &imareg1, imaregions &imareg2) : imasites(imareg1.nlig(),imareg1.ncol()) {
	nblayer=1; s_homog=-1;
	Tima=new long int[nbpix*nblayer];
	valnul=/*-1*/0;
	int i,j;
	for (i=0; i<nbpix; i++) {
		for (j=0; j<nblayer; j++) Tima[i*nblayer+j]=valnul;
		Tsites[i]=&(Tima[i*nblayer]);
	}
//	imareg1.verif_1cc_per_region(); imareg2.verif_1cc_per_region();
	int nreg1=imareg1.nbregions, nreg2=imareg2.nbregions, n=0, k, ncc,ii,jj,kk;
	cout<<" nreg 1 = "<<nreg1<<", nreg2 = "<<nreg2<<"\n";
	int **T_corresp=new int*[nreg1+1];
	for (i=0; i<=nreg1; i++) {T_corresp[i]=new int[nreg2+1]; for (j=0; j<=nreg2; j++) T_corresp[i][j]=-1;}
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) if (T_corresp[imareg1(i,j)][imareg2(i,j)]<0) T_corresp[imareg1(i,j)][imareg2(i,j)]=++n;
	cout<<n<<" intersections de regions non vides\n";
	nbregions=n+1;
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) (*this)(i,j)=T_corresp[imareg1(i,j)][imareg2(i,j)];
	imabin imab(nblig,nbcol);
	imadata<int> imacc(nblig,nbcol);
	for (k=0; k<=n; k++) {
		imab.mise_a_zero();
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) if ((*this)(i,j)==k) imab(i,j)=1;
		if (imab.norm()>0) {
			imacc=imab.composantes_connexes(ncc);
			if (ncc>1) {
				for (kk=2; kk<=ncc; kk++) {
					nbregions++;
					for (ii=0; ii<nblig; ii++)
						for (jj=0; jj<nbcol; jj++) if (imacc(ii,jj)==kk) (*this)(ii,jj)=nbregions;
				}
			}
		}
	}
	cout<<" finalement "<<nbregions<<" regions generees apres etiquettage en composantes connexes des intersections non vides\n";
	verif_1cc_per_region();
	for (i=0; i<nreg1; i++) if (T_corresp[i]!=NULL) delete[] T_corresp[i];
	if (T_corresp!=NULL) delete[] T_corresp;
	cout<<nbregions<<" regions construites a partir des "<<nreg1<<" regions de l'image 1 et des "<<nreg2<<" regions de l'image 2\n";
}

bool imaregions::verif_homog(float x1, float x2) {
	bool b=0;
	if (abs(x1-x2)<=s_homog) b=1;
	return b;
}

imaregions::imaregions (imadata<float> &ima, float s, char *selectgermes, bool iaf) : imasites(ima.nlig(),ima.ncol()) {
	nblayer=1;
	Tima=new long int[nbpix*nblayer];
	valnul=-1;
	int i,j,i0,i2,j0,j2;
	long int l;
	for (l=0; l<nbpix; l++) {
		for (j=0; j<nblayer; j++) Tima[l*nblayer+j]=valnul;
		Tsites[l]=&(Tima[l*nblayer]);
	}
	s_homog=s;
	bool fini=0;
	imabin imadeja(nblig,nbcol), imaLpix(nblig,nbcol); imadeja.mise_a_zero();
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) if (ima(i,j)==valnul) imadeja(i,j)=1; //cout<<imadeja.norm()<<"\n";
	int x0, y0, nreg=1;
	long int npixrest=nbpix;
	float x;
	char *selectH="histo", *selectA="aleat";
	cout<<" segmentation par croissance de regions\n";
	while (npixrest>0) {
		if (iaf) cout<<" npixrest "<<npixrest<<"\n";
		x0=-1; y0=-1;
		if (strcmp(selectgermes,selectH)==0) {
			imadata<float> imabis(ima);
			for (i=0; i<nblig; i++)
				for (j=0; j<nbcol; j++) if (imadeja(i,j)) imabis(i,j)=imabis.v_nulle();
			double prec;
			imabis.histogramme(imadeja.negatif(),(int)(256/s_homog),1);
			x=(float)imabis.modeI(prec);
			if (iaf) cout<<" mode des pixels restants = "<<x<<" a "<<prec<<" pres\n";
			float vmin=1.e+9, v;
			int ii, jj, nv;
			for (i=0; i<nblig; i++)
				for (j=0; j<nbcol; j++)
					if (imabis(i,j)!=imabis.v_nulle() && abs(imabis(i,j)-x)<prec) {
						v=0; nv=0;
						for (ii=maxi(i-1,0); ii<=mini(i+1,nblig-1); ii++)
							for (jj=maxi(j-1,0); jj<=mini(j+1,nbcol-1); jj++)
								if (imabis(ii,jj)!=imabis.v_nulle()) { 
									v+=abs(ima(ii,jj)-x);
									nv++;
								}
						if (nv>1) v=v/(nv-1);
						if (v<vmin) {
							x0=i; y0=j;
							vmin=v;
						}
				}
			if (iaf) cout<<" nouveau germe en ("<<x0<<","<<y0<<")\n";
		}
		if (strcmp(selectgermes,selectA)==0 || x0==-1 || y0==-1) {
			int nmaxtirages=1000000, ntirages=0;
			do {
				ntirages++;
				x0=(int)((float)rand()/RAND_MAX*nblig); x0=maxi(0,mini(x0,nblig-1));
				y0=(int)((float)rand()/RAND_MAX*nbcol); y0=maxi(0,mini(y0,nbcol-1));
			} while (imadeja(x0,y0) && ntirages<nmaxtirages);
			if (iaf) cout<<" nouveau germe en ("<<x0<<","<<y0<<")\n";
			if (iaf && ntirages>=nmaxtirages) {cout<<" pb dans le tirage aleatoire "<<imadeja.norm()<<"\n"; npixrest=0;}
		}
//		liste_pixels Lpix; imaLpix.mise_a_zero(); Lpix.insere(x0,y0); 
		int *Tpix=new int[npixrest]; Tpix[0]=x0*nbcol+y0;
		float xmoy=ima(x0,y0); imaLpix(x0,y0)=1;
		unsigned int n=1, ind_extrait=0, ind_insere=1; 
		while (ind_extrait<ind_insere) {//cout<<imadeja.norm()<<"\n";
//		while (Lpix.nb_elts()>0) {//cout<<imadeja.norm()<<"\n";
//			elt_liste E=Lpix.extrait(); int ix=E.x, iy=E.y;
			int ix=Tpix[ind_extrait]/nbcol, iy=Tpix[ind_extrait]%nbcol; ind_extrait++;
			imaLpix(ix,iy)=0; //if (nreg>=2) {cout<<" pixel "<<ix<<" "<<iy<<" extrait\n"; }
			(*this)(ix,iy)=nreg; imadeja(ix,iy)=1; npixrest--;
			x=ima(ix,iy);
			i0=maxi(0,ix-1); i2=mini(ix+1,nblig-1); j0=maxi(0,iy-1); j2=mini(iy+1,nbcol-1);
			for (i=i0; i<=i2; i++)
				for (j=j0; j<=j2; j++) { 
//					if (nreg==2) cout<<i<<" "<<j<<" "<<ima(i,j)<<" "<<xmoy<<" "<<verif_homog(ima(i,j),xmoy)<<" "<<!imaLpix(i,j)<<" "<<!imadeja(i,j)<<" "<<ima(i,j)<<"\n";
					if (!imaLpix(i,j) && !imadeja(i,j) && verif_homog(ima(i,j),xmoy)) {
//						Lpix.insere(i,j); 
						Tpix[ind_insere]=i*nbcol+j; ind_insere++;
						imaLpix(i,j)=1;
						xmoy=(xmoy*n+ima(i,j))/(n+1); n++;
					}
				}
		}
		if (Tpix!=NULL) delete[] Tpix;
		if (iaf) cout<<" nreg "<<nreg<<" moyenne "<<xmoy<<" contient "<<n<<" pixels\n"; //char aa; cin>>aa;
		nreg++;
	}
	nbregions=nreg;
}

imaregions::imaregions (imadata<float> &ima, float s, unsigned short int l0, bool ifus) : imasites(ima.nlig(),ima.ncol()) {
	nblayer=1;
	Tima=new long int[nbpix*nblayer];
	valnul=-1;
	nbregions=0;
	unsigned int jj;
  int i,j,i2,j2,layer=0;
	long int l;
	for (l=0; l<nbpix; l++) {
		for (j=0; j<nblayer; j++) Tima[l*nblayer+j]=valnul;
		Tsites[l]=&(Tima[l*nblayer]);
		}
	s_homog=s;
	cout<<" segmentation par quadtree\n";
	double x=log((double)nblig)/log(2.);
	if (nblig!=nbcol || abs(x-(int)x)>1.e-6) {
		cout<<" #lig = "<<nblig<<" # col = "<<nbcol<<" => pb dans dimensions l'image pour appliquer le quadtree\n";
		}
	else {
		const unsigned short int N=nblig, L=(unsigned short int)(x+1.e-6), t0=(unsigned short int)(N/pow((float)2,l0));
		const unsigned int N2=N*N;
		unsigned short int *taille=new unsigned short int[N2];
		cout<<" N = "<<N<<" N2 = "<<N2<<" L = "<<L<<" t0 = "<<t0<<"\n";
		for (jj=0; jj<N2; jj++) taille[jj]=0;
		const float valmin0=FLT_MAX, valmax0=-FLT_MAX;
		float *Tmin=new float[N2];
		for (jj=0; jj<N2; jj++) Tmin[jj]=valmin0;
		float *Tmax=new float[N2];
		for (jj=0; jj<N2; jj++) Tmax[jj]=valmax0;
		float *Tmoy=new float[N2];
		for (jj=0; jj<N2; jj++) Tmoy[jj]=0.;
		unsigned int *lig0=new unsigned int[N2];
		for (jj=0; jj<N2; jj++) lig0[jj]=nblig;
		unsigned int *col0=new unsigned int[N2];
		for (jj=0; jj<N2; jj++) col0[jj]=nbcol;
// initialisation
		int k, n=(int)pow((float)4,l0), ires=(int)pow((float)4,L-l0);
    unsigned int t, cle, un=(unsigned int)n;
		cout<<" n = "<<n<<" ires = "<<ires<<"\n";
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				i2=i; j2=j; cle=0; t=1;
				for (k=0; k<L; k++) {
					cle+=t*(j2%2); j2=j2/2;
					t*=2;
					cle+=t*(i2%2); i2=i2/2;
					t*=2;
					}
				cle=cle/ires;
				if (cle>=0 && cle<un) {
					float x=ima(i,j,layer);
					if (x<Tmin[cle]) Tmin[cle]=x;
					if (x>Tmax[cle]) Tmax[cle]=x;
					Tmoy[cle]+=x;
					if ((unsigned int)i<lig0[cle]) lig0[cle]=(unsigned int)i;
					if ((unsigned int)j<col0[cle]) col0[cle]=(unsigned int)j;
					Tima[i*nbcol+j]=cle;
					}
				else cout<<" Pb dans calcul de la cle = "<<cle<<" en ("<<i<<","<<j<<")\n";
				}
		for (i=0; i<n; i++) {
			taille[i]=t0;
			Tmoy[i]=Tmoy[i]/t0/t0;
//			cout<<" bloc "<<i<<" taille "<<taille[i]<<" lig,col "<<lig0[i]<<" "<<col0[i];
//			cout<<" Tmin "<<Tmin[i]<<" Tmax "<<Tmax[i]<<" Tmoy "<<Tmoy[i]<<"\n";
			}
		nbregions=n;
		j=l0;
		t=t0;
		k=1;
		cout<<" fin initialisation\n";
// fusion
		while (j>0) {
			int incre=(int)pow((float)4,l0-j+1);
			for (i=0; i<n; i+=incre) {
				if (t==taille[i] && t==taille[i+k] && t==taille[i+2*k] && t==taille[i+3*k]) {
					float xmax=maxi(maxi(maxi(Tmax[i],Tmax[i+k]),Tmax[i+2*k]),Tmax[i+3*k]);
					float xmin=mini(mini(mini(Tmin[i],Tmin[i+k]),Tmin[i+2*k]),Tmin[i+3*k]);
					if ((xmax-xmin)<s_homog) {
//						cout<<xmax-xmin<<" => fusion des blocs "<<i<<" "<<i+k<<" "<<i+2*k<<" "<<i+3*k<<"\n";
						taille[i]=2*t;
						taille[i+k]=taille[i+2*k]=taille[i+3*k]=0;
						Tmax[i]=xmax;
						Tmin[i]=xmin;
						Tmoy[i]=(Tmoy[i]+Tmoy[i+k]+Tmoy[i+2*k]+Tmoy[i+3*k])/4;
						for (i2=0; i2<(int)t; i2++)
							for (j2=0; j2<(int)t; j2++) {
								Tima[(lig0[i+k]+i2)*nbcol+col0[i+k]+j2]=i;
								Tima[(lig0[i+2*k]+i2)*nbcol+col0[i+2*k]+j2]=i;
								Tima[(lig0[i+3*k]+i2)*nbcol+col0[i+3*k]+j2]=i;
								}
						}
//					else
//						cout<<xmax-xmin<<" => NON fusion des blocs "<<i<<" "<<i+k<<" "<<i+2*k<<" "<<i+3*k<<"\n";
					}
				}
			j--;
			t*=2;
			k*=4;
			nbregions-=3;
			}
		cout<<" fin fusion\n";
		j=l0;
// division
		for (i=0; i<n; i++) {
			if (taille[i]>0 && taille[i]<=t0 && taille[i]>1) {
				while (Tmax[i]-Tmin[i]>=s_homog) {
					t=taille[i]/2;
					taille[i]=t;
					taille[n]=taille[n+1]=taille[n+2]=t;
					i2=lig0[i]; j2=col0[i];
					lig0[n]  =i2;   col0[n]  =j2+t;
					lig0[n+1]=i2+t; col0[n+1]=j2;
					lig0[n+2]=i2+t; col0[n+2]=j2+t;
					Tmin[i]  =valmin0; Tmax[i]  =valmax0;
					Tmin[n]  =valmin0; Tmax[n]  =valmax0;
					Tmin[n+1]=valmin0; Tmax[n+1]=valmax0;
					Tmin[n+2]=valmin0; Tmax[n+2]=valmax0;
					for (i2=0; i2<(int)t; i2++)
						for (j2=0; j2<(int)t; j2++) {
							x=ima(lig0[i]+i2,col0[i]+j2,layer);
							if (x<Tmin[i]) Tmin[i]=(float)x;
							if (x>Tmax[i]) Tmax[i]=(float)x;
							x=ima(lig0[n]+i2,col0[n]+j2);
							if (x<Tmin[n]) Tmin[n]=(float)x;
							if (x>Tmax[n]) Tmax[n]=(float)x;
							Tima[(lig0[n]+i2)*nbcol+col0[n]+j2]=n;
							x=ima(lig0[n+1]+i2,col0[n+1]+j2);
							if (x<Tmin[n+1]) Tmin[n+1]=(float)x;
							if (x>Tmax[n+1]) Tmax[n+1]=(float)x;
							Tima[(lig0[n+1]+i2)*nbcol+col0[n+1]+j2]=n+1;
							x=ima(lig0[n+2]+i2,col0[n+2]+j2);
							if (x<Tmin[n+2]) Tmin[n+2]=(float)x;
							if (x>Tmax[n+2]) Tmax[n+2]=(float)x;
							Tima[(lig0[n+2]+i2)*nbcol+col0[n+2]+j2]=n+2;
							}
					n+=3;
					nbregions+=3;
					}
				}
			}
// fin division
		cout<<" # regions final = "<<nbregions<<"\n";
// creation de l'image des regions
		long int * Tireg=new long int[n];
		int ireg=0;
		for (k=0; k<n; k++)
			if (taille[k]>0) Tireg[k]=ireg++;
		for (l=0; l<nbpix; l++) {
			k=Tima[l];
			if (k>=0 && k<n) Tima[l]=Tireg[k];
			}
		nbregions=ireg;
		cout<<" # regions final = "<<nbregions<<"\n";
/*		for (k=0; k<n; k++)
			if (taille[k]>0) cout<<" region "<<k<<" taille "<<taille[k]<<" label "<<Tireg[k]<<"\n";
			else cout<<" region "<<k<<" taille "<<taille[k]<<" pas de label\n";*/
		unsigned short int *taillef=new unsigned short int[nbregions];
		float *Tminf=new float[nbregions];
		float *Tmaxf=new float[nbregions];
		float *Tmoyf=new float[nbregions];
		unsigned int *lig0f=new unsigned int[nbregions];
		unsigned int *col0f=new unsigned int[nbregions];
		for (k=0; k<n; k++) {
			if (taille[k]>0) {
				j=Tireg[k];
				if (j>=0 && j<(int)nbregions) {
					taillef[j]=taille[k]; lig0f[j]=lig0[k]; col0f[j]=col0[k];
					Tminf[j]=Tmin[k]; Tmaxf[j]=Tmax[k]; Tmoyf[j]=Tmoy[k];
					}
				else cout<<" Pb ds etiquette de region "<<k<<" : "<<j<<" <0 ou >="<<nbregions<<"\n";
				}
			}
		if (Tmin!=NULL) delete[] Tmin;
		if (Tmax!=NULL) delete[] Tmax;
		if (Tmoy!=NULL) delete[] Tmoy;
		if (lig0!=NULL) delete[] lig0;
		if (col0!=NULL) delete[] col0;
		if (taille!=NULL) delete[] taille;
		if (Tireg!=NULL) delete[] Tireg;
// fusion des regions
		if (ifus) {
			n=nbregions;
			bool fini=0;
			unsigned int l0,k0,l2,k2,ll0,ll2,kk0,kk2,iter=0;
			while (!fini) {
				fini=1;cout<<" iteration "<<iter++<<" => #reg = "<<nbregions<<"\n";
				for (j=0; j<n; j++) {
					t=taillef[j];
					if (t>0) {
						l0=lig0f[j]; l2=lig0f[j]+t;
						k0=col0f[j]; k2=col0f[j]+t;
						for (jj=j+1; jj<(unsigned int)n; jj++) {
							t=taillef[jj];
							if (t>0 && Tima[l0*nbcol+k0]<Tima[lig0f[jj]*nbcol+col0f[jj]]) {
								ll0=lig0f[jj]; ll2=lig0f[jj]+t;
								kk0=col0f[jj]; kk2=col0f[jj]+t;
								if ((l2==ll0 && ((k0>=kk0 && k2<=kk2) || (kk0>=k0 && kk2<=k2))) ||
									(l0==ll2 && ((k0>=kk0 && k2<=kk2) || (kk0>=k0 && kk2<=k2))) ||
									(k2==kk0 && ((l0>=ll0 && l2<=ll2) || (ll0>=l0 && ll2<=l2))) ||
									(k0==kk2 && ((l0>=ll0 && l2<=ll2) || (ll0>=l0 && ll2<=l2)))) {
									float xmax=maxi(Tmaxf[j],Tmaxf[jj]), xmin=mini(Tminf[j],Tminf[jj]);
									if ((xmax-xmin)<s_homog/1.0) {
										ireg=Tima[l0*nbcol+k0];
										i=Tima[ll0*nbcol+kk0];
//										for (i2=0; i2<t; i2++)
//											for (j2=0; j2<t; j2++) Tima[(ll0+i2)*nbcol+kk0+j2]=i;
										for (i2=0; i2<nblig; i2++)
											for (j2=0; j2<nbcol; j2++)
												if (Tima[i2*nbcol+j2]==i) Tima[i2*nbcol+j2]=ireg;
										nbregions--;
										Tmaxf[j]=xmax;
										Tminf[j]=xmin;
										Tmaxf[jj]=xmax;
										Tminf[jj]=xmin;
//										cout<<" fusion des regions "<<j<<" et "<<jj<<" => i="<<i<<", #reg = "<<nbregions<<"\n";
//										if (nbregions==0) {char aa; cin>>aa;}
										if (fini) fini=0;
										}
									}
								}
							}
						}
					}
				}
			cout<<" # regions final = "<<nbregions<<"\n";
			long int *Tireg=new long int[n];
			for (k=0; k<n; k++) taillef[k]=0;
			for (l=0; l<nbpix; l++) {k=Tima[l]; if (k>=0 && k<n) taillef[k]++;}
			ireg=0;
			for (k=0; k<n; k++)
				if (taillef[k]>0) Tireg[k]=ireg++;
			for (l=0; l<nbpix; l++) {
				k=Tima[l];
				if (k>=0 && k<n) Tima[l]=Tireg[k];
				}
			nbregions=ireg;
			cout<<" # regions final = "<<nbregions<<"\n";
			for (k=0; k<n; k++) taillef[k]=0;
			for (l=0; l<nbpix; l++) {k=Tima[l]; if (k>=0 && k<n) taillef[k]++;}
//			for (k=0; k<n; k++)
//				if (taillef[k]>0) cout<<" region "<<k<<" contient "<<taillef[k]<<" pixels\n";
			if (Tireg!=NULL) delete[] Tireg;
			}
// fin fusion des regions
		if (Tminf!=NULL) delete[] Tminf;
		if (Tmaxf!=NULL) delete[] Tmaxf;
		if (Tmoyf!=NULL) delete[] Tmoyf;
		if (lig0f!=NULL) delete[] lig0f;
		if (col0f!=NULL) delete[] col0f;
		if (taillef!=NULL) delete[] taillef;
		}
	}

imaregions::imaregions (imadata<float> &ima, double _s_homog) : imasites(ima.nlig(),ima.ncol()) {
	nblayer=1;
	Tima=new long int[nbpix*nblayer];
	valnul=/*-1*/0;
	int k=0, i,j;
  unsigned int ui,uj;
  unsigned long jj;
	long int l;
	for (l=0; l<nbpix; l++) {
		for (j=0; j<nblayer; j++) Tima[l*nblayer+j]=valnul;
		Tsites[l]=&(Tima[l*nblayer]);
	}
	s_homog=(float)_s_homog;
	cout<<" fusion de regions dans graphe\n";

// initialisations
	const bool init=1;
	sommet *Tsommets=NULL; arrete *Tarretes=NULL;
	unsigned long int n=0,m=0,nbsommets,nbarretes,nbarretes_max,nv,mmin,mv,nbarrete2,nv2,n1,n2;
	double valmin, valminv, val1, val2;
	const double s_homog2=s_homog*s_homog;
	if (init==1) {
		const int nmax_reg_init=4096;
		const float epsilon=(float)1.e-6, delta_pas_quantif=1.f;
		float pas_quantif=1.f, xx, x_it;
		int nreg_init=nblig*nbcol, ireg, kreg, np;
		bool fini_init=0;
		imadata<float> ima2(ima), imaHomoge(nblig,nbcol); for (i=0; i<nblig; i++) for (j=0; j<nbcol; j++) imaHomoge(i,j)=(float)(i*nbcol+j);
		int *tab_np=new int[nblig*nbcol];
		do {
			x_it=0;
			do {
				l=0; x_it+=1;
				for (i=0; i<nblig; i++) 
					for (j=0; j<nbcol-1; j++) 
						if (fabs(ima2(i,j)-ima2(i,j+1))<=epsilon && imaHomoge(i,j+1)!=imaHomoge(i,j)) {xx=mini(imaHomoge(i,j+1),imaHomoge(i,j)); imaHomoge(i,j+1)=imaHomoge(i,j)=xx; l++;}
				for (i=0; i<nblig-1; i++) 
					for (j=0; j<nbcol; j++) 
						if (fabs(ima2(i,j)-ima2(i+1,j))<=epsilon && imaHomoge(i+1,j)!=imaHomoge(i,j)) {xx=mini(imaHomoge(i+1,j),imaHomoge(i,j)); imaHomoge(i+1,j)=imaHomoge(i,j)=xx; l++;}
				if (((int)x_it)%100==0) cout<<" it "<<(int)x_it<<" l = "<<l<<"\n";
			} while (l>0 && x_it<10000);
			for (l=0; l<nblig*nbcol; l++) tab_np[l]=0;
			for (i=0; i<nblig; i++) for (j=0; j<nbcol; j++) tab_np[(int)imaHomoge(i,j)]++;
			k=0; for (l=0; l<nblig*nbcol; l++) if (tab_np[l]>0) k++;
			cout<<" a l'initialisation "<<k<<" regions de contraste <= "<<pas_quantif<<"\n";
			if (k>nmax_reg_init) {
				pas_quantif+=delta_pas_quantif; if (pas_quantif<1.f) cout<<" attention pas quatification < 1 ???????????\n";
				for (i=0; i<nblig; i++) for (j=0; j<nbcol; j++) {ima2(i,j)=floor(ima(i,j)/pas_quantif);}
			} else {
				fini_init=1; nreg_init=k;
				Tsommets=new sommet[nreg_init];
				n=0;
				int *T_regcree=new int[nblig*nbcol]; for (l=0; l<nblig*nbcol; l++) T_regcree[l]=-1;
				for (i=0; i<nblig; i++) 
					for (j=0; j<nbcol; j++) {
						ireg=(int)imaHomoge(i,j);
						if (T_regcree[ireg]==-1) {
							Tsommets[n].Lpix.insere(i,j); Tsommets[n].val=ima(i,j); Tsommets[n].valid=1; Tsommets[n].n_arretes=0; 
							T_regcree[ireg]=n; n++; //cout<<" creation sommet "<<n-1<<"\n";
						} else {
							kreg=T_regcree[ireg]; np=Tsommets[kreg].Lpix.nb_elts();
							Tsommets[kreg].Lpix.insere(i,j); 
							Tsommets[kreg].val=(np*Tsommets[kreg].val+ima(i,j))/(np+1);
						}
					}
				nbsommets=n; cout<<" "<<nbsommets<<" sommets crees\n";
				m=0; // création des arretes (non orientées) APRES création de toutes les régions
				arrete *Tarretes_bis=new arrete[nreg_init*nreg_init/2];
				bool *T_arretecree=new bool[nreg_init*nreg_init]; for (l=0; l<nreg_init*nreg_init; l++) T_arretecree[l]=false;
				int i0,i2,j0,j2;
				for (i=0; i<nblig; i++) 
					for (j=0; j<nbcol; j++) { //cout<<" i = "<<i<<" j = "<<j<<"\n";
						i0=maxi(0,i-1); i2=mini(i+1,nblig-1); j0=maxi(0,j-1); j2=mini(j+1,nbcol-1);
						ireg=(int)imaHomoge(i,j); if (T_regcree[ireg]<0) {cout<<" region de label ireg = "<<ireg<<" ne correspond pas a 1 region de T_regcree ????????\n"; char aa; cin>>aa;}
						kreg=(int)imaHomoge(i0,j); if (T_regcree[kreg]<0) {cout<<" region de label kreg = "<<kreg<<" ne correspond pas a 1 region de T_regcree ????????\n"; char aa; cin>>aa;}
						if (ireg!=kreg && T_regcree[ireg]>=0 && T_regcree[kreg]>=0 && !T_arretecree[T_regcree[ireg]*nreg_init+T_regcree[kreg]]) {//cout<<" creation arrete i0 "<<m<<"\n";
							T_arretecree[T_regcree[ireg]*nreg_init+T_regcree[kreg]]=true; 
							T_arretecree[T_regcree[kreg]*nreg_init+T_regcree[ireg]]=true; 
							Tarretes_bis[m].pred=T_regcree[ireg]; Tarretes_bis[m].succ=T_regcree[kreg]; 
							Tarretes_bis[m].val=pow(Tsommets[Tarretes_bis[m].succ].val-Tsommets[Tarretes_bis[m].pred].val,2);
							Tarretes_bis[m].valid=1; m++;}
						kreg=(int)imaHomoge(i2,j); if (T_regcree[kreg]<0) {cout<<" region de label kreg = "<<kreg<<" ne correspond pas a 1 region de T_regcree ????????\n"; char aa; cin>>aa;}
						if (ireg!=kreg && T_regcree[ireg]>=0 && T_regcree[kreg]>=0 && !T_arretecree[T_regcree[ireg]*nreg_init+T_regcree[kreg]]) {//cout<<" creation arrete i2 "<<m<<"\n";
							T_arretecree[T_regcree[ireg]*nreg_init+T_regcree[kreg]]=true; 
							T_arretecree[T_regcree[kreg]*nreg_init+T_regcree[ireg]]=true; 
							Tarretes_bis[m].pred=T_regcree[ireg]; Tarretes_bis[m].succ=T_regcree[kreg]; 
							Tarretes_bis[m].val=pow(Tsommets[Tarretes_bis[m].succ].val-Tsommets[Tarretes_bis[m].pred].val,2);
							Tarretes_bis[m].valid=1; m++;}
						kreg=(int)imaHomoge(i,j0); if (T_regcree[kreg]<0) {cout<<" region de label kreg = "<<kreg<<" ne correspond pas a 1 region de T_regcree ????????\n"; char aa; cin>>aa;}
						if (ireg!=kreg && T_regcree[ireg]>=0 && T_regcree[kreg]>=0 && !T_arretecree[T_regcree[ireg]*nreg_init+T_regcree[kreg]]) {//cout<<" creation arrete j0 "<<m<<"\n";
							T_arretecree[T_regcree[kreg]*nreg_init+T_regcree[ireg]]=true; 
							T_arretecree[T_regcree[ireg]*nreg_init+T_regcree[kreg]]=true; 
							Tarretes_bis[m].pred=T_regcree[ireg]; Tarretes_bis[m].succ=T_regcree[kreg]; 
							Tarretes_bis[m].val=pow(Tsommets[Tarretes_bis[m].succ].val-Tsommets[Tarretes_bis[m].pred].val,2);
							Tarretes_bis[m].valid=1; m++;}
						kreg=(int)imaHomoge(i,j2); if (T_regcree[kreg]<0) {cout<<" region de label kreg = "<<kreg<<" ne correspond pas a 1 region de T_regcree ????????\n"; char aa; cin>>aa;}
						if (ireg!=kreg && T_regcree[ireg]>=0 && T_regcree[kreg]>=0 && !T_arretecree[T_regcree[ireg]*nreg_init+T_regcree[kreg]]) {//cout<<" creation arrete j2 "<<m<<"\n";
							T_arretecree[T_regcree[ireg]*nreg_init+T_regcree[kreg]]=true; 
							T_arretecree[T_regcree[kreg]*nreg_init+T_regcree[ireg]]=true; 
							Tarretes_bis[m].pred=T_regcree[ireg]; Tarretes_bis[m].succ=T_regcree[kreg]; 
							Tarretes_bis[m].val=pow(Tsommets[Tarretes_bis[m].succ].val-Tsommets[Tarretes_bis[m].pred].val,2);
							Tarretes_bis[m].valid=1; m++;}
					}
				if (T_regcree!=NULL) delete[] T_regcree; if (T_arretecree!=NULL) delete[] T_arretecree;
				nbarretes=m; cout<<" "<<nbarretes<<" arretes creees\n"; nbarretes_max=nbarretes;
				Tarretes=new arrete[nbarretes_max]; 
				for (m=0; m<nbarretes_max; m++) Tarretes[m]=Tarretes_bis[m];
				if (Tarretes_bis!=NULL) delete[] Tarretes_bis;
			}
			imaHomoge.sauve_Ima("imaHomoge.dat"); ima2.sauve_ImaPGM("ima_quantif.pgm"); //char aa; cin>>aa;
		}	while (!fini_init);
		if (tab_np!=NULL) delete[] tab_np;
	}
	if (init==0) {
		n=m=0;
		Tsommets=new sommet[nbpix]; nbarretes_max=nblig*(nbcol-1)+(nblig-1)*nbcol;
		Tarretes=new arrete[nbarretes_max];    // au max en 4-connexite
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				Tsommets[n].Lpix.insere(i,j); Tsommets[n].val=ima(i,j,k);         // ou autre valeur representative de la region
				Tsommets[n].valid=1; Tsommets[n].n_arretes=0;                           // initialisation juste après
				if (i>0) {
					Tarretes[m].pred=n; Tarretes[m].succ=n-nbcol;
					Tarretes[m].val=pow(Tsommets[Tarretes[m].succ].val-Tsommets[Tarretes[m].pred].val,2);        // ou autre cout
					Tarretes[m].valid=1; m++;
				}
				if (j>0) {
					Tarretes[m].pred=n; Tarretes[m].succ=n-1;
					Tarretes[m].val=pow(Tsommets[Tarretes[m].succ].val-Tsommets[Tarretes[m].pred].val,2);        // ou autre cout
					Tarretes[m].valid=1; m++;
				}
				n++;
			}
		nbsommets=n;
		nbarretes=m;
	}
	const unsigned long int nmaxsommets=nbsommets, nmaxarretes=nbarretes;
	for (m=0; m<nbarretes; m++)
		if (Tarretes[m].valid && Tarretes[m].val>s_homog2) Tarretes[m].valid=0;
	for (m=0; m<nbarretes; m++)
		if (Tarretes[m].valid) {
			n=Tarretes[m].pred;
			Tsommets[n].T_arretes[Tsommets[n].n_arretes]=m; Tsommets[n].n_arretes++;
			n=Tarretes[m].succ;
			Tsommets[n].T_arretes[Tsommets[n].n_arretes]=m; Tsommets[n].n_arretes++;
		}
	k=0;
	for (n=0; n<nbsommets; n++) if ((int)Tsommets[n].n_arretes>k) k=Tsommets[n].n_arretes;
	cout<<" nombre max d'arretes par sommet trouve = "<<k<<"\n";

// fusion des regions
	bool fini=0, iOK;
	unsigned int iter=0;
	while (!fini) {
		cout<<" iteration "<<iter++<<" : # regions = "<<nbsommets<<"\n";
		unsigned long int *Tarrete2=new unsigned long int[nbarretes];
		nbarrete2=0;
		bool *T_isommets=new bool[nmaxsommets];
		for (n=0; n<nmaxsommets; n++) T_isommets[n]=Tsommets[n].valid;
// verification de la consistance du graphe
		for (n=0; n<nmaxsommets; n++)
			if (T_isommets[n])
				for (ui=0; ui<Tsommets[n].n_arretes; ui++) {m=Tsommets[n].T_arretes[ui];
					if (m<0 && m>=nbarretes_max) {
						cout<<" indice d'arrete non valide : sommet "<<n<<" a "<<Tsommets[n].n_arretes<<" arretes; la "<<ui<<"eme arrete a pour numero "<<m<<" <0 ou >="<<nbarretes_max<<"\n"; {char aa; cin>>aa;}}
					else if (Tarretes[m].valid && n!=Tarretes[m].pred && n!=Tarretes[m].succ) {
							cout<<" * Pb : sommet "<<n<<" touche ? arrete "<<m<<" de "<<Tarretes[m].pred<<" vers "<<Tarretes[m].succ<<"\n"; {char aa; cin>>aa;}}
				}
// selection des arretes de moindre cout par accord mutuel
		for (n=0; n<nmaxsommets; n++) {
			if (T_isommets[n]) {
				valmin=s_homog2;
				for (ui=0; ui<Tsommets[n].n_arretes; ui++) {
					m=Tsommets[n].T_arretes[ui];
					if (m>=0 && m<nbarretes_max && Tarretes[m].valid && Tarretes[m].val<valmin) {
							valmin=Tarretes[m].val;
							mmin=m;
					}
				}
				if (valmin<s_homog2) {
					m=mmin;
					if (n==Tarretes[m].pred) nv=Tarretes[m].succ;
					else {
						if (n==Tarretes[m].succ) nv=Tarretes[m].pred;
						else cout<<" * Pb : sommet "<<n<<" touche ? arrete "<<m<<" de "<<Tarretes[m].pred<<" vers "<<Tarretes[m].succ<<"\n";
					}
					if (T_isommets[nv]) {
						valminv=valmin;
						iOK=1;
						for (ui=0; ui<Tsommets[nv].n_arretes; ui++) {
							mv=Tsommets[nv].T_arretes[ui];
							if (m!=mv && Tarretes[mv].valid && Tarretes[mv].val<valminv) iOK=0;
						}
						if (iOK) {
							Tarrete2[nbarrete2++]=mmin;	T_isommets[n]=0; T_isommets[nv]=0;
						}
					}
				}
			}
		}
// fusion des regions associees aux arretes de moindre cout
		if (nbarrete2>0) {
			for (jj=0; jj<nbarrete2; jj++) {
				m=Tarrete2[jj];
				n=Tarretes[m].pred;
				nv=Tarretes[m].succ;
				n1=Tsommets[n].Lpix.nb_elts(); val1=Tsommets[n].val;
				n2=Tsommets[nv].Lpix.nb_elts(); val2=Tsommets[nv].val;
				while (Tsommets[nv].Lpix.nb_elts()>0) {Tsommets[n].Lpix.insere(Tsommets[nv].Lpix.extrait());}
				Tsommets[n].val=(n1*val1+n2*val2)/(n1+n2);                         // ou autre...
				Tarretes[m].valid=0;
				nbarretes--;
				for (ui=0; ui<Tsommets[nv].n_arretes; ui++) {
					mv=Tsommets[nv].T_arretes[ui];
					if (Tarretes[mv].valid) {
						if (nv==Tarretes[mv].pred) nv2=Tarretes[mv].succ;
						else {
							if (nv==Tarretes[mv].succ) nv2=Tarretes[mv].pred;
							else cout<<" + Pb : sommet "<<nv<<" touche ? arrete "<<mv<<" de "<<Tarretes[mv].pred<<" vers "<<Tarretes[mv].succ<<"\n";
						}
						if (nv2!=n) {
							bool nv2OK=1;
							for (uj=0; uj<Tsommets[n].n_arretes; uj++)
								if (Tarretes[Tsommets[n].T_arretes[uj]].pred==nv2 ||	Tarretes[Tsommets[n].T_arretes[uj]].succ==nv2) nv2OK=0;
							if (nv2OK) {
								if (nv==Tarretes[mv].pred) Tarretes[mv].pred=n;
								else Tarretes[mv].succ=n;
//								Tarretes[mv].val=pow(Tsommets[n].val-Tsommets[n2].val,2); // ou autre...
								Tsommets[n].T_arretes[Tsommets[n].n_arretes++]=mv;
								if (Tsommets[n].n_arretes>=nmaxarreteparsommet) {
									cout<<" au sommet "<<n<<", pb de dim de tableau Tsommets[n].T_arretes[.] -> elimination des 'trous'\n";
									cout<<" on va passer de "<<Tsommets[n].n_arretes<<" arretes a "; 
									int ii, jj;
									for (ii=0; ii<(int)Tsommets[n].n_arretes; ii++) {
										if (!Tarretes[Tsommets[n].T_arretes[ii]].valid) {cout<<" 1 trou d'elimine...\n";
											for (jj=ii; jj<(int)Tsommets[n].n_arretes-1; jj++) Tsommets[n].T_arretes[jj]=Tsommets[n].T_arretes[jj+1];
											Tsommets[n].n_arretes--;
										}
									} cout<<Tsommets[n].n_arretes<<" arretes\n";
									if (Tsommets[n].n_arretes>=nmaxarreteparsommet) {cout<<"Pb de dim de tableau Tsommets[n].T_arretes[.] persiste !!!\n"; char aa; cin>>aa;}
								}
							}
							else Tarretes[mv].valid=0;
						}
					}
				}
				Tsommets[nv].valid=0;
				nbsommets--;
			}
			for (m=0; m<nmaxarretes; m++)
				if (Tarretes[m].valid) {
					Tarretes[m].val=pow(Tsommets[Tarretes[m].succ].val-Tsommets[Tarretes[m].pred].val,2);      // ou autre cout
					if (Tarretes[m].val>s_homog2) Tarretes[m].valid=0;
				}
		}
		else fini=1;
		if (Tarrete2!=NULL) delete[] Tarrete2;
		if (T_isommets!=NULL) delete[] T_isommets;
	}
// Creation de l'image des regions
	int ireg=valnul;
	for (n=0; n<nmaxsommets; n++) {
		if (Tsommets[n].Lpix.nb_elts()>0) ireg++;
		while (Tsommets[n].Lpix.nb_elts()>0) {
			elt_liste P=Tsommets[n].Lpix.extrait();
			if (P.x>=0 && P.x<nblig && P.y>=0 && P.y<nbcol) Tima[P.x*nbcol+P.y]=ireg; // couche 0
		}
	}
	nbregions=ireg-valnul;
	if (Tsommets!=NULL) delete[] Tsommets;
	if (Tarretes!=NULL) delete[] Tarretes;
	cout<<" fin fusion de regions dans graphe\n";
}

imaregions::imaregions (imadata<float> &ima, unsigned int nreg) : imasites(ima.nlig(),ima.ncol()) {
	nblayer=1; valnul=/*-1*/0;
	Tima=new long int[nbpix*nblayer];
	int k=0, i,j;
	unsigned int ui,uj;
	unsigned long jj;
	for (i=0; i<nbpix; i++) {
		for (j=0; j<nblayer; j++) Tima[i*nblayer+j]=valnul;
		Tsites[i]=&(Tima[i*nblayer]);
  }
	s_homog=-1;
	cout<<" segmentation par graphe\n";

// initialisations
	const bool init=1, i_reg1CC=0;
	sommet *Tsommets=NULL; arrete *Tarretes=NULL;
	unsigned long int n=0,m=0,nbsommets,nbarretes,nbarretes_max,nv,mmin,mv,nbarrete2,nv2,n1,n2;
	double valmin, valminv, val1, val2;
	const double valmin0=FLT_MAX;
	imabin ima1reg(nblig,nbcol);
	liste_pixels Lpix2; elt_liste P;
	int ncc, kreg, np;
	if (init==1) {
		const int nmax_reg_init=4096;
		const float delta_pas_quantif=1.f;
		int nreg_init, ireg;
		imadata<float> imaHomoge=ima.init1_segmentationregions(nmax_reg_init,delta_pas_quantif,nreg_init);
		Tsommets=new sommet[nreg_init];
		n=0;
		int *T_regcree=new int[nblig*nbcol]; for (i=0; i<nblig*nbcol; i++) T_regcree[i]=-1;
		for (i=0; i<nblig; i++) 
			for (j=0; j<nbcol; j++) {
				ireg=(int)imaHomoge(i,j);
				if (T_regcree[ireg]==-1) {
					Tsommets[n].Lpix.insere(i,j); Tsommets[n].val=ima(i,j); Tsommets[n].valid=1; Tsommets[n].n_arretes=0; 
					T_regcree[ireg]=n; n++; //cout<<" creation sommet "<<n-1<<"\n";
				} else {
					kreg=T_regcree[ireg]; np=Tsommets[kreg].Lpix.nb_elts();
					Tsommets[kreg].Lpix.insere(i,j); 
					Tsommets[kreg].val=(np*Tsommets[kreg].val+ima(i,j))/(np+1);
				}
			}
		nbsommets=n; cout<<" "<<nbsommets<<" sommets crees\n";
// verification de la consistance des regions
		if (i_reg1CC) {
			for (n=0; n<nbsommets; n++)
				if (Tsommets[n].valid) {if (n%100==0) cout<<" region "<<n<<" verifiee\n";
					ima1reg.mise_a_zero(); Lpix2=Tsommets[n].Lpix;
					while (Lpix2.nb_elts()>0) {
						P=Lpix2.extrait();
						if (P.x>=0 && P.x<nblig && P.y>=0 && P.y<nbcol) ima1reg(P.x,P.y)=1;
					}
					ima1reg.composantes_connexes(ncc,8,0);
					if (ncc!=1) {cout<<" a l'initialisation, PB : region "<<n<<" a plusieurs CC : ncc = "<<ncc<<"\n"; 
						ima1reg.imaunsignedchar().sauve_ImaPGM("ima1reg.pgm"); char aa; cin>>aa;}
//					else cout<<" a l'initialisation, region "<<n<<" Ok : ncc = "<<ncc<<"\n";
				}
		}
// création des arretes (non orientées) APRES création de toutes les régions
		m=0; // création des arretes (non orientées) APRES création de toutes les régions
		arrete *Tarretes_bis=new arrete[nreg_init*nreg_init/2];
		bool *T_arretecree=new bool[nreg_init*nreg_init]; for (i=0; i<nreg_init*nreg_init; i++) T_arretecree[i]=false;
		int i0,i2,j0,j2;
		for (i=0; i<nblig; i++) 
			for (j=0; j<nbcol; j++) { //cout<<" i = "<<i<<" j = "<<j<<"\n";
				i0=maxi(0,i-1); i2=mini(i+1,nblig-1); j0=maxi(0,j-1); j2=mini(j+1,nbcol-1);
				ireg=(int)imaHomoge(i,j); if (T_regcree[ireg]<0) {cout<<" region de label ireg = "<<ireg<<" ne correspond pas a 1 region de T_regcree ????????\n"; char aa; cin>>aa;}
				kreg=(int)imaHomoge(i0,j); if (T_regcree[kreg]<0) {cout<<" region de label kreg = "<<kreg<<" ne correspond pas a 1 region de T_regcree ????????\n"; char aa; cin>>aa;}
				if (ireg!=kreg && T_regcree[ireg]>=0 && T_regcree[kreg]>=0 && !T_arretecree[T_regcree[ireg]*nreg_init+T_regcree[kreg]]) {//cout<<" creation arrete i0 "<<m<<"\n";
					T_arretecree[T_regcree[ireg]*nreg_init+T_regcree[kreg]]=true; 
					T_arretecree[T_regcree[kreg]*nreg_init+T_regcree[ireg]]=true; 
					Tarretes_bis[m].pred=T_regcree[ireg]; Tarretes_bis[m].succ=T_regcree[kreg]; 
					Tarretes_bis[m].val=pow(Tsommets[Tarretes_bis[m].succ].val-Tsommets[Tarretes_bis[m].pred].val,2);
					Tarretes_bis[m].valid=1; m++;}
				kreg=(int)imaHomoge(i2,j); if (T_regcree[kreg]<0) {cout<<" region de label kreg = "<<kreg<<" ne correspond pas a 1 region de T_regcree ????????\n"; char aa; cin>>aa;}
				if (ireg!=kreg && T_regcree[ireg]>=0 && T_regcree[kreg]>=0 && !T_arretecree[T_regcree[ireg]*nreg_init+T_regcree[kreg]]) {//cout<<" creation arrete i2 "<<m<<"\n";
					T_arretecree[T_regcree[ireg]*nreg_init+T_regcree[kreg]]=true; 
					T_arretecree[T_regcree[kreg]*nreg_init+T_regcree[ireg]]=true; 
					Tarretes_bis[m].pred=T_regcree[ireg]; Tarretes_bis[m].succ=T_regcree[kreg]; 
					Tarretes_bis[m].val=pow(Tsommets[Tarretes_bis[m].succ].val-Tsommets[Tarretes_bis[m].pred].val,2);
					Tarretes_bis[m].valid=1; m++;}
				kreg=(int)imaHomoge(i,j0); if (T_regcree[kreg]<0) {cout<<" region de label kreg = "<<kreg<<" ne correspond pas a 1 region de T_regcree ????????\n"; char aa; cin>>aa;}
				if (ireg!=kreg && T_regcree[ireg]>=0 && T_regcree[kreg]>=0 && !T_arretecree[T_regcree[ireg]*nreg_init+T_regcree[kreg]]) {//cout<<" creation arrete j0 "<<m<<"\n";
					T_arretecree[T_regcree[kreg]*nreg_init+T_regcree[ireg]]=true; 
					T_arretecree[T_regcree[ireg]*nreg_init+T_regcree[kreg]]=true; 
					Tarretes_bis[m].pred=T_regcree[ireg]; Tarretes_bis[m].succ=T_regcree[kreg]; 
					Tarretes_bis[m].val=pow(Tsommets[Tarretes_bis[m].succ].val-Tsommets[Tarretes_bis[m].pred].val,2);
					Tarretes_bis[m].valid=1; m++;}
				kreg=(int)imaHomoge(i,j2); if (T_regcree[kreg]<0) {cout<<" region de label kreg = "<<kreg<<" ne correspond pas a 1 region de T_regcree ????????\n"; char aa; cin>>aa;}
				if (ireg!=kreg && T_regcree[ireg]>=0 && T_regcree[kreg]>=0 && !T_arretecree[T_regcree[ireg]*nreg_init+T_regcree[kreg]]) {//cout<<" creation arrete j2 "<<m<<"\n";
					T_arretecree[T_regcree[ireg]*nreg_init+T_regcree[kreg]]=true; 
					T_arretecree[T_regcree[kreg]*nreg_init+T_regcree[ireg]]=true; 
					Tarretes_bis[m].pred=T_regcree[ireg]; Tarretes_bis[m].succ=T_regcree[kreg]; 
					Tarretes_bis[m].val=pow(Tsommets[Tarretes_bis[m].succ].val-Tsommets[Tarretes_bis[m].pred].val,2);
					Tarretes_bis[m].valid=1; m++;}
			}
		if (T_regcree!=NULL) delete[] T_regcree; if (T_arretecree!=NULL) delete[] T_arretecree;
		nbarretes=m; cout<<" "<<nbarretes<<" arretes creees\n"; nbarretes_max=nbarretes;
		Tarretes=new arrete[nbarretes_max]; 
		for (m=0; m<nbarretes_max; m++) Tarretes[m]=Tarretes_bis[m];
		if (Tarretes_bis!=NULL) delete[] Tarretes_bis;
	}
	if (init==0) {
		n=m=0;
		Tsommets=new sommet[nbpix]; nbarretes_max=nblig*(nbcol-1)+(nblig-1)*nbcol;
		Tarretes=new arrete[nbarretes_max];    // au max en 4-connexite
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				Tsommets[n].Lpix.insere(i,j); Tsommets[n].val=ima(i,j,k);         // ou autre valeur representative de la region
				Tsommets[n].valid=1; Tsommets[n].n_arretes=0;                           // initialisation juste après
				if (i>0) {
					Tarretes[m].pred=n; Tarretes[m].succ=n-nbcol;
					Tarretes[m].val=pow(Tsommets[Tarretes[m].succ].val-Tsommets[Tarretes[m].pred].val,2);        // ou autre cout
					Tarretes[m].valid=1; m++;
				}
				if (j>0) {
					Tarretes[m].pred=n; Tarretes[m].succ=n-1;
					Tarretes[m].val=pow(Tsommets[Tarretes[m].succ].val-Tsommets[Tarretes[m].pred].val,2);        // ou autre cout
					Tarretes[m].valid=1; m++;
				}
				n++;
			}
		nbsommets=n;
		nbarretes=m;
	}
// debut de la segmentation par graphe
	const unsigned long int nmaxsommets=nbsommets, nmaxarretes=nbarretes;
	for (m=0; m<nbarretes; m++)
		if (Tarretes[m].valid) {
			n=Tarretes[m].pred;
			Tsommets[n].T_arretes[Tsommets[n].n_arretes]=m; Tsommets[n].n_arretes++;
			n=Tarretes[m].succ;
			Tsommets[n].T_arretes[Tsommets[n].n_arretes]=m; Tsommets[n].n_arretes++;
		}
	k=0;
	for (n=0; n<nbsommets; n++) if ((int)Tsommets[n].n_arretes>k) k=Tsommets[n].n_arretes;
	cout<<" nombre max d'arretes par sommet trouve = "<<k<<"\n";
// fusion des regions
	bool fini=0, iOK;
	unsigned int iter=0;
	while (!fini) {
		cout<<" iteration "<<iter++<<" : # regions = "<<nbsommets<<"\n";
		unsigned long int *Tarrete2=new unsigned long int[nbarretes];
		nbarrete2=0;
		bool *T_isommets=new bool[nmaxsommets];
		for (n=0; n<nmaxsommets; n++) T_isommets[n]=Tsommets[n].valid;
// verification de la consistance des regions
		if (i_reg1CC) {
			for (n=0; n<nmaxsommets; n++)
				if (T_isommets[n]) {
					ima1reg.mise_a_zero(); Lpix2=Tsommets[n].Lpix;
					while (Lpix2.nb_elts()>0) {
						P=Lpix2.extrait();
						if (P.x>=0 && P.x<nblig && P.y>=0 && P.y<nbcol) ima1reg(P.x,P.y)=1;
					}
					ima1reg.composantes_connexes(ncc,8,0);
					if (ncc!=1) {cout<<" PB : region "<<n<<" a plusieurs CC : ncc = "<<ncc<<"\n"; 
						ima1reg.imaunsignedchar().sauve_ImaPGM("ima1reg.pgm"); char aa; cin>>aa;}
				}
		}
// verification de la consistance du graphe
		for (n=0; n<nmaxsommets; n++)
			if (T_isommets[n])
				for (ui=0; ui<Tsommets[n].n_arretes; ui++) {m=Tsommets[n].T_arretes[ui];
					if (m<0 && m>=nbarretes_max) {
						cout<<" indice d'arrete non valide : sommet "<<n<<" a "<<Tsommets[n].n_arretes<<" arretes; la "<<ui<<"eme arrete a pour numero "<<m<<" <0 ou >="<<nbarretes_max<<"\n"; {char aa; cin>>aa;}}
					else if (Tarretes[m].valid && n!=Tarretes[m].pred && n!=Tarretes[m].succ) {
							cout<<" * Pb : sommet "<<n<<" touche ? arrete "<<m<<" de "<<Tarretes[m].pred<<" vers "<<Tarretes[m].succ<<"\n"; {char aa; cin>>aa;}}
				}
// selection des arretes de moindre cout par accord mutuel
		for (n=0; n<nmaxsommets; n++) {
			if (T_isommets[n]) {
				valmin=valmin0;
				for (ui=0; ui<Tsommets[n].n_arretes; ui++) {
					m=Tsommets[n].T_arretes[ui];
					if (m>=0 && m<nbarretes_max && Tarretes[m].valid && Tarretes[m].val<valmin) {
							valmin=Tarretes[m].val;
							mmin=m;
					}
				}
				if (valmin<valmin0) {
					m=mmin;
					if (n==Tarretes[m].pred) nv=Tarretes[m].succ;
					else {
						if (n==Tarretes[m].succ) nv=Tarretes[m].pred;
						else cout<<" * Pb : sommet "<<n<<" touche ? arrete "<<m<<" de "<<Tarretes[m].pred<<" vers "<<Tarretes[m].succ<<"\n";
					}
					if (T_isommets[nv]) {
						valminv=valmin;
						iOK=1;
						for (ui=0; ui<Tsommets[nv].n_arretes; ui++) {
							mv=Tsommets[nv].T_arretes[ui];
							if (m!=mv && Tarretes[mv].valid && Tarretes[mv].val<valminv) iOK=0;
						}
						if (iOK) {
							Tarrete2[nbarrete2++]=mmin;	T_isommets[n]=0; T_isommets[nv]=0;
						}
					}
				}
			}
		}
// fusion des regions associees aux arretes de moindre cout
		if (nbarrete2>0) {
			for (jj=0; jj<nbarrete2; jj++) {
				m=Tarrete2[jj];
				n=Tarretes[m].pred;
				nv=Tarretes[m].succ;
				n1=Tsommets[n].Lpix.nb_elts(); val1=Tsommets[n].val;
				n2=Tsommets[nv].Lpix.nb_elts(); val2=Tsommets[nv].val;
				while (Tsommets[nv].Lpix.nb_elts()>0) {Tsommets[n].Lpix.insere(Tsommets[nv].Lpix.extrait());}
				Tsommets[n].val=(n1*val1+n2*val2)/(n1+n2);                         // ou autre...
				Tarretes[m].valid=0;
				nbarretes--;
				for (ui=0; ui<Tsommets[nv].n_arretes; ui++) {
					mv=Tsommets[nv].T_arretes[ui];
					if (Tarretes[mv].valid) {
						if (nv==Tarretes[mv].pred) nv2=Tarretes[mv].succ;
						else {
							if (nv==Tarretes[mv].succ) nv2=Tarretes[mv].pred;
							else cout<<" + Pb : sommet "<<nv<<" touche ? arrete "<<mv<<" de "<<Tarretes[mv].pred<<" vers "<<Tarretes[mv].succ<<"\n";
						}
						if (nv2!=n) {
							bool nv2OK=1;
							for (uj=0; uj<Tsommets[n].n_arretes; uj++)
								if (Tarretes[Tsommets[n].T_arretes[uj]].pred==nv2 || Tarretes[Tsommets[n].T_arretes[uj]].succ==nv2) nv2OK=0;
							if (nv2OK) {
								if (nv==Tarretes[mv].pred) Tarretes[mv].pred=n;
								else Tarretes[mv].succ=n;
//								Tarretes[mv].val=pow(Tsommets[n].val-Tsommets[n2].val,2); // ou autre...
								Tsommets[n].T_arretes[Tsommets[n].n_arretes++]=mv;
								if (Tsommets[n].n_arretes>=nmaxarreteparsommet) {
									cout<<" au sommet "<<n<<", pb de dim de tableau Tsommets[n].T_arretes[.] -> elimination des 'trous'\n";
									cout<<" on va passer de "<<Tsommets[n].n_arretes<<" arretes a "; 
									int ii, jj;
									for (ii=0; ii<(int)Tsommets[n].n_arretes; ii++) {
										if (!Tarretes[Tsommets[n].T_arretes[ii]].valid) {cout<<" 1 trou d'elimine...\n";
											for (jj=ii; jj<(int)Tsommets[n].n_arretes-1; jj++) Tsommets[n].T_arretes[jj]=Tsommets[n].T_arretes[jj+1];
											Tsommets[n].n_arretes--;
										}
									} cout<<Tsommets[n].n_arretes<<" arretes\n";
									if (Tsommets[n].n_arretes>=nmaxarreteparsommet) {cout<<"Pb de dim de tableau Tsommets[n].T_arretes[.] persiste !!!\n"; char aa; cin>>aa;}
								}
							}
							else Tarretes[mv].valid=0;
						}
//						else {
//							for (uj=0; uj<Tsommets[n].n_arretes; uj++) {
//								if (Tarretes[Tsommets[n].T_arretes[uj]].pred==nv ||
//									Tarretes[Tsommets[n].T_arretes[uj]].succ==nv) {
//										for (k=uj; k<Tsommets[n].n_arretes-1; k++)
//											Tsommets[n].T_arretes[k]=Tsommets[n].T_arretes[k+1];
//										Tsommets[n].n_arretes--;
//										uj--;
//									}
//								}
//							nbarretes--;
//							}
					}
				}
				Tsommets[nv].valid=0;
				nbsommets--;
			}
			for (m=0; m<nmaxarretes; m++)
				if (Tarretes[m].valid)
					Tarretes[m].val=pow(Tsommets[Tarretes[m].succ].val-Tsommets[Tarretes[m].pred].val,2);      // ou autre cout
// Sauvegarde intermediaire de l'image des regions
			if ((nbsommets>1000 && nbsommets%500==0) || (nbsommets<=1000 && nbsommets%100==0) || (nbsommets<=100 && nbsommets%10==0) || (nbsommets<=20)) {
				char nomficI[80], aa[5]; strcpy(nomficI,"imasauv0000.dat"); strcpy(aa,"0000"); 
				int k=nbsommets; for (int i=4; i>0; i--) {aa[i-1]='0'+k%10; k=k/10;}
				strncpy(nomficI+7,aa,4);
				imadata<int> imasauv(nblig,nbcol); imasauv.mise_a_zero();
				int ireg=valnul; 
				for (n=0; n<nmaxsommets; n++) {
					Lpix2=Tsommets[n].Lpix;
					if (Lpix2.nb_elts()>0) {
						ireg++;
						while (Lpix2.nb_elts()>0) {
							P=Lpix2.extrait();
							if (P.x>=0 && P.x<nblig && P.y>=0 && P.y<nbcol) imasauv(P.x,P.y)=ireg;
						}
					}
				}
				imasauv.sauve_ImaBSQ(nomficI);
			}
		}
		else fini=1;
		if (nbsommets<=nreg) fini=1;
		if (Tarrete2!=NULL) delete[] Tarrete2;
		if (T_isommets!=NULL) delete[] T_isommets;
	}
// Creation de l'image des regions
	int ireg=valnul;
	for (n=0; n<nmaxsommets; n++) {
		if (Tsommets[n].Lpix.nb_elts()>0) ireg++;
		while (Tsommets[n].Lpix.nb_elts()>0) {
			elt_liste P=Tsommets[n].Lpix.extrait();
			if (P.x>=0 && P.x<nblig && P.y>=0 && P.y<nbcol) Tima[P.x*nbcol+P.y]=ireg; // couche 0
		}
	}
	nbregions=ireg-valnul;
	if (Tsommets!=NULL) delete[] Tsommets;
	if (Tarretes!=NULL) delete[] Tarretes;
	cout<<" fin segmentation par graphe\n";
}

unsigned long int imaregions::comptearrete (imadata<int> &imaRg, unsigned long int nbsommets) {
	imabin Iarrete(nbsommets,nbsommets); 
	Iarrete.mise_a_zero();
	int nblig=imaRg.nlig(), nbcol=imaRg.ncol(), i,j,k,n,i0,i2,j0,j2,ii,jj;
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) {
			n=imaRg(i,j);
			if (n>=0 && n<(int)nbsommets) {
				i0=maxi(i-1,0); i2=mini(i+1,nblig-1); j0=maxi(j-1,0); j2=mini(j+1,nbcol-1);
				for (ii=i0; ii<=i2; ii++)
					for (jj=j0; jj<=j2; jj++) {
						k=imaRg(ii,jj);
						if ((ii==i || jj==j) && k>=0 && k<(int)nbsommets) Iarrete(k,n)=Iarrete(n,k)=1; // cas 4-connexite
//						if (k>=0 && k<nreginit) Iarrete(k,n)=Iarrete(n,k)=1;										 // cas 8-connexite
					}
			}
		}
	return (unsigned long int)(Iarrete.norm()/2);
}

void imaregions::cree_graph_init (imadata<int> &imaRg, imadata<float> &ima, regionV* Tregions, arrete* Tarretes, 
																	unsigned long int &nbarretes, unsigned long int &nbsommets) {
	cout<<" entree dans cree_graph_init\n";
	int nblig=imaRg.nlig(), nbcol=imaRg.ncol(), nreginit=0, narrinit=0,i,j,k,m,n,i0,i2,j0,j2,ii,jj, n1,n2;
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) if (imaRg(i,j)>nreginit) nreginit=imaRg(i,j);
	nreginit++;
	cout<<" initialement "<<nreginit<<" regions\n"; 
	if (nreginit!=nbsommets) {
		cout<<" attention "<<nreginit<<" <> "<<nbsommets<<"\n";
		if (nreginit>(int)nbsommets) {
			cout<<" pb Tableau de regions trop petit !!!\n"; cout<<" ancienne @ Tregions"<<Tregions<<"\n";
			if (Tregions!=NULL) delete[] Tregions;
			Tregions=new regionV[nreginit]; cout<<" nouvelle @ Tregions"<<Tregions<<"\n"; char aa; cin>>aa;
		}
		nbsommets=nreginit;
	}
	unsigned int **Treg_vois=new unsigned int*[nreginit], *T_npix=new unsigned int[nreginit];
	float **T_valr=new float*[nreginit];
	for (n=0; n<nreginit; n++) {Treg_vois[n]=new unsigned int[nreginit]; T_valr[n]=new float[ima.ncanaux()];}
	for (n=0; n<nreginit; n++) {
		T_npix[n]=0; for (k=0; k<ima.ncanaux(); k++) T_valr[n][k]=0.f;
		for (j=0; j<nreginit; j++) Treg_vois[n][j]=0;
	}
	double lg_a, xx;
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) {
			n=imaRg(i,j);
			if (n>=0 && n<nreginit) {
				Tregions[n].Lpix.insere(i,j);
				T_npix[n]++; for (k=0; k<ima.ncanaux(); k++) T_valr[n][k]+=ima(i,j,k);	
				i0=maxi(i-1,0); i2=mini(i+1,nblig-1); j0=maxi(j-1,0); j2=mini(j+1,nbcol-1);
				for (ii=i0; ii<=i2; ii++)
					for (jj=j0; jj<=j2; jj++) {
						k=imaRg(ii,jj);
						if ((ii==i || jj==j) && k>=0 && k<nreginit) Treg_vois[n][k]++; // cas 4-connexite
//						if (k>=0 && k<nreginit) Treg_vois[n][k]++;										 // cas 8-connexite
					}
			}
		}
	narrinit=0;
	for (n=0; n<nreginit; n++) {
		if (T_npix[n]>0) {
			Tregions[n].val=new double[ima.ncanaux()];// allocation du tableau des valeurs des canaux sur la region
			for (k=0; k<ima.ncanaux(); k++) Tregions[n].val[k]=T_valr[n][k]/T_npix[n];  // sa valeur d'un canal est la moyenne des pixels de la CC
			Tregions[n].valid=1;                     // elle est valide
			Tregions[n].n_arretes=0;                 // aucune arete n'a ete definie a ce stade
			m=0;
			for (j=0; j<nreginit; j++) 
				if (j!=n && Treg_vois[n][j]>0) m++;
			Tregions[n].T_arretes=new unsigned int[m];   // allocation du tableau des arrêtes
			narrinit+=m;
//			Tregions[n].n_arretes=m;
		} else Tregions[n].valid=0;
	}
	narrinit/=2;
	if (narrinit!=nbarretes) {
		cout<<" attention "<<narrinit<<" <> "<<nbarretes<<"\n";
		if (narrinit>(int)nbarretes) {
			cout<<" pb Tableau des arretes trop petit !!!\n"; cout<<" ancienne @ Tarretes"<<Tarretes<<"\n";
			if (Tarretes!=NULL) delete[] Tarretes;
			Tarretes=new arrete[narrinit]; cout<<" nouvelle @ Tarretes"<<Tarretes<<"\n"; char aa; cin>>aa;
		}
		nbarretes=narrinit;
	}
	cout<<" le graphe des regions a "<<nbsommets<<" sommets et "<<nbarretes<<" arretes en tout\n"; 
	m=0;
	for (n=0; n<nreginit; n++)
		for (j=n+1; j<nreginit; j++) {
			lg_a=Treg_vois[n][j];
			if (lg_a>0) {
				Tarretes[m].length=lg_a;		          // longueur initiale de l'arrete
				Tarretes[m].pred=n;			              // region de depart de l'arrete
				Tarretes[m].succ=j;	                  // region d'arrivee de l'arrete
				n1=Tregions[n].Lpix.nb_elts(); n2=Tregions[j].Lpix.nb_elts();
				xx=0.; for (k=0; k<ima.ncanaux(); k++) xx+=pow(Tregions[n].val[k]-Tregions[j].val[k],2);
				Tarretes[m].val=xx*n1*n2/(n1+n2)/lg_a; // valeur de l'arrete
				Tarretes[m].valid=1;
				m++;
			}
		}
// boucle calculant le nombre d'arretes par region
	for (m=0; m<narrinit; m++) { // boucle sur le # max d'arretes car les arretes valides peuvent etre n'importe ou dans le tableau
		if (Tarretes[m].valid) {							  // si l'arrete existe, pour chaque sommet
			n1=Tarretes[m].pred;							  // on determine le numero du sommet
			Tregions[n1].T_arretes[Tregions[n1].n_arretes]=m; // l'arrete est associee au sommet
			Tregions[n1].n_arretes++;						  // le # d'arretes du sommet est incremente
			n2=Tarretes[m].succ;
			Tregions[n2].T_arretes[Tregions[n2].n_arretes]=m;
			Tregions[n2].n_arretes++;
		}
	}
	if (T_npix!=NULL) delete[] T_npix;
	for (n=0; n<nreginit; n++) {
		if (Treg_vois[n]!=NULL) delete[] Treg_vois[n]; if (T_valr[n]!=NULL) delete[] T_valr[n];}
	if (Treg_vois!=NULL) delete[] Treg_vois;
	if (T_valr!=NULL) delete[] T_valr;
}

imaregions::imaregions (imadata<float> &ima, imadata<int> &imaRg, unsigned int nreg, char *nomfics) : imasites(ima.nlig(),ima.ncol()) { // segmentation Koepfler (fonctionnelle Mumford & Sha simplifiee)
	nblayer=1; valnul=/*-1*/0;
	Tima=new long int[nbpix*nblayer];
	int i,j,k=0,np,dim=ima.ncanaux();
  unsigned int ui,uj;
  unsigned long jj;
	for (i=0; i<nbpix; i++) {
		for (j=0; j<nblayer; j++) Tima[i*nblayer+j]=valnul;
		Tsites[i]=&(Tima[i*nblayer]);
	}
	unsigned long int n=0,m=0,nbsommets=(unsigned long int)imaRg.maxI()+1,nbarretes,nv,mv,nbarrete2,nv2,n1,n2,n3,n4,m1,m2; cout<<" nbsommets = "<<nbsommets<<"\n";
	regionV *Tregions=new regionV[nbsommets]; cout<<" @ Tregions "<<Tregions<<"\n";
	nbarretes=comptearrete (imaRg,nbsommets); cout<<" nbarretes = "<<nbarretes<<"\n";
	arrete *Tarretes=new arrete[nbarretes]; cout<<" @ Tarretes "<<Tarretes<<"\n";
	cree_graph_init (imaRg,ima,Tregions,Tarretes,nbarretes,nbsommets); cout<<" @ Tarretes "<<Tarretes<<"\n"; cout<<" @ Tregions "<<Tregions<<"\n";
//	for (m=0; m<nbsommets; m++) cout<<" m "<<m<<" Tregions[m].valid "<<(int)Tregions[m].valid<<"\n";
//	for (m=0; m<nbarretes; m++) cout<<" m "<<m<<" Tarretes[m].valid "<<(int)Tarretes[m].valid<<" Tarretes[m].val "<<(int)Tarretes[m].val<<"\n";

	cout<<" segmentation selon la fonctionnelle simplifiee de Mumford et Shah\n";
	const bool i_reg1CC=0;
	liste_pixels Lpix2; elt_liste P;

// Variables et initialisation
	double *val1=new double[dim], *val2=new double[dim]; 
	const double valmin0=FLT_MAX;
	double lambda=valmin0, xx, eps=/*0.*/1.e-2;
	const unsigned long int nmaxsommets=nbsommets, nmaxarretes=nbarretes;
	cout<<" nmaxsommets "<<nmaxsommets<<" nmaxarretes "<<nmaxarretes<<"\n";
	unsigned long int *Tarrete2=new unsigned long int[nmaxarretes]; // tableau des arretes a supprimer (taille = # courant d'arretes)
	cout<<" + Tarrete2 cree en @ "<<Tarrete2<<"\n"; //char aa; cin>>aa;
// boucle calculant le nombre d'arretes par region
	for (m=0; m<nmaxarretes; m++) { // boucle sur le # max d'arretes car les arretes valides peuvent etre n'importe ou dans le tableau
		if (Tarretes[m].valid) {
			xx=Tarretes[m].val;
			if (xx<lambda) {lambda=xx; cout<<" fusion regions "<<Tarretes[m].pred<<" & "<<Tarretes[m].succ<<" si lambda>="<<lambda<<"\n";}
		}
	}
	lambda+=eps;
	cout<<" valeur du prochain parametre d'echelle lambda = "<<lambda<<"\n";
//	for (m=0; m<nbarretes; m++) cout<<" +++ m "<<m<<" Tarretes[m].valid "<<(int)Tarretes[m].valid<<" Tarretes[m].val "<<(int)Tarretes[m].val<<"\n";
// iterations sur la fusion des regions
	bool fini=0;
	unsigned int iter=0;
	int iter_aff=2;
	while (!fini) { //cout<<iter<<" "<<iter_aff<<" "<<iter%iter_aff<<" ";
//		if ((iter%iter_aff)==0)
		cout<<" iteration "<<iter++<<" : # regions = "<<nbsommets<<"\n";
//		unsigned long int *Tarrete2=new unsigned long int[nbarretes]; // tableau des arretes a supprimer (taille = # courant d'arretes)
//		if (Tarrete2!NULL) delete[] Tarrete2; Tarrete2=new unsigned long int[nbarretes]; // tableau des arretes a supprimer (taille = # courant d'arretes)
//		cout<<" Tarrete2 cree en @ "<<Tarrete2<<"\n"; char aa; cin>>aa;
		nbarrete2=0;
// determination des arretes a supprimer
		for (m=0; m<nmaxarretes; m++) {//cout<<m<<" "<<(int)Tarretes[m].valid<<" "<<lambda<<" >= ? "<<Tarretes[m].val;
			if (Tarretes[m].valid && Tarretes[m].val<=lambda) {Tarrete2[nbarrete2++]=m; 
			cout<<" $$$$$ "<<m<<" "<<(int)Tarretes[m].valid<<" "<<lambda<<" >= "<<Tarretes[m].val<<"\n";}}
// fusion des regions
		if (nbarrete2>0) {                                // si il y a des fusions de regions
//			if (iter%iter_aff==0)
				cout<<" # de fusions a realiser = "<<nbarrete2<<"\n";
			for (jj=0; jj<nbarrete2; jj++) {		      // boucle sur le # de fusions
				m=Tarrete2[jj]; 
				if (Tarretes[m].valid) { // il se peut qu'une arrete ait ete invalidee par les fusions precedentes
					n=Tarretes[m].pred;
					nv=Tarretes[m].succ;
					n1=Tregions[n].Lpix.nb_elts(); for (k=0; k<dim; k++) val1[k]=Tregions[n].val[k];
					n2=Tregions[nv].Lpix.nb_elts(); for (k=0; k<dim; k++) val2[k]=Tregions[nv].val[k];
// sommation des pixels frontieres s'il existe une arrete de chaque region n et nv vers une meme autre region
					for (ui=0; ui<Tregions[n].n_arretes; ui++) { 
						m1=Tregions[n].T_arretes[ui]; 
						if (Tarretes[m1].valid) { 
							n3=Tarretes[m1].pred; 
							if (n3==n) n3=Tarretes[m1].succ; 
							for (uj=0; uj<Tregions[nv].n_arretes; uj++) { 
								m2=Tregions[nv].T_arretes[uj]; 
								if (Tarretes[m2].valid) { 
									n4=Tarretes[m2].pred; 
									if (n4==nv) n4=Tarretes[m2].succ; 
									if (n3!=nv && n4!=n && n3==n4)
										Tarretes[m1].length+=Tarretes[m2].length; 
								}
							}
						}
					} 
					Tarretes[m].valid=0;
					nbarretes--; 
// actualisation de la liste des pixels et de la valeur de la region subsistante
					Tregions[n].Lpix+=Tregions[nv].Lpix;//cout<<" operateur += entre listes pixels : a la fin "<<Tregions[n].Lpix.nb_elts()<<" pixels\n";
					/*Tregions[nv].Lpix.vide();*/ Tregions[nv].valid=0; //cout<<" operateur vide de liste pixels : fin\n";
//					Tregions[nv].Lpix.affiche();
					for (k=0; k<dim; k++) Tregions[n].val[k]=(n1*val1[k]+n2*val2[k])/(n1+n2);
// actualisation de la liste des arretes de la region subsistante :
// pour toutes les arretes de la region absorbee : determination de la tierse region cible,
// si cette region n'est pas la region absorbante, et s'il n'y a pas d'arrete la reliant a
// la region absorbante, l'arrete est rajoutee a la liste des arretes de la region absorbante
					n3=0;
					for (ui=0; ui<Tregions[nv].n_arretes; ui++) {
						mv=Tregions[nv].T_arretes[ui];
						bool nv2OK=0;
						if (Tarretes[mv].valid) {
							if (nv==Tarretes[mv].pred) nv2=Tarretes[mv].succ;
							else {
								if (nv==Tarretes[mv].succ) nv2=Tarretes[mv].pred;
								else cout<<" Pb : sommet "<<nv<<" touche ? arrete "<<mv<<" de "<<Tarretes[mv].pred<<" vers "<<Tarretes[mv].succ<<"\n";
							}
							if (nv2!=n) {
								nv2OK=1;
								for (uj=0; uj<Tregions[n].n_arretes; uj++) {
									if (Tarretes[Tregions[n].T_arretes[uj]].valid &&
									    (Tarretes[Tregions[n].T_arretes[uj]].pred==nv2 || Tarretes[Tregions[n].T_arretes[uj]].succ==nv2)) nv2OK=0;
									}
								}
							}
							Tarretes[mv].valid=nv2OK;
							n3+=nv2OK; 
//							cout<<" arrete "<<Tarretes[mv].pred<<" vers "<<Tarretes[mv].succ<<" a transferer ? "<<nv2OK<<"\n";
					}
					if (n3>0) { //cout<<" en tout "<<n3<<" arretes a transferer\n";
						n4=0;
						for (uj=0; uj<Tregions[n].n_arretes; uj++)
							if (Tarretes[Tregions[n].T_arretes[uj]].valid) n4++;
						unsigned long int *T_arretesFus=new unsigned long int[n4+n3];
						jj=0;
						for (uj=0; uj<Tregions[n].n_arretes; uj++) {
							mv=Tregions[n].T_arretes[uj];
							if (Tarretes[mv].valid) T_arretesFus[jj++]=mv;
							}
						if (Tregions[n].T_arretes!=NULL) delete[] Tregions[n].T_arretes;
						for (ui=0; ui<Tregions[nv].n_arretes; ui++) {
							mv=Tregions[nv].T_arretes[ui];
							if (Tarretes[mv].valid) {
								if (nv==Tarretes[mv].pred) Tarretes[mv].pred=n;
								else
									if (nv==Tarretes[mv].succ) Tarretes[mv].succ=n;
									else cout<<"pb!!!!!!!!!!!!!";
								T_arretesFus[jj++]=mv;
							}
						}
						if (jj!=n3+n4) {
							cout<<" Pb ds comptage des arretes valides\n";
							char aa; cin>>aa;
						}
						Tregions[n].T_arretes=new unsigned int[jj];
						for (uj=0; uj<jj; uj++)
							Tregions[n].T_arretes[uj]=T_arretesFus[uj];
						Tregions[n].n_arretes=jj;
						if (T_arretesFus!=NULL) delete[] T_arretesFus;
					}
// destruction de la region absorbee
					if (Tregions[nv].T_arretes!=NULL) delete[] Tregions[nv].T_arretes;
					Tregions[nv].valid=0;
					nbsommets--;
				}
			} // fin de la fusion de toutes les regions qui devaient l'etre
//			cout<<" fin de la fusion de toutes les regions qui devaient l'etre\n";
// actualisation des couts des arretes
			for (m=0; m<nmaxarretes; m++)
				if (Tarretes[m].valid) {
					n1=Tarretes[m].succ; n2=Tarretes[m].pred;
					if (Tregions[n1].valid==0 || Tregions[n2].valid==0) Tarretes[m].valid=0;
					else {
						n3=Tregions[n1].Lpix.nb_elts();	n4=Tregions[n2].Lpix.nb_elts();
						xx=0.; for (k=0; k<dim; k++) xx+=pow(Tregions[n1].val[k]-Tregions[n2].val[k],2);
						Tarretes[m].val=(xx*n3*n4/(n3+n4))/Tarretes[m].length;
					}
				}
//			cout<<" fin actualisation des arretes\n";
		}
		else {
			if (eps==0) eps=1.e-9;
			else eps*=10;
//			cout<<nbarrete2<<" aretes a supprimer ????????? eps = "<<eps<<"\n";
		}
// elimination des arretes non valides (non necessaire a chaque iteration) et recalcul de la longueur de frontieres entre regions
		if (iter%1==0 || nbsommets<=nreg) { //cout<<" elimination des arretes non valides & recalcul de la longueur de frontieres\n";
			for (n=0; n<nmaxsommets; n++)
				if (Tregions[n].valid) {
					n4=0;
					for (uj=0; uj<Tregions[n].n_arretes; uj++)
						if (Tarretes[Tregions[n].T_arretes[uj]].valid) n4++;
					unsigned long int *T_arretesFus=new unsigned long int[n4];
					jj=0;
					for (uj=0; uj<Tregions[n].n_arretes; uj++) {
						mv=Tregions[n].T_arretes[uj];
						if (Tarretes[mv].valid) T_arretesFus[jj++]=mv;
					}
					if (Tregions[n].T_arretes!=NULL) delete[] Tregions[n].T_arretes;
					Tregions[n].T_arretes=new unsigned int[n4];
					for (uj=0; uj<n4; uj++)
						Tregions[n].T_arretes[uj]=T_arretesFus[uj];
					Tregions[n].n_arretes=n4;
					if (T_arretesFus!=NULL) delete[] T_arretesFus;
				}
//			cout<<" fin elimination des arretes non valides\n";
			for (n=0; n<nmaxsommets; n++) {
				if (Tregions[n].valid) {//cout<<" region "<<n<<" valide, contient "<<Tregions[n].Lpix.nb_elts()<<" pixels";
					liste_pixels Lpix2=Tregions[n].Lpix;
					for (k=0; k<dim; k++) val1[k]=0.;
					int npix_n=0;
					while (Lpix2.nb_elts()>0) {
						elt_liste P=Lpix2.extrait();
						if (P.x>=0 && P.x<nblig && P.y>=0 && P.y<nbcol) {
							Tima[P.x*nbcol+P.y]=n;
							for (k=0; k<dim; k++) val1[k]+=ima(P.x,P.y,k);
							npix_n++;
						}
					}
					for (k=0; k<dim; k++) val1[k]/=npix_n;
					xx=0.; for (k=0; k<dim; k++) xx+=pow(val1[k]-Tregions[n].val[k],2.);
					if (npix_n!=Tregions[n].Lpix.nb_elts() || xx>1.e-6) {
						cout<<" region "<<n<<" contient "<<Tregions[n].Lpix.nb_elts()<<" pixels = "<<npix_n<<" ???\n";
						cout<<" Valeur de la region ????????????? (";	for (k=0; k<dim; k++) cout<<Tregions[n].val[k]<<" ";
						cout<<") actualisee a "; for (k=0; k<dim; k++) cout<<val1[k]<<" "; cout<<")\n";
						for (k=0; k<dim; k++) Tregions[n].val[k]=val1[k];
					}
				}
			}
//			cout<<" fin verification valeur et #pixels des regions encore valides\n";
			for (m=0; m<nmaxarretes; m++) {
				if(Tarretes[m].valid) {
					n1=Tarretes[m].pred;
					n2=Tarretes[m].succ;
					if (Tregions[n1].Lpix.nb_elts()<=Tregions[n2].Lpix.nb_elts()) {
						n=n1;	nv=n2;
					}	else {
						n=n2;	nv=n1;
					}
					n3=0;
					for (i=0; i<nblig; i++)
						for (j=0; j<nbcol; j++)
							if (Tima[i*nbcol+j]==n)
								if ((i>0 && Tima[(i-1)*nbcol+j]==nv) || (i<nblig-1 && Tima[(i+1)*nbcol+j]==nv) ||
									(j>0 && Tima[i*nbcol+j-1]==nv) ||	(j<nbcol-1 && Tima[i*nbcol+j+1]==nv) ) n3++;
					if (n3!=Tarretes[m].length) {
//						cout<<" longueur arrete "<<m<<" recalculee : "<<n3;
//						cout<<" au lieu de "<<Tarretes[m].length<<"\n";
						Tarretes[m].val*=Tarretes[m].length/n3;
						Tarretes[m].length=n3;
					}
				}
			}
		}
//		cout<<" fin recalcul de la longueur de frontieres\n";
// Sauvegarde intermediaire de l'image des regions
		if ((nbsommets>1000 && nbsommets%500==0) || (nbsommets<=1000 && nbsommets%100==0) || (nbsommets<=100 && nbsommets%10==0) || (nbsommets<=20)) {
			char nomficI[80], aa[5]; strcpy(aa,"0000"); 
			if (nbsommets<256) strcpy(nomficI,"imasauv0000.pgm"); 
			else strcpy(nomficI,"imasauv0000.dat");
			np=nbsommets; for (int i=4; i>0; i--) {aa[i-1]='0'+np%10; np=np/10;}
			strncpy(nomficI+7,aa,4);
			imadata<int> imasauv(nblig,nbcol); imasauv.mise_a_zero();
			int ireg=valnul; 
			for (n=0; n<nmaxsommets; n++) {
				if (Tregions[n].valid) {
					ireg++;
					Lpix2=Tregions[n].Lpix;
					while (Lpix2.nb_elts()>0) {
						P=Lpix2.extrait();
						if (P.x>=0 && P.x<nblig && P.y>=0 && P.y<nbcol) imasauv(P.x,P.y)=ireg;
					}
				}
			}
			if (nbsommets<256) {imasauv.sauve_ImaPGM(nomficI); cout<<" sauve image "<<nomficI<<"\n";}
			else imasauv.sauve_ImaBSQ(nomficI);
		}
// determination de la prochaine valeur du parametre d'echelle lambda
		lambda=valmin0; mv=0;
		for (m=0; m<nmaxarretes; m++) {
			if(Tarretes[m].valid) {
				n1=Tarretes[m].pred; n2=Tarretes[m].succ;
				n3=Tregions[n1].Lpix.nb_elts();	n4=Tregions[n2].Lpix.nb_elts();
				xx=0.; for (k=0; k<dim; k++) xx+=pow(Tregions[n2].val[k]-Tregions[n1].val[k],2);
				xx=xx*n4*n3/(n3+n4)/Tarretes[m].length;
				if (xx!=Tarretes[m].val) {
					cout<<"####### mise a jour de la valeur d'arrete : "<<xx<<" au lieu de "<<Tarretes[m].val<<"\n";
					Tarretes[m].val=xx;
				}
				if(xx<lambda) lambda=xx;
				mv++;
			}
		}
		lambda+=eps;
		if (iter%iter_aff==0)
			cout<<mv<<" arretes => valeur du prochain parametre d'echelle lambda = "<<lambda<<"\n";
		if (mv==0) {char aa; cin>>aa;}
		if (nbsommets<=nreg) fini=1;	                                       // critere d'arret
//		if (Tarrete2!=NULL) delete[] Tarrete2;
	}
	if (Tarrete2!=NULL) delete[] Tarrete2;
// affichage des caracteristiques des regions
	int ireg=valnul;
	for (n=0; n<nmaxsommets; n++) {
		if (Tregions[n].valid) {
			ireg++;
			if (iter%iter_aff==0) {
				cout<<" region "<<ireg<<" no "<<n<<" a pour valeur ("; for (k=0; k<dim; k++) cout<<Tregions[n].val[k]<<" "; cout<<"), contient ";
				cout<<Tregions[n].Lpix.nb_elts()<<" pixels et a "<<Tregions[n].n_arretes<<" arretes\n";
			}
			for (uj=0; uj<Tregions[n].n_arretes; uj++) {
				m=Tregions[n].T_arretes[uj];
				if (iter%iter_aff==0) {
					cout<<" arrete "<<m<<" : "<<Tarretes[m].pred<<" -> "<<Tarretes[m].succ;
					cout<<" valide? "<<Tarretes[m].valid<<", longueur = "<<Tarretes[m].length<<" cout = "<<Tarretes[m].val<<"\n";
				}
				n1=Tarretes[m].succ; n2=Tarretes[m].pred;
				n3=Tregions[n1].Lpix.nb_elts(); n4=Tregions[n2].Lpix.nb_elts();
			}
			if (iter%iter_aff==0) cout<<"\n";
		}
	}
// creation de l'image des regions
	ireg=valnul;
	for (n=0; n<nmaxsommets; n++) {
		if (Tregions[n].valid) {
			ireg++;
			while (Tregions[n].Lpix.nb_elts()>0) {
				elt_liste P=Tregions[n].Lpix.extrait();
				if (P.x>=0 && P.x<nblig && P.y>=0 && P.y<nbcol) Tima[P.x*nbcol+P.y]=ireg;
			}
			if (Tregions[n].T_arretes!=NULL) delete[] Tregions[n].T_arretes;
		}
	}
	nbregions=ireg-valnul;
	while (!verif_1cc_per_region()) verif_1cc_per_region();
	if (val1!=NULL) delete[] val1; if (val2!=NULL) delete[] val2;
	if (Tregions!=NULL) delete[] Tregions;
	if (Tarretes!=NULL) delete[] Tarretes;
}

imaregions::imaregions (imadata<float> &ima, unsigned int nreg, char *nomfics) : imasites(ima.nlig(),ima.ncol()) { // segmentation Koepfler (fonctionnelle Mumford & Sha simplifiee)
	nblayer=1;
	Tima=new long int[nbpix*nblayer];
	valnul=/*-1*/0;
	int i,j,k=0,np;
	unsigned int ui,uj;
	unsigned long jj;
	for (i=0; i<nbpix; i++) {
		for (j=0; j<nblayer; j++) Tima[i*nblayer+j]=valnul;
		Tsites[i]=&(Tima[i*nblayer]);
	}
	unsigned long int n=0,m=0,nbsommets,nbarretes,nv,mv,nbarrete2,nv2,n1,n2,n3,n4,m1,m2;
	regionV *Tregions=NULL;
	arrete *Tarretes=NULL;
	cout<<" segmentation selon la fonctionnelle simplifiee de Mumford et Shah\n";
	const int init=1, dim=ima.ncanaux(); const bool i_reg1CC=0;
	liste_pixels Lpix2; elt_liste P;
	double xx;
// Initialement=segmentation chaque groupe de pixels de valeur constante est une région
	if (init==1) {
		const int nmax_reg_init=4096;
		const float delta_pas_quantif=1.f;
		int nreg_init, ireg, kreg;
		imadata<float> imaHomoge=ima.init1_segmentationregions(nmax_reg_init,delta_pas_quantif,nreg_init);
		Tregions=new regionV[nreg_init];
		n=0;
		int *T_regcree=new int[nblig*nbcol]; for (i=0; i<nblig*nbcol; i++) T_regcree[i]=-1;
		for (i=0; i<nblig; i++) 
			for (j=0; j<nbcol; j++) {
				ireg=(int)imaHomoge(i,j);
				if (T_regcree[ireg]==-1) {
					Tregions[n].Lpix.insere(i,j); for (k=0; k<dim; k++) Tregions[n].val[k]=ima(i,j,k); Tregions[n].valid=1; Tregions[n].n_arretes=0; 
					T_regcree[ireg]=n; n++; //cout<<" creation sommet "<<n-1<<"\n";
				} else {
					kreg=T_regcree[ireg]; np=Tregions[kreg].Lpix.nb_elts();
					Tregions[kreg].Lpix.insere(i,j); 
					for (k=0; k<dim; k++) Tregions[kreg].val[k]=(np*Tregions[kreg].val[k]+ima(i,j,k))/(np+1);
				}
			}
		nbsommets=n; cout<<" "<<nbsommets<<" sommets crees\n";
// verification de la consistance des regions
		if (i_reg1CC) {
			imabin ima1reg(nblig,nbcol);
			int ncc;
			for (n=0; n<nbsommets; n++)
				if (Tregions[n].valid) {if (n%100==0) cout<<" region "<<n<<" verifiee\n";
					ima1reg.mise_a_zero(); Lpix2=Tregions[n].Lpix;
					while (Lpix2.nb_elts()>0) {
						P=Lpix2.extrait();
						if (P.x>=0 && P.x<nblig && P.y>=0 && P.y<nbcol) ima1reg(P.x,P.y)=1;
					}
					ima1reg.composantes_connexes(ncc,8,0);
					if (ncc!=1) {cout<<" a l'initialisation, PB : region "<<n<<" a plusieurs CC : ncc = "<<ncc<<"\n"; 
						ima1reg.imaunsignedchar().sauve_ImaPGM("ima1reg.pgm"); char aa; cin>>aa;}
//					else cout<<" a l'initialisation, region "<<n<<" Ok : ncc = "<<ncc<<"\n";
				}
		}
// création des arretes (non orientées) APRES création de toutes les régions
		arrete *Tarretes_bis=new arrete[nreg_init*nreg_init/2];
		bool *T_arretecree=new bool[nreg_init*nreg_init]; for (i=0; i<nreg_init*nreg_init; i++) T_arretecree[i]=false;
		unsigned int *Treg_vois=new unsigned int[nreg_init*nreg_init]; for (i=0; i<nreg_init*nreg_init; i++) Treg_vois[i]=false;
		int i0,i2,j0,j2; m=0;
		for (i=0; i<nblig; i++) 
			for (j=0; j<nbcol; j++) { //cout<<" i = "<<i<<" j = "<<j<<"\n";
				i0=maxi(0,i-1); i2=mini(i+1,nblig-1); j0=maxi(0,j-1); j2=mini(j+1,nbcol-1);
				ireg=(int)imaHomoge(i,j); if (T_regcree[ireg]<0) {cout<<" region de label ireg = "<<ireg<<" ne correspond pas a 1 region de T_regcree ????????\n"; char aa; cin>>aa;}
				kreg=(int)imaHomoge(i0,j); if (T_regcree[kreg]<0) {cout<<" region de label kreg = "<<kreg<<" ne correspond pas a 1 region de T_regcree ????????\n"; char aa; cin>>aa;}
				if (ireg!=kreg && T_regcree[ireg]>=0 && T_regcree[kreg]>=0 && !T_arretecree[T_regcree[ireg]*nreg_init+T_regcree[kreg]]) {
					T_arretecree[T_regcree[ireg]*nreg_init+T_regcree[kreg]]=T_arretecree[T_regcree[kreg]*nreg_init+T_regcree[ireg]]=true; m++; //cout<<" creation arrete i0 "<<m<<"\n";
					Treg_vois[T_regcree[ireg]*nreg_init+T_regcree[kreg]]++; Treg_vois[T_regcree[kreg]*nreg_init+T_regcree[ireg]]++;
				}
				kreg=(int)imaHomoge(i2,j); if (T_regcree[kreg]<0) {cout<<" region de label kreg = "<<kreg<<" ne correspond pas a 1 region de T_regcree ????????\n"; char aa; cin>>aa;}
				if (ireg!=kreg && T_regcree[ireg]>=0 && T_regcree[kreg]>=0 && !T_arretecree[T_regcree[ireg]*nreg_init+T_regcree[kreg]]) {
					T_arretecree[T_regcree[ireg]*nreg_init+T_regcree[kreg]]=T_arretecree[T_regcree[kreg]*nreg_init+T_regcree[ireg]]=true; m++; //cout<<" creation arrete i2 "<<m<<"\n";
					Treg_vois[T_regcree[ireg]*nreg_init+T_regcree[kreg]]++; Treg_vois[T_regcree[kreg]*nreg_init+T_regcree[ireg]]++;
				}
				kreg=(int)imaHomoge(i,j0); if (T_regcree[kreg]<0) {cout<<" region de label kreg = "<<kreg<<" ne correspond pas a 1 region de T_regcree ????????\n"; char aa; cin>>aa;}
				if (ireg!=kreg && T_regcree[ireg]>=0 && T_regcree[kreg]>=0 && !T_arretecree[T_regcree[ireg]*nreg_init+T_regcree[kreg]]) {
					T_arretecree[T_regcree[kreg]*nreg_init+T_regcree[ireg]]=T_arretecree[T_regcree[ireg]*nreg_init+T_regcree[kreg]]=true; m++; //cout<<" creation arrete j0 "<<m<<"\n";
					Treg_vois[T_regcree[ireg]*nreg_init+T_regcree[kreg]]++; Treg_vois[T_regcree[kreg]*nreg_init+T_regcree[ireg]]++;
				}
				kreg=(int)imaHomoge(i,j2); if (T_regcree[kreg]<0) {cout<<" region de label kreg = "<<kreg<<" ne correspond pas a 1 region de T_regcree ????????\n"; char aa; cin>>aa;}
				if (ireg!=kreg && T_regcree[ireg]>=0 && T_regcree[kreg]>=0 && !T_arretecree[T_regcree[ireg]*nreg_init+T_regcree[kreg]]) {
					T_arretecree[T_regcree[ireg]*nreg_init+T_regcree[kreg]]=T_arretecree[T_regcree[kreg]*nreg_init+T_regcree[ireg]]=true; m++; //cout<<" creation arrete j2 "<<m<<"\n";
					Treg_vois[T_regcree[ireg]*nreg_init+T_regcree[kreg]]++; Treg_vois[T_regcree[kreg]*nreg_init+T_regcree[ireg]]++;
				}
			}
		nbarretes=m; cout<<" "<<nbarretes<<" arretes creees\n";
		m=0; 
		for (i=0; i<nreg_init; i++)
			for (j=i+1; j<nreg_init; j++) 
				if (T_arretecree[i*nreg_init+j]) {
					Tarretes_bis[m].pred=i; 																									// region de depart de l'arrete
					Tarretes_bis[m].succ=j; 																									// region d'arrivee de l'arrete
					n1=Tregions[i].Lpix.nb_elts(); n2=Tregions[j].Lpix.nb_elts();
					xx=0.; for (k=0; k<dim; k++) xx+=pow(Tregions[i].val[k]-Tregions[j].val[k],2);
					Tarretes_bis[m].val=xx*n1*n2/(n1+n2);  // valeur de l'arrete
					Tarretes_bis[m].valid=1;
					Tarretes_bis[m].length=Treg_vois[i*nreg_init+j];												  // longueur initiale de l'arrete
					m++;
				}
		nbarretes=m; cout<<" "<<nbarretes<<" arretes creees\n";
		for (i=0; i<nreg_init; i++) {
			m=0; for (j=0; j<nreg_init; j++) if (T_arretecree[i*nreg_init+j]) m++;
			Tregions[i].T_arretes=new unsigned int[m];   // allocation dynamique du tableau des arrêtes pour chaque région
		}
		for (m=0; m<nbarretes; m++)  // boucle sur le # max d'arretes car les arretes valides peuvent etre n'importe ou dans le tableau
			if (Tarretes_bis[m].valid) {							  // si l'arrete existe, pour chaque sommet
				n1=Tarretes_bis[m].pred;							  // on determine le numero du sommet
				Tregions[n1].T_arretes[Tregions[n1].n_arretes]=m; // l'arrete est associee au sommet
				Tregions[n1].n_arretes++;						  // le # d'arretes du sommet est incremente
				n2=Tarretes_bis[m].succ;
				Tregions[n2].T_arretes[Tregions[n2].n_arretes]=m;
				Tregions[n2].n_arretes++;
			}
		if (T_regcree!=NULL) delete[] T_regcree; if (T_arretecree!=NULL) delete[] T_arretecree; 
		if (Treg_vois!=NULL) delete[] Treg_vois;
		Tarretes=new arrete[nbarretes]; for (m=0; m<nbarretes; m++) Tarretes[m]=Tarretes_bis[m];
		if (Tarretes_bis!=NULL) delete[] Tarretes_bis;
	}
	if (init==0) {
		int ii,jj,i0,i2,j0,j2,ncc,ntcc;
		const float eps=0.f;
		float x, dx, dxmin=(float)1.e+9, xmin=(float)ima.minI(), xmax=(float)ima.maxI(), xx;
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				x=ima(i,j);
				i0=maxi(i-1,0); i2=mini(i+1,nblig-1); j0=maxi(j-1,0); j2=mini(j+1,nbcol-1);
				for (ii=i0; ii<=i2; ii++)
					for (jj=j0; jj<=j2; jj++) {
						dx=fabs(x-ima(ii,jj));
						if (dx>0 && dx<dxmin) dxmin=dx;
					}
			}
		cout<<" dxmin="<<dxmin<<"\n";
		imabin ima_x(nblig,nbcol);
		imadata<int> ima_xcc(nblig,nbcol), ima_cc(nblig,nbcol);
		ima_cc.mise_a_zero(); ntcc=0;
		for (xx=xmin; xx<=xmax; xx+=dxmin) {
			ima_x.mise_a_zero();
			for (i=0; i<nblig; i++)
				for (j=0; j<nbcol; j++) {
					x=ima(i,j);
					if (x>=xx && x<=xx+eps) ima_x(i,j)=1;
				}
			ima_xcc.mise_a_zero();
			ima_xcc=ima_x.composantes_connexes(ncc,8,0);
			for (i=0; i<nblig; i++)
				for (j=0; j<nbcol; j++) 
					if (ima_xcc(i,j)>0) ima_cc(i,j)=ima_xcc(i,j)+ntcc;
			ntcc+=ncc;
		}
		cout<<" initialement "<<ntcc<<" regions\n"; (imadata<float>(ima_cc)).sauve_Ima("ima_cc.dat");
		cree_graph_init (ima_cc,ima,Tregions,Tarretes,nbarretes,nbsommets);
/*		Tregions=new region[ntcc];
		unsigned int **Treg_vois=new unsigned int*[ntcc];
		for (n=0; n<ntcc; n++) Treg_vois[n]=new unsigned int[ntcc];
		for (n=0; n<ntcc; n++) 
			for (j=0; j<ntcc; j++) Treg_vois[n][j]=0;
		nbarretes=0;
		for (n=0; n<ntcc; n++) {
			val_cc=0.; ncc=0; m=0;
			for (i=0; i<nblig; i++)
				for (j=0; j<nbcol; j++)
					if (ima_cc(i,j)==n+1) {
						Tregions[n].Lpix.insere(i,j);            // la region a les pixels de la CC
						val_cc+=ima(i,j); ncc++;
						i0=maxi(i-1,0); i2=mini(i+1,nblig-1); j0=maxi(j-1,0); j2=mini(j+1,nbcol-1);
//						for (ii=i0; ii<=i2; ii++)
//							for (jj=j0; jj<=j2; jj++) 
//								if (ima_cc(ii,jj)>=1) Treg_vois[n][ima_cc(ii,jj)-1]+=1; // cas 8-connexite
						for (ii=i0; ii<=i2; ii++)
							for (jj=j0; jj<=j2; jj++) 
								if ((ii==i || jj==j) && ima_cc(ii,jj)>=1) Treg_vois[n][ima_cc(ii,jj)-1]+=1; // cas 4-connexite
					}
			Tregions[n].val=val_cc/ncc;              // sa valeur est la moyenne des pixels de la CC
			Tregions[n].valid=1;                     // elle est valide
			Tregions[n].n_arretes=0;                 // aucune arete n'a ete definie a ce stade
			for (j=0; j<ntcc; j++) 
				if (j!=n && Treg_vois[n][j]>0) m+=1; 
//			cout<<" la region "<<n+1<<" a "<<m<<" regions voisines\n";
			Tregions[n].T_arretes=new unsigned long int[m];   // allocation du tableau des arrêtes
			nbarretes+=m;
		}
		cout<<" le graphe des regions a "<<nbarretes/2<<" arretes en tout\n"; 
		nbsommets=ntcc;
		nbarretes/=2;
		Tarretes=new arrete[nbarretes];
		m=0;
		for (n=0; n<ntcc; n++)
			for (j=n+1; j<ntcc; j++)
				if (Treg_vois[n][j]>0) {
					Tarretes[m].pred=n;			                        // region de depart de l'arrete
					Tarretes[m].succ=j;	                            // region d'arrivee de l'arrete
					n1=Tregions[n].Lpix.nb_elts(); n2=Tregions[j].Lpix.nb_elts();
					Tarretes[m].val=pow(Tregions[n].val-Tregions[j].val,2)*n1*n2/(n1+n2); // valeur de l'arrete
					Tarretes[m].valid=1;
					Tarretes[m].length=Treg_vois[n][j];		          // longueur initiale de l'arrete
					m++;
				}
		for (m=0; m<nbarretes; m++)  // boucle sur le # max d'arretes car les arretes valides peuvent etre n'importe ou dans le tableau
			if (Tarretes[m].valid) {							  // si l'arrete existe, pour chaque sommet
				n1=Tarretes[m].pred;							  // on determine le numero du sommet
				Tregions[n1].T_arretes[Tregions[n1].n_arretes]=m; // l'arrete est associee au sommet
				Tregions[n1].n_arretes++;						  // le # d'arretes du sommet est incremente
				n2=Tarretes[m].succ;
				Tregions[n2].T_arretes[Tregions[n2].n_arretes]=m;
				Tregions[n2].n_arretes++;
			}
		for (n=0; n<ntcc; n++)
			if (Treg_vois[n]!=NULL) delete[] Treg_vois[n];
		if (Treg_vois!=NULL) delete[] Treg_vois;*/
	}
	if (init==-1) {
// Initialement=segmentation triviale, i.e. chaque pixel est une région
		Tregions=new regionV[nbpix]; n=0;
// Initialement, chaque region a autant d'arrêtes que de voisins en k-connexité (ici k=4)
		Tarretes=new arrete[nblig*(nbcol-1)+(nblig-1)*nbcol]; m=0;
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				Tregions[n].Lpix.insere(i,j);            // la region a 1 seul pixel : (i,j)
				for (k=0; k<dim; k++) Tregions[n].val[k]=ima(i,j,k);              // sa valeur est celle de son unique pixel
				Tregions[n].valid=1;                     // elle est valide
				Tregions[n].n_arretes=0;                 // aucune arete n'a ete definie a ce stade
				if (i>0) {                                        // creation des arretes verticales
					Tarretes[m].pred=n;			                        // region de depart de l'arrete
					Tarretes[m].succ=n-nbcol;	                      // region d'arrivee de l'arrete
					xx=0.; for (k=0; k<dim; k++) xx+=pow(ima(i-1,j,k)-ima(i,j,k),2);
					Tarretes[m].val=xx/2; // valeur de l'arrete
					Tarretes[m].valid=1;
					Tarretes[m].length=1;		                        // longueur initiale de l'arrete
					m++;
				}
				if (j>0) {                                        // creation des arretes horizontales
					Tarretes[m].pred=n; Tarretes[m].succ=n-1;
					Tarretes[m].val=pow(ima(i,j-1,k)-ima(i,j,k),2)/2;
					Tarretes[m].valid=1; Tarretes[m].length=1;
					m++;
				}
				Tregions[n].T_arretes=new unsigned int[4];   // allocation du tableau des arrêtes
				n++;
			}
		nbsommets=n; nbarretes=m;
		for (m=0; m<nbarretes; m++)  // boucle sur le # max d'arretes car les arretes valides peuvent etre n'importe ou dans le tableau
			if (Tarretes[m].valid) {							  // si l'arrete existe, pour chaque sommet
				n1=Tarretes[m].pred;							  // on determine le numero du sommet
				Tregions[n1].T_arretes[Tregions[n1].n_arretes]=m; // l'arrete est associee au sommet
				Tregions[n1].n_arretes++;						  // le # d'arretes du sommet est incremente
				n2=Tarretes[m].succ;
				Tregions[n2].T_arretes[Tregions[n2].n_arretes]=m;
				Tregions[n2].n_arretes++;
			}
}
// Variables et initialisation
	double *val1=new double[dim], *val2=new double[dim]; 
	const double valmin0=FLT_MAX;
	double lambda=valmin0;
	const unsigned long int nmaxsommets=nbsommets, nmaxarretes=nbarretes;
// boucle calculant le nombre d'arretes par region
	for (m=0; m<nmaxarretes; m++)  // boucle sur le # max d'arretes car les arretes valides peuvent etre n'importe ou dans le tableau
		if (Tarretes[m].valid) {							  // si l'arrete existe, pour chaque sommet
//			n1=Tarretes[m].pred;							  // on determine le numero du sommet
//			Tregions[n1].T_arretes[Tregions[n1].n_arretes]=m; // l'arrete est associee au sommet
//			Tregions[n1].n_arretes++;						  // le # d'arretes du sommet est incremente
//			n2=Tarretes[m].succ;
//			Tregions[n2].T_arretes[Tregions[n2].n_arretes]=m;
//			Tregions[n2].n_arretes++;
//			n3=Tregions[n1].Lpix.nb_elts();           // lambda minimal pour avoir au moins 1 fusion
//			n4=Tregions[n2].Lpix.nb_elts();
//			xx=(pow(Tregions[n2].val - Tregions[n1].val,2)*n4*n3/(n3+n4))/Tarretes[m].length;
			xx=Tarretes[m].val;
			if (xx<lambda) {
				lambda=xx; 
//				cout<<" fusion regions "<<n1<<"("<<Tregions[n1].val<<","<<n3<<") & "<<n2<<"("<<Tregions[n2].val<<","<<n4<<") si lambda>="<<lambda<<" ("<<Tarretes[m].length<<")\n";
			}
		}
//	{cout<<" fin init "; char aa; cin>>aa;}
	double eps=/*0.*/1.e-2;
	lambda+=eps;
	cout<<" valeur du prochain parametre d'echelle lambda = "<<lambda<<"\n";
// iterations sur la fusion des regions
	bool fini=0;
	unsigned int iter=0;
	int iter_aff=2;
	while (!fini) { //cout<<iter<<" "<<iter_aff<<" "<<iter%iter_aff<<" ";
//		if ((iter%iter_aff)==0)
		cout<<" iteration "<<iter++<<" : # regions = "<<nbsommets<<"\n";
		unsigned long int * Tarrete2=new unsigned long int[nbarretes]; // tableau des arretes a supprimer (taille = # courant d'arretes)
		nbarrete2=0;
// determination des arretes a supprimer
		for (m=0; m<nmaxarretes; m++)
			if (Tarretes[m].valid && Tarretes[m].val<=lambda) Tarrete2[nbarrete2++]=m;
// fusion des regions
		if (nbarrete2>0) {                                // si il y a des fusions de regions
//			if (iter%iter_aff==0)
				cout<<" # de fusions a realiser = "<<nbarrete2<<"\n";
			for (jj=0; jj<nbarrete2; jj++) {		      // boucle sur le # de fusions
				m=Tarrete2[jj];
				if (Tarretes[m].valid) { // il se peut qu'une arrete ait ete invalidee par les fusions precedentes
					n=Tarretes[m].pred;
					nv=Tarretes[m].succ;
					n1=Tregions[n].Lpix.nb_elts(); for (k=0; k<dim; k++) val1[k]=Tregions[n].val[k];
					n2=Tregions[nv].Lpix.nb_elts(); for (k=0; k<dim; k++) val2[k]=Tregions[nv].val[k];
// sommation des pixels frontieres s'il existe une arrete de chaque region n et nv vers une meme autre region
					for (ui=0; ui<Tregions[n].n_arretes; ui++) {
						m1=Tregions[n].T_arretes[ui];
						if (Tarretes[m1].valid) {
							n3=Tarretes[m1].pred;
							if (n3==n) n3=Tarretes[m1].succ;
							for (uj=0; uj<Tregions[nv].n_arretes; uj++) {
								m2=Tregions[nv].T_arretes[uj];
								if (Tarretes[m2].valid) {
									n4=Tarretes[m2].pred;
									if (n4==nv) n4=Tarretes[m2].succ;
									if (n3!=nv && n4!=n && n3==n4)
										Tarretes[m1].length+=Tarretes[m2].length;
								}
							}
						}
					}
					Tarretes[m].valid=0;
					nbarretes--;
// actualisation de la liste des pixels et de la valeur de la region subsistante
//					while (Tregions[nv].Lpix.nb_elts()>0)
//						Tregions[n].Lpix.insere(Tregions[nv].Lpix.extrait());
//					cout<<" regions "<<n<<"&"<<nv<<" operateur += entre listes pixels de "<<Tregions[n].Lpix.nb_elts()<<"+"<<Tregions[nv].Lpix.nb_elts()<<" pixels\n";
//					Tregions[n].Lpix.affiche(); cout<<"***\n";Tregions[nv].Lpix.affiche();
					Tregions[n].Lpix+=Tregions[nv].Lpix;//cout<<" operateur += entre listes pixels : a la fin "<<Tregions[n].Lpix.nb_elts()<<" pixels\n";
//					Tregions[n].Lpix.affiche();
					/*Tregions[nv].Lpix.vide();*/ Tregions[nv].valid=0; //cout<<" operateur vide de liste pixels : fin\n";
//					Tregions[nv].Lpix.affiche();
					for (k=0; k<dim; k++) Tregions[n].val[k]=(n1*val1[k]+n2*val2[k])/(n1+n2);
//					{char aa; cin>>aa;}
// actualisation de la liste des arretes de la region subsistante :
// pour toutes les arretes de la region absorbee : determination de la tierse region cible,
// si cette region n'est pas la region absorbante, et s'il n'y a pas d'arrete la reliant a
// la region absorbante, l'arrete est rajoutee a la liste des arretes de la region absorbante
					n3=0;
					for (ui=0; ui<Tregions[nv].n_arretes; ui++) {
						mv=Tregions[nv].T_arretes[ui];
						bool nv2OK=0;
						if (Tarretes[mv].valid) {
							if (nv==Tarretes[mv].pred) nv2=Tarretes[mv].succ;
							else {
								if (nv==Tarretes[mv].succ) nv2=Tarretes[mv].pred;
								else cout<<" Pb : sommet "<<nv<<" touche ? arrete "<<mv<<" de "<<Tarretes[mv].pred<<" vers "<<Tarretes[mv].succ<<"\n";
							}
							if (nv2!=n) {
								nv2OK=1;
								for (uj=0; uj<Tregions[n].n_arretes; uj++) {
									if (Tarretes[Tregions[n].T_arretes[uj]].valid &&
									    (Tarretes[Tregions[n].T_arretes[uj]].pred==nv2 || Tarretes[Tregions[n].T_arretes[uj]].succ==nv2)) nv2OK=0;
									}
								}
							}
							Tarretes[mv].valid=nv2OK;
							n3+=nv2OK; 
//							cout<<" arrete "<<Tarretes[mv].pred<<" vers "<<Tarretes[mv].succ<<" a transferer ? "<<nv2OK<<"\n";
					}
					if (n3>0) { //cout<<" en tout "<<n3<<" arretes a transferer\n";
						n4=0;
						for (uj=0; uj<Tregions[n].n_arretes; uj++)
							if (Tarretes[Tregions[n].T_arretes[uj]].valid) n4++;
						unsigned long int *T_arretesFus=new unsigned long int[n4+n3];
						jj=0;
						for (uj=0; uj<Tregions[n].n_arretes; uj++) {
							mv=Tregions[n].T_arretes[uj];
							if (Tarretes[mv].valid) T_arretesFus[jj++]=mv;
							}
						if (Tregions[n].T_arretes!=NULL) delete[] Tregions[n].T_arretes;
						for (ui=0; ui<Tregions[nv].n_arretes; ui++) {
							mv=Tregions[nv].T_arretes[ui];
							if (Tarretes[mv].valid) {
								if (nv==Tarretes[mv].pred) Tarretes[mv].pred=n;
								else
									if (nv==Tarretes[mv].succ) Tarretes[mv].succ=n;
									else cout<<"pb!!!!!!!!!!!!!";
								T_arretesFus[jj++]=mv;
							}
						}
						if (jj!=n3+n4) {
							cout<<" Pb ds comptage des arretes valides\n";
							char aa; cin>>aa;
						}
						Tregions[n].T_arretes=new unsigned int[jj];
						for (uj=0; uj<jj; uj++)
							Tregions[n].T_arretes[uj]=T_arretesFus[uj];
						Tregions[n].n_arretes=jj;
						if (T_arretesFus!=NULL) delete[] T_arretesFus;
					}
// destruction de la region absorbee
					if (Tregions[nv].T_arretes!=NULL) delete[] Tregions[nv].T_arretes;
					Tregions[nv].valid=0;
					nbsommets--;
				}
			} // fin de la fusion de toutes les regions qui devaient l'etre
//			cout<<" fin de la fusion de toutes les regions qui devaient l'etre\n";
// actualisation des couts des arretes
			for (m=0; m<nmaxarretes; m++)
				if (Tarretes[m].valid) {
					n1=Tarretes[m].succ; n2=Tarretes[m].pred;
					if (Tregions[n1].valid==0 || Tregions[n2].valid==0) Tarretes[m].valid=0;
					else {
						n3=Tregions[n1].Lpix.nb_elts();	n4=Tregions[n2].Lpix.nb_elts();
						xx=0.; for (k=0; k<dim; k++) xx+=pow(Tregions[n1].val[k]-Tregions[n2].val[k],2);
						Tarretes[m].val=(xx*n3*n4/(n3+n4))/Tarretes[m].length;
					}
				}
//			cout<<" fin actualisation des arretes\n";
		}
		else {
			if (eps==0) eps=1.e-9;
			else eps*=10;
//			cout<<nbarrete2<<" aretes a supprimer ????????? eps = "<<eps<<"\n";
		}
// elimination des arretes non valides (non necessaire a chaque iteration)
//  et recalcul de la longueur de frontieres entre regions
		if (iter%1==0 || nbsommets<=nreg) { //cout<<" elimination des arretes non valides & recalcul de la longueur de frontieres\n";
			for (n=0; n<nmaxsommets; n++)
				if (Tregions[n].valid) {
					n4=0;
					for (uj=0; uj<Tregions[n].n_arretes; uj++)
						if (Tarretes[Tregions[n].T_arretes[uj]].valid) n4++;
					unsigned long int *T_arretesFus=new unsigned long int[n4];
					jj=0;
					for (uj=0; uj<Tregions[n].n_arretes; uj++) {
						mv=Tregions[n].T_arretes[uj];
						if (Tarretes[mv].valid) T_arretesFus[jj++]=mv;
					}
					if (Tregions[n].T_arretes!=NULL) delete[] Tregions[n].T_arretes;
					Tregions[n].T_arretes=new unsigned int[n4];
					for (uj=0; uj<n4; uj++)
						Tregions[n].T_arretes[uj]=T_arretesFus[uj];
					Tregions[n].n_arretes=n4;
					if (T_arretesFus!=NULL) delete[] T_arretesFus;
				}
//			cout<<" fin elimination des arretes non valides\n";
			for (n=0; n<nmaxsommets; n++) {
				if (Tregions[n].valid) {//cout<<" region "<<n<<" valide, contient "<<Tregions[n].Lpix.nb_elts()<<" pixels";
					liste_pixels Lpix2=Tregions[n].Lpix;
					for (k=0; k<dim; k++) val1[k]=0.;
					int npix_n=0;
					while (Lpix2.nb_elts()>0) {
						elt_liste P=Lpix2.extrait();
						if (P.x>=0 && P.x<nblig && P.y>=0 && P.y<nbcol) {
							Tima[P.x*nbcol+P.y]=n;
							val1[k]+=ima(P.x,P.y,k);
							npix_n++;
						}
					}
					for (k=0; k<dim; k++) val1[k]/=npix_n; 
					xx=0.; for (k=0; k<dim; k++) xx+=fabs(val1[k]-Tregions[n].val[k]); 
					if (npix_n!=Tregions[n].Lpix.nb_elts() || xx>1.e-6) {
						cout<<" region "<<n<<" contient "<<Tregions[n].Lpix.nb_elts()<<" pixels = "<<npix_n<<" ???\n";
						cout<<" Valeur de la region ????????????? ("; for (k=0; k<dim; k++) cout<<Tregions[n].val[k]<<" "; 
						cout<<") actualisee a "; for (k=0; k<dim; k++) cout<<val1[k]<<" "; cout<<"\n";
						for (k=0; k<dim; k++) Tregions[n].val[k]=val1[k];
					}
				}
			}
//			cout<<" fin verification valeur et #pixels des regions encore valides\n";
			for (m=0; m<nmaxarretes; m++) {
				if(Tarretes[m].valid) {
					n1=Tarretes[m].pred;
					n2=Tarretes[m].succ;
					if (Tregions[n1].Lpix.nb_elts()<=Tregions[n2].Lpix.nb_elts()) {
						n=n1;
						nv=n2;
					}
					else {
						n=n2;
						nv=n1;
					}
					n3=0;
					for (i=0; i<nblig; i++)
						for (j=0; j<nbcol; j++)
							if (Tima[i*nbcol+j]==n)
								if ((i>0 && Tima[(i-1)*nbcol+j]==nv) ||
									(i<nblig-1 && Tima[(i+1)*nbcol+j]==nv) ||
									(j>0 && Tima[i*nbcol+j-1]==nv) ||
									(j<nbcol-1 && Tima[i*nbcol+j+1]==nv) ) n3++;
					if (n3!=Tarretes[m].length) {
//						cout<<" longueur arrete "<<m<<" recalculee : "<<n3;
//						cout<<" au lieu de "<<Tarretes[m].length<<"\n";
						Tarretes[m].val*=Tarretes[m].length/n3;
						Tarretes[m].length=n3;
					}
				}
			}
		}
//		cout<<" fin recalcul de la longueur de frontieres\n";
// Sauvegarde intermediaire de l'image des regions
		if ((nbsommets>1000 && nbsommets%500==0) || (nbsommets<=1000 && nbsommets%100==0) || (nbsommets<=100 && nbsommets%10==0) || (nbsommets<=20)) {
			char nomficI[80], aa[5]; strcpy(nomficI,"imasauv0000.dat"); strcpy(aa,"0000"); 
			np=nbsommets; for (int i=4; i>0; i--) {aa[i-1]='0'+np%10; np=np/10;}
			strncpy(nomficI+7,aa,4);
			imadata<int> imasauv(nblig,nbcol); imasauv.mise_a_zero();
			int ireg=valnul; 
			for (n=0; n<nmaxsommets; n++) {
				if (Tregions[n].valid) {
					ireg++;
					Lpix2=Tregions[n].Lpix;
					while (Lpix2.nb_elts()>0) {
						P=Lpix2.extrait();
						if (P.x>=0 && P.x<nblig && P.y>=0 && P.y<nbcol) imasauv(P.x,P.y)=ireg;
					}
				}
			}
			imasauv.sauve_ImaBSQ(nomficI);
		}
// determination de la prochaine valeur du parametre d'echelle lambda
		lambda=valmin0; mv=0;
		for (m=0; m<nmaxarretes; m++) {
			if(Tarretes[m].valid) {
				n1=Tarretes[m].pred; n2=Tarretes[m].succ;
				n3=Tregions[n1].Lpix.nb_elts();	n4=Tregions[n2].Lpix.nb_elts();
				xx=0.; for (k=0; k<dim; k++) xx+=pow(Tregions[n2].val[k]-Tregions[n1].val[k],2);
				xx=(xx*n4*n3/(n3+n4))/Tarretes[m].length;
				if (xx!=Tarretes[m].val) {
					cout<<"####### mise a jour de la valeur d'arrete : "<<xx<<" au lieu de "<<Tarretes[m].val<<"\n";
					Tarretes[m].val=xx;
				}
				if(xx<lambda) lambda=xx;
				mv++;
			}
		}
		lambda+=eps;
		if (iter%iter_aff==0)
			cout<<mv<<" arretes => valeur du prochain parametre d'echelle lambda = "<<lambda<<"\n";
		if (mv==0) {char aa; cin>>aa;}
		if (nbsommets<=nreg) fini=1;	                                       // critere d'arret
		if (Tarrete2!=NULL) delete[] Tarrete2;
	}
// affichage des caracteristiques des regions
	int ireg=valnul;
	for (n=0; n<nmaxsommets; n++) {
		if (Tregions[n].valid) {
			ireg++;
			if (iter%iter_aff==0) {
				cout<<" region "<<ireg<<" no "<<n<<" a pour valeur "<<Tregions[n].val<<", contient ";
				cout<<Tregions[n].Lpix.nb_elts()<<" pixels et a "<<Tregions[n].n_arretes<<" arretes\n";
			}
			for (uj=0; uj<Tregions[n].n_arretes; uj++) {
				m=Tregions[n].T_arretes[uj];
				if (iter%iter_aff==0) {
					cout<<" arrete "<<m<<" : "<<Tarretes[m].pred<<" -> "<<Tarretes[m].succ;
					cout<<" valide? "<<Tarretes[m].valid<<", longueur = "<<Tarretes[m].length<<" cout = "<<Tarretes[m].val<<"\n";
				}
				n1=Tarretes[m].succ; n2=Tarretes[m].pred;
				n3=Tregions[n1].Lpix.nb_elts(); n4=Tregions[n2].Lpix.nb_elts();
			}
			if (iter%iter_aff==0) cout<<"\n";
		}
	}
// creation de l'image des regions
	ireg=valnul;
	for (n=0; n<nmaxsommets; n++) {
		if (Tregions[n].valid) {
			ireg++;
			while (Tregions[n].Lpix.nb_elts()>0) {
				elt_liste P=Tregions[n].Lpix.extrait();
				if (P.x>=0 && P.x<nblig && P.y>=0 && P.y<nbcol) Tima[P.x*nbcol+P.y]=ireg;
			}
			if (Tregions[n].T_arretes!=NULL) delete[] Tregions[n].T_arretes;
		}
	}
	nbregions=ireg-valnul;
	verif_1cc_per_region();
	if (val1!=NULL) delete[] val1; if (val2!=NULL) delete[] val2; 
	if (Tregions!=NULL) delete[] Tregions;
	if (Tarretes!=NULL) delete[] Tarretes;
}

imaregions::imaregions (imadata<float> &ima, int iconnex, bool ielimLPE) : imasites(ima.nlig(),ima.ncol()) { // LPE : ligne de partage des eaux
	nblayer=1;
	Tima=new long int[nbpix*nblayer];
	valnul=-1;
	int i,j,k=0,l,n;
	for (long int ii=0; ii<nbpix; ii++) {
		for (j=0; j<nblayer; j++) Tima[ii*nblayer+j]=valnul;
		Tsites[ii]=&(Tima[ii*nblayer]);
	}
	s_homog=-1;
// calcul de l'image des bassins versants
	const float coef255=(float)mini(255./(ima.maxI()-ima.minI()),1.), coef_nblevels=1.f/*0.05f*/;
	imadata<BYTE> imaX((ima+(float)(-ima.minI()))*coef255); imaX.statbasic(1);
	if (fabs(coef_nblevels-1.f)>1.e-6) {imaX=ima.imaunsignedchar(1)*coef_nblevels; imaX.statbasic(1);}
	imabin imaXi(nblig,nbcol), imaXi_1(nblig,nbcol), imaXi_rgeo(nblig,nbcol), imaIZi(nblig,nbcol), imaWi(nblig,nbcol);
	imaXi.mise_a_zero();
	bool icheck=1;
	for (i=0; i<256*coef_nblevels; i++) {
		n=0;
		for (l=0; l<nblig; l++)
			for (k=0; k<nbcol; k++)
				if (imaX(l,k)==i) {imaXi(l,k)=1; n++;}	//image des niveaux inférieurs à i
		if (n>0) {cout<<" niveau "<<i<<" : "<<n<<" pixels ; ";
			if (i>0) {
				imaXi_rgeo=imaXi.reconstruction_geodesique(imaXi_1,iconnex);
//				imaXi_rgeo.imaunsignedchar().sauve_ImaPGM("imaXi_rgeo.pgm"); cout<<" reconstruction geodesique OK\n";
				imaIZi=imaXi.zones_influence_geodesique(imaWi,icheck,iconnex); 
//				imaIZi.imaunsignedchar().sauve_ImaPGM("imaIZi.pgm"); cout<<" zones d'influence geodesique OK\n";
				for (l=0; l<nblig; l++)
					for (k=0; k<nbcol; k++) if ((imaXi(l,k)&&!imaXi_rgeo(l,k))||imaIZi(l,k)) imaWi(l,k)=1;
			} else imaWi=imaXi;
			imaXi_1=imaXi;
//			imaWi.imaunsignedchar().sauve_ImaPGM("imaWi.pgm");
		}
	}
// image des regions deduite
	int nregions=0;
	imadata<int> imacc=imaWi.composantes_connexes(nregions,iconnex);
//	imadata<int> imacc=imaWi.CompoConnexes_MM(nregions);
	imabin imab(nblig+2,nbcol+2);
	for (l=0; l<nblig+2; l++) imab(l,0)=imab(l,nbcol+1)=1;
	for (k=0; k<nbcol+2; k++) imab(0,k)=imab(nblig+1,k)=1;
	for (l=0; l<nblig; l++)
		for (k=0; k<nbcol; k++) {
			if (imacc(l,k)>0) (*this)(l,k)=imacc(l,k);
			else {(*this)(l,k)=0; imab(l+1,k+1)=1;}
		}
	nbregions=nregions;
	cout<<" au final "<<nbregions<<" regions construites\n";
	imabin imabel=imab.elagage(40); imabel.imaunsignedchar().sauve_ImaPGM("imaWiEl.pgm");
	cout<<" avant elagage "<<imab.norm()<<" apres elagage "<<imabel.norm()<<"\n";
	long int T[8][2];
	int l0,l2,k0,k2,m;
	bool trouve;
	for (l=0; l<nblig; l++)
		for (k=0; k<nbcol; k++) {
			if (imab(l,k) && !imabel(l,k)) {
				for (i=0; i<8; i++)
					for (j=0; j<2; j++) T[i][j]=0;
				l0=maxi(l-1,0); l2=mini(l+1,nblig-1); k0=maxi(k-1,0); k2=mini(k+1,nbcol-1);
				m=0;
				for (i=l0; i<=l2; i++)
					for (j=k0; j<=k2; j++)
						if (i!=l || j!=k) {
							trouve=0;
							for (n=0; n<m; n++)
								if (!trouve && T[n][0]==(*this)(i,j)) {T[n][1]++; trouve=1;}
							if (!trouve && m<8) {T[m][0]=(*this)(i,j); T[m][1]=1; m++;}
						}
//				for (n=0; n<m; n++) if (T[n][1]>0) cout<<T[n][0]<<" "<<T[n][1]<<", "; cout<<"\n";
				i=-1;
				for (n=0; n<m; n++)
					if (T[n][1]>i || ((*this)(l,k)==0 && T[n][1]>=i)) {i=T[n][1]; (*this)(l,k)=T[n][0];}
				i=(*this)(l,k);
				if (i!=0) 
					for (n=0; n<m; n++)
						if (T[n][0]!=0 && T[n][0]!=i && T[n][1]>0) (*this)(l,k)=0;
			}
		}
	if (ielimLPE) {
		for (l=0; l<nblig; l++)
			for (k=0; k<nbcol; k++) {
				if ((*this)(l,k)==0) {
					for (i=0; i<8; i++)
						for (j=0; j<2; j++) T[i][j]=0;
					l0=maxi(l-1,0); l2=mini(l+1,nblig-1); k0=maxi(k-1,0); k2=mini(k+1,nbcol-1);
					m=0;
					for (i=l0; i<=l2; i++)
						for (j=k0; j<=k2; j++)
							if (i!=l || j!=k && (*this)(i,j)>0) {
								trouve=0;
								for (n=0; n<m; n++)
									if (!trouve && T[n][0]==(*this)(i,j)) {T[n][1]++; trouve=1;}
								if (!trouve && m<8) {T[m][0]=(*this)(i,j); T[m][1]=1; m++;}
							}
					i=-1;
					for (n=0; n<m; n++)
						if (T[n][1]>i || ((*this)(l,k)==0 && T[n][1]>=i)) {i=T[n][1]; (*this)(l,k)=T[n][0];}
/*					i=(*this)(l,k);
					if (i!=0) 
						for (n=0; n<m; n++)
							if (T[n][0]!=0 && T[n][0]!=i && T[n][1]>0) (*this)(l,k)=0;*/
				}
			}
	}
	cout<<" au final "<<nbregions<<" regions construites\n";
}

void imaregions::suppression_trop_petites_regions (int npixmin) { cout<<" suppression_trop_petites_regions : initialement "<<nbregions<<"\n";
	if (npixmin>0) {
		int *npix_reg=new int[nbregions+1],i,j;
		long int k,l;
		for (k=0; k<=nbregions; k++) npix_reg[k]=0;
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				k=(*this)(i,j);
				if (k>=0 && k<=nbregions) npix_reg[k]++;
			}
		for (k=0; k<=nbregions; k++) if (npix_reg[k]>0) cout<<" region "<<k<<" contient "<<npix_reg[k]<<" pixels\n";
		for (k=nbregions; k>0; k--) {
			if (npix_reg[k]<npixmin && k<=nbregions) {
				for (i=0; i<nblig; i++)
					for (j=0; j<nbcol; j++) if ((*this)(i,j)==k) (*this)(i,j)=0;
				for (i=0; i<nblig; i++)
					for (j=0; j<nbcol; j++) if ((*this)(i,j)>(int)k) (*this)(i,j)-=1;
				for (l=k; l<nbregions; l++) npix_reg[l]=npix_reg[l+1];
				npix_reg[nbregions]=0;
				nbregions--; k=(k+1,nbregions+1);
			}
		}
		cout<<" apres elimimination des regions de moins de "<<npixmin<<" pixels, il reste "<<nbregions<<" regions :\n";
		for (k=0; k<=nbregions; k++) cout<<" region "<<k<<" contient "<<npix_reg[k]<<" pixels\n";
		if (npix_reg!=NULL) delete[] npix_reg;
	}
}

void imaregions::suppression_plus_petites_regions (int nregmax) { cout<<" suppression_plus_petites_regions : initialement "<<nbregions<<"\n";
	if (nbregions>nregmax) {
//		int nblig=imareg.nlig(),nbcol=imareg.ncol(),nreg=imareg.nregions(), ;
		int n,i,j,k,ireg,*T_nbpix=new int[nbregions+1];
		double *T_nbpix_bis=new double[nbregions+1];
		for (n=0; n<nbregions+1; n++) T_nbpix[n]=0;
		for (i=0; i<nblig; i++) for (j=0; j<nbcol; j++) T_nbpix[(*this)(i,j)]++;
		for (n=0; n<nbregions+1; n++) cout<<" "<<T_nbpix[n]; cout<<"\n";
		cout<<" segmentation comprend "<<nbregions<<" regions -> trop de regions -> on supprime les "<<nbregions-nregmax<<" regions plus petites en nombre de pixels\n";
		for (n=0; n<nbregions+1; n++) T_nbpix_bis[n]=(double)T_nbpix[n];
		tri_rapide(T_nbpix_bis,nbregions+1);
		for (n=0; n<nbregions+1; n++) cout<<" "<<(int)T_nbpix_bis[n]; cout<<"\n";
		n=0;
		while (n<=nbregions-nregmax) {
			ireg=0; while (ireg<nbregions+1 && T_nbpix[ireg]!=T_nbpix_bis[n]) ireg++;
			if (ireg>=nbregions+1) cout<<" Pb....\n";
			else {
//				cout<<" elimination de la region "<<ireg<<" qui contient "<<T_nbpix[ireg]<<" pixels\n";
				for (i=0; i<nblig; i++) for (j=0; j<nbcol; j++) if ((*this)(i,j)==ireg) (*this)(i,j)=0;
				T_nbpix[ireg]=-1;
			}
			n++;
		}
//		for (n=1; n<nreg+1; n++) {
		for (n=1; n<=nregmax; n++)
			if (T_nbpix[n]<=0) {
				k=maxi(n+1,nregmax+1); while (k<nbregions+1 && T_nbpix[k]<=0) k++; 
				if (k<nbregions+1) {
//					cout<<" region "<<k<<" devient "<<n<<"\n";
					for (i=0; i<nblig; i++) for (j=0; j<nbcol; j++) if ((*this)(i,j)==k) (*this)(i,j)=n;
					T_nbpix[n]=T_nbpix[k]; T_nbpix[k]=-1;
				}
			}
//		for (n=0; n<nbregions+1; n++) T_nbpix[n]=0; for (i=0; i<nblig; i++) for (j=0; j<nbcol; j++) T_nbpix[imareg(i,j)]++;
//		for (n=0; n<nbregions+1; n++) cout<<" region "<<n<<" contient "<<T_nbpix[n]<<" pixels\n";
		n=0; for (i=0; i<nblig; i++) for (j=0; j<nbcol; j++) if ((*this)(i,j)>n) n=(*this)(i,j);
		cout<<" il reste donc "<<n<<" regions\n"; //imareg.nregions()=n;
		if (T_nbpix!=NULL) delete[] T_nbpix; if (T_nbpix_bis!=NULL) delete[] T_nbpix_bis;
	}
}

void imaregions::affiche() const {
	int i,j,k;
	long int *adval;
	long int ireg;
	for (k=0; k<nblayer; k++) {
		cout<<"*************** image : 'couche' "<<k+1<<" ***************\n";
		for (i=0; i<nblig; i++) {
			for (j=0; j<nbcol; j++) {
				adval=(long int *)Tsites[i*nbcol+j];
				ireg=*(adval+k);
				cout<<(int)ireg<<" ";
				}
			cout<<"\n";
			}
		cout<<"----------------------------------------------------------\n";
		}
	}

void imaregions::sauve_Ima(string nomfich, int icanal) const {
	char *nomfic=(char *)nomfich.c_str();
	sauve_Ima(nomfic,icanal);
}

void imaregions::sauve_Ima(char *nomfich, int icanal) const {
	ofstream sortie;
	sortie.open(nomfich,ios::out|ios::binary);
	if (!sortie) {
		cout<<" ouverture de "<<nomfich<<" impossible\n";
		exit (-1);
		}
	long int* adval=(long int*)Tsites[0];
	int itype=2;
	cout<<" image "<<nblig<<" lig.& "<<nbcol<<" col., de type bytes => itype = "<<itype<<"\n";
	sortie.write((char *)&nblig,sizeof(int));
	sortie.write((char *)&nbcol,sizeof(int));
	sortie.write((char *)&itype,sizeof(int));
	for (long int i=0; i<nbpix; i++) {
		adval=(long int*)Tsites[i]+icanal;
		sortie.write((char *)adval,sizeof(long int));
		}
	sortie.close();
	}

void imaregions::sauve_ImaUc(string nomfich, int icanal) const {
	char *nomfic=(char *)nomfich.c_str();
	sauve_ImaUc(nomfic,icanal);
}

void imaregions::sauve_ImaUc(char *nomfich,int icanal) const {
	int i,j,k,labmax=0;
	imadata<BYTE> I1c(nblig,nbcol,1);
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) {
//			k=(*this)(i,j,icanal);
			k=Tima[(i*nbcol+j)*nblayer+icanal];
			if (k>labmax) labmax=k;
			}
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) {
			I1c(i,j,0)=0;
//			k=(*this)(i,j,icanal);
			k=Tima[(i*nbcol+j)*nblayer+icanal];
			if (k>0) {
				k=k*255/labmax;
				I1c(i,j,0)=k;
				}
			}
	I1c.sauve_ImaBSQ(nomfich);
	}

