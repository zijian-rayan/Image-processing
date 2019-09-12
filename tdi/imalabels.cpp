#include "imalabels.h"
// #include "graph.h"
#include "imabin.h"

void imalabels::cmptvs4connex (int i, int j, int k, int nclas, int *Nvs) {
	int i0,i2,j0,j2,ii,jj,icl;
	i0=maxi(i-1,0);
	i2=mini(i+2,nblig);
	j0=maxi(j-1,0);
	j2=mini(j+2,nbcol);
	for (ii=0; ii<=nclas; ii++) Nvs[ii]=0;
	for (ii=i0; ii<i2; ii++)
		for (jj=j0; jj<j2; jj++) {
			if (ii==i || jj==j) {
				icl=(int)Tima[(ii*nbcol+jj)*nblayer+k];
				if (icl>0 && icl<=nclas) Nvs[icl]++;
				}
			}
	for (ii=1; ii<=nclas; ii++) Nvs[0]+=Nvs[ii];
	}

void imalabels::cmptvs8connex (int i, int j, int k, int nclas, int *Nvs) {
	int i0,i2,j0,j2,ii,jj,icl;
	i0=maxi(i-1,0);
	i2=mini(i+2,nblig);
	j0=maxi(j-1,0);
	j2=mini(j+2,nbcol);
	for (ii=0; ii<=nclas; ii++) Nvs[ii]=0;
	for (ii=i0; ii<i2; ii++)
		for (jj=j0; jj<j2; jj++) {
			icl=(int)Tima[(ii*nbcol+jj)*nblayer+k];
			if (icl>0 && icl<=nclas) Nvs[icl]++;
			}
	Nvs[(int)Tima[(i*nbcol+j)*nblayer+k]]--;
	for (ii=1; ii<=nclas; ii++) Nvs[0]+=Nvs[ii];
	}

double imalabels::U_lineprocess (int j0, int j2, int i0, int i2, float a0, float a1, float a2) {
	double x=0;
	int n=0;
	if (j0!=0) n++; if (j2!=0) n++; if (i0!=0) n++; if (i2!=0) n++;
	if (n>0) {
		if (n==2) {
			if (j0!=0 && j2!=0 || i0!=0 && i2!=0) x+=a0;
			else x+=a1;
		}
		if (n==1 || n==4) x+=a2;
		if (n==3) x+=a1;
	}
	return x;
}

void imalabels::cmptvs4connexLP (int i, int j, int nclas, int *Nvs) {
	int ipixl,ipix,icl;
	for (icl=0; icl<=nclas; icl++) Nvs[icl]=0;
	ipixl=i*nbcol;
	if (i>0) {
		ipix=(ipixl-nbcol+j)*nblayer;                   //pixel (i-1,j)
		if (Tima[ipix+3]==0) {
			icl=(int)Tima[ipix+1];
			if (icl>0 && icl<=nclas) Nvs[icl]++;
			}
		}
	if (i<nblig-1) {
		ipix=(ipixl+j)*nblayer;                         //pixel (i,j)
		if (Tima[ipix+3]==0) {
			icl=(int)Tima[ipix+nbcol*nblayer+1];         //pixel (i+1,j)
			if (icl>0 && icl<=nclas) Nvs[icl]++;
			}
		}
	if (j>0) {
		ipix=(ipixl+j-1)*nblayer;                       //pixel (i,j-1)
		if (Tima[ipix+2]==0) {
			icl=(int)Tima[ipix+1];
			if (icl>0 && icl<=nclas) Nvs[icl]++;
			}
		}
	if (i<nblig-1) {
		ipix=(ipixl+j)*nblayer;                         //pixel (i,j)
		if (Tima[ipix+2]==0) {
			icl=(int)Tima[ipix+nblayer+1];               //pixel (i,j+1)
			if (icl>0 && icl<=nclas) Nvs[icl]++;
			}
		}
	for (icl=1; icl<=nclas; icl++) Nvs[0]+=Nvs[icl];
   }

imadata<float> imalabels::imaU0condClas(imadata<float> &ima, statclasses Tclas) {
	const int nclas=Tclas.nclas();
	int i,j,icl;
	long int ipixl;
	pixel<float> pix, dist;
	imadata<float> imaU0(nblig,nbcol,nclas+1);
	for (i=0; i<nblig; i++) {
		ipixl=i*nbcol;
		for (j=0; j<nbcol; j++) {
			pix=ima.pix(i,j);
			dist=Tclas.distclasses(pix);
			icl=(int)dist[0];
			imaU0.pix(i,j,dist);
			}
		}
	return imaU0;
}

imalabels::imalabels(imadata<float> &ima, statclasses Tclas, float bPotts, char *algo, float pct_ch) : imasites(ima.nlig(),ima.ncol()) {
	Tima=NULL;
	const int nclas=Tclas.nclas(), itermax=10000, nit_raznch=1, nchmax=(int)(pct_ch*nbpix);
	const float rTemp=0.999f;//0.998f;
	float nch=0.f, temp=2000.f;//10.f;
	int i,j,icl,ii, *Nvs=new int[nclas+1], iter=0;
	long int l,ipixl;
	double *Papriori=new double[nclas+1], xx, xcumul, xmin;
	char *algoICM="ICM", *algoSAG="SAG", *algoSAM="SAM";
	cout<<" classification MRF "<<nclas<<" classes, parametre a priori beta "<<bPotts<<", algorithme "<<algo<<"\n";
	nblayer=2;
	Tima=new BYTE[nbpix*nblayer];
	pixel<float> pix, dist;
	for (l=0; l<nbpix; l++) {
		for (j=0; j<nblayer; j++) Tima[l*nblayer+j]=0; Tsites[l]=&(Tima[l*nblayer]);
	}
	imadata<float> imaU0(nblig,nbcol,nclas+1);
	imadata<BYTE> imaSauv(nblig,nbcol); int factD=(int)floor(255.f/nclas);
	for (i=0; i<nblig; i++) {               // calcul du terme d'attache aux donnees
		ipixl=i*nbcol;                                   // stocke sur la couche 0
		for (j=0; j<nbcol; j++) {
			pix=ima.pix(i,j); dist=Tclas.distclasses(pix);
			icl=(int)dist[0]; imaU0.pix(i,j,dist);
			imaSauv(i,j)=Tima[(ipixl+j)*nblayer]=Tima[(ipixl+j)*nblayer+1]=(BYTE)icl;
		}
	}
	char nit[5]="0000", nomfic[100]="classif_iter0000.pgm"; (imaSauv*(float)factD).sauve_ImaPGM(nomfic);
	do {
		iter=iter+1; temp=temp*rTemp;
		if (iter%nit_raznch==0) nch=0.f;
		for (i=0; i<nblig; i++) {
			ipixl=i*nbcol;
			for (j=0; j<nbcol; j++) {
				for (ii=0; ii<=nclas; ii++) Papriori[ii]=0.;
				cmptvs8connex (i,j,1,nclas,Nvs); //cmptvs4connex (i,j,1,nclas,Nvs);
				if (strcmp(algo,algoSAM)==0) {
					icl=(int)Tima[(i*nbcol+j)*nblayer+1];
					do ii=(rand()%nclas)+1; while (ii==icl);
					Papriori[icl]=imaU0(i,j,icl)+bPotts*(Nvs[0]-Nvs[icl]);
					Papriori[ii]=imaU0(i,j,ii)+bPotts*(Nvs[0]-Nvs[ii]);
					xx=Papriori[ii]-Papriori[icl];
					if (xx<0) icl=ii;
					else {
						xcumul=exp(-xx/temp);
						if ((double)rand()/RAND_MAX<xcumul) icl=ii;
					}
				}
				else {
					for (ii=1; ii<=nclas; ii++) Papriori[ii]=(imaU0(i,j,ii)+bPotts*(Nvs[0]-Nvs[ii]))/temp;
					xx=0.; for (ii=1; ii<=nclas; ii++) xx+=exp(-Papriori[ii]);
					if (xx!=0. && strcmp(algo,algoSAG)==0) {
						for (ii=1; ii<=nclas; ii++) Papriori[ii]=exp(-Papriori[ii])/xx;
						xx=(double)rand()/RAND_MAX;
						xcumul=0.; icl=0;
						do {icl++; xcumul+=Papriori[icl]; } while (xcumul<xx && icl<nclas);
						if (xcumul<xx) icl=(int)Tima[(ipixl+j)*nblayer+1];
					}
					else {                                      //notamment algo="ICM
						xmin=1.e+9;
						icl=0;
						for (ii=1; ii<=nclas; ii++) {
							if (Papriori[ii]<xmin) {icl=ii; xmin=Papriori[ii];}
						}
					}
				}
				if (icl<=0 || icl>nclas) cout<<"Pb icl "<<i<<" "<<j<<" : "<<icl<<" "<<nclas<<"\n";
				if ((int)Tima[(ipixl+j)*nblayer+1] != icl) {
					imaSauv(i,j)=Tima[(ipixl+j)*nblayer+1]=(BYTE)icl; nch=nch+1.f;}
			}
		}
		if (strcmp(algo,algoICM)==0 || iter%100==0) {
			cout<<" iter. "<<setw(3)<<iter<<" temp. "<<setw(6)<<setprecision(3)<<temp;
			cout<<" => #changements = "<<setw(6)<<nch<<"\n";
			num_2_char_blabla(iter,nit,4); cout<<" iteration "<<nit; strncpy(nomfic+12,nit,4); cout<<" sauvegarde dans "<<nomfic<<"\n";;
			(imaSauv*(float)factD).sauve_ImaPGM(nomfic);
		}
	} while (iter<itermax && nch>nchmax);
	cout<<" fin de la convergence pour "<<iter<<" iterations, derniere it : "<<nch<<" pixels changes\n";
	if (Nvs!=NULL) delete[] Nvs;
	if (Papriori!=NULL) delete[] Papriori;
	}

/*

imalabels::imalabels(imadata<float> & ima, float bPotts, statclasses Tclas) : imasites(ima.nlig(),ima.ncol()) {
	Tima=NULL;
	const int nclas=Tclas.nclas();
	int i,j,icl,i0,i2,ii,j0,j2,jj;
	long int l,ipixl;
	bool fin=0;
	cout<<" classification MRF "<<nclas<<" classes, parametre a priori beta "<<bPotts<<", optimisation par graph-cut\n";
	nblayer=2;
	Tima=new BYTE[nbpix*nblayer];
	pixel<float> pix, dist;
	for (l=0; l<nbpix; l++) {
		for (j=0; j<nblayer; j++) Tima[l*nblayer+j]=0;
		Tsites[l]=&(Tima[l*nblayer]);
		}
	imadata<float> imaU0(nblig,nbcol,nclas+1);
	for (i=0; i<nblig; i++) {               // calcul du terme d'attache aux donnees
		ipixl=i*nbcol;
		for (j=0; j<nbcol; j++) {
			pix=ima.pix(i,j);
			dist=Tclas.distclasses(pix);
			icl=(int)dist[0];
			Tima[(ipixl+j)*nblayer]=(BYTE)icl;
			Tima[(ipixl+j)*nblayer+1]=(BYTE)icl;
			imaU0.pix(i,j,dist);
		}
	}
	graph G;
	cout<<" debut creation du graphe\n";
	G.ajoute(node(nbpix+1)); //cout<<" noeud "<<nbpix+1<<"\n";
	G.ajoute(node(nbpix+2)); //cout<<" noeud "<<nbpix+2<<"\n";
	for (i=0; i<nblig; i++) {
		ipixl=i*nbcol+1;
		if (i%10==0) {cout<<" ligne "<<i<<"\n";}
		for (j=0; j<nbcol; j++) {
			G.ajoute(node(ipixl+j)); 
			G.ajoute(arc(nbpix+1,ipixl+j,imaU0(i,j,1))); 
			G.ajoute(arc(ipixl+j,nbpix+2,imaU0(i,j,2))); 
		}
	}
	cout<<" fin de creation des noeuds et aretes vers S et T => "<<G.n_nodes()<<" noeuds\n";
	for (i=0; i<nblig; i++) {
		if (i%10==0) cout<<" ligne "<<i<<"\n";
		for (j=0; j<nbcol; j++) {
			ipixl=i*nbcol+j+1;
			i0=maxi(i-1,0); i2=mini(i+1,nblig-1);
			j0=maxi(j-1,0); j2=mini(j+1,nbcol-1);
			for (ii=i0; ii<=i2; ii++)
				for (jj=j0; jj<=j2; jj++)
//					if (i!=ii || j!=jj) G.ajoute(arc(ipixl,ii*nbcol+jj+1,bPotts));
					if (abs(i-ii)+abs(j-jj)==1) G.ajoute(arc(ipixl,ii*nbcol+jj+1,bPotts)); // 4-connexite
		}
	}
	cout<<" fin de creation des aretes d'adjacence des pixels-noeuds => "<<G.n_arcs()<<" arcs\n";
//	G.affiche(); char aa; cin>>aa;
	while (!fin) {
		fin=G.flotmax_FordFukerson(nbpix+1,nbpix+2);
		int n1=0, n2=0, n0=0;
		for (i=0; i<nblig; i++) {
			ipixl=i*nbcol+1;
			for (j=0; j<nbcol; j++) {
				if (G.sature(nbpix+1,ipixl+j) && !G.sature(ipixl+j,nbpix+2)) {
					Tima[(ipixl+j-1)*nblayer+1]=(BYTE)1;
					n1++;
				}
				else  {
					if (G.sature(ipixl+j,nbpix+2) && !G.sature(nbpix+1,ipixl+j)) {
						Tima[(ipixl+j-1)*nblayer+1]=(BYTE)2;
						n2++;
					}
					else {
						n0++; 
						Tima[(ipixl+j-1)*nblayer+1]=(BYTE)0;
					}
				}
			}
		}
		cout<<n1<<" pixels dans la classes 0, "<<n2<<" pixels dans la classe 1, "<<n0<<" pixels non encore labelises\n";
		sauve3d_Ima("./classifGraphCut_test0",1);
		int k, n=1, iter=0, naffect, *Tn=new int[3], *T_n=new int[3];
		bool iaf=0;
		while (n>0 && iter<50) {
			n=0; iter++;
			for (i=0; i<nblig; i++) {
				ipixl=i*nbcol;
				for (j=0; j<nbcol; j++) {
					if (Tima[(ipixl+j)*nblayer+1]==(BYTE)0) {
						if (iaf) {
							cout<<" pixel de coordonnees "<<i<<","<<j<<" non labelise "<<Tima[(ipixl+j)*nblayer+1]<<" a pour voisins :\n";
							if (j>0) cout<<i<<","<<j-1<<" de label "<<(int)Tima[(ipixl+j-1)*nblayer+1]<<"\n";
							if (j<nbcol-1) cout<<i<<","<<j+1<<" de label "<<(int)Tima[(ipixl+j+1)*nblayer+1]<<"\n";
							if (i>0) cout<<i-1<<","<<j<<" de label "<<(int)Tima[((i-1)*nbcol+j)*nblayer+1]<<"\n";
							if (i<nblig-1) cout<<i+1<<","<<j<<" de label "<<(int)Tima[((i+1)*nbcol+j)*nblayer+1]<<"\n";
						}
						for (k=0; k<3; k++) Tn[k]=0;
						if (j>0 && Tima[(ipixl+j-1)*nblayer+1]!=(BYTE)0)
							if (!G.sature(ipixl+j,ipixl+j+1) && !G.sature(ipixl+j+1,ipixl+j)) {
								Tima[(ipixl+j)*nblayer+1]=Tima[(ipixl+j-1)*nblayer+1];
								Tn[Tima[(ipixl+j)*nblayer+1]]++;
							}
						if (j<nbcol-1 && Tima[(ipixl+j+1)*nblayer+1]!=(BYTE)0)
							if (!G.sature(ipixl+j+1,ipixl+j+2) && !G.sature(ipixl+j+2,ipixl+j+1)) {
								Tima[(ipixl+j)*nblayer+1]=Tima[(ipixl+j+1)*nblayer+1];
								Tn[Tima[(ipixl+j)*nblayer+1]]++;
							}
						if (i>0 && Tima[((i-1)*nbcol+j)*nblayer+1]!=(BYTE)0)
							if (!G.sature((i-1)*nbcol+j+1,ipixl+j+1) && !G.sature(ipixl+j+1,(i-1)*nbcol+j+1)) {
								Tima[(ipixl+j)*nblayer+1]=Tima[((i-1)*nbcol+j)*nblayer+1];
								Tn[Tima[(ipixl+j)*nblayer+1]]++;
							}
						if (i<nblig-1 && Tima[((i+1)*nbcol+j)*nblayer+1]!=(BYTE)0)
							if (!G.sature((i+1)*nbcol+j+1,ipixl+j+1) && !G.sature(ipixl+j+1,(i+1)*nbcol+j+1)) {
								Tima[(ipixl+j)*nblayer+1]=Tima[((i+1)*nbcol+j)*nblayer+1];
								Tn[Tima[(ipixl+j)*nblayer+1]]++;
							}
						naffect=0;
						for (k=0; k<3; k++) {if (iaf) cout<<k<<" "<<Tn[k]<<" ; "; if (Tn[k]>0) naffect++;}
						if (naffect==1) n++;
						else
							if (naffect>1)
								Tima[(ipixl+j)*nblayer+1]=(BYTE)0;
						if (iaf) cout<<" => affectation finale = "<<(int)Tima[(ipixl+j)*nblayer+1]<<"\n";
					}
				}
			}
			cout<<n<<" pixels supplementaires labelises\n";
		}
		if (fin) {
			n=1;
			while (n>0) {
				n=0;
				for (i=0; i<nblig; i++)
					for (j=0; j<nbcol; j++)
						if (Tima[(i*nbcol+j)*nblayer+1]==(BYTE)0) {
							if (iaf) {
								cout<<" pixel ("<<i<<","<<j<<") labelise "<<(int)Tima[(i*nbcol+j)*nblayer+1]<<" correspond au noeud "<<i*nbcol+j+1<<"\n";
								G.affichenoeud(i*nbcol+j+1);
								cout<<" derniere tentative de labelisation\n";
							}
							ipixl=i*nbcol;
							for (k=0; k<3; k++) Tn[k]=0;
							if (j>0 && Tima[(ipixl+j-1)*nblayer+1]!=(BYTE)0)
								if (!G.sature(ipixl+j,ipixl+j+1) || !G.sature(ipixl+j+1,ipixl+j)) {
									Tima[(ipixl+j)*nblayer+1]=Tima[(ipixl+j-1)*nblayer+1];
									Tn[Tima[(ipixl+j)*nblayer+1]]++;
								}
							if (j<nbcol-1 && Tima[(ipixl+j+1)*nblayer+1]!=(BYTE)0)
								if (!G.sature(ipixl+j+1,ipixl+j+2) || !G.sature(ipixl+j+2,ipixl+j+1)) {
									Tima[(ipixl+j)*nblayer+1]=Tima[(ipixl+j+1)*nblayer+1];
									Tn[Tima[(ipixl+j)*nblayer+1]]++;
								}
							if (i>0 && Tima[((i-1)*nbcol+j)*nblayer+1]!=(BYTE)0)
								if (!G.sature((i-1)*nbcol+j+1,ipixl+j+1) || !G.sature(ipixl+j+1,(i-1)*nbcol+j+1)) {
									Tima[(ipixl+j)*nblayer+1]=Tima[((i-1)*nbcol+j)*nblayer+1];
									Tn[Tima[(ipixl+j)*nblayer+1]]++;
								}
							if (i<nblig-1 && Tima[((i+1)*nbcol+j)*nblayer+1]!=(BYTE)0)
								if (!G.sature((i+1)*nbcol+j+1,ipixl+j+1) || !G.sature(ipixl+j+1,(i+1)*nbcol+j+1)) {
									Tima[(ipixl+j)*nblayer+1]=Tima[((i+1)*nbcol+j)*nblayer+1];
									Tn[Tima[(ipixl+j)*nblayer+1]]++;
								}
							naffect=0;
							for (k=0; k<3; k++) {if (iaf) cout<<k<<" "<<Tn[k]<<" ; "; if (Tn[k]>0) naffect++;}
							if (naffect==1) n++;
							else
								if (naffect>1)
									Tima[(ipixl+j)*nblayer+1]=(BYTE)0;
						}
				cout<<n<<" pixels supplementaires labelises\n";
			}
		}
		if (Tn!=NULL) delete[] Tn;
		if (T_n!=NULL) delete[] T_n;
//		cout<<n1<<" pixels dans la classes 0, "<<n2<<" pixels dans la classe 1, "<<n0<<" pixels non encore labelises\n";
		sauve3d_Ima("./classifGraphCut_test1",1);
	}
}
*/


void imalabels::affiche() {
	int i,j,k;
	BYTE *adval;
	BYTE icl;
	for (k=0; k<nblayer; k++) {
		cout<<"*************** image : 'couche' "<<k+1<<" ***************\n";
		for (i=0; i<nblig; i++) {
			for (j=0; j<nbcol; j++) {
				adval=(BYTE *)Tsites[i*nbcol+j];
				icl=*(adval+k);
				cout<<(int)icl<<" ";
				}
			cout<<"\n";
			}
		cout<<"----------------------------------------------------------\n";
		}
	}

imalabels imalabels::operator - (imalabels & ima2) {
	int nlig=mini(nblig,ima2.nblig), ncol=mini(nbcol,ima2.nbcol), 
		nlayer=mini(nblayer,ima2.nblayer), i, j, k;
	imalabels imaRes(nlig,ncol,nlayer);
	for (i=0; i<nlig; i++)
		for (j=0; j<ncol; j++)
			for (k=0; k<nlayer; k++) {
				if ((*this)(i,j,k)==ima2(i,j,k)) imaRes(i,j,k)=1;
				else imaRes(i,j,k)=0;
				}
	return imaRes;
	}

void imalabels::correspondances_labels (imalabels &ima2, bool **Tcorresp, bool affich) {
	int c1=0,c2=0,cc=0,nbrlig=mini(nblig,ima2.nlig()),nbrcol=mini(nbcol,ima2.ncol()),k=0;
	int i,j,ii,jj,l1,l2,n_iter,nch,ntirages,n1,n2,n3,nOK;
	bool accept_chg, iOK;
	double coef_iter0=20.,coef_beta=1.-0.75,ntiragesmax=1000,U2maxim=0,U,rTemp=0.9,Temp,xx,c;
	for (i=0; i<nbrlig; i++)
		for (j=0; j<nbrcol; j++) {
			if ((*this)(i,j,k)>c1) c1=(*this)(i,j,k);
			if (ima2(i,j,k)>c2) c2=ima2(i,j,k);
			}
	if (affich) cout<<" image 1 : "<<c1<<" classes, comparee a image 2 : "<<c2<<" classes\n";
	int n_iter0=(int)(coef_iter0*c1*c2),n_itermax=100*n_iter0,beta=(int)(coef_beta*nbrlig*nbrcol/maxi(c1,c2));
	double Temp0=0.*nbrlig*nbrcol;
	bool **Tcorres2=NULL;
	Tcorres2=new bool*[c1+1];
	for (i=0; i<=c1; i++)
		Tcorres2[i]=new bool[c2+1];
	int** H2d=NULL;
	H2d=new int*[c1+1];
	for (i=0; i<=c1; i++) H2d[i]=new int[c2+1];
	bool* T_E=new bool[maxi(c1,c2)+1];
// initialisations
	for (i=0; i<=c1; i++)
		for (j=0; j<=c2; j++) H2d[i][j]=0;
	for (i=0; i<nbrlig; i++)
		for (j=0; j<nbrcol; j++) {
			l1=(*this)(i,j,k);
			l2=ima2(i,j,k);
			if (l1>=0 && l1<=c1 && l2>=0 && l2<=c2) H2d[l1][l2]++;
			}
	for (i=0; i<=c1; i++)
		for (j=0; j<=c2; j++) Tcorresp[i][j]=0;
	if (c1<c2) {
		for (i=1; i<=c1; i++) Tcorresp[i][i]=1;
		for (i=c1+1; i<=c2; i++) Tcorresp[c1][i]=1;
		}
	else {
		for (i=1; i<=c2; i++) Tcorresp[i][i]=1;
		for (i=c2+1; i<=c1; i++) Tcorresp[i][c2]=1;
		}
	U2maxim=0;
	for (i=1; i<=c1; i++)
		for (j=1; j<=c2; j++) 
			if (Tcorresp[i][j]) U2maxim+=H2d[i][j];
	U2maxim+=beta*mini(c1,c2);
	cc=mini(c1,c2);
	Temp=Temp0;
	n_iter=0;
// test d'une nouvelle configuration : 
	do {
		n_iter++;
		Temp*=rTemp;
		if (n_iter%n_iter0==1) nch=0;
		for (i=1; i<=c1; i++)
			for (j=1; j<=c2; j++) Tcorres2[i][j]=Tcorresp[i][j];
// tirage d'une nouvelle correspondance
		ntirages=0;
		do {
			l1=(rand()%c1)+1;
			l2=(rand()%c2)+1;
			xx=(double)rand()/(double)RAND_MAX;
			ntirages++;
			} while (ntirages<ntiragesmax && (l1<1  || l1>c1 || l2<1 || l2>c2 || Tcorresp[l1][l2] || xx==0.5) );
		Tcorresp[l1][l2]=1;
		accept_chg=1;
// verification de la contrainte (3)
		n1=0;
		for (i=1; i<=c1; i++) n1+=Tcorresp[i][l2];
		n2=0;
		for (j=1; j<=c2; j++) n2+=Tcorresp[l1][j];
		if (n1>1 && n2>1) {
//			cout<<" condition (3) non respectee en ("<<l1<<","<<l2<<") : \n";
			if (xx>0.5) {//cout<<" => forcage de la somme des tij a 1 en colonne\n";
				for (i=1; i<=c1; i++) 
					if (i!=l1 && Tcorresp[i][l2]) {
						Tcorresp[i][l2]=0;
						n1=0;
						for (j=1; j<=c2; j++) n1+=Tcorresp[i][j];
						if (n1==0) {
							nOK=0;
							for (j=1; j<=c2; j++) {
								iOK=0;
								n2=0;
								for (ii=1; ii<=c1; ii++) {
									n2+=Tcorresp[ii][j];
									if (n2==1) {
										n3=0;
										for (jj=1; jj<=c2; jj++) n3+=Tcorresp[ii][jj];
										if (n3==1) iOK=1;
										}
									}
								T_E[j]=iOK;
								nOK+=iOK;
								}
							if (nOK>0) {
								jj=(rand()%nOK)+1;
								j=0;
								do {
									j++;
									jj-=T_E[j];
									} while (jj>0 && j<=c2);
								Tcorresp[i][j]=1;
								}
							else accept_chg=0;
							}
						}
				}
			else {//cout<<" => forcage de la somme des tij a 1 en ligne\n";
				for (j=1; j<=c2; j++) 
					if (j!=l2 && Tcorresp[l1][j]) {
						Tcorresp[l1][j]=0;
						n2=0;
						for (i=1; i<=c1; i++) n2+=Tcorresp[i][j];
						if (n2==0) {
							nOK=0;
							for (i=1; i<=c1; i++) {
								iOK=0;
								n1=0;
								for (jj=1; jj<=c2; jj++) {
									n1+=Tcorresp[i][jj];
									if (n1==1) {
										n3=0;
										for (ii=1; ii<=c1; ii++) n3+=Tcorresp[ii][jj];
										if (n3==1) iOK=1;
										}
									}
								T_E[i]=iOK;
								nOK+=iOK;
								}
							if (nOK>0) {
								ii=(rand()%nOK)+1;
								i=0;
								do {
									i++;
									ii-=T_E[i];
									} while (ii>0 && i<=c1);
								Tcorresp[i][j]=1;
								}
							else accept_chg=0;
							}
						}
				}
			}
		if (accept_chg) {
			for (l1=1; l1<=c1; l1++)
				for (l2=1; l2<=c2; l2++) {
					if (Tcorresp[l1][l2]) {
						n1=0;
						for (i=1; i<=c1; i++) n1+=Tcorresp[i][l2];
						n2=0;
						for (j=1; j<=c2; j++) n2+=Tcorresp[l1][j];
						if ((n1>1 && n2>1) || (n1==0 && n2==0)) {
//							cout<<" condition (3) non respectee en ("<<l1<<","<<l2<<") : \n";
//							cout<<" => mise a "<<!Tcorresp[l1][l2]<<" du terme ("<<l1<<","<<l2<<")\n";
							Tcorresp[l1][l2]=!Tcorresp[l1][l2];
							}
						}
					}
			c=0.;
			for (i=1; i<=c1; i++) {
				n1=0;
				jj=0;
				for (j=1; j<=c2; j++) {
					n1+=Tcorresp[i][j];
					if (jj==0 && n1==1) jj=j;
					}
				if (n1>1) c+=1.;
				else {
					n2=0;
					for (ii=1; ii<=c1; ii++) n2+=Tcorresp[ii][jj];
					if (n2==0) {
						cout<<" somme des termes de la colonne "<<jj<<" nulle\n";
						cout<<" somme des termes de la ligne "<<i<<" = "<<n1<<"\n";
						}
					else c+=1./(double)n2;
/*					xx=0.;
					for (j=1; j<=c2; j++) {
						if (Tcorresp[i][j]) {
							n2=0;
							for (ii=1; ii<=c1; ii++) n2+=Tcorresp[ii][j];
							xx+=(double)Tcorresp[i][j]/(double)n2;
							}
						}
					c+=xx;*/
					}
				}
			if (fabs(c-around(c))>1.e-6) 
				cout<<" pb dans le calcul du # de classes 'communes' = "<<c<<"!="<<around(c)<<"\n";
			U=0;
			for (i=1; i<=c1; i++)
				for (j=1; j<=c2; j++) 
					if (Tcorresp[i][j]) U+=H2d[i][j];
			U+=beta*c;
			if (U<=U2maxim && (Temp==0 ||
				(double)rand()/(double)RAND_MAX>=exp((U-U2maxim)/Temp)) ) accept_chg=0;
			}
		if (!accept_chg) {
			for (i=1; i<=c1; i++)
				for (j=1; j<=c2; j++) Tcorresp[i][j]=Tcorres2[i][j];
			}
		else {
			nch++;
			U2maxim=U;
			cc=around(c);
			}
		if (affich && n_iter%n_iter0==0) {
			cout<<" iteration "<<n_iter<<" temp = "<<Temp<<" : # changements = "<<nch<<"\n";
			}
	} while (n_iter%n_iter0!=0 || (nch>0 && n_iter<n_itermax));
// affichage des resultats et fin
//	cout<<" iteration "<<n_iter<<" temp = "<<Temp<<" : # changements = "<<nch<<"\n";
	if (affich) {
		cout<<"\n******************************************\n matrice des correspondances finale :\n";
		for (i=1; i<=c1; i++) {
			for (j=1; j<=c2; j++) cout<<" "<<Tcorresp[i][j];
			cout<<"\n";
			}
		cout<<" # de classes 'communes' = "<<cc<<"\n";
		cout<<" % de pixels de 'meme' label = "<<100.*(U2maxim-beta*cc)/nbrlig/nbrcol<<" %\n";
		}
	for (i=0; i<=c1; i++) 
		if (H2d[i]!=NULL) delete[] H2d[i];
	if (H2d!=NULL) delete[] H2d;
	for (i=0; i<=c1; i++)
		if (Tcorres2[i]!=NULL) delete[] Tcorres2[i];
	if (Tcorres2!=NULL) delete[] Tcorres2;
	if (T_E!=NULL) delete[] T_E;
	}
 
void imalabels::sauve_Ima(char *nomfich, int icanal) const {
	ofstream sortie;
	sortie.open(nomfich,ios::out|ios::binary);
	if (!sortie) {
		cout<<" ouverture de "<<nomfich<<" impossible\n";
		exit (-1);
		}
	BYTE* adval=(BYTE*)Tsites[0];
	int itype, sizeval;
	itype=0;
	sizeval=sizeof(BYTE);
	cout<<" image "<<nblig<<" lig.& "<<nbcol<<" col., de type ";
	cout<<" bytes => itype = "<<itype<<"\n";
	sortie.write((char *)&nblig,sizeof(int));
	sortie.write((char *)&nbcol,sizeof(int));
	sortie.write((char *)&itype,sizeof(int));
	for (long int i=0; i<nbpix; i++) {
		adval=(BYTE*)Tsites[i]+icanal;
		sortie.write((char *)adval,sizeval);
		}
	sortie.close();
	}

void imalabels::sauve3d_Ima(char *nomfich, int icanal) const {
	int i,j,k;
	imadata<BYTE> I3c(nblig,nbcol,3);
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) {
			I3c(i,j,0)=I3c(i,j,1)=I3c(i,j,2)=0;
//			k=(*this)(i,j,icanal);
			k=Tima[((i*nbcol)+j)*nblayer+icanal];
			if (k>0) {
				if (k>=1 && k<=3) I3c(i,j,k-1)=255;
				else {
					if (k>=4 && k<=7) {
						if (k==4) {I3c(i,j,0)=I3c(i,j,1)=255;}
//						if (k==5) {I3c(i,j,0)=I3c(i,j,2)=255;}
						if (k==5) {I3c(i,j,0)=I3c(i,j,2)=155;}
						if (k==6) {I3c(i,j,1)=I3c(i,j,2)=255;}
						if (k==7) {I3c(i,j,0)=I3c(i,j,1)=I3c(i,j,2)=255;}
						}
					else cout<<" classe "<<k<<" => rajouter couleur\n";
					}
				}
			}
//	I3c.sauve_ImaBSQ(nomfich);
	I3c.sauve_ImaBIP(nomfich);
	}

void imalabels::sauve1d_Ima(char *nomfich, int icanal) const {
	int i,j,k,labmax=0;
	imadata<BYTE> I1c(nblig,nbcol,1);
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) {
//			k=(*this)(i,j,icanal);
			k=Tima[((i*nbcol)+j)*nblayer+icanal];
			if (k>labmax) labmax=k;
		}
	cout<<" classe de 1 a "<<labmax<<"\n";
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) {
			I1c(i,j,0)=0;
//			k=(*this)(i,j,icanal);
			k=Tima[((i*nbcol)+j)*nblayer+icanal];
			if (k>0) {
				k=k*255/labmax;
				I1c(i,j,0)=k;
			}
		}
	I1c.sauve_ImaBSQ(nomfich);
	}

imalabels::imalabels(imadata<float> &ima, statclasses Tclas, char cline, float bPotts, float a0, float a1, float a2)
                     : imasites(ima.nlig(),ima.ncol()) {
	if (cline=='0') {*this=imalabels(ima,Tclas,bPotts);}
	else {
		Tima=NULL;
		const int nclas=Tclas.nclas(), itermax=1000, ncasLP=16;
		const float rTemp=0.99f;
		int i,j,icl,nch,ii,jj,i0,i2,j0,j2,nlinesV,nlinesH,il0,nchV,nchH;
		unsigned int iter=0;
		int *Nvs=new int[nclas+1];
		float temp=100.f;
		double *Papriori=new double[nclas+1], xx, xcumul, xmin;
		int *jj_Uiimin=new int[nclas+1];
		long int l,ipixl;
		nblayer=4;
		Tima=new BYTE[nbpix*nblayer];
		pixel<float> pix, dist, dist0;
		for (l=0; l<nbpix; l++) {
			for (j=0; j<nblayer; j++) Tima[l*nblayer+j]=0;
			Tsites[l]=&(Tima[l*nblayer]);
		}
		imadata<float> imaU0(nblig,nbcol,nclas+1);
		for (i=0; i<nblig; i++) {         // calcul du terme d'attache aux donnees
			ipixl=i*nbcol;                                // stocke sur la couche 0
			for (j=0; j<nbcol; j++) {
				pix=ima.pix(i,j);
				dist0=Tclas.distclasses(pix);	dist=Tclas.U0gaus_classes(pix);
				icl=(int)dist[0];
				Tima[(ipixl+j)*nblayer]=(BYTE)icl; Tima[(ipixl+j)*nblayer+1]=(BYTE)icl;
				imaU0.pix(i,j,dist);
			}
		}
		nlinesV=nlinesH=0;
		for (i=0; i<nblig-1; i++) {        // initialisation des champs des lignes
			ipixl=i*nbcol;               // stockes sur les couches 2 (VL) et 3 (HL)
			for (j=0; j<nbcol-1; j++) {
				icl=Tima[(ipixl+j)*nblayer];
				if (icl==Tima[(ipixl+j+1)*nblayer]) Tima[(ipixl+j)*nblayer+2]=0;
				else {Tima[(ipixl+j)*nblayer+2]=1; nlinesV++;}
				if (icl==Tima[(ipixl+nbcol+j)*nblayer]) Tima[(ipixl+j)*nblayer+3]=0;
				else {Tima[(ipixl+j)*nblayer+3]=1; nlinesH++;}
//				Tima[(ipixl+j)*nblayer+2]=rand()%2; Tima[(ipixl+j)*nblayer+3]=rand()%2;
				Tima[(ipixl+j)*nblayer+2]=Tima[(ipixl+j)*nblayer+3]=0;
			}
		}
		cout<<" initialement "<<nlinesH<<" lig.H & "<<nlinesV<<" lig.V\n";
		char *algoG="glob";
//		char *algoG="sequ";
		if (strcmp(algoG,"glob")==0) {
			temp=1;
			do {
				iter=iter+1; temp=temp*rTemp; nch=0;
				for (i=0; i<nblig; i++) {
					ipixl=i*nbcol;
					for (j=0; j<nbcol; j++) {
						icl=(int)Tima[(i*nbcol+j)*nblayer+1];
						for (ii=0; ii<=nclas; ii++) Papriori[ii]=imaU0(i,j,ii);
						for (ii=1; ii<=nclas; ii++) {
							Tima[(i*nbcol+j)*nblayer+1]=ii;	jj_Uiimin[ii]=0;	xmin=1.e+9;
							for (jj=0; jj<ncasLP; jj++) {
								il0=jj; Tima[(i*nbcol+j)*nblayer+3]=(BYTE)(il0%2);
								il0=il0/2; Tima[(i*nbcol+j)*nblayer+2]=(BYTE)(il0%2);
								il0=il0/2; if (i>0) Tima[((i-1)*nbcol+j)*nblayer+3]=(BYTE)(il0%2);
								il0=il0/2; if (j>0) Tima[(i*nbcol+j-1)*nblayer+2]=(BYTE)(il0%2);
								xx=0.; i0=0; i2=0; j0=0; j2=0;
								if (i>0 && j>0) {i0=(int)Tima[((i-1)*nbcol+j-1)*nblayer+3]; j0=(int)Tima[((i-1)*nbcol+j-1)*nblayer+2];}
								if (i>0) i2=(int)Tima[((i-1)*nbcol+j)*nblayer+3];
								if (j>0) j2=(int)Tima[(i*nbcol+j-1)*nblayer+2];
								xx+=U_lineprocess (i0,i2,j0,j2,a0,a1,a2);
								i0=0; i2=0; j0=0; j2=0;
								if (i>0) i0=(int)Tima[((i-1)*nbcol+j)*nblayer+3];
								if (i>0 && j<nbcol-2) i2=(int)Tima[((i-1)*nbcol+j+1)*nblayer+3];
								if (i>0) j0=(int)Tima[((i-1)*nbcol+j)*nblayer+2];
								j2=(int)Tima[(i*nbcol+j)*nblayer+2];
								xx+=U_lineprocess (i0,i2,j0,j2,a0,a1,a2);
								i0=0; i2=0; j0=0; j2=0;
								i0=(int)Tima[(i*nbcol+j)*nblayer+3];
								if (j<nbcol-2) i2=(int)Tima[(i*nbcol+j+1)*nblayer+3];
								j0=(int)Tima[(i*nbcol+j)*nblayer+2];
								if (i<nblig-2) j2=(int)Tima[((i+1)*nbcol+j)*nblayer+2];
								xx+=U_lineprocess (i0,i2,j0,j2,a0,a1,a2);
								i0=0; i2=0; j0=0; j2=0;
								if (j>0) i0=(int)Tima[(i*nbcol+j-1)*nblayer+3];
								i2=(int)Tima[(i*nbcol+j)*nblayer+3];
								if (j>0) j0=(int)Tima[(i*nbcol+j-1)*nblayer+2];
								if (j>0 && i<nblig-2) j2=(int)Tima[((i+1)*nbcol+j-1)*nblayer+2];
								xx+=U_lineprocess (i0,i2,j0,j2,a0,a1,a2);
								cmptvs4connexLP (i,j,nclas,Nvs);
								xx+=bPotts*(Nvs[0]-2*Nvs[ii]);
								if (xx<xmin) {
									jj_Uiimin[ii]=jj;
									xmin=xx;
								}
							}
							Papriori[ii]+=xmin;
						}
						Tima[(ipixl+j)*nblayer+1]=icl;
						xmin=1.e+9;                                   // algo="ICM
						icl=0;
						for (ii=1; ii<=nclas; ii++)
							if (Papriori[ii]<xmin) {icl=ii; xmin=Papriori[ii];}
						if (icl<=0 || icl>nclas) cout<<"Pb icl "<<ii<<" : "<<icl<<" "<<nclas<<"\n";
						else {
							if ((int)Tima[(ipixl+j)*nblayer+1]!=icl) nch++;
							Tima[(ipixl+j)*nblayer+1]=(BYTE)icl;
							il0=jj_Uiimin[icl];
							Tima[(ipixl+j)*nblayer+3]=il0%2;
							il0=il0/2;
							Tima[(ipixl+j)*nblayer+2]=il0%2;
							il0=il0/2;
							if (i>0) Tima[(ipixl-nbcol+j)*nblayer+3]=il0%2;
							il0=il0/2;
							if (j>0) Tima[(ipixl+j-1)*nblayer+2]=il0%2;
						}
					}
				}
				cout<<" iter. "<<setw(3)<<iter<<" temp. "<<setw(6)<<setprecision(3)<<temp;
				cout<<" => #changements = "<<setw(6)<<nch<<"\n";
				nlinesV=0; nlinesH=0;
				for (l=0; l<nbpix; l++) {nlinesV+=Tima[l*nblayer+2]; nlinesH+=Tima[l*nblayer+3];}
				cout<<" #linesV = "<<nlinesV<<", #linesH = "<<nlinesH<<"\n";
			} while (iter<itermax && nch>0 || temp>1);
		}
		if (strcmp(algoG,"sequ")==0) {
			char *algo="SAM";
			do {
				iter=iter+1;
				temp=temp*rTemp;
				nch=0; nchV=0; nchH=0; nlinesV=0; nlinesH=0;
				for (i=0; i<nblig; i++) {
					ipixl=i*nbcol;
					for (j=0; j<nbcol; j++) {
						icl=(int)Tima[(i*nbcol+j)*nblayer+1];
						for (ii=0; ii<2; ii++) Papriori[ii]=0.;
						il0=Tima[(ipixl+j)*nblayer+2];
						for (ii=0; ii<2; ii++) {                  // lignes verticales
							Tima[(ipixl+j)*nblayer+2]=ii;
							if (i>0 || j<nbcol-2) {
								i0=Tima[(ipixl-nbcol+j)*nblayer+2];
								j0=Tima[(ipixl-nbcol+j)*nblayer+3];
								j2=Tima[(ipixl-nbcol+j+1)*nblayer+3];
								Papriori[ii]+=U_lineprocess (ii,i0,j0,j2,a0,a1,a2);
							}
							if (i<nblig-2 || j<nbcol-2) {
								i2=Tima[(ipixl+nbcol+j)*nblayer+2];
								j0=Tima[(ipixl+j)*nblayer+3];
								j2=Tima[(ipixl+j+1)*nblayer+3];
								Papriori[ii]+=U_lineprocess (ii,i2,j0,j2,a0,a1,a2);
							}
							cmptvs4connexLP (i,j,nclas,Nvs);
							Papriori[ii]+=bPotts*(Nvs[0]-Nvs[icl])-bPotts*Nvs[icl];
						}
						xx=Papriori[1-il0]-Papriori[il0];
						xcumul=xx+1.;
						if (xx>=0) {
							xcumul=exp(-xx/temp);
							xx=(double)rand()/RAND_MAX;
						}
						if (xx<xcumul) {
							nchV++;
							Tima[(ipixl+j)*nblayer+2]=1-il0;
						}
						else Tima[(ipixl+j)*nblayer+2]=il0;
						if (Tima[(ipixl+j)*nblayer+2]==1) nlinesV++;
						for (jj=0; jj<2; jj++) Papriori[jj]=0.;
						il0=Tima[(ipixl+j)*nblayer+3];
						for (jj=0; jj<2; jj++) {                  // lignes horizontales
							Tima[(ipixl+j)*nblayer+3]=jj;
							if (j>0 || i<nblig-2) {
								j0=Tima[(ipixl+j-1)*nblayer+3];
								i0=Tima[(ipixl+j-1)*nblayer+2];
								i2=Tima[(ipixl+nbcol+j-1)*nblayer+2];
								Papriori[jj]+=U_lineprocess (jj,j0,i0,i2,a0,a1,a2);
							}
							if (i<nblig-2 || j<nbcol-2) {
								j2=Tima[(ipixl+j+1)*nblayer+3];
								i0=Tima[(ipixl+j)*nblayer+2];
								i2=Tima[(ipixl+nbcol+j)*nblayer+2];
								Papriori[jj]+=U_lineprocess (jj,j2,i0,i2,a0,a1,a2);
							}
							cmptvs4connexLP (i,j,nclas,Nvs);
							Papriori[jj]+=bPotts*(Nvs[0]-Nvs[icl])-bPotts*Nvs[icl];
						}
						xx=Papriori[1-il0]-Papriori[il0];
						xcumul=xx+1.;
						if (xx>=0) {
							xcumul=exp(-xx/temp);
							xx=(double)rand()/RAND_MAX;
						}
						if (xx<xcumul) {
							nchH++;
							Tima[(ipixl+j)*nblayer+3]=1-il0;
						}
						else Tima[(ipixl+j)*nblayer+3]=il0;
						if (Tima[(ipixl+j)*nblayer+3]==1) nlinesH++;
						cmptvs4connexLP (i,j,nclas,Nvs);
						for (ii=0; ii<=nclas; ii++) Papriori[ii]=0.;
						if (algo=="SAM") {
							icl=(int)Tima[(i*nbcol+j)*nblayer+1];
							do {
								ii=(rand()%nclas)+1;
							} while (ii==icl);
							Papriori[icl]=imaU0(i,j,icl)+bPotts*(Nvs[0]-2*Nvs[icl]);
							Papriori[ii]=imaU0(i,j,ii)+bPotts*(Nvs[0]-2*Nvs[ii]);
							xx=Papriori[ii]-Papriori[icl];
							if (xx<0) icl=ii;
							else {
								xcumul=exp(-xx/temp);
								xx=(double)rand()/RAND_MAX;
								if (xx<xcumul) icl=ii;
							}
						}
						else {                                   // non teste
							for (ii=1; ii<=nclas; ii++)
								Papriori[ii]=imaU0(i,j,ii)+bPotts*(Nvs[0]-2*Nvs[ii]);
							for (ii=1; ii<=nclas; ii++) Papriori[ii]=Papriori[ii]/temp;
							xx=0.;
							for (ii=1; ii<=nclas; ii++) xx+=exp(-Papriori[ii]);
							if (xx!=0. && algo=="SAG") {
								for (ii=1; ii<=nclas; ii++) Papriori[ii]=exp(-Papriori[ii])/xx;
								xx=(double)rand()/RAND_MAX;
								xcumul=0.;
								icl=0;
								do {
									icl++;
									xcumul+=Papriori[icl];
								} while (xcumul<xx && icl<nclas);
								if (xcumul<xx) icl=(int)Tima[(ipixl+j)*nblayer+1];
							}
							else {                                //notamment algo="ICM
								xmin=1.e+9;
								icl=0;
								for (ii=1; ii<=nclas; ii++) {
									if (Papriori[ii]<xmin) {
										icl=ii;
										xmin=Papriori[ii];
									}
								}
							}
						}
						if (icl<=0 || icl>nclas)
							cout<<"Pb icl "<<ii<<" "<<jj<<" : "<<icl<<" "<<nclas<<"\n";
						if ((int)Tima[(ipixl+j)*nblayer+1] != icl) {
							Tima[(ipixl+j)*nblayer+1]=(BYTE)icl;
							nch++;
						}
					}
				}
				cout<<" iter. "<<setw(3)<<iter<<" temp. "<<setw(6)<<setprecision(3)<<temp;
				cout<<" => #changements = "<<setw(6)<<nch<<"\n";
				cout<<" #linesV = "<<nlinesV<<" "<<nchV<<", #linesH = "<<nlinesH<<" "<<nchH<<"\n";
			} while (iter<itermax && nch>0 || temp>1);
		}
		cout<<" fin de la convergence pour "<<iter<<" iterations, ";
		cout<<nch<<" pixels encore changes\n";
		if (Nvs!=NULL) delete[] Nvs;
		if (Papriori!=NULL) delete[] Papriori;
		if (jj_Uiimin!=NULL) delete[] jj_Uiimin;
		}
	}

imalabels::imalabels(int nl, int nc, int nclas, float bPotts, float Bdat) : imasites(nl,nc) {
	Tima=NULL;
	const int itermax=10000;
	const float temp0=nclas*bPotts, rtemp=0.9f;/*0.99f/*, alpha=0.2f*/;
	int i,j,k,nbclas,icl,iter,nch,ii, *Nvs=new int[nclas+1];
	long int l,ipixl;
	float temp;
	double *Papriori=new double[nclas+1], xx, xcumul, xmin;
	nblayer=nclas;
	Tima=new BYTE[nbpix*nblayer];
	for (l=0; l<nbpix; l++) {
		for (j=0; j<nblayer; j++) Tima[l*nblayer+j]=0;
		Tsites[l]=&(Tima[l*nblayer]);
	}
	unsigned int    number;
  errno_t         err;
	for (k=1; k<nclas; k++) {
		nbclas=k+1;
		for (i=0; i<nblig; i++) {          // simulation du terme d'attache aux donnees
			ipixl=i*nbcol;                                // stocke sur la couche 0
			for (j=0; j<nbcol; j++) {
//				icl=(rand()%nbclas)+1;

        err = rand_s(&number); if (err!= 0) cout<<"The rand_s function failed!\n";
        icl=(int)((double)number/(double)UINT_MAX*nclas+1);
				icl=maxi(mini(icl,nclas+1),0);
				Tima[(ipixl+j)*nblayer]=(BYTE)icl;
				Tima[(ipixl+j)*nblayer+k]=(BYTE)icl;
			}
		}
		iter=0; temp=temp0;
		do {
			iter=iter+1;
			temp=temp*rtemp;
			nch=0;
			for (i=0; i<nblig; i++) {
				ipixl=i*nbcol;
				for (j=0; j<nbcol; j++) {
					for (ii=0; ii<=nbclas; ii++) Papriori[ii]=0.;
					cmptvs8connex (i,j,k,nbclas,Nvs);
//					cmptvs4connex (i,j,k,nbclas,Nvs);
					for (ii=1; ii<=nbclas; ii++) Papriori[ii]=Bdat+bPotts*(Nvs[0]-Nvs[ii]);
//					for (ii=1; ii<=nbclas; ii++) {
//						Papriori[ii]=Bdat;
//						for (int jj=1; jj<=nbclas; jj++) {
//							if (ii!=jj) Papriori[ii]+=bPotts*Nvs[jj];
//							if (ii!=jj) Papriori[ii]+=bPotts*(1+alpha*(2*abs(ii-jj)/(nbclas-1)-1))*Nvs[jj];
//						}
//					}
					Papriori[(int)Tima[(i*nbcol+j)*nblayer+k]]-=Bdat;
					for (ii=1; ii<=nbclas; ii++) Papriori[ii]=Papriori[ii]/temp;
					xx=0.;
					for (ii=1; ii<=nbclas; ii++) xx+=exp(-Papriori[ii]);
					if (xx==0.) {
						xmin=1.e+9;
						icl=0;
						for (ii=1; ii<=nbclas; ii++) {
							if (Papriori[ii]<xmin) {
								icl=ii;
								xmin=Papriori[ii];
							}
						}
					}
					else {
						for (ii=1; ii<=nbclas; ii++) Papriori[ii]=exp(-Papriori[ii])/xx;
//						xx=(double)rand()/RAND_MAX;
						err=rand_s(&number);
						if (err!=0) cout<<"The rand_s function failed!\n";
						xx=(double)number/(double)UINT_MAX*1.;
						xcumul=0.;
						icl=0;
						do {
							icl++;
							xcumul+=Papriori[icl];
						} while (xcumul<xx && icl<nbclas);
						if (xcumul<xx) {
							icl=(int)Tima[(ipixl+j)*nblayer+k];
//							cout<<" un coup pour rien "<<xcumul<<" < "<<xx<<"\n";
						}
					}
					if (icl<=0 || icl>nbclas)
						cout<<"Pb icl "<<i<<" "<<j<<" : "<<icl<<" "<<nbclas<<"\n";
					if ((int)Tima[(ipixl+j)*nblayer+k] != icl) {
						Tima[(ipixl+j)*nblayer+k]=(BYTE)icl;
						nch++;
					}
				}
			}
			if (iter<10 || iter%10==0) cout<<" iter. "<<iter<<" #changements = "<<nch<<"\n";
		} while (iter<itermax && nch>0);
		cout<<" fin de la convergence pour "<<iter<<" iterations, "<<nch<<" pixels encore changes\n";
	} cout<<" fin du constructeur par simulation d'une realisation MRF\n";
	if (Nvs!=NULL) delete[] Nvs;
  if (Papriori!=NULL) delete[] Papriori;
}

imalabels::imalabels(imadata<float> &ima, int nbcl, char *algo) : imasites(ima.nlig(),ima.ncol()) {
	Tima=NULL;
	const unsigned int nbsimG=10*nbcol*nbcl, msimG=nbsimG/4;
	const float rTfT0=0.02f, rTemp=exp(log(rTfT0)/nbsimG);
	nblayer=nbsimG+1;
	Tima=new BYTE[nbpix*nblayer];
	int i,j,k,icl,iconfigV,n;
	long int l;
	for (l=0; l<nbpix; l++) {
		for (j=0; j<nblayer; j++) Tima[l*nblayer+j]=0;
		Tsites[l]=&(Tima[l*nblayer]);
		}
	statclasses Tclas;
	int dim=ima.ncanaux();
	stat1class cl(1,dim);
	cl.dcova()=1.;
	for (j=0; j<dim; j++) {
		float x=(float)ima.varI(j);
		cl.cova(j,j)=x;
		cl.icova(j,j)=1.f/x;
		cl.dcova()*=x;
		}
	cl.proba()=1.f/nbcl;
	for (k=0; k<nbcl; k++) {
		cl.label()=k+1;
		for (j=0; j<dim; j++)
			cl.mean(j)=(float)(ima.minI(j)+(0.5+k)*(ima.maxI(j)-ima.minI(j))/nbcl);
		Tclas.ajoute(cl);
		}
	cout<<" initialisation des classes a :\n"; Tclas.affiche();
	const int nconnex=4;
	const int nconfigV=(int)pow(nbcl+1.f,nconnex);
	double *Taij=new double[nbcl*nconfigV];
	for (i=0; i<nconfigV*nbcl; i++) Taij[i]=1.f/nbcl;
	int* T_Nvs=new int[nconfigV*(nbcl+1)];
	for (i=0; i<nconfigV*(nbcl+1); i++) T_Nvs[i]=0;
	for (i=0; i<nconfigV; i++) {
		iconfigV=i;
		j=i*(nbcl+1);
		for (k=nconnex-1; k>=0; k--) {
			icl=iconfigV/(int)pow(nbcl+1.f,k);
			if (icl>0 && icl<=nbcl) {
				T_Nvs[j+icl]++;
				iconfigV-=icl*(int)pow(nbcl+1.f,k);
				}
			else {
				if (icl<0 || icl>nbcl)
					cout<<" Pb ds calcul du voisinage associe a la configuration "<<i<<"\n";
				}
			}
		for (icl=1; icl<=nbcl; icl++) T_Nvs[j]+=T_Nvs[j+icl];
		}
	int **T_Uicl_config=new int*[nbcl*(nconnex+1)];
	for (icl=1; icl<=nbcl; icl++) {
		for (k=0; k<=nconnex; k++) {
			n=0;
			for (i=0; i<nconfigV; i++)
				if (T_Nvs[i*(nbcl+1)+icl]==k) n++;
			j=(icl-1)*(nconnex+1)+k;
			T_Uicl_config[j]=new int[n+1];
			T_Uicl_config[j][0]=n;
			n=0;
			for (i=0; i<nconfigV; i++)
				if (T_Nvs[i*(nbcl+1)+icl]==k) T_Uicl_config[j][++n]=i;
			}
		}
/*	for (icl=1; icl<=nbcl; icl++)
		for (k=0; k<=nconnex; k++) {
			j=(icl-1)*(nconnex+1)+k;
			n=T_Uicl_config[j][0];
			for (i=1; i<=n; i++) {
				iconfigV=T_Uicl_config[j][i];
				cout<<" icl."<<icl<<" "<<k<<"-k : "<<iconfigV<<" => #(lv=icl) = ";
				cout<<T_Nvs[iconfigV*(nbcl+1)+icl]<<"\n";
				}
			}
	char aa; cin>>aa;*/
	unsigned int iterEM=0, nchEM=0, iterG;
	const unsigned int itermax=100, nchmax=nbpix*0;
	const float Temp0=1.0f;
	do {
		iterEM++; iterG=0;
// tirage des echantillons a partir de l'echatillonneur de Gibbs
		pixel<float> pix, dist;
		imadata<float> imaU0(nblig,nbcol,nbcl+1);
		int ipixl, icl0, nchG=0, ii, iconfigV;
		for (i=0; i<nblig; i++) {               // calcul du terme d'attache aux donnees
			ipixl=i*nbcol;                                   // stocke sur la couche 0
			for (j=0; j<nbcol; j++) {
				pix=ima.pix(i,j);
				dist=Tclas.U0gaus_classes(pix);
				icl=(int)dist[0];
				for (k=0; k<(int)nbsimG; k++) Tima[(ipixl+j)*nblayer+k]=(BYTE)icl;
				imaU0.pix(i,j,dist);
				}
			}
		double *Papriori=new double[nbcl+1];
		int *Nvs=new int[nbcl+1];
		double xx, sxx, xcumul, bPotts=0.;
		float Temp=Temp0;
		do {
			iterG++;
			for (i=0; i<nblig; i++) {
				ipixl=i*nbcol;
				for (j=0; j<nbcol; j++) {
					icl0=(int)Tima[(ipixl+j)*nblayer+iterG-1];
					cmptvs4connex (i,j,iterG-1,nbcl,Nvs);
					for (ii=1; ii<=nbcl; ii++)
						Papriori[ii]=exp(-(imaU0(i,j,ii)+bPotts*(Nvs[0]+1-Nvs[ii]))/Temp);
					sxx=0.;
					for (ii=1; ii<=nbcl; ii++) sxx+=Papriori[ii];
					if (sxx!=0) {
						for (ii=1; ii<=nbcl; ii++) Papriori[ii]/=sxx;
						xx=(double)rand()/RAND_MAX;
						xcumul=0.; icl=0;
						do {
							icl++;
							xcumul+=Papriori[icl];
							} while (xcumul<xx && icl<nbcl);
						if (xcumul<xx) icl=icl0;
						}
					else {
						cout<<" dans MPM, sxx="<<sxx<<" en "<<i<<" "<<j<<"\n";
						for (ii=1; ii<=nbcl; ii++) {
							cout<<Papriori[ii]<<" "<<exp(-imaU0(i,j,ii))<<" "<<imaU0(i,j,ii);
							cout<<"\n";
						}
						char bb; cin>>bb;
						}
					if (icl<1 || icl>nbcl)
						cout<<"Pb icl "<<i<<" "<<j<<" : "<<icl<<" "<<nbcl<<"\n";
					Tima[(ipixl+j)*nblayer+iterG]=(BYTE)icl;
					if (icl0 != icl) nchG=nchG+1;
					}
				}
			if (iterG%(nbsimG/10)==0) {
				cout<<" iter. "<<setw(3)<<iterG-nbsimG/10<<" a "<<setw(3)<<iterG;
				cout<<" => #changements = "<<setw(6)<<nchG<<"\n";
				nchG=0;
				}
			Temp*=rTemp;
/*			for (i=0; i<nblig; i++) {
				for (j=0; j<nbcol; j++) cout<<" "<<(int)Tima[(i*nbcol+j)*nblayer+iterG];
				cout<<"\n";
				}*/
			} while (iterG<nbsimG-1);
// approximation des esperances
		double *Pxs=new double[nbpix*nbcl];
		for (icl=1; icl<=nbcl; icl++) {
			Papriori[icl]=0;
			for (ii=0; ii<nconfigV; ii++) Taij[(icl-1)*nconfigV+ii]=0;
			for (l=0; l<nbpix; l++) Pxs[l*nbcl+icl-1]=0;
			}
		sxx=0;
		for (iterG=msimG; iterG<nbsimG; iterG++) 
			for (i=0; i<nblig; i++)
				for (j=0; j<nbcol; j++) {
					icl=(int)Tima[(i*nbcol+j)*nblayer+iterG];
					Pxs[(i*nbcol+j)*nbcl+icl-1]+=1.;
					Papriori[icl]+=1.;
					iconfigV=0;
					if (i>0) {
						icl=(int)Tima[((i-1)*nbcol+j)*nblayer+iterG];
						iconfigV+=icl*(int)pow(nbcl+1.f,3);
						}
					if (j>0) {
						icl=(int)Tima[(i*nbcol+j-1)*nblayer+iterG];
						iconfigV+=icl*(int)pow(nbcl+1.f,2);
						}
					if (i<nblig-1) {
						icl=(int)Tima[((i+1)*nbcol+j)*nblayer+iterG];
						iconfigV+=icl*(nbcl+1);
						}
					if (j<nbcol-1) {
						icl=(int)Tima[(i*nbcol+j+1)*nblayer+iterG];
						iconfigV+=icl;
						}
					Taij[(icl-1)*nconfigV+iconfigV]+=1;
					sxx+=1;
/*					cout<<" voisinage : ";
					if (i>0) cout<<(int)Tima[((i-1)*nbcol+j)*nblayer+iterG]<<" ";
					if (j>0) cout<<(int)Tima[(i*nbcol+j-1)*nblayer+iterG]<<" ";
					if (i<nblig-1) cout<<(int)Tima[((i+1)*nbcol+j)*nblayer+iterG]<<" ";
					if (j<nbcol-1) cout<<(int)Tima[(i*nbcol+j+1)*nblayer+iterG]<<"\n";
					for (icl=0; icl<=nbcl; icl++) Nvs[icl]=0;
					cout<<" config."<<iconfigV<<" => ";
					for (k=nconnex-1; k>=0; k--) {
						icl=iconfigV/(int)pow(nbcl+1,k);//cout<<" cl."<<icl;
						if (icl>0 && icl<=nbcl) {
							Nvs[icl]++;
							iconfigV-=icl*(int)pow(nbcl+1,k);
							}
						}
					for (icl=1; icl<=nbcl; icl++) cout<<" cl."<<icl<<" #= "<<Nvs[icl]; cout<<"\n";
					char bb; cin>>bb;*/
					}
		for (icl=1; icl<=nbcl; icl++) {
			Papriori[icl]/=sxx;
			for (l=0; l<nbpix; l++) Pxs[l*nbcl+icl-1]/=(nbsimG-msimG);
			}
/*		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				cout<<" "<<i<<" "<<j<<" => Pxs : ";
				for (icl=1; icl<=nbcl; icl++) cout<<" "<<Pxs[(i*nbcol+j)*nbcl+icl-1];
				cout<<"\n";
				}*/
		for (ii=0; ii<nconfigV; ii++) {
			sxx=0;
			for (icl=1; icl<=nbcl; icl++) sxx+=Taij[(icl-1)*nconfigV+ii];
			if (sxx>0)
				for (icl=1; icl<=nbcl; icl++) Taij[(icl-1)*nconfigV+ii]/=sxx;
//			else 
//				for (icl=1; icl<=nbcl; icl++) Taij[(icl-1)*nconfigV+ii]=1./nbcl;
			}
// remise a jour des parametres
		int n_bPotts=0;
		bPotts=0;
		for (icl=1; icl<=nbcl; icl++) {
			for (k=0; k<=nconnex; k++) {
				j=(icl-1)*(nconnex+1)+k;
				sxx=0;
				n=T_Uicl_config[j][0];
				if (n>0) {
					for (ii=1; ii<=n; ii++) {
						iconfigV=T_Uicl_config[j][ii];
						sxx+=Taij[(icl-1)*nconfigV+iconfigV];
						}
					sxx/=n;
					}
				if (sxx>0) {
					j=iconfigV*(nbcl+1);
					bPotts-=log(sxx)/(T_Nvs[j]+1-T_Nvs[j+icl]);
					n_bPotts++;
					}
				}
			}
/*		bool inul;cout<<" calcul de beta\n";
		for (ii=0; ii<nconfigV; ii++) {
			inul=1;
			for (icl=1; icl<=nbcl; icl++) 
				if (Taij[(icl-1)*nconfigV+ii]>0) inul=0;
			if (!inul) {
				for (icl=0; icl<=nbcl; icl++) Nvs[icl]=T_Nvs[ii*(nbcl+1)+icl];
				for (icl=1; icl<=nbcl; icl++) {
					xx=Taij[(icl-1)*nconfigV+ii];
					cout<<" icl "<<icl<<" "<<xx<<" ";
					if (xx>0) {cout<<-log(xx)<<" "<<(Nvs[0]+1-Nvs[icl]);
						bPotts=-log(xx)/(Nvs[0]+1-Nvs[icl]);
						n_bPotts++;
						}
					}
				cout<<"\n";
				}
			}*/
		if (n_bPotts>0) bPotts/=n_bPotts;
		else cout<<" # de configurations d'estimation de bPotts nul ???? "<<n_bPotts<<"\n";
		cout<<" estimation de bPotts : "<<bPotts<<"\n";
		for (icl=1; icl<=nbcl; icl++) {
			cl.label()=icl;
			sxx=0;
			for (l=0; l<nbpix; l++) sxx+=Pxs[l*nbcl+icl-1];
			if (sxx>0 && sxx<=nbpix) cl.proba()=(float)(sxx/nbpix);
			else {
				cout<<" Pb ds calcul de la proba a priori de "<<icl<<" : "<<sxx/nbpix<<"\n";
				for (l=0; l<nbpix; l++) cout<<" "<<Pxs[l*nbcl+icl-1];cout<<"\n";
				for (i=0; i<nblig; i++) {
					for (j=0; j<nbcol; j++) cout<<" "<<icl;
					cout<<"\n";
					}
				char bb; cin>>bb;
				}
			cl.dcova()=1.;
			for (k=0; k<ima.ncanaux(); k++) {
				sxx=0;
				for (i=0; i<nblig; i++) 
					for (j=0; j<nbcol; j++)
						sxx+=ima(i,j,k)*Pxs[(i*nbcol+j)*nbcl+icl-1];
				sxx=sxx/nbpix/cl.proba();
				cl.mean(k)=(float)sxx;
				xx=sxx;
				sxx=0;
				for (i=0; i<nblig; i++) 
					for (j=0; j<nbcol; j++)
						sxx+=pow(ima(i,j,k)-xx,2)*Pxs[(i*nbcol+j)*nbcl+icl-1];
				sxx=sxx/(nbpix-1)/cl.proba();
				cl.cova(k,k)=(float)sxx;
				cl.dcova()*=(float)sxx;
				}
			Tclas.modif(cl);
			}
		Tclas.affiche();
		if (Pxs!=NULL) delete[] Pxs;
		if (Nvs!=NULL) delete[] Nvs;
		if (Papriori!=NULL) delete[] Papriori;
// critere d'arret
		nchEM=0;
		for (i=0; i<nblig; i++) {
			ipixl=i*nbcol;
			for (j=0; j<nbcol; j++) {
				ii=(ipixl+j)*nblayer+nbsimG;
				if (Tima[ii-1]!=Tima[ii]) {
					nchEM++;
					Tima[ii]=Tima[ii-1];
					}
				}
			}
		cout<<" iteration "<<iterEM<<" de l'EM Gibbsien => # de pixels changes = "<<nchEM<<"\n";
//		char aa; cin>>aa;
		} while ((nchEM>nchmax || iterEM<20) && iterEM<itermax);
	if (Taij!=NULL) delete[] Taij;
	if (T_Nvs!=NULL) delete[] T_Nvs;
	for (j=0; j<nbcl*(nconnex+1); j++) 
		if (T_Uicl_config[j]!=NULL) delete[] T_Uicl_config[j];
	if (T_Uicl_config!=NULL) delete[] T_Uicl_config;
	}

int imalabels::HorizLength (int i, int j, int icl_m, int &j_debut, int k) {
	int icl, iclij, Maxj=j,Minj=j, horiz_length=-1;
	if (i>=0 && i<nblig && j>=0 && j<nbcol) {
		iclij=(int)Tima[(i*nbcol+j)*nblayer+k];
		if (iclij==icl_m) {
			do icl=(int)Tima[(i*nbcol+mini(++Maxj,nbcol-1))*nblayer+k];
			while (icl==icl_m && Maxj<nbcol);		
			do icl=(int)Tima[(i*nbcol+maxi(--Minj,0))*nblayer+k];
			while (icl==icl_m && Minj>=0);
			j_debut=++Minj;
			horiz_length=Maxj-Minj;
		} else {
			do icl=(int)Tima[(i*nbcol+mini(++Maxj,nbcol-1))*nblayer+k];
			while (icl==icl_m && Maxj<nbcol);
			Maxj-=(j+1);
			do icl=(int)Tima[(i*nbcol+maxi(--Minj,0))*nblayer+k];
			while (icl==icl_m && Minj>=0);
			Minj=j-1-Minj;
			if (Maxj>Minj) {horiz_length=Maxj; j_debut=j+1;}
			else {horiz_length=Minj; j_debut=j-1-Minj;}
		}
	}
	return horiz_length;
}

imalabels::imalabels(imadata<float> & ima, double s, double delta_s, float bPotts, 
                     /*char *algo, */float pct_ch) : imasites(ima.nlig(),ima.ncol()) {
	Tima=NULL;
	const int nclas=2, icl_m=1, icl_b=2, itermax=50, nit_raznch=1, nchmax=(int)(pct_ch*nbpix);
//	const float rTemp=0.998f;
	const double alpha=10./delta_s;
	float /*temp=10.f,*/nch=0.f;
	int i,j,icl,ii,jj,j_d,lg0,lg2,lg1,lenghtV,sens/*, *Nvs=new int[nclas+1]*/,iter=0;
	long int l,ipixl;
	double *Papriori=new double[nclas+1], xx/*, xcumul*/, xmin;
//	char *algoICM="ICM", *algoSAG="SAG", *algoSAM="SAM";
	bool i_connex_OK;
	cout<<" classification MRF "<<nclas<<" classes, favorise classe "<<icl_m<<" de largeur (en colonne) reguliere, parametre a priori beta "<<bPotts/*<<", algorithme "<<algo*/<<"\n";
	nblayer=2;
	Tima=new BYTE[nbpix*nblayer];
//	pixel<float> pix, dist;
	for (l=0; l<nbpix; l++) {
		for (j=0; j<nblayer; j++) Tima[l*nblayer+j]=0;
		Tsites[l]=&(Tima[l*nblayer]);
	}
	imadata<float> imaU0(nblig,nbcol,nclas+1);
//	imabin imab(nblig,nbcol);
	for (i=0; i<nblig; i++) {               // calcul du terme d'attache aux donnees
		ipixl=i*nbcol;                                   // stocke sur la couche 0
		for (j=0; j<nbcol; j++) {
			xx=ima(i,j);
			imaU0(i,j,icl_m)=(float)log(1+exp(-alpha*(xx-s)));
			imaU0(i,j,icl_b)=(float)log(1+exp(alpha*(xx-s)));
			if (imaU0(i,j,icl_m)>imaU0(i,j,icl_b)) icl=icl_b;
			else {
				icl=icl_m; 
//				imab(i,j)=1;
			}
			Tima[(ipixl+j)*nblayer]=(BYTE)icl;
			Tima[(ipixl+j)*nblayer+1]=(BYTE)icl;
		}
	}
	sauve1d_Ima("classif_iter0",1);
/*	eltstruct S3(3,3);
	imab=imab.dilate(S3);					// dilatation de la classification aveugle
	for (i=0; i<nblig; i++) {               
		ipixl=i*nbcol;
		for (j=0; j<nbcol; j++) {
			imab(i,j)?icl=icl_m:icl=icl_b;
			Tima[(ipixl+j)*nblayer]=(BYTE)icl;
			Tima[(ipixl+j)*nblayer+1]=(BYTE)icl;
		}
	}*/
	do {
		iter=iter+1;
//		temp=temp*rTemp;
		if (iter%nit_raznch==0) nch=0.f;
		for (i=0; i<nblig; i++) {
			ipixl=i*nbcol;
			for (j=0; j<nbcol; j++) {
				for (ii=0; ii<=nclas; ii++) Papriori[ii]=0.;
				lg0=HorizLength(i-1,j,icl_m,j_d,1);
				lg2=HorizLength(i+1,j,icl_m,j_d,1);
				if (lg0>lg2) {lenghtV=lg0; sens=-1;}
				else
					if (lg2>lg0 || lg2>=0) {lenghtV=lg2; sens=+1;}
					else {lenghtV=0; sens=0;}
				icl=Tima[(i*nbcol+j)*nblayer+1];
				for (ii=1; ii<=nclas; ii++) {
					Tima[(i*nbcol+j)*nblayer+1]=ii;
					lg1=HorizLength(i,j,icl_m,j_d,1);
					i_connex_OK=0; jj=j_d;
					do if (Tima[((i+sens)*nbcol+jj++)*nblayer+1]==icl_m) i_connex_OK=1;
					while (!i_connex_OK && jj<j_d+lg1);
					Papriori[ii]=imaU0(i,j,ii);
					if (lenghtV*i_connex_OK>0)
						Papriori[ii]+=bPotts*abs(lg1-lenghtV);			
					else
						if(ii==icl_m)
							Papriori[ii]+=bPotts*100;  // passer 100 en parametre constant
				}
				Tima[(i*nbcol+j)*nblayer+1]=icl;
//				for (ii=1; ii<=nclas; ii++) Papriori[ii]=Papriori[ii]/temp;
				xx=0.;
				for (ii=1; ii<=nclas; ii++) xx+=exp(-Papriori[ii]);
/*				if (xx!=0. && stricmp(algo,algoSAG)==0) {
					for (ii=1; ii<=nclas; ii++) Papriori[ii]=exp(-Papriori[ii])/xx;
					xx=(double)rand()/RAND_MAX;
					xcumul=0.;
					icl=0;
					do xcumul+=Papriori[++icl];
					while (xcumul<xx && icl<nclas);
					if (xcumul<xx) icl=(int)Tima[(ipixl+j)*nblayer+1];
				}
				else { */                                     //notamment algo="ICM
					xmin=1.e+9;
					icl=0;
					for (ii=1; ii<=nclas; ii++) {
						if (Papriori[ii]<xmin) {
							icl=ii;
							xmin=Papriori[ii];
						}
					}
/*				}*/
				if (icl<=0 || icl>nclas)
					cout<<"Pb icl "<<i<<" "<<j<<" : "<<icl<<" "<<nclas<<"\n";
				if ((int)Tima[(ipixl+j)*nblayer+1] != icl) {
					Tima[(ipixl+j)*nblayer+1]=(BYTE)icl;
					nch=nch+1.f;
//					cout<<" chgt en "<<i<<" "<<j<<" : icl = "<<icl<<"\n";
				}
			}
		}
/*		if (stricmp(algo,algoICM)==0 || iter%100==0) {*/
			cout<<" iter. "<<setw(3)<<iter/*<<" temp. "<<setw(6)<<setprecision(3)<<temp*/;
			cout<<" => #changements = "<<setw(6)<<nch<<"\n"; //char aa; cin>>aa;
/*		}*/
		if (iter==1)    sauve1d_Ima("classif_iter1",1);
		if (iter==10)   sauve1d_Ima("classif_iter10",1);
//		if (iter==100)  sauve1d_Ima("classif_iter100",1);
//		if (iter==1000) sauve1d_Ima("classif_iter1000",1);
//		if (iter==2000) sauve1d_Ima("classif_iter2000",1);
//		if (iter==5000) sauve1d_Ima("classif_iter5000",1);
	} while (iter<itermax && nch>nchmax);
	cout<<" fin de la convergence pour "<<iter<<" iterations, ";
	cout<<nch<<" pixels encore changes\n";
/*	if (Nvs!=NULL) delete[] Nvs;*/
	if (Papriori!=NULL) delete[] Papriori;
}

imalabels::imalabels(int nHRl, int nHRc, imadata<float> & ima_BR, statclasses Tclas, float bPotts, 
                     char * algo) : imasites(ima_BR.nlig()*nHRl,ima_BR.ncol()*nHRc) {
	Tima=NULL;
	const int nclas=Tclas.nclas();
	cout<<" classification de l'image BR en "<<nclas<<" classes\n";
	const int itermax=50;
	const float rTemp=0.998f;
	float temp=10.f,nch=0.f;
	int i,j,icl,ii,iHR,jHR,k,icl0,i2,j2,i2min,i2max,j2min,j2max,i2HR,j2HR;
	long int l,ipixl;
	unsigned int iter=0,it,nchBR;
	int* Nvs=new int[nclas+1];
	double* Papriori=new double[nclas+1];
	double xx;
	nblayer=2;
	Tima=new BYTE[nbpix*nblayer];
	int nk_BR=ima_BR.ncanaux();
	pixel<float> pixBR, simBR(nk_BR), newBR(nk_BR);
	stat1class cl;
	for (l=0; l<nbpix; l++) {
		for (j=0; j<nblayer; j++) Tima[l*nblayer+j]=0;
		Tsites[l]=&(Tima[l*nblayer]);
		}
// initialisation des parametres
	int nl_BR=ima_BR.nlig(), nc_BR=ima_BR.ncol(), nl_HR=nblig, nc_HR=nbcol;
	cout<<" #lig HR/BR = "<<nl_HR<<"/"<<nl_BR<<" = "<<nHRl<<" ?\n";
	cout<<" #col HR/BR = "<<nc_HR<<"/"<<nc_BR<<" = "<<nHRc<<" ?\n";
	if (nl_HR/nl_BR!=nHRl) {
		nHRl=nl_HR/nl_BR;
		cout<<" rectification du ratio des resolutions en  lignes : "<<nHRl<<"\n";
		}
	if (nc_HR/nc_BR!=nHRc) {
		nHRc=nc_HR/nc_BR;
		cout<<" rectification du ratio des resolutions en colonnes : "<<nHRc<<"\n";
		}
// initialisation aleatoire de l'image des labels HR
	for (l=0; l<nbpix; l++) {
		icl=(rand()%nclas)+1;
		for (j=0; j<nblayer; j++) Tima[l*nblayer+j]=icl;
		}
	(*this).sauve3d_Ima("./classifBR_init",0);
// iterations de la classification
	for (ii=1; ii<=nclas; ii++) {
		cl=Tclas.extract(ii);
		cout<<" classe "<<ii<<"\n";
		for (k=0; k<nk_BR; k++) cout<<" "<<cl.mean(k);
		cout<<"\n";
		}
//	{cout<<" on en est ici "; char aa; cin>>aa;}
	do {
		iter=iter+1;
		if (iter%1==0) nch=0.f;
		for (i=0; i<nl_BR; i++) { // pour chaque ligne BR
			for (j=0; j<nc_BR; j++) { // pour chaque colonne BR
				pixBR=ima_BR.pix(i,j); // valeur observee en basse resolution
/*				cout<<" pixel BR "<<i<<" "<<j<<": ";
				for (ii=0; ii<nk_BR; ii++) cout<<pixBR[ii]<<" ";
				cout<<"\n";
				char aa; cin>>aa;*/
// calcul de la valeur BR simulee d'apres l'image des labels et les centres des classes
				for (ii=0; ii<=nclas; ii++) Nvs[ii]=0;
				for (iHR=0; iHR<nHRl; iHR++) {
					ipixl=(i*nHRl+iHR)*nc_HR;
					for (jHR=0; jHR<nHRc; jHR++) {
						icl=Tima[(ipixl+j*nHRc+jHR)*nblayer]; // label considerant uniquement le terme d'attache aux donnees
						Nvs[icl]++;
						Nvs[0]++;
						}
					}
				for (k=0; k<nk_BR; k++) simBR[k]=0;
				for (ii=1; ii<=nclas; ii++) {
					cl=Tclas.extract(ii);
					for (k=0; k<nk_BR; k++) simBR[k]+=Nvs[ii]*cl.mean(k);
					}
				for (k=0; k<nk_BR; k++) simBR[k]/=Nvs[0];
				xx=0.;
				for (k=0; k<nk_BR; k++) xx+=pow(pixBR[k]-simBR[k],2);
				Papriori[0]=xx;
				it=0;
				do {
					nchBR=0;
					for (iHR=0; iHR<nHRl; iHR++) { // pour chaque ligne HR
						ipixl=(i*nHRl+iHR)*nc_HR;
						for (jHR=0; jHR<nHRc; jHR++) { // pour chaque colonne HR
							for (ii=1; ii<=nclas; ii++) Papriori[ii]=Papriori[0];
							icl0=Tima[(ipixl+j*nHRc+jHR)*nblayer];
							do {
								icl=(rand()%nclas)+1;
							} while (icl==icl0);
							for (k=0; k<nk_BR; k++) newBR[k]=simBR[k]*Nvs[0];
							cl=Tclas.extract(icl);
							for (k=0; k<nk_BR; k++) newBR[k]+=cl.mean(k);
							cl=Tclas.extract(icl0);
							for (k=0; k<nk_BR; k++) newBR[k]-=cl.mean(k);
							for (k=0; k<nk_BR; k++) newBR[k]/=Nvs[0];
							xx=0.;
							for (k=0; k<nk_BR; k++) xx+=pow(pixBR[k]-newBR[k],2);
							Papriori[icl]=xx;
							if (Papriori[icl]<Papriori[0]) {
								Tima[(ipixl+j*nHRc+jHR)*nblayer]=icl;
								for (k=0; k<nk_BR; k++) simBR[k]=newBR[k];
								Papriori[0]=Papriori[icl];
								nchBR++;
								}
							}
						}
					it++;
					} while (it<10 && nchBR>0);
				nch+=nchBR;
				}
			}
			if (iter%1==0) {
				cout<<" iter. "<<setw(3)<<iter<<" => #changements = "<<setw(6)<<nch<<"\n";
			}
		} while (iter<itermax && nch>0);
// regularisation : 
// pour chaque pixel BR
//		pour chaque couple de pixels HR du pixel BR courant
//			si les labels sont differets
//				permuter les labels
//				calculer l'energie de voisinage associee au couple de pixels HR
//				decider ou non de conserver la permutation
//				mise a jour de l'image des labels selon decision precedente
	for (l=0; l<nbpix; l++) Tima[l*nblayer+1]=Tima[l*nblayer];
	int ***TnvsHR=NULL;
	TnvsHR=new int** [nHRl+2];
	for (iHR=-1; iHR<=nHRl; iHR++) {
		TnvsHR[iHR+1]=new int*[nHRc+2];
		for (jHR=-1; jHR<=nHRc; jHR++) TnvsHR[iHR+1][jHR+1]=new int[nclas+1];
		}
	imadata<float> imaU(nl_BR,nc_BR);
	iter=0;
	do {
		iter++;
		if (iter%1==0) nch=0.f;
		double Ur=0.;
		for (i=0; i<nl_BR; i++) { // pour chaque ligne BR
			for (j=0; j<nc_BR; j++) { // pour chaque colonne BR
/*				pixBR=ima_BR.pix(i,j); // valeur observee en basse resolution
// calcul de la valeur BR simulee d'apres l'image des labels et les centres des classes
				for (ii=0; ii<=nclas; ii++) Nvs[ii]=0;
				for (iHR=0; iHR<nHRl; iHR++) {
					ipixl=(i*nHRl+iHR)*nc_HR;
					for (jHR=0; jHR<nHRc; jHR++) {
						icl=Tima[(ipixl+j*nHRc+jHR)*nblayer+1]; // label considerant uniquement le terme d'attache aux donnees
						Nvs[icl]++;
						Nvs[0]++;
						}
					}
				for (k=0; k<nk_BR; k++) simBR[k]=0;
				for (ii=1; ii<=nclas; ii++) {
					cl=Tclas.extract(ii);
					for (k=0; k<nk_BR; k++) simBR[k]+=Nvs[ii]*cl.mean(k);
					}
				for (k=0; k<nk_BR; k++) simBR[k]/=Nvs[0];
				xx=0.;
				for (k=0; k<nk_BR; k++) xx+=pow(pixBR[k]-simBR[k],2);*/
//				cout<<" pixel BR "<<i<<" "<<j<<" ";
				for (iHR=-1; iHR<=nHRl; iHR++) 
					for (jHR=-1; jHR<=nHRc; jHR++)
						for (ii=0; ii<=nclas; ii++) TnvsHR[iHR+1][jHR+1][ii]=0;
				xx=0;
				for (iHR=-1; iHR<=nHRl; iHR++)
					if (i*nHRl+iHR>=0 && i*nHRl+iHR<nl_HR) {
						for (jHR=-1; jHR<=nHRc; jHR++) {
							if (j*nHRc+jHR>=0 && j*nHRc+jHR<nc_HR) {
								i2min=maxi(i*nHRl+iHR-1,0);
								i2max=mini(i*nHRl+iHR+1,nl_HR-1);
								j2min=maxi(j*nHRc+jHR-1,0);
								j2max=mini(j*nHRc+jHR+1,nc_HR-1);
								for (i2=i2min; i2<=i2max; i2++)
									for (j2=j2min; j2<=j2max; j2++) {
										TnvsHR[iHR+1][jHR+1][0]++;
										icl=Tima[(i2*nc_HR+j2)*nblayer+1];
										TnvsHR[iHR+1][jHR+1][icl]++;
										}
								TnvsHR[iHR+1][jHR+1][0]--;
								icl=Tima[((i*nHRl+iHR)*nc_HR+j*nHRc+jHR)*nblayer+1];
								TnvsHR[iHR+1][jHR+1][icl]--;
								xx+=0.5*(TnvsHR[iHR+1][jHR+1][0]-TnvsHR[iHR+1][jHR+1][icl]);
								}
							}
						}
				imaU(i,j)=(float)xx;
				Papriori[0]=xx;
				it=0;
				do {
					nchBR=0;
					for (iHR=0; iHR<nHRl; iHR++) { // pour chaque ligne HR
						ipixl=(i*nHRl+iHR)*nc_HR;
						for (jHR=0; jHR<nHRc; jHR++) { // pour chaque colonne HR
							icl0=Tima[(ipixl+j*nHRc+jHR)*nblayer+1];

							for (i2HR=iHR; i2HR<nHRl; i2HR++) 
								for (j2HR=0; j2HR<nHRc; j2HR++)
									if (i2HR>iHR || j2HR>jHR) {
										icl=Tima[((i*nHRl+i2HR)*nc_HR+j*nHRc+j2HR)*nblayer+1];
										if (icl0!=icl) {
											xx=Papriori[0];
											xx=xx+TnvsHR[i2HR+1][j2HR+1][icl]-TnvsHR[i2HR+1][j2HR+1][icl0];
											i2min=maxi(i2HR-1,0);
											i2max=mini(i2HR+1,nHRl-1);
											j2min=maxi(j2HR-1,0);
											j2max=mini(j2HR+1,nHRc-1);
											for (i2=i2min; i2<=i2max; i2++)
												for (j2=j2min; j2<=j2max; j2++) {
													TnvsHR[i2+1][j2+1][icl]--;
													TnvsHR[i2+1][j2+1][icl0]++;
													}
											xx=xx+TnvsHR[iHR+1][jHR+1][icl0]-TnvsHR[iHR+1][jHR+1][icl];
											i2min=maxi(iHR-1,0);
											i2max=mini(iHR+1,nHRl-1);
											j2min=maxi(jHR-1,0);
											j2max=mini(jHR+1,nHRc-1);
											for (i2=i2min; i2<=i2max; i2++)
												for (j2=j2min; j2<=j2max; j2++) {
													TnvsHR[i2+1][j2+1][icl0]--;
													TnvsHR[i2+1][j2+1][icl]++;
													}
											if (xx<Papriori[0]) {
												nchBR+=2;
												Tima[(ipixl+j*nHRc+jHR)*nblayer+1]=icl;
												Tima[((i*nHRl+i2HR)*nc_HR+j*nHRc+j2HR)*nblayer+1]=icl0;
												imaU(i,j)=(float)xx;
												Papriori[0]=xx;
												icl0=icl;
												}
											else {
												i2min=maxi(i2HR-1,0);
												i2max=mini(i2HR+1,nHRl-1);
												j2min=maxi(j2HR-1,0);
												j2max=mini(j2HR+1,nHRc-1);
												for (i2=i2min; i2<=i2max; i2++)
													for (j2=j2min; j2<=j2max; j2++) {
														TnvsHR[i2+1][j2+1][icl0]--;
														TnvsHR[i2+1][j2+1][icl]++;
														}
												i2min=maxi(iHR-1,0);
												i2max=mini(iHR+1,nHRl-1);
												j2min=maxi(jHR-1,0);
												j2max=mini(jHR+1,nHRc-1);
												for (i2=i2min; i2<=i2max; i2++)
													for (j2=j2min; j2<=j2max; j2++) {
														TnvsHR[i2+1][j2+1][icl]--;
														TnvsHR[i2+1][j2+1][icl0]++;
														}
												}
/*				for (int iiHR=-1; iiHR<=nHRl; iiHR++) 
					for (int jjHR=-1; jjHR<=nHRc; jjHR++)
						for (ii=0; ii<=nclas; ii++) TnvsHR[iiHR+1][jjHR+1][ii]=0;
				xx=0;
				for (iiHR=-1; iiHR<=nHRl; iiHR++)
					if (i*nHRl+iiHR>=0 && i*nHRl+iiHR<nl_HR) {
						for (int jjHR=-1; jjHR<=nHRc; jjHR++) {
							if (j*nHRc+jjHR>=0 && j*nHRc+jjHR<nc_HR) {
								i2min=maxi(i*nHRl+iiHR-1,0);
								i2max=mini(i*nHRl+iiHR+1,nl_HR-1);
								j2min=maxi(j*nHRc+jjHR-1,0);
								j2max=mini(j*nHRc+jjHR+1,nc_HR-1);
								for (int ii2=i2min; ii2<=i2max; ii2++)
									for (int jj2=j2min; jj2<=j2max; jj2++) {
										TnvsHR[iiHR+1][jjHR+1][0]++;
										int icl2=Tima[(ii2*nc_HR+jj2)*nblayer+1];
										TnvsHR[iiHR+1][jjHR+1][icl2]++;
										}
								TnvsHR[iiHR+1][jjHR+1][0]--;
								int icl2=Tima[((i*nHRl+iiHR)*nc_HR+j*nHRc+jjHR)*nblayer+1];
								TnvsHR[iiHR+1][jjHR+1][icl2]--;
								xx+=0.5*(TnvsHR[iiHR+1][jjHR+1][0]-TnvsHR[iiHR+1][jjHR+1][icl2]);
								}
							}
						}
					if (xx!=imaU(i,j) || xx>250 || xx<0 || imaU(i,j)>250 || imaU(i,j)<0) {
						cout<<" calcul de l'energie xx = "<<xx<<" imaU = "<<imaU(i,j)<<"\n";
						char cc; cin>>cc;
						}*/
											}
										}
							}
						}
					it++;
					Ur+=xx;
					nch+=nchBR;
//					cout<<" # ch = "<<nchBR;
//					char aa; cin>>aa;
					} while (it<10 && nchBR>0);
//				cout<<"\n";
				nch+=nchBR;
				}
			}
			if (iter%1==0) {
				cout<<" iter. "<<setw(3)<<iter<<" U = "<<setw(9)<<setprecision(6)<<Ur;
				cout<<" => #changements = "<<setw(6)<<nch<<"\n";
			}
		} while (iter<itermax && nch>0);
		for (iHR=-1; iHR<=nHRl; iHR++) {
			for (jHR=-1; jHR<=nHRc; jHR++) if (TnvsHR[iHR+1][jHR+1]!=NULL) delete TnvsHR[iHR+1][jHR+1];
			if (TnvsHR[iHR+1]!=NULL) delete[] TnvsHR[iHR+1];
			}
		if (TnvsHR!=NULL) delete[] TnvsHR;
	cout<<" fin de la convergence pour "<<iter<<" iterations, ";
	cout<<nch<<" pixels encore changes\n";
	if (Nvs!=NULL) delete[] Nvs;
	if (Papriori!=NULL) delete[] Papriori;
	}
