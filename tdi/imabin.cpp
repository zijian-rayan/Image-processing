#include "imabin.h"

void imabin::affiche(const int ip) const {
	cout<<" image binaire de "<<nblig<<" lignes et "<<nbcol<<" colonnes \n";
	long int n=0;
	for (int i=0; i<nblig; i+=ip) {
		n=i*nbcol;
		for (int j=0; j<nbcol; j+=ip) {
			if (n>=0 && n<nbpix) cout<<setw(2)<<Tima[n];
			n+=ip;
		}
		cout<<"\n";
	}
}

imabin::operator imadata<float>() {
	int i,j;
	imadata<float> ima_f(nblig,nbcol,1);
	float x;
	for (i=0; i<nblig; i++) {
		for (j=0; j<nbcol; j++) {
			bool * adval=(bool*)Tsites[i*nbcol+j];
			x=*(adval);
			ima_f(i,j,0)=x;
			ima_f[i*nbcol+j]=&(ima_f(i,j,0));
			}
		}
	return ima_f;
}

imabin::operator imadata<BYTE>() {
	int i,j;
	imadata<BYTE> ima_u(nblig,nbcol,1);
	int x;
	for (i=0; i<nblig; i++) {
		for (j=0; j<nbcol; j++) {
			bool * adval=(bool*)Tsites[i*nbcol+j];
			x=*(adval);
			ima_u(i,j,0)=x*255;
			ima_u[i*nbcol+j]=&(ima_u(i,j,0));
			}
		}
	return ima_u;
}

imadata<BYTE> imabin::imaunsignedchar(const bool etaldyn) {
	int i,j,x;
	imadata<BYTE> ima_u(nblig,nbcol,1);
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) {
			x=(int)(*this)(i,j);
			if (etaldyn) x*=255;
			ima_u(i,j,0)=(BYTE)x;
		}
	return ima_u;
}

imabin imabin::operator+ (const imabin &ima2) {
	int nlig=mini(nblig,ima2.nlig()),ncol=mini(nbcol,ima2.ncol()),i,j;
	imabin imaRes(nlig,ncol);
	bool b;
	for (i=0; i<nlig; i++)
		for (j=0; j<ncol; j++) {
			b=((int)(*this)(i,j)!=(int)ima2(i,j));
			imaRes(i,j)=b;
		}
	return imaRes;
}

imabin imabin::operator- (const imabin &ima2) {
	return ((*this)+ima2);
}

imabin imabin::operator|| (const imabin &ima2) {
	int nlig=mini(nblig,ima2.nlig()),ncol=mini(nbcol,ima2.ncol()),i,j;
	imabin imaRes(nlig,ncol);
	bool b;
	for (i=0; i<nlig; i++)
		for (j=0; j<ncol; j++) {
			b=((int)(*this)(i,j)||(int)ima2(i,j));
			imaRes(i,j)=b;
		}
	return imaRes;
}

imabin imabin::operator&& (const imabin &ima2) {
	int nlig=mini(nblig,ima2.nlig()),ncol=mini(nbcol,ima2.ncol()),i,j;
	imabin imaRes(nlig,ncol);
	bool b;
	for (i=0; i<nlig; i++)
		for (j=0; j<ncol; j++) {
			b=((int)(*this)(i,j)&&(int)ima2(i,j));
			imaRes(i,j)=b;
		}
	return imaRes;
}

bool imabin::operator== (const imabin &ima2) {
	bool b;
	int i,j;
	b=1; i=0; j=0;
	while (b && i<nblig && j<nbcol) {
		if ((*this)(i,j)!=ima2(i,j)) b=0;
		else {
			if (j<nbcol-1) j++;
			else {i++; j=0;}
		}
	}
	return b;
}

bool imabin::operator!= (const imabin &ima2) {
	return !((*this)==ima2);
}

void imabin::mise_a_zero() {
	int i,j;
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) (*this)(i,j)=0;
}

void imabin::mise_a_un() {
	int i,j;
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) (*this)(i,j)=1;
}
imabin imabin::negatif() {
	int i,j;
	imabin imaRes(nblig,nbcol);
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) imaRes(i,j)=!(*this)(i,j);
	return imaRes;
}

imabin imabin::ajoutbords() {
	int i,j;
	imabin imaRes(nblig+2,nbcol+2);
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++)
			imaRes(i+1,j+1)=(*this)(i,j);
	return imaRes;
}

imabin imabin::retirebords() {
	int i,j;
	imabin imaRes(nblig-2,nbcol-2);
	for (i=1; i<nblig-1; i++)
		for (j=1; j<nbcol-1; j++)
			imaRes(i-1,j-1)=(*this)(i,j);
	return imaRes;
}

float imabin::norm() const {
	int i,j;
	float n=0.f;
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++)
			n+=(*this)(i,j);
	return n;
}

void imabin::sauve_Ima(char *nomfich, bool bin) const {
	ofstream sortie;
	sortie.open(nomfich,ios::out|ios::binary);
	if (!sortie) {
		cout<<" ouverture de "<<nomfich<<" impossible\n";
		exit (-1);
		}
	bool * adval=(bool*)Tsites[0];
	BYTE val;
	int itype=0, sizeval=sizeof(BYTE);
	cout<<" image "<<nblig<<" lig.& "<<nbcol<<" col., de type ";
	cout<<typeid(*adval).name()<<" => itype = "<<itype<<"\n";
	sortie.write((char *)&nblig,sizeof(int));
	sortie.write((char *)&nbcol,sizeof(int));
	sortie.write((char *)&itype,sizeof(int));
	for (long int i=0; i<nbpix; i++) {
		val=Tima[i];
		if (!bin) val=val*255;
		sortie.write((char *)&val,sizeval);
		}
	sortie.close();
}

imadata<int> imabin::composantes_connexes (int &ncc, int iconnex, bool iaf) {
	if (iaf) cout<<" calcul des composantes connexes\n";
	imadata<int> imaRes(nblig,nbcol,1);
	imaRes.mise_a_zero();
	const int nT_eq=100000;
	int n,i,j,jj,k,l,m,p,*T_eq=new int[nT_eq],jdec=(iconnex==8?1:0);
	for (i=0; i<nT_eq; i++) T_eq[i]=0;
	bool icfl, ichg;
	n=0;
	if ((*this)(0,0)) imaRes(0,0)=++n; // traitement premiere ligne
	for (j=1; j<nbcol; j++)
		if ((*this)(0,j)) {
			k=imaRes(0,j-1); imaRes(0,j)=(k>0?k:++n); //if (k>0) imaRes(0,j)=k; else imaRes(0,j)=++n;
//			if (iaf) {cout<<" pixel (0,"<<j<<") "<<imaRes(0,j)<<"\n";}
		}
	for (i=1; i<nblig; i++) {					// traitement lignes suivantes
		if ((*this)(i,0)) {
			k=imaRes(i-1,0); imaRes(i,0)=(k>0?k:++n); //if (k>0) imaRes(i,0)=k; else imaRes(i,0)=++n;
//			if (iaf) {cout<<" pixel ("<<i<<",0) "<<imaRes(i,0)<<"\n";}
		}
		for (j=1; j<nbcol; j++) {
			if ((*this)(i,j)) {
				icfl=0; k=imaRes(i,j-1);
				for (jj=maxi(j-jdec,0); jj<=mini(j+jdec,nbcol-1); jj++) {
					l=imaRes(i-1,jj);
					if (l>0)
						if (k==0) k=l;
						else
							if (k!=l && T_eq[k]!=l && T_eq[l]!=k) {
								if (icfl && iaf) {
									cout<<" potentiel conflit labels ! configuration : connexite = "<<iconnex<<", "<<imaRes(i,j-1)<<"\n";
									cout<<" "<<imaRes(i-1,j-1)<<" "<<imaRes(i-1,j)<<" "<<imaRes(i-1,mini(j+1,nbcol-1))<<"\n";
									cout<<" label dans T_eq "<<T_eq[imaRes(i,j-1)]<<" "<<T_eq[imaRes(i-1,j-1)]<<" "<<T_eq[imaRes(i-1,j)]<<" "<<T_eq[imaRes(i-1,mini(j+1,nbcol-1))]<<"\n";
								}
								icfl=1; //if (iaf) cout<<" T_eq "<<k<<" "<<l<<" "<<mini(k,l)<<" "<<T_eq[mini(k,l)]<<" ";
								m=mini(k,l); p=maxi(k,l);
								if (T_eq[m]!=0) m=T_eq[m];
								if (T_eq[p]!=m) {
									bool iok=0;
									if (!iok && T_eq[p]==0) {T_eq[p]=m; iok=1;}
									if (!iok && T_eq[p]!=0 && T_eq[p]>m && T_eq[T_eq[p]]==0) {T_eq[T_eq[p]]=m; iok=1;}
									if (!iok && T_eq[m]==0 && T_eq[p]!=0) {T_eq[m]=T_eq[p]; iok=1;}
									if (!iok && T_eq[m]==0) {T_eq[m]=p; iok=1;}
									if (!iok) { 
//										if (iaf) cout<<" et aussi T_eq "<<maxi(k,l)<<" "<<m<<" "<<mini(k,l)<<" "<<T_eq[mini(k,l)]<<" ";
										if (m>T_eq[maxi(k,l)]) T_eq[m]=T_eq[maxi(k,l)]; 
										else {T_eq[T_eq[maxi(k,l)]]=m; T_eq[maxi(k,l)]=m;}
//										T_eq[maxi(k,l)]=T_eq[m]=mini(T_eq[maxi(k,l)],m);
									}
								}
//								if (iaf) {for (int ii=0; ii<nT_eq; ii++) if (T_eq[ii]!=0) cout<<" T_eq["<<ii<<"]="<<(int)T_eq[ii]<<"\n";}
								k=m;
							}
					}
				if (k>0) imaRes(i,j)=k;
				else {
					imaRes(i,j)=++n; //cout<<" increment de n en "<<i<<" "<<j<<" k="<<k<<" => n = "<<n<<"\n";
					if (n>=nT_eq) {
						cout<<"\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n Table d'equivalence trop petite !!!\n";
						cout<<"\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
						char aa; cin>>aa;
					}
				}
//				if (iaf) {cout<<" pixel ("<<i<<","<<j<<") "<<imaRes(i,j)<<"\n";}
			}
		}
	}
	do {
		ichg=0;
		for (i=0; i<nT_eq; i++) {
			l=T_eq[i];
			if (l!=0 && T_eq[l]!=0) {
				m=mini(l,T_eq[l]);
				if (m!=l) {
					T_eq[i]=T_eq[l];
					ichg=1;
				}
			}
		}
		if (iaf) {for (int ii=0; ii<nT_eq; ii++) if (T_eq[ii]!=0) cout<<" # T_eq["<<ii<<"]="<<(int)T_eq[ii]<<"\n"; char aa; cin>>aa;}
	} while (ichg);
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) {
			l=imaRes(i,j);
			if (T_eq[l]!=0) imaRes(i,j)=T_eq[l];
		}
// verification etiquettage en composantes connexes
	bool lastverif=1;
	do {
		for (i=0; i<nT_eq; i++) T_eq[i]=0;
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				l=imaRes(i,j);
				if (l>0) {
					if (i>0) {k=imaRes(i-1,j); if (k>0 && l!=k) {T_eq[maxi(k,l)]=mini(k,l); }}
					if (i<nblig-1) {k=imaRes(i+1,j); if (k>0 && l!=k) {T_eq[maxi(k,l)]=mini(k,l); }}
					if (j>0) {k=imaRes(i,j-1); if (k>0 && l!=k) {T_eq[maxi(k,l)]=mini(k,l); }}
					if (j<nbcol-1) {k=imaRes(i,j+1); if (k>0 && l!=k) {T_eq[maxi(k,l)]=mini(k,l); }}
//					if (iconnex==8) 
				}
			}
		lastverif=0;
		for (i=0; i<nT_eq; i++) 
			if (T_eq[i]>0) {cout<<"\n !!!! equivalence trouvee a posteriori entre "<<i<<" et "<<T_eq[i]<<" !!!!\n";	lastverif=1; /*char aa; cin>>aa;*/}
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {l=imaRes(i,j);	if (T_eq[l]!=0) imaRes(i,j)=T_eq[l];}
	} while (lastverif);
// stat sur composantes connexes
	for (i=0; i<nT_eq; i++) T_eq[i]=0;
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) T_eq[imaRes(i,j)]++;
	if (iaf) {for (i=0; i<nT_eq; i++) if (T_eq[i]!=0) cout<<" composante "<<i<<" : # pixels = "<<T_eq[i]<<"\n";}
	n=0;
	for (i=1; i<nT_eq; i++) 
		if (T_eq[i]!=0) T_eq[i]=++n;
	if (iaf) cout<<" # composantes connexes trouvees = "<<n<<" labelisees entre 1 et "<<n<<"\n";
	ncc=0;
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) {
			l=imaRes(i,j);
			if (iaf && l!=0 && T_eq[l]!=0 && l!=T_eq[l]) cout<<" l = "<<l<<" et T_eq[l] = "<<T_eq[l]<<"\n";
			if (l!=0 && T_eq[l]!=0) imaRes(i,j)=T_eq[l];
			if (imaRes(i,j)>ncc) ncc=imaRes(i,j);
		}
	if (iaf) cout<<" # composantes connexes trouvees = "<<ncc<<" labelisees entre 1 et "<<ncc<<"\n";
//	if (ncc>1) {int ncc2; cout<<" # composantes connexes trouvees = "<<ncc<<" labelisees entre 1 et "<<ncc<<"\n"; char aa; cin>>aa; composantes_connexes (ncc2,iconnex,1);}
	for (i=0; i<nT_eq; i++) T_eq[i]=0;
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) T_eq[imaRes(i,j)]++;
	if (iaf) {for (i=0; i<=ncc; i++) if (T_eq[i]!=0) cout<<" composante "<<i<<" : # pixels = "<<T_eq[i]<<"\n";}
	if (T_eq!=NULL) delete[] T_eq; if (iaf) cout<<" fin composantes_connexes\n";
	return imaRes;
}

imadata<int> imabin::boxes_rectangulaires (int &ncc, bool iaf) {
	if (iaf) cout<<" decomposition sous forme de rectangles\n";
  imadata<int> imaRes(nblig,nbcol,1); imaRes.mise_a_zero();
	int n,i,j,j0,j2,jj,k;
	bool iok;
  n=0;
	for (j=0; j<nbcol; j++) { // traitement de la premiere ligne a part car pas de ligne precedente
    k=0; 
    if ((*this)(0,j)) {
			if (j>0) k=imaRes(0,j-1);
			if (k>0) imaRes(0,j)=k;
			else {imaRes(0,j)=++n; /*if (iaf) cout<<" box "<<n<<" commence en (0,"<<j<<")\n";*/}
		}
  }
	for (i=1; i<nblig; i++) { //cout<<" i = "<<i<<" n = "<<n<<"\n"; 
		for (j=0; j<nbcol; j++) {
			if ((*this)(i,j)) {
				k=imaRes(i-1,j); iok=1; j0=j; j2=j;
				if (k==0) iok=0;
				else {
					if (j>0 && imaRes(i-1,j-1)==k) iok=0;
					else {
						do {if (!(*this)(i,j2++)) iok=0;} while (j2<nbcol && imaRes(i-1,j2)==k);
					}
				}
				if (iok) {
					for (jj=j0; jj<j2; jj++) imaRes(i,jj)=k;
					j=j2-1;
				} else {
					n++; jj=j; if (iaf) cout<<" box "<<n<<" commence en ("<<i<<","<<j<<")\n";
					do {imaRes(i,jj++)=n;} while (jj<nbcol && (*this)(i,jj));
					j=jj-1;
				}
			}
		}
	}
	ncc=n;
	if (iaf) {for (i=0; i<nblig; i++) {for (j=0; j<nbcol; j++) cout<<(int)imaRes(i,j)<<" "; cout<<"\n";}}
	if (iaf) cout<<" # rectangles trouves = "<<ncc<<" labelisees entre 1 et "<<ncc<<"\n";
	return imaRes;
}

imabin imabin::dilate(eltstruct B) {
	int i,j,i0,i2,j0,j2,ii,jj/*,nl2,nc2*/,iB;
	int nl=B.nl,nc=B.nc, x0=B.x0, y0=B.y0;
	bool z, id;
	imabin imaRes(*this);
	id=0;
/*	if ((nl>4 || nc>4) && (nl%2)==1 && (nc%2)==1 &&
	    (nl==nc || nl==1 || nc==1) ) {
		nl2=maxi(1,(nl-3)/2); // cas des elements struct non centrés????
		nc2=maxi(1,(nc-3)/2);
		imabin imaB(nl,nc), imaC(nl,nc);
		for (i=0; i<nl; i++)
			for (j=0; j<nc; j++) imaC(i,j)=B.aT[i*nc+j];
		if (!id && nl==1) {
			eltstruct B1(1,3);
			for (i=0; i<3; i++) B1.aT[i]=1;
			for (j=nc/2-1; j<=nc/2+1; j++) imaB(0,j)=1;
			if (imaB.dilate(B1.transpose()????,nc2)==imaC) {
				id=1;
				imaRes=dilate(B1,nc2);
				}
			}
		if (!id && nc==1) {
			eltstruct B1(3,1);
			for (i=0; i<3; i++) B1.aT[i]=1;
			for (i=nl/2-1; i<=nl/2+1; i++) imaB(i,0)=1;
			if (imaB.dilate(B1.transpose()????,nc2)==imaC) {
				id=1;
				imaRes=dilate(B1,nc2);
				}
			}
		if (!id && nl==nc) {
			eltstruct B1(3,3);
			for (i=0; i<9; i++) B1.aT[i]=1;
			for (i=nl/2-1; i<=nl/2+1; i++)
				for (j=nc/2-1; j<=nc/2+1; j++) imaB(i,j)=1;
			if (imaB.dilate(B1.transpose()????,nl2)==imaC) {
				id=1;
				imaRes=dilate(B1,nl2);
				}
			if (!id) {
				for (i=0; i<9; i++) B1.aT[i]=1;
				B1.aT[0]=0; B1.aT[2]=0; B1.aT[6]=0; B1.aT[8]=0;
				for (i=0; i<nl; i++)
					for (j=0; j<nc; j++) imaB(i,j)=0;
			for (i=nl/2-1; i<=nl/2+1; i++)
				for (j=nc/2-1; j<=nc/2+1; j++)
					if (i==nl/2 || j==nc/2) imaB(i,j)=1;
				if (imaB.dilate(B1,nl2)==imaC) {
					id=1;
					imaRes=dilate(B1,nl2);
					}
				}
			if (!id) {
				for (i=0; i<9; i++) B1.aT[i]=0;
				B1.aT[0]=1; B1.aT[4]=1; B1.aT[8]=1;
				for (i=0; i<nl; i++)
					for (j=0; j<nc; j++) imaB(i,j)=0;
				imaB(nl/2-1,nc/2-1)=1; imaB(nl/2,nc/2)=1;
				imaB(nl/2+1,nc/2+1)=1;
				if (imaB.dilate(B1,nl2)==imaC) {
					id=1;
					imaRes=dilate(B1,nl2);
					}
				}
			if (!id) {
				for (i=0; i<9; i++) B1.aT[i]=0;
				B1.aT[2]=1; B1.aT[4]=1; B1.aT[6]=1;
				for (i=0; i<nl; i++)
					for (j=0; j<nc; j++) imaB(i,j)=0;
				imaB(nl/2-1,nc/2+1)=1; imaB(nl/2,nc/2)=1;
				imaB(nl/2+1,nc/2-1)=1;
				if (imaB.dilate(B1,nl2)==imaC) {
					id=1;
					imaRes=dilate(B1,nl2);
					}
				}
			}
		}*/
	if (id==0) {
		for (i=0; i<nblig; i++) {
			i0=maxi(i-y0,0);
			i2=mini(i+nl-1-y0,nblig-1);
			for (j=0; j<nbcol; j++) {
				j0=maxi(j-x0,0);
				j2=mini(j+nc-1-x0,nbcol-1);
				z=0;
				for (ii=i0; ii<=i2; ii++) {
					if (!z) {
						iB=(ii-i+y0)*nc;
						for (jj=j0; jj<=j2; jj++)
							if (!z && (*this)(ii,jj) && B.aT[iB+jj-j+x0]) z=1;
						}
					}
				imaRes(i,j)=z;
				}
			}
		}
	return imaRes;
}

imabin imabin::erode(eltstruct B) {
	int i,j,i0,i2,j0,j2,ii,jj,iB;
	int nl=B.nl,nc=B.nc, x0=B.x0, y0=B.y0;
	bool z;
	imabin imaRes(*this);
	for (i=0; i<nblig; i++) {
		i0=maxi(i-y0,0);
		i2=mini(i+nl-1-y0,nblig-1);
		for (j=0; j<nbcol; j++) {
			j0=maxi(j-x0,0);
			j2=mini(j+nc-1-x0,nbcol-1);
			z=1;
			for (ii=i0; ii<=i2; ii++) {
				if (z) {
					iB=(ii-i+y0)*nc;
					for (jj=j0; jj<=j2; jj++)
						if (z && !(*this)(ii,jj) && B.aT[iB+jj-j+x0]) z=0;
					}
				}
			imaRes(i,j)=z;
			}
		}
	return imaRes;
}

imabin imabin::dilate(eltstruct B1, int n) { // calcul de B1 a partir de B ?
	int i,j,i0,i2,j0,j2,ii,jj,iB,nn;
	int nl=B1.nl,nc=B1.nc, x0=B1.x0, y0=B1.y0;
	bool z;
	imabin imaR(*this), imaC(*this);
	for (nn=0; nn<n; nn++) {
		for (i=0; i<nblig; i++) {
			i0=maxi(i-y0,0);
			i2=mini(i+nl-1-y0,nblig-1);
			for (j=0; j<nbcol; j++) {
				j0=maxi(j-x0,0);
				j2=mini(j+nc-1-x0,nbcol-1);
				z=0;
				for (ii=i0; ii<=i2; ii++) {
					if (!z) {
						iB=(ii-i+y0)*nc;
						for (jj=j0; jj<=j2; jj++)
							if (!z && B1.aT[iB+jj-j+x0] && imaC(ii,jj)) z=1;
						}
					}
				imaR(i,j)=z;
				}
			}
		imaC=imaR;
		}
	return imaR;
}

imabin imabin::erode(eltstruct B1, int n) {
	int i,j,i0,i2,j0,j2,ii,jj,iB,nn;
	int nl=B1.nl, nc=B1.nc, x0=B1.x0, y0=B1.y0;
	bool z;
	imabin imaR(*this), imaC(*this);
	for (nn=0; nn<n; nn++) {
		for (i=0; i<nblig; i++) {
			i0=maxi(i-y0,0);
			i2=mini(i+nl-1-y0,nblig-1);
			for (j=0; j<nbcol; j++) {
				j0=maxi(j-x0,0);
				j2=mini(j+nc-1-x0,nbcol-1);
				z=1;
				for (ii=i0; ii<=i2; ii++) {
					if (z) {
						iB=(ii-i+y0)*nc;
						for (jj=j0; jj<=j2; jj++)
							if (z && !imaC(ii,jj) && B1.aT[iB+jj-j+x0]) z=0;
						}
					}
				imaR(i,j)=z;
				}
			}
		imaC=imaR;
		}
	return imaR;
}

imabin imabin::dilate(eltstruct B1, eltstruct B2) { // calcul de B1,B2 a partir de B ?
	imabin imaR(nblig,nbcol), imaC(nblig,nbcol);
	imaC=dilate(B1);
	imaR=imaC.dilate(B2);
	return imaR;
}

imabin imabin::erode(eltstruct B1, eltstruct B2) { // calcul de B1,B2 a partir de B ?
	imabin imaR(nblig,nbcol), imaC(nblig,nbcol);
	imaC=erode(B1);
	imaR=imaC.erode(B2);
	return imaR;
}

imadata<float> imabin::Tr_dist (const imabin &Masque, int inoyau) const {
	float *noyau=NULL;
	int dim,dim1,dim2,i,j,i0,i2,j0,j2,ii,jj,k;
	imadata<float> ima_d(nblig,nbcol,1);
	float dist, x;
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) ima_d(i,j)=FLT_MAX*(1-(*this)(i,j));
	if (inoyau<1 || inoyau>4) {
		cout<<" parametre inoyau "<<inoyau<<" hors domaine => valeur par defaut = 1\n";
		inoyau=1;
	}
	switch (inoyau) {
		case 2: noyau=new float[dim=9];
				for (i=0; i<dim; i++) noyau[i]=1.f;
				noyau[(dim-1)/2]=0.f; dim2=3; dim1=1; break;
		case 3: noyau=new float[dim=9];
				for (i=0; i<dim; i++) {
					if (i%2==0) noyau[i]=4.f/3.f;
					else noyau[i]=1.f; }
				noyau[(dim-1)/2]=0.f; dim2=3; dim1=1; break;
		case 4: noyau=new float[dim=25];
				noyau[1]=noyau[3]=noyau[5]=noyau[9]=11.f/5.f;
				noyau[15]=noyau[19]=noyau[21]=noyau[23]=11.f/5.f;
				noyau[6]=noyau[8]=noyau[16]=noyau[18]=7.f/5.f;
				noyau[7]=noyau[11]=noyau[13]=noyau[17]=1.f;
				noyau[(dim-1)/2]=0.f;
				noyau[0]=noyau[2]=noyau[4]=noyau[10]=FLT_MAX;
				noyau[14]=noyau[20]=noyau[22]=noyau[24]=FLT_MAX;
				dim2=5; dim1=2; break;
	}
	int nch, it=0;
	do {cout<<" iteration "<<++it<<"\n";
		nch=0;
		if (inoyau==1) {                                              // cas 4-connexite
			for (i=0; i<nblig; i++)
				for (j=0; j<nbcol; j++)
					if (Masque(i,j)) {
						dist=ima_d(i,j);
						if (i>0 && Masque(i-1,j)) {x=ima_d(i-1,j)+1; if (x<dist) dist=x;}
						if (j>0 && Masque(i,j-1)) {x=ima_d(i,j-1)+1; if (x<dist) dist=x;}
						if (ima_d(i,j)>dist) {ima_d(i,j)=dist; nch++;}
					}
			for (i=nblig-1; i>=0; i--)
				for (j=nbcol-1; j>=0; j--) 
					if (Masque(i,j)) {
						dist=ima_d(i,j);
						if (i<nblig-1 && Masque(i+1,j)) {x=ima_d(i+1,j)+1; if (x<dist) dist=x;}
						if (j<nbcol-1 && Masque(i,j+1)) {x=ima_d(i,j+1)+1; if (x<dist) dist=x;}
						if (ima_d(i,j)>dist) {ima_d(i,j)=dist; nch++;}
					}
		} else {
			for (i=0; i<nblig; i++) {
				i0=maxi(0,i-dim1);
				for (j=0; j<nbcol; j++) {
					j0=maxi(0,j-dim1); j2=mini(j+dim1,nbcol-1);
					if (Masque(i,j)) {
						dist=ima_d(i,j);
						for (ii=i0; ii<=i; ii++)
							for (jj=j0; jj<=j2; jj++)
								if (Masque(ii,jj)) {
									k=(ii-i+dim1)*dim2+(jj-j+dim1);
									if (k<(dim-1)/2) {
										x=ima_d(ii,jj)+noyau[k];
										if (x<dist) dist=x;
									}
								}
						if (ima_d(i,j)>dist) {ima_d(i,j)=dist; nch++;}
					}
				}
			}
			for (i=nblig-1; i>=0; i--) {
				i2=mini(i+dim1,nblig-1);
				for (j=nbcol-1; j>=0; j--) {
					j0=maxi(0,j-dim1); j2=mini(j+dim1,nbcol-1);
					if (Masque(i,j)) {
						dist=ima_d(i,j);
						for (ii=i; ii<=i2; ii++)
							for (jj=j0; jj<=j2; jj++)
								if (Masque(ii,jj)) {
									k=(ii-i+dim1)*dim2+(jj-j+dim1);
									if (k>(dim-1)/2) {
										x=ima_d(ii,jj)+noyau[k];
										if (x<dist) dist=x;
									}
								}
						if (ima_d(i,j)>dist) {ima_d(i,j)=dist; nch++;}
					}
				}
			}
		}
	} while (nch>0);
	if (noyau!=NULL) delete[] noyau;
	float xinfini=FLT_MAX;
	if (ima_d.maxI()>=xinfini) {
		float xmax=0.f;
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++)
				if (ima_d(i,j)<xinfini && ima_d(i,j)>xmax) xmax=ima_d(i,j);
		cout<<" Max distance hors infini = "<<xmax<<"\n";
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++)
				if (ima_d(i,j)>=xinfini) ima_d(i,j)=2*xmax;
	} else cout<<" Max distance = "<<ima_d.maxI()<<" < +Inf "<<xinfini<<"\n";
	return ima_d;
}

imadata<float> imabin::Tr_dist (int inoyau) const {
	int *noyau=NULL;
	int dim,dim1,dim2,i,j,i0,i2,j0,j2,ii,jj,k;
	imadata<float> ima_d(nblig,nbcol,1);
	float dist,x, coef;
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++)	ima_d(i,j)=FLT_MAX*(1-(*this)(i,j));
	if (inoyau<1 || inoyau>4) {
		cout<<" parametre inoyau "<<inoyau<<" hors domaine => valeur par defaut = 1\n";
		inoyau=1;
	}
	if (inoyau==1) { //cout<<" noyau 4 connexite\n";                    // cas 4-connexite
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				dist=ima_d(i,j);
				if (i>0) {x=ima_d(i-1,j)+1; if (x<dist) dist=x;}
				if (j>0) {x=ima_d(i,j-1)+1; if (x<dist) dist=x;}
				ima_d(i,j)=dist;
			}
		for (i=nblig-1; i>=0; i--)
			for (j=nbcol-1; j>=0; j--) {
				dist=ima_d(i,j);
				if (i<nblig-1) {x=ima_d(i+1,j)+1; if (x<dist) dist=x;}
				if (j<nbcol-1) {x=ima_d(i,j+1)+1; if (x<dist) dist=x;}
				ima_d(i,j)=dist;
			}
	} else {
		switch (inoyau) {
			case 2: noyau=new int[dim=9]; //cout<<" noyau 8 connexite\n";
				for (i=0; i<dim; i++) noyau[i]=1;
				noyau[(dim-1)/2]=0; dim2=3; dim1=1; coef=1.; break;
			case 3: noyau=new int[dim=9]; //cout<<" noyau Euclidien 3x3\n";
				for (i=0; i<dim; i++) {
					if (i%2==0) noyau[i]=4;
					else noyau[i]=3;
				}
				noyau[(dim-1)/2]=0;	dim2=3; dim1=1; coef=3.; break;
			case 4: noyau=new int[dim=25]; //cout<<" noyau Euclidien 5x5\n";
				noyau[1]=noyau[3]=noyau[5]=noyau[9]=noyau[15]=noyau[19]=noyau[21]=noyau[23]=11;
				noyau[6]=noyau[8]=noyau[16]=noyau[18]=7;
				noyau[7]=noyau[11]=noyau[13]=noyau[17]=5;
				noyau[(dim-1)/2]=0;
				noyau[0]=noyau[2]=noyau[4]=noyau[10]=32768;
				noyau[14]=noyau[20]=noyau[22]=noyau[24]=32768;
				dim2=5; dim1=2; coef=5.; break;
		}
		for (i=0; i<nblig; i++) {
			i0=maxi(0,i-dim1);
			for (j=0; j<nbcol; j++) {
				j0=maxi(0,j-dim1); j2=mini(j+dim1,nbcol-1);
				dist=ima_d(i,j);
				for (ii=i0; ii<=i; ii++)
					for (jj=j0; jj<=j2; jj++) {
						k=(ii-i+dim1)*dim2+(jj-j+dim1);
						if (k<(dim-1)/2) {
							x=ima_d(ii,jj)+noyau[k];
							if (x<dist) dist=x;
						}
					}
				ima_d(i,j)=dist;
			}
		}
		for (i=nblig-1; i>=0; i--) {
			i2=mini(i+dim1,nblig-1);
			for (j=nbcol-1; j>=0; j--) {
				j0=maxi(0,j-dim1); j2=mini(j+dim1,nbcol-1);
				dist=ima_d(i,j);
				for (ii=i; ii<=i2; ii++)
					for (jj=j0; jj<=j2; jj++) {
						k=(ii-i+dim1)*dim2+(jj-j+dim1);
						if (k>(dim-1)/2) {
							x=ima_d(ii,jj)+noyau[k];
							if (x<dist) dist=x;
						}
					}
				ima_d(i,j)=dist;
			}
		}
		if (inoyau!=2) {
			for (i=0; i<nblig; i++)
				for (j=0; j<nbcol; j++) ima_d(i,j)=ima_d(i,j)/coef;
		}
	}
	if (noyau!=NULL) delete[] noyau;
	return ima_d;
}

imabin imabin::dilate(const imadata<float> &ima_d, float seuil) {
	imabin imaRes(nblig,nbcol);
	int i,j;
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) {
			if (ima_d(i,j)>seuil) imaRes(i,j)=0;
			else imaRes(i,j)=1;
			}
	return imaRes;
}

imabin imabin::dilate(int inoyau, float seuil) {
	imabin imaRes(nblig,nbcol);
	imadata<float> ima_d;
	ima_d=Tr_dist(inoyau);
	imaRes=dilate(ima_d,seuil);
	return imaRes;
}

imabin imabin::erode(const imadata<float> &ima_d, float seuil) {
	imabin imaRes(nblig,nbcol);
	int i,j;
	for (i=0; i<nblig; i++) {
		for (j=0; j<nbcol; j++) {
			if (ima_d(i,j)>seuil) imaRes(i,j)=1;
			else imaRes(i,j)=0;
			}
		}
	return imaRes;
}

imabin imabin::erode(int inoyau, float seuil) {
	imabin imaRes(nblig,nbcol);
	imadata<float> ima_d=negatif().Tr_dist(inoyau);
	imaRes=erode(ima_d,seuil);
	return imaRes;
}

imabin imabin::ouverture(eltstruct B) {
	imabin imaRes(nblig,nbcol);
	imaRes=erode(B).dilate(B.transpose());
	return imaRes;
}

imabin imabin::ouverture(eltstruct B, int k) {
	imabin imaRes(nblig,nbcol);
	imaRes=erode(B,k).dilate(B.transpose(),k);
	return imaRes;
}

imabin imabin::fermeture(eltstruct B) {
	int i,j,dl0,dl2,dc0,dc2;
	dl0=B.nl-1-B.y0; dl2=B.y0; dc0=B.nc-1-B.x0; dc2=B.x0;
	imabin imaC(nblig+dl0+dl2,nbcol+dc0+dc2), imaR(nblig+dl0+dl2,nbcol+dc0+dc2);
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) imaC(i+dl0,j+dc0)=(*this)(i,j);
	imaR=imaC.dilate(B).erode(B.transpose());
	imabin imaRes(nblig,nbcol);
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) imaRes(i,j)=imaR(i+dl0,j+dc0);
	return imaRes;
}

imabin imabin::fermeture(eltstruct B, int k) {
	int i,j,dl0,dl2,dc0,dc2;
	dl0=B.nl-1-B.y0; dl2=B.y0; dc0=B.nc-1-B.x0; dc2=B.x0;
	imabin imaC(nblig+dl0+dl2,nbcol+dc0+dc2), imaR(nblig+dl0+dl2,nbcol+dc0+dc2);
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) imaC(i+dl0,j+dc0)=(*this)(i,j);
	imaR=imaC.dilate(B,k).erode(B.transpose(),k);
	imabin imaRes(nblig,nbcol);
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) imaRes(i,j)=imaR(i+dl0,j+dc0);
	return imaRes;
}

imabin imabin::tophat(eltstruct B) {
	imabin imaC=ouverture(B), imaRes=(*this)-imaC;
	return imaRes;
}

imabin imabin::tophat_c(eltstruct B) {
	imabin imaRes=fermeture(B)-(*this);
	return imaRes;
}

imabin imabin::reconstruction_geodesique(imabin &marq, int iconnex) {
	imabin imaRes(nblig,nbcol);
	int i,j,i0,i2,j0,j2,ii,jj;
	if (marq.nblig<nblig || marq.nbcol<nbcol) {
		imabin marqE(nblig,nbcol);
		int nlig=mini(nblig,marq.nblig), ncol=mini(nbcol,marq.nbcol);
		for (i=0; i<nlig; i++)
			for (j=0; j<ncol; j++) marqE(i,j)=marq(i,j);
		marq=marqE;
	}
	imabin imadejadsLpix(nblig,nbcol);
	liste_pixels Lpix;
	elt_liste E;
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) if ((*this)(i,j) && marq(i,j)) {Lpix.insere(i,j); imadejadsLpix(i,j)=1;}
	while (Lpix.nb_elts()>0) {
		E=Lpix.extrait();
		i=E.x; j=E.y;
		imaRes(i,j)=1;
		i0=maxi(i-1,0); i2=mini(i+1,nblig-1);
		j0=maxi(j-1,0); j2=mini(j+1,nbcol-1);
		if (iconnex==4) {
			for (ii=i0; ii<=i2; ii++)
				for (jj=j0; jj<=j2; jj++)
//					if (imaRes(ii,jj)!=1 && (ii==i || jj==j) && (*this)(ii,jj) && !Lpix.existe(ii,jj)) Lpix.insere(ii,jj);
					if (imaRes(ii,jj)!=1 && (ii==i || jj==j) && (*this)(ii,jj) && !imadejadsLpix(ii,jj)) {
						Lpix.insere(ii,jj);
						imadejadsLpix(ii,jj)=1;
					}
		}
		else {
			for (ii=i0; ii<=i2; ii++)
				for (jj=j0; jj<=j2; jj++)
//					if (imaRes(ii,jj)!=1 && (*this)(ii,jj) && !Lpix.existe(ii,jj)) Lpix.insere(ii,jj);
					if (imaRes(ii,jj)!=1 && (*this)(ii,jj) && !imadejadsLpix(ii,jj)) {
						Lpix.insere(ii,jj);
						imadejadsLpix(ii,jj)=1;
					}
			}
		}
//	cout<<"\n";
	return imaRes;
}

imabin imabin::reconstruction_geodesique(int y0, int x0, int iconnex) {
	imabin imaRes(nblig,nbcol);
	liste_pixels Lpix;
	elt_liste E;
	int i,j,i0,i2,j0,j2,ii,jj;
	if ((*this)(y0,x0)) Lpix.insere(y0,x0);
	while (Lpix.nb_elts()>0) {
		E=Lpix.extrait();
		i=E.x; j=E.y;
		imaRes(i,j)=1;
		i0=maxi(i-1,0); i2=mini(i+1,nblig-1);
		j0=maxi(j-1,0); j2=mini(j+1,nbcol-1);
		if (iconnex==4) {
			for (ii=i0; ii<=i2; ii++)
				for (jj=j0; jj<=j2; jj++)
					if (imaRes(ii,jj)!=1 && (ii==i || jj==j) &&
					    (*this)(ii,jj) && !Lpix.existe(ii,jj)) Lpix.insere(ii,jj);
			}
		else {
			for (ii=i0; ii<=i2; ii++)
				for (jj=j0; jj<=j2; jj++)
					if (imaRes(ii,jj)!=1 && (*this)(ii,jj) && !Lpix.existe(ii,jj))
						Lpix.insere(ii,jj);
			}
		}
	return imaRes;
}

imabin imabin::bouche_trou(int iconnex) {
    imabin imadon(*this), imaRes(nblig,nbcol);
    int i,j,n, ncc;
		imadata<int> imacc=imadon.negatif().composantes_connexes(ncc,iconnex); 
		int *T_n=new int[ncc];
		for (i=0; i<ncc; i++) T_n[i]=0;
    n=nbcol-1;
    for (i=0; i<nblig; i++) {
      if (imacc(i,0)>0) T_n[imacc(i,0)-1]++;
      if (imacc(i,n)>0) T_n[imacc(i,n)-1]++;
    }
    n=nblig-1;
    for (i=1; i<nbcol-1; i++) {
      if (imacc(0,i)>0) T_n[imacc(0,i)-1]++;
      if (imacc(n,i)>0) T_n[imacc(n,i)-1]++;
    }
    int nmax=0, iccfond=0;
		for (i=0; i<ncc; i++) 
			if (T_n[i]>nmax) {nmax=T_n[i]; iccfond=i+1;}; cout<<" composante connexe representant le fond = "<<iccfond<<"\n";
		if (T_n!=NULL) delete[] T_n;
		if (iccfond!=0) {
			imaRes.mise_a_zero();
			for (i=0; i<nblig; i++) 
				for (j=0; j<nbcol; j++) imaRes(i,j)=(imacc(i,j)==iccfond?1:0);
		}
    return imaRes.negatif();
  }

imabin imabin::bouche_trou2(int iconnex) {
    imabin imadon2(nblig+2,nbcol+2), imaRes(nblig,nbcol);
    int i,j,n, ncc, iccfond=0;
		for (i=1; i<=nblig; i++) for (j=1; j<=nbcol; j++) imadon2(i,j)=(*this)(i-1,j-1);
		imadata<int> imacc=imadon2.negatif().composantes_connexes(ncc,iconnex); 
		iccfond=imacc(0,0); 
		if (iccfond!=imacc(0,nbcol+1) || iccfond!=imacc(nblig+1,nbcol+1) || iccfond!=imacc(nblig+1,0) ) {
			cout<<" Pb : plusieurs CC susceptibles de representer le fond...\n"; iccfond=0;
			int *T_n=new int[ncc]; for (i=0; i<ncc; i++) T_n[i]=0;
		  n=nbcol+1;
			for (i=0; i<nblig+2; i++) {if (imacc(i,0)>0) T_n[imacc(i,0)-1]++; if (imacc(i,n)>0) T_n[imacc(i,n)-1]++;}
			n=nblig+1;
			for (i=1; i<nbcol+2; i++) {if (imacc(0,i)>0) T_n[imacc(0,i)-1]++; if (imacc(n,i)>0) T_n[imacc(n,i)-1]++;}
			int nmax=0; for (i=0; i<ncc; i++) if (T_n[i]>nmax) {nmax=T_n[i]; iccfond=i+1;}; cout<<" composante connexe representant le fond = "<<iccfond<<"\n";
			if (T_n!=NULL) delete[] T_n;
		}
		if (iccfond!=0) {
			for (i=0; i<nblig; i++) 
				for (j=0; j<nbcol; j++) imaRes(i,j)=(imacc(i+1,j+1)==iccfond?1:0);
		}
    return imaRes.negatif();
  }

imadata<int> imabin::CompoConnexes_MM (int &ncc, int iconnex, bool iaf) {
	if (iaf) cout<<" calcul des composantes connexes\n";
	imadata<int> imaRes(nblig,nbcol,1);
	imaRes.mise_a_zero();
	imabin imaC(nblig,nbcol);
	int n,i,j,ii,jj;
	n=0;
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) {
			if (imaRes(i,j)==0 && (*this)(i,j)) {
				n++;
//				cout<<" composante "<<n<<"\n";
				imaC=reconstruction_geodesique(i,j,iconnex);
				for (ii=0; ii<nblig; ii++)
					for (jj=0; jj<nbcol; jj++)
						if (imaC(ii,jj)) imaRes(ii,jj)=n;
				}
			}
	ncc=n;
	if (iaf) cout<<" # composantes connexes trouvees = "<<ncc<<" labelisees entre 1 et "<<ncc<<"\n";
	return imaRes;
}

imabin imabin::erode_ultime(const imadata<float> &ima_d, int iconnex) {
	imabin imaMaxLoc(nblig,nbcol), imavaldist(nblig,nbcol), imaRegion(nblig,nbcol), imaRes(nblig,nbcol);
	int i,j,i0,i2,j0,j2,ii,jj;
	float d;
	bool quatrecnx=0, imaxloc, imaxreg;
	if (iconnex==4) quatrecnx=1;
	for (i=0; i<nblig; i++) {
		i0=maxi(i-1,0); i2=mini(i+1,nblig-1);
		for (j=0; j<nbcol; j++) {
			d=ima_d(i,j);
			if (d>0) {
				j0=maxi(j-1,0); j2=mini(j+1,nbcol-1);
				imaxloc=1; ii=i0; jj=j0;
				while (imaxloc && ii<=i2 && jj<=j2) {
					if ((!quatrecnx || ii==i || jj==j) && d<ima_d(ii,jj)) imaxloc=0;
					else {
						if (jj<j2) jj++;
						else {jj=j0; ii++;}
					}
				}
				imaMaxLoc(i,j)=imaxloc;
			}
			else imaMaxLoc(i,j)=0;
		}
	}
	for (i=0; i<nblig; i++) {
		for (j=0; j<nbcol; j++) {
			if (imaMaxLoc(i,j)) {
				d=ima_d(i,j);
				for (ii=0; ii<nblig; ii++)
					for (jj=0; jj<nbcol; jj++)
						if (ima_d(ii,jj)==d) imavaldist(ii,jj)=1;
						else imavaldist(ii,jj)=0;
				imaRegion=imavaldist.reconstruction_geodesique(i,j,iconnex);
				imaxreg=1; ii=0, jj=0;
				while (imaxreg && ii<nblig && jj<nbcol) {
					if (imaRegion(ii,jj)) {
            i0=maxi(ii-1,0); i2=mini(ii+1,nblig-1);	j0=maxi(jj-1,0); j2=mini(jj+1,nbcol-1);
						d=ima_d(ii,jj);
						if (d<ima_d(i0,jj) || d<ima_d(i2,jj) || d<ima_d(ii,j0) || d<ima_d(ii,j2) ) imaxreg=0;
						if (imaxreg && !quatrecnx && (d<ima_d(i0,j0) || d<ima_d(i0,j2) || d<ima_d(i2,j0) || d<ima_d(i2,j2)) ) imaxreg=0;
					}
					jj++;
					if (jj==nbcol) {ii++; jj=0;}
				}
				for (ii=0; ii<nblig; ii++)
					for (jj=0; jj<nbcol; jj++)
						if (imaRegion(ii,jj)) {
							if (!imaxreg) imaRes(ii,jj)=0;
							else imaRes(ii,jj)=1;
							imaMaxLoc(ii,jj)=0;
						}
			}
		}
	}
	return imaRes;
}

imabin imabin::erode_ultime(int inoyau) {
	imadata<float> ima_d=negatif().Tr_dist(inoyau);
	int iconnex=8;
	if (inoyau==1) iconnex=4;
	imabin imaRes=erode_ultime(ima_d,iconnex);
//	imaRes.affiche();
	return imaRes;
}

imabin imabin::transfo_tout_ou_rien (eltstruct B1, eltstruct B2, int n) const { // on teste l'appartenance à B1 et la non appartenance à B2
	int i,j,i0,i2,j0,j2,ii,jj,iB,nn;
	int nl1=B1.nl, nc1=B1.nc, x01=B1.x0, y01=B1.y0, nl2=B2.nl, nc2=B2.nc, x02=B2.x0, y02=B2.y0;
//	cout<<" elt structurant 1 centre en "<<y01<<" "<<x01<<" elt structurant 2 centre en "<<y02<<" "<<x02<<"\n";
	bool z;
	imabin imaR(nblig,nbcol), imaC(*this);
	for (nn=0; nn<n; nn++) {
		for (i=0; i<nblig; i++) {
			i0=maxi(i-y01,0);
			i2=mini(i+nl1-1-y01,nblig-1);
			for (j=0; j<nbcol; j++) {
				j0=maxi(j-x01,0);
				j2=mini(j+nc1-1-x01,nbcol-1);
				z=1;
				for (ii=i0; ii<=i2; ii++) {
					if (z) {
						iB=(ii-i+y01)*nc1;
						for (jj=j0; jj<=j2; jj++)
							if (z && B1.aT[iB+jj-j+x01] && !imaC(ii,jj)) z=0;
					}
				}
				imaR(i,j)=z;
			}
		}
		for (i=0; i<nblig; i++) {
			i0=maxi(i-y02,0);
			i2=mini(i+nl2-1-y02,nblig-1);
			for (j=0; j<nbcol; j++) {
				j0=maxi(j-x02,0);
				j2=mini(j+nc2-1-x02,nbcol-1);
				z=imaR(i,j);
				for (ii=i0; ii<=i2; ii++) {
					if (z) {
						iB=(ii-i+y02)*nc2;
						for (jj=j0; jj<=j2; jj++)
							if (z && B2.aT[iB+jj-j+x02] && imaC(ii,jj)) z=0;
					}
				}
				imaR(i,j)=z;
			}
		}
	}
	return imaR;
}

imabin imabin::zones_influence_geodesique (const imabin &marq, bool &icheck, int iconnex) {
	imabin imaRes(marq), imaRes1(nblig+2,nbcol+2);
	imaRes1.mise_a_zero(); // imaRes1 est une image plus grande pour eviter les pbs de bords
	int i,j,i1,j1,iB,ncc,k,n,nE,ii,d_i=0,d_j=0,iii,jjj;
	for (i=0; i<nblig; i++) for (j=0; j<nbcol; j++) imaRes1(i+1,j+1)=marq(i,j);
	for (i=0; i<nblig+2; i++) { 
		imaRes1(i,0)=imaRes1(i,1); imaRes1(i,nbcol+1)=imaRes1(i,nbcol);
//		if (i>1 && imaRes1(i,1) && !imaRes1(i-1,1) && !imaRes1(i-2,1)) imaRes1(i-1,0)=imaRes1(i,1); // ??????????????????????
//		if (i<nblig && imaRes1(i,1) && !imaRes1(i+1,1) && !imaRes1(i+2,1)) imaRes1(i+1,0)=imaRes1(i,1); // ??????????????????
	}
	for (j=0; j<nbcol+2; j++) {
		imaRes1(0,j)=imaRes1(1,j); imaRes1(nblig+1,j)=imaRes1(nblig,j);
//		if (j>1 && imaRes1(1,j) && !imaRes1(1,j-1) && !imaRes1(1,j-2)) imaRes1(0,j-1)=imaRes1(1,j);
//		if (j<nbcol && imaRes1(1,j) && !imaRes1(1,j+1) && !imaRes1(1,j+2)) imaRes1(0,j+1)=imaRes1(1,j);
	}
	imaRes1(0,0)=imaRes1(1,1); imaRes1(nblig+1,nbcol+1)=imaRes1(nblig,nbcol); 
	imaRes1(0,nbcol+1)=imaRes1(1,nbcol); imaRes1(nblig+1,0)=imaRes1(nblig,1);
	if (icheck) {
		imadata<int> imacc=imaRes1.composantes_connexes(ncc,iconnex); // pour le cas où il y ait des composantes connexes
//		imacc.sauve_ImaPGM2("imacc_marqueurs.pgm");
		icheck=0;
		for (k=1; k<=ncc; k++) {
			n=0;
			for (i=0; i<nblig+2; i++) 
				for (j=0; j<nbcol+2; j++) if (imacc(i,j)==k) n++; // indices décalés de 1 pour imacc
			if (n<3) {
				icheck=1;
				int i0=nblig+2, i2=-1, j0=nbcol+2, j2=-1, ii, jj;
				for (i=0; i<nblig+2; i++) 
					for (j=0; j<nbcol+2; j++) if (imacc(i,j)==k) {if (i<i0) i0=i; if (i>i2) i2=i; if (j<j0) j0=j; if (j>j2) j2=j; }
				bool trouve=0; iB=0;
				do {
					switch (iB) {                          // éventuellement jusqu'à 4 directions à considérer successivement
						case 0: d_i=-1; d_j=0; break; case 1: d_i=0; d_j=-1; break; 
						case 2: d_i=+1; d_j=0; break; case 3: d_i=0; d_j=+1; break; }
					ii=maxi(i0,1); jj=maxi(j0,1);
					do {
						if (imacc(ii,jj)==k) {
							if (ii+d_i>0 && ii+d_i<nblig+1 && jj+d_j>0 && jj+d_j<nbcol+1 && !imaRes1(ii+d_i,jj+d_j) 
									&& (*this)(ii+d_i-1,jj+d_j-1)) { // indices décalés de 1 pour *this par rapport à imaRes1
								trouve=1;
								for (iii=maxi(ii+d_i-1,1); iii<=mini(ii+d_i+1,nblig); iii++)
									for (jjj=maxi(jj+d_j-1,1); jjj<=mini(jj+d_j+1,nbcol); jjj++) 
										if (imaRes1(iii,jjj) && imacc(iii,jjj)!=k) trouve=0;
								if (trouve) {
/*									for (iii=maxi(ii+d_i-1,1); iii<=mini(ii+d_i+1,nblig); iii++) {
										for (jjj=maxi(jj+d_j-1,1); jjj<=mini(jj+d_j+1,nbcol); jjj++) cout<<(int)imaRes1(iii,jjj)<<" ";
										cout<<"\n";}*/
									imaRes1(ii+d_i,jj+d_j)=1; imacc(ii+d_i,jj+d_j)=k; 
/*									cout<<" on ajoute le pixel ("<<ii+d_i-1<<","<<jj+d_j-1<<") aux marqueurs\n"; 
									for (iii=maxi(ii+d_i-1,1); iii<=mini(ii+d_i+1,nblig); iii++) {
										for (jjj=maxi(jj+d_j-1,1); jjj<=mini(jj+d_j+1,nbcol); jjj++) cout<<(int)imaRes1(iii,jjj)<<" ";
										cout<<"\n";}*/ //char aa; cin>>aa;
									n++; if (n<3) trouve=0; 
								}
							}
						}
						jj++; if (jj>j2) {ii++; jj=j0;}
					} while (!trouve && ii<=i2 && jj<=j2);
					iB++;
				} while (!trouve && iB<=4);
			}
		}
	}
	liste_pixels Lpix, L2, Lretrieve;
	elt_liste E;
	bool xv;
	for (i=0; i<nblig; i++)
//		for (j=0; j<nbcol; j++) if ((*this)(i,j) && !marq(i,j)) Lpix.insere(i,j);
		for (j=0; j<nbcol; j++) if ((*this)(i,j) && !imaRes1(i+1,j+1)) Lpix.insere(i,j); // indices décalés de 1 pour imaRes1
	cout<<" initialement "<<Lpix.nb_elts()<<" dans ZI, iconnex = "<<iconnex<<"\n";
	int nch=1;
	while (Lpix.nb_elts()>0 && nch>0) {
		nch=0;
		L2=Lpix; elt_liste *T_E=new elt_liste[L2.nb_elts()]; nE=0;
		while (L2.nb_elts()>0) {E=L2.extrait(); T_E[nE].x=E.x; T_E[nE].y=E.y; nE++;}
		for (iB=0; iB<4; iB++) {
			switch (iB) {                          // 4 directions à considérer successivement
				case 0: d_i=-1; d_j=0; break; case 1: d_i=0; d_j=-1; break; 
				case 2: d_i=+1; d_j=0; break; case 3: d_i=0; d_j=+1; break; }
			n=0; 
			while (n<nE) {
				i=T_E[n].x; j=T_E[n].y; n++; i1=i+1; j1=j+1; xv=imaRes1(i1+d_i,j1+d_j);
				if (xv) {
					ii=0;
					if (iconnex==4) {
						ii+=(int)(imaRes1(i1,j1+1)&&(!imaRes1(i1+1,j1+1)||!imaRes1(i1+1,j1)));
						ii+=(int)(imaRes1(i1+1,j1)&&(!imaRes1(i1+1,j1-1)||!imaRes1(i1,j1-1)));
						ii+=(int)(imaRes1(i1,j1-1)&&(!imaRes1(i1-1,j1-1)||!imaRes1(i1-1,j1)));
						ii+=(int)(imaRes1(i1-1,j1)&&(!imaRes1(i1-1,j1+1)||!imaRes1(i1,j1+1)));
					} 
					if (iconnex==8) {
						ii+=(int)(!imaRes1(i1,j1+1)&&(imaRes1(i1+1,j1+1)||imaRes1(i1+1,j1)));
						ii+=(int)(!imaRes1(i1+1,j1)&&(imaRes1(i1+1,j1-1)||imaRes1(i1,j1-1)));
						ii+=(int)(!imaRes1(i1,j1-1)&&(imaRes1(i1-1,j1-1)||imaRes1(i1-1,j1)));
						ii+=(int)(!imaRes1(i1-1,j1)&&(imaRes1(i1-1,j1+1)||imaRes1(i1,j1+1)));
					}
					if (ii==1) Lretrieve.insere(i,j);
				}
			}
			while (Lretrieve.nb_elts()>0) {
				E=Lretrieve.extrait();
				imaRes1(E.x+1,E.y+1)=1;
				Lpix.supprime(E.x,E.y);
				nch++;
			}
		}
		if (T_E!=NULL) delete[] T_E;
/*		for (iB=0; iB<8; iB++) {                           // 8 masques à utiliser successivement
			n=0; //cout<<" iB = "<<iB<<"\n";
			while (n<nE) {
				i=T_E[n].x; j=T_E[n].y; n++; epais=0; xv=0; i1=i+1; j1=j+1;
				if (iB==0||iB==1) xv=imaRes1(i1-1,j1);
				else if (iB==2||iB==3) xv=imaRes1(i1,j1-1);
					 else if (iB==4||iB==5) xv=imaRes1(i1+1,j1);
						  else if (iB==6||iB==7) xv=imaRes1(i1,j1+1);
				if (xv) {
					if (iconnex==4) {
						i1=i+1; j1=j+1;
						switch (iB) {
							case 0: if (!imaRes1(i1+1,j1-1)&&!imaRes1(i1+1,j1)&&!imaRes1(i1+1,j1+1)&&
												imaRes1(i1-1,j1-1)&& imaRes1(i1-1,j1)&& imaRes1(i1-1,j1+1)) epais=1; break;
							case 1: if (!imaRes1(i1,j1-1)&&!imaRes1(i1+1,j1)&&!imaRes1(i1+1,j1-1)&&
												imaRes1(i1,j1+1)&&imaRes1(i1-1,j1)&&imaRes1(i1-1,j1+1)) epais=1; break;
							case 2: if (!imaRes1(i1-1,j1+1)&&!imaRes1(i1,j1+1)&&!imaRes1(i1+1,j1+1)&&
												imaRes1(i1-1,j1-1)&& imaRes1(i1,j1-1)&& imaRes1(i1+1,j1-1)) epais=1; break;
							case 3: if (!imaRes1(i1,j1+1)&&!imaRes1(i1+1,j1)&&!imaRes1(i1+1,j1+1)&&
												imaRes1(i1,j1-1)&&imaRes1(i1-1,j1)&&imaRes1(i1-1,j1-1)) epais=1; break;
							case 4: if (!imaRes1(i1-1,j1-1)&&!imaRes1(i1-1,j1)&&!imaRes1(i1-1,j1+1)&&
												imaRes1(i1+1,j1-1)&& imaRes1(i1+1,j1)&& imaRes1(i1+1,j1+1)) epais=1; break;
							case 5: if (!imaRes1(i1,j1+1)&&!imaRes1(i1-1,j1)&&!imaRes1(i1-1,j1+1)&&
												imaRes1(i1,j1-1)&&imaRes1(i1+1,j1)&&imaRes1(i1+1,j1-1)) epais=1; break;
							case 6: if (!imaRes1(i1-1,j1-1)&&!imaRes1(i1,j1-1)&&!imaRes1(i1+1,j1-1)&&
												imaRes1(i1-1,j1+1)&& imaRes1(i1,j1+1)&& imaRes1(i1+1,j1+1)) epais=1; break;
							case 7: if (!imaRes1(i1,j1-1)&&!imaRes1(i1-1,j1)&&!imaRes1(i1-1,j1-1)&&
												imaRes1(i1,j1+1)&&imaRes1(i1+1,j1)&&imaRes1(i1+1,j1+1)) epais=1; break;
						}
					}
					if (iconnex==8) {
						i1=i+1; j1=j+1;
						switch (iB) {
							case 0: if (!imaRes1(i1+1,j1-1)&&!imaRes1(i1+1,j1)&&!imaRes1(i1+1,j1+1)&&imaRes1(i1-1,j1)) epais=1; break;
							case 1: if (!imaRes1(i1+1,j1-1)&&imaRes1(i1,j1+1)&& imaRes1(i1-1,j1)) epais=1; break;
							case 2: if (!imaRes1(i1-1,j1+1)&&!imaRes1(i1,j1+1)&&!imaRes1(i1+1,j1+1)&&imaRes1(i1,j1-1)) epais=1; break;
							case 3: if (!imaRes1(i1+1,j1+1)&&imaRes1(i1,j1-1)&& imaRes1(i1-1,j1)) epais=1; break;
							case 4: if (!imaRes1(i1-1,j1-1)&&!imaRes1(i1-1,j1)&&!imaRes1(i1-1,j1+1)&&imaRes1(i1+1,j1)) epais=1; break;
							case 5: if (!imaRes1(i1-1,j1+1)&&imaRes1(i1,j1-1)&& imaRes1(i1+1,j1)) epais=1; break;
							case 6: if (!imaRes1(i1-1,j1-1)&&!imaRes1(i1,j1-1)&&!imaRes1(i1+1,j1-1)&&imaRes1(i1,j1+1)) epais=1; break;
							case 7: if (!imaRes1(i1-1,j1-1)&&imaRes1(i1,j1+1)&& imaRes1(i1+1,j1)) epais=1; break;
						}
					}
				}
				if (epais) Lretrieve.insere(i,j);
			}
			while (Lretrieve.nb_elts()>0) {
				E=Lretrieve.extrait();
				imaRes1(E.x+1,E.y+1)=1;
				Lpix.supprime(E.x,E.y);
				nch++;
			}
		}
		if (T_E!=NULL) delete[] T_E;*/
	}
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) imaRes(i,j)=imaRes1(i+1,j+1);
	return imaRes;
}

imadata<BYTE> imabin::detect_coin (int l, bool c_s, bool c_e) {
	const int noUL=1, noUR=2, noLL=3, noLR=4, moUL=5, moUR=6, moLL=7, moLR=8;
	eltstruct B1_UL(l+1,l+1,1,1), B2_UL(l+1,l+1,1,1), B1_UR(l+1,l+1,1,l-1), B2_UR(l+1,l+1,1,l-1),
			  B1_LL(l+1,l+1,l-1,1), B2_LL(l+1,l+1,l-1,1), B1_LR(l+1,l+1,l-1,l-1), B2_LR(l+1,l+1,l-1,l-1);
	int i,j,nUL=0,nUR=0,nLL=0,nLR=0,n=0,mUL=0,mUR=0,mLL=0,mLR=0,m=0;
	imadata<BYTE> imaRes(nblig,nbcol);
	imabin imaCoin(nblig,nbcol);

/************************** debut a programmer ****************************/




/***************************  fin a programmer ****************************/

	cout<<n<<" coin(s) saillant(s) = ";
	cout<<nUL<<" coin(s) UL + "<<nUR<<" coin(s) UR + "<<nLL<<" coin(s) LL + "<<nLR<<" coin(s) LR\n";
	cout<<m<<" coin(s) entrant(s) = ";
	cout<<mUL<<" coin(s) UL + "<<mUR<<" coin(s) UR + "<<mLL<<" coin(s) LL + "<<mLR<<" coin(s) LR\n";
	return imaRes;
}

/*imabin imabin::elagage (int nb_it) {
	imabin imaRes(*this);
	liste_pixels Lpix, L2, Lretrieve;
	elt_liste E;
	int i,j;
	bool iretrieve;
	for (i=1; i<nblig-1; i++)                                       // on ne traite pas les bords
		for (j=1; j<nbcol-1; j++) if ((*this)(i,j)) Lpix.insere(i,j);
	int nch=1, it=0;
	while (Lpix.nb_elts()>0 && nch>0 && it<nb_it) {
		nch=0;
		for (int iB=0; iB<8; iB++) {                           // 8 masques a utiliser successivement
			L2=Lpix;
			while (L2.nb_elts()>0) {
				E=L2.extrait(); i=E.x; j=E.y; iretrieve=0;
				switch (iB) {
					case 0: iretrieve=imaRes(i-1,j)&&!imaRes(i,j-1)&&!imaRes(i,j+1)&&!imaRes(i+1,j-1)&&!imaRes(i+1,j)&&!imaRes(i+1,j+1); break;
					case 1: iretrieve=imaRes(i+1,j)&&!imaRes(i,j-1)&&!imaRes(i,j+1)&&!imaRes(i-1,j-1)&&!imaRes(i-1,j)&&!imaRes(i-1,j+1); break;
					case 2: iretrieve=imaRes(i,j-1)&&!imaRes(i-1,j)&&!imaRes(i+1,j)&&!imaRes(i-1,j+1)&&!imaRes(i,j+1)&&!imaRes(i+1,j+1); break;
					case 3: iretrieve=imaRes(i,j+1)&&!imaRes(i-1,j)&&!imaRes(i+1,j)&&!imaRes(i-1,j-1)&&!imaRes(i,j-1)&&!imaRes(i+1,j-1); break;
					case 4: iretrieve=imaRes(i-1,j-1)&&!imaRes(i-1,j+1)&&!imaRes(i,j+1)&&!imaRes(i+1,j-1)&&!imaRes(i+1,j)&&!imaRes(i+1,j+1); break;
					case 5: iretrieve=imaRes(i-1,j+1)&&!imaRes(i-1,j-1)&&!imaRes(i,j-1)&&!imaRes(i+1,j-1)&&!imaRes(i+1,j)&&!imaRes(i+1,j+1); break;
					case 6: iretrieve=imaRes(i+1,j-1)&&!imaRes(i-1,j-1)&&!imaRes(i-1,j)&&!imaRes(i-1,j+1)&&!imaRes(i,j+1)&&!imaRes(i+1,j+1); break;
					case 7: iretrieve=imaRes(i+1,j+1)&&!imaRes(i-1,j-1)&&!imaRes(i-1,j)&&!imaRes(i-1,j+1)&&!imaRes(i,j-1)&&!imaRes(i+1,j-1); break;
				}
				if (iretrieve) Lretrieve.insere(i,j);
			}
			while (Lretrieve.nb_elts()>0) {
				E=Lretrieve.extrait();
				imaRes(E.x,E.y)=0;
				Lpix.supprime(E.x,E.y);
				nch++;
			}
		}
		it++;
	}
	return imaRes;
}*/

imabin imabin::elagage (int nb_it, bool itriple) {
	imabin imaRes=ajoutbords(), imaMarq; //imaRes.imaunsignedchar().sauve_ImaPGM("ima_a_elaguer.pgm");
	liste_pixels Lpix, L2, Lretrieve;
	elt_liste E;
	int i,j;
	bool iretrieve;
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) if ((*this)(i,j)) Lpix.insere(i+1,j+1); // +1 sur coordonnées à cause des bords ajoutés !
	int nch=1, it=0;
	if (itriple) {
		imaMarq=imabin(nblig+2,nbcol+2);
		int ii,jj,i0,i2,j0,j2,n;
		for (i=1; i<nblig; i++)					// detection des points triples
			for (j=1; j<nbcol; j++) {
				n=0; i0=maxi(0,i-1); i2=mini(i+1,nblig-1); j0=maxi(0,j-1); j2=mini(j+1,nbcol-1); 
				for (ii=i0; ii<=i2; ii++) for (jj=j0; jj<=j2; jj++) n+=(int)imaRes(ii,jj);
				if (n>3) imaMarq(i,j)=1;
			}
	}
	while (Lpix.nb_elts()>0 && nch>0 && it<nb_it) {
		nch=0;
		for (int iB=0; iB<8; iB++) {                           // 8 masques a utiliser successivement
			L2=Lpix;
			while (L2.nb_elts()>0) {
				E=L2.extrait(); i=E.x; j=E.y; iretrieve=0;
				if (!itriple || !imaMarq(i,j))
					switch (iB) {
						case 0: iretrieve=imaRes(i-1,j)&&!imaRes(i,j-1)&&!imaRes(i,j+1)&&!imaRes(i+1,j-1)&&!imaRes(i+1,j)&&!imaRes(i+1,j+1); break;
						case 1: iretrieve=imaRes(i+1,j)&&!imaRes(i,j-1)&&!imaRes(i,j+1)&&!imaRes(i-1,j-1)&&!imaRes(i-1,j)&&!imaRes(i-1,j+1); break;
						case 2: iretrieve=imaRes(i,j-1)&&!imaRes(i-1,j)&&!imaRes(i+1,j)&&!imaRes(i-1,j+1)&&!imaRes(i,j+1)&&!imaRes(i+1,j+1); break;
						case 3: iretrieve=imaRes(i,j+1)&&!imaRes(i-1,j)&&!imaRes(i+1,j)&&!imaRes(i-1,j-1)&&!imaRes(i,j-1)&&!imaRes(i+1,j-1); break;
						case 4: iretrieve=imaRes(i-1,j-1)&&!imaRes(i-1,j+1)&&!imaRes(i,j+1)&&!imaRes(i+1,j-1)&&!imaRes(i+1,j)&&!imaRes(i+1,j+1); break;
						case 5: iretrieve=imaRes(i-1,j+1)&&!imaRes(i-1,j-1)&&!imaRes(i,j-1)&&!imaRes(i+1,j-1)&&!imaRes(i+1,j)&&!imaRes(i+1,j+1); break;
						case 6: iretrieve=imaRes(i+1,j-1)&&!imaRes(i-1,j-1)&&!imaRes(i-1,j)&&!imaRes(i-1,j+1)&&!imaRes(i,j+1)&&!imaRes(i+1,j+1); break;
						case 7: iretrieve=imaRes(i+1,j+1)&&!imaRes(i-1,j-1)&&!imaRes(i-1,j)&&!imaRes(i-1,j+1)&&!imaRes(i,j-1)&&!imaRes(i+1,j-1); break;
					}
				if (iretrieve) Lretrieve.insere(i,j);
			}
			while (Lretrieve.nb_elts()>0) {
				E=Lretrieve.extrait();
				imaRes(E.x,E.y)=0;
				Lpix.supprime(E.x,E.y);
				nch++;
			}
		}
		it++;
	} //imaRes.imaunsignedchar().sauve_ImaPGM("ima_elaguee.pgm");
	return imaRes.retirebords();
}

imabin imabin::elagage_saufextrem (int nb_it, int iconnex) {
// attention : fonction à tester !!!!!!!!! pourquoi ii dans imaRes.elagage(ii,1) ?
	int ncc, i,j,k,ii;
	imabin imaRes(*this), imaDif; 
	for (ii=2; ii<=nb_it; ++ii) {
//	for (ii=nb_it; ii<=nb_it; ++ii) {
		imaDif=imaRes-imaRes.elagage(ii,1); //imaDif.imaunsignedchar().sauve_ImaPGM("./imaDif.pgm");
		imadata<int> imaDif_cc=imaDif.composantes_connexes(ncc,iconnex);
		int *np_cc=new int[ncc]; for (k=0; k<ncc; k++) np_cc[k]=0;
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				k=imaDif_cc(i,j)-1; 
				if (k>=0 && k<ncc) np_cc[k]++; 
			}
//		cout<<" elagage "<<ii<<" : \n"; for (k=0; k<ncc; k++) cout<<" comp. "<<k<<" nbpix = "<<np_cc[k]<<"; \n"; cout<<"\n";
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++) {
				k=imaDif_cc(i,j)-1;
				if (k>=0 && k<ncc && np_cc[k]<ii) imaRes(i,j)=0;
			}
			if (np_cc!=NULL) delete[] np_cc; //{imaRes.imaunsignedchar().sauve_ImaPGM("./imaRes.pgm"); char aa; cin>>aa;}
	}
	return imaRes;
}

imabin imabin::detect_extremites_sq () {// s'applique à une image de squelette
// detection grossiere des extremites par différence avec élagage de 1px
	imabin ima_extr=(*this)-(*this).elagage(1); //ima_extr.imaunsignedchar().sauve_ImaPGM("./ima_extr1px.pgm");
	int i,j,k,ncc_e;
	imadata<int> ima_extr_cc=ima_extr.composantes_connexes(ncc_e,8); //ima_extr_cc.imaunsignedchar().sauve_ImaPGM("./ima_extr1px_cc.pgm");
	int **tab_cc_e=new int*[ncc_e+1]; for (k=0; k<=ncc_e; k++) tab_cc_e[k]=new int[3]; // pour chaque extremité: nbre de points, coordonnées du dernier point
	for (k=0; k<=ncc_e; k++) for (j=0; j<3; j++) tab_cc_e[k][j]=0; 
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) {k=ima_extr_cc(i,j); if (k>0) {tab_cc_e[k][0]++; tab_cc_e[k][1]=i; tab_cc_e[k][2]=j;} }
// raffinement des extremités non ponctuelles à 1px
	float xk,yk,nk; int ik,jk,dmin,l,ii,jj;
	for (k=1; k<=ncc_e; k++) 
		if (tab_cc_e[k][0]>=2) { //si extremité non ponctuelle, extremité approximée par point inclu dans l'extremité non ponctuelle et le + proche du centroide
			xk=yk=nk=0.f; i=tab_cc_e[k][1]; j=tab_cc_e[k][2]; l=(tab_cc_e[k][0]+1)/2;
			for (ii=maxi(0,i-l); ii<=mini(nblig-1,i+l); ii++)
				for (jj=maxi(0,j-l); jj<=mini(nbcol-1,j+l); jj++)
					if (ima_extr(ii,jj) && ima_extr_cc(ii,jj)==k) {xk+=ii; yk+=jj; nk+=1;}
			xk/=nk; yk/=nk;
			ik=(int)xk; jk=(int)yk; if (ima_extr(ik,jk) && ima_extr_cc(ik,jk)==k) dmin=0; else dmin=INT_MAX; 
			for (ii=maxi(0,i-l); ii<=mini(nblig-1,i+l); ii++)
				for (jj=maxi(0,j-l); jj<=mini(nbcol-1,j+l); jj++) {
					if (ima_extr(ii,jj) && ima_extr_cc(ii,jj)==k) {
						ima_extr(ii,jj)=0; 
						if (dmin>0 && abs(ii-ik)+abs(jj-jk)<dmin) {dmin=abs(ii-ik)+abs(jj-jk); ik=ii; jk=jj;}
					}
				}
			ima_extr(ik,jk)=1; 
		}
	for (k=0; k<=ncc_e; k++) if (tab_cc_e[k]!=NULL) delete[] tab_cc_e[k]; if (tab_cc_e!=NULL) delete[] tab_cc_e;
// élimination des extrémités trop proches du squelette principal
	imabin ima_extr2=(*this)-(*this).elagage(2); //ima_extr2.imaunsignedchar().sauve_ImaPGM("./ima_extr2px.pgm");
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++)
			if (ima_extr2(i,j)) {
				k=0; 
				for (ii=maxi(0,i-1); ii<=mini(nblig-1,i+1); ii++)
					for (jj=maxi(0,j-1); jj<=mini(nbcol-1,j+1); jj++) {
						k+=(int)ima_extr2(ii,jj);
					}
				if (k<2) ima_extr2(i,j)=0; 
			}
	ima_extr=ima_extr&&ima_extr2; ima_extr.imaunsignedchar().sauve_ImaPGM("./ima_extr1px_b.pgm");
	cout<<" nb extremites = "<<ima_extr.norm()<<"\n";
	return ima_extr;
}

imadata<BYTE> imabin::calc_orient_sq (bool ireg, bool iverif) {cout<<" entree dans calc_orient_sq\n"; // s'applique à une image de squelette
	imabin ima_sq8(*this);
	if (iverif) {ima_sq8=squelette(8,2); /*cout<<" ima_sq8 calculee : "<<ima_sq8.norm()-norm()<<"\n"; ima_sq8.imaunsignedchar().sauve_ImaPGM("ima_sq8.pgm");*/}
	imadata<BYTE> ima_orient(nblig,nbcol);
	const int sz_tab_o=256; const float res_o_deg=22.5f;
	int tab_o[sz_tab_o], i,j,k; // table codant en dur les orientations (0,22.5,45,67.5,90,112.5,135,158.5) selon la config du voisinage
	for (k=0; k<sz_tab_o; k++) tab_o[k]=-999;
	tab_o[1]=tab_o[16]=tab_o[17]=(int)(0.f/res_o_deg); tab_o[4]=tab_o[64]=tab_o[68]=(int)(90.f/res_o_deg); 
	tab_o[2]=tab_o[32]=tab_o[34]=(int)(45.f/res_o_deg); tab_o[8]=tab_o[128]=tab_o[136]=(int)(135.f/res_o_deg); 
	tab_o[32+1]=tab_o[16+2]=(int)(22.5f/res_o_deg); tab_o[32+4]=tab_o[64+2]=(int)(67.5f/res_o_deg); 
	tab_o[64+8]=tab_o[128+4]=(int)(112.5f/res_o_deg); tab_o[1+8]=tab_o[128+16]=(int)(157.5f/res_o_deg); 
	tab_o[2+8]=tab_o[128+32]=(int)(0.f/res_o_deg); tab_o[2+128]=tab_o[8+32]=(int)(90.f/res_o_deg);

	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++)
			if (ima_sq8(i,j)) {
				k=0; if (i>0 && ima_sq8(i-1,j)) k+=4; if (i+1<nblig && ima_sq8(i+1,j)) k+=64;
				if (j+1<nbcol) {if (ima_sq8(i,j+1)) k+=1; if (i>0 && ima_sq8(i-1,j+1)) k+=2; if (i+1<nblig && ima_sq8(i+1,j+1)) k+=128;}
				if (j>0) {if (ima_sq8(i,j-1)) k+=16; if (i>0 && ima_sq8(i-1,j-1)) k+=8; if (i+1<nblig && ima_sq8(i+1,j-1)) k+=32;}
//				cout<<i<<" "<<j<<" k "<<k<<" tab_o[k] "<<tab_o[k]<<"\n";
				if (tab_o[k]>=0) ima_orient(i,j)=tab_o[k]+1; else ima_orient(i,j)=255;
			}
	ima_orient.sauve_ImaPGM("ima_orient.pgm"); cout<<" sauve ima_orient.pgm\n"; 
	if (ireg) {cout<<" etape de regularisation\n";
		eltstruct S3(3,3);
		imabin ima_unknown(ima_orient,255), ima_unknown_dil=ima_unknown.dilate(S3);
		int ncc,ncc2,n,iorient; bool iaff=1;
//		imadata<int> ima_unknown_dil_cc=ima_unknown_dil.composantes_connexes(ncc2,8,iaff), ima_unknown_cc=ima_unknown.composantes_connexes(ncc,8,iaff);
		imadata<int> ima_unknown_dil_cc=ima_unknown_dil.CompoConnexes_MM(ncc2,8,iaff), ima_unknown_cc=ima_unknown.CompoConnexes_MM(ncc,8,iaff);
		ima_unknown_dil_cc.imaunsignedchar().sauve_ImaPGM("./ima_unknown_dil_cc.pgm"); ima_unknown_cc.imaunsignedchar().sauve_ImaPGM("./ima_unknown_cc.pgm");
		cout<<" ncc = "<<ncc<<", ncc2 = "<<ncc2<<"\n";
		int **tab=new int*[ncc2]; for (k=0; k<ncc2; k++) {tab[k]=new int[10]; for (i=1; i<10; i++) tab[k][i]=0; }
		for (i=0; i<nblig; i++) 
			for (j=0; j<nbcol; j++) {
				if (ima_unknown_dil(i,j) && !ima_unknown(i,j)) {
					k=ima_unknown_dil_cc(i,j)-1;
					if (ima_orient(i,j)>0 && ima_orient(i,j)<9) tab[k][ima_orient(i,j)]++;
				}
		}
//		cout<<" ***********\n ";
		for (k=0; k<ncc2; k++) {
			tab[k][9]=iorient=-1; n=0;
			for (i=1; i<9; i++) if (tab[k][i]>n) {n=tab[k][i]; iorient=i;}
			tab[k][9]=iorient; //cout<<" cc "<<k+1<<" orientation "<<tab[k][9]<<"\n";
		}
//		cout<<" %%%%%%%%%%%%\n ";
		for (i=0; i<nblig; i++)
			for (j=0; j<nbcol; j++)
				if (ima_unknown(i,j)) {
					k=ima_unknown_dil_cc(i,j)-1; if (k>=0 && k<ncc2) iorient=tab[k][9];
					if (iorient>=0 && iorient<9) ima_orient(i,j)=iorient; //cout<<" pixel ("<<i<<","<<j<<") recoit orientation "<<ima_orient(i,j)<<"\n";
				}
//		cout<<" ;;;;;;;;;;;;\n ";

		for (k=0; k<ncc2; k++) if (tab[k]!=NULL) delete[] tab[k]; if (tab!=NULL) delete[] tab;
		ima_orient.sauve_ImaPGM("ima_orient2.pgm"); cout<<" sauve ima_orient2.pgm\n"; 
	}cout<<" fin calc_orient_sq\n"; 
	return ima_orient;
}

imabin imabin::enveloppe_convexe () {
	imabin imaRes(*this);
	liste_pixels Lpix, Lretrieve, L2;
	elt_liste E;
	int i,j;
	bool epais;
	for (i=1; i<nblig-1; i++)                                       // on ne traite pas les bords
		for (j=1; j<nbcol-1; j++) if (!(*this)(i,j)) Lpix.insere(i,j);
	cout<<" initialement "<<Lpix.nb_elts()<<" pixels a traiter\n";
	int nch=1, nb_iB=12, it=0;
	while (Lpix.nb_elts()>0 && nch>0 && it<1000) { it++;
		nch=0;
		for (int iB=0; iB<nb_iB; iB++) {                           // 12 masques a utiliser successivement en 4-connexite
			L2=Lpix;
			while (L2.nb_elts()>0) {
				E=L2.extrait(); i=E.x; j=E.y; epais=0;
				switch (iB) {
					case 0: epais=!imaRes(i-1,j)&&imaRes(i,j-1)&&imaRes(i,j+1); break;
					case 1: epais=!imaRes(i+1,j)&&imaRes(i,j-1)&&imaRes(i,j+1); break;
					case 2: epais=!imaRes(i,j-1)&&imaRes(i-1,j)&&imaRes(i+1,j); break;
					case 3: epais=!imaRes(i,j+1)&&imaRes(i-1,j)&&imaRes(i+1,j); break;
					case 4: epais=imaRes(i+1,j)&&imaRes(i,j-1)&&imaRes(i-1,j-1); break;
					case 5: epais=imaRes(i+1,j)&&imaRes(i,j+1)&&imaRes(i-1,j+1); break;
					case 6: epais=imaRes(i,j-1)&&imaRes(i-1,j)&&imaRes(i-1,j+1); break;
					case 7: epais=imaRes(i,j-1)&&imaRes(i+1,j)&&imaRes(i+1,j+1); break;
					case 8: epais=imaRes(i-1,j)&&imaRes(i,j+1)&&imaRes(i+1,j+1); break;
					case 9: epais=imaRes(i-1,j)&&imaRes(i,j-1)&&imaRes(i+1,j-1); break;
					case 10: epais=imaRes(i,j+1)&&imaRes(i+1,j)&&imaRes(i+1,j-1); break;
					case 11: epais=imaRes(i,j+1)&&imaRes(i-1,j)&&imaRes(i-1,j-1); break;
				}
				if (epais) Lretrieve.insere(i,j);
			}
			while (Lretrieve.nb_elts()>0) {
				E=Lretrieve.extrait();
				imaRes(E.x,E.y)=1;
				Lpix.supprime(E.x,E.y);
				nch++;
			}
		}
		cout<<" iteration "<<it<<" : nch "<<nch<<" sur "<<Lpix.nb_elts()<<"\n";
		imaRes.imaunsignedchar().sauve_Ima("env_convexe.dat");
	}
	return imaRes;
}

imabin imabin::squelette (int iconnex, BYTE imarq, bool iaf) {
	imabin marq(nblig,nbcol), imaRes(nblig,nbcol);
	liste_pixels Lpix, L2, Lretrieve;
	elt_liste E;
	int n,i,j,i0,i2,j0,j2,icard;
	bool xv,x0,x1,x2,x3,x4,x5,x6,x7,elimine;
	char nomfich[14]="squelette_it0";
	imaRes=(*this);
	switch (imarq) {
		case 0 : {
			if (iaf) cout<<" points d'ancrage = erode ultime\n";
			marq=erode_ultime();
			break; }
		case 1 : {
			if (iaf) cout<<" points d'ancrage = maxima locaux de la transformee en distance\n";
			int inoyau=1; if (iconnex==8) inoyau=2;
			imabin imatemp=negatif();
			imadata<float> imdist=imatemp.Tr_dist(inoyau);
//			imdist.imaunsignedchar().sauve_ImaPGM("imdist.pgm");
			for (i=1; i<nblig-1; i++)                                       // on ne traite pas les bords
				for (j=1; j<nbcol-1; j++)
					if ((imdist(i,j)>0) && (imdist(i,j)>=imdist(i-1,j)) && (imdist(i,j)>=imdist(i+1,j)) &&
						(imdist(i,j)>=imdist(i,j-1)) && (imdist(i,j)>=imdist(i,j+1))) marq(i,j)=1;
			if (iconnex==8)
				for (i=1; i<nblig-1; i++)                                       // on ne traite pas les bords
					for (j=1; j<nbcol-1; j++)
						if (marq(i,j) && (imdist(i,j)<imdist(i-1,j-1) || imdist(i,j)<imdist(i-1,j+1) ||
							 imdist(i,j)<imdist(i+1,j-1) || imdist(i,j)<imdist(i+1,j+1))) marq(i,j)=0;
			break; }
		case 2 : {
			if (iaf) cout<<" points d'ancrage = maxima locaux de la transformee en distance puis affinement\n";
			imabin imatemp=negatif();
			imadata<float> imdist=imatemp.Tr_dist(4);
			for (i=1; i<nblig-1; i++)                                       // on ne traite pas les bords
				for (j=1; j<nbcol-1; j++)
					if ((imdist(i,j)>0) && (imdist(i,j)>=imdist(i-1,j)) && (imdist(i,j)>=imdist(i+1,j)) &&
						(imdist(i,j)>=imdist(i,j-1)) && (imdist(i,j)>=imdist(i,j+1))) marq(i,j)=1;
			int ii,jj,i0,i2,j0,j2,n;
			for (i=1; i<nblig-1; i++)																				// on elimine les points triple des points d'ancrage
				for (j=1; j<nbcol-1; j++) {
					n=0; i0=maxi(0,i-1); i2=mini(i+1,nblig-1); j0=maxi(0,j-1); j2=mini(j+1,nbcol-1); 
					for (ii=i0; ii<=i2; ii++) for (jj=j0; jj<=j2; jj++) n+=(int)marq(ii,jj);
					if (n>3) marq(i,j)=0;
				}
			break; }
		default : if (iaf) cout<<" pas de point d'ancrage\n";
	}
	marq.imaunsignedchar().sauve_ImaPGM("points_ancrage.pgm");
	for (i=1; i<nblig-1; i++)                                       // on ne traite pas les bords
		for (j=1; j<nbcol-1; j++) if ((*this)(i,j)) Lpix.insere(i,j);
	int nch=1;
	while (Lpix.nb_elts()>0 && nch>0) {
		nch=0;
		for (icard=0; icard<4; icard++) { // 4 cas de voisins dans Xc : Nord, Est, Sud, Ouest
			L2=Lpix;
			while (L2.nb_elts()>0) {
				E=L2.extrait(); i=E.x; j=E.y; elimine=0;
				if (!marq(i,j)) {
					if (icard==0) xv=imaRes(i-1,j);
					else if (icard==1) xv=imaRes(i,j+1);
						else if (icard==2) xv=imaRes(i+1,j);
							else if (icard==3) xv=imaRes(i,j-1);
					if (!xv) {
						i0=i-1; j0=j-1; i2=i+1; j2=j+1;
						x0=imaRes(i,j2); x1=imaRes(i0,j2); x2=imaRes(i0,j); x3=imaRes(i0,j0);
						x4=imaRes(i,j0); x5=imaRes(i2,j0); x6=imaRes(i2,j); x7=imaRes(i2,j2);
						n=0;
						if (iconnex==8) {
							if (!x0 && (x1 || x2)) n++; if (!x2 && (x3 || x4)) n++;
							if (!x4 && (x5 || x6)) n++;	if (!x6 && (x7 || x0)) n++;
						} else {
							if (x0 && (!x1 || !x2)) n++; if (x2 && (!x3 || !x4)) n++;
							if (x4 && (!x5 || !x6)) n++; if (x6 && (!x7 || !x0)) n++;
						}
						if (n==1) elimine=1;
					}
				}
				if (elimine) Lretrieve.insere(i,j);
			}
			while (Lretrieve.nb_elts()>0) {
				E=Lretrieve.extrait(); imaRes(E.x,E.y)=0;
				Lpix.supprime(E.x,E.y);                                // a verifier
				nch++;                                                 // a verifier
			}
		}
	}
/*	bool iretrieve;                                                // elagage
	while (Lpix.nb_elts()>0) {E=Lpix.extrait();	Lpix.supprime(E.x,E.y);}
	for (i=1; i<nblig-1; i++)                                      // on ne traite pas les bords ni les points d'ancrage
		for (j=1; j<nbcol-1; j++) if (imaRes(i,j)&&!marq(i,j)) Lpix.insere(i,j);
	cout<<" # de pixels du squelette hors points d'ancrage = "<<Lpix.nb_elts()<<"\n";
	int it=0, nb_it=50;
	nch=1;
	while (Lpix.nb_elts()>0 && nch>0 && it<nb_it) {
		nch=0;
		for (int iB=0; iB<8; iB++) {                           // 8 masques a utiliser successivement
			L2=Lpix;
			while (L2.nb_elts()>0) {
				E=L2.extrait(); i=E.x; j=E.y; iretrieve=0;
				if (marq(i,j)) cout<<" it "<<it<<" point d'ancrage\n";
				switch (iB) {
					case 0: iretrieve=imaRes(i-1,j)&&!imaRes(i,j-1)&&!imaRes(i,j+1)&&!imaRes(i+1,j-1)&&!imaRes(i+1,j)&&!imaRes(i+1,j+1); break;
					case 1: iretrieve=imaRes(i+1,j)&&!imaRes(i,j-1)&&!imaRes(i,j+1)&&!imaRes(i-1,j-1)&&!imaRes(i-1,j)&&!imaRes(i-1,j+1); break;
					case 2: iretrieve=imaRes(i,j-1)&&!imaRes(i-1,j)&&!imaRes(i+1,j)&&!imaRes(i-1,j+1)&&!imaRes(i,j+1)&&!imaRes(i+1,j+1); break;
					case 3: iretrieve=imaRes(i,j+1)&&!imaRes(i-1,j)&&!imaRes(i+1,j)&&!imaRes(i-1,j-1)&&!imaRes(i,j-1)&&!imaRes(i+1,j-1); break;
					case 4: iretrieve=imaRes(i-1,j-1)&&!imaRes(i-1,j+1)&&!imaRes(i,j+1)&&!imaRes(i+1,j-1)&&!imaRes(i+1,j)&&!imaRes(i+1,j+1); break;
					case 5: iretrieve=imaRes(i-1,j+1)&&!imaRes(i-1,j-1)&&!imaRes(i,j-1)&&!imaRes(i+1,j-1)&&!imaRes(i+1,j)&&!imaRes(i+1,j+1); break;
					case 6: iretrieve=imaRes(i+1,j-1)&&!imaRes(i-1,j-1)&&!imaRes(i-1,j)&&!imaRes(i-1,j+1)&&!imaRes(i,j+1)&&!imaRes(i+1,j+1); break;
					case 7: iretrieve=imaRes(i+1,j+1)&&!imaRes(i-1,j-1)&&!imaRes(i-1,j)&&!imaRes(i-1,j+1)&&!imaRes(i,j-1)&&!imaRes(i+1,j-1); break;
				}
				if (iretrieve) Lretrieve.insere(i,j);
			}
			while (Lretrieve.nb_elts()>0) {
				E=Lretrieve.extrait();
				imaRes(E.x,E.y)=0;
				Lpix.supprime(E.x,E.y);
				nch++;
			}
		}
		it++; cout<<" iteration "<<it<<" : "<<nch<<" pixels elagues\n";
	}*/
	return imaRes;
}

imabin imabin::squelette (imabin &ima_anchor, int iconnex, bool iaf) {
	imabin imaRes(nblig,nbcol);
	liste_pixels Lpix, L2, Lretrieve;
	elt_liste E;
	int n,i,j,i0,i2,j0,j2,icard;
	bool xv,x0,x1,x2,x3,x4,x5,x6,x7,elimine;
	imaRes=(*this);
	for (i=1; i<nblig-1; i++)                                       // on ne traite pas les bords
		for (j=1; j<nbcol-1; j++) if ((*this)(i,j)) Lpix.insere(i,j);
	int nch=1;
	while (Lpix.nb_elts()>0 && nch>0) {
		nch=0;
		for (icard=0; icard<4; icard++) { // 4 cas de voisins dans Xc : Nord, Est, Sud, Ouest
			L2=Lpix;
			while (L2.nb_elts()>0) {
				E=L2.extrait(); i=E.x; j=E.y; elimine=0;
				if (!ima_anchor(i,j)) {
					if (icard==0) xv=imaRes(i-1,j);
					else if (icard==1) xv=imaRes(i,j+1);
						else if (icard==2) xv=imaRes(i+1,j);
							else if (icard==3) xv=imaRes(i,j-1);
					if (!xv) {
						i0=i-1; j0=j-1; i2=i+1; j2=j+1;
						x0=imaRes(i,j2); x1=imaRes(i0,j2); x2=imaRes(i0,j); x3=imaRes(i0,j0);
						x4=imaRes(i,j0); x5=imaRes(i2,j0); x6=imaRes(i2,j); x7=imaRes(i2,j2);
						n=0;
						if (iconnex==8) {
							if (!x0 && (x1 || x2)) n++; if (!x2 && (x3 || x4)) n++; if (!x4 && (x5 || x6)) n++;	if (!x6 && (x7 || x0)) n++;
						}	else {
							if (x0 && (!x1 || !x2)) n++; if (x2 && (!x3 || !x4)) n++; if (x4 && (!x5 || !x6)) n++; if (x6 && (!x7 || !x0)) n++;
						}
						if (n==1) elimine=1;
					}
				}
				if (elimine) Lretrieve.insere(i,j);
			}
			while (Lretrieve.nb_elts()>0) {
				E=Lretrieve.extrait(); imaRes(E.x,E.y)=0;
				Lpix.supprime(E.x,E.y); nch++;                                // a verifier
			}
		}
	}
	return imaRes;
}

/* ---------------------------------------------------------------------*/
/* Transformee de Fourier                                               */
imadata<float> imabin::transformee_Fourier() {
/* ---------------------------------------------------------------------*/
	const int yc=nbcol/2, xc=nblig/2;
	const float xnli=(float)nblig, xnco=(float)nbcol, PI_2=2.f*(float)PI;
	int		i,j,u,v;
	double   z,xx,yy,cosz,sinz;
	imadata<float> imaFourier(nblig,nbcol,3);
	cout<<" image de la transformee de Fourier : "<<imaFourier.nlig()<<" lignes et "<<imaFourier.ncol()<<" colonnes\n";
	imaFourier.mise_a_zero ();
	imadata<double> imaF_1D(nblig,nbcol,2);
	for (j=0; j<nbcol; j++) { //cout<<" colonne "<<j<<" ";
		for (i=0; i<nblig; i++) {
			xx=0; yy=0;
			for (u=0; u<nblig; u++) {
				z=-PI_2*u*i/xnli;
				xx+=(*this)(u,j)*cos(z);
				yy+=(*this)(u,j)*sin(z);
			}
			imaF_1D(i,j,0)=xx;
			imaF_1D(i,j,1)=yy;
		}
	}
	for (i=0; i<nblig; i++) { //cout<<" ligne "<<i<<" ";
		for (j=0; j<nbcol; j++) {
			xx=0; yy=0;
			for (v=0; v<nbcol; v++) {
				z=-PI_2*v*j/xnco; cosz=cos(z); sinz=sin(z);
				xx+=imaF_1D(i,v,0)*cosz-imaF_1D(i,v,1)*sinz;
				yy+=imaF_1D(i,v,0)*sinz+imaF_1D(i,v,1)*cosz;
			}
			imaFourier(i,j,0)=(float)xx;
			imaFourier(i,j,1)=(float)yy;
		}
	}
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) {
			imaFourier(i,j,0)=imaFourier(i,j,0)/nblig/nbcol;
			imaFourier(i,j,1)=imaFourier(i,j,1)/nblig/nbcol;
		}
/*	for (i=0; i<nblig; i++) { if (i%1==0) cout<<" ligne "<<i<<"\n";
		for (j=0; j<nbcol; j++) {
			xx=0; yy=0;
			for (u=0; u<nblig; u++)
				for (v=0; v<nbcol; v++) {
					z=-PI_2*(u*i/xnli+v*j/xnco);
					xx+=(*this)(u,v)*cos(z);
					yy+=(*this)(u,v)*sin(z);
				}
			imaFourier(i,j,0)=xx;
			imaFourier(i,j,1)=yy;
		}
	}*/
	for (i=0; i<nblig/2; i++)
		for (j=0; j<nbcol/2; j++) {
			xx=imaFourier(i,j,0); yy=imaFourier(i,j,1);
			imaFourier(i,j,0)=imaFourier(i+xc,j+yc,0); imaFourier(i,j,1)=imaFourier(i+xc,j+yc,1);
			imaFourier(i+xc,j+yc,0)=(float)xx; imaFourier(i+xc,j+yc,1)=(float)yy;
			xx=imaFourier(i+xc,j,0); yy=imaFourier(i+xc,j,1);
			imaFourier(i+xc,j,0)=imaFourier(i,j+yc,0); imaFourier(i+xc,j,1)=imaFourier(i,j+yc,1);
			imaFourier(i,j+yc,0)=(float)xx; imaFourier(i,j+yc,1)=(float)yy;
		}
	for (i=0; i<imaFourier.nlig(); i++)
		for (j=0; j<imaFourier.ncol(); j++) imaFourier(i,j,2)=pow(pow(imaFourier(i,j,0),2)+pow(imaFourier(i,j,1),2),0.5f);
	imaFourier.statbasic(1);
	return imaFourier;
}

/* ---------------------------------------------------------------------*/
/* Transformee de Fourier                                               */
imadata<float> imabin::transformee_Fourier_Inv() {
/* ---------------------------------------------------------------------*/
	const int yc=nbcol/2, xc=nblig/2;
	const float xnli=(float)nblig, xnco=(float)nbcol, PI_2=2.f*(float)PI;
	int		i,j,u,v;
	double   z,xx,yy,cosz,sinz;
	imadata<float> imaFourier(*this), imaFourierInv(nblig,nbcol,3);
	cout<<" image de la transformee de Fourier inverse : "<<imaFourierInv.nlig()<<" lignes et "<<imaFourierInv.ncol()<<" colonnes\n";
	imaFourierInv.mise_a_zero ();
	for (i=0; i<nblig/2; i++)
		for (j=0; j<nbcol/2; j++) {
			xx=imaFourier(i,j); 
			imaFourier(i,j)=imaFourier(i+xc,j+yc);
			imaFourier(i+xc,j+yc)=(float)xx; 
			xx=imaFourier(i+xc,j);
			imaFourier(i+xc,j)=imaFourier(i,j+yc); 
			imaFourier(i,j+yc)=(float)xx;
/*			xx=imaFourier(i,j,0); yy=imaFourier(i,j,1);
			imaFourier(i,j,0)=imaFourier(i+xc,j+yc,0); imaFourier(i,j,1)=imaFourier(i+xc,j+yc,1);
			imaFourier(i+xc,j+yc,0)=xx; imaFourier(i+xc,j+yc,1)=yy;
			xx=imaFourier(i+xc,j,0); yy=imaFourier(i+xc,j,1);
			imaFourier(i+xc,j,0)=imaFourier(i,j+yc,0); imaFourier(i+xc,j,1)=imaFourier(i,j+yc,1);
			imaFourier(i,j+yc,0)=xx; imaFourier(i,j+yc,1)=yy;*/
		}
	imadata<double> imaF_1D(nblig,nbcol,2);
	for (j=0; j<nbcol; j++) { //cout<<" colonne "<<j<<" ";
		for (i=0; i<nblig; i++) {
			xx=0; yy=0;
			for (u=0; u<nblig; u++) {
				z=+PI_2*u*i/xnli;
				xx+=imaFourier(u,j)*cos(z);
				yy+=imaFourier(u,j)*sin(z);
			}
			imaF_1D(i,j,0)=xx;
			imaF_1D(i,j,1)=yy;
		}
	}
	for (i=0; i<nblig; i++) { //cout<<" ligne "<<i<<" ";
		for (j=0; j<nbcol; j++) {
			xx=0; yy=0;
			for (v=0; v<nbcol; v++) {
				z=+PI_2*v*j/xnco; cosz=cos(z); sinz=sin(z);
				xx+=imaF_1D(i,v,0)*cosz-imaF_1D(i,v,1)*sinz;
				yy+=imaF_1D(i,v,0)*sinz+imaF_1D(i,v,1)*cosz;
			}
			imaFourierInv(i,j,0)=(float)xx;
			imaFourierInv(i,j,1)=(float)yy;
		}
	}
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++) {
			imaFourierInv(i,j,0)=imaFourierInv(i,j,0)/nblig/nbcol;
			imaFourierInv(i,j,1)=imaFourierInv(i,j,1)/nblig/nbcol;
		}
	for (i=0; i<imaFourierInv.nlig(); i++)
		for (j=0; j<imaFourierInv.ncol(); j++) imaFourierInv(i,j,2)=pow(pow(imaFourierInv(i,j,0),2)+pow(imaFourierInv(i,j,1),2),0.5f);
	imaFourierInv.statbasic(1);
	return imaFourierInv;
}

/* ---------------------------------------------------------------------*/
/* Transformee de Hough : precision d'angle de 1 degre                  */
imadata<int> imabin::transformee_Hough_1pt(bool iaffich)
/* ---------------------------------------------------------------------*/ {
	const int yc=nbcol/2, xc=nblig/2, N_Rho=(int)sqrt((float)nblig*nblig+nbcol*nbcol);
	int		i,j,iRho,iTheta;
	const float pas_Theta_deg=1.f;
	float   Theta_deg, Theta_rad, Rho/*, r, phi*/;
	imadata<int> imaHough(N_Rho+1,360);
	cout<<" image de la transformee de Hough : "<<imaHough.nlig()<<" lignes et "<<imaHough.ncol()<<" colonnes\n";
	imaHough.mise_a_zero();
	imabin Ima_deja(N_Rho+1,360);
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++)
			if ((*this)(i,j)) {
				Ima_deja.mise_a_zero(); // Ima_deja sert à ce qu'un pixel de l'image d'origine ne vote pas plusieurs fois pour le même pixel de la transforméee de Hough
//				r=pow((float)((j-yc)*(j-yc)+(i-xc)*(i-xc)),0.5f); phi=atan2((float)(i-xc),(float)(j-yc));
				for (Theta_deg=-180.f; Theta_deg<180.f; Theta_deg+=pas_Theta_deg) {
					Theta_rad=Theta_deg*(float)PI/180.f;
					Rho=(j-yc)*sin(Theta_rad)+(i-xc)*cos(Theta_rad);
					iRho=(int)(Rho+N_Rho/2); iTheta=(int)Theta_deg+180;
					Ima_deja(iRho,iTheta)=1;
				}
				for (iTheta=0; iTheta<360; iTheta++) 
					for (iRho=0; iRho<=N_Rho; iRho++)
						if (Ima_deja(iRho,iTheta)) imaHough(iRho,iTheta)++;
			}
	if (iaffich) {
		cout<<" transformee de Hough : axe des abscisses = angle des droites, de -180 a +180 degres";
		cout<<", axe des ordonnees = distance de la droite au centre du domaine image, 0 (droites radiales au milieu de l'axe";
		cout<<", valeurs croisantes de rho (de -rhomax a +rhomax) avec l'indice de ligne";
		cout<<", valeurs positives/negatives fct de |rho| et theta: droite intersecte la colonne du milieu du domaine de Hough";
		cout<<" en rho/cos(theta), et la ligne du milieu du domaine de Hough en rho/sin(theta)\n";
		imaHough.statbasic(1); 
	} else 
		imaHough.statbasic(0); 
	return imaHough;
}

/* ---------------------------------------------------------------------*/
/* Transformee de Hough : precision d'angle de 1 degre, de rho : 1 pixel*/
imadata<int> imabin::transformee_Hough()
/* ---------------------------------------------------------------------*/ {
	const int yc=nbcol/2, xc=nblig/2, N_Rho=(int)sqrt((float)nblig*nblig+nbcol*nbcol);
	int		i,j,iRho,iRho2,iTheta;
//	const float pas_Theta_deg=0.25f, pas_Rho=0.025f;
	const float pas_Theta_deg=1.f, pas_Rho=0.1f;
	float   Theta_deg, Theta_rad, Rho, r, phi;
	imadata<int> imaHough(N_Rho+1,360);
	cout<<" image de la transformee de Hough : "<<imaHough.nlig()<<" lignes et "<<imaHough.ncol()<<" colonnes\n";
	imaHough.mise_a_zero();
	imabin Ima_deja(N_Rho+1,360);
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++)
			if ((*this)(i,j)) {
				Ima_deja.mise_a_zero(); // Ima_deja sert à ce qu'un pixel de l'image d'origine ne vote pas plusieurs fois pour le même pixel de la transforméee de Hough
				r=pow((float)((j-yc)*(j-yc)+(i-xc)*(i-xc)),0.5f); phi=atan2((float)(i-xc),(float)(j-yc));
				for (Theta_deg=-180.f; Theta_deg<180.f; Theta_deg+=pas_Theta_deg) {
					Theta_rad=Theta_deg*(float)PI/180.f;
					Rho=(j-yc)*sin(Theta_rad)+(i-xc)*cos(Theta_rad);
					iRho=(int)(Rho+N_Rho/2); iTheta=(int)Theta_deg+180;
					if (iRho>=0 && iRho<=N_Rho)
						if (!Ima_deja(iRho,iTheta)) {imaHough(iRho,iTheta)++; Ima_deja(iRho,iTheta)=1;}
				}
				if (r>0) {
					for (Rho=0.f; Rho<=mini(r,N_Rho/2.f); Rho+=pas_Rho) {
						iRho2=(int)(Rho+N_Rho/2); iRho=(int)(N_Rho/2-Rho);
						Theta_rad=modulo2PI((float)PI/2-(-acos(Rho/r)+phi)); 
						iTheta=(int)(Theta_rad/(float)PI*180.f); if (iTheta==360) iTheta=0;
						if (!Ima_deja(iRho,iTheta)) {imaHough(iRho,iTheta)++; Ima_deja(iRho,iTheta)=1;}
						iTheta=(int)(modulo2PI(Theta_rad+PI)/PI*180.f); if (iTheta==360) iTheta=0;
						if (!Ima_deja(iRho2,iTheta)) {imaHough(iRho2,iTheta)++; Ima_deja(iRho2,iTheta)=1;}
						Theta_rad=modulo2PI((float)PI/2-(acos(Rho/r)+phi-2*(float)PI)); 
						iTheta=(int)(Theta_rad/(float)PI*180.f); if (iTheta==360) iTheta=0;
						if (!Ima_deja(iRho,iTheta)) {imaHough(iRho,iTheta)++; Ima_deja(iRho,iTheta)=1;}
						iTheta=(int)(modulo2PI(Theta_rad+PI)/PI*180.f); if (iTheta==360) iTheta=0;
						if (!Ima_deja(iRho2,iTheta)) {imaHough(iRho2,iTheta)++; Ima_deja(iRho2,iTheta)=1;}
					}
				}
			}
	cout<<" transformee de Hough : axe des abscisses = angle des droites, de -180 a +180 degres";
	cout<<", axe des ordonnees = distance de la droite au centre du domaine image, 0 (droites radiales au milieu de l'axe";
	cout<<", valeurs croisantes de rho (de -rhomax a +rhomax) avec l'indice de ligne";
	cout<<", valeurs positives/negatives fct de |rho| et theta: droite intersecte la colonne du milieu du domaine de Hough";
	cout<<" en rho/cos(theta), et la ligne du milieu du domaine de Hough en rho/sin(theta)\n";
	imaHough.statbasic(1); 
	return imaHough;
}

/* ---------------------------------------------------------------------*/
/* Transformee de Hough : precision d'angle de 1 degre, de rho : 1 pixel*/
void imabin::reconst_hough_transform (imadata<float> &imaHough, int nblig_ima, int nbcol_ima)
/* ---------------------------------------------------------------------*/ {
	const double /*eps=1.e-6*/ eps=0.1;
	const int N_Rho=imaHough.nlig(), N_Theta=imaHough.ncol(), n_dx=10;
	cout<<" image de la transformee de Hough de "<<N_Rho<<" lignes sur "<<N_Theta<<" colonnes\n";
	int		xc, yc, iTheta, iRho, i,j,k;
	float x, y, fTheta, fRho;
	imaHough.statbasic(1);
//	float s=imaHough.moyI()+15*pow(imaHough.varI(),0.5);
//	float s=(float)(imaHough.moyI()+10*pow(imaHough.varI(),0.5));
//	float s=imaHough.moyI()+5*pow(imaHough.varI(),0.5);
	cout<<" seuil a 10% = "<<imaHough.seuil_percentile(1.f-0.1f)<<" seuil a 5% = "<<imaHough.seuil_percentile(1.f-0.05f);
	cout<<" seuil a 1% = "<<imaHough.seuil_percentile(1.f-0.01f)<<" seuil a 0.5% = "<<imaHough.seuil_percentile(1.f-0.005f);
	cout<<" seuil a 0.1% = "<<imaHough.seuil_percentile(1.f-0.001f)<<" seuil a 0.05% = "<<imaHough.seuil_percentile(1.f-0.0005f);
	float s=imaHough.seuil_percentile(1.f-0.01f);
	cout<<" valeurs depassant le seuil "<<s<<"\n";
	imabin ima_max(imaHough,s); ima_max.imaunsignedchar().sauve_ImaPGM("./max_imaTHough.pgm");
//	ima_max=ima_max.erode_ultime(); ima_max.imaunsignedchar().sauve_ImaPGM("./max_imaTHough_EU.pgm");
	eltstruct S3(5,5); ima_max=ima_max.dilate(S3);
	imabin ima_marq(N_Rho,N_Theta); ima_marq.mise_a_zero();
	for (i=0; i<N_Rho; i++) ima_marq(i,0)=ima_max(i,N_Theta-1); // à cause du modulo sur les angles
	ima_max=ima_max+ima_max.reconstruction_geodesique(ima_marq,8); ima_max.imaunsignedchar().sauve_ImaPGM("./max_imaTHough_D5.pgm");
	int ncc, n;
	imadata<int> ima_max_cc=ima_max.composantes_connexes(ncc,8,0);
	cout<<ncc<<" composantes connexes\n";
	float val_max=-999;
	for (k=1; k<=ncc; k++) {
		n=xc=yc=0; val_max=-999;
		for (iRho=0; iRho<N_Rho; iRho++)
			for (iTheta=0; iTheta<N_Theta; iTheta++)
				if (ima_max_cc(iRho,iTheta)==k) {
//					xc+=iRho; yc+=iTheta; n++; 
					if (imaHough(iRho,iTheta)>val_max) {val_max=imaHough(iRho,iTheta); xc=iRho; yc=iTheta;}
					ima_max(iRho,iTheta)=0;
				}
//		if (n>0) ima_max(xc/=n,yc/=n)=1;
		ima_max(xc,yc)=1;
	} ima_max.imaunsignedchar().sauve_ImaPGM("./max_imaTHough_1pixperCC.pgm");
	n=0;
/*	for (iRho=0; iRho<N_Rho; iRho++)
		for (iTheta=0; iTheta<N_Theta; iTheta++)
			if (imaHough(iRho,iTheta)>=s) { //cout<<iRho<<" "<<iTheta<<" ";
				cout<<" presence d'un max en Rho="<<iRho-N_Rho/2<<" vs centre image, et theta="<<iTheta-N_Theta/2<<" deg. de valeur "<<imaHough(iRho,iTheta)<<"\n";
				n++;
			}*/
	for (iRho=0; iRho<N_Rho; iRho++)
		for (iTheta=0; iTheta<N_Theta; iTheta++)
			if (ima_max(iRho,iTheta)) {
				cout<<" presence d'un max en Rho="<<iRho-N_Rho/2<<" vs centre image, et theta="<<iTheta-N_Theta/2<<" deg. de valeur "<<imaHough(iRho,iTheta)<<"\n";
				n++;
			}
	val_pos *T=new val_pos[n];
	n=0;
	for (iRho=0; iRho<N_Rho; iRho++)
		for (iTheta=0; iTheta<N_Theta; iTheta++)
//		for (iTheta=170; iTheta<190; iTheta++)
//			if (imaHough(iRho,iTheta)>s) { 
			if (ima_max(iRho,iTheta)) { 
				T[n].pos=iRho*N_Theta+iTheta; T[n].val=imaHough(iRho,iTheta); n++;
			}
	tri_rapide(T,n);
//	imabin imaReconsH(nblig_ima,nbcol_ima);
	nblig=nblig_ima; nbcol=nbcol_ima; yc=nbcol/2; xc=nblig/2; 
	(*this)=imabin(nblig,nbcol); mise_a_zero(); cout<<" en sortie image reconstruite de "<<xc<<" lig. et "<<yc<<" col.\n";
	for (k=0; k<n; k++) {
		cout<<" max en Rho="<<(int)(T[k].pos)/N_Theta-N_Rho/2<<" vs centre im., et theta="<<(int)(T[k].pos)%N_Theta-N_Theta/2<<" deg. de valeur "<<T[k].val<<"\n";
//		cout<<n<<" "<<T[k].pos<<" presence d'un max en Rho="<<(int)(T[k].pos)/360-N_Rho/2<<" vs centre image, et theta="<<(int)(T[k].pos)%360-180<<" deg. de valeur "<<T[k].val<<"\n";
		fTheta=((int)(T[k].pos)%N_Theta-N_Theta/2)*(float)PI/180.f;
		fRho=(float)(T[k].pos)/N_Theta-N_Rho/2;
//		if (fabs(sin(fTheta))>eps && fabs(cos(fTheta))>eps) {
		cout<<" reconstruction dans l'image\n";
		if (fabs(sin(fTheta))>eps) {
			for (i=0; i<nblig; i++) {
				for (j=0; j<n_dx; j++) {
					x=i+(float)j/n_dx;
					y=(fRho-(x-xc)*cos(fTheta))/sin(fTheta)+yc;
//					if (y>0 && y<nbcol) imaReconsH(mini(around(x),nblig-1),mini(maxi(0,around(y)),nbcol-1))=1;
					if (y>0 && y<nbcol) (*this)(mini(around(x),nblig-1),mini(maxi(0,around(y)),nbcol-1))=1;
				}
			}
		} else {
			if (fabs(sin(fTheta))<=eps) {
				for (j=0; j<nbcol; j++) {
					for (i=0; i<n_dx; i++) {
						y=j+(float)i/n_dx;
						x=(fRho-(y-yc)*sin(fTheta))/cos(fTheta)+xc;
						if (x>0 && x<nblig) (*this)(mini(around(x),nblig-1),mini(maxi(0,around(y)),nbcol-1))=1;
					}
				}
			}
		}
	}
	if (T!=NULL) delete[] T;
//	imaReconsH.imaunsignedchar().sauve_ImaPGM("reconstructionHough.pgm");
//	return imaReconsH;
}

void imabin::reconst_hough_transform (imadata<float> &imaHough, int nblig_ima, int nbcol_ima, imabin ima_max)
/* ---------------------------------------------------------------------*/ {
	const double /*eps=1.e-6*/ eps=0.1;
	const int N_Rho=imaHough.nlig(), N_Theta=imaHough.ncol(), n_dx=10;
	cout<<" image de la transformee de Hough de "<<N_Rho<<" lignes sur "<<N_Theta<<" colonnes\n";
	int		xc, yc, iTheta, iRho, i,j,k;
	float x, y, fTheta, fRho;
/*	imaHough.statbasic(1);
	cout<<" seuil a 10% = "<<imaHough.seuil_percentile(1.-0.1)<<" seuil a 5% = "<<imaHough.seuil_percentile(1.-0.05);
	cout<<" seuil a 1% = "<<imaHough.seuil_percentile(1.-0.01)<<" seuil a 0.5% = "<<imaHough.seuil_percentile(1.-0.005);
	cout<<" seuil a 0.1% = "<<imaHough.seuil_percentile(1.-0.001)<<" seuil a 0.05% = "<<imaHough.seuil_percentile(1.-0.0005);
	float s=imaHough.seuil_percentile(1.-0.01);
	cout<<" valeurs depassant le seuil "<<s<<"\n";
	imabin ima_max(imaHough,s); ima_max.imaunsignedchar().sauve_ImaPGM("./max_imaTHough.pgm");*/
	eltstruct S3(5,5); ima_max=ima_max.dilate(S3);
	imabin ima_marq(N_Rho,N_Theta); ima_marq.mise_a_zero();
	for (i=0; i<N_Rho; i++) ima_marq(i,0)=ima_max(i,N_Theta-1); // à cause du modulo sur les angles
	ima_max=ima_max+ima_max.reconstruction_geodesique(ima_marq,8); ima_max.imaunsignedchar().sauve_ImaPGM("./max_imaTHough_D5.pgm");
	int ncc, n;
	imadata<int> ima_max_cc=ima_max.composantes_connexes(ncc,8,0);
	cout<<ncc<<" composantes connexes\n";
	float val_max=-999;
	for (k=1; k<=ncc; k++) {
		n=xc=yc=0; val_max=-999;
		for (iRho=0; iRho<N_Rho; iRho++)
			for (iTheta=0; iTheta<N_Theta; iTheta++)
				if (ima_max_cc(iRho,iTheta)==k) {
					if (imaHough(iRho,iTheta)>val_max) {val_max=imaHough(iRho,iTheta); xc=iRho; yc=iTheta;}
					ima_max(iRho,iTheta)=0;
				}
		ima_max(xc,yc)=1;
	} ima_max.imaunsignedchar().sauve_ImaPGM("./max_imaTHough_1pixperCC.pgm");
	n=(int)ima_max.norm();
	val_pos *T=new val_pos[n];
	n=0;
	for (iRho=0; iRho<N_Rho; iRho++)
		for (iTheta=0; iTheta<N_Theta; iTheta++)
			if (ima_max(iRho,iTheta)) { 
				T[n].pos=iRho*N_Theta+iTheta; T[n].val=imaHough(iRho,iTheta); n++;
			}
	tri_rapide(T,n);
	nblig=nblig_ima; nbcol=nbcol_ima; yc=nbcol/2; xc=nblig/2; 
	(*this)=imabin(nblig,nbcol); mise_a_zero(); cout<<" en sortie image reconstruite de "<<xc<<" lig. et "<<yc<<" col.\n";
	for (k=0; k<n; k++) {
		cout<<" max en Rho="<<(int)(T[k].pos)/N_Theta-N_Rho/2<<" vs centre im., et theta="<<(int)(T[k].pos)%N_Theta-N_Theta/2<<" deg. de valeur "<<T[k].val<<"\n";
		fTheta=((int)(T[k].pos)%N_Theta-N_Theta/2)*(float)PI/180.f;
		fRho=(float)(T[k].pos)/N_Theta-N_Rho/2;
//		if (fabs(sin(fTheta))>eps && fabs(cos(fTheta))>eps) {
		if (fabs(sin(fTheta))>eps) {
			cout<<" reconstruction dans l'image\n";
			for (i=0; i<nblig; i++) {
				for (j=0; j<n_dx; j++) {
					x=i+(float)j/n_dx;
					y=(fRho-(x-xc)*cos(fTheta))/sin(fTheta)+yc;
//					if (y>0 && y<nbcol) imaReconsH(mini(around(x),nblig-1),mini(maxi(0,around(y)),nbcol-1))=1;
					if (y>0 && y<nbcol) (*this)(mini(around(x),nblig-1),mini(maxi(0,around(y)),nbcol-1))=1;
				}
			}
		} else {
			if (fabs(sin(fTheta))<=eps) {
				for (j=0; j<nbcol; j++) {
					for (i=0; i<n_dx; i++) {
						y=j+(float)i/n_dx;
						x=(fRho-(y-yc)*sin(fTheta))/cos(fTheta)+xc;
						if (x>0 && x<nblig) (*this)(mini(around(x),nblig-1),mini(maxi(0,around(y)),nbcol-1))=1;
					}
				}
			}
		}
	}
	if (T!=NULL) delete[] T;
}

/* ---------------------------------------------------------------------*/
/* Transformee de Hough : precision d'angle de 1 degre, de rho : 1 pixel*/
imadata<int> imabin::norme_transformee_Hough(int ligTH0, int colTH0, int d_ligTH, int d_colTH)
/* ---------------------------------------------------------------------*/ {
	const int yc=nbcol/2, xc=nblig/2, N_Rho=(int)sqrt((float)nblig*nblig+nbcol*nbcol), nvois=13;
	const int nbligTH=N_Rho+1, nbcolTH=360, nbcanTH=nvois;
	int		i,j,k,x,y,iphi,iRho,iRho2;
	const float pas_deg=1.f/*0.1f*/, pas_Rho=0.1f;
	float   Theta_rad,Rho,Rho2,r,r2,phi,phi2,Rho2min,Rho2max;

	imadata<int> imaHough_Norm(nbligTH,nbcolTH,nbcanTH);
	cout<<" image-norme de la transformee de Hough : "<<nbligTH<<" lignes, "<<nbcolTH<<" colonnes et "<<nbcanTH<<" canaux\n";
	imaHough_Norm.mise_a_zero();
	imabin Ima_deja(nblig,nbcol),ImaTH_deja(nbligTH,nbcolTH); // diminuer les dimensions de ImaTH_deja au voisinage
	for (i=maxi(ligTH0,0); i<=mini(ligTH0+d_ligTH,nbligTH-1); i++)
		for (j=maxi(colTH0,0); j<=mini(colTH0+d_colTH,nbcolTH-1); j++) {
			Rho=(float)i-N_Rho/2; phi=modulo2PI((float)(j-180)*(float)PI/180.f); 
			if (i%10==0 && j%10==0) cout<<i<<" "<<j<<" Rho = "<<Rho<<", phi = "<<phi/PI*180<<"\n";
			Ima_deja.mise_a_zero();
			short int sign_Rho=(Rho>=0?1:-1);
			for (r=Rho+sign_Rho*pas_Rho; fabs(r)<N_Rho; r+=sign_Rho*pas_Rho) {r2=r;
				if (r2!=0) {
					for (k=0; k<2; k++) { 
						if (k==0) Theta_rad=modulo2PI(-phi+asin(Rho/r2)); 
						if (k==1) Theta_rad=modulo2PI((float)PI-phi-asin(Rho/r2)); 
						x=(int)(r2*sin(Theta_rad))+xc; y=(int)(r2*cos(Theta_rad))+yc;
						if (x>=0 && x<nblig && y>=0 && y<nbcol && !Ima_deja(x,y)) {
							ImaTH_deja.mise_a_zero();
							r=sign_Rho*pow((float)((y-yc)*(y-yc)+(x-xc)*(x-xc)),0.5f); Theta_rad=modulo2PI(atan2((float)(x-xc),(float)(y-yc)));
							for (phi2=phi-2*(float)PI/180.f; phi2<=phi+2*(float)PI/180.f; phi2+=pas_deg*(float)PI/180.f) {
								Rho2=(y-yc)*sin(phi2)+(x-xc)*cos(phi2);
								iRho=(int)(Rho2+N_Rho/2); iphi=(int)(phi2/PI*180)+180; if (iphi>=360) iphi-=360;
								if (iRho>=0 && iRho<=N_Rho && !ImaTH_deja(iRho,iphi)) {
//									imaHough_Norm(iRho,iphi,0)++;
									if (iRho==i && iphi==j) imaHough_Norm(iRho,iphi,0)++; 
									else if (iRho==i && iphi==j+1) imaHough_Norm(iRho,iphi,1)++;
										else if (iRho==i+1 && iphi==j) imaHough_Norm(iRho,iphi,2)++;
											else if (iRho==i && iphi==j-1) imaHough_Norm(iRho,iphi,3)++;
												else if (iRho==i-1 && iphi==j) imaHough_Norm(iRho,iphi,4)++;
													else if (iRho==i+1 && iphi==j+1) imaHough_Norm(iRho,iphi,5)++;
														else if (iRho==i+1 && iphi==j-1) imaHough_Norm(iRho,iphi,6)++;
															else if (iRho==i-1 && iphi==j-1) imaHough_Norm(iRho,iphi,7)++;
																else if (iRho==i-1 && iphi==j+1) imaHough_Norm(iRho,iphi,8)++;
																	else if (iRho==i && iphi==j+2) imaHough_Norm(iRho,iphi,9)++;
																		else if (iRho==i+2 && iphi==j) imaHough_Norm(iRho,iphi,10)++;
																			else if (iRho==i && iphi==j-2) imaHough_Norm(iRho,iphi,11)++;
																				else if (iRho==i-2 && iphi==j) imaHough_Norm(iRho,iphi,12)++;
									ImaTH_deja(iRho,iphi)=1;
								}
							}
							if (r!=0) {
								if (Rho<0) {Rho2min=maxi(Rho-2.f,r); Rho2max=mini(Rho+2.f,0.f);}
								else {Rho2min=maxi(Rho-2.f,0.f); Rho2max=mini(Rho+2.f,r);} 
								for (Rho2=Rho2min; Rho2<=Rho2max; Rho2+=pas_Rho) {
									iRho2=(int)(N_Rho/2+sign_Rho*Rho2); iRho=(int)(N_Rho/2-sign_Rho*Rho2); 
									phi2=modulo2PI((float)PI/2-(-acos(Rho2/r)+Theta_rad)); iphi=(int)(phi2/PI*180.);
									if (!ImaTH_deja(iRho,iphi)) {
//										imaHough_Norm(iRho,iphi,0)++;
										if (iRho==i && iphi==j) imaHough_Norm(iRho,iphi,0)++; 
										else if (iRho==i && iphi==j+1) imaHough_Norm(iRho,iphi,1)++;
											else if (iRho==i+1 && iphi==j) imaHough_Norm(iRho,iphi,2)++;
												else if (iRho==i && iphi==j-1) imaHough_Norm(iRho,iphi,3)++;
													else if (iRho==i-1 && iphi==j) imaHough_Norm(iRho,iphi,4)++;
														else if (iRho==i+1 && iphi==j+1) imaHough_Norm(iRho,iphi,5)++;
															else if (iRho==i+1 && iphi==j-1) imaHough_Norm(iRho,iphi,6)++;
																else if (iRho==i-1 && iphi==j-1) imaHough_Norm(iRho,iphi,7)++;
																	else if (iRho==i-1 && iphi==j+1) imaHough_Norm(iRho,iphi,8)++;
																		else if (iRho==i && iphi==j+2) imaHough_Norm(iRho,iphi,9)++;
																			else if (iRho==i+2 && iphi==j) imaHough_Norm(iRho,iphi,10)++;
																				else if (iRho==i && iphi==j-2) imaHough_Norm(iRho,iphi,11)++;
																					else if (iRho==i-2 && iphi==j) imaHough_Norm(iRho,iphi,12)++;
										ImaTH_deja(iRho,iphi)=1;}
									iphi=(int)(modulo2PI(phi2+PI)/PI*180.f);
									if (!ImaTH_deja(iRho2,iphi)) {
//										imaHough_Norm(iRho2,iphi,0)++;
										if (iRho2==i && iphi==j) imaHough_Norm(iRho2,iphi,0)++; 
										else if (iRho2==i && iphi==j+1) imaHough_Norm(iRho2,iphi,1)++;
											else if (iRho2==i+1 && iphi==j) imaHough_Norm(iRho2,iphi,2)++;
												else if (iRho2==i && iphi==j-1) imaHough_Norm(iRho2,iphi,3)++;
													else if (iRho2==i-1 && iphi==j) imaHough_Norm(iRho2,iphi,4)++;
														else if (iRho2==i+1 && iphi==j+1) imaHough_Norm(iRho2,iphi,5)++;
															else if (iRho2==i+1 && iphi==j-1) imaHough_Norm(iRho2,iphi,6)++;
																else if (iRho2==i-1 && iphi==j-1) imaHough_Norm(iRho2,iphi,7)++;
																	else if (iRho2==i-1 && iphi==j+1) imaHough_Norm(iRho2,iphi,8)++;
																		else if (iRho2==i && iphi==j+2) imaHough_Norm(iRho2,iphi,9)++;
																			else if (iRho2==i+2 && iphi==j) imaHough_Norm(iRho2,iphi,10)++;
																				else if (iRho2==i && iphi==j-2) imaHough_Norm(iRho2,iphi,11)++;
																					else if (iRho2==i-2 && iphi==j) imaHough_Norm(iRho2,iphi,12)++;
										ImaTH_deja(iRho2,iphi)=1;}
									phi2=(float)modulo2PI(PI/2-(acos(Rho2/r)+Theta_rad-2*PI)); iphi=(int)(phi2/PI*180.);
									if (!ImaTH_deja(iRho,iphi)) {
//										imaHough_Norm(iRho,iphi,0)++;
										if (iRho==i && iphi==j) imaHough_Norm(iRho,iphi,0)++; 
										else if (iRho==i && iphi==j+1) imaHough_Norm(iRho,iphi,1)++;
											else if (iRho==i+1 && iphi==j) imaHough_Norm(iRho,iphi,2)++;
												else if (iRho==i && iphi==j-1) imaHough_Norm(iRho,iphi,3)++;
													else if (iRho==i-1 && iphi==j) imaHough_Norm(iRho,iphi,4)++;
														else if (iRho==i+1 && iphi==j+1) imaHough_Norm(iRho,iphi,5)++;
															else if (iRho==i+1 && iphi==j-1) imaHough_Norm(iRho,iphi,6)++;
																else if (iRho==i-1 && iphi==j-1) imaHough_Norm(iRho,iphi,7)++;
																	else if (iRho==i-1 && iphi==j+1) imaHough_Norm(iRho,iphi,8)++;
																		else if (iRho==i && iphi==j+2) imaHough_Norm(iRho,iphi,9)++;
																			else if (iRho==i+2 && iphi==j) imaHough_Norm(iRho,iphi,10)++;
																				else if (iRho==i && iphi==j-2) imaHough_Norm(iRho,iphi,11)++;
																					else if (iRho==i-2 && iphi==j) imaHough_Norm(iRho,iphi,12)++;
										ImaTH_deja(iRho,iphi)=1;}
									iphi=(int)(modulo2PI(phi2+PI)/PI*180.f);
									if (!ImaTH_deja(iRho2,iphi)) {
//										imaHough_Norm(iRho2,iphi,0)++;
										if (iRho2==i && iphi==j) imaHough_Norm(iRho2,iphi,0)++; 
										else if (iRho2==i && iphi==j+1) imaHough_Norm(iRho2,iphi,1)++;
											else if (iRho2==i+1 && iphi==j) imaHough_Norm(iRho2,iphi,2)++;
												else if (iRho2==i && iphi==j-1) imaHough_Norm(iRho2,iphi,3)++;
													else if (iRho2==i-1 && iphi==j) imaHough_Norm(iRho2,iphi,4)++;
														else if (iRho2==i+1 && iphi==j+1) imaHough_Norm(iRho2,iphi,5)++;
															else if (iRho2==i+1 && iphi==j-1) imaHough_Norm(iRho2,iphi,6)++;
																else if (iRho2==i-1 && iphi==j-1) imaHough_Norm(iRho2,iphi,7)++;
																	else if (iRho2==i-1 && iphi==j+1) imaHough_Norm(iRho2,iphi,8)++;
																		else if (iRho2==i && iphi==j+2) imaHough_Norm(iRho2,iphi,9)++;
																			else if (iRho2==i+2 && iphi==j) imaHough_Norm(iRho2,iphi,10)++;
																				else if (iRho2==i && iphi==j-2) imaHough_Norm(iRho2,iphi,11)++;
																					else if (iRho2==i-2 && iphi==j) imaHough_Norm(iRho2,iphi,12)++;
										ImaTH_deja(iRho2,iphi)=1;}
								}
							}
							Ima_deja(x,y)=1;
						}
					}
				}
				r=r2;
			}
//			Ima_deja.imaunsignedchar(1).sauve_ImaPGM("testlig_Hough_Norm.pgm");
//			(Ima_deja.transformee_Hough()).imaunsignedchar().sauve_ImaPGM("testTH_Hough_Norm.pgm");
//			imadata<float>(imaHough_Norm,0).imaunsignedchar().sauve_ImaPGM("testTH2_Hough_Norm.pgm");
//}
		}
//	cout<<" transformee de Hough : axe des abscisses = angle des droites, de -180 a +180 degres";
//	cout<<", axe des ordonnees = distance de la droite au centre du domaine image, 0 (droites radiales au milieu de l'axe";
//	cout<<", valeurs croisantes de rho (de -rhomax a +rhomax) avec l'indice de ligne";
//	cout<<", valeurs positives/negatives fct de |rho| et theta: droite intersecte la colonne du milieu du domaine de Hough";
//	cout<<" en rho/cos(theta), et la ligne du milieu du domaine de Hough en rho/sin(theta)\n";
	imaHough_Norm.sauve_ImaBSQ("./norme13_lignes_Hough.bsq");
	imaHough_Norm.statbasic(1); 
	return imaHough_Norm;
}

/* ---------------------------------------------------------------------*/
/* Transformee de Hough : precision d'angle de 1 degre, de rho : 1 pixel*/
imadata<int> imabin::norme_transformee_Hough(const imabin &masq)
/* ---------------------------------------------------------------------*/ {
	const int yc=nbcol/2, xc=nblig/2, N_Rho=(int)sqrt((float)nblig*nblig+nbcol*nbcol), nvois=13;
	const int nbligTH=N_Rho+1, nbcolTH=360, nbcanTH=nvois;
	int		i,j,k,n,x,y,iphi,iRho,iRho2;
	const float pas_deg=1.f/*0.1f*/, pas_Rho=0.1f;
	float   Theta_rad,Rho,Rho2,r,r2,phi,phi2,Rho2min,Rho2max;

	imadata<int> imaHough_Norm(nbligTH,nbcolTH,nbcanTH);
	cout<<" image-norme de la transformee de Hough : "<<nbligTH<<" lignes, "<<nbcolTH<<" colonnes et "<<nbcanTH<<" canaux\n";
	cout<<" calculee seulement sur "<<masq.norm()<<" points (donnes par un masque binaire)\n";
	imaHough_Norm.mise_a_zero();
	imabin Ima_deja(nblig,nbcol),ImaTH_deja(nbligTH,nbcolTH); // diminuer les dimensions de ImaTH_deja au voisinage
	for (i=0; i<nbligTH; i++)
		for (j=0; j<nbcolTH; j++) {
			if (masq(i,j)==1) {
				Rho=(float)i-N_Rho/2; phi=(float)modulo2PI((j-180)*PI/180.); 
				if (i%10==0 && j%10==0) cout<<i<<" "<<j<<" Rho = "<<Rho<<", phi = "<<phi/PI*180<<"\n";
				Ima_deja.mise_a_zero();
				short int sign_Rho=(Rho>=0?1:-1);
				for	(r=Rho+sign_Rho*pas_Rho; fabs(r)<N_Rho; r+=sign_Rho*pas_Rho) {r2=r;
					if (r2!=0) {
						for (k=0; k<2; k++) { 
							if (k==0) Theta_rad=modulo2PI(-phi+asin(Rho/r2)); 
							if (k==1) Theta_rad=modulo2PI((float)PI-phi-asin(Rho/r2)); 
							x=(int)(r2*sin(Theta_rad))+xc; y=(int)(r2*cos(Theta_rad))+yc;
							if (x>=0 && x<nblig && y>=0 && y<nbcol && !Ima_deja(x,y)) {
								ImaTH_deja.mise_a_zero();
								r=sign_Rho*pow((float)((y-yc)*(y-yc)+(x-xc)*(x-xc)),0.5f); Theta_rad=modulo2PI(atan2((float)(x-xc),(float)(y-yc)));
								for (phi2=phi-2*(float)PI/180.f; phi2<=phi+2*(float)PI/180.f; phi2+=pas_deg*(float)PI/180.f) {
									Rho2=(y-yc)*sin(phi2)+(x-xc)*cos(phi2);
									iRho=(int)(Rho2+N_Rho/2); iphi=(int)(phi2/PI*180)+180; if (iphi>=360) iphi-=360;
									if (iRho>=0 && iRho<=N_Rho && !ImaTH_deja(iRho,iphi)) {
//										imaHough_Norm(iRho,iphi,0)++;
										if (iRho==i && iphi==j) n=0; 
										else if (iRho==i && iphi==j+1) n=1;
											else if (iRho==i+1 && iphi==j) n=2;
												else if (iRho==i && iphi==j-1) n=3;
													else if (iRho==i-1 && iphi==j) n=4;
														else if (iRho==i+1 && iphi==j+1) n=5;
															else if (iRho==i+1 && iphi==j-1) n=6;
																else if (iRho==i-1 && iphi==j-1) n=7;
																	else if (iRho==i-1 && iphi==j+1) n=8;
																		else if (iRho==i && iphi==j+2) n=9;
																			else if (iRho==i+2 && iphi==j) n=10;
																				else if (iRho==i && iphi==j-2) n=11;
																					else if (iRho==i-2 && iphi==j) n=12;
																						else n=-1;
										if (n>=0) imaHough_Norm(i,j,n)++; 
										if (i==30 && j==31) {cout<<" point ("<<i<<","<<j<<") -> en ("<<iRho<<","<<iphi<<") -> n = "<<n<<"\n"; char aa; cin>>aa;}
										if (n>=0 && iRho==30 && iphi==31) {cout<<" point ("<<iRho<<","<<iphi<<") renforce par ("<<x<<","<<y<<") -> en ("<<i<<","<<j<<","<<n<<") : "<<imaHough_Norm(i,j,n)<<"\n"; }
										ImaTH_deja(iRho,iphi)=1;
									}
								}
								if (r!=0) {
									if (Rho<0) {Rho2min=maxi(Rho-2.f,r); Rho2max=mini(Rho+2.f,0.f);}
									else {Rho2min=maxi(Rho-2.f,0.f); Rho2max=mini(Rho+2.f,r);} 
									for (Rho2=Rho2min; Rho2<=Rho2max; Rho2+=pas_Rho) {
										iRho2=(int)(N_Rho/2+sign_Rho*Rho2); iRho=(int)(N_Rho/2-sign_Rho*Rho2); 
										phi2=(float)modulo2PI(PI/2-(-acos(Rho2/r)+Theta_rad)); iphi=(int)(phi2/PI*180.);
										if (!ImaTH_deja(iRho,iphi)) {
//											imaHough_Norm(iRho,iphi,0)++;
											if (iRho==i && iphi==j) n=0; 
											else if (iRho==i && iphi==j+1) n=1;
												else if (iRho==i+1 && iphi==j) n=2;
													else if (iRho==i && iphi==j-1) n=3;
														else if (iRho==i-1 && iphi==j) n=4;
															else if (iRho==i+1 && iphi==j+1) n=5;
																else if (iRho==i+1 && iphi==j-1) n=6;
																	else if (iRho==i-1 && iphi==j-1) n=7;
																		else if (iRho==i-1 && iphi==j+1) n=8;
																			else if (iRho==i && iphi==j+2) n=9;
																				else if (iRho==i+2 && iphi==j) n=10;
																					else if (iRho==i && iphi==j-2) n=11;
																						else if (iRho==i-2 && iphi==j) n=12;
																							else n=-1;
											if (n>=0) imaHough_Norm(i,j,n)++; 
if (n>=0 && iRho==30 && iphi==31) cout<<" point ("<<iRho<<","<<iphi<<") renforce par ("<<x<<","<<y<<") -> en ("<<i<<","<<j<<","<<n<<") : "<<imaHough_Norm(i,j,n)<<"\n";
											ImaTH_deja(iRho,iphi)=1;}
										iphi=(int)(modulo2PI(phi2+PI)/PI*180.f);
										if (!ImaTH_deja(iRho2,iphi)) {
//											imaHough_Norm(iRho2,iphi,0)++;
											if (iRho2==i && iphi==j) n=0; 
											else if (iRho2==i && iphi==j+1) n=1;
												else if (iRho2==i+1 && iphi==j) n=2;
													else if (iRho2==i && iphi==j-1) n=3;
														else if (iRho2==i-1 && iphi==j) n=4;
															else if (iRho2==i+1 && iphi==j+1) n=5;
																else if (iRho2==i+1 && iphi==j-1) n=6;
																	else if (iRho2==i-1 && iphi==j-1) n=7;
																		else if (iRho2==i-1 && iphi==j+1) n=8;
																			else if (iRho2==i && iphi==j+2) n=9;
																				else if (iRho2==i+2 && iphi==j) n=10;
																					else if (iRho2==i && iphi==j-2) n=11;
																						else if (iRho2==i-2 && iphi==j) n=12;
																							else n=-1;
											if (n>=0) imaHough_Norm(i,j,n)++; 
if (n>=0 && iRho2==30 && iphi==31) cout<<" point ("<<iRho2<<","<<iphi<<") renforce par ("<<x<<","<<y<<") -> en ("<<i<<","<<j<<","<<n<<") : "<<imaHough_Norm(i,j,n)<<"\n";
											ImaTH_deja(iRho2,iphi)=1;}
										phi2=(float)modulo2PI(PI/2-(acos(Rho2/r)+Theta_rad-2*PI)); iphi=(int)(phi2/PI*180.);
										if (!ImaTH_deja(iRho,iphi)) {
//											imaHough_Norm(iRho,iphi,0)++;
											if (iRho==i && iphi==j) n=0; 
											else if (iRho==i && iphi==j+1) n=1;
												else if (iRho==i+1 && iphi==j) n=2;
													else if (iRho==i && iphi==j-1) n=3;
														else if (iRho==i-1 && iphi==j) n=4;
															else if (iRho==i+1 && iphi==j+1) n=5;
																else if (iRho==i+1 && iphi==j-1) n=6;
																	else if (iRho==i-1 && iphi==j-1) n=7;
																		else if (iRho==i-1 && iphi==j+1) n=8;
																			else if (iRho==i && iphi==j+2) n=9;
																				else if (iRho==i+2 && iphi==j) n=10;
																					else if (iRho==i && iphi==j-2) n=11;
																						else if (iRho==i-2 && iphi==j) n=12;
																							else n=-1;
											if (n>=0) imaHough_Norm(i,j,n)++; 
if (n>=0 && iRho==30 && iphi==31) cout<<" point ("<<iRho<<","<<iphi<<") renforce par ("<<x<<","<<y<<") -> en ("<<i<<","<<j<<","<<n<<") : "<<imaHough_Norm(i,j,n)<<"\n";
											ImaTH_deja(iRho,iphi)=1;}
										iphi=(int)(modulo2PI(phi2+PI)/PI*180.f);
										if (!ImaTH_deja(iRho2,iphi)) {
//											imaHough_Norm(iRho2,iphi,0)++;
											if (iRho2==i && iphi==j) n=0; 
											else if (iRho2==i && iphi==j+1) n=1;
												else if (iRho2==i+1 && iphi==j) n=2;
													else if (iRho2==i && iphi==j-1) n=3;
														else if (iRho2==i-1 && iphi==j) n=4;
															else if (iRho2==i+1 && iphi==j+1) n=5;
																else if (iRho2==i+1 && iphi==j-1) n=6;
																	else if (iRho2==i-1 && iphi==j-1) n=7;
																		else if (iRho2==i-1 && iphi==j+1) n=8;
																			else if (iRho2==i && iphi==j+2) n=9;
																				else if (iRho2==i+2 && iphi==j) n=10;
																					else if (iRho2==i && iphi==j-2) n=11;
																						else if (iRho2==i-2 && iphi==j) n=12;
																							else n=-1;
											if (n>=0) imaHough_Norm(i,j,n)++; 
if (n>=0 && iRho2==30 && iphi==31) cout<<" point ("<<iRho2<<","<<iphi<<") renforce par ("<<x<<","<<y<<") -> en ("<<i<<","<<j<<","<<n<<") : "<<imaHough_Norm(i,j,n)<<"\n";
											ImaTH_deja(iRho2,iphi)=1;}
									}
								}
								Ima_deja(x,y)=1;
							}
						}
					}
					r=r2;
				}
			}
		}
//	cout<<" transformee de Hough : axe des abscisses = angle des droites, de -180 a +180 degres";
//	cout<<", axe des ordonnees = distance de la droite au centre du domaine image, 0 (droites radiales au milieu de l'axe";
//	cout<<", valeurs croisantes de rho (de -rhomax a +rhomax) avec l'indice de ligne";
//	cout<<", valeurs positives/negatives fct de |rho| et theta: droite intersecte la colonne du milieu du domaine de Hough";
//	cout<<" en rho/cos(theta), et la ligne du milieu du domaine de Hough en rho/sin(theta)\n";
	imaHough_Norm.sauve_ImaBSQ("./norme13_lignes_Hough.bsq");
	imaHough_Norm.statbasic(1); 
	return imaHough_Norm;
}

/* ---------------------------------------------------------------------*/
/* Transformee de Hough : precision d'angle de 1 degre, de rho : 1 pixel*/
imadata<int> imabin::hough_transform ()
/* ---------------------------------------------------------------------*/ {
	const double /*eps=1.e-6*/ eps=0.1;
	const int yc=nbcol/2, xc=nblig/2, N_Rho=(int)sqrt((float)nblig*nblig+nbcol*nbcol), n_dx=10;
	int		iTheta, iRho, i,j,k;
	float   x,y,fTheta, fRho;

	imadata<int> imaHough((int)(N_Rho+1),360);

	cout<<" image de la transformee de Hough : "<<imaHough.nlig()<<" lignes et "<<imaHough.ncol()<<" colonnes\n";

	for (iRho=0; iRho<=N_Rho; iRho++)
		for (iTheta=0; iTheta<360; iTheta++)
			imaHough(iRho,iTheta)=0;
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++)
			if ((*this)(i,j)) {
				for (iTheta=-180; iTheta<180; iTheta++) {
					fTheta=(float)(iTheta*PI/180.);
					fRho=(j-yc)*sin(fTheta)+(i-xc)*cos(fTheta);
					if ((fRho>=-N_Rho/2) &&  (fRho<=N_Rho/2)) {
						imaHough((int)(fRho+N_Rho/2),(iTheta+180))++;
						}
					}
				}
	imaHough.statbasic();
//	float s=imaHough.moyI()+15*pow(imaHough.varI(),0.5);
	float s=(float)(imaHough.moyI()+10*pow(imaHough.varI(),0.5));
//	float s=imaHough.moyI()+5*pow(imaHough.varI(),0.5);
	cout<<" valeurs depassant le seuil "<<s<<"\n";
	int n=0;
	for (iRho=0; iRho<=N_Rho; iRho++)
		for (iTheta=0; iTheta<360; iTheta++)
			if (imaHough(iRho,iTheta)>s) { //cout<<iRho<<" "<<iTheta<<" ";
//				cout<<" presence d'un max en Rho="<<iRho-N_Rho/2<<" vs centre image, et theta="<<iTheta-180<<" deg. de valeur "<<imaHough(iRho,iTheta)<<"\n";
				n++;
			}
	val_pos *T=new val_pos[n];
	n=0;
	for (iRho=0; iRho<=N_Rho; iRho++)
		for (iTheta=0; iTheta<360; iTheta++)
			if (imaHough(iRho,iTheta)>s) { 
				T[n].pos=iRho*360+iTheta; T[n].val=imaHough(iRho,iTheta); n++;
			}
	tri_rapide(T,n);
	imabin imaReconsH(nblig,nbcol);
	for (k=0; k<n; k++) {
		cout<<" max en Rho="<<(int)(T[k].pos)/360-N_Rho/2<<" vs centre im., et theta="<<(int)(T[k].pos)%360-180<<" deg. de valeur "<<T[k].val<<"\n";
//		cout<<n<<" "<<T[k].pos<<" presence d'un max en Rho="<<(int)(T[k].pos)/360-N_Rho/2<<" vs centre image, et theta="<<(int)(T[k].pos)%360-180<<" deg. de valeur "<<T[k].val<<"\n";
		fTheta=((int)(T[k].pos)%360-180)*(float)PI/180.f;
		fRho=(float)(T[k].pos)/360-N_Rho/2;
		if (fabs(sin(fTheta))>eps && fabs(cos(fTheta))>eps) {
			cout<<" reconstruction dans l'image\n";
			for (i=0; i<nblig; i++) {
				for (j=0; j<n_dx; j++) {
					x=i+(float)j/n_dx;
					y=(fRho-(x-xc)*cos(fTheta))/sin(fTheta)+yc;
					if (y>0 && y<nbcol) imaReconsH(mini(around(x),nblig-1),mini(maxi(0,around(y)),nbcol-1))=1;
				}
			}
		}
	}
	if (T!=NULL) delete[] T;
	imaReconsH.imaunsignedchar().sauve_ImaPGM("reconstructionHough.pgm");
/*	s=imaHough.maxI()-(imaHough.maxI()-imaHough.minI())/4;
//	s=imaHough.maxI()-(imaHough.maxI()-imaHough.minI())/2;
	cout<<" valeurs depassant le seuil "<<s<<"\n";
	for (iRho=0; iRho<=N_Rho; iRho++)
		for (iTheta=0; iTheta<360; iTheta++)
			if (imaHough(iRho,iTheta)>s) { cout<<iRho<<" "<<iTheta<<" ";
				cout<<" presence d'un max en Rho="<<iRho-N_Rho/2<<" vs centre image, et theta="<<iTheta-180<<" deg. de valeur "<<imaHough(iRho,iTheta)<<"\n";
			}*/
	return imaHough;
}

/* ---------------------------------------------------------------------*/
/* Analyse de la Transformee de Hough et extraction des segments*/
//imadata<int> imabin::analyse_hough (imadata<int> hough) {
/* ---------------------------------------------------------------------*/
/*	int x, y, i,j,iRho, iTheta, rmax, tmax, maxAct=0, maxPrec, iter=0;
	const int xc=nblig/2, yc=nbcol/2, N_Rho=(int)sqrt(nblig*nblig+nbcol*nbcol);
	float fRho, fTheta, R, T;
	imadata<int> imres(nblig, nbcol);
	for (x=0; x<nblig; x++)
		for(y=0; y<nbcol; y++) imres(x,y)=0;
	hough.sauve_ImaPGM("houghINIT.pgm");

	do {
		if (iter>0) maxPrec=maxAct;
		maxAct=0 ;
		for (iRho=0; iRho<=N_Rho; iRho++)
			for (iTheta=0; iTheta<=360; iTheta++)
				if (hough(iRho, iTheta)>maxAct) {
					maxAct=hough(iRho,iTheta);
					rmax = iRho-N_Rho/2;
					tmax= iTheta-180;
				}
				cout<<"max="<<maxAct;
		if (maxAct>20) {
			cout<<"DROITE "<<iter << " : ";
			cout<<"Angle ="<<(tmax) <<" RHO :"<<(rmax)<<"\n";
		    fTheta	=    (tmax) *  PI		/180;
		   	fRho	=	(rmax);
			if ((tmax==90)|| (tmax==-90)) {
				if (tmax==90)  y=(int)(rmax+yc);
				if (tmax==-90) y=(int)(nbcol-(rmax+yc));
				if (y>=0 &&  y<nbcol) {
					for(x=0; x<nblig; x++) {
					 imres(x,y)=255;
					 for (iTheta=-180; iTheta<180; iTheta++) {
						T= (float)iTheta * PI/180;
						R= (y-yc) * sin(T)+ (x-xc) * cos(T);
							if( ( R>=-N_Rho/2) &&  (R<=N_Rho/2) ) {
								i=(int)(R+N_Rho/2);
								j=(iTheta+180);
								if (hough(i,j)>0) hough(i,j)--;
							}
					}
				 }
			}
		}
		else
			if (fabs(sin(fTheta))<0.05) {
				x=(int)(fRho+xc);
				if (x>=0 &&  x<nblig) {
					for(y=0; y<nbcol; y++) {
						imres(x,y)=255;
						for (iTheta=-180; iTheta<180; iTheta++) {
							T= (float)iTheta * PI/180;
							R= (y-yc) * sin(T)+ (x-xc) * cos(T);
							if ((R>=-N_Rho/2) &&  (R<=N_Rho/2)) {
								i=(int)(R+N_Rho/2);
								j=(iTheta+180);
								if (hough(i,j)>0) hough(i,j)--;
							}
						}
					}
				}
			}
			else {
				for (x=0; x<nblig; x++) {
					y=(int)(((-cos(fTheta)/sin(fTheta))*(x-xc)+ (fRho)/sin(fTheta)))+yc+0.5;
    				if ((y>=0) && (y<nbcol)) {
						imres(x,y)=255;
  						for (iTheta=-180; iTheta<180; iTheta++) {
							T= (float)iTheta * PI/180;
							R= (y-yc) * sin(T)+ (x-xc) * cos(T);
							if( ( R>=-N_Rho/2) &&  (R<=N_Rho/2) ) {
								i=(int)(R+N_Rho/2);
								j=(iTheta+180);
								if (hough(i,j)>0) hough(i,j)--;
							}
						}
					}
				}
			}
		}
		iter++;
	} while((iter<100) && (maxAct>20));
	cout<<"imres.ncol&lig "<<imres.ncol()<<" "<<imres.nlig()<<"\n";
	hough.sauve_ImaPGM("houghFIN.pgm");
	return(imres);
}*/

imadata<int> imabin::transformee_Hough_cercles(const int N_Rho) {
	const int yc=nbcol/2, xc=nblig/2;
	int		i,j,iRho,ii,jj;
	float   xRho;
	imadata<int> imaHough(nblig,nbcol,N_Rho+1);
	cout<<" image de la transformee de Hough : "<<imaHough.nlig()<<" lignes et "<<imaHough.ncol()<<" colonnes et "<<imaHough.ncanaux()<<" canaux\n";
	imaHough.mise_a_zero();
//	imabin Ima_deja(nblig,nbcol);
	float pas_Rho=1.f/*10.f*/;
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++)
			if ((*this)(i,j)) { //cout<<"("<<i<<","<<j<<") ";
				for (ii=0; ii<nblig; ii++)
					for (jj=0; jj<nbcol; jj++) {
						xRho=pow(pow((float)ii-i,2.f)+pow((float)jj-j,2.f),0.5f);
						iRho=(int)(xRho/pas_Rho);
						if (iRho>=0 && iRho<=N_Rho) imaHough(ii,jj,iRho)++;
					}
/*				iRhomaxi=mini(mini(maxi(i,nblig-i),maxi(j,nbcol-j)),N_Rho);
				for (iRho=1; iRho<iRhomaxi; iRho++) {
					xRho=iRho*pas_Rho;
					imin=maxi(i-(int)xRho,0); imax=mini(i+(int)xRho,nblig-1); jmin=maxi(j-(int)xRho,0); jmax=mini(j+(int)xRho,nbcol-1);
					Rho2=xRho*xRho;
/*					Ima_deja.mise_a_zero();
					for (ic=imin; ic<imax; ic++) {
						xx=Rho2-(float)(i-ic)*(float)(i-ic);
						if (xx>=0) {
							xx=pow(xx,0.5f); 
							jc=(int)xx+j;  // jc-j=sqrt(Rho2-(float)(i-ic)*(float)(i-ic))
							if (jc>=0 && jc<nbcol && !Ima_deja(ic,jc)) {imaHough(ic,jc,iRho)++; Ima_deja(ic,jc)=1;}
							jc=-(int)xx+j; // j-jc=sqrt(Rho2-(float)(i-ic)*(float)(i-ic))
							if (jc>=0 && jc<nbcol && !Ima_deja(ic,jc)) {imaHough(ic,jc,iRho)++; Ima_deja(ic,jc)=1;}
						}
					}
					for (jc=jmin; jc<jmax; jc++) {
						xx=Rho2-(float)(j-jc)*(float)(j-jc);
						if (xx>=0) {
							xx=pow(xx,0.5f); 
							ic=(int)xx+i;  // ic-i=sqrt(Rho2-(float)(j-jc)*(float)(j-jc))
							if (ic>=0 && ic<nblig && !Ima_deja(ic,jc)) {imaHough(ic,jc,iRho)++; Ima_deja(ic,jc)=1;}
							ic=-(int)xx+i; // i-ic=sqrt(Rho2-(float)(j-jc)*(float)(j-jc))
							if (ic>=0 && ic<nblig && !Ima_deja(ic,jc)) {imaHough(ic,jc,iRho)++; Ima_deja(ic,jc)=1;}
						}
					}*/
/*          for (ii=imin; ii<=imax; ii++)
            for (jj=jmin; jj<=jmax; jj++) {
              xx=(float)(i-ii)*(float)(i-ii)+(float)(j-jj)*(float)(j-jj);
//              if (fabs(xx-Rho2)<100) imaHough(ii,jj,iRho)++;
              if (fabs(xx-Rho2)<1.f) imaHough(ii,jj,iRho)++;
            }
				}*/
			}
	cout<<" transformee de Hough : axe des abscisses = ligne du centre du cercle";
	cout<<", axe des ordonnees = colonne du centre du cercle, axe des canaux = rayon du cercle\n";
	imaHough.statbasic(0);
	return imaHough;
}

imadata<int> imabin::transformee_Hough_squares(const int l_side_max) {
	int		i,j,ii,jj,iRho,iRho_min,iRho_max;
	imadata<int> imaHough(nblig,nbcol,l_side_max);
	cout<<" image de la transformee de Hough : "<<imaHough.nlig()<<" lignes et "<<imaHough.ncol()<<" colonnes et "<<imaHough.ncanaux()<<" canaux\n";
	imaHough.mise_a_zero();
//	imabin Ima_deja(nblig,nbcol);
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++)
			if ((*this)(i,j)) { cout<<"("<<i<<","<<j<<") ";
				for (ii=0; ii<=i; ii++)
					for (jj=0; jj<=j; jj++) {
						iRho_min=maxi(i-ii,j-jj); iRho_max=mini(mini(nblig-ii,nbcol-jj),l_side_max); 
						for (iRho=iRho_min; iRho<iRho_max; iRho++) imaHough(ii,jj,iRho)++;
					}
			}
	cout<<" transformee de Hough : axe des abscisses = ligne du coin Upper Left du rectangle";
	cout<<", axe des ordonnees = colonne du du coin UL, axe des canaux = longueur du cote du carre\n";
	imaHough.statbasic(0);
	return imaHough;
}

/* ---------------------------------------------------------------------
 Parametres coocurrence binaire
 K : demi largeur de la fenetre de calcul
 o : orientation : 0 : - , 1 : \ , 2: | , 3: /
 d : distance
---------------------------------------------------------------------*/
imadata<float> imabin::coocurrence (int K, int o, int d, unsigned short int p) {
	int i,j,ii,jj,i0,i2,j0,j2,p1, p2;
//	const int a=around(d/sqrt(2.));
	const int a=d;                              // cas de la 8-connexite
	imadata<float> ima(nblig,nbcol);
	matrice2D<float> h(2,2);
//	cout<<" matrice de coocurrence calculee pour d="<<d<<" et o="<<o<<"\n";
	if (d>K) {cout<<" distance trop élevee! => d = "<<K<<"\n "; d=K;}
	switch(o) {
		case 0 : for (i=0; i<nblig; i++)		// orientation 0 (-)
					for (j=0; j<nbcol; j++) {
						for (ii=0; ii<2; ii++)
							for (jj=0; jj<2; jj++) h(ii,jj)=0;
						i0=maxi(0,i-K); i2=mini(nblig-1,i+K);
						j0=maxi(0,j-K); j2=mini(nbcol-1-d,j+K);
						for (ii=i0; ii<=i2; ii++)
							for(jj=j0; jj<=j2; jj++) {
								p1=(*this)(ii,jj);
								p2=(*this)(ii,jj+d);
								h(p1,p2)+=1.;
							}
						switch (p) {
							case 0: ima(i,j)=(float)h.entropy(); break;
							case 1: ima(i,j)=(float)h.contrast(); break;
							case 2: ima(i,j)=(float)h.dissimilarity(); break;
							case 3: ima(i,j)=(float)h.homogeneity(); break;
							default: break;
						}
					}
				break;
		case 1 : for (i=0; i<nblig; i++)		   // orientation 1 (\)
					for (j=0; j<nbcol; j++) {
						for (ii=0; ii<2; ii++)
							for (jj=0; jj<2; jj++) h(ii,jj)=0;
						i0=maxi(0,i-K); i2=mini(nblig-1-a,i+K);
						j0=maxi(0,j-K); j2=mini(nbcol-1-a,j+K);
						for (ii=i0; ii<=i2; ii++)
							for(jj=j0; jj<=j2; jj++) {
								p1=(*this)(ii,jj);
								p2=(*this)(ii+a,jj+a);
								h(p1,p2)+=1.;
							}
						switch (p) {
							case 0: ima(i,j)=(float)h.entropy(); break;
							case 1: ima(i,j)=(float)h.contrast(); break;
							case 2: ima(i,j)=(float)h.dissimilarity(); break;
							case 3: ima(i,j)=(float)h.homogeneity(); break;
							default: break;
						}
					}
				break;
		case 2 : for (i=0; i<nblig; i++)		// orientation 2 (|)
					for (j=0; j<nbcol; j++) {
						for (ii=0; ii<2; ii++)
							for (jj=0; jj<2; jj++) h(ii,jj)=0;
						i0=maxi(0,i-K); i2=mini(nblig-1-d,i+K);
						j0=maxi(0,j-K); j2=mini(nbcol-1,j+K);
						for (ii=i0; ii<=i2; ii++)
							for(jj=j0; jj<=j2; jj++) {
								p1=(*this)(ii,jj);
								p2=(*this)(ii+d,jj);
								h(p1,p2)+=1.;
							}
						switch (p) {
							case 0: ima(i,j)=(float)h.entropy(); break;
							case 1: ima(i,j)=(float)h.contrast(); break;
							case 2: ima(i,j)=(float)h.dissimilarity(); break;
							case 3: ima(i,j)=(float)h.homogeneity(); break;
							default: break;
						}
					}
				break;
		case 3 : for (i=0; i<nblig; i++)		   // orientation 3 (/)
					for (j=0; j<nbcol; j++) {
						for (ii=0; ii<2; ii++)
							for (jj=0; jj<2; jj++) h(ii,jj)=0;
						i0=maxi(0,i-K); i2=mini(nblig-1-a,i+K);
						j0=maxi(0+a,j-K); j2=mini(nbcol-1,j+K);
						for (ii=i0; ii<=i2; ii++)
							for(jj=j0; jj<=j2; jj++) {
								p1=(*this)(ii,jj);
								p2=(*this)(ii+a,jj-a);
								h(p1,p2)+=1.;
							}
						switch (p) {
							case 0: ima(i,j)=(float)h.entropy(); break;
							case 1: ima(i,j)=(float)h.contrast(); break;
							case 2: ima(i,j)=(float)h.dissimilarity(); break;
							case 3: ima(i,j)=(float)h.homogeneity(); break;
							default: break;
						}
					}
				break;
		default: break;
	}
	return ima;
}

void imabin::trace_line (float fTheta, float fRho, float fThetaP1, float fThetaP2, bool iaff) {
	if (iaff) cout<<" trace_line fTheta = "<<fTheta*180/PI<<", fRho = "<<fRho<<", entre "<<fThetaP1*180/PI<<" et "<<fThetaP2*180/PI<<" deg.\n";
	if (fThetaP1>PI) fThetaP1-=(float)(2*PI); if (fThetaP2>PI) fThetaP2-=(float)(2*PI);
	if (iaff) cout<<" trace_line fTheta = "<<fTheta*180/PI<<", fRho = "<<fRho<<", entre "<<fThetaP1*180/PI<<" et "<<fThetaP2*180/PI<<" deg.\n";
	bool add2PI=0;
/*	if (maxi(fThetaP1,fThetaP2)-mini(fThetaP1,fThetaP2)>PI) {
		if (fThetaP2==mini(fThetaP1,fThetaP2)) fThetaP2+=2*PI; else fThetaP1+=2*PI; 
		add2PI=1; }*/
	const float fThetaMin=mini(fThetaP1,fThetaP2), fThetaMax=maxi(fThetaP1,fThetaP2); if (iaff) cout<<", i.e. entre "<<fThetaMin*180/PI<<" et "<<fThetaMax*180/PI<<" deg.\n";
	const int n_dx=10, yc=nblig/2, xc=nbcol/2;
	int i,j,ii,jj;
	float x,y,phi; //double cos_phi,sin_phi;
	if (fabs(cos(fTheta))>0.5) {if (iaff) cout<<" fTheta "<<fTheta/PI*180<<" col. vs lig. / pente = "<<sin(fTheta)/cos(fTheta)<<"\n";
		for (i=0; i<nblig; i++)
			for (j=0; j<n_dx; j++) { // y est la coord. lig. (orientée vers le bas) et x la coord. col.
				y=-(i+(float)j/n_dx); x=(fRho-(y+yc)*sin(fTheta))/cos(fTheta)+xc; 
				phi=atan2(y+yc,x-xc); //cos_phi=cos(phi); sin_phi=sin(phi); 
				if (phi<0 && add2PI) phi+=(float)(2*PI); 
				ii=mini(maxi(0,around(-y)),nblig-1); jj=mini(maxi(0,around(x)),nbcol-1); 
				if (around(x)>=0 && around(x)<nbcol && phi>=fThetaMin && phi<=fThetaMax) (*this)(ii,jj)=true;
//				if (around(x)>=0 && around(x)<nbcol && j==0) cout<<" col = "<<x<<" lig = "<<-y<<" phi = "<<phi*180/PI<<"\n";
			}
	} else {if (iaff) cout<<" fTheta "<<fTheta/PI*180<<" lig. vs col. / pente = "<<cos(fTheta)/sin(fTheta)<<"\n";
		for (j=0; j<nbcol; j++)
			for (i=0; i<n_dx; i++) {
				x=j+(float)i/n_dx; y=(fRho-(x-xc)*cos(fTheta))/sin(fTheta)-yc; // x la coord. col. et y est la coord. lig. (orientée vers le bas)
				phi=atan2(y+yc,x-xc); //cos_phi=cos(phi); sin_phi=sin(phi);
				if (phi<0 && add2PI) phi+=(float)(2*PI); 
				ii=mini(maxi(0,around(-y)),nblig-1); jj=mini(maxi(0,around(x)),nbcol-1);
				if (around(-y)>=0 && around(-y)<nblig && phi>=fThetaMin && phi<=fThetaMax) (*this)(ii,jj)=true;
			}
	} if (iaff) cout<<" fin trace_line : segment de longueur "<<norm()<<" pixels\n";
}

void imabin::trace_line (float fTheta, float fRho, int lmin, int lmax, int cmin, int cmax) {
//	cout<<" trace_line2 fTheta = "<<fTheta*180/PI<<", fRho = "<<fRho<<", entre ("<<lmin<<","<<cmin<<") et ("<<lmax<<","<<cmax<<")\n";
	const int n_dx=10, yc=nblig/2, xc=nbcol/2;
	int i,j,ii,jj;
	float x,y;
	if (lmin>lmax) {i=lmax; lmax=lmin; lmin=i;}
	if (cmin>cmax) {j=cmax; cmax=cmin; cmin=j;}
	if (fabs(cos(fTheta))>0.5) {//cout<<" fTheta "<<fTheta/PI*180<<" col. vs lig / pente = "<<sin(fTheta)/cos(fTheta)<<"\n";
		for (i=lmin; i<=lmax; i++)
			for (j=0; j<n_dx; j++) { // y est la coord. lig. (orientée vers le bas) et x la coord. col.
				y=-(i+(float)j/n_dx); x=(fRho-(y+yc)*sin(fTheta))/cos(fTheta)+xc; 
				ii=mini(maxi(0,around(-y)),nblig-1); jj=mini(maxi(0,around(x)),nbcol-1); 
				if (around(x)>=cmin && around(x)<=cmax) (*this)(ii,jj)=true;
			}
	} else {//cout<<" fTheta "<<fTheta/PI*180<<" col. vs lig / pente = "<<cos(fTheta)/sin(fTheta)<<"\n";
		for (j=cmin; j<=cmax; j++) {//cout<<" j = "<<j;
			for (i=0; i<n_dx; i++) {
				x=j+(float)i/n_dx; y=(fRho-(x-xc)*cos(fTheta))/sin(fTheta)-yc; // x la coord. col. et y est la coord. lig. (orientée vers le bas)
				ii=mini(maxi(0,around(-y)),nblig-1); jj=mini(maxi(0,around(x)),nbcol-1);
				if (around(-y)>=lmin && around(-y)<=lmax) (*this)(ii,jj)=true;
			}
		}
	} //cout<<" fin trace_line : segment de longueur "<<norm()<<" pixels\n";
}

void imabin::trace_line (int l1, int c1, int l2, int c2) {//cout<<" entree dans trace_line avec deja "<<norm()<<" pixels\n"; imaunsignedchar().sauve_ImaPGM("bid0.ppgm");
	const int n_dx=10;
	int xc=nbcol/2, yc=nblig/2, x1=c1-xc, y1=yc-l1, x2=c2-xc, y2=yc-l2, lmin=mini(l1,l2),lmax=maxi(l1,l2),cmin=mini(c1,c2),cmax=maxi(c1,c2),i,j,ii,jj,n=(int)norm();
	float x,y, fRho=FLT_MAX, fTheta=atan2((float)x2-x1,(float)y1-y2);	// y est la coord. ligne (axe des y vers le haut) et x la coord. colonne, l'angle est dans le sens trigo à partir de l'axe x

	if (y1!=y2 || x1!=x2) fRho=(float)fabs((float)x1*y2-x2*y1)/pow(pow((float)y1-y2,2.f)+pow((float)x1-x2,2.f),0.5f);
	else {cout<<" points P1 et P2 confondus !!!\n";}
	if (x1*y2-x2*y1>0) fTheta+=(float)PI;

	if (fabs(cos(fTheta))>0.5) {//cout<<" fTheta "<<fTheta/PI*180<<" col. vs lig / pente = "<<sin(fTheta)/cos(fTheta)<<"\n";
		for (i=lmin; i<=lmax; i++)
			for (j=0; j<n_dx; j++) { // y est la coord. lig. (orientée vers le bas) et x la coord. col.
				y=-(i+(float)j/n_dx); x=(fRho-(y+yc)*sin(fTheta))/cos(fTheta)+xc; 
				ii=mini(maxi(0,around(-y)),nblig-1); jj=mini(maxi(0,around(x)),nbcol-1); 
				if (around(x)>=cmin && around(x)<=cmax) {/*cout<<((*this)(ii,jj)==true?"alr ":"new ")<<ii<<" "<<jj<<" ; ";*/ (*this)(ii,jj)=true; }
			}
	} else {//cout<<" fTheta "<<fTheta/PI*180<<" col. vs lig / pente = "<<cos(fTheta)/sin(fTheta)<<"\n";
		for (j=cmin; j<=cmax; j++) {//cout<<" j = "<<j;
			for (i=0; i<n_dx; i++) {
				x=j+(float)i/n_dx; y=(fRho-(x-xc)*cos(fTheta))/sin(fTheta)-yc; // x la coord. col. et y est la coord. lig. (orientée vers le bas)
				ii=mini(maxi(0,around(-y)),nblig-1); jj=mini(maxi(0,around(x)),nbcol-1);
				if (around(-y)>=lmin && around(-y)<=lmax) {/*cout<<((*this)(ii,jj)==true?"alr ":"new ")<<ii<<" "<<jj<<" ; ";*/ (*this)(ii,jj)=true; }
			}
		}
	} cout<<" fin trace_line : segment ("<<c1<<","<<l1<<") -> ("<<c2<<","<<l2<<") de longueur "<<norm()-n<<" pixels\n"; //imaunsignedchar().sauve_ImaPGM("bid1.ppgm"); {char aa; cin>>aa;}
}