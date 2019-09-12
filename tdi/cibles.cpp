#include "cibles.h"

const int cible::valnul=-999;

void cible::affiche() const {
	cout<<" cible "<<nocib<<"\n";
	if (!pos_valid)
		cout<<" centre de la cible inconnu\n";
	else
		cout<<" coordonnees du centre : ("<<x_centre<<","<<y_centre<<")\n";
	if (!dim_valid)
		cout<<" rectangle englobant de la cible inconnu\n";
	else {
		cout<<" dimensions du rectangle englobant : "<<dx_rect<<"x"<<dy_rect<<"\n";
		cout<<" position du coin en haut a gauche : ("<<x_ul_rect<<","<<y_ul_rect<<")\n";
		}
	if (!vit_valid)
		cout<<" vitesse de la cible inconnue\n";
	else
		cout<<" vitesse : v=("<<x_dir_vit<<","<<y_dir_vit<<"), ||v||="<<mod_vit<<"\n";
	}

double cible::distance_cibles (cible &cib) {
	double xx1=0., xx2=0.;
	if (pos_valid && cib.centre_valide()) xx1+=pow((double)xcentre()-cib.xcentre(),2.)+pow((double)ycentre()-cib.ycentre(),2.);
	if (dim_valid && cib.rectangle_valide()) {
		xx2+=pow((double)dx_rectangle()-cib.dx_rectangle(),2.)+pow((double)dy_rectangle()-cib.dy_rectangle(),2.);
//		xx+=pow(x_ul_rectangle()-cib.x_ul_rectangle(),2)+pow(y_ul_rectangle()-cib.y_ul_rectangle(),2);
		}
	return xx1+xx2*10;
	}


void liste_cibles::insere (cible &cib) {
//	cout<<" insertion nouvel element => #elts = ";
	elt_cible *nouveau;
	nouveau=new elt_cible;
	nouveau->cib=cib;
	nouveau->suivant=debut;
	debut=nouveau;
	nbelt++;
//	cout<<nbelt<<"\n";
	}

void liste_cibles::insere (cible *cib) {
//	cout<<" insertion nouvel element => #elts = ";
	elt_cible *nouveau;
	nouveau=new elt_cible;
	nouveau->cib=*cib;
	nouveau->suivant=debut;
	debut=nouveau;
	nbelt++;
//	cout<<nbelt<<"\n";
	}

cible liste_cibles::extrait () {
	cible cib0;
	elt_cible res;
	res.cib=cib0; res.suivant=NULL;
	if (debut!=NULL) {
		res=*debut;
		delete debut;
		debut=res.suivant;
		nbelt--;
		}
	return res.cib;
	}

void liste_cibles::supprime (int icib) {
	courant=debut;
	precedent=NULL;
	bool trouve=0;
	cible cibc;
	while (!trouve && courant!=NULL) {
		cibc=courant->cib;
		if (cibc.no_cible()==icib) {
			trouve=1;
			if (precedent!=NULL) {
				precedent->suivant=courant->suivant;
				delete courant;
			}
			else
				if (courant==debut) {
					debut=courant->suivant;
					delete courant;
				}
				else
					cout<<" Pb pas de precedent\n";
		}
		else {
			precedent=courant;
			courant=courant->suivant;
		}
	}
	if (trouve) nbelt--;
}

bool liste_cibles::existe (cible &cib) {
	courant=debut;
	bool trouve=0;
	cible cibc;
	while (!trouve && courant!=NULL) {
		cibc=courant->cib;
		if (cibc==cib) {
			cout<<" trouve : \n"; cibc.affiche(); cout<<" = \n"; cib.affiche();
			trouve=1;
			}
		else courant=courant->suivant;
		}
	return trouve;
	}

bool liste_cibles::existe (cible *cib) {
	courant=debut;
	bool trouve=0;
	cible cibc;
	while (!trouve && courant!=NULL) {
		cibc=courant->cib;
		if (cibc==*cib) {
			cout<<" trouve : \n"; cibc.affiche(); cout<<" = \n"; cib->affiche();
			trouve=1;
			}
		else courant=courant->suivant;
		}
	return trouve;
	}

bool liste_cibles::existe (int icib) {
	courant=debut;
	bool trouve=0;
	cible cibc;
	while (!trouve && courant!=NULL) {
		cibc=courant->cib;
		if (cibc.no_cible()==icib) trouve=1;
		else courant=courant->suivant;
		}
	return trouve;
	}

void liste_cibles::affiche () {
	courant=debut;
	cible cibc;
	while (courant!=NULL) {
		cibc=courant->cib;
		cibc.affiche();
		courant=courant->suivant;
		}
	}

int liste_cibles::label_max () {
	int n=0;
	courant=debut;
	cible cibc;
	while (courant!=NULL) {
		cibc=courant->cib;
		if (cibc.no_cible()>n) n=cibc.no_cible();
		courant=courant->suivant;
		}
	return n;
	}

void liste_cibles::correspondances_cibles (liste_cibles &L2, bool **Tcorresp, bool affich) {
//	double coef_beta=1.-0.75;
//	int ,beta=coef_beta*nbrlig*nbrcol/maxi(c1,c2);
	double coef_iter0=20., pmvt=0.8, pdeath=(1.-pmvt)/2., pbirth=pdeath+(1.-pmvt)/2.;
	int c1=nb_elts(), c2=L2.nb_elts(), cc=0, n_iter0=(int)(coef_iter0*c1*c2), n_itermax=100*n_iter0, ntiragesmax=1000;
	int i,j,ii,jj,n_iter,ndeath,nbirth,nswap,nch=0,ntirages,n,n1,n2,l1,l2;
	bool accept_chg;
	double beta=0.,Temp0=0.,rTemp=0.9,Temp,U2minim,xx,x,U,c;
	if (affich) cout<<" liste 1 : "<<c1<<" cibles, comparee a liste 2 : "<<c2<<" cibles\n";
	bool **Tcorres2=NULL;
	Tcorres2=new bool*[c1+1];
	for (i=0; i<=c1; i++)
		Tcorres2[i]=new bool[c2+1];
	double** D2d=NULL;
	D2d=new double*[c1+1];
	for (i=0; i<=c1; i++) D2d[i]=new double[c2+1];
	cible *T1cib=new cible[nb_elts()], *T2cib=new cible[L2.nb_elts()];
	int *Tdeath[2], *Tbirth[2], *Tswap[2];
	for (i=0; i<2; i++) {
		Tdeath[0]=new int[c1*c2];
		Tdeath[1]=new int[c1*c2];
		Tbirth[0]=new int[c1*c2];
		Tbirth[1]=new int[c1*c2];
		Tswap[0]=new int[c1*c2];
		Tswap[1]=new int[c1*c2];
		}
	if (affich) cout<<" liste 1 :\n"; (*this).affiche();
	if (affich) cout<<" liste 2 :\n"; L2.affiche();
// initialisations
	i=0;
	while (nb_elts()>0) T1cib[i++]=extrait();
	j=0;
	while (L2.nb_elts()>0) T2cib[j++]=L2.extrait();
	for (i=0; i<=c1; i++)
		for (j=0; j<=c2; j++) D2d[i][j]=0.;
	for (i=1; i<=c1; i++) {
		xx=FLT_MAX;
		for (j=1; j<=c2; j++) {
			D2d[i][j]=T1cib[i-1].distance_cibles(T2cib[j-1]);
			if (D2d[i][j]<xx) xx=D2d[i][j];
			}
		x=FLT_MAX;
		for (j=1; j<=c2; j++)
			if (D2d[i][j]<x && D2d[i][j]>xx) x=D2d[i][j];
		beta+=x;
		}
	for (j=1; j<=c2; j++) {
		xx=FLT_MAX;
		for (i=1; i<=c1; i++)
			if (D2d[i][j]<xx) xx=D2d[i][j];
		x=FLT_MAX;
		for (i=1; i<=c1; i++)
			if (D2d[i][j]<x && D2d[i][j]>xx) x=D2d[i][j];
		beta+=x;
		}
	beta/=(c1+c2);
	cout<<" beta = "<<beta<<"\n";
/*	beta=0;
	for (i=1; i<=c1; i++) {
		xx=FLT_MAX;
		for (j=1; j<=c2; j++)
			if (D2d[i][j]<xx) xx=D2d[i][j];
		beta+=xx;
		}
	for (j=1; j<=c2; j++) {
		xx=FLT_MAX;
		for (i=1; i<=c1; i++)
			if (D2d[i][j]<xx) xx=D2d[i][j];
		beta+=xx;
		}
	beta/=(c1+c2);
	cout<<" beta = "<<beta<<"\n";
	beta=0;
	for (j=1; j<=c2; j++)
		for (i=1; i<=c1; i++) beta+=D2d[i][j];
	beta/=(c1*c2);
	cout<<" beta = "<<beta<<"\n";*/

	cout<<"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n matrice des distances :\n";
	for (i=1; i<=c1; i++) {
		for (j=1; j<=c2; j++) cout<<" "<<setw(8)<<D2d[i][j];
		cout<<"\n";
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
	cout<<"\n******************************************\n matrice des correspondances initiale :\n";
	for (i=1; i<=c1; i++) {
		for (j=1; j<=c2; j++) cout<<" "<<Tcorresp[i][j];
		cout<<"\n";
		}
//	char aa; cin>>aa;
	U2minim=0.;
	for (i=1; i<=c1; i++)
		for (j=1; j<=c2; j++) 
			if (Tcorresp[i][j]) U2minim+=D2d[i][j];
	cc=mini(c1,c2);
	U2minim-=beta*cc;
	cout<<" energie = "<<U2minim<<"\n";
	Temp=Temp0;
	n_iter=0;
// test d'une nouvelle configuration : 
	do {
		n_iter++;
		Temp*=rTemp;
		if (n_iter%n_iter0==1) nch=0;
		for (i=1; i<=c1; i++)
			for (j=1; j<=c2; j++) Tcorres2[i][j]=Tcorresp[i][j];
// tirage d'un mouvement
		ndeath=0; nbirth=0; nswap=0;
		for (i=1; i<=c1; i++)
			for (j=1; j<=c2; j++) {
				if (Tcorresp[i][j]) {
					Tdeath[0][ndeath]=i;
					Tdeath[1][ndeath]=j;
					ndeath++;
					n1=0;
					for (ii=1; ii<=c1; ii++) n1+=Tcorresp[ii][j];
					n2=0;
					for (jj=1; jj<=c2; jj++) n2+=Tcorresp[i][jj];
					if (n1==1 && n2==1) {
						Tswap[0][nswap]=i;
						Tswap[1][nswap]=j;
						nswap++;
						}
					}
				else {
					n1=0;
					for (ii=1; ii<=c1; ii++) n1+=Tcorresp[ii][j];
					n2=0;
					for (jj=1; jj<=c2; jj++) n2+=Tcorresp[i][jj];
					if (n1==0 || n2==0) {
						Tbirth[0][nbirth]=i;
						Tbirth[1][nbirth]=j;
						nbirth++;
						}
					}
				}
//		cout<<" ndeath "<<ndeath<<" nbirth "<<nbirth<<" nswap "<<nswap<<"\n";
		ntirages=0;
		do {
			xx=(double)rand()/(double)RAND_MAX;
			ntirages++;
			} while ( !((xx<pdeath && ndeath>0) || (xx>=pdeath && xx<pbirth && nbirth>0) || 
						(xx>=pbirth && nswap>1)) && ntirages<ntiragesmax);
		if (ntirages>=ntiragesmax) {
			cout<<" nbre de tirages max atteint dans tirage du mouvement\n";
			accept_chg=0;
			}
		else {
			accept_chg=1;
			if (xx>=0 && xx<pdeath) { // suppression d'une relation de correspondance
				n=(rand()%ndeath);
				i=Tdeath[0][n];
				j=Tdeath[1][n];
				if (Tcorresp[i][j]) Tcorresp[i][j]=0;
				else cout<<" suppression de la correspondance "<<i<<" "<<j<<" mais Tcorresp[i][j] deja nul ????\n";
//				cout<<" suppression de la correspondance "<<i<<" "<<j<<"\n";
			}
			if (xx>=pdeath && xx<pbirth) { // ajout d'une relation de correspondance
				n=(rand()%nbirth);
				i=Tbirth[0][n];
				j=Tbirth[1][n];
				if (!Tcorresp[i][j]) Tcorresp[i][j]=1;
				else cout<<" ajout de la correspondance "<<i<<" "<<j<<" mais Tcorresp[i][j] deja =1 ????\n";
//				cout<<" ajout de la correspondance "<<i<<" "<<j<<"\n";
				}
			if (xx>=pbirth) {
				n1=(rand()%nswap);
				do {
					n2=(rand()%nswap);
					ntirages++;
				} while ( n1==n2 && ntirages<ntiragesmax);
				if (ntirages>=ntiragesmax) {
					cout<<" nbre de tirages max atteint dans cas swap\n";
					accept_chg=0;
					}
				else {
					i=Tswap[0][n1];
					j=Tswap[1][n1];
					ii=Tswap[0][n2];
					jj=Tswap[1][n2];
					if (Tcorresp[i][j] && Tcorresp[ii][jj] && !Tcorresp[i][jj] && !Tcorresp[ii][j]) {
						Tcorresp[i][j]=0;
						Tcorresp[ii][jj]=0;
						Tcorresp[ii][j]=1;
						Tcorresp[i][jj]=1;
						}
					else cout<<" swap entre les correspondances "<<i<<" "<<j<<" et "<<ii<<" "<<jj<<" mais ????\n";
					}
//				cout<<" swap entre les correspondances "<<i<<" "<<j<<" et "<<ii<<" "<<jj<<" \n";
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
						if (n1>1 && n2>1) {
							cout<<" condition (3) non respectee en ("<<l1<<","<<l2<<") : \n";
							cout<<" => mise a "<<!Tcorresp[l1][l2]<<" du terme ("<<l1<<","<<l2<<")\n";
							Tcorresp[l1][l2]=!Tcorresp[l1][l2];
							}
						}
					}
			U=0;
			for (i=1; i<=c1; i++)
				for (j=1; j<=c2; j++) 
					if (Tcorresp[i][j]) U+=D2d[i][j];
			c=0;
/*			for (i=1; i<=c1; i++) {
				n1=0;
				for (j=1; j<=c2; j++) n1+=Tcorresp[i][j];
				if (n1==0) c+=1.;
				}
			for (j=1; j<=c2; j++) {
				n2=0;
				for (i=1; i<=c1; i++) n2+=Tcorresp[i][j];
				if (n2==0) c+=1.;
				}*/
			for (i=1; i<=c1; i++) {
				n1=0;
				jj=0;
				for (j=1; j<=c2; j++) {
					n1+=Tcorresp[i][j];
					if (jj==0 && n1==1) jj=j;
					}
				if (n1>1) c+=1.;
				else {
					if (n1==1) {
						n2=0;
						for (ii=1; ii<=c1; ii++) n2+=Tcorresp[ii][jj];
						if (n2==0) {
							cout<<" somme des termes de la ligne "<<i<<" = "<<n1<<"\n";
							cout<<" somme des termes de la colonne "<<jj<<" nulle\n";
	for (int iii=1; iii<=c1; iii++) {
		for (int jjj=1; jjj<=c2; jjj++) cout<<" "<<Tcorresp[iii][jjj];
		cout<<"\n";
		}
		char aa; cin>>aa;
						}
						else c+=1./(double)n2;
						}
					}
				}
//			cout<<" # de classes 'communes' = "<<c<<"\n";
			if (fabs(c-around(c))>1.e-6) 
				cout<<" pb dans le calcul du # de classes 'communes' = "<<c<<"!="<<around(c)<<"\n";
			U-=beta*c;
			if (U>=U2minim && (Temp==0 ||
				(double)rand()/(double)RAND_MAX>=exp((U2minim-U)/Temp)) ) accept_chg=0;
			}
		if (!accept_chg) {
			for (i=1; i<=c1; i++)
				for (j=1; j<=c2; j++) Tcorresp[i][j]=Tcorres2[i][j];
//			cout<<" changement refuse\n";
			}
		else {
			nch++;
			U2minim=U;
			cc=around(c);
//			cout<<" changement accepte\n";
			}
/*	cout<<"\n******************************************\n matrice des correspondances courante :\n";
	for (i=1; i<=c1; i++) {
		for (j=1; j<=c2; j++) cout<<" "<<Tcorresp[i][j];
		cout<<"\n";
		}
	cout<<" energie = "<<U2minim<<"\n";*/
//	char aa; cin>>aa;
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
		cout<<" energie = "<<U2minim<<"\n";
//		cout<<" % de pixels de 'meme' label = "<<100.*(U2minim-beta*cc)/nbrlig/nbrcol<<" %\n";
		}
	for (i=0; i<=c1; i++) 
		if (D2d[i]!=NULL) delete[] D2d[i];
	if (D2d!=NULL) delete[] D2d;
	for (i=0; i<=c1; i++)
		if (Tcorres2[i]!=NULL) delete[] Tcorres2[i];
	if (Tcorres2!=NULL) delete[] Tcorres2;
	if (T1cib!=NULL) delete[] T1cib;
	if (T2cib!=NULL) delete[] T2cib;
	for (i=0; i<2; i++) {
		if (Tdeath[i]!=NULL) delete[] Tdeath[i];
		if (Tbirth[i]!=NULL) delete[] Tbirth[i];
		if (Tswap[i]!=NULL) delete[] Tswap[i];
		}
	}
 
