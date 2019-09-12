#include "liste_pixels.h"

liste_pixels& liste_pixels::operator += (liste_pixels& L) {
//	cout<<" fin a @ "<<fin<<" et suivant @ "<<fin->suivant;
	fin->suivant=L.debut; //cout<<" devient @ "<<fin->suivant<<"\n";
	fin=L.fin; 	
	nbelt+=L.nbelt;
	irect=irect&&L.irect;
	if (irect) {
		xmin=mini(xmin,L.xmin); xmax=maxi(xmax,L.xmax); 
		ymin=mini(ymin,L.ymin); ymax=maxi(ymax,L.ymax); 
	}
	L.nbelt=0; L.debut=L.courant=L.precedent=L.fin=NULL; L.irect=0; // les elements de L ont changes de liste proprietaire
	return *this;
}

void liste_pixels::insere (const elt_liste &E) {
//	cout<<" insertion nouvel element => #elts = ";
	elt_liste *nouveau=new elt_liste;
	nouveau->x=E.x;
	nouveau->y=E.y;
	nouveau->suivant=debut;
	debut=nouveau;
	nbelt++;
	if (irect) {
		if (xmin<0 || E.x<xmin) xmin=E.x;
		if (E.x>xmax) xmax=E.x;
		if (ymin<0 || E.y<ymin) ymin=E.y;
		if (E.y>ymax) ymax=E.y;
	}
	if (nbelt==1) fin=debut;
}

void liste_pixels::insere (const elt_liste *E) {
//	cout<<" insertion nouvel element => #elts = ";
	elt_liste *nouveau=new elt_liste;
	nouveau->x=E->x;
	nouveau->y=E->y;
	nouveau->suivant=debut;
	debut=nouveau;
	nbelt++;
	if (irect) {
		if (xmin<0 || E->x<xmin) xmin=E->x;
		if (E->x>xmax) xmax=E->x;
		if (ymin<0 || E->y<ymin) ymin=E->y;
		if (E->y>ymax) ymax=E->y;
	}
	if (nbelt==1) fin=debut;
}

void liste_pixels::insere (const int x0, const int y0) {
	elt_liste *nouveau=new elt_liste;
	nouveau->x=x0;
	nouveau->y=y0;
	nouveau->suivant=debut;
	debut=nouveau;
	nbelt++;
	if (irect) {
		if (xmin<0 || x0<xmin) xmin=x0;
		if (x0>xmax) xmax=x0;
		if (ymin<0 || y0<ymin) ymin=y0;
		if (y0>ymax) ymax=y0;
	}
	if (nbelt==1) fin=debut;
}

elt_liste liste_pixels::extrait () {
	elt_liste res={-1,-1,NULL};
	if (debut!=NULL) {
		res=*debut;
		delete debut;
		debut=res.suivant;
		nbelt--;
	}
//	rectangle_englobant();
	irect=0;
	return res;
}

void liste_pixels::supprime (int x0, int y0) {
	courant=debut;
	precedent=NULL;
	bool trouve;
	trouve=0;
	while (!trouve && courant!=NULL) {
		if (courant->x==x0 && courant->y==y0) {
			trouve=1;
			if (precedent!=NULL) {
				if (courant==fin) fin=precedent;
				precedent->suivant=courant->suivant;
				delete courant;
			} else
				if (courant==debut) {
					debut=courant->suivant;
					if (courant==fin) fin=debut;
					delete courant;
				}
				else cout<<" Pb pas de precedent\n";
		} else {
			precedent=courant;
			courant=courant->suivant;
		}
	}
	if (trouve) {
		nbelt--;
//		rectangle_englobant();
		irect=0;
	}
}

bool liste_pixels::existe (elt_liste &E) {
	courant=debut;
	bool trouve;
	trouve=0;
	while (!trouve && courant!=NULL) {
		if (courant->x==E.x && courant->y==E.y) {
			cout<<" trouve : "<<courant->x<<" = "<<E.x;
			cout<<" "<<courant->y<<" = "<<E.y;
			trouve=1;
		} else courant=courant->suivant;
	}
	return trouve;
}

bool liste_pixels::existe (elt_liste *E) { // ou elt_liste& E
	courant=debut;
	bool trouve;
	trouve=0;
	while (!trouve && courant!=NULL) {
		if (courant->x==E->x && courant->y==E->y) {
			cout<<" trouve : "<<courant->x<<" = "<<E->x;
			cout<<" "<<courant->y<<" = "<<E->y;
			trouve=1;
		} else courant=courant->suivant;
	}
	return trouve;
}

bool liste_pixels::existe (int x0, int y0) {
	courant=debut;
	bool trouve;
	trouve=0;
	while (!trouve && courant!=NULL) {
		if (courant->x==x0 && courant->y==y0) trouve=1;
		else courant=courant->suivant;
	}
	return trouve;
}

void liste_pixels::affiche () const {
	elt_liste *courant2=debut;
	while (courant2!=NULL) {
		cout<<" ("<<courant2->x<<","<<courant2->y<<")";
		courant2=courant2->suivant;
	}
}

void liste_pixels::rectangle_englobant () {
	courant=debut;
	xmax=ymax=INT_MIN;
	xmin=ymin=INT_MAX;
	while (courant!=NULL) {
		if (courant->x<xmin) xmin=courant->x;
		if (courant->x>xmax) xmax=courant->x;
		if (courant->y<ymin) ymin=courant->y;
		if (courant->y>ymax) ymax=courant->y;
		courant=courant->suivant;
	}
}

unsigned int liste_pixels::nb_cc (liste_pixels* L, unsigned short int icc, bool iaf) {
	int i,j,k,i_cc,j_cc,n,ncc=0;
	rectangle_englobant ();
	int nl=xmax-xmin+1, nc=ymax-ymin+1, i_cc0=nl, j_cc0=nc;
//	if (iaf) cout<<" xmin= "<<xmin<<" xmax= "<<xmax<<" ymin= "<<ymin<<" ymax= "<<ymax<<"\n";
	imabin imaL(nl,nc);
	courant=debut;
	while (courant!=NULL) {
		imaL(courant->x-xmin,courant->y-ymin)=1;
		courant=courant->suivant;
	}
	imadata<int> imaCC_L=imaL.composantes_connexes(ncc,icc,iaf);
	if (iaf) imaCC_L.affiche();
	if (L->nb_elts()>0) L->~liste_pixels();
	if (ncc>=1) {
		for (k=1; k<=ncc; k++) {
			i_cc=nl; j_cc=nc; n=i_cc*nc+j_cc;
			for (i=0; i<nl; i++)
				for (j=0; j<nc; j++)
					if (imaCC_L(i,j)==k) {
						if (i*nc+j<n) {i_cc=i; j_cc=j; n=i_cc*nc+j_cc;}
					}
			L->insere(i_cc+xmin,j_cc+ymin);
			if (n<i_cc0*nc+j_cc0) {i_cc0=i_cc; j_cc0=j_cc;}
		}
		L->supprime(i_cc0+xmin,j_cc0+ymin);
		L->insere(i_cc0+xmin,j_cc0+ymin); // de façon a avoir la composante connexe 'exterieure' en position 0
	}
	if (iaf) {
		cout<<" # de composantes connexes de N = "<<ncc<<" = "<<L->nb_elts()<<"\n";
		cout<<" affichage de la liste des pixels 1 representant par composante connexe de N :\n";
		L->affiche();
	}
	return (unsigned int)ncc;
}
