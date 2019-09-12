#ifndef _TREE_H
#define _TREE_H

#include <iostream>
using namespace std;

#include "liste_pixels.h"

class elt_tree {
public:
	int ident;
	float val;
	liste_pixels L, L_ini_hole;
//	int nb_freres;
	elt_tree *fils, *frere, *pere;
	BYTE etage;
//	bool valid;
public:
	elt_tree () {
		ident=0; val=0.f; etage=0;
		fils=frere=pere=NULL;
	}
	~elt_tree () {
//		cout<<" destructeur de elt_tree @ "<<this<<"\n";
//		if (fils!=NULL) delete fils;
//		if (frere!=NULL) delete frere;
//		if (pere!=NULL) delete pere;
	}
	elt_tree (const elt_tree &E) {
		ident=E.ident; val=E.val; etage=E.etage; L=E.L; 
		fils=frere=pere=NULL; L_ini_hole=E.L_ini_hole;
	}
	elt_tree& operator= (const elt_tree &E) {
		if (this!=&E) {
			if (&E!=NULL) {
				ident=E.ident; val=E.val; etage=E.etage; L=E.L; L_ini_hole=E.L_ini_hole;}
			else {
				ident=0; val=0.f; etage=0; L=L_ini_hole=liste_pixels();}
			fils=frere=pere=NULL; 
		}
		return *this;
	}
};

class ad_elt_tree {
public:
	BYTE type_origin;
	elt_tree *ad_origine, *adres_elt;
public:
	ad_elt_tree () {
		type_origin=0;
		ad_origine=adres_elt=NULL;
	}
	~ad_elt_tree () {
//		if (ad_origine!=NULL) delete ad_origine;
//		if (adres_elt!=NULL) delete adres_elt;
	}
};

class liste_ad_elt_tree {
public:
	int n_adres;
	ad_elt_tree *L;
public:
	liste_ad_elt_tree (int n=0) {
		n_adres=0; L=NULL; 
		if (n>0) L=new ad_elt_tree[n];
	}
	~liste_ad_elt_tree () {
		if (L!=NULL) delete[] L;
	}
	void add_elt2liste (elt_tree *courant, elt_tree *entreeT) { //cout<<" on rajoute les fils et freres de "<<courant<<" avec comme entree "<<entreeT<<"\n";
		if (courant->frere!=NULL) {
			L[n_adres].ad_origine=entreeT;
			L[n_adres].adres_elt=courant->frere;
			L[n_adres].type_origin=2;
			n_adres++;
		}
		if (courant->fils!=NULL) {
			L[n_adres].ad_origine=entreeT;
			L[n_adres].adres_elt=courant->fils;
			L[n_adres].type_origin=1;
			n_adres++;
		}
	}
	BYTE sub_elt2liste (elt_tree* &courant, elt_tree* &entreeT) { 
		n_adres--;
		courant=L[n_adres].adres_elt;
		entreeT=L[n_adres].ad_origine; //cout<<" on reccupere "<<courant<<" avec comme entree "<<entreeT<<"\n";
		return L[n_adres].type_origin;
	}
};

class tree {
	int nb_etages, nb_elts;
	elt_tree *racine, *courant, *precedent;
public:
	tree () {
		racine=courant=precedent=NULL;
		nb_etages=nb_elts=0;
	}
	~tree () {
		int n_adres;
		if (nb_elts>0) {
			elt_tree **liste_adres=new elt_tree*[nb_elts];
			liste_adres[0]=racine;
			courant=liste_adres[0];
			n_adres=1;
			while (n_adres>0) {
				n_adres--;
				courant=liste_adres[n_adres];
				if (courant->frere!=NULL) liste_adres[n_adres++]=courant->frere;
				if (courant->fils!=NULL) liste_adres[n_adres++]=courant->fils;
//				cout<<" destruction de elt_tree a l'adresse "<<courant<<"\n";
				delete courant; 
				nb_elts--; //cout<<" reste "<<nb_elts<<"\n";
			}
			if (liste_adres!=NULL) delete[] liste_adres;
		}
	}
	tree (const tree &T) {
		int i, n_adres;
		elt_tree **liste_adres=new elt_tree*[nb_elts];
		for (i=0; i<nb_elts; i++) liste_adres[i]=NULL;
		n_adres=0;
		liste_adres[n_adres++]=racine;
		courant=liste_adres[n_adres-1];
		while (courant!=NULL) {
			n_adres--;
			if (courant->frere!=NULL) liste_adres[n_adres++]=courant->frere;
			if (courant->fils!=NULL) liste_adres[n_adres++]=courant->fils;
			delete courant;
			nb_elts--;
			courant=liste_adres[n_adres-1];
		}
		if (liste_adres!=NULL) delete[] liste_adres;
		liste_ad_elt_tree liste_elt_tree(T.nb_elts);
		courant=T.racine;
		if (courant==NULL) {
			racine=courant=NULL;
			nb_etages=nb_elts=0;
		}
		else {
			elt_tree *nouveau=new elt_tree; *nouveau=*courant;
			racine=nouveau;
			nb_etages=nb_elts=1;
			liste_elt_tree.add_elt2liste (courant,nouveau);	
			while (liste_elt_tree.n_adres>0) {
				BYTE itype=liste_elt_tree.sub_elt2liste (courant,precedent);
				elt_tree *nouveau2=new elt_tree; *nouveau2=*courant;
				if (itype==2 && precedent->frere==NULL) {
					precedent->frere=nouveau2;
					nouveau2->pere=precedent->pere;
				} else
					if (precedent->frere!=NULL) cout<<" place du frere deja occupee !!!\n";
				if (itype==1 && precedent->fils==NULL) {
					precedent->fils=nouveau2;
					nouveau2->pere=precedent;
				} else
					if (precedent->fils!=NULL) cout<<" place du fils deja occupee !!!\n";
				nb_elts++;
				liste_elt_tree.add_elt2liste (courant,nouveau2);	
			}
		}
		nb_etages=T.nb_etages;
		if (nb_elts!=T.nb_elts) cout<<" probleme dans recopie arbre\n";
		courant=racine;
	}
	tree& operator= (const tree& T) { //cout<<" operateur d'affectation de la classe tree\n";
		if (this!=&T) { //cout<<" les arbres sont effectivement differents donc on ecrase le 1er : "<<nb_elts<<" elements a detruire\n";
			if (nb_elts>0) {
				int i, n_adres;
				elt_tree **liste_adres=new elt_tree*[nb_elts];
				for (i=0; i<nb_elts; i++) liste_adres[i]=NULL;
				n_adres=0;
				liste_adres[n_adres++]=racine;
				courant=liste_adres[n_adres-1];
				while (courant!=NULL) {
					n_adres--;
					if (courant->frere!=NULL) liste_adres[n_adres++]=courant->frere; //cout<<" on stocke adresse du frere "<<courant->frere<<"\n";}
					if (courant->fils!=NULL) liste_adres[n_adres++]=courant->fils; //cout<<" on stocke adresse du fils "<<courant->frere<<"\n";}
					delete courant; 
					nb_elts--;
					courant=liste_adres[n_adres-1];
				}
				if (liste_adres!=NULL) delete[] liste_adres;
			}
			liste_ad_elt_tree liste_elt_tree(T.nb_elts);
			courant=T.racine;
			if (courant==NULL) {
				racine=courant=NULL;
				nb_etages=nb_elts=0;
			}
			else {
				elt_tree* nouveau=new elt_tree; *nouveau=*courant;
				racine=nouveau;
				nb_etages=nb_elts=1;
				liste_elt_tree.add_elt2liste (courant,racine); //affiche(courant);	
				while (liste_elt_tree.n_adres>0) {
					BYTE itype=liste_elt_tree.sub_elt2liste (courant,precedent); 
//					cout<<" @ courant recupere = "<<courant<<" @ precedent recupere = "<<precedent<<"\n";affiche(courant); affiche(precedent);
					elt_tree *nouveau2=new elt_tree; *nouveau2=*courant; 
					if (itype==2 && precedent->frere==NULL) {
						precedent->frere=nouveau2;
						nouveau2->pere=precedent->pere;
					} else
						if (precedent->frere!=NULL) cout<<" place du frere deja occupee !!!\n";
					if (itype==1 && precedent->fils==NULL) {
						precedent->fils=nouveau2;
						nouveau2->pere=precedent;
					}
					else
						if (precedent->fils!=NULL) cout<<" place du fils deja occupee !!!\n";
					nb_elts++;
					liste_elt_tree.add_elt2liste (courant,nouveau2);	
				}
			}
			nb_etages=T.nb_etages;
			if (nb_elts!=T.nb_elts) cout<<" probleme dans affectation arbre\n";
			courant=racine;
		}
		return *this;
	}
	int nbre_elts () {
		return nb_elts;
	}
	int nbre_etages () {
		return nb_etages;
	}
	elt_tree* cherche (int no, bool iaf=0) {
		bool trouve=0;
		elt_tree *ad_no=NULL, *ad_courant; ad_courant=courant;
		if (nb_elts>0) {
			int i, n_adres;
			elt_tree **liste_adres=new elt_tree*[nb_elts];
			for (i=0; i<nb_elts; i++) liste_adres[i]=NULL;
			n_adres=0;
			liste_adres[n_adres++]=racine;
			courant=liste_adres[n_adres-1];
			while (!trouve && courant!=NULL) {
				n_adres--;
				if (courant->frere!=NULL) liste_adres[n_adres++]=courant->frere;
				if (courant->fils!=NULL) liste_adres[n_adres++]=courant->fils;
				if (courant->ident==no) {
					trouve=1;
					ad_no=courant; if (iaf) cout<<" courant @"<<courant<<" a le bon ident "<<courant->ident<<"="<<no<<"\n";
				}
				else {
					if (iaf) cout<<" courant @"<<courant<<" n'a pas le bon ident "<<courant->ident<<"\n";
					if (n_adres>0) courant=liste_adres[n_adres-1]; 
					else courant=NULL;
				}
			}
			if (liste_adres!=NULL) delete[] liste_adres;
			courant=ad_courant;
		}
		return ad_no;
	}
	void insere_pere (elt_tree* E, const int no_lien, bool iaf=0) {
		if (iaf) {cout<<" element a inserer : \n"; affiche(E); cout<<" en tant que pere de element d'ident "<<no_lien<<"\n";}
		elt_tree *nouveau=new elt_tree; *nouveau=*E;
		courant=cherche(no_lien); if (iaf) cout<<" insere a l'@ "<<nouveau<<" en tant que pere de : "<<courant<<"\n";
		if (courant!=NULL) {
			nouveau->fils=courant;
			courant->pere=nouveau;
			courant=courant->frere;
			while (courant!=NULL) {
				courant->pere=nouveau;
				courant=courant->frere;
			}
		}
		nb_elts++;
		racine=nouveau; 
	}
	void insere (elt_tree* E, const int no_lien, const BYTE type_lien, bool iaf=0) {
		if (iaf) {cout<<" element a inserer : \n"; affiche(E);}
		elt_tree *nouveau=new elt_tree; *nouveau=*E;
		precedent=cherche(no_lien); if (iaf) cout<<" insere a l'@ "<<nouveau<<" a partir de : "<<precedent<<"\n";
		if (precedent==NULL) {
			racine=nouveau; if (iaf) cout<<" insertion a la racine\n";
		}
		else {
			if (type_lien==2) {
				if (precedent->frere!=NULL) {
					cout<<" place du frere deja occupee !!!\n";
					while (precedent->frere!=NULL) {precedent=precedent->frere;}
				}
				if (precedent->frere==NULL) {
					precedent->frere=nouveau; if (iaf) cout<<" insertion en tant que frere de "<<precedent<<"\n";
					nouveau->pere=precedent->pere; //if (iaf) affiche(nouveau);
				} else cout<<" echec dans la recherche du frere benjamin\n";	
			}
			if (type_lien==1) {
				if (precedent->fils==NULL) {
					precedent->fils=nouveau; if (iaf) cout<<" insertion en tant que fils de "<<precedent<<"\n";
					nouveau->pere=precedent; //if (iaf) affiche(nouveau);
				} else cout<<" place du fils deja occupee !!!\n";
			}
		}
		nb_elts++;
		if (E->frere!=NULL || E->fils!=NULL) {if (iaf) cout<<" insertion des fils et freres de "<<E<<"\n";
			const int nmax_elt=100;
			liste_ad_elt_tree liste_elt_tree(nmax_elt);
			liste_elt_tree.add_elt2liste(E,nouveau);
			while (liste_elt_tree.n_adres>0) {
				BYTE itype=liste_elt_tree.sub_elt2liste(courant,precedent);
				elt_tree *nouveau2=new elt_tree; *nouveau2=*courant;
				if (itype==2) {
					if (precedent->frere==NULL) {
						precedent->frere=nouveau2;
						nouveau2->pere=precedent->pere;
					} else cout<<" place du frere deja occupee !!!\n";
				}
				if (itype==1) {
					if (precedent->fils==NULL) {
						precedent->fils=nouveau2;
						nouveau2->pere=precedent;
					} else cout<<" place du fils deja occupee !!!\n";
				}
				nb_elts++;
//				liste_elt_tree.add_elt2liste(courant,nouveau2);
				if (liste_elt_tree.n_adres>=nmax_elt) cout<<" pb : nombre d'elements du sous-arbre a inserer trop grand_n";
				else liste_elt_tree.add_elt2liste (courant,nouveau2);
			}
		}
		else 
			if (iaf) cout<<" pas d'insertion de fils ou freres pour "<<E<<"\n";
	}
	void bouche_trou (elt_tree* E, const int no_lien, const BYTE type_lien, int x, int y) {
		insere(E,no_lien,type_lien);
		cherche(no_lien)->L_ini_hole.supprime (x,y);
	}
	imadata<float> ima_shape (int nlig, int ncol) { cout<<"\n\n calcul de l'image des shapes\n";
		imadata<float> I_Shape(nlig,ncol,3);
		int i,n;
		float v;
		nb_etages=0;
		liste_ad_elt_tree liste_adres(nb_elts);
		liste_adres.L[(liste_adres.n_adres)++].adres_elt=racine;
		courant=racine;
		courant->etage=0; nb_etages++; 
		elt_liste E; 
		liste_pixels L2;
		L2=(courant->L); 
		n=courant->etage; i=courant->ident; v=courant->val;
		while (L2.nb_elts()>0) {E=L2.extrait(); I_Shape(E.x,E.y,0)=(float)i; I_Shape(E.x,E.y,1)=(float)n; I_Shape(E.x,E.y,2)=v;}
		while (liste_adres.n_adres>0 && courant!=NULL) {
			(liste_adres.n_adres)--;
			if (courant->frere!=NULL) { cout<<" trouve frere\n";
				n=(liste_adres.n_adres)++;
				liste_adres.L[n].adres_elt=courant->frere;
				liste_adres.L[n].ad_origine=courant;
				liste_adres.L[n].type_origin=2;
			}
			if (courant->fils!=NULL) { cout<<" trouve fils\n";
				n=(liste_adres.n_adres)++;
				liste_adres.L[n].adres_elt=courant->fils;
				liste_adres.L[n].ad_origine=courant;
				liste_adres.L[n].type_origin=1;
			}
			n=liste_adres.n_adres-1; //cout<<" nb adresses = "<<n+1<<"\n";
			courant=liste_adres.L[n].adres_elt;
			if (courant!=NULL && n>=0) {
				courant->etage=(liste_adres.L[n].ad_origine)->etage+(2-liste_adres.L[n].type_origin);
				cout<<" ############## element \n"; affiche(courant); cout<<" a l'etage "<<(int)(courant->etage)<<"\n";
				n=(int)(courant->etage); i=courant->ident; v=courant->val; L2=(courant->L); 
				while (L2.nb_elts()>0) {
					E=L2.extrait(); 
					if (n>I_Shape(E.x,E.y,1)) {I_Shape(E.x,E.y,0)=(float)i; I_Shape(E.x,E.y,1)=(float)n; I_Shape(E.x,E.y,2)=v;}
				}
				if (n+1>nb_etages) nb_etages=n+1;
			}
		}
		courant=racine;
		return I_Shape;
	}
	elt_tree* search_hole () {
		elt_tree *elt_with_hole=NULL;
		int i, n_adres;
		elt_tree **liste_adres=new elt_tree*[nb_elts];
		for (i=0; i<nb_elts; i++) liste_adres[i]=NULL;
		n_adres=0;
		liste_adres[n_adres++]=racine;
		courant=liste_adres[n_adres-1];
		while (n_adres>0 && courant!=NULL && elt_with_hole==NULL) {
			n_adres--;
			if (courant->frere!=NULL) liste_adres[n_adres++]=courant->frere;
			if (courant->fils!=NULL) liste_adres[n_adres++]=courant->fils;
			if ((courant->L_ini_hole).nb_elts()>0) elt_with_hole=courant;
			courant=liste_adres[n_adres-1];
		}
		if (nb_elts>0 && liste_adres!=NULL) delete[] liste_adres;
		courant=racine;
		return elt_with_hole;
	}
	void affiche (elt_tree* aE) {
		cout<<" Element de l'arbre a l'adresse "<<aE<<"\n";
		cout<<"	 a pour identifiant "<<aE->ident<<"\n";
		cout<<"	 a pour valeur "<<aE->val<<"\n";
		cout<<"	 contient "<<(aE->L).nb_elts()<<" pixel(s) "; (aE->L).affiche(); cout<<"\n";
		if ((aE->L_ini_hole).nb_elts()>0) {
			cout<<"	 a "<<(aE->L_ini_hole).nb_elts()<<" trou(s) represente(s) par le(s) pixel(s) "; (aE->L_ini_hole).affiche(); cout<<"\n";
		} else cout<<"	 n'a pas de trou\n";
		if (aE->pere!=NULL) cout<<"	 a un pere a l'adresse "<<aE->pere;
		else cout<<"	 n'a pas de pere,";
		if (aE->fils!=NULL) cout<<" a un fils a l'adresse "<<aE->fils;
		else cout<<" n'a pas de fils,";
		if (aE->frere!=NULL) cout<<" et un frere a l'adresse "<<aE->frere<<"\n";
		else cout<<" et n'a pas de frere\n";
	}
	void affiche () {
		int i, n_adres;
		elt_tree **liste_adres=new elt_tree*[nb_elts];
		for (i=0; i<nb_elts; i++) liste_adres[i]=NULL;
		n_adres=0;
		liste_adres[n_adres++]=racine;
		courant=liste_adres[n_adres-1];
		while (n_adres>0 && courant!=NULL) {
			n_adres--;
			if (courant->frere!=NULL) liste_adres[n_adres++]=courant->frere;
			if (courant->fils!=NULL) liste_adres[n_adres++]=courant->fils;
			affiche(courant);
			courant=liste_adres[n_adres-1];
		}
		if (nb_elts>0 && liste_adres!=NULL) delete[] liste_adres;
		courant=racine;
	}

	elt_tree* cherche (int ilig, int icol) {
		bool trouve=0;
		elt_tree *ad_no=NULL, *ad_courant; ad_courant=courant;
		if (nb_elts>0) {
			int i, n_adres;
			elt_tree **liste_adres=new elt_tree*[nb_elts];
			for (i=0; i<nb_elts; i++) liste_adres[i]=NULL;
			n_adres=0;
			liste_adres[n_adres++]=racine; 
			courant=liste_adres[n_adres-1];
			while (!trouve && courant!=NULL) {
				n_adres--;
				cout<<" n_adres = "<<n_adres<<" courant "<<courant<<"\n";
				if (courant->frere!=NULL) liste_adres[n_adres++]=courant->frere;
				if (courant->fils!=NULL) liste_adres[n_adres++]=courant->fils;
				if ((courant->L).existe(ilig,icol)) {
					trouve=1;
					ad_no=courant; cout<<" courant @"<<courant<<" contient le pixel recherche dans sa liste : "; (courant->L).affiche(); cout<<"\n";
				}
				else {
					if (n_adres>0) courant=liste_adres[n_adres-1]; 
					else courant=NULL;
				}
			}
			if (liste_adres!=NULL) delete[] liste_adres;
			courant=ad_courant;
		}
		return ad_no;
	}

	void Arbre_Min (imadata<int>,bool=0); // calcul du Min_Tree
	void Arbre_Max (imadata<int>,bool=0); // calcul du Max_Tree
	int calc_etages (bool=0);
	imadata<int> image_Tree (int,int,bool=0);
};


#endif