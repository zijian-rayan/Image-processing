#ifndef _FICHIERS_H
#define _FICHIERS_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
using namespace std;
#include "def.h"
#include "fonctions.h"
#include "pixel.h"
#include "imasites.h"
#include "imadata.h"
#include "imabin.h"
#include "statclass.h"
#include "imalabels.h"
#include "imacontours.h"

// conversion chaine caracteres
/*void num_2_char_blabla (int nb, char *no, int n) {
	int i=around(pow((double)10,n-1.)),j;
	for (j=0; j<n; j++) {
		*(no+j)=0X30+((nb/i)%10); 
		i/=10;
	}
}*/
void num_2_char_blabla (int nb, char *no, int n);

template <class T> class imadata;
class imalabels;
class imabin;

class fichimage_sortie
{	ofstream sortie;
	char *nomfich;
	int nlig, ncol, itype, sizeval, offsetlig;
	bool iO;
//	int offsetfich, sizelig, sizecanal;
 public:
	fichimage_sortie(string nomfile, bool itest=0) {
		iO=0;
		sortie.open(nomfile.c_str(),ios::out|ios::binary);
		if (!sortie) {
			if (!itest) {cout<<" ouverture de "<<nomfile<<" impossible\n"; exit (-1);}
		}
		else iO=1;
		nomfich=(char*)nomfile.c_str();
	}
	fichimage_sortie(char *nomfile) {
		sortie.open(nomfile,ios::out|ios::binary);
		if (!sortie) {
			cout<<" ouverture de "<<nomfile<<" impossible\n";
			exit (-1);
		}
		nomfich=nomfile;
	}
	~fichimage_sortie() {
		if (sortie) sortie.close();
	}
//    void ecrit_enteteIma (int,char * formatfile="svm");
    void ecrit_ImaLab (imalabels&,const int=0);
    void ecrit_ImaRGB (imadata<BYTE>&);
//    template <class T> imadata<T> ecrit_NrawIma (T,int=1);
//    template <class T> imadata<T> ecrit_bsqIma (T,int=1);
//    template <class T> imadata<T> lit_bilIma (T,int=1);       //non verifie
//    template <class T> imadata<T> lit_bipIma (T,int=1);       //non verifie
};

class fichimage_entree
{	ifstream entree;
	char *nomfich;
	int nlig, ncol, itype, sizeval, offsetfich, offsetlig, sizelig, sizecanal;
	bool iO;
 public:
	fichimage_entree(string nomfile, bool itest=0) {//cout<<" ouverture de "<<nomfile<<"\n";
		iO=0;
		entree.open(nomfile.c_str(),ios::in|ios::binary);
		if (!entree) {
			if (!itest) {cout<<" ouverture de "<<nomfile<<" impossible\n"; exit (-1);}
		}
		else {iO=1; cout<<" ouverture de "<<nomfile<<" reussie\n";}
		nomfich=new char[nomfile.length()+1];
		strcpy (nomfich, nomfile.c_str()); //nomfich=(char*)nomfile.c_str(); 
//		cout<<" nomfile = "<<nomfile.c_str()<<" "<<nomfile.length()+1<<" -> nomfich = "<<nomfich<<"\n";
	}
	fichimage_entree(char *nomfile) {
		entree.open(nomfile,ios::in|ios::binary);
		if (!entree) {
			cout<<" ouverture de "<<nomfile<<" impossible\n";
			exit (-1);
		}
		nomfich=new char[strlen(nomfile)+1];
		strcpy(nomfich,nomfile);
	}
	~fichimage_entree() {if (entree) entree.close(); if (nomfich!=NULL) delete[] nomfich;}

	fichimage_entree (const fichimage_entree &f) {
		cout<<" constructeur copie\n";
		nomfich=new char[strlen(f.nomfich)+1];
		strcpy(nomfich,f.nomfich);
		nlig=f.nlig; ncol=f.ncol; itype=f.itype; sizeval=f.sizeval; offsetfich=f.offsetfich;
		offsetlig=f.offsetlig; sizelig=f.sizelig; sizecanal=f.sizecanal; 
//		f.entree.close();
		entree.open(nomfich,ios::in|ios::binary);
		if (!entree) {
			cout<<" ouverture de "<<nomfich<<" impossible"<<entree<<"\n";
			exit (-1);
		}
	}

	fichimage_entree& operator=(const fichimage_entree &f) {
		cout<<" operateur affectation "<<f.nomfich<<"\n";
		if (this != &f) {
			if (nomfich!=NULL) delete[] nomfich;
			nomfich=new char[strlen(f.nomfich)+1];
			strcpy(nomfich,f.nomfich);
			nlig=f.nlig; ncol=f.ncol; itype=f.itype; sizeval=f.sizeval; offsetfich=f.offsetfich;
			offsetlig=f.offsetlig; sizelig=f.sizelig; sizecanal=f.sizecanal; 
//			f.entree.close();
			entree.close();
			entree.open(nomfich,ios::in|ios::binary);
			if (!entree) {
				cout<<" ouverture de "<<nomfich<<" impossible "<<entree<<"\n";
				exit (-1);
			}
		}
		return *this;
	}

	int type () const {return itype;}
	bool ok () const {return iO;}
	void changefichier(char*);
	int lit_enteteIma (char* ="svm");
	void init_enteteIma (const int, const int, const int, const int=1, const int=0, ostream& =cout);
	void ini_enteteIma (const int, const int, const int, const int=1, const int=0);
	int recopie_enteteIma (const fichimage_entree&, ostream& =cout);
	void saut_debut_bsq (const int=0, const int=0);
	void goto_imaN (const int n) {offsetfich=sizecanal*n; /*cout<<" offsetfich = "<<offsetfich<<"\n";*/}
	imalabels lit_ImaLab (const int=0);
	template <class T> imadata<T> lit_rawIma (T,const int=0);
	imadata<unsigned __int16> lit_rawIma(const int=0);
	template <class T> imadata<T> lit_NrawIma (T,const int=1);
	imadata<BYTE> lit_NrawIma_ui (const int=1);
	imadata<BYTE> lit_1rawIma_rgb (const int, const int, const int=3);
	imadata<BYTE> LoadPGM ();
	imadata<unsigned __int16> LoadPGM2 ();
	imadata<int32> lit_NrawIma_i32 (const int=1);
	imadata<float> lit_NrawIma_fl (const int=1);
	imadata<float> lit_ImaBSQ_fl (const int=1);
	template <class T> imadata<T> lit_bsqIma (T,const int=1);
	template <class T> imadata<T> lit_bilIma (T,const int=1);       //non verifie
	imadata<BYTE> lit_bipIma_ui (const int=1);
	imabin lit_Imabinaire ();
};

#endif