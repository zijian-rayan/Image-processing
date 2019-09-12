#include "fichiers.h"

void num_2_char_blabla (int nb, char *no, int n) {
	int i=around(pow((double)10,n-1.)),j;
	for (j=0; j<n; j++) {
		*(no+j)=0X30+((nb/i)%10);
		i/=10;
	}
}

void fichimage_sortie::ecrit_ImaLab (imalabels &ima, const int k) {
	int i,j;
	BYTE val;
	nlig=ima.nlig(); ncol=ima.ncol();
	offsetlig=ima.ncol();
	sortie.write((char*)&offsetlig,sizeof(int)); // entete SLH + acces indirect fortran !!!
	sortie.write((char*)&nlig,sizeof(int));
	sortie.write((char*)&ncol,sizeof(int));
	itype=0;
	sizeval=sizeof(BYTE);
	sortie.write((char*)&i,sizeof(int));
	sortie.write((char*)&offsetlig,sizeof(int));
	for (i=0; i<nlig; i++) {
		sortie.write((char*)&offsetlig,sizeof(int));
		for (j=0; j<ncol; j++) {
			val=ima(i,j,k);
			sortie.write((char*)&val,sizeval);
		}
		sortie.write((char*)&offsetlig,sizeof(int));
	}
}

void fichimage_sortie::ecrit_ImaRGB (imadata<BYTE> &ima) {
	int i,j,k;
	BYTE val;
	nlig=ima.nlig(); ncol=ima.ncol(); itype=0; 
	sizeval=sizeof(BYTE); offsetlig=ima.ncol(); 
	for (i=0; i<nlig; i++)
		for (j=0; j<ncol; j++)
			for (k=0; k<3; k++) {
				val=ima(i,j,k);
				sortie.write((char*)&val,sizeval);
			}
}

void fichimage_entree::changefichier(char *nomfile) {
	entree.close();
	entree.open(nomfile,ios::in|ios::binary);
	if (!entree) {
		cout<<" ouverture de "<<nomfile<<" impossible\n";
		exit (-1);
	}
	nomfich=nomfile;
}

int fichimage_entree::lit_enteteIma (char *formatfile) {
	offsetfich=0;
	offsetlig=0;
	if (strcmp(formatfile,"svm")==0) {
		entree.seekg(sizeof(int),ios::beg);
		entree.read((char*)&nlig,sizeof(int));
		entree.read((char*)&ncol,sizeof(int));
		entree.read((char*)&itype,sizeof(int));
		cout<<" fichier de "<<nlig<<" lig. et "<<ncol;
		cout<<" col. de donnees de type "<<itype<<"\n";
		offsetfich=5*sizeof(int);
		offsetlig=sizeof(int);
	}
	else
		if (strcmp(formatfile,"dat")==0) {
			entree.read((char*)&nlig,sizeof(int));
			entree.read((char*)&ncol,sizeof(int));
			entree.read((char*)&itype,sizeof(int));
			cout<<" fichier de "<<nlig<<" lig. et "<<ncol;
			cout<<" col. de donnees de type "<<itype<<"\n";
			offsetfich=3*sizeof(int);
		}
		else {
			cout<<" donner les # de lignes et de colonnes et le type de l'image : ";
			cin>>nlig>>ncol>>itype;
		}
	init_enteteIma (nlig,ncol,itype,1,offsetfich);
/*	switch (itype) {
		case 0: sizeval=sizeof(BYTE);
			break;
		case 1: sizeval=sizeof(short int);
			break;
		case 2: sizeval=sizeof(long int);
			break;
		case 3: sizeval=sizeof(float);
			break;
		case 4: sizeval=sizeof(double);
			break;
		default: cout<<" type de donnees inconnu\n";
			break;
	}
	sizelig=sizeval*ncol+2*offsetlig;
	sizecanal=sizelig*nlig;
	cout<<" image(s) de "<<nlig<<" lignes, "<<ncol<<" colonnes de type "<<itype<<"\n";
	cout<<" taille en-tete fichier = "<<offsetfich<<", taille entete ligne = "<<offsetlig<<"\n";
	cout<<" taille enregistrement = "<<sizelig<<", taille canal = "<<sizecanal<<"\n";*/
	return itype;
}

void fichimage_entree::ini_enteteIma (const int _nlig, const int _ncol, const int _itype, const int nbip, const int _offsetfich) {
	nlig=_nlig; ncol=_ncol; itype=_itype; offsetfich=_offsetfich;
	offsetlig=0;
	switch (itype) {
		case 0: sizeval=sizeof(BYTE); break;
		case 1: sizeval=sizeof(short int); break;
		case 2: sizeval=sizeof(long int); break;
		case 3: sizeval=sizeof(float); break;
		case 4: sizeval=sizeof(double); break;
		default: cout<<" type de donnees inconnu\n"; sizeval=0; break;
	}
	sizelig=sizeval*nbip*ncol+2*offsetlig;
	sizecanal=sizelig*nlig;
	cout<<" fichier-image de "<<nlig<<" lig. x "<<ncol<<" col. de type "<<itype<<" code sur "<<sizeval<<" octets\n";
	cout<<" offsetfich = "<<offsetfich<<" sizecanal = "<<sizecanal<<" sizelig = "<<sizelig<<"\n";
}

void fichimage_entree::init_enteteIma (const int _nlig, const int _ncol, const int _itype, const int nbip, const int _offsetfich, ostream &os) {
	nlig=_nlig; ncol=_ncol; itype=_itype; offsetfich=_offsetfich;
	offsetlig=0;
	switch (itype) {
		case 0: sizeval=sizeof(BYTE); break;
		case 1: sizeval=sizeof(unsigned __int16); break;
		case 2: sizeval=sizeof(long int); break;
		case 3: sizeval=sizeof(float); break;
		case 4: sizeval=sizeof(double); break;
		default: cout<<" type de donnees inconnu\n"; sizeval=0; break;
	}
	sizelig=sizeval*nbip*ncol+2*offsetlig;
	sizecanal=sizelig*nlig;
	os<<" image(s) de "<<nlig<<" lignes, "<<ncol<<" colonnes de type "<<itype<<"\n";
	os<<" taille en-tete fichier = "<<offsetfich<<", taille entete ligne = "<<offsetlig<<"\n";
	os<<" taille enregistrement = "<<sizelig<<", taille canal bsq = "<<sizecanal<<"\n";
}

int fichimage_entree::recopie_enteteIma (const fichimage_entree &f, ostream &os) {
	nlig=f.nlig;
	ncol=f.ncol;
	itype=f.itype;
	offsetfich=f.offsetfich;
	offsetlig=f.offsetlig;
	sizeval=f.sizeval;
	sizelig=f.sizelig;
	sizecanal=f.sizecanal;
	os<<" image(s) de "<<nlig<<" lignes, "<<ncol<<" colonnes de type "<<itype<<"\n";
	os<<" taille en-tete fichier = "<<offsetfich<<", taille entete ligne = "<<offsetlig<<"\n";
	os<<" taille enregistrement = "<<sizelig<<", taille canal = "<<sizecanal<<"\n";
	return itype;
}

void fichimage_entree::saut_debut_bsq (const int ilig0, const int icol0) {
	offsetfich=offsetfich+ilig0*sizelig;
	offsetlig=offsetlig+icol0*sizeval;
	nlig=nlig-ilig0;
	ncol=ncol-icol0;
}

imalabels fichimage_entree::lit_ImaLab (const int icanal) {
	int i,j;
	imalabels ima(nlig,ncol);
	BYTE val;
	for (i=0; i<nlig; i++) {
		entree.seekg(offsetfich+sizecanal*icanal+sizelig*i+offsetlig,ios::beg);
		for (j=0; j<ncol; j++) {
			entree.read((char*)&val,sizeval);
			ima(i,j,0)=val;
		}
	}
	return ima;
}

imabin fichimage_entree::lit_Imabinaire () {
	int i,j;
	imabin ima(nlig,ncol);
	BYTE val;
	for (i=0; i<nlig; i++) {
		entree.seekg(offsetfich+sizelig*i+offsetlig,ios::beg);
		for (j=0; j<ncol; j++) {
			entree.read((char*)&val,sizeval);
			ima(i,j)=!!val;
		}
	}
	return ima;
}

template <class T> imadata<T> fichimage_entree::lit_rawIma (T u, const int icanal) {
	int i,j;
	imadata<T> ima(nlig,ncol);
	T val;
	for (i=0; i<nlig; i++) {
		entree.seekg(offsetfich+sizecanal*icanal+sizelig*i+offsetlig,ios::beg);
		for (j=0; j<ncol; j++) {
			entree.read((char*)&val,sizeval);
			ima(i,j,0)=val;
		}
	}
	return ima;
}

imadata<unsigned __int16> fichimage_entree::lit_rawIma (const int icanal) {
	int i,j;
	imadata<unsigned __int16> ima(nlig,ncol);
	unsigned __int16 val;
	for (i=0; i<nlig; i++) {
		entree.seekg(offsetfich+sizecanal*icanal+sizelig*i+offsetlig,ios::beg);
		for (j=0; j<ncol; j++) {
			entree.read((char*)&val,sizeval);
			ima(i,j,0)=val;
		}
	}
	return ima;
}

template <class T> imadata<T> fichimage_entree::lit_NrawIma (T u, const int ncanaux) {
	int i,j,k;
	imadata<T> ima(nlig,ncol,ncanaux);
	T val;
	for (i=0; i<nlig; i++) {
		entree.seekg(offsetfich+sizelig*i+offsetlig,ios::beg);
		for (j=0; j<ncol; j++) {
			entree.read((char*)&val,sizeval);
			ima(i,j,0)=val;
		}
	}
	for (k=1; k<ncanaux; k++) {
		cout<<" repertoire+nom du fichier image #"<<k+1<<" en entree : ";
//		cin>>setw(LGMAX_NOM_FICH)>>nomfich;
		cin>>nomfich;
		ifstream ific;
		ific.open(nomfich,ios::in|ios::binary);
		if (!ific) {
			cout<<" ouverture de "<<nomfich<<" impossible\n";
			exit (-1);
		}
		for (i=0; i<nlig; i++) {
			ific.seekg(offsetfich+sizelig*i+offsetlig,ios::beg);
			for (j=0; j<ncol; j++) {
				ific.read((char*)&val,sizeval);
				ima(i,j,k)=val;
			}
		}
	}
	return ima;
}

imadata<BYTE> fichimage_entree::lit_NrawIma_ui (const int ncanaux) {
	int i,j,k;
	imadata<BYTE> ima(nlig,ncol,ncanaux);
	BYTE val;
	for (i=0; i<nlig; i++) {
		entree.seekg(offsetfich+sizelig*i+offsetlig,ios::beg);
		for (j=0; j<ncol; j++) {
			entree.read((char*)&val,sizeval);
			ima(i,j,0)=val;
		}
	}
	const int LGMAX_NOM_FICH=120;
	char nomfich[LGMAX_NOM_FICH+1];
	for (k=1; k<ncanaux; k++) {
		cout<<" repertoire+nom du fichier image #"<<k+1<<" en entree : ";
		cin>>setw(LGMAX_NOM_FICH)>>nomfich;
		ifstream ific;
		ific.open(nomfich,ios::in|ios::binary);
		if (!ific) {
			cout<<" ouverture de "<<nomfich<<" impossible\n";
			exit (-1);
		}
		for (i=0; i<nlig; i++) {
			ific.seekg(offsetfich+sizelig*i+offsetlig,ios::beg);
			for (j=0; j<ncol; j++) {
				ific.read((char*)&val,sizeval);
				ima(i,j,k)=val;
			}
		}
	}
	return ima;
}

imadata<BYTE> fichimage_entree::lit_1rawIma_rgb (const int nlig, const int ncol, const int ncan) {
	imadata<BYTE> ima(nlig,ncol,ncan);
	int i,j,k;
	BYTE val;
	for (i=0; i<nlig; i++) {
		entree.seekg(offsetfich+sizelig*i+offsetlig,ios::beg);
		for (j=0; j<ncol; j++)
			for (k=0; k<ncan; k++) {
				entree.read((char*)&val,sizeval);
				ima(i,j,k)=val;
			}
	}
	offsetfich+=sizelig*nlig;
	return ima;
}

imadata<BYTE> fichimage_entree::LoadPGM () { /* cette version ne lit que le type P5 */
	cout<<" chargement de l'image au format pgm P5 ou P6\n";
	char* buffer=new char[2];
	int i, j, k, ngr, ncanaux;
	offsetlig=0; offsetfich=0;
	entree.read(buffer,2*sizeof(char)); offsetfich+=2;
	char aa;
	entree.read(&aa,sizeof(char)); cout<<aa; offsetfich++;
	if (buffer[0]!='P' || (buffer[1]!='5' && buffer[1]!='6')) {
		cout<<" entete "<<buffer[0]<<buffer[1]<<" du fichier "<<nomfich<<" invalide\n";
	} else {
		cout<<" entete "<<buffer[0]<<buffer[1]<<" du fichier "<<nomfich<<" valide\n";
		ncanaux=1;
		if (buffer[1]=='6') ncanaux=3;
		bool fin=1;                                   // lecture des lignes de commentaires (commençant par #)
		entree.read(&aa,sizeof(char)); cout<<aa; offsetfich++;
		if (aa=='#') {
			fin=0;
			while (!fin) {
				entree.read(&aa,sizeof(char)); cout<<aa; offsetfich++;
				if (aa=='\n') {
					entree.read(&aa,sizeof(char)); cout<<aa; offsetfich++;
					if (aa!='#') fin=1;
				}
			}
		}
		cout<<" fin des lignes de commentaires\n";
		offsetfich-=1; entree.seekg(offsetfich,ios::beg); // on revient 1 car. en arrière pour lire nbre de col, lig, niv gris
		fin=0; ncol=0;
		while (!fin) {
			entree.read(&aa,sizeof(char)); offsetfich++;
			if (aa!='\n' && aa!=' ') ncol=ncol*10+atoi(&aa);
			else fin=1;
		}
		fin=0; nlig=0;
		while (!fin) {
			entree.read(&aa,sizeof(char)); offsetfich++;
			if (aa!='\n' && aa!=' ') nlig=nlig*10+atoi(&aa);
			else fin=1;
		}
		cout<<" fichier de "<<nlig<<" lig. et "<<ncol<<" col.\n";
		fin=0; ngr=0;
		while (!fin) {
			entree.read(&aa,sizeof(char)); offsetfich++;
			if (aa!='\n' && aa!=' ') ngr=ngr*10+atoi(&aa);
			else fin=1;
		}
		if (ngr==255) itype=0;
		if (ngr==65535) itype=1;
		cout<<" "<<ngr<<" niv. de gris => donnees de type "<<itype<<"\n";
		offsetfich*=sizeof(char);
		switch (itype) {
			case 0: sizeval=sizeof(BYTE); break;
			case 1: sizeval=sizeof(unsigned __int16); break;
			case 2: sizeval=sizeof(long int); break;
			case 3: sizeval=sizeof(float); break;
			case 4: sizeval=sizeof(double); break;
			default: cout<<" type de donnees inconnu\n"; break;
		}
		sizelig=sizeval*ncol+2*offsetlig;
		sizecanal=sizelig*nlig;
	}
	if (buffer!=NULL) delete[] buffer; 
	cout<<" nb canaux = "<<ncanaux<<"\n";
	imadata<BYTE> ima(nlig,ncol,ncanaux);
	if (itype==0) {
		BYTE val;
		for (i=0; i<nlig; i++) {
//			entree.seekg(offsetfich+sizelig*i+offsetlig,ios::beg);
			for (j=0; j<ncol; j++) {
				for (k=0; k<ncanaux; k++) {
					entree.read((char*)&val,sizeval);
					ima(i,j,k)=val;
				}
			}
		}
	}
	if (itype==1) {
		BYTE val1,val2;
		for (i=0; i<nlig; i++) {
			for (j=0; j<ncol; j++) {
				for (k=0; k<ncanaux; k++) {
					entree.read((char*)&val1,sizeof(BYTE));
					entree.read((char*)&val2,sizeof(BYTE));
					if (val1>0) cout<<" attention troncature des valeurs au dela de 255 !!!\n";
					ima(i,j,k)=(BYTE)val2;
				}
			}
		}
	}
	ima.statbasic(1);
	entree.seekg(0,ios::beg);
	return ima;
}

imadata<unsigned __int16> fichimage_entree::LoadPGM2 () { /* cette version ne lit que le type P5 */
	cout<<" chargement de l'image au format pgm P5 ou P6\n";
	char* buffer=new char[2];
	int i, j, k, ngr, ncanaux;
	offsetlig=0; offsetfich=0;
	entree.read(buffer,2*sizeof(char)); offsetfich+=2;
	char aa;
	entree.read(&aa,sizeof(char)); cout<<aa; offsetfich++;
	if (buffer[0]!='P' || (buffer[1]!='5' && buffer[1]!='6')) {
		cout<<" entete "<<buffer[0]<<buffer[1]<<" du fichier "<<nomfich<<" invalide\n";
	} else {
		cout<<" entete "<<buffer[0]<<buffer[1]<<" du fichier "<<nomfich<<" valide\n";
		ncanaux=1;
		if (buffer[1]=='6') ncanaux=3;
		bool fin=1;                                   // lecture des lignes de commentaires (commençant par #)
		entree.read(&aa,sizeof(char)); cout<<aa; offsetfich++;
		if (aa=='#') {
			fin=0;
			while (!fin) {
				entree.read(&aa,sizeof(char)); cout<<aa; offsetfich++;
				if (aa=='\n') {
					entree.read(&aa,sizeof(char)); cout<<aa; offsetfich++;
					if (aa!='#') fin=1;
				}
			}
		}
		cout<<" fin des lignes de commentaires\n";
		offsetfich-=1; entree.seekg(offsetfich,ios::beg); // on revient 1 car. en arrière pour lire nbre de col, lig, niv gris
		fin=0; ncol=0;
		while (!fin) {
			entree.read(&aa,sizeof(char)); offsetfich++;
			if (aa!='\n' && aa!=' ') ncol=ncol*10+atoi(&aa);
			else fin=1;
		}
		fin=0; nlig=0;
		while (!fin) {
			entree.read(&aa,sizeof(char)); offsetfich++;
			if (aa!='\n' && aa!=' ') nlig=nlig*10+atoi(&aa);
			else fin=1;
		}
		cout<<" fichier de "<<nlig<<" lig. et "<<ncol<<" col.\n";
		fin=0; ngr=0;
		while (!fin) {
			entree.read(&aa,sizeof(char)); offsetfich++;
			if (aa!='\n' && aa!=' ') ngr=ngr*10+atoi(&aa);
			else fin=1;
		}
		if (ngr==255) itype=0;
		if (ngr==65535) itype=1;
		cout<<" "<<ngr<<" niv. de gris => donnees de type "<<itype<<"\n";
		offsetfich*=sizeof(char);
		switch (itype) {
			case 0: sizeval=sizeof(BYTE); break;
			case 1: sizeval=sizeof(unsigned __int16); break;
			case 2: sizeval=sizeof(long int); break;
			case 3: sizeval=sizeof(float); break;
			case 4: sizeval=sizeof(double); break;
			default: cout<<" type de donnees inconnu\n"; break;
		}
		sizelig=sizeval*ncol+2*offsetlig;
		sizecanal=sizelig*nlig;
	}
	if (buffer!=NULL) delete[] buffer; 
	cout<<" nb canaux = "<<ncanaux<<"\n";
	imadata<unsigned __int16> ima(nlig,ncol,ncanaux);
	if (itype==0) {
		BYTE val;
		for (i=0; i<nlig; i++) {
//			entree.seekg(offsetfich+sizelig*i+offsetlig,ios::beg);
			for (j=0; j<ncol; j++) {
				for (k=0; k<ncanaux; k++) {
					entree.read((char*)&val,sizeval);
					ima(i,j,k)=(unsigned __int16)val;
				}
			}
		}
	}
	if (itype==1) {
		ima=imadata<unsigned __int16>(nlig,ncol,ncanaux);
		BYTE val1,val2;
		for (i=0; i<nlig; i++) {
			for (j=0; j<ncol; j++) {
				for (k=0; k<ncanaux; k++) {
					entree.read((char*)&val1,sizeof(BYTE));
					entree.read((char*)&val2,sizeof(BYTE));
					ima(i,j,k)=(unsigned __int16)val1*256+val2;
				}
			}
		}
	}
	ima.statbasic(1);
	entree.seekg(0,ios::beg);
	return ima;
}

imadata<int32> fichimage_entree::lit_NrawIma_i32 (const int ncanaux) {
	int i,j,k;
	imadata<int32> ima(nlig,ncol,ncanaux);
	int32 val;
	for (i=0; i<nlig; i++) {
		entree.seekg(offsetfich+sizelig*i+offsetlig,ios::beg);
		for (j=0; j<ncol; j++) {
			entree.read((char*)&val,sizeval);
			ima(i,j,0)=val;
		}
	}
	const int LGMAX_NOM_FICH=120;
	char nomfich[LGMAX_NOM_FICH+1];
	for (k=1; k<ncanaux; k++) {
		cout<<" repertoire+nom du fichier image #"<<k+1<<" en entree : ";
		cin>>setw(LGMAX_NOM_FICH)>>nomfich;
		ifstream ific;
		ific.open(nomfich,ios::in|ios::binary);
		if (!ific) {
			cout<<" ouverture de "<<nomfich<<" impossible\n";
			exit (-1);
		}
		for (i=0; i<nlig; i++) {
			ific.seekg(offsetfich+sizelig*i+offsetlig,ios::beg);
			for (j=0; j<ncol; j++) {
				ific.read((char*)&val,sizeval);
				ima(i,j,k)=val;
			}
		}
	}
	return ima;
}

imadata<float> fichimage_entree::lit_NrawIma_fl (const int ncanaux) {
	int i,j,k;
	imadata<float> ima(nlig,ncol,ncanaux);
	float val;
	for (i=0; i<nlig; i++) {
		entree.seekg(offsetfich+sizelig*i+offsetlig,ios::beg);
		for (j=0; j<ncol; j++) {
			entree.read((char*)&val,sizeval);
			ima(i,j,0)=val;
			}
		}
	const int LGMAX_NOM_FICH=120;
	char nomfich[LGMAX_NOM_FICH+1];
	for (k=1; k<ncanaux; k++) {
		cout<<" repertoire+nom du fichier image #"<<k+1<<" en entree : ";
		cin>>setw(LGMAX_NOM_FICH)>>nomfich;
//		cin>>nomfich;
		ifstream ific;
		ific.open(nomfich,ios::in|ios::binary);
		if (!ific) {
			cout<<" ouverture de "<<nomfich<<" impossible\n";
			exit (-1);
			}
		for (i=0; i<nlig; i++) {
			ific.seekg(offsetfich+sizelig*i+offsetlig,ios::beg);
			for (j=0; j<ncol; j++) {
				ific.read((char*)&val,sizeval);
				ima(i,j,k)=val;
				}
			}
		}
	return ima;
	}

imadata<float> fichimage_entree::lit_ImaBSQ_fl (const int ncanaux) {
	int i,j,k;
	imadata<float> ima(nlig,ncol,ncanaux);
	cout<<" lecture d'1 image "<<nlig<<" lig. x "<<ncol<<" col. x "<<ncanaux<<" canaux\n";
	cout<<" offsetfich = "<<offsetfich<<" sizecanal = "<<sizecanal<<" sizelig = "<<sizelig<<"\n";
	float val;
	for (k=0; k<ncanaux; k++)
		for (i=0; i<nlig; i++) {
			entree.seekg(offsetfich+sizecanal*k+sizelig*i+offsetlig,ios::beg);
			for (j=0; j<ncol; j++) {
				entree.read((char*)&val,sizeval);
				ima(i,j,k)=val;
			}
		}
	return ima;
	}

template <class T> imadata<T> fichimage_entree::lit_bsqIma (T u, const int ncanaux) {
	int i,j,k;
	imadata<T> ima(nlig,ncol,ncanaux);
	T val;
	for (k=0; k<ncanaux; k++)
		for (i=0; i<nlig; i++) {
			entree.seekg(offsetfich+sizecanal*k+sizelig*i+offsetlig,ios::beg);
			for (j=0; j<ncol; j++) {
				entree.read((char*)&val,sizeval);
				ima(i,j,k)=val;
				}
			}
	return ima;
	}

template <class T> imadata<T> fichimage_entree::lit_bilIma (T u, const int ncanaux) {
   int i,j,k;
   imadata<T> ima(nlig,ncol,ncanaux);
   T val;
   sizelig=(sizeval*ncol+2*offsetlig)*ncanaux;
   for (i=0; i<nlig; i++) {
      for (k=0; k<ncanaux; k++) {
         entree.seekg(offsetfich+sizelig*i+offsetlig,ios::beg);
         for (j=0; j<ncol; j++) {
             entree.read((char*)&val,sizeval);
             ima(i,j,k)=val;
            }
         }
      }
   return ima;
   }

imadata<BYTE> fichimage_entree::lit_bipIma_ui (const int ncanaux) {
   int i,j,k;
   imadata<BYTE> ima(nlig,ncol,ncanaux);
   BYTE val;
   sizelig=sizeval*ncanaux*ncol+2*offsetlig;
   for (i=0; i<nlig; i++) {
      entree.seekg(offsetfich+sizelig*i+offsetlig,ios::beg);
      for (j=0; j<ncol; j++) {
         for (k=0; k<ncanaux; k++) {
            entree.read((char*)&val,sizeval);
            ima(i,j,k)=val;
            }
         }
      }
   return ima;
   }

