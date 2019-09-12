#include <iostream>
using namespace std;

#include "menu.h"
#include "fichiers.h"
#include "sampleset.h"
#include "statclass.h"
#include "imaregions.h"

/* -------------------------------- */
int main_inv_ngr(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<4) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<"\n";
		exit(-1);
	}
    char *filename = argv[2];
	fichimage_entree ifich(filename);
	imadata<float> imadon=(imadata<float>)ifich.LoadPGM();

	for (int k=0; k<imadon.ncanaux(); k++) 
		imadon.inverse_ngr(255, k);
    
	imadon.sauve_ImaPGM(argv[3]);
    return(-1);
}

/* -------------------------------- */
int main_bin_sup(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<5) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<"\n";
		exit(-1);
	}
    char *filename = argv[3];
	fichimage_entree ifich(filename);
	imadata<float> imadon=(imadata<float>)ifich.LoadPGM();

	imabin imab(imadon, atof(argv[2]));
    
	imadata<BYTE>(imab).sauve_ImaPGM(argv[4]);
    return(-1);
}

/* -------------------------------- */
int main_bin_inf(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<5) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<"\n";
		exit(-1);
	}
    char *filename = argv[3];
	fichimage_entree ifich(filename);
	imadata<float> imadon=(imadata<float>)ifich.LoadPGM();

	const eps=0.01;
	imabin imab(imadon, atof(argv[2])+eps);
	imab=imab.negatif();
    
	imadata<BYTE>(imab).sauve_ImaPGM(argv[4]);
    return(-1);
}

/* -------------------------------- */
int main_seuil_hysteresis(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<6) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<"\n";
		exit(-1);
	}
    char *filename = argv[4];
	fichimage_entree ifich(filename);
	imadata<float> imadon=(imadata<float>)ifich.LoadPGM();

	imadata<BYTE> imaR;
	imadon.seuil_hysteresis (atof(argv[3]),atof(argv[2]),imaR);
    
	imaR.sauve_ImaPGM(argv[5]);
    return(-1);
}

/* -------------------------------- */
int main_gauss_noise(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<5) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<"\n";
		exit(-1);
	}
    char *filename = argv[3];
	fichimage_entree ifich(filename);
	imadata<float> imadon=(imadata<float>)ifich.LoadPGM();

	imadon.add_gauss_noise (atof(argv[2]));
    
	imadon.sauve_ImaPGM(argv[4]);
    return(-1);
}

/* -------------------------------- */
int main_impul_noise(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<5) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<"\n";
		exit(-1);
	}
    char *filename = argv[3];
	fichimage_entree ifich(filename);
	imadata<float> imadon=(imadata<float>)ifich.LoadPGM();

	imadon.add_impul_noise (atof(argv[2]));
    
	imadon.sauve_ImaPGM(argv[4]);
    return(-1);
}

/* -------------------------------- */
int main_filtr_Gauss(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<5) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<"\n";
		exit(-1);
	}
    char *filename = argv[3];
	fichimage_entree ifich(filename);
	imadata<float> imadon=(imadata<float>)ifich.LoadPGM();

	imadon=imadon.filtregaussienne (pow(atof(argv[2]),2),1);
    
	imadon.sauve_ImaPGM(argv[4]);
    return(-1);
}

/* -------------------------------- */
int main_filtr_moyen(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<5) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<"\n";
		exit(-1);
	}
    char *filename = argv[3];
	fichimage_entree ifich(filename);
	imadata<float> imadon=(imadata<float>)ifich.LoadPGM();

	imadon=imadon.filtremoyenne (atoi(argv[2]), atoi(argv[2]));
    
	imadon.sauve_ImaPGM(argv[4]);
    return(-1);
}

/* -------------------------------- */
int main_filtr_median(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<5) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<"\n";
		exit(-1);
	}
    char *filename = argv[3];
	fichimage_entree ifich(filename);
	imadata<float> imadon=(imadata<float>)ifich.LoadPGM();

	imadon=imadon.filtremedian (atoi(argv[2]), atoi(argv[2]));
    
	imadon.sauve_ImaPGM(argv[4]);
    return(-1);
}

/* -------------------------------- */
int main_filtr_Nagao(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<4) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<"\n";
		exit(-1);
	}
    char *filename = argv[2];
	fichimage_entree ifich(filename);
	imadata<float> imadon=(imadata<float>)ifich.LoadPGM();

	imadon=imadon.filtreNagao ();
    
	imadon.sauve_ImaPGM(argv[3]);
    return(-1);
}

/* -------------------------------- */
int main_filtr_SNN(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<5) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<"\n";
		exit(-1);
	}
    char *filename = argv[3];
	fichimage_entree ifich(filename);
	imadata<float> imadon=(imadata<float>)ifich.LoadPGM();

	imadon=imadon.filtreSymNearNeigh (atoi(argv[2]));
    
	imadon.sauve_ImaPGM(argv[4]);
    return(-1);
}

/* -------------------------------- */
int main_filtr_FAS(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<7) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<" "<<argv[5]<<" "<<argv[6]<<"\n";
		exit(-1);
	}
    char *filename = argv[5];
	fichimage_entree ifich(filename);
	imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
	eltstruct S3(3,3); 
	if (atoi(argv[2])==4) {S3(0,0)=0; S3(0,2)=0; S3(2,0)=0; S3(2,2)=0;}
	bool sens=(atoi(argv[4])==0);

	imadon=imadon.filtrealterne(S3,atoi(argv[3]),sens);
    
	imadon.sauve_ImaPGM(argv[6]);
    return(-1);
}

/* -------------------------------- */
int main_psnr(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<7) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<" "<<argv[5]<<" "<<argv[6]<<"\n";
		exit(-1);
	}
    char *filename1 = argv[2], *filename2 = argv[3];
	fichimage_entree ific1(filename1), ific2(filename2);
	imadata<float> imadon=(imadata<float>)ific2.LoadPGM();
	imadata<float> imarec=(imadata<float>)ific1.LoadPGM();

	imarec.psnr(imadon);
    
    return(-1);
}
/* -------------------------------- */
int main_negatif(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<4) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<"\n";
		exit(-1);
	}
    char *filename = argv[2];
	fichimage_entree ifich(filename);
	imabin imadon((imadata<float>)ifich.LoadPGM(),1);

	imadon=imadon.negatif ();
    
	imadon.imaunsignedchar().sauve_ImaPGM(argv[3]);
    return(-1);
}

/* -------------------------------- */
int main_plus_bin(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<5) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<"\n";
		exit(-1);
	}
    char *filename1 = argv[2], *filename2 = argv[3];
	fichimage_entree ific1(filename1), ific2(filename2);
	imabin imadon1((imadata<float>)ific1.LoadPGM(),1), imadon2((imadata<float>)ific2.LoadPGM(),1);

	imadon1=imadon1+imadon2;
    
	imadon1.imaunsignedchar().sauve_ImaPGM(argv[4]);
    return(-1);
}

/* -------------------------------- */
int main_et(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<5) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<"\n";
		exit(-1);
	}
    char *filename1 = argv[2], *filename2 = argv[3];
	fichimage_entree ific1(filename1), ific2(filename2);
	imabin imadon1((imadata<float>)ific1.LoadPGM(),1), imadon2((imadata<float>)ific2.LoadPGM(),1);

	imadon1=imadon1&&imadon2;
    
	imadon1.imaunsignedchar().sauve_ImaPGM(argv[4]);
    return(-1);
}

/* -------------------------------- */
int main_ou(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<5) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<"\n";
		exit(-1);
	}
    char *filename1 = argv[2], *filename2 = argv[3];
	fichimage_entree ific1(filename1), ific2(filename2);
	imabin imadon1((imadata<float>)ific1.LoadPGM(),1), imadon2((imadata<float>)ific2.LoadPGM(),1);

	imadon1=imadon1||imadon2;
    
	imadon1.imaunsignedchar().sauve_ImaPGM(argv[4]);
    return(-1);
}

/* -------------------------------- */
int main_tr_dist(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<5) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<"\n";
		exit(-1);
	}
    char *filename = argv[3];
	fichimage_entree ifich(filename);
	imabin imadon((imadata<float>)ifich.LoadPGM(),1);

	imadata<float> imadist=imadon.Tr_dist (atoi(argv[2]));
    imadist.statbasic(1);

	imadist.sauve_ImaPGM(argv[4]);
    return(-1);
}

/* -------------------------------- */
int main_erode_bin(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<6) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<" "<<argv[5]<<"\n";
		exit(-1);
	}
    char *filename = argv[4];
	fichimage_entree ifich(filename);
	imabin imadon((imadata<float>)ifich.LoadPGM(),1);

	eltstruct S3(3,3); 
	if (atoi(argv[3])==4) {S3(0,0)=0; S3(0,2)=0; S3(2,0)=0; S3(2,2)=0;}

	imadon=imadon.erode (S3, (atoi(argv[2])-1)/2);
    
	imadon.imaunsignedchar().sauve_ImaPGM(argv[5]);
    return(-1);
}

/* -------------------------------- */
int main_dilate_bin(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<6) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<" "<<argv[5]<<"\n";
		exit(-1);
	}
    char *filename = argv[4];
	fichimage_entree ifich(filename);
	imabin imadon((imadata<float>)ifich.LoadPGM(),1);

	eltstruct S3(3,3); 
	if (atoi(argv[3])==4) {S3(0,0)=0; S3(0,2)=0; S3(2,0)=0; S3(2,2)=0;}

	imadon=imadon.dilate (S3, (atoi(argv[2])-1)/2);
    
	imadon.imaunsignedchar().sauve_ImaPGM(argv[5]);
    return(-1);
}

/* -------------------------------- */
int main_ouvre_bin(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<6) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<" "<<argv[5]<<"\n";
		exit(-1);
	}
    char *filename = argv[4];
	fichimage_entree ifich(filename);
	imabin imadon((imadata<float>)ifich.LoadPGM(),1);

	eltstruct S3(3,3); 
	if (atoi(argv[3])==4) {S3(0,0)=0; S3(0,2)=0; S3(2,0)=0; S3(2,2)=0;}

	int k=(atoi(argv[2])-1)/2;
	imadon=imadon.erode (S3,k).dilate(S3,k);
    
	imadon.imaunsignedchar().sauve_ImaPGM(argv[5]);
    return(-1);
}

/* -------------------------------- */
int main_ferme_bin(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<6) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<" "<<argv[5]<<"\n";
		exit(-1);
	}
    char *filename = argv[4];
	fichimage_entree ifich(filename);
	imabin imadon((imadata<float>)ifich.LoadPGM(),1);

	eltstruct S3(3,3); 
	if (atoi(argv[3])==4) {S3(0,0)=0; S3(0,2)=0; S3(2,0)=0; S3(2,2)=0;}

	int k=(atoi(argv[2])-1)/2;
	imadon=imadon.dilate(S3,k).erode (S3,k);
    
	imadon.imaunsignedchar().sauve_ImaPGM(argv[5]);
    return(-1);
}

/* -------------------------------- */
int main_tophat_bin(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<7) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<" "<<argv[5]<<" "<<argv[6]<<"\n";
		exit(-1);
	}
    char *filename = argv[5];
	fichimage_entree ifich(filename);
	imabin imadon((imadata<float>)ifich.LoadPGM(),1), imares;

	int k=(atoi(argv[2])-1)/2, i=atoi(argv[4]);
	eltstruct S3(3,3); 
	if (atoi(argv[3])==4) {S3(0,0)=0; S3(0,2)=0; S3(2,0)=0; S3(2,2)=0;}

	if (i==0) imares=imadon.dilate(S3,k).erode(S3,k)-imadon;
	else      imares=imadon-imadon.erode(S3,k).dilate(S3,k);
    
	imares.imaunsignedchar().sauve_ImaPGM(argv[6]);
    return(-1);
}

/* -------------------------------- */
int main_reconstr_geod(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<6) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<" "<<argv[5]<<"\n";
		exit(-1);
	}
    char *filename1 = argv[4], *filename2 = argv[3];
	fichimage_entree ific1(filename1), ific2(filename2);
	imabin imadon((imadata<float>)ific1.LoadPGM(),1), imamarq((imadata<float>)ific2.LoadPGM(),1);

	int k=atoi(argv[2]);
	imabin imares=imadon.reconstruction_geodesique (imamarq, k);

	imares.imaunsignedchar().sauve_ImaPGM(argv[5]);
    return(-1);
}

/* -------------------------------- */
int main_etiquette_cc(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<5) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<"\n";
		exit(-1);
	}
    char *filename = argv[3];
	fichimage_entree ifich(filename);
	imabin imadon((imadata<float>)ifich.LoadPGM(),1);

	int ncc=0, k=atoi(argv[2]);
	imadata<int> imacc=imadon.composantes_connexes (ncc, k);

	imacc.sauve_ImaPGM(argv[4]);
    return(-1);
}

/* -------------------------------- */
int main_erode_ultime(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<5) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<"\n";
		exit(-1);
	}
    char *filename = argv[3];
	fichimage_entree ifich(filename);
	imabin imadon((imadata<float>)ifich.LoadPGM(),1);

	imadon=imadon.erode_ultime (atoi(argv[2])/4);  // 4-connexite => inoyau=1, 8-connexite => inoyau=2 

	imadon.imaunsignedchar().sauve_ImaPGM(argv[4]);
    return(-1);
}

/* -------------------------------- */
int main_detect_coin(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<5) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<"\n";
		exit(-1);
	}
    char *filename = argv[3];
	fichimage_entree ifich(filename);
	imabin imadon((imadata<float>)ifich.LoadPGM(),1);

	imadata<BYTE> imacoin=imadon.detect_coin (atoi(argv[2]),1,1); 

	imacoin.sauve_ImaPGM(argv[4]);
    return(-1);
}

/* -------------------------------- */
int main_env_convexe(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<4) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<"\n";
		exit(-1);
	}
    char *filename = argv[2];
	fichimage_entree ifich(filename);
	imabin imadon((imadata<float>)ifich.LoadPGM(),1);

	imadon=imadon.enveloppe_convexe (); 

	imadon.imaunsignedchar().sauve_ImaPGM(argv[3]);
    return(-1);
}

/* -------------------------------- */
int main_squelette(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<5) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<"\n";
		exit(-1);
	}
    char *filename = argv[3];
	fichimage_entree ifich(filename);
	imabin imadon((imadata<float>)ifich.LoadPGM(),1);

	imadon=imadon.squelette (atoi(argv[2])); 

	imadon.imaunsignedchar().sauve_ImaPGM(argv[4]);
    return(-1);
}

/* -------------------------------- */
int main_elagage(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<5) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<"\n";
		exit(-1);
	}
    char *filename = argv[3];
	fichimage_entree ifich(filename);
	imabin imadon((imadata<float>)ifich.LoadPGM(),1);

	imadon=imadon.elagage (atoi(argv[2])); 

	imadon.imaunsignedchar().sauve_ImaPGM(argv[4]);
    return(-1);
}

/* -------------------------------- */
int main_zon_infl_geod(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<6) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<" "<<argv[5]<<"\n";
		exit(-1);
	}
    char *filename1 = argv[3], *filename2 = argv[4];
	fichimage_entree ific1(filename1), ific2(filename2);
	imabin imamarq((imadata<float>)ific1.LoadPGM(),1), imadon((imadata<float>)ific2.LoadPGM(),1);

	imadon=imadon.zones_influence_geodesique(imamarq, atoi(argv[2])); 

	imadon.imaunsignedchar().sauve_ImaPGM(argv[5]);
    return(-1);
}

/* -------------------------------- */
int main_class_cmeans(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<5) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<"\n";
		exit(-1);
	}
    char *filename = argv[3];
	fichimage_entree ifich(filename);
	imadata<float> imadon=(imadata<float>)ifich.LoadPGM();

//	pixel<float> pix;
//	sample<pixel<float> > s_pix;
	sampleset<pixel<float> > datset=imadon.dataset();
	statclasses clkmeans=datset.k_means(atoi(argv[2]));
	clkmeans.affiche();
	imalabels imacl(imadon,clkmeans,0.,"ICM");

	imacl.conv2imBYTE().sauve_ImaPGM(argv[4]);
    return(-1);
}

/* -------------------------------- */
int main_class_kppv(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<6) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<" "<<argv[5]<<"\n";
		exit(-1);
	}
    char *filename = argv[4];
	fichimage_entree ifich(filename);
	imadata<float> imadon=(imadata<float>)ifich.LoadPGM();

	sampleset<pixel<float> > datlearn;
	int nt=0,k,ncl=atoi(argv[2]),i,n=0,l,j;
	for (k=0; k<ncl; k++) {
		cout<<" Donner le nombre d'echantillons de la classe "<<k+1<<" : "; 
		cin>>n;
		for (i=0; i<n; i++) {
			cout<<" coordonnees pixel (lig,col) de l'echantillon "<<i<<" : ";
			cin>>l; cin>>j;
			if (l>=0 && j>=0 && l<imadon.nlig() && j<imadon.ncol()) datlearn.ajoute(imadon.pix(l,j),nt++,k);
		}
	}
	datlearn.affiche();
	sampleset<pixel<float> > datset=imadon.dataset();
	datset.k_ppv (datlearn, atoi(argv[3]));
	int nblig=imadon.nlig(), nbcol=imadon.ncol();
	imalabels imacl(nblig,nbcol);
	for (i=0; i<nblig; i++)
		for (j=0; j<nbcol; j++)
			imacl(i,j)=datset.remove(i*nbcol+j).label();

	imacl.conv2imBYTE().sauve_ImaPGM(argv[5]);
    return(-1);
}

/* -------------------------------- */
int main_class_mrf(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<7) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<" "<<argv[5]<<" "<<argv[6]<<"\n";
		exit(-1);
	}
    char *filename = argv[5];
	fichimage_entree ifich(filename);
	imadata<float> imadon=(imadata<float>)ifich.LoadPGM();

	int d=imadon.ncanaux(),k,ncl=atoi(argv[2]),j;
	statclasses Tcl;
	stat1class cl(1,d);
	for (k=0; k<ncl; k++) {
		cl.label()=k+1; 
		cout<<" Donner le vecteur moyenne de la classe "<<k+1<<" : ";
		for (j=0; j<d; j++) cin>>cl.mean(j);
		cout<<" Donner le vecteur variance de la classe "<<k+1<<" : ";
		for (j=0; j<d; j++) cin>>cl.cova(j,j);
		for (j=0; j<d; j++) cl.icova(j,j)=1./cl.cova(j,j);
		cl.dcova()=1.; for (j=0; j<d; j++) cl.dcova()*=cl.cova(j,j);
		Tcl.ajoute(cl);
	}
	Tcl.affiche();
	imalabels imacl(imadon, Tcl, atof(argv[3]), argv[4]);

	imacl.conv2imBYTE(1).sauve_ImaPGM(argv[6]);
    return(-1);
}

/* -------------------------------- */
int main_class_emgibbs(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<5) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<"\n";
		exit(-1);
	}
    char *filename = argv[3];
	fichimage_entree ifich(filename);
	imadata<float> imadon=(imadata<float>)ifich.LoadPGM();

	imalabels imacl(imadon, atoi(argv[2]));

	imacl.conv2imBYTE(1).sauve_ImaPGM(argv[4]);
    return(-1);
}

/* -------------------------------- */
int main_segm_classif(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<4) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<"\n";
		exit(-1);
	}
    char *filename = argv[2];
	fichimage_entree ifich(filename);
	imalabels imalab((imadata<BYTE>)ifich.LoadPGM());

	imaregions imareg(imalab);

	if (imareg.nregions()<256) imareg.conv2imBYTE().sauve_ImaPGM(argv[3]);
    return(-1);
}

/* -------------------------------- */
int main_image_gradient(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<6) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<" "<<argv[5]<<"\n";
		exit(-1);
	}
    char *filename = argv[3];
	fichimage_entree ifich(filename);
	imadata<float> imadon=(imadata<float>)ifich.LoadPGM();

	int m=atoi(argv[2]);
	char *masq;
	switch (m) {
		case 0:  masq="Prewitt"; break;
		case 2:  masq="MDIF"; break;
		default: masq="Sobel";
	}
	imadata<float> imadir(imadon);
	imadata<float> imagrd=imadon.gradient(imadir,masq);

	imagrd.sauve_ImaPGM(argv[4]);
	float coef=90/(float)PI;
	imadir=(imadir*coef)+90.;
	cout<<" les valeurs de l'image de la direction du gradient sont 1<->2deg.\n";
	imadir.sauve_ImaPGM(argv[5]);
    return(-1);
}

/* -------------------------------- */
int main_cont_gradient(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<6) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<" "<<argv[5]<<"\n";
		exit(-1);
	}
    char *filename = argv[4];
	fichimage_entree ifich(filename);
	imadata<float> imadon=(imadata<float>)ifich.LoadPGM();

	int m=atoi(argv[2]);
	char *masq;
	switch (m) {
		case 0:  masq="Prewitt"; break;
		case 2:  masq="MDIF"; break;
		default: masq="Sobel";
	}
	unsigned short int lgmax=atoi(argv[3]);
	imacontours imacont(imadon,masq,lgmax,1);

	imacont.conv2im1dBYTE(0).sauve_ImaPGM(argv[5]);
    return(-1);
}

/* -------------------------------- */
int main_cont_laplacien(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<6) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<" "<<argv[5]<<"\n";
		exit(-1);
	}
    char *filename = argv[4];
	fichimage_entree ifich(filename);
	imadata<float> imadon=(imadata<float>)ifich.LoadPGM();

	int mG=atoi(argv[2]), mL=atoi(argv[3]);
	char *masqG, *masqL;
	switch (mG) {
		case 0:  masqG="Prewitt"; break;
		case 2:  masqG="MDIF"; break;
		default: masqG="Sobel";
	}
	switch (mL) {
		case 4:  masqL="4connex"; break;
		default: masqL="8connex";
	}
	imacontours imacont(imadon,masqL,masqG,1);

	imacont.conv2im1dBYTE(0).sauve_ImaPGM(argv[5]);
    return(-1);
}

/* -------------------------------- */
int main_cont_optimal(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<5) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<"\n";
		exit(-1);
	}
    char *filename = argv[4];
	fichimage_entree ifich(filename);
	imadata<float> imadon=(imadata<float>)ifich.LoadPGM();

	int f=atoi(argv[2]);
	char *filtre;
	switch (f) {
	case 0: filtre="Deriche"; break;
	case 1: filtre="Shen"; break;
	}
	imacontours imacont(imadon,atof(argv[3]),filtre,1);

	imacont.conv2im1dBYTE(0).sauve_ImaPGM(argv[4]);
    return(-1);
}

/* -------------------------------- */
int main_affine_contours(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<5) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<"\n";
		exit(-1);
	}
    char *filename1 = argv[2], *filename2 = argv[3];
	fichimage_entree ific1(filename1), ific2(filename2);
	imadata<float> imadon=(imadata<float>)ific2.LoadPGM();

	char *masq="Sobel";
	imadata<float> imadir(imadon);
	imadata<float> imagrd=imadon.gradient(imadir,masq);
	imacontours imacont((imadata<BYTE>)ific1.LoadPGM());

	imacont.imagrd_maxloc(imagrd,imadir);

	imacont.conv2im1dBYTE(0).sauve_ImaPGM(argv[4]);
    return(-1);
}

/* -------------------------------- */
int main_prolonge_contours(int argc, char *argv[])
/* -------------------------------- */
{   
    if(argc<6) {
		cout<<"Il manque un argument!: "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<" "<<argv[5]<<"\n";
		exit(-1);
	}
    char *filename1 = argv[3], *filename2 = argv[4];
	fichimage_entree ific1(filename1), ific2(filename2);
	imadata<float> imagrd=(imadata<float>)ific2.LoadPGM();
	imacontours imacont((imadata<BYTE>)ific1.LoadPGM());

	imacont.prolongecontours(imagrd,atoi(argv[2]));

	imacont.conv2im1dBYTE(0).sauve_ImaPGM(argv[5]);
    return(-1);
}

/* --------------------------- */
//void main(int argc, char *argv[])
void main_TdI(int argc, char *argv[])
/* --------------------------- */
{           
    if ((argc==1)|| ((argc==2) && (!strcmp(argv[1],"-Menu")))) {
        menu();       
        return;
    } else {
		if(!strcmp(argv[1],"-inv_Ngr"))			{ main_inv_ngr				(argc, argv); return;}
		if(!strcmp(argv[1],"-bin_Sup"))			{ main_bin_sup				(argc, argv); return;}
		if(!strcmp(argv[1],"-bin_Inf"))			{ main_bin_inf				(argc, argv); return;}
		if(!strcmp(argv[1],"-seuil_hyst"))		{ main_seuil_hysteresis 	(argc, argv); return;}
		if(!strcmp(argv[1],"-gauss_noise"))		{ main_gauss_noise			(argc, argv); return;}
		if(!strcmp(argv[1],"-impul_noise"))		{ main_impul_noise			(argc, argv); return;}
		if(!strcmp(argv[1],"-filtr_Gauss"))		{ main_filtr_Gauss			(argc, argv); return;}
		if(!strcmp(argv[1],"-filtr_moyen"))		{ main_filtr_moyen			(argc, argv); return;}
		if(!strcmp(argv[1],"-filtr_median"))	{ main_filtr_median			(argc, argv); return;}
		if(!strcmp(argv[1],"-filtr_Nagao"))		{ main_filtr_Nagao			(argc, argv); return;}
		if(!strcmp(argv[1],"-filtr_SNN"))		{ main_filtr_SNN			(argc, argv); return;}
		if(!strcmp(argv[1],"-filtr_FAS"))		{ main_filtr_FAS			(argc, argv); return;}
		if(!strcmp(argv[1],"-psnr"))			{ main_psnr					(argc, argv); return;}
		if(!strcmp(argv[1],"-negatif"))			{ main_negatif				(argc, argv); return;}
		if(!strcmp(argv[1],"-plus"))			{ main_plus_bin				(argc, argv); return;}
		if(!strcmp(argv[1],"-et"))				{ main_et					(argc, argv); return;}
		if(!strcmp(argv[1],"-ou"))				{ main_ou					(argc, argv); return;}
		if(!strcmp(argv[1],"-trf_dist"))		{ main_tr_dist				(argc, argv); return;}
		if(!strcmp(argv[1],"-erode_bin"))		{ main_erode_bin			(argc, argv); return;}
		if(!strcmp(argv[1],"-dilate_bin"))		{ main_dilate_bin			(argc, argv); return;}
		if(!strcmp(argv[1],"-ouvre_bin"))		{ main_ouvre_bin			(argc, argv); return;}
		if(!strcmp(argv[1],"-ferme_bin"))		{ main_ferme_bin			(argc, argv); return;}
		if(!strcmp(argv[1],"-tophat_bin"))		{ main_tophat_bin			(argc, argv); return;}
		if(!strcmp(argv[1],"-recon_geod"))		{ main_reconstr_geod		(argc, argv); return;}
		if(!strcmp(argv[1],"-etiq_cc"))			{ main_etiquette_cc			(argc, argv); return;}
		if(!strcmp(argv[1],"-erode_ult"))		{ main_erode_ultime			(argc, argv); return;}
		if(!strcmp(argv[1],"-detect_coin"))		{ main_detect_coin			(argc, argv); return;}
		if(!strcmp(argv[1],"-env_convexe"))		{ main_env_convexe			(argc, argv); return;}
		if(!strcmp(argv[1],"-squelette"))		{ main_squelette			(argc, argv); return;}
		if(!strcmp(argv[1],"-elagage"))			{ main_elagage				(argc, argv); return;}
		if(!strcmp(argv[1],"-zig"))				{ main_zon_infl_geod		(argc, argv); return;}
		if(!strcmp(argv[1],"-class_cmeans"))	{ main_class_cmeans			(argc, argv); return;}
		if(!strcmp(argv[1],"-class_kppv"))		{ main_class_kppv			(argc, argv); return;}
		if(!strcmp(argv[1],"-class_MRF"))		{ main_class_mrf			(argc, argv); return;}
		if(!strcmp(argv[1],"-class_EMG"))		{ main_class_emgibbs		(argc, argv); return;}
		if(!strcmp(argv[1],"-seg_class"))		{ main_segm_classif			(argc, argv); return;}
		if(!strcmp(argv[1],"-ima_grd"))			{ main_image_gradient		(argc, argv); return;}
		if(!strcmp(argv[1],"-cont_grd"))		{ main_cont_gradient		(argc, argv); return;}
		if(!strcmp(argv[1],"-cont_lpl"))		{ main_cont_laplacien		(argc, argv); return;}
		if(!strcmp(argv[1],"-cont_opt"))		{ main_cont_optimal			(argc, argv); return;}
		if(!strcmp(argv[1],"-aff_cont"))		{ main_affine_contours		(argc, argv); return;}
		if(!strcmp(argv[1],"-prol_cont"))		{ main_prolonge_contours	(argc, argv); return;}
//		if(!strcmp(argv[1],"-1D"))				{ main_texture1D			(argc, argv); return;}   
//		if(!strcmp(argv[1],"-2D"))				{ main_texture2D			(argc, argv); return;} 
//		if(!strcmp(argv[1],"-Corr"))			{ main_correlation			(argc, argv); return;}
//		if(!strcmp(argv[1],"-segm_classif"))	{ main_segm_classif			(argc, argv); return;}
//		if(!strcmp(argv[1],"-croissance"))		{ main_croissance			(argc, argv); return;}
//		if(!strcmp(argv[1],"-segmentation"))	{ main_segmentation			(argc, argv); return;}
//		if(!strcmp(argv[1],"-rgb2ist"))			{ main_RGB2IST				(argc, argv); return;}
//		if(!strcmp(argv[1],"-rgb2lab"))			{ main_RGB2LAB				(argc, argv); return;}
//		if(!strcmp(argv[1],"-rgb2l1"))			{ main_RGB2L1				(argc, argv); return;}
//		if(!strcmp(argv[1],"-sobel"))			{ main_sobel				(argc, argv); return;}
//		if(!strcmp(argv[1],"-dizenzo"))			{ main_diZenzo				(argc, argv); return;}
//		if(!strcmp(argv[1],"-median_marg"))		{ main_MedianMarginal		(argc, argv); return;}
//		if(!strcmp(argv[1],"-median_vect"))		{ main_MedianVectoriel		(argc, argv); return;}
//		if(!strcmp(argv[1],"-reconG_bin"))		{ main_Reconstr_Geo_bin		(argc, argv); return;}
//		if(!strcmp(argv[1],"-erode_gr"))		{ main_Erosion_gr			(argc, argv); return;}
//		if(!strcmp(argv[1],"-dilate_gr"))		{ main_Dilatation_gr		(argc, argv); return;}
//		if(!strcmp(argv[1],"-gmorph"))			{ main_gmorph				(argc, argv); return;}
//		if(!strcmp(argv[1],"-lmorph"))			{ main_lmorph				(argc, argv); return;}
//		if(!strcmp(argv[1],"-FAS"))				{ main_FAS					(argc, argv); return;}
//		if(!strcmp(argv[1],"-reconG_gr"))		{ main_reconstrGeodgr		(argc, argv); return;}
		cout<<"Cette option n'existe pas!\n";           
    }
    return;
}
