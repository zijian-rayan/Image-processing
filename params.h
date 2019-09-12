#ifndef _PARAMS_H
#define _PARAMS_H

#include "stdafx.h"
#include "tdi/tdi.h"
#include "SharedData.h"

// classe pour rediriger cout vers la console actuelle ou vers n'importe quoi d'autre
//class MyStr: public streambuf
//{
//  typedef char_traits<char> _Tr;
//public:
//  virtual int_type overflow(int_type c= _Tr::eof())
//  {
//    //    ::MessageBox(0,"Overflow","Message",MB_OK);
//    char v=_Tr::to_char_type(c);
//    DWORD nbc;
//    WriteConsole(hcons,&v,1,&nbc,0);
//    return _Tr::not_eof(c);
//  }
//  virtual int sync()
//  {
//    //    ::MessageBox(0,"Sync","Message",MB_OK);
//    return 0;
//  }
//  
//  static bool ChangeCout()
//  {
//    static MyStr my_str;
//    //cout=ostream(&my_str);
//    return true;
//  }
//  static bool done;
//  static HANDLE hcons;
//};

// parametres utilisés dans l'appel aux fonctions de traitement d'images
struct PARAM_FNAME
{
  char* fname;
};

struct PARAM_INT
{
  int val, val_min, val_max, val_def;
};

struct PARAM_FLOAT
{
  float val, val_min, val_max, val_def;
};

struct PARAM_SEL
{
  BYTE sel_cur, sel_def;
  const char* sel_str;
};

// classe qui encapsule n'importe quel paramètre, son type et sa description
class Param_Type
{
public:
  BYTE type;
  const char* desc;
  union
  {
    PARAM_FNAME  p_fname;
    PARAM_INT    p_int;
    PARAM_FLOAT  p_float;
    PARAM_SEL    p_sel;
  };
  Param_Type(): type(255) {}
  Param_Type(const char* _desc, bool _is_inp): type(_is_inp?0:1), desc(_desc)
  { p_fname.fname=0; }
  Param_Type(const char* _desc, int _val_min, int _val_max, int _val_def): type(2), desc(_desc)
  { p_int.val=p_int.val_def=_val_def; p_int.val_min=_val_min; p_int.val_max=_val_max; }
  Param_Type(const char* _desc, float _val_min, float _val_max, float _val_def): type(3), desc(_desc)
  { p_float.val=p_float.val_def=_val_def; p_float.val_min=_val_min; p_float.val_max=_val_max; }
  Param_Type(const char* _desc, BYTE _sel_def, const char* _sel_str): type(4), desc(_desc)
  { p_sel.sel_cur=p_sel.sel_def=_sel_def; p_sel.sel_str=_sel_str; }
  ~Param_Type()
  {
    if((type==0 || type==1) && p_fname.fname!=0) delete[] p_fname.fname;
  }
  static Param_Type para[];
private:
  Param_Type(const Param_Type& p) {} // copie publique C++ interdite
};

// classe qui encapsule la fonction à appeler avec paramètres, sa description
// elle a comme membres statiques toutes les fonctions à lancer en 2ème thread par Launch
class Func_Type
{
public:
  BYTE cat;           // categorie (groupe) de fonctions
  const char* name;   // nom de la fonction
  const char* desc;   // description a afficher
  DWORD nb_param;     // nb. de parametres
  Param_Type* params; // tableau de parametres de type polymorphe Param_Type
  bool (*pfct)(int nbp, Param_Type* p); // ptr. de fonction a executer 
  // dans _sparam il faut mettre toujours les sorties à la fin !
  Func_Type(BYTE _cat, const char* _name, const char* _desc, const char* _sparam, bool (*_pfct)(int nbp, Param_Type* p)):
  name(_name), desc(_desc), pfct(_pfct), cat(_cat)
  {
    nb_param=strlen(_sparam);
    // à verifier que _sparam n'a pas des valeurs > taille de para
    params=new Param_Type[nb_param];
    for(DWORD i=0; i<nb_param; i++)
      memcpy(params+i,Param_Type::para+_sparam[i],sizeof(Param_Type));
  }
  static DWORD WINAPI Launch(void* arg) // pas de MFC à partir de là !
  {
    //AllocConsole();
    //MyStr::hcons=GetStdHandle(STD_OUTPUT_HANDLE);
    Func_Type* pfunc=(Func_Type*)arg;
    char title[256]; ::CharToOem(pfunc->name,title);
    SetConsoleTitle(title);
    pfunc->pfct(pfunc->nb_param,pfunc->params);
    ::MessageBox(0,pfunc->name,"Fin d'execution",MB_OK);
    //FreeConsole();
    return 0;
  }
  static bool op_operation_ima(int nbp, Param_Type* p)
  {
		bool ajust_dyn=1;
    fichimage_entree ific1(p[0].p_fname.fname), ific2(p[1].p_fname.fname);
    imadata<float> imadon1=(imadata<float>)ific1.LoadPGM(), imadon2=(imadata<float>)ific2.LoadPGM();
    if (p[2].p_sel.sel_cur==0) imadon1=imadon1+imadon2;
    else {
      if (p[2].p_sel.sel_cur<=2) {  
        imadon1=imadon1-imadon2;
        if (p[2].p_sel.sel_cur==2) imadon1=imadon1.absI();
      } else {
				float eps=1.f, pc=0.9f;
        if (imadon2.minI()>eps) imadon1=imadon1/imadon2;
        else imadon1=(imadon1+eps)/(imadon2+eps);
        imadon1.statbasic(1); float xmax=(float)imadon1.maxI();
        if (xmax<=1) imadon1=imadon1*100.f;
        else 
          if (xmax<=25) imadon1=imadon1*10.f;
				if (ajust_dyn) {
					float s=imadon1.seuil_percentile (pc), fact=255.f/maxi(s,1.f); 
					imadon1=imadon1*fact;
				}
      }
    }
    imadon1.imaunsignedchar().sauve_ImaPGM(p[3].p_fname.fname);
    return true;
  }
  static bool op_operation_scal(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
    float x=p[1].p_float.val;
    if (p[2].p_sel.sel_cur==1) imadon=imadon*x;
    else imadon=imadon+x;
		int i,j,nl=imadon.nlig(),nc=imadon.ncol();
		for (i=0; i<nl; i++) for (j=0; j<nc; j++) imadon(i,j)=mini(maxi(0.f,imadon(i,j)),255.f);
    imadon.sauve_ImaPGM(p[3].p_fname.fname);
    return true;
  }
  static bool op_inv_ngr(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
    for (int k=0; k<imadon.ncanaux(); k++) imadon.inverse_ngr(255,k);
    imadon.sauve_ImaPGM(p[1].p_fname.fname);
    return true;
  }
  static bool op_etir_dyn(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
    imadon.statbasic();
    for (int i=0; i<imadon.ncanaux(); i++) {
      float a=255/(float)(imadon.maxI(i)-imadon.minI(i)), b=-a*(float)imadon.minI(i);
      imadon.copiecanal(i,(imadon.mult_canal(a,i)).add_canal(b,i),i);
    }
    imadon.sauve_ImaPGM(p[1].p_fname.fname);
    return true;
  }
  static bool op_egal_histo(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
    cout<<" statistiques image origine :\n"; imadon.statbasic(1);
    int i,n=imadon.ncanaux();
    float a,b;
    for (i=0; i<n; i++) {
      a=(float)(255./(imadon.maxI(i)-imadon.minI(i))); b=-(float)(a*imadon.minI(i));
      imadon.copiecanal(i,(imadon.mult_canal(a,i)).add_canal(b,i),i);
    }
    /*cout<<" stat avant egalisation :\n"; */imadon.statbasic();
    imadon.egalisehisto(128,0);
    /*cout<<" stat apres egalisation :\n"; */imadon.statbasic();
    for (i=0; i<n; i++) {
      a=(float)(255./(imadon.maxI(i)-imadon.minI(i))); b=-(float)(a*imadon.minI(i));
      imadon.copiecanal(i,(imadon.mult_canal(a,i)).add_canal(b,i),i);
    }
    cout<<" statistiques apres re-etirement dynamique :\n"; imadon.statbasic(1);
    imadon.sauve_ImaPGM(p[1].p_fname.fname);
    return true;
  }
  static bool op_masq_ima(int nbp, Param_Type* p)
  {
    fichimage_entree ific1(p[0].p_fname.fname), ific2(p[1].p_fname.fname);
    imadata<float> imadon1=(imadata<float>)ific1.LoadPGM(), imadon2=(imadata<float>)ific2.LoadPGM();
    int i,j,k,nblig=mini(imadon1.nlig(),imadon2.nlig()),nbcol=mini(imadon1.ncol(),imadon2.ncol()),dim=imadon1.ncanaux(); 
    for (i=0; i<nblig; i++)
      for (j=0; j<nbcol; j++) 
        if (imadon2(i,j)==0) for (k=0; k<dim; k++) imadon1(i,j,k)=0;
    imadon1.sauve_ImaPGM(p[2].p_fname.fname);
    return true;
  }
  static bool op_bin_sup(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
    imabin imab(imadon, p[1].p_float.val);
    imadata<BYTE>(imab).sauve_ImaPGM(p[2].p_fname.fname);
    return true;
  }
  static bool op_bin_inf(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
    imabin imab(imadon, p[1].p_float.val+0.01f);
    imab=imab.negatif();
    imadata<BYTE>(imab).sauve_ImaPGM(p[2].p_fname.fname);
    return true;
  }
  static bool op_seuil_hysteresis(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
    imadata<BYTE> imaR;
    imadon.seuil_hysteresis(p[2].p_float.val,p[1].p_float.val,imaR);
    (imaR*255).sauve_ImaPGM(p[3].p_fname.fname);
    return true;
  }
  static bool op_bin_inter(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
    imabin imab1(imadon,p[1].p_float.val+0.01f), imab2(imadon,p[2].p_float.val);
    imab1=imab1.negatif();
    imadata<BYTE>(imab1&&imab2).sauve_ImaPGM(p[3].p_fname.fname);
    return true;
  }
  static bool op_bin_pct_sup(int nbp, Param_Type* p)
  {
    fichimage_entree ific1(p[0].p_fname.fname), ific2(p[1].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ific1.LoadPGM();
    imabin imaROI((imadata<float>)ific2.LoadPGM(),1);
    float s=imadon.seuil_percentile(imaROI,1.f-p[2].p_float.val);
    imabin imab(imadon,s);
    imadata<BYTE>(imab).sauve_ImaPGM(p[3].p_fname.fname);
    return true;
  }
  static bool op_bin_pct_inf(int nbp, Param_Type* p)
  {
    fichimage_entree ific1(p[0].p_fname.fname), ific2(p[1].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ific1.LoadPGM();
    imabin imaROI((imadata<float>)ific2.LoadPGM(),1);
    float s=imadon.seuil_percentile(imaROI,p[2].p_float.val);
    imabin imab(imadon,s);
    imab=imab.negatif()+imaROI;
    imadata<BYTE>(imab).sauve_ImaPGM(p[3].p_fname.fname);
    return true;
  }
  static bool op_bin_otsu(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
    imabin imab; imadon.seuil_otsu(imab);
    imadata<BYTE>(imab).sauve_ImaPGM(p[1].p_fname.fname);
    return true;
  }
  static bool op_bin_maxAbs(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
    imabin imab(imadon,(const float)imadon.maxI());
    imadata<BYTE>(imab).sauve_ImaPGM(p[1].p_fname.fname);
    return true;
  }
  static bool op_bin_maxReg(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
		imabin imadb(imadon,0);
		imabin imab=imadb.erode_ultime(imadon,8);
    imadata<BYTE>(imab).sauve_ImaPGM(p[1].p_fname.fname);
    return true;
  }
  static bool op_gauss_noise(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
    imadon.add_gauss_noise(p[1].p_float.val*p[1].p_float.val);
    imadon.imaunsignedchar().sauve_ImaPGM(p[2].p_fname.fname);
    return true;
  }
  static bool op_impul_noise(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
    imadon.add_impul_noise(p[1].p_float.val);
    imadon.imaunsignedchar().sauve_ImaPGM(p[2].p_fname.fname);
    return true;
  }
  static bool op_filtr_Gauss(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
    imadon=imadon.filtregaussienne(p[1].p_float.val,1);
    imadon.sauve_ImaPGM(p[2].p_fname.fname);
    return true;
  }
  static bool op_filtr_moyen(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
    int w=p[1].p_sel.sel_cur*2+3;
    imadon=imadon.filtremoyenne(w,w);
    imadon.sauve_ImaPGM(p[2].p_fname.fname);
    return true;
  }
  static bool op_filtr_median(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
    int w=p[1].p_sel.sel_cur*2+3;
    imadon=imadon.filtremedian(w,w);
    imadon.sauve_ImaPGM(p[2].p_fname.fname);
    return true;
  }
  static bool op_filtr_Nagao(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
    int w=p[1].p_sel.sel_cur*2+3;
    imadon=imadon.filtreNagao(w,0,7,1);
    imadon.sauve_ImaPGM(p[2].p_fname.fname);
    return true;
  }
	  static bool op_filtr_NagaoMed(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
    int w=p[1].p_sel.sel_cur*2+3;
    imadon=imadon.filtreNagao(w,1,11,1);
    imadon.sauve_ImaPGM(p[2].p_fname.fname);
    return true;
  }
  static bool op_filtr_SNN(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
    int w=p[1].p_sel.sel_cur*2+3;
    imadon=imadon.filtreSymNearNeigh(w);
    imadon.sauve_ImaPGM(p[2].p_fname.fname);
    return true;
  }
  static bool op_filtr_FAS(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
    eltstruct S3(3,3); 
    if(p[2].p_sel.sel_cur==0) {S3(0,0)=0; S3(0,2)=0; S3(2,0)=0; S3(2,2)=0;}
    imadon=imadon.filtrealterne(S3,p[3].p_int.val,p[1].p_sel.sel_cur==0);
    imadon.sauve_ImaPGM(p[4].p_fname.fname);
    return true;
  }
  static bool op_psnr(int nbp, Param_Type* p)
  {
    fichimage_entree ific1(p[0].p_fname.fname), ific2(p[1].p_fname.fname);
    imadata<float> imarec=(imadata<float>)ific1.LoadPGM(), imadon=(imadata<float>)ific2.LoadPGM();
    imarec.psnr(imadon);
    return true;
  }
  static bool op_distHist(int nbp, Param_Type* p)
  {
    fichimage_entree ific1(p[0].p_fname.fname), ific2(p[1].p_fname.fname);
    imadata<float> imadon1=(imadata<float>)ific1.LoadPGM(), imadon2=(imadata<float>)ific2.LoadPGM();
		int nbin=64;
		bool affich=1, ioptim=0;
		imadon1.histogramme (nbin,affich,ioptim); imadon2.histogramme (nbin,affich,ioptim);
		cout<<imadon1.distHisto_Bhattacharyya(imadon2);
    return true;
  }
  static bool op_rehaus_contra(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
    eltstruct S3(3,3); 
    if(p[1].p_sel.sel_cur==0) {S3(0,0)=0; S3(0,2)=0; S3(2,0)=0; S3(2,2)=0;}
    int n=p[2].p_sel.sel_cur+1;
    float alpha=p[3].p_float.val;
/*	imadata<float> imaerode=imadon.erode(S3,n);
    imadata<float> imadilate=imadon.dilate(S3,n);
    int nbli=imadon.nlig(), nbco=imadon.ncol(), i,j;
    float delta;
    for (i=0; i<nbli; i++)
      for (j=0; j<nbco; j++) {
        delta=imadilate(i,j)-imaerode(i,j);
        if (imadon(i,j)<imaerode(i,j)+alpha*delta) imadon(i,j)=imaerode(i,j);
        else
          if (imadon(i,j)>imadilate(i,j)-alpha*delta) imadon(i,j)=imadilate(i,j);
		}*/
    imadon.rehausse_contraste(S3,n,alpha,alpha).sauve_ImaPGM(p[4].p_fname.fname);
    return true;
  }
  static bool op_xor_bin(int nbp, Param_Type* p)
  {
    fichimage_entree ific1(p[0].p_fname.fname), ific2(p[1].p_fname.fname);
    imabin imadon1((imadata<float>)ific1.LoadPGM(),1), imadon2((imadata<float>)ific2.LoadPGM(),1);
    imadon1=imadon1+imadon2;
    imadon1.imaunsignedchar().sauve_ImaPGM(p[2].p_fname.fname);
    return true;
  }
  static bool op_and_bin(int nbp, Param_Type* p)
  {
    fichimage_entree ific1(p[0].p_fname.fname), ific2(p[1].p_fname.fname);
    imabin imadon1((imadata<float>)ific1.LoadPGM(),1), imadon2((imadata<float>)ific2.LoadPGM(),1);
    (imadon1&&imadon2).imaunsignedchar().sauve_ImaPGM(p[2].p_fname.fname);
    return true;
  }
  static bool op_or_bin(int nbp, Param_Type* p)
  {
    fichimage_entree ific1(p[0].p_fname.fname), ific2(p[1].p_fname.fname);
    imabin imadon1((imadata<float>)ific1.LoadPGM(),1), imadon2((imadata<float>)ific2.LoadPGM(),1);
    imadon1=imadon1||imadon2;
    imadon1.imaunsignedchar().sauve_ImaPGM(p[2].p_fname.fname);
    return true;
  }
  static bool op_not_bin(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imabin imadon((imadata<float>)ifich.LoadPGM(),1);
    imadon=imadon.negatif();
    imadon.imaunsignedchar().sauve_ImaPGM(p[1].p_fname.fname);
    return true;
  }
  static bool op_tr_dist(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imabin imadon((imadata<float>)ifich.LoadPGM(),1);
    imadata<float> imadist=imadon.Tr_dist(p[1].p_int.val+1);
    imadist.statbasic(1);
    imadist.sauve_ImaPGM(p[2].p_fname.fname);
    return true;
  }
  static bool op_tr_dist_geo(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imabin imadon((imadata<float>)ifich.LoadPGM(),1);
    fichimage_entree ificm(p[1].p_fname.fname);
    imabin imamsk((imadata<float>)ificm.LoadPGM(),1);
    imadata<float> imadist=imadon.Tr_dist(imamsk,p[2].p_int.val+1);
    imadist.statbasic(1);
    imadist.sauve_ImaPGM(p[3].p_fname.fname);
    return true;
  }
	static bool op_erode_bin(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imabin imadon((imadata<float>)ifich.LoadPGM(),1);
    eltstruct S3(3,3); 
    if(p[2].p_sel.sel_cur==0) {S3(0,0)=0; S3(0,2)=0; S3(2,0)=0; S3(2,2)=0;}
    imadon=imadon.erode(S3, p[1].p_sel.sel_cur+1);
    imadon.imaunsignedchar().sauve_ImaPGM(p[3].p_fname.fname);
    return true;
  }
  static bool op_dilate_bin(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imabin imadon((imadata<float>)ifich.LoadPGM(),1);
    eltstruct S3(3,3); 
    if(p[2].p_sel.sel_cur==0) {S3(0,0)=0; S3(0,2)=0; S3(2,0)=0; S3(2,2)=0;}
    imadon=imadon.dilate(S3, p[1].p_sel.sel_cur+1);
    imadon.imaunsignedchar().sauve_ImaPGM(p[3].p_fname.fname);
    return true;
  }
  static bool op_ouvre_bin(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imabin imadon((imadata<float>)ifich.LoadPGM(),1);
    eltstruct S3(3,3); 
    if(p[2].p_sel.sel_cur==0) {S3(0,0)=0; S3(0,2)=0; S3(2,0)=0; S3(2,2)=0;}
    imadon=imadon.erode(S3, p[1].p_sel.sel_cur+1).dilate(S3,p[1].p_sel.sel_cur+1);
    imadon.imaunsignedchar().sauve_ImaPGM(p[3].p_fname.fname);
    return true;
  }
  static bool op_ferme_bin(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imabin imadon((imadata<float>)ifich.LoadPGM(),1);
    eltstruct S3(3,3); 
    if(p[2].p_sel.sel_cur==0) {S3(0,0)=0; S3(0,2)=0; S3(2,0)=0; S3(2,2)=0;}
    imadon=imadon.dilate(S3, p[1].p_sel.sel_cur+1).erode(S3,p[1].p_sel.sel_cur+1);
    imadon.imaunsignedchar().sauve_ImaPGM(p[3].p_fname.fname);
    return true;
  }
  static bool op_tophat_bin(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imabin imadon((imadata<float>)ifich.LoadPGM(),1), imares;
    int k=p[2].p_sel.sel_cur+1;
    eltstruct S3(3,3); 
    if(p[3].p_sel.sel_cur==0) {S3(0,0)=0; S3(0,2)=0; S3(2,0)=0; S3(2,2)=0;}
    if(p[1].p_sel.sel_cur==0) imares=imadon.dilate(S3,k).erode(S3,k)-imadon;
    else imares=imadon-imadon.erode(S3,k).dilate(S3,k);
    imares.imaunsignedchar().sauve_ImaPGM(p[4].p_fname.fname);
    return true;
  }
  static bool op_reconst_geod_bin(int nbp, Param_Type* p)
  {
    fichimage_entree ific1(p[0].p_fname.fname), ific2(p[1].p_fname.fname);
    imabin imadon((imadata<float>)ific1.LoadPGM(),1), imamarq((imadata<float>)ific2.LoadPGM(),1);
    int k=4;
    if (p[2].p_sel.sel_cur==1) k=8;
    imadon.reconstruction_geodesique (imamarq,k).imaunsignedchar().sauve_ImaPGM(p[3].p_fname.fname);
    return true;
  }
/*  static bool op_etiquette_cc(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imabin imadon((imadata<float>)ifich.LoadPGM(),1);
    int ncc=0, k=4;
    if (p[1].p_sel.sel_cur==1) k=8;
    imadata<int> imacc=imadon.composantes_connexes (ncc, k);
    imacc.sauve_ImaPGM(p[2].p_fname.fname);
    return true;
  }*/
	static bool op_etiquette_cc(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imabin imadon((imadata<float>)ifich.LoadPGM(),1);
    int ncc=0, k=4;
    if (p[1].p_sel.sel_cur==1) k=8;
    imadata<int> imacc;
    if (p[2].p_sel.sel_cur==0) {imacc=imadon.composantes_connexes (ncc, k); cout<<" algo rapide\n";}
    else {imacc=imadon.CompoConnexes_MM (ncc, k);cout<<" algo Morpho Math\n";}
		cout<<" ******** nbre de composantes connexes = "<<ncc<<" ********\n";
    imacc.sauve_ImaPGM(p[3].p_fname.fname);
    return true;
  }
	static bool op_elimin_cc(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imabin imadon((imadata<float>)ifich.LoadPGM(),1);
    imaregions imareg(imadon);
    imareg.suppression_trop_petites_regions(p[1].p_int.val);
    for (int i=0; i<imadon.nlig(); i++)
      for (int j=0; j<imadon.ncol(); j++) if (imareg(i,j)==0) imadon(i,j)=0;
    imadon.imaunsignedchar().sauve_ImaPGM(p[2].p_fname.fname);
    return true;
  }
	static bool op_etiquette_boxes(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imabin imadon((imadata<float>)ifich.LoadPGM(),1);
    int ncc=0;
    imadata<int> imabox=imadon.boxes_rectangulaires (ncc); 
    imabox.sauve_ImaPGM(p[1].p_fname.fname);
    return true;
  }
  static bool op_erode_ultime(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imabin imadon((imadata<float>)ifich.LoadPGM(),1);  
    imadon=imadon.erode_ultime (p[1].p_sel.sel_cur+1);  // 4-connexite => inoyau=1, 8-connexite => inoyau=2   
    imadon.imaunsignedchar().sauve_ImaPGM(p[2].p_fname.fname);
    return true;
  }
  static bool op_detect_coin(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imabin imadon((imadata<float>)ifich.LoadPGM(),1);
    imadata<BYTE> imacoin=imadon.detect_coin (p[1].p_int.val,1,1);
    imacoin.sauve_ImaPGM(p[2].p_fname.fname);
    return true;
  }
  static bool op_env_convexe(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imabin imadon((imadata<float>)ifich.LoadPGM(),1);
    imadon=imadon.enveloppe_convexe (); 
    imadon.imaunsignedchar().sauve_ImaPGM(p[1].p_fname.fname);
    return true;
  }
  static bool op_squelette(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imabin imadon((imadata<float>)ifich.LoadPGM(),1);
    int k=8;
    if (p[1].p_sel.sel_cur==0) k=4;
		BYTE imarq=2/*1*/;
    imadon=imadon.squelette(k,imarq); 
    imadon.imaunsignedchar().sauve_ImaPGM(p[2].p_fname.fname);
    return true;
  }
  static bool op_elagage(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imabin imadon((imadata<float>)ifich.LoadPGM(),1);
    imadon=imadon.elagage (p[1].p_int.val); 
    imadon.imaunsignedchar().sauve_ImaPGM(p[2].p_fname.fname);
    return true;
  }
  static bool op_zon_infl_geod(int nbp, Param_Type* p)
  {
    fichimage_entree ific1(p[0].p_fname.fname), ific2(p[1].p_fname.fname);
    imabin imamarq((imadata<float>)ific1.LoadPGM(),1), imadon((imadata<float>)ific2.LoadPGM(),1);
    int k=4;
    if (p[2].p_sel.sel_cur==1) k=8;
		bool icheck=1;
    imadon=imadon.zones_influence_geodesique(imamarq,icheck,k); 
		imabin imab=imadon.negatif(); 
		imab=imab.elagage(40);
		imadon=imab.negatif();
    imadon.imaunsignedchar().sauve_ImaPGM(p[3].p_fname.fname);
    return true;
  }
  static bool op_elim_objet_bord(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imabin imadon((imadata<float>)ifich.LoadPGM(),1);
    imabin imamarq(imadon.nlig(),imadon.ncol());
    int i,n;
    n=imamarq.ncol()-1;
    for (i=0; i<imamarq.nlig(); i++) {
      imamarq(i,0)=imamarq(i,n)=1;
    }
    n=imamarq.nlig()-1;
    for (i=0; i<imamarq.ncol(); i++) {
      imamarq(0,i)=imamarq(n,i)=1;
    }
    int k=4;
    if (p[1].p_sel.sel_cur==1) k=8;
    (imadon-imadon.reconstruction_geodesique(imamarq,k)).imaunsignedchar().sauve_ImaPGM(p[2].p_fname.fname);
    return true;
  }
  static bool op_bouche_trou(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imabin imadon((imadata<float>)ifich.LoadPGM(),1), imares(imadon.nlig(),imadon.ncol());
    int i,j,n, ncc, nmax=0, iccfond=0, k=4;
    if (p[1].p_sel.sel_cur==1) k=8;
		imadata<int> imacc=imadon.negatif().composantes_connexes(ncc,k); //imacc.imaunsignedchar().sauve_ImaPGM(p[2].p_fname.fname);
		int *T_n=new int[ncc];
		for (i=0; i<ncc; i++) T_n[i]=0;
    n=imacc.ncol()-1;
    for (i=0; i<imacc.nlig(); i++) {
      if (imacc(i,0)>0) T_n[imacc(i,0)-1]++;
      if (imacc(i,n)>0) T_n[imacc(i,n)-1]++;
    }
    n=imacc.nlig()-1;
    for (i=1; i<imacc.ncol()-1; i++) {
      if (imacc(0,i)>0) T_n[imacc(0,i)-1]++;
      if (imacc(n,i)>0) T_n[imacc(n,i)-1]++;
    }
		for (i=0; i<ncc; i++) 
			if (T_n[i]>nmax) {nmax=T_n[i]; iccfond=i+1;}; cout<<" composante connexe representant le fond = "<<iccfond<<"\n";
		if (T_n!=NULL) delete[] T_n;
		if (iccfond!=0) {
			imares.mise_a_zero();
			for (i=0; i<imacc.nlig(); i++) 
				for (j=0; j<imacc.ncol(); j++) imares(i,j)=(imacc(i,j)==iccfond?1:0);
		}
    imares.negatif().imaunsignedchar().sauve_ImaPGM(p[2].p_fname.fname);
    return true;
  }
  static bool op_erode_fct(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
    eltstruct S3(3,3); 
    if(p[2].p_sel.sel_cur==0) {S3(0,0)=0; S3(0,2)=0; S3(2,0)=0; S3(2,2)=0;}
    imadon=imadon.erode(S3, p[1].p_sel.sel_cur+1);
    imadon.imaunsignedchar().sauve_ImaPGM(p[3].p_fname.fname);
    return true;
  }
  static bool op_dilate_fct(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
    eltstruct S3(3,3); 
    if(p[2].p_sel.sel_cur==0) {S3(0,0)=0; S3(0,2)=0; S3(2,0)=0; S3(2,2)=0;}
    imadon=imadon.dilate(S3, p[1].p_sel.sel_cur+1);
    imadon.imaunsignedchar().sauve_ImaPGM(p[3].p_fname.fname);
    return true;
  }
  static bool op_ouvre_fct(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
    eltstruct S3(3,3); 
    if(p[2].p_sel.sel_cur==0) {S3(0,0)=0; S3(0,2)=0; S3(2,0)=0; S3(2,2)=0;}
    imadon=imadon.erode(S3,p[1].p_sel.sel_cur+1).dilate(S3,p[1].p_sel.sel_cur+1);
    imadon.imaunsignedchar().sauve_ImaPGM(p[3].p_fname.fname);
    return true;
  }
  static bool op_ferme_fct(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
    eltstruct S3(3,3); 
    if(p[2].p_sel.sel_cur==0) {S3(0,0)=0; S3(0,2)=0; S3(2,0)=0; S3(2,2)=0;}
    imadon=imadon.dilate(S3,p[1].p_sel.sel_cur+1).erode(S3,p[1].p_sel.sel_cur+1);
    imadon.imaunsignedchar().sauve_ImaPGM(p[3].p_fname.fname);
    return true;
  }
	static bool op_top_hat(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
    eltstruct S3(3,3); 
    if(p[2].p_sel.sel_cur==0) {S3(0,0)=0; S3(0,2)=0; S3(2,0)=0; S3(2,2)=0;}
		if (p[3].p_sel.sel_cur==0) 
			imadon=imadon-imadon.erode(S3,p[1].p_sel.sel_cur+1).dilate(S3,p[1].p_sel.sel_cur+1);
		else imadon=imadon.dilate(S3,p[1].p_sel.sel_cur+1).erode(S3,p[1].p_sel.sel_cur+1)-imadon;
		imadon.statbasic(1);
    imadon.imaunsignedchar().sauve_ImaPGM(p[4].p_fname.fname);
    return true;
  }
	static bool op_top_hat_2(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
    eltstruct S3(3,3); 
    if(p[2].p_sel.sel_cur==0) {S3(0,0)=0; S3(0,2)=0; S3(2,0)=0; S3(2,2)=0;}
//		if (p[3].p_sel.sel_cur==0) 
			imadon=imadon-imadon.reconst_geod(imadon.erode(S3,p[1].p_sel.sel_cur+1),S3);
//		else imadon=imadon.dilate(S3,p[1].p_sel.sel_cur+1).erode(S3,p[1].p_sel.sel_cur+1)-imadon;
		imadon.statbasic(1);
    imadon.imaunsignedchar().sauve_ImaPGM(p[4].p_fname.fname);
    return true;
  }
  static bool op_reconst_geod_fct(int nbp, Param_Type* p)
  {
    fichimage_entree ific1(p[0].p_fname.fname), ific2(p[1].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ific1.LoadPGM(), imamarq=(imadata<float>)ific2.LoadPGM();
    eltstruct S3(3,3); 
    if(p[2].p_sel.sel_cur==0) {S3(0,0)=0; S3(0,2)=0; S3(2,0)=0; S3(2,2)=0;}
    imadon=imadon.reconst_geod(imamarq,S3);
    imadon.imaunsignedchar().sauve_ImaPGM(p[3].p_fname.fname);
    return true;
  }
	static bool op_h_max_min(int nbp, Param_Type* p) {
    fichimage_entree ific(p[0].p_fname.fname);
		int isens=1; if (p[1].p_sel.sel_cur==1) isens=-1;
		float h=p[2].p_float.val;
    imadata<float> imadon=((imadata<float>)ific.LoadPGM())*(float)isens, imamarq=imadon+(-h);
    eltstruct S3(3,3); 
    if(p[3].p_sel.sel_cur==0) {S3(0,0)=0; S3(0,2)=0; S3(2,0)=0; S3(2,2)=0;}
    imadon=(imadon.reconst_geod(imamarq,S3))*(float)isens;
    imadon.imaunsignedchar().sauve_ImaPGM(p[4].p_fname.fname);
    return true;
	}
	static bool op_h_max(int nbp, Param_Type* p) {
    fichimage_entree ific(p[0].p_fname.fname);
		float h=p[1].p_float.val;
    imadata<float> imadon=(imadata<float>)ific.LoadPGM(), imamarq=imadon+(-h);
    eltstruct S3(3,3); 
    if(p[2].p_sel.sel_cur==0) {S3(0,0)=0; S3(0,2)=0; S3(2,0)=0; S3(2,2)=0;}
    imadon=imadon.reconst_geod(imamarq,S3);
    imadon.imaunsignedchar().sauve_ImaPGM(p[3].p_fname.fname);
    return true;
	}
	static bool op_h_min(int nbp, Param_Type* p) {
    fichimage_entree ific(p[0].p_fname.fname);
		float h=p[1].p_float.val;
    imadata<float> imadon=((imadata<float>)ific.LoadPGM())*(-1), imamarq=imadon+(-h);
    eltstruct S3(3,3); 
    if(p[2].p_sel.sel_cur==0) {S3(0,0)=0; S3(0,2)=0; S3(2,0)=0; S3(2,2)=0;}
    imadon=(imadon.reconst_geod(imamarq,S3))*(-1);
    imadon.imaunsignedchar().sauve_ImaPGM(p[3].p_fname.fname);
    return true;
	}
	static bool op_max_reg(int nbp, Param_Type* p) {
    fichimage_entree ific(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ific.LoadPGM();
		float h=p[1].p_float.val;
    eltstruct S3(3,3); 
    if(p[2].p_sel.sel_cur==0) {S3(0,0)=0; S3(0,2)=0; S3(2,0)=0; S3(2,2)=0;}
    imadon=imadon.marq_max_regionaux(S3,h);
    imadon.imaunsignedchar().sauve_ImaPGM(p[3].p_fname.fname);
    return true;
	}
	static bool op_grad_mm(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
    eltstruct ES(3,3); 
    if(p[1].p_sel.sel_cur==0) {ES(0,0)=0; ES(0,2)=0; ES(2,0)=0; ES(2,2)=0;}
    imadon=imadon.gradient_m(ES); 
    imadon.imaunsignedchar().sauve_ImaPGM(p[2].p_fname.fname);
    return true;
  }
  static bool op_grad_multiechmm(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
    eltstruct ES(3,3); 
    if(p[2].p_sel.sel_cur==0) {ES(0,0)=0; ES(0,2)=0; ES(2,0)=0; ES(2,2)=0;}
		int ordre=p[1].p_sel.sel_cur+1;
		imadon=imadon.dilate(ES,ordre)-imadon.erode(ES,ordre);
		imadon=(imadon-imadon.erode(ES,ordre).dilate(ES,ordre)).erode(ES,ordre-1);
    imadon.imaunsignedchar().sauve_ImaPGM(p[3].p_fname.fname);
    return true;
  }
  static bool op_lapl_mm(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
    eltstruct ES(3,3); 
    if(p[1].p_sel.sel_cur==0) {ES(0,0)=0; ES(0,2)=0; ES(2,0)=0; ES(2,2)=0;}
    imadon=imadon.laplacien_m(ES); 
    imadon.imaunsignedchar().sauve_ImaPGM(p[2].p_fname.fname);
    return true;
  }
  static bool op_class_ppv(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
    ifstream fechant;
    fechant.open(p[1].p_fname.fname,ios::in);
    if (!fechant) {
      cout<<" ouverture de "<<fechant<<" impossible\n";
      return false;
    }
    else {
      sampleset<pixel<float> > datlearn;
      int k,ncl=0,i,n=0,l,j;
      fechant>>n;
      for (i=0; i<n; i++) {
        fechant>>l>>j>>k;
        if (l>=0 && j>=0 && l<imadon.nlig() && j<imadon.ncol()) datlearn.ajoute(imadon.pix(l,j),i,k);
        cout<<" coordonnees pixel (lig,col) de l'echantillon "<<i<<" : "<<" classe "<<k<<"\n";
        if (k>ncl) ncl=k;
      }
      fechant.close();
      datlearn.affiche();
      sampleset<pixel<float> > datset=imadon.dataset();
      datset.k_ppv (datlearn,(p[2].p_sel.sel_cur)*2+1);
      int nblig=imadon.nlig(), nbcol=imadon.ncol();
      imalabels imacl(nblig,nbcol);
      for (i=0; i<nblig; i++)
        for (j=0; j<nbcol; j++)
          imacl(i,j)=datset.remove(i*nbcol+j).label();
      imacl.conv2imBYTE().sauve_ImaPGM(p[3].p_fname.fname);
      return true;
    }
  }
  static bool op_class_cmeans(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
    pixel<float> pix;
    sample<pixel<float> > s_pix;
    sampleset<pixel<float> > datset=imadon.dataset();
		int ncl=p[1].p_int.val;
    statclasses clkmeans=datset.k_means(ncl);
    clkmeans.affiche();
    imalabels imacl(imadon,clkmeans,0.,"ICM");
    imacl.conv2imBYTE().sauve_ImaPGM(p[2].p_fname.fname);
		ofstream fclasses; fclasses.open("classes.txt",ios::out);
    if (!fclasses) {
      cout<<" ouverture de "<<fclasses<<" impossible\n"; return false;
    } else {
      int d=imadon.ncanaux(),k,j;
      stat1class cl(1,d); 
			fclasses<<ncl<<" "<<d<<"\n";
      for (k=0; k<ncl; k++) {
        cl=clkmeans.extract(k);
        for (j=0; j<d; j++) fclasses<<cl.mean(j);
        fclasses<<"\n";
        for (j=0; j<d; j++) fclasses<<cl.cova(j,j);
        fclasses<<"\n";
      }
      fclasses.close();
    }
    return true;
  }
  static bool op_class_cmeans_mask(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname), ificm(p[1].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
    imabin imaROI((imadata<float>)ificm.LoadPGM(),1);
    pixel<float> pix;
    sample<pixel<float> > s_pix;
    sampleset<pixel<float> > datset=imadon.dataset(imaROI);
    statclasses clkmeans=datset.k_means(p[2].p_int.val);
    clkmeans.affiche();
    imalabels imacl(imadon,clkmeans,0.,"ICM");
    int i,j,nl=imadon.nlig(),nc=imadon.ncol();
    for (i=0; i<nl; i++)
      for (j=0; j<nc; j++) if (!imaROI(i,j)) imacl(i,j)=0;
    imacl.conv2imBYTE().sauve_ImaPGM(p[3].p_fname.fname);
    return true;
  }
  static bool op_class_mrf(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
    ifstream fclasses;
    fclasses.open(p[1].p_fname.fname,ios::in);
    if (!fclasses) {
      cout<<" ouverture de "<<fclasses<<" impossible\n";
      return false;
    }
    else {
      int d=imadon.ncanaux(),k,ncl,j;
      statclasses Tcl;
      stat1class cl(1,d);
      fclasses>>ncl>>d;
      for (k=0; k<ncl; k++) {
        cl.label()=k+1; 
        for (j=0; j<d; j++) fclasses>>cl.mean(j);
        for (j=0; j<d; j++) fclasses>>cl.cova(j,j);
        for (j=0; j<d; j++) cl.icova(j,j)=1.f/cl.cova(j,j);
        cl.dcova()=1.; for (j=0; j<d; j++) cl.dcova()*=cl.cova(j,j);
        Tcl.ajoute(cl);
      }
      fclasses.close();
      Tcl.affiche();
      char* algo="ICM";
      if (p[3].p_sel.sel_cur==1) algo="SAG";
      if (p[3].p_sel.sel_cur==2) algo="SAM";
      imalabels imacl(imadon,Tcl,p[2].p_float.val,algo);
//      imacl.conv2imBYTE(1).sauve_ImaPGM(p[4].p_fname.fname);
      imacl.conv2imBYTE(0).sauve_ImaPGM("./classif.pgm");
			imacl.conv2imFACOL(0).sauve_ImaPGM(p[4].p_fname.fname);
      return true;
    }
  }
  static bool op_class_mrf_lines(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
    ifstream fclasses;
    fclasses.open(p[1].p_fname.fname,ios::in);
    if (!fclasses) {
      cout<<" ouverture de "<<fclasses<<" impossible\n";
      return false;
    }
    else {
      int d=imadon.ncanaux(),k,ncl,j;
      statclasses Tcl;
      stat1class cl(1,d);
      fclasses>>ncl>>d;
      for (k=0; k<ncl; k++) {
        cl.label()=k+1; 
        for (j=0; j<d; j++) fclasses>>cl.mean(j);
        for (j=0; j<d; j++) fclasses>>cl.cova(j,j);
        for (j=0; j<d; j++) cl.icova(j,j)=1.f/cl.cova(j,j);
        cl.dcova()=1.; for (j=0; j<d; j++) cl.dcova()*=cl.cova(j,j);
        Tcl.ajoute(cl);
      }
      fclasses.close();
      Tcl.affiche();
      char* algo="ICM";
      if (p[3].p_sel.sel_cur==1) algo="SAG";
      if (p[3].p_sel.sel_cur==2) algo="SAM";
//      imalabels imacl(imadon,Tcl,p[2].p_float.val,algo);
			float a0=0.9f, a1=1.8f, a2=2.7f;
			imalabels imacl(imadon,Tcl,'1',p[2].p_float.val,a0,a1,a2);
      imacl.conv2imBYTE(1).sauve_ImaPGM(p[4].p_fname.fname);
      return true;
    }
  }
	static bool op_class_emgibbs(int nbp, Param_Type* p)
  {
    return true;
  }
  static bool op_bord_class(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<BYTE> imalab=ifich.LoadPGM();
		int nblig=imalab.nlig(), nbcol=imalab.ncol(), i0,i2,j0,j2,l;
		imabin imabords(nblig,nbcol);
		if(p[1].p_sel.sel_cur==0) { // 4-connexite
			for (int i=0; i<nblig; i++) {
				i0=maxi(0,i-1); i2=mini(i+1,nblig-1);
				for (int j=0; j<nbcol; j++) {
					j0=maxi(0,j-1); j2=mini(j+1,nbcol-1); l=imalab(i,j);
					imabords(i,j)=(l!=imalab(i0,j)||l!=imalab(i2,j)||l!=imalab(i,j0)||l!=imalab(i,j2));
				}
			}
		} else { // 8-connexite
			for (int i=0; i<nblig; i++) {
				i0=maxi(0,i-1); i2=mini(i+1,nblig-1);
				for (int j=0; j<nbcol; j++) {
					j0=maxi(0,j-1); j2=mini(j+1,nbcol-1); l=imalab(i,j);
					imabords(i,j)=(l!=imalab(i0,j)||l!=imalab(i2,j)||l!=imalab(i,j0)||l!=imalab(i,j2)||
												 l!=imalab(i0,j0)||l!=imalab(i0,j2)||l!=imalab(i2,j0)||l!=imalab(i2,j2));
				}
			}
		}
		cout<<" longueur des bords des classes = "<<imabords.norm()<<"\n";
    imabords.imaunsignedchar().sauve_ImaPGM(p[2].p_fname.fname);
    return true;
  }
  static bool op_image_gradient(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
    char *masq;
    switch (p[1].p_sel.sel_cur) {
      case 0:  masq="Prewitt"; break;
      case 2:  masq="MDIF"; break;
      default: masq="Sobel";
    }
    imadata<float> imadir(imadon);
    imadata<float> imagrd=imadon.gradient(imadir,masq);
    imagrd.sauve_ImaPGM(p[2].p_fname.fname);
    /*	float coef=90/(float)PI;
    imadir=(imadir*coef)+90.;
    cout<<" les valeurs de l'image de la direction du gradient sont 1<->2deg.\n";
    imadir.sauve_ImaPGM(p[3].p_fname.fname);*/
    return true;
  }
  static bool op_cont_gradient(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
    char *masq;
    switch (p[1].p_sel.sel_cur) {
      case 0:  masq="Prewitt"; break;
      case 2:  masq="MDIF"; break;
      default: masq="Sobel";
    }
    imacontours imacont(imadon,masq,p[2].p_int.val,1);
    imacont.conv2imBYTE().sauve_ImaPGM(p[3].p_fname.fname);
    return true;
  }
  static bool op_cont_laplacien(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
    char *masqG, *masqL;
    switch (p[1].p_sel.sel_cur) {
      case 0:  masqG="Prewitt"; break;
      case 2:  masqG="MDIF"; break;
      default: masqG="Sobel";
    }
    switch (p[2].p_sel.sel_cur) {
      case 0:  masqL="4connex"; break;
      default: masqL="8connex";
    }
    imacontours imacont(imadon,masqL,masqG,1);
    imacont.conv2imBYTE().sauve_ImaPGM(p[3].p_fname.fname);
    return true;
  }
  static bool op_cont_optimal(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
    char *filtre;
    switch (p[1].p_sel.sel_cur) {
      case 0: filtre="Deriche"; break;
      case 1: filtre="Shen"; break;
    }
    imacontours imacont(imadon,p[2].p_float.val,filtre,1);
    imacont.conv2imBYTE().sauve_ImaPGM(p[3].p_fname.fname);
    return true;
  }
  static bool op_affine_contours(int nbp, Param_Type* p)
  {
    fichimage_entree ific1(p[0].p_fname.fname), ific2(p[1].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ific2.LoadPGM();
    char *masq="Sobel";
    imadata<float> imadir(imadon);
    imadata<float> imagrd=imadon.gradient(imadir,masq);
    imacontours imacont((imadata<BYTE>)ific1.LoadPGM());
    imacont.imagrd_maxloc(imagrd,imadir);
    imacont.conv2imBYTE().sauve_ImaPGM(p[2].p_fname.fname);
    return true;
  }
  static bool op_prolonge_contours(int nbp, Param_Type* p)
  {
    fichimage_entree ific1(p[0].p_fname.fname), ific2(p[1].p_fname.fname);
    imadata<float> imagrd=(imadata<float>)ific2.LoadPGM();
    imacontours imacont((imadata<BYTE>)ific1.LoadPGM());
//    imacont.prolongecontours(imagrd,p[2].p_int.val,p[3].p_float.val);
    imacont.edgedrawing(imagrd,p[2].p_int.val,p[4].p_int.val,p[3].p_float.val,1);
    imacont.conv2imBYTE().sauve_ImaPGM(p[5].p_fname.fname);
    return true;
  }
  static bool op_edge_drawing(int nbp, Param_Type* p)
  {
    fichimage_entree ific1(p[0].p_fname.fname), ific2(p[1].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ific2.LoadPGM();
    imacontours imacont((imadata<BYTE>)ific1.LoadPGM());
    imacont.edgedrawing(imadon,p[2].p_int.val,p[3].p_int.val);
//    imacont.edgedrawing(imadon,p[2].p_int.val,p[3].p_int.val,1);
    imacont.conv2imBYTE().sauve_ImaPGM(p[4].p_fname.fname);
    return true;
  }
/*  static bool op_test_2_2(int nbp, Param_Type* p)
  {
    cout<<"Entree #1: "<<p[0].p_fname.fname<<'\n';
    cout<<"Entree #2: "<<p[1].p_fname.fname<<'\n';
    cout<<"Sortie #1: "<<p[2].p_fname.fname<<'\n';
    cout<<"Sortie #2: "<<p[3].p_fname.fname<<'\n';
    return true;
  }*/
  static bool op_hough(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imabin imadon((imadata<float>)ifich.LoadPGM(),1);
    imadata<int> imhough=imadon.transformee_Hough();
//  imadata<int> imhough=imadon.hough_transform();
    imhough.statbasic(1);
    imhough.sauve_ImaPGM(p[1].p_fname.fname);
    return true;
  }
  static bool op_rec_hough(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname), ifich2(p[1].p_fname.fname);
    imadata<float> imaTHough=(imadata<float>)ifich.LoadPGM(), ima=(imadata<float>)ifich2.LoadPGM();
		int nblig=ima.nlig(), nbcol=ima.ncol();
    imabin im_rec; im_rec.reconst_hough_transform(imaTHough,nblig,nbcol);
	  im_rec.imaunsignedchar().sauve_ImaPGM(p[2].p_fname.fname);
    return true;
  }
  static bool op_profil_bin(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imabin imab((imadata<float>)ifich.LoadPGM(),1);
    int i,j,nblig=imab.nlig(),nbcol=imab.ncol(),n=16; 
    if (p[1].p_sel.sel_cur==0) {
      imadata<float> imadon(n,nbcol);
      for (j=0; j<nbcol; j++) {
        for (i=0; i<nblig; i++) imadon(0,j)+=imab(i,j);
        for (i=1; i<n; i++) imadon(i,j)=imadon(0,j);
      }
      imadon.sauve_ImaPGM(p[2].p_fname.fname);
    } else {
      imadata<float> imadon(n,nblig);
      for (i=0; i<nblig; i++) {
        for (j=0; j<nbcol; j++) imadon(0,i)+=imab(i,j);
        for (j=1; j<n; j++) imadon(j,i)=imadon(0,i);
      }
      imadon.sauve_ImaPGM(p[2].p_fname.fname);
    }
    return true;
  }
  static bool op_analyse_hough(int nbp, Param_Type* p)
  {
    fichimage_entree ific1(p[0].p_fname.fname); 
    imabin imadon((imadata<float>)ific1.LoadPGM(),1);
    imadata<int> imhough=imadon.hough_transform();
//	imadata<int> imsegments=imadon.analyse_hough(imhough);
	//imhough.statbasic(1);
//	imsegments.sauve_ImaPGM(p[1].p_fname.fname);
    return true;
  }
  static bool op_hough_cercles(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imabin imadon((imadata<float>)ifich.LoadPGM(),1);
    imadata<int> imhough=imadon.transformee_Hough_cercles();
    imhough.statbasic(1);
		float xx=FLT_MIN;
		for (int k=0; k<10; k++) {
			imadata<int> imahough_r(imhough,k);
			imahough_r.statbasic(1);
//			if (imahough_r.maxI()>=2*PI*k*10) {
			if (imahough_r.maxI()>xx) {
				cout<<" max a "<<imahough_r.maxI()<<" pour R="<<10*k<<"\n"; 
				imahough_r.sauve_ImaPGM(p[1].p_fname.fname);
				xx=(float)imahough_r.maxI();
			}
		}
    return true;
  }
  static bool op_points_interet(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ifich.LoadPGM();
		imadata<BYTE> imapointinteret(imadon.nlig(),imadon.ncol());
		keypoints(imadon,p[1].p_float.val,p[2].p_float.val,(bool)1,(bool)1).conv2imBYTE(imapointinteret); 
		imapointinteret.sauve_ImaPGM(p[3].p_fname.fname); 
    return true;
  }
  static bool op_corresp_points_interet(int nbp, Param_Type* p)
  {
    fichimage_entree ific1(p[0].p_fname.fname), ific2(p[1].p_fname.fname);
    imadata<BYTE> imakpts1=ific1.LoadPGM(), imakpts2=ific2.LoadPGM();
		float dx,dy,rapport,angle;
		imadata<float> imaFM=imakpts1.transformee_FourierMellin(imakpts2,dx,dy,rapport,angle,1);
//		imaFM.imaunsignedchar(1).sauve_ImaPGM(p[2].p_fname.fname);
		imakpts2.projette_ima(dx,dy,rapport,angle).imaunsignedchar(0).sauve_ImaPGM(p[2].p_fname.fname);
    return true;
  }
  static bool op_segm_classif(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imalabels imalab((imadata<BYTE>)ifich.LoadPGM());
    imaregions imareg(imalab);
    if (imareg.nregions()<256) imareg.conv2imBYTE().sauve_ImaPGM(p[1].p_fname.fname);
    return true;
  }
  static bool op_croise_2segment(int nbp, Param_Type* p)
  {
    fichimage_entree ifich0(p[0].p_fname.fname), ifich1(p[1].p_fname.fname);
    imaregions imareg1((imadata<BYTE>)ifich0.LoadPGM()), imareg2((imadata<BYTE>)ifich1.LoadPGM());
    imaregions imareg(imareg1,imareg2);
    if (imareg.nregions()<256) imareg.conv2imBYTE().sauve_ImaPGM(p[2].p_fname.fname);
		else if (imareg.nregions()<65536) imareg.conv2imUI().sauve_ImaPGM2(p[2].p_fname.fname);
    return true;
  }	
	static bool op_reg_growing(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon((imadata<float>)ifich.LoadPGM());
    char* select="histo";
    if (p[2].p_sel.sel_cur==1) select="aleat";
    imaregions imareg(imadon,p[1].p_float.val,select,1);
		imareg.suppression_trop_petites_regions(5);
		imareg.suppression_plus_petites_regions (255);
    imareg.conv2imBYTE().sauve_ImaPGM(p[3].p_fname.fname);
    return true;
  }
  static bool op_reg_quadtree(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon((imadata<float>)ifich.LoadPGM());
    imaregions imareg(imadon,p[1].p_float.val,p[2].p_int.val,1);
		imareg.suppression_trop_petites_regions(5);
		imareg.suppression_plus_petites_regions (255);
    imareg.conv2imBYTE().sauve_ImaPGM(p[3].p_fname.fname);
    return true;
  }
  static bool op_fusion_reg(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon((imadata<float>)ifich.LoadPGM());
    imaregions imareg(imadon,(double)p[1].p_float.val);
		imareg.suppression_trop_petites_regions(5);
		imareg.suppression_plus_petites_regions (255);
    imareg.conv2imBYTE().sauve_ImaPGM(p[2].p_fname.fname);
    return true;
  }
	static bool op_reg_lpe(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon((imadata<float>)ifich.LoadPGM());
    imaregions imareg(imadon);
		int n=(unsigned int)p[1].p_int.val;
		imareg.suppression_trop_petites_regions((int)n);
		cout<<" nombre de regions de taille >= "<<n<<" : "<<imareg.nregions()<<"\n";
/*		while (imareg.nregions()>255) {
			n++;
			imareg.suppression_trop_petites_regions(n);
			cout<<" nombre de regions de taille >= "<<n<<" : "<<imareg.nregions()<<"\n";
		}*/
//    imareg.conv2imBYTE().sauve_ImaPGM(p[2].p_fname.fname);
    if (imareg.nregions()<256) imareg.conv2imBYTE().sauve_ImaPGM(p[2].p_fname.fname);
		else if (imareg.nregions()<65536) imareg.conv2imUI().sauve_ImaPGM2(p[2].p_fname.fname);
		int i,j,nl=imadon.nlig(),nc=imadon.ncol();
		imabin imacont(nl,nc); imacont.mise_a_zero();
		for (i=0; i<nl; i++)
			for (j=0; j<nc; j++) if ((imareg(i,j)==0) || (i>0 && imareg(i,j)!=imareg(i-1,j)) || (j>0 && imareg(i,j)!=imareg(i,j-1))) imacont(i,j)=1;
		for (i=0; i<nl; i++) imacont(i,0)=imacont(i,nc-1)=1;
		for (j=0; j<nc; j++) imacont(0,j)=imacont(nl-1,j)=1;
		imacont.imaunsignedchar().sauve_ImaPGM("imaCont_LPE.pgm");
		imacont=imacont.elagage(40);
		imacont.imaunsignedchar().sauve_ImaPGM("imaCont_LPE_eb40.pgm");
		imacouleur<BYTE> imacol(imareg.conv2imBYTE(0,0),0);
		char nomfic[100];
		strncpy(nomfic,p[2].p_fname.fname,strlen(p[2].p_fname.fname)+1); cout<<nomfic;
		strncpy(nomfic+strlen(p[2].p_fname.fname)-4,".ppm",strlen(".ppm")+1); cout<<nomfic;
		imacol.sauve_ImaPGM(nomfic);
    return true;
  }
  static bool op_reg_graphe(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon((imadata<float>)ifich.LoadPGM());
    imaregions imareg(imadon,(unsigned int)p[1].p_int.val);
    imareg.conv2imBYTE().sauve_ImaPGM(p[2].p_fname.fname);
    return true;
  }
  static bool op_reg_mumford(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon((imadata<float>)ifich.LoadPGM());
    imaregions imareg(imadon,(unsigned int)p[1].p_int.val,p[2].p_fname.fname);
    imareg.conv2imBYTE().sauve_ImaPGM(p[2].p_fname.fname);
    return true;
  }
/*   static bool op_fus_reg_mumford(int nbp, Param_Type* p)
  {
    fichimage_entree ifich0(p[0].p_fname.fname), ifich1(p[1].p_fname.fname);
    imadata<float> imadon((imadata<float>)ifich0.LoadPGM());
    imadata<int> imareginit(ifich1.LoadPGM2());
    imaregions imareg(imadon,imareginit,(unsigned int)p[2].p_int.val,p[3].p_fname.fname);
//    imareg.conv2imBYTE().sauve_ImaPGM(p[3].p_fname.fname);
    if (imareg.nregions()<256) imareg.conv2imBYTE().sauve_ImaPGM(p[3].p_fname.fname);
		else if (imareg.nregions()<65536) imareg.conv2imUI().sauve_ImaPGM2(p[3].p_fname.fname);
    return true;
  }*/
		static bool op_suppix_SLIC(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon((imadata<float>)ifich.LoadPGM());
		long int nreg=(long int)p[1].p_int.val;
		double fact_m=p[2].p_float.val;
    imaregions imareg(imadon,nreg,fact_m);
		float s=(float)imareg.ncol()*imareg.nlig()/nreg;
		const float fact_szreg_aborbe=1/20.f, fact_szreg_ok=1/4.f;
		imareg.reconnecte_reg((int)(s*fact_szreg_ok),(int)(s*fact_szreg_aborbe));
/*		imareg.suppression_trop_petites_regions((int)n);
		cout<<" nombre de regions de taille >= "<<n<<" : "<<imareg.nregions()<<"\n";*/
		cout<<" image superpixels avec "<<imareg.nregions()<<" regions\n";
		if (imareg.nregions()<256) {cout<<" **** "<<p[3].p_fname.fname; imareg.conv2imBYTE().sauve_ImaPGM(p[3].p_fname.fname);}
		else if (imareg.nregions()<65536) imareg.conv2imUI().sauve_ImaPGM2(p[3].p_fname.fname);
		int i,j,nl=imadon.nlig(),nc=imadon.ncol();
		imabin imacont(nl,nc); imacont.mise_a_zero();
		for (i=0; i<nl; i++)
			for (j=0; j<nc; j++) if ((imareg(i,j)==0) || (i>0 && imareg(i,j)!=imareg(i-1,j)) || (j>0 && imareg(i,j)!=imareg(i,j-1))) imacont(i,j)=1;
		for (i=0; i<nl; i++) imacont(i,0)=imacont(i,nc-1)=1;
		for (j=0; j<nc; j++) imacont(0,j)=imacont(nl-1,j)=1;
		imacont.imaunsignedchar().sauve_ImaPGM("imaCont_SLIC.pgm");
/*		imacont=imacont.elagage(40);
		imacont.imaunsignedchar().sauve_ImaPGM("imaCont_LPE_eb40.pgm");*/
    return true;
  }
	static bool op_waterpix(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon((imadata<float>)ifich.LoadPGM());
		long int nreg=(long int)p[1].p_int.val;
		double param_k=p[2].p_float.val;
		int k=4;
    if (p[3].p_sel.sel_cur==1) k=8;
		imaregions imareg(imadon,param_k,nreg,k); cout<<" # superpiels = "<<imareg.nregions()<<" "<<(unsigned int)imareg.nregions()<<"\n";
/*		float s=(float)imareg.ncol()*imareg.nlig()/nreg;
		const float fact_szreg_aborbe=1/20.f, fact_szreg_ok=1/4.f;
		imareg.reconnecte_reg(s*fact_szreg_ok,s*fact_szreg_aborbe);*/
		if ((unsigned int)imareg.nregions()<256) {cout<<imareg.nregions()<<" **** "<<p[4].p_fname.fname; imareg.conv2imBYTE().sauve_ImaPGM(p[4].p_fname.fname);}
		else //if ((unsigned int)imareg.nregions()<65536) 
			{cout<<imareg.nregions()<<" $$$$ "<<p[4].p_fname.fname; imareg.conv2imUI().sauve_ImaPGM2(p[4].p_fname.fname);}
		int i,j,nl=imadon.nlig(),nc=imadon.ncol();
		imabin imacont(nl,nc); imacont.mise_a_zero();
		for (i=0; i<nl; i++)
			for (j=0; j<nc; j++) if ((imareg(i,j)==0) || (i>0 && imareg(i,j)!=imareg(i-1,j)) || (j>0 && imareg(i,j)!=imareg(i,j-1))) imacont(i,j)=1;
		for (i=0; i<nl; i++) imacont(i,0)=imacont(i,nc-1)=1;
		for (j=0; j<nc; j++) imacont(0,j)=imacont(nl-1,j)=1;
		imacont.imaunsignedchar().sauve_ImaPGM("imaCont_waterpix.pgm");
    return true;
  }
	 static bool op_param_objets(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imaobjets ima_obj((imadata<BYTE>)ifich.LoadPGM());
    ima_obj.basic_param();
    ima_obj.affiche();
    return true;
  }
  static bool op_detect_disques(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imaobjets imaobj((imadata<BYTE>)ifich.LoadPGM());
    imaobj.detection_disques(p[1].p_float.val,p[2].p_float.val).imaunsignedchar().sauve_ImaPGM(p[3].p_fname.fname);
    return true;
  }
	static bool op_class_tailleCC(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imabin imadon((imadata<float>)ifich.LoadPGM(),1);
    int ncc=0, k=8;
    imadata<int> imacc=imadon.composantes_connexes (ncc, k);
    imaobjets imaobj((imadata<BYTE>)imacc);
    imaobj.basic_param(0); imaobj.affiche();
		int nbO=imaobj.nb_obj(); float *T_n=new float[nbO];
    sample<pixel<float> > s_pix;
		sampleset<pixel<float> > datset;
		for (int k=1; k<nbO; k++) {
			T_n[k]=imaobj.objet_i(k).surf_obj(); 
			s_pix=pixel<float>(1,(int)T_n[k]); datset.ajoute(s_pix);} 
		datset.affiche();
    statclasses clkmeans=datset.k_means(2); clkmeans.affiche();
		float seuil=(clkmeans.extract(1).mean(0)+clkmeans.extract(2).mean(0))*0.5f; cout<<" seuil de classification = "<<seuil<<"\n";
		BYTE *T_l=new BYTE[nbO]; T_l[0]=0;
		for (int k=1; k<nbO; k++) if (T_n[k]<=seuil) T_l[k]=1; else T_l[k]=2;
		int nblig=imacc.nlig(), nbcol=imacc.ncol();
		imadata<BYTE> imalab(nblig,nbcol);
		for (int i=0; i<nblig; i++) for (int j=0; j<nbcol; j++) imalab(i,j)=T_l[imacc(i,j)];
		imalab.sauve_ImaPGM(p[1].p_fname.fname);
		if (T_l!=NULL) delete[] T_l; if (T_n!=NULL) delete[] T_n;
    return true;
  }
  static bool op_conv_couleur(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imacouleur<BYTE> imacol((imadata<BYTE>)ifich.LoadPGM());
    switch (p[1].p_sel.sel_cur) {
      case 0: imacol.RVB2IST();
              imacol.sauve_ImaPGM(p[2].p_fname.fname);break;
      case 1: imacol.RVB2Lab();
              imacol.sauve_ImaPGM(p[2].p_fname.fname);break;
      case 2: imacol.RVB2I1I2I3();
              imacol.sauve_ImaPGM(p[2].p_fname.fname); break;
      case 3: imacol.RVB2NormeL1();
              imacol.sauve_ImaPGM(p[2].p_fname.fname); break;
    }
    return true;
  }
  static bool op_conv_false_color(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<BYTE> ima(ifich.LoadPGM());
		imacouleur<BYTE> imacol(ima,0);
		imacol.sauve_ImaPGM(p[1].p_fname.fname);
    return true;
  }
  static bool op_split_components(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imacouleur<BYTE> imacol((imadata<BYTE>)ifich.LoadPGM());
    imadata<BYTE> im0=imacol.select_composante(0);
    imadata<BYTE> im1=imacol.select_composante(1);
    imadata<BYTE> im2=imacol.select_composante(2);
    im0.sauve_ImaPGM(p[1].p_fname.fname);
    im1.sauve_ImaPGM(p[2].p_fname.fname);
    im2.sauve_ImaPGM(p[3].p_fname.fname);
    return true;
  }
  static bool op_superpos_can(int nbp, Param_Type* p)
  {
    fichimage_entree ific0(p[0].p_fname.fname), ific1(p[1].p_fname.fname), ific2(p[2].p_fname.fname);
		imadata<BYTE> ima0(ific0.LoadPGM());
    imacouleur<BYTE> imacol(ima0.nlig(),ima0.ncol());
		imacol.copiecanal(0,ima0); imacol.copiecanal(1,ific1.LoadPGM()); imacol.copiecanal(2,ific2.LoadPGM());
    imacol.sauve_ImaPGM(p[3].p_fname.fname);
    return true;
  }
/*  static bool op_contraste(int nbp, Param_Type* p)
  {
	fichimage_entree ifich(p[0].p_fname.fname);
    imadata<BYTE> ima((imadata<BYTE>)ifich.LoadPGM());
//	imadata<float> imares=ima.ima_contrast(p[1].p_int.val,p[2].p_int.val);
//	imares.sauve_ImaPGM(p[3].p_fname.fname);
	imadata<float> imares=ima.ima_contrast(p[1].p_int.val);
	imares.sauve_ImaPGM(p[2].p_fname.fname);
   return true;
  }
  static bool op_entropy(int nbp, Param_Type* p)
  {
	fichimage_entree ifich(p[0].p_fname.fname);
    imadata<BYTE> ima((imadata<BYTE>)ifich.LoadPGM());
//	imadata<float> imares=ima.ima_entropy(p[1].p_int.val,p[2].p_int.val);
//	imares.sauve_ImaPGM(p[3].p_fname.fname);
	imadata<float> imares=ima.ima_entropy(p[1].p_int.val);
	imares.imaunsignedchar(1).sauve_ImaPGM(p[2].p_fname.fname);
   return true;
  }
  static bool op_energy(int nbp, Param_Type* p)
  {
	fichimage_entree ifich(p[0].p_fname.fname);
    imadata<BYTE> ima((imadata<BYTE>)ifich.LoadPGM());
//    imadata<float> imares=ima.ima_energy(p[1].p_int.val,p[2].p_int.val);
//	imares.sauve_ImaPGM(p[3].p_fname.fname);
    imadata<float> imares=ima.ima_energy(p[1].p_int.val);
	imares.imaunsignedchar(1).sauve_ImaPGM(p[2].p_fname.fname);
   return true;
  }*/
  static bool op_texture1(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<BYTE> ima((imadata<BYTE>)ifich.LoadPGM());
    imadata<float> imares(ima);
    switch (p[2].p_sel.sel_cur) {
       case 0: imares=ima.ima_contrast(p[1].p_int.val); break;
       case 1: imares=ima.ima_entropy(p[1].p_int.val); break;
       case 2: imares=ima.ima_energy(p[1].p_int.val); break;
       default: break;
    }
		imares.statbasic(1);
		if (p[2].p_sel.sel_cur==1) {
			cout<<" multiplication des valeurs d'entropie par 256/log2(64) = 4\n";
			imares=imares*4.f; imares.statbasic(1);
			imares.imaunsignedchar(0).sauve_ImaPGM(p[3].p_fname.fname);
		}
    else imares.imaunsignedchar(1).sauve_ImaPGM(p[3].p_fname.fname);
    return true;
  }
  static bool op_moment(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<BYTE> ima((imadata<BYTE>)ifich.LoadPGM());
//    imadata<float> imares=ima.ima_moment(p[3].p_int.val,p[1].p_int.val,p[2].p_int.val);
//	imares.sauve_ImaPGM(p[4].p_fname.fname);
    imadata<float> imares=ima.ima_moment(p[2].p_int.val,p[1].p_int.val);
    imares.imaunsignedchar(1).sauve_ImaPGM(p[3].p_fname.fname);
    return true;
  }
  static bool op_cooc_gris(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<BYTE> ima((imadata<BYTE>)ifich.LoadPGM());
		int dx, dy;
		switch (p[2].p_sel.sel_cur) {
			case 0: dx=0; dy=1; break;
			case 1: dx=1; dy=1; break;
			case 2: dx=1; dy=0; break;
			case 3: dx=-1; dy=1; break;
		}
		dx*=p[3].p_int.val; dy*=p[3].p_int.val; cout<<" decalage en ligne et colonne = "<<dx<<", "<<dy<<"\n";
    imadata<float> imares=ima.ima_carac_cooccurence(p[1].p_int.val,dx,dy,p[4].p_sel.sel_cur);
		imares.statbasic(1);
    imares.imaunsignedchar(1).sauve_ImaPGM(p[5].p_fname.fname);
    return true;
  }
  static bool op_cooc_all_gris(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<BYTE> ima((imadata<BYTE>)ifich.LoadPGM());
		int W=p[1].p_int.val, dx, dy, i,j,k;
		char nomfich0[20], nomfich[30];
		const int ndir=4, nech=3;
		imadata<float> imares[ndir*nech];
		for (i=0; i<ndir; i++) {
			switch (i) {
				case 0: dx=0; dy=1; strcpy(nomfich0,"ima_cooc_h"); break;
				case 1: dx=1; dy=1; strcpy(nomfich0,"ima_cooc_d"); break;
				case 2: dx=1; dy=0; strcpy(nomfich0,"ima_cooc_v"); break;
				case 3: dx=-1; dy=1; strcpy(nomfich0,"ima_cooc_b"); break;
			}
			dx*=W; dy*=W; 
			for (j=0; j<nech; j++) {cout<<" decalage en ligne et colonne = "<<dx<<", "<<dy<<"\n";
				k=i*nech+j;
		    imares[k]=ima.ima_carac_cooccurence(W,dx,dy,p[2].p_sel.sel_cur);
				imares[k].statbasic(1);
				strcpy(nomfich,nomfich0); 
				switch (j) {
					case 0: strcat(nomfich,"_W1.pgm"); break;
					case 1: strcat(nomfich,"_W2.pgm"); break;
					case 2: strcat(nomfich,"_W4.pgm"); break;
				}
				imares[k].imaunsignedchar().sauve_ImaPGM(nomfich);
				dx/=2; dy/=2; 
			}
		}
		imadata<float> imares2;
		for (j=0; j<nech; j++) {cout<<" ratio h/v \n";
			imares2=imares[j]/imares[2*nech+j];
			imares2.statbasic(1);
			strcpy(nomfich,"ima_cooc_ratiohv"); 
			switch (j) {
				case 0: strcat(nomfich,"_W1.pgm"); break;
				case 1: strcat(nomfich,"_W2.pgm"); break;
				case 2: strcat(nomfich,"_W4.pgm"); break;
			}
			imares2.imaunsignedchar(1).sauve_ImaPGM(nomfich);
		}
		imares2.mise_a_zero();
		for (i=0; i<imares2.nlig(); i++)
			for (j=0; j<imares2.ncol(); j++) {
				float xmax=imares[0](i,j);
				for (k=1; k<ndir*nech; k++)	if (imares[k](i,j)>xmax) {xmax=imares[k](i,j); imares2(i,j)=(float)k;}
			}
		cout<<" stat imares2 :\n"; imares2.statbasic(1); //imares2=imares2*20; imares2.statbasic(1);
		(imares2*20).imaunsignedchar().sauve_ImaPGM(p[3].p_fname.fname);
    return true;
  }
  static bool op_coocurence(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<BYTE> ima((imadata<BYTE>)ifich.LoadPGM());
		int W=p[1].p_int.val, dx, dy, i,j,k;
		const int ndir=4, nech=4;
		imadata<float> imares[ndir*nech];
		for (i=0; i<ndir; i++) {
			switch (i) {
				case 0: dx=0; dy=1; break;
				case 1: dx=1; dy=1; break;
				case 2: dx=1; dy=0; break;
				case 3: dx=-1; dy=1; break;
			}
			dx*=W; dy*=W; 
			for (j=0; j<nech; j++) {cout<<" decalage en ligne et colonne = "<<dx<<", "<<dy<<"\n";
				k=j*ndir+i;
		    imares[k]=ima.cooccurence(dx,dy);// imares[k].statbasic(1);
				dx/=2; dy/=2; 
			}
		}
		int nblig2=imares[0].nlig()*nech, nbcol2=imares[0].ncol()*ndir,ii,jj,i0,j0;
		imadata<float> imares2(nblig2,nbcol2);
		for (i=0; i<ndir; i++)
			for (j=0; j<nech; j++) {
				k=j*ndir+i; i0=j*nblig2/nech; j0=i*nbcol2/ndir;
				for (ii=0; ii<nblig2/nech; ii++)
					for (jj=0; jj<nbcol2/ndir; jj++) imares2(ii+i0,jj+j0)=imares[k](ii,jj);
			}
		imares2.statbasic(1);
		imares2.imaunsignedchar(1).sauve_ImaPGM(p[2].p_fname.fname);
    return true;
  }
  static bool op_cooc_binaire(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imabin ima((imadata<BYTE>)ifich.LoadPGM(),1);
    imadata<float> imares=ima.coocurrence(p[1].p_int.val,p[2].p_sel.sel_cur,p[3].p_int.val,p[4].p_sel.sel_cur);
		imares.statbasic(1);
    imares.imaunsignedchar(1).sauve_ImaPGM(p[5].p_fname.fname);
    return true;
  }
	static bool op_transfo_geom(int nbp, Param_Type* p)
	{
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon((imadata<float>)ifich.LoadPGM());
		float dx=p[1].p_float.val,dy=p[2].p_float.val,phi=p[3].p_float.val,hfact=p[4].p_float.val;
		bool resize=1;
		unsigned char fill=0;
		imadata<float> imares=imadon.projette_ima((double)dx,(double)dy,(double)hfact,(double)phi,resize,fill);
    imares.imaunsignedchar(1).sauve_ImaPGM(p[5].p_fname.fname);
		return true;
	}
	static bool op_sous_ima(int nbp, Param_Type* p)
	{
    fichimage_entree ifich(p[0].p_fname.fname);
    imadata<float> imadon((imadata<float>)ifich.LoadPGM());
		imadon.sousIma(p[1].p_int.val+p[3].p_int.val,p[2].p_int.val+p[4].p_int.val,p[3].p_int.val,p[4].p_int.val,p[5].p_int.val,p[6].p_int.val);
    imadon.imaunsignedchar(0).sauve_ImaPGM(p[7].p_fname.fname);
		return true;
	}
	static bool op_rectif_corrphase(int nbp, Param_Type* p)
	{
    fichimage_entree ific1(p[0].p_fname.fname), ific2(p[1].p_fname.fname);
    imadata<float> imadon1((imadata<float>)ific1.LoadPGM()), imadon2((imadata<float>)ific2.LoadPGM());
		int dx,dy; 
		float corrmax,corr,var1,var2;
		imadata<float> log_corrphase_TF=imadon1.correlation_phase (imadon2,dx,dy,corrmax,1,0);
		int nblig=imadon1.nlig(), nbcol=imadon1.ncol(), nblig2, nbcol2, i,j,ii,jj, i0=dx-1,i2=dx+1,j0=dy-1,j2=dy+1, n;
		imadata<float> imadon2_to_1;
		corrmax=0.;
		bool resize=1; 
		unsigned char fill=0;
		for (ii=i0; ii<=i2; ii++) 
			for (jj=j0; jj<=j2; jj++) {
				imadon2_to_1=imadon2.projette_ima((double)(-ii),(double)(-jj),1.,0.,resize,fill);
				corr=var1=var2=0.; n=0;
				nblig2=imadon2_to_1.nlig(); nbcol2=imadon2_to_1.ncol();
				for (i=0; i<mini(nblig,nblig2); i++)
					for (j=0; j<mini(nbcol,nbcol2); j++) 
						if (imadon2_to_1(i,j)!=0 && imadon1(i,j)!=0) {corr+=imadon2_to_1(i,j)*imadon1(i,j); var1+=imadon1(i,j)*imadon1(i,j); var2+=imadon2_to_1(i,j)*imadon2_to_1(i,j); n++;}
				if (n>0) {corr/=n; var1/=n; var2/=n;}
				if (corr>corrmax) {corrmax=corr; dx=ii; dy=jj;}
				cout<<" Translation testee entre les 2 images : ("<<ii<<","<<jj<<") => correlation = "<<corr/pow(var1*var2,0.5f)<<"\n";
			}
		cout<<"\n => Translation trouvee entre les 2 images : ("<<dx<<","<<dy<<")\n";
		imadon2_to_1=imadon2.projette_ima((double)(-dx),(double)(-dy),1.,0.,resize,fill);
    imadon2_to_1.imaunsignedchar(1).sauve_ImaPGM(p[2].p_fname.fname);
		return true;
	}
	static bool op_rectif_FourrierMellin(int nbp, Param_Type* p)
	{
    fichimage_entree ific1(p[0].p_fname.fname), ific2(p[1].p_fname.fname);
    imadata<float> imadon1((imadata<float>)ific1.LoadPGM()), imadon2((imadata<float>)ific2.LoadPGM());
		float dx,dy,rapport,angle;
		bool resize=1;
		unsigned char fill=2;
		imadata<float> log_corrphase_TF, imadon2_to_1(imadon2), imadon2_to_1b(imadon2);
		const float eps=(float)(1.e-3);
		const int n_iter=4;
		int k;
		for (k=0; k<n_iter; k++) {
			if (k>0) cout<<"\n\n On affine recalage en recalculant ("<<k+1<<"eme iteration) transformee :\n";
			log_corrphase_TF=imadon1.transformee_FourierMellin (imadon2_to_1,dx,dy,rapport,angle);
			if ((fabs(rapport-1)>=eps || fabs(angle)>=eps) && k<n_iter-1) {
				imadon2_to_1=imadon2_to_1.projette_ima(0.,0.,(double)rapport,(double)angle,resize,fill);
				imadon2_to_1b=imadon2_to_1b.projette_ima(0.,0.,(double)rapport,(double)angle,resize,0);
			} else {
				imadon2_to_1=imadon2_to_1.projette_ima((double)dx,(double)dy,(double)rapport,(double)angle,resize,fill);
				imadon2_to_1b=imadon2_to_1b.projette_ima((double)dx,(double)dy,(double)rapport,(double)angle,resize,0);
				if (fabs(rapport-1)<eps && fabs(angle)<eps) k=n_iter;
			}
			imadon2_to_1b.imaunsignedchar(1).sauve_ImaPGM(p[2].p_fname.fname);
			switch (k) {
				case 0: imadon2_to_1.imaunsignedchar(1).sauve_ImaPGM("imadon2_to_1_iter1.pgm"); imadon2_to_1b.imaunsignedchar(1).sauve_ImaPGM("imadon2_to_1b_iter1.pgm"); break;
				case 1: imadon2_to_1.imaunsignedchar(1).sauve_ImaPGM("imadon2_to_1_iter2.pgm"); imadon2_to_1b.imaunsignedchar(1).sauve_ImaPGM("imadon2_to_1b_iter2.pgm"); break;
				case 2: imadon2_to_1.imaunsignedchar(1).sauve_ImaPGM("imadon2_to_1_iter3.pgm"); imadon2_to_1b.imaunsignedchar(1).sauve_ImaPGM("imadon2_to_1b_iter3.pgm"); break;
				case 3: imadon2_to_1.imaunsignedchar(1).sauve_ImaPGM("imadon2_to_1_iter4.pgm"); imadon2_to_1b.imaunsignedchar(1).sauve_ImaPGM("imadon2_to_1b_iter4.pgm"); break;
				case 4: imadon2_to_1.imaunsignedchar(1).sauve_ImaPGM("imadon2_to_1_iter5.pgm"); imadon2_to_1b.imaunsignedchar(1).sauve_ImaPGM("imadon2_to_1b_iter5.pgm"); break;
				default: imadon2_to_1.imaunsignedchar(1).sauve_ImaPGM("imadon2_to_1_iterN.pgm"); imadon2_to_1b.imaunsignedchar(1).sauve_ImaPGM("imadon2_to_1b_iterN.pgm"); 
			}
		}
		imadon2_to_1b.sousIma(imadon1.nlig(),imadon1.ncol(),0,0,1,1);
		imadon2_to_1b.imaunsignedchar(1).sauve_ImaPGM(p[2].p_fname.fname);
		return true;
	}
	static bool op_rectif_carre_imacol(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imacouleur<BYTE> imacol((imadata<BYTE>)ifich.LoadPGM()), imacol_bis(imacol);
		imadata<float> imV=imacol.select_composante(1); imV.imaunsignedchar().sauve_ImaPGM("imaV.pgm");
		imacol_bis.RVB2IST(); imadata<float> imI=imacol_bis.select_composante(0); imI.imaunsignedchar().sauve_ImaPGM("imaI.pgm");
		int i,j, nbli=imI.nlig(), nbco=imI.ncol();
		imabin imab(nbli,nbco); imab.mise_a_zero();
		for (i=0; i<nbli; i++)
			for (j=0; j<nbco; j++)
				if (imV(i,j)>127 && (float)imV(i,j)/(float)imI(i,j)>1.5) imab(i,j)=1;
		imab.imaunsignedchar(1).sauve_ImaPGM("imabin.pgm");
		eltstruct S3(3,3);
		const int npoints=4, d=3;
		int **coord=new int*[npoints], k, nbli2=240, nbco2=240;
		matrice2D<double> X(npoints,d), Y(npoints,2);
		for (k=0; k<npoints; k++) coord[k]=new int[4];
		coord[0][0]=0; coord[0][1]=0;  
		coord[1][0]=nbco2; coord[1][1]=0;  
		coord[2][0]=nbco2; coord[2][1]=nbli2; 
		coord[3][0]=0; coord[3][1]=nbli2; 
		bool itrouve=0;
		for (i=0; i<nbli; i++)
			for (j=0; j<=i; j++)
				if (!imab(i-j,j)) {coord[0][2]=j; coord[0][3]=i-j; i=nbli; j=nbco; itrouve=1;}
		if (itrouve) cout<<" coin UL : "<<coord[0][3]<<" "<<coord[0][2]<<"\n";
		else cout<<" coin UL non trouve !!!!!\n";
		itrouve=0;
		for (i=0; i<nbli; i++)
			for (j=nbco-1; j>=maxi(nbco-1-i,0); j--) 
				if (!imab(i+j-nbco+1,j)) {coord[1][2]=j; coord[1][3]=i+j-nbco+1; i=nbli; j=-1; itrouve=1;}
		if (itrouve) cout<<" coin UR : "<<coord[1][3]<<" "<<coord[1][2]<<"\n";
		else cout<<" coin UR non trouve !!!!!\n";
		itrouve=0;
		for (i=nbli-1; i>=0; i--)
			for (j=nbco-1; j>=nbco-1-(nbli-1-i); j--)//{cout<<i+j-nbco+1<<" "<<j<<"\n";
				if (!imab(i-(j-nbco+1),j)) {coord[2][2]=j; coord[2][3]=i-(j-nbco+1); i=-1; j=-1; itrouve=1;}
			//}
		if (itrouve) cout<<" coin LR : "<<coord[2][3]<<" "<<coord[2][2]<<"\n";
		else cout<<" coin LR non trouve !!!!!\n";
		itrouve=0;
		for (i=nbli-1; i>=0; i--)
			for (j=0; j<=nbli-1-i; j++)
				if (!imab(i+j,j)) {coord[3][2]=j; coord[3][3]=i+j; i=-1; j=nbco; itrouve=1;}
		if (itrouve) cout<<" coin LL : "<<coord[3][3]<<" "<<coord[3][2]<<"\n";
		else cout<<" coin LL non trouve !!!!!\n";
		for (k=0; k<npoints; k++) {
			X(k,0)=1; X(k,1)=coord[k][0]; X(k,2)=coord[k][1];
			if (d>5) {X(k,3)=pow(X(k,1),2); X(k,4)=pow(X(k,2),2); X(k,5)=X(k,1)*X(k,2);}
			Y(k,0)=coord[k][2]; Y(k,1)=coord[k][3];
		}
		for (k=0; k<npoints; k++) if (coord[k]!=NULL) delete[] coord[k];
		if (coord!=NULL) delete[] coord;
		matrice2D<double> Xt=X.transpo();
		matrice2D<double> XXt(Xt,X,d,npoints,d); 
		bool invOK;
		matrice2D<double> XXt_1=XXt.inverse(invOK);
		cout<<" verification de l'inversion de X.Xt : inversion ";
		matrice2D<double> Id(d,d); for (i=0; i<d; i++) Id(i,i)=1; 
		(matrice2D<double>(XXt,XXt_1,d,d,d)==Id && matrice2D<double>(XXt_1,XXt,d,d,d)==Id)?cout<<"OK":cout<<"Not OK";cout<<"\n";
		matrice2D<double> A(d,2);
		for (k=0; k<d; k++) A(k,0)=A(k,1)=0; A(1,0)=A(2,1)=1;
		A=matrice2D<double>(matrice2D<double>(XXt_1,Xt,d,d,npoints),Y,d,npoints,2);
		matrice2D<double> D=Y-matrice2D<double>(X,A,npoints,d,2); 
		cout<<" erreur d'estimation par point d'amer\n"; D.affiche();
		cout<<" coefficients de la transformation\n"; A.affiche();
// Application de la transformation et projection de l'image
		imadata<float> imacol_proj(imacol);
		imacol_proj=imacol.projette_ima(A,nbli2,nbco2,imab,1);
		imacol_proj.imaunsignedchar().sauve_ImaPGM(p[1].p_fname.fname);
    return true;
  }
  static bool op_flot_optique(int nbp, Param_Type* p)
  {
    fichimage_entree ific1(p[0].p_fname.fname), ific2(p[1].p_fname.fname);
    imadata<float> ima1((imadata<float>)ific1.LoadPGM()), ima2((imadata<float>)ific2.LoadPGM());
    imaflot<float> imares;
	/* constructeur Horn and Schunck */
    switch (p[3].p_sel.sel_cur) {
      case 0: imares=imaflot<float>(ima1,ima2,p[2].p_float.val); break;
    }
    imares.imaunsignedchar(1).sauve_ImaPGM(p[4].p_fname.fname);
    return true;
  }
  static bool op_pb1(int nbp, Param_Type* p)
  { 
		float haut=2.75f, K1=248.f, K2=19.1f, focale=K1/haut, nblig_ima=39.f, alpha=(nblig_ima/2.f-K2)/focale, nbcol_ima=52.f, colcent_ima=(float)nbcol_ima/2.f, gamma=2*atan(colcent_ima/focale); // parametres de la caméra
		cout<<" par defaut : focale = "<<focale<<", ouverture = "<<gamma*180/PI<<"\n";
		float entreAxe_Roues=5.2f, rayon_Roue=4.0f/2, V_rot=40.f, V_transl=200.f, pas_Z=10.f, Dmin=K1/(nblig_ima-K2), prec=0.1f, Rplot=5.f; // parametres du déplacement du robot
		ifstream fcontexte; fcontexte.open(p[1].p_fname.fname,ios::in);
    if (!fcontexte) {
      cout<<" ouverture de du fichier de contexte "<<fcontexte<<" impossible -> on prend les parametres par defaut.\n";
    }
    else {
      fcontexte>>K1>>K2>>V_rot>>V_transl>>pas_Z>>Dmin;
			focale=K1/haut; alpha=(nblig_ima/2.f-K2)/focale; gamma=2*atan(colcent_ima/focale); // parametres de la caméra
			cout<<" Pour les parametres specifies : focale = "<<focale<<", ouverture = "<<gamma*180/PI<<", angle d'inclinaison = "<<alpha*180/PI<<"\n";
			cout<<" Vitesse dans les tournants sur place = "<<V_rot<<", vitesse en translation = "<<V_transl<<", pas de deplacement en translation = "<<pas_Z<<", distance minimale au plot = "<<Dmin<<"\n";
      fcontexte.close();
		}
		imacouleur<BYTE> imacol;
		imadata<float> imaR, imaG, imaB, imaRG, imaI;
		imabin imabR, imabG, imabB, imab;
    eltstruct S3(3,3); 
		int ncc, i, n, lmin, cmin, cmax;
		float lbary, cbary, surf, surf_max, compac, compac_max, theta, VG, VD, delta_t, dist;
		imadata<int> imaCC;
		bool continuer, ecrire;
		unsigned char nbcouleurs=3, couleur=0, etat=0;  // etats : 0=recherche_plot, 1=avance_vers_plot, 2&3=contourne_plot; couleurs : 0=Rouge, 1=Jaune, 2=Vert;
		cout<<" debut de la boucle infinie (taper 0 pour en sortir)\n";
		do { 
			VG=VD=delta_t=0.f; ecrire=0;
			if (etat==0 || etat==1) {
		    fichimage_entree ifich(p[0].p_fname.fname);
				imacol=(imadata<BYTE>)ifich.LoadPGM();
				imaR=imacol.select_composante(0); imaG=imacol.select_composante(1); imaB=imacol.select_composante(2); imaI=(imaR+imaG+imaB); 
				switch (couleur) {
					case 0: imaRG=(imaR/(imaG+1.f)); imaRG.seuil_hysteresis(2.f,1.5f,imabR); imabR.imaunsignedchar().sauve_ImaPGM("./imabRsurG.pgm");
									imabG=imabin(imaG/(imaI+1.f),0.30f+0.000001f).negatif(); imabG.imaunsignedchar().sauve_ImaPGM("./imabG.pgm");
									imab=imabR&&imabG; imab.imaunsignedchar().sauve_ImaPGM("./imabRetG.pgm");
									imab=imab.ouverture(S3,+1).fermeture(S3,+1); break;
					case 1: (imaR/(imaB+1.f)).seuil_hysteresis(1.5f,1.f,imabR); imabR.imaunsignedchar().sauve_ImaPGM("./imabR.pgm");
									(imaG/(imaB+1.f)).seuil_hysteresis(1.5f,1.f,imabG); imabG.imaunsignedchar().sauve_ImaPGM("./imabG.pgm");
									imabB=imabin(imaB/(imaI+1.f),0.30f+0.000001f).negatif(); imabB.imaunsignedchar().sauve_ImaPGM("./imabB.pgm");
									imab=imabR&&imabG&&imabB; imab.imaunsignedchar().sauve_ImaPGM("./imabRG.pgm");
									imab=imab.ouverture(S3,+1).fermeture(S3,+1); break;
					case 2: imaRG=(imaG/(imaR+1.f)); imaRG.seuil_hysteresis(2.f,1.5f,imabG); imabG.imaunsignedchar().sauve_ImaPGM("./imabGsurR.pgm");
									imabR=imabin(imaR/(imaI+1.f),0.30f+0.000001f).negatif(); imabR.imaunsignedchar().sauve_ImaPGM("./imabR.pgm");
									imab=imabR&&imabG; imab.imaunsignedchar().sauve_ImaPGM("./imabRetG.pgm");
									imab=imab.ouverture(S3,+1).fermeture(S3,+1); break;
					default: break;
				}
				imab.imaunsignedchar().sauve_ImaPGM(p[2].p_fname.fname);
				ncc=0; imaCC=imab.composantes_connexes(ncc,8); cout<<" nombre de composantes connexes = "<<ncc<<"\n";
				if (ncc>1) {
					imabin imamarq(imab.nlig(),imab.ncol());
					for (i=0; i<imamarq.nlig(); i++) {
						n=imamarq.ncol()-1;
						imamarq(i,0)=imamarq(i,n)=1;
					}
					for (i=0; i<imamarq.ncol(); i++) {
						n=imamarq.nlig()-1;
						imamarq(0,i)=imamarq(n,i)=1;
					}
					imab=imab-imab.reconstruction_geodesique(imamarq,8);
					imaCC=imab.composantes_connexes(ncc,8); cout<<" apres elimination des objets touchant le bord : nombre de composantes connexes = "<<ncc<<"\n";
				}
				if (ncc>1) {
					imadata<float> imaBsurI=imaB/(imaI+1.f);
					n=imaBsurI.nlig()-1;
					for (i=0; i<imaBsurI.ncol(); i++) imaBsurI(n,i)=1.f;
					imabB=imabin(imaB,50.f);
					imaBsurI.seuil_hysteresis(0.9f,0.33f,imabB); 
					imabB=imabB&&imabin(imaB,50.f); imabB.imaunsignedchar().sauve_ImaPGM("./imabB.pgm");
					imabB=imabB.dilate(S3,+2);
					imab=imab.reconstruction_geodesique(imabB,8);
					imaCC=imab.composantes_connexes(ncc,8); cout<<" apres elimination des objets ne touchant pas le tapis : nombre de composantes connexes = "<<ncc<<"\n";
				}
				imaobjets ima_obj(imab.imaunsignedchar()); ima_obj.basic_param(); ima_obj.affiche();
				objet obj=ima_obj.objet_i(1);
				if (ncc>1) {
					surf_max=obj.surf_obj(); compac_max=surf_max;
					if (obj.l_max_obj()-obj.l_min_obj()>0) compac_max/=(obj.l_max_obj()-obj.l_min_obj()); if (obj.c_max_obj()-obj.c_min_obj()>0) compac_max/=(obj.c_max_obj()-obj.c_min_obj());
					for (i=2; i<=ncc; i++) {
						surf=ima_obj.objet_i(i).surf_obj(); compac=surf;
						if (ima_obj.objet_i(i).l_max_obj()-ima_obj.objet_i(i).l_min_obj()>0) compac/=(ima_obj.objet_i(i).l_max_obj()-ima_obj.objet_i(i).l_min_obj()); 
						if (ima_obj.objet_i(i).c_max_obj()-ima_obj.objet_i(i).l_min_obj()>0) compac/=(ima_obj.objet_i(i).c_max_obj()-ima_obj.objet_i(i).c_min_obj());
						if (compac>compac_max && (surf>surf_max || surf>surf_max*0.8)) {compac_max=compac; surf_max=maxi(surf,surf_max); obj=ima_obj.objet_i(i);}
					}
				}
				lmin=obj.l_max_obj(); lbary=obj.l_baryc_obj(); cbary=obj.c_baryc_obj(); cmin=obj.c_min_obj(); cmax=obj.c_max_obj(); 
				if (ncc<1) etat=0;
				else etat=1;
			}
			cout<<" etat "<<(int)etat<<" cible couleur "<<(int)couleur<<"\n";
			switch (etat) {
				case 0 : // recherche de la cible
					theta=gamma*(1.f-prec);
					cout<<"\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
					cout<<" etat "<<(int)etat<<" cible couleur "<<(int)couleur<<" non visible\n -> tourner sur place de l'angle "<<theta*180/PI<<" grace a la commande \n";
					VG=-V_rot; VD=V_rot; delta_t=1000*(float)(theta*1000/PI/4.*entreAxe_Roues/rayon_Roue/V_rot); ecrire=1;
					cout<<" VG="<<VG<<", VD="<<VD<<" pendant "<<delta_t<<" milliemes de seconde";
					cout<<"\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
					break;
				case 1 : // approche de la cible
					if (cmax<colcent_ima || colcent_ima<cmin) { // on recentre
//					if (cbary<colcent_ima*(1-2*prec)||cbary>colcent_ima*(1+2*prec)) { // on recentre
						theta=atan((colcent_ima-cbary)/focale);
						cout<<"\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
						cout<<" etat "<<(int)etat<<" pour avoir la cible couleur "<<(int)couleur<<" dans l'axe,\n tourner sur place de l'angle "<<theta*180/PI<<" grace a la commande \n";
						VG=-V_rot; VD=+V_rot; delta_t=1000*(float)(theta*1000/PI/4.*entreAxe_Roues/rayon_Roue/V_rot); if (delta_t<0) {VG*=-1; VD*=-1; delta_t*=-1;} ecrire=1;
						cout<<" VG="<<VG<<", VD="<<VD<<" pendant "<<delta_t<<" milliemes de seconde";
						cout<<"\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
					} else {
						if (lmin<(nblig_ima-1)) { // on avance
							dist=K1/(lmin-K2);
							cout<<"\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
							cout<<" etat "<<(int)etat<<" obstacle couleur "<<(int)couleur<<" a la distance "<<dist<<"\n -> avancer d'une distance "<<min(pas_Z,dist-Dmin*(1+prec))<<" soit la commande \n";
							VG=V_transl; VD=V_transl; delta_t=1000*(float)(min(pas_Z,dist-Dmin*(1+prec))*1000./2./PI/rayon_Roue/V_transl); if (delta_t<0) {VG*=-1; VD*=-1; delta_t*=-1;} ecrire=1;
							cout<<" VG="<<VG<<", VD="<<VD<<" pendant "<<delta_t<<" milliemes de seconde";
							cout<<"\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
						} else etat=2;
					} break;
				case 2 : // contournement de la cible
					theta=(float)PI/2.f;
					cout<<"\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
					cout<<" etat "<<(int)etat<<" pour contourner la cible couleur "<<(int)couleur<<",\n d'abord tourner sur place de l'angle "<<theta*180/PI<<" grace a la commande \n";
					VG=-V_rot; VD=V_rot; delta_t=1000*(float)(theta*1000/PI/4.*entreAxe_Roues/rayon_Roue/V_rot); ecrire=1;
					cout<<" VG="<<VG<<", VD="<<VD<<" pendant "<<delta_t<<" milliemes de seconde";
					cout<<"\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
					etat=3; break;
				case 3 : 
					theta=3.f*(float)PI/2.f*(1.f+prec);
					cout<<"\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
					cout<<" puis il faut contourner le plot avec l'angle "<<theta*180/PI<<" grace a la commande \n";
					VG=V_transl; VD=V_transl*((Dmin+Rplot)-entreAxe_Roues/2.f)/((Dmin+Rplot)+entreAxe_Roues/2.f); delta_t=1000*(float)(theta*1000/(float)PI/2.f*((Dmin+Rplot)+entreAxe_Roues/2.f)/rayon_Roue/V_transl); ecrire=1;
//					VG=V_transl; VD=V_transl*(Dmin-entreAxe_Roues/2.)/(Dmin+entreAxe_Roues/2.); delta_t=1000*(float)(theta*1000/PI/2.*(Dmin+entreAxe_Roues/2.)/rayon_Roue/V_transl); ecrire=1;
					cout<<" VG="<<VG<<", VD="<<VD<<" pendant "<<delta_t<<" milliemes de seconde";
					cout<<"\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
					etat=4; break;
				case 4 : 
					couleur++; couleur=couleur%nbcouleurs;
					etat=0; break;
			}
			if (ecrire) {
				ofstream fcommande; fcommande.open(p[3].p_fname.fname,ios::out);
				if (!fcommande) {cout<<" ouverture de du fichier de commande "<<fcommande<<" impossible !!!\n";}
				else {
					fcommande<<(int)VG<<" "<<(int)VD<<" "<<(int)delta_t<<"\n";
					fcommande.close();
				}
				continuer=(getchar()!=0);
			}
		} while (continuer);

    return true;
  }
  static bool op_pb1bis(int nbp, Param_Type* p)
  { 
		float haut=2.75f, K1=248.f, K2=19.1f, focale=K1/haut, nblig_ima=39.f, alpha=(nblig_ima/2.f-K2)/focale, nbcol_ima=52.f, colcent_ima=(float)nbcol_ima/2.f, gamma=2*atan(colcent_ima/focale); // parametres de la caméra
		cout<<" par defaut : focale = "<<focale<<", ouverture = "<<gamma*180/PI<<"\n";
		float entreAxe_Roues=5.2f, rayon_Roue=4.0f/2, V_rot=40.f, V_transl=200.f, pas_Z=10.f, Dmin=K1/(nblig_ima-K2), prec=0.1f, Rplot=5.f; // parametres du déplacement du robot
		ifstream fcontexte; fcontexte.open(p[0].p_fname.fname,ios::in);
    if (!fcontexte) {
      cout<<" ouverture de du fichier de contexte "<<fcontexte<<" impossible -> on prend les parametres par defaut.\n";
    }
    else {
      fcontexte>>K1>>K2>>V_rot>>V_transl>>pas_Z>>Dmin;
			focale=K1/haut; alpha=(nblig_ima/2.f-K2)/focale; gamma=2*atan(colcent_ima/focale); // parametres de la caméra
			cout<<" Pour les parametres specifies : focale = "<<focale<<", ouverture = "<<gamma*180/PI<<", angle d'inclinaison = "<<alpha*180/PI<<"\n";
			cout<<" Vitesse dans les tournants sur place = "<<V_rot<<", vitesse en translation = "<<V_transl<<", pas de deplacement en translation = "<<pas_Z<<", distance minimale au plot = "<<Dmin<<"\n";
      fcontexte.close();
		}
		imacouleur<BYTE> imacol((int)nblig_ima,(int)nbcol_ima);
		imadata<float> imaR, imaG, imaB, imaRG, imaI;
		imabin imabR, imabG, imabB, imab;
    eltstruct S3(3,3); 
		int size_ligima=(int)nbcol_ima*3, ncc, i, n, lmin, cmin, cmax;
		float lbary, cbary, surf, surf_max, compac, compac_max, theta, VG, VD, delta_t, dist;
		imadata<int> imaCC;
		bool continuer=1, ecrire;
		unsigned char nbcouleurs=3, couleur=0, etat=0;  // etats : 0=recherche_plot, 1=avance_vers_plot, 2&3=contourne_plot; couleurs : 0=Rouge, 1=Jaune, 2=Vert;
		cout<<" debut de la boucle infinie (taper 0 pour en sortir)\n";
		CSharedData<RobotState> Data(SH_NAME,false);
		do { 
			VG=VD=delta_t=0.f; ecrire=0;
			if (etat==0 || etat==1) {
				Data->read_cam=1;
				for (i=0;i<Data->IMG_SZ; i++) imacol(i/size_ligima,(i%size_ligima)/3,i%3)=Data->img[i];
				Data->read_cam=0;
				imaR=imacol.select_composante(0); imaG=imacol.select_composante(1); imaB=imacol.select_composante(2); imaI=(imaR+imaG+imaB); 
				switch (couleur) {
					case 0: imaRG=(imaR/(imaG+1.f)); imaRG.seuil_hysteresis(2.f,1.5f,imabR); imabR.imaunsignedchar().sauve_ImaPGM("./imabRsurG.pgm");
									imabG=imabin(imaG/(imaI+1.f),0.30f+0.000001f).negatif(); imabG.imaunsignedchar().sauve_ImaPGM("./imabG.pgm");
									imab=imabR&&imabG; imab.imaunsignedchar().sauve_ImaPGM("./imabRetG.pgm");
									imab=imab.ouverture(S3,+1).fermeture(S3,+1); break;
					case 1: (imaR/(imaB+1.f)).seuil_hysteresis(1.5f,1.f,imabR); imabR.imaunsignedchar().sauve_ImaPGM("./imabR.pgm");
									(imaG/(imaB+1.f)).seuil_hysteresis(1.5f,1.f,imabG); imabG.imaunsignedchar().sauve_ImaPGM("./imabG.pgm");
									imabB=imabin(imaB/(imaI+1.f),0.30f+0.000001f).negatif(); imabB.imaunsignedchar().sauve_ImaPGM("./imabB.pgm");
									imab=imabR&&imabG&&imabB; imab.imaunsignedchar().sauve_ImaPGM("./imabRG.pgm");
									imab=imab.ouverture(S3,+1).fermeture(S3,+1); break;
					case 2: imaRG=(imaG/(imaR+1.f)); imaRG.seuil_hysteresis(2.f,1.5f,imabG); imabG.imaunsignedchar().sauve_ImaPGM("./imabGsurR.pgm");
									imabR=imabin(imaR/(imaI+1.f),0.30f+0.000001f).negatif(); imabR.imaunsignedchar().sauve_ImaPGM("./imabR.pgm");
									imab=imabR&&imabG; imab.imaunsignedchar().sauve_ImaPGM("./imabRetG.pgm");
									imab=imab.ouverture(S3,+1).fermeture(S3,+1); break;
					default: break;
				}
				imab.imaunsignedchar().sauve_ImaPGM(p[1].p_fname.fname);
				ncc=0; imaCC=imab.composantes_connexes(ncc,8); cout<<" nombre de composantes connexes = "<<ncc<<"\n";
				if (ncc>1) {
					imabin imamarq(imab.nlig(),imab.ncol());
					for (i=0; i<imamarq.nlig(); i++) {
						n=imamarq.ncol()-1;
						imamarq(i,0)=imamarq(i,n)=1;
					}
					for (i=0; i<imamarq.ncol(); i++) {
						n=imamarq.nlig()-1;
						imamarq(0,i)=imamarq(n,i)=1;
					}
					imab=imab-imab.reconstruction_geodesique(imamarq,8);
					imaCC=imab.composantes_connexes(ncc,8); cout<<" apres elimination des objets touchant le bord : nombre de composantes connexes = "<<ncc<<"\n";
				}
				if (ncc>1) {
					imadata<float> imaBsurI=imaB/(imaI+1.f);
					n=imaBsurI.nlig()-1;
					for (i=0; i<imaBsurI.ncol(); i++) imaBsurI(n,i)=1.f;
					imabB=imabin(imaB,50.f);
					imaBsurI.seuil_hysteresis(0.9f,0.33f,imabB); 
					imabB=imabB&&imabin(imaB,50.f); imabB.imaunsignedchar().sauve_ImaPGM("./imabB.pgm");
					imabB=imabB.dilate(S3,+2);
					imab=imab.reconstruction_geodesique(imabB,8);
					imaCC=imab.composantes_connexes(ncc,8); cout<<" apres elimination des objets ne touchant pas le tapis : nombre de composantes connexes = "<<ncc<<"\n";
				}
				imaobjets ima_obj(imab.imaunsignedchar()); ima_obj.basic_param(); ima_obj.affiche();
				objet obj=ima_obj.objet_i(1);
				if (ncc>1) {
					surf_max=obj.surf_obj(); compac_max=surf_max;
					if (obj.l_max_obj()-obj.l_min_obj()>0) compac_max/=(obj.l_max_obj()-obj.l_min_obj()); if (obj.c_max_obj()-obj.c_min_obj()>0) compac_max/=(obj.c_max_obj()-obj.c_min_obj());
					for (i=2; i<=ncc; i++) {
						surf=ima_obj.objet_i(i).surf_obj(); compac=surf;
						if (ima_obj.objet_i(i).l_max_obj()-ima_obj.objet_i(i).l_min_obj()>0) compac/=(ima_obj.objet_i(i).l_max_obj()-ima_obj.objet_i(i).l_min_obj()); 
						if (ima_obj.objet_i(i).c_max_obj()-ima_obj.objet_i(i).l_min_obj()>0) compac/=(ima_obj.objet_i(i).c_max_obj()-ima_obj.objet_i(i).c_min_obj());
						if (compac>compac_max && (surf>surf_max || surf>surf_max*0.8)) {compac_max=compac; surf_max=maxi(surf,surf_max); obj=ima_obj.objet_i(i);}
					}
				}
				lmin=obj.l_max_obj(); lbary=obj.l_baryc_obj(); cbary=obj.c_baryc_obj(); cmin=obj.c_min_obj(); cmax=obj.c_max_obj(); 
				if (ncc<1) etat=0;
				else etat=1;
			}
			cout<<" etat "<<(int)etat<<" cible couleur "<<(int)couleur<<"\n";
			switch (etat) {
				case 0 : // recherche de la cible
					theta=gamma*(1.f-prec);
					cout<<"\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
					cout<<" etat "<<(int)etat<<" cible couleur "<<(int)couleur<<" non visible\n -> tourner sur place de l'angle "<<theta*180/PI<<" grace a la commande \n";
					VG=-V_rot; VD=V_rot; delta_t=1000*(float)(theta*1000/PI/4.*entreAxe_Roues/rayon_Roue/V_rot); ecrire=1;
					cout<<" VG="<<VG<<", VD="<<VD<<" pendant "<<delta_t<<" milliemes de seconde";
					cout<<"\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
					break;
				case 1 : // approche de la cible
					if (cmax<colcent_ima || colcent_ima<cmin) { // on recentre
//					if (cbary<colcent_ima*(1-2*prec)||cbary>colcent_ima*(1+2*prec)) { // on recentre
						theta=atan((colcent_ima-cbary)/focale);
						cout<<"\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
						cout<<" etat "<<(int)etat<<" pour avoir la cible couleur "<<(int)couleur<<" dans l'axe,\n tourner sur place de l'angle "<<theta*180/PI<<" grace a la commande \n";
						VG=-V_rot; VD=+V_rot; delta_t=1000*(float)(theta*1000/PI/4.*entreAxe_Roues/rayon_Roue/V_rot); if (delta_t<0) {VG*=-1; VD*=-1; delta_t*=-1;} ecrire=1;
						cout<<" VG="<<VG<<", VD="<<VD<<" pendant "<<delta_t<<" milliemes de seconde";
						cout<<"\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
					} else {
						if (lmin<(nblig_ima-1)) { // on avance
							dist=K1/(lmin-K2);
							cout<<"\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
							cout<<" etat "<<(int)etat<<" obstacle couleur "<<(int)couleur<<" a la distance "<<dist<<"\n -> avancer d'une distance "<<min(pas_Z,dist-Dmin*(1+prec))<<" soit la commande \n";
							VG=V_transl; VD=V_transl; delta_t=1000*(float)(min(pas_Z,dist-Dmin*(1+prec))*1000./2./PI/rayon_Roue/V_transl); if (delta_t<0) {VG*=-1; VD*=-1; delta_t*=-1;} ecrire=1;
							cout<<" VG="<<VG<<", VD="<<VD<<" pendant "<<delta_t<<" milliemes de seconde";
							cout<<"\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
						} else etat=2;
					} break;
				case 2 : // contournement de la cible
					theta=(float)PI/2.f;
					cout<<"\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
					cout<<" etat "<<(int)etat<<" pour contourner la cible couleur "<<(int)couleur<<",\n d'abord tourner sur place de l'angle "<<theta*180/PI<<" grace a la commande \n";
					VG=-V_rot; VD=V_rot; delta_t=1000*(float)(theta*1000/PI/4.*entreAxe_Roues/rayon_Roue/V_rot); ecrire=1;
					cout<<" VG="<<VG<<", VD="<<VD<<" pendant "<<delta_t<<" milliemes de seconde";
					cout<<"\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
					etat=3; break;
				case 3 : 
					theta=3.f*(float)PI/2.f*(1.f+prec);
					cout<<"\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
					cout<<" puis il faut contourner le plot avec l'angle "<<theta*180/PI<<" grace a la commande \n";
					VG=V_transl; VD=V_transl*((Dmin+Rplot)-entreAxe_Roues/2.f)/((Dmin+Rplot)+entreAxe_Roues/2.f); delta_t=1000*(float)(theta*1000/(float)PI/2.f*((Dmin+Rplot)+entreAxe_Roues/2.f)/rayon_Roue/V_transl); ecrire=1;
//					VG=V_transl; VD=V_transl*(Dmin-entreAxe_Roues/2.)/(Dmin+entreAxe_Roues/2.); delta_t=1000*(float)(theta*1000/PI/2.*(Dmin+entreAxe_Roues/2.)/rayon_Roue/V_transl); ecrire=1;
					cout<<" VG="<<VG<<", VD="<<VD<<" pendant "<<delta_t<<" milliemes de seconde";
					cout<<"\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
					etat=4; break;
				case 4 : 
					couleur++; couleur=couleur%nbcouleurs;
					etat=0; break;
			}
			Data->speed_l=(double)VG; Data->speed_r=(double)VD;
			Data->write_speed=mini((int)delta_t,255);
			do {;} while (Data->write_speed>0);
		} while (continuer);

    return true;
  }
  static bool op_pb2(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
    imacouleur<BYTE> imacol((imadata<BYTE>)ifich.LoadPGM());
    imadata<float> imaR=imacol.select_composante(0), imaG=imacol.select_composante(1), imaB=imacol.select_composante(2);
 /************************ debut a completer *************************/
	
 /************************* fin a completer *************************/
    return true;
  }
  static bool op_pb3(int nbp, Param_Type* p)
  {
    char *nomfich=p[0].p_fname.fname, *extention=".pgm", no[4]="000";
    int nbima=maxi(p[1].p_int.val,0);
    fichimage_entree ifich(nomfich); cout<<" ouverture du fichier "<<nomfich<<"\n";
    imadata<float> ima1((imadata<float>)ifich.LoadPGM());
    int nblig=ima1.nlig(),nbcol=ima1.ncol(), ad_debut_no=strlen(nomfich)-strlen(no)-strlen(extention), no0=0;
    int i,j,k,ii,jj,i0,i2,j0,j2,n,sn;
//	strncpy(nomfich+ad_debut_no,no,3);
    strncpy(no,nomfich+ad_debut_no,3); cout<<" numero initial des fichiers "<<no<<"\n";
    i=100; for (j=0; j<3; j++) {no0+=(*(no+j)-0X30)*i; i/=10;} cout<<" no0 = "<<no0<<"\n";
    imadata<float> ima2, imadif0, ima_dir1, ima_dir2, imadif1, imatemp(nblig,nbcol);
    imadata<BYTE> imadec0(nblig,nbcol), imadec1(nblig,nbcol), imaRes(nblig,nbcol);
    const float x_abaisse_ng=10.f;
    const int nl=5, nc=5, nl2=nl/2, nc2=nc/2;
    eltstruct ES(nl,nc);
    double s0_b, s0_h, s1_b, s1_h, xx;
//	ima1.egalisehisto(256,0);
    for (k=no0+1; k<no0+nbima; k++) {
      i=100; for (j=0; j<3; j++) {*(no+j)=0X30+((k/i)%10); i/=10;}
      strncpy(nomfich+ad_debut_no,no,3);
      ifich.changefichier(nomfich);
      cout<<" ouverture du fichier "<<nomfich<<"\n";
      ima2=(imadata<float>)ifich.LoadPGM();
/************************ debut a completer *************************/
//		ima2.egalisehisto(256,0);
      imadif0=(ima2-ima1).absI(); //imadif0.imaunsignedchar().sauve_ImaPGM("imadif00.pgm");
      imadif0=imadif0.reconst_geod(imadif0+(-x_abaisse_ng),ES); imadif0.imaunsignedchar().sauve_ImaPGM("imadif0.pgm");
      imadif1=(ima2.gradient(ima_dir2)-ima1.gradient(ima_dir1)).absI(); imadif1.imaunsignedchar().sauve_ImaPGM("imadif01.pgm");
//		imatemp=(ima_dir2-ima_dir1).abs()*(180.f/PI); imatemp.imaunsignedchar().sauve_ImaPGM("imatemp.pgm");
//		imadif0=imadif0.fermeture(ES).ouverture(ES); imadif0.imaunsignedchar().sauve_ImaPGM("imadif0.pgm");
//		imadif1=imadif1.fermeture(ES).ouverture(ES); imadif1.imaunsignedchar().sauve_ImaPGM("imadif1.pgm");
      imadif0.statbasic(0);
      s0_b=imadif0.moyI()+2*pow(imadif0.varI(),0.5); s0_h=imadif0.moyI()+3*pow(imadif0.varI(),0.5); 
      s1_b=imadif1.moyI()+2*pow(imadif1.varI(),0.5); s1_h=imadif1.moyI()+3*pow(imadif1.varI(),0.5); 
      cout<<" seuils : s0="<<s0_b<<" "<<s0_h<<" s1="<<s1_b<<" "<<s1_h<<"\n";
      for (i=0; i<nblig; i++)
        for (j=0; j<nbcol; j++) {
          if (imadif0(i,j)<s0_b) imadec0(i,j)=0;
          else {
            if (imadif0(i,j)<=s0_h) imadec0(i,j)=1;
            else imadec0(i,j)=2;
          }
          if (imadif1(i,j)<s1_b) imadec1(i,j)=0;
          else {
            if (imadif1(i,j)<=s1_h) imadec1(i,j)=1;
            else imadec1(i,j)=2;
          }
        }
      imadec0.sauve_ImaPGM("imadec0.pgm"); imadec1.sauve_ImaPGM("imadec1.pgm");
      imatemp.mise_a_zero();
      for (i=0; i<nblig; i++) {
        i0=maxi(i-nl2,0); i2=mini(i+nl2,nblig-1);
        for (j=0; j<nbcol; j++) {
          j0=maxi(j-nc2,0); j2=mini(j+nc2,nbcol-1);
          sn=n=0;
          for (ii=i0; ii<=i2; ii++)
            for (jj=j0; jj<=j2; jj++) {
//						n+=imadec0(ii,jj)+imadec1(ii,jj); 
              n+=maxi(imadec0(ii,jj),imadec1(ii,jj)); 
              sn++;
            }
          xx=(double)n/2.0/sn;
          for (ii=i0; ii<=i2; ii++)
            for (jj=j0; jj<=j2; jj++) imatemp(ii,jj)=maxi((float)xx,imatemp(ii,jj));
          if (xx>0.5)
            for (ii=i0; ii<=i2; ii++)
              for (jj=j0; jj<=j2; jj++) imaRes(ii,jj)=1;
        }
      }
      imadec0.sauve_ImaPGM("imadec0.pgm"); imadec1.sauve_ImaPGM("imadec1.pgm");
      imatemp.imaunsignedchar(1).sauve_ImaPGM("imatemp.pgm"); imaRes.sauve_ImaPGM("imaRes.pgm");
/************************* fin a completer *************************/
      ima1=ima2;
    }
    return true;
  }
  static bool op_pb5(int nbp, Param_Type* p)
  {
    const int rayonS3=3;
    const BOOL itps=0;
    const BYTE ncl=2;
    const double beta=2000.;
    fichimage_entree ifich(p[0].p_fname.fname);
    cout<<" ouverture de l'image couleur\n";
    imacouleur<BYTE> imacol((imadata<BYTE>)ifich.LoadPGM());
    cout<<" conversion dans l'espace IST\n";
    imacol.RVB2IST();
    cout<<" selection de la composante I\n";
    imadata<float> imaI=(imadata<float>)imacol.select_composante(0); imaI.imaunsignedchar().sauve_ImaPGM("imaI.pgm");
    eltstruct S3(3,3);
    cout<<" erosion fonctionnelle avec elt structurant 8-connexe de rayon "<<rayonS3<<"\n";
    imadata<float> imaE=imaI.erode(S3,rayonS3); imaE.imaunsignedchar().sauve_ImaPGM("imaE.pgm");
    imadata<float> imaR;
    if (itps==0) {
      cout<<" reconstruction geodesique fonctionnelle\n";
      imaR=imaI.reconst_geod(imaE,S3); imaR.imaunsignedchar().sauve_ImaPGM("imaR.pgm");
    }
    else {
      cout<<" dilatation fonctionnelle\n";
      imaR=imaE.dilate(S3,rayonS3);
    }
    imadata<float> imaD=imaI-imaR; imaD.imaunsignedchar().sauve_ImaPGM("imaD.pgm");
    statclasses clkmeans;
    if (itps==0) {
      cout<<" classification c-means en "<<(int)ncl<<" classes\n";
      sampleset<pixel<float> > datset=imaD.dataset();
      clkmeans=datset.k_means(ncl); clkmeans.affiche();
    }
    else {
      stat1class cl(1,1);
      cl.cova(0,0)=1.f; cl.icova(0,0)=1.f/cl.cova(0,0);
      cl.label()=1; cl.mean(0)=10.f; clkmeans.ajoute(cl);
      cl.label()=2; cl.mean(0)=150.f; clkmeans.ajoute(cl);
      if (ncl==3) {cl.label()=3; cl.mean(0)=50.f; clkmeans.ajoute(cl);}
    }
    clkmeans.affiche();
    cout<<" classification markovienne en "<<(int)ncl<<" classes avec beta = "<<beta<<"\n";
    imalabels imaDcl(imaD,clkmeans,(float)beta,"ICM",0.01f);
    imadata<BYTE> imaDc=imaDcl.conv2imBYTE(1); imaDc.sauve_ImaPGM("imaDc.pgm");
    stat1class clas;
    BYTE lab=1, k;
    double max=clkmeans.extract(lab).mean(0);
    cout<<" classe "<<(int)lab<<" moy = "<<max<<"\n";
    for (k=2; k<=ncl; k++) {
      clas=clkmeans.extract(k); cout<<" classe "<<(int)k<<" moy = "<<clas.mean(0)<<"\n";
      if (clas.mean(0)>max) {
        max=clas.mean(0); lab=k;
      }
    }
    float val=(int)lab*255.f/ncl;
    cout<<" classe des lignes blanches de label "<<(int)lab<<" val = "<<val<<"\n";
    imabin imab=(imabin(imaDc,val-1.5f))&&(imabin(imaDc,val+1.5f).negatif()); imab.imaunsignedchar().sauve_ImaPGM("imabin_Dc.pgm");
    if (itps==0) {
      int ncc=0;
      imadata<int> imacc=imab.composantes_connexes(ncc,8); imacc.imaunsignedchar().sauve_ImaPGM("ima_Dc_cc.pgm");
//		imaregions imacc2(imab); imacc2.conv2imBYTE().sauve_ImaPGM("ima_Dc_cc2.pgm");
      imaobjets imaobj(imacc.imaunsignedchar(),0,0); imaobj.affiche();
      imaobj.regroupe_alignements(0); imaobj.affiche();
      imaobj.conv2imBYTE().sauve_ImaPGM("ima_Dc_cc2.pgm");
    } else {
      imadata<float> ima_c=imadata<float>(imab);
      cout<<" detection de contours (Sobel)\n";
      imacontours imacont(ima_c,"Sobel",1,1);
      imadata<BYTE> ima_Dcc=imacont.conv2imBYTE(); ima_Dcc.sauve_ImaPGM("imaDcc.pgm");
      cout<<" transformee de Hough\n";
      imadata<int> imhough=imabin(ima_Dcc,1).hough_transform();
      imhough.statbasic(1);
      fichimage_entree ific1("reconstructionHough.pgm");
      imabin imaReconsH((imadata<float>)ific1.LoadPGM(),1);
      cout<<" fermeture sur l'image reconstruite de Hough\n";
      imaReconsH=imaReconsH.dilate(S3,2).erode(S3,3);
      imaReconsH.imaunsignedchar().sauve_ImaPGM("reconstructionHough_F.pgm");
//		(imab&&imaReconsH).imaunsignedchar().sauve_ImaPGM(p[1].p_fname.fname);
/*		cout<<" squelette de la fermeture sur l'image reconstruite de Hough\n";
	    imaReconsH=imaReconsH.squelette (8); 
      cout<<" reconstruction geodesique de la classe des lignes blanches a partir du squelette precedent\n";
      imab.reconstruction_geodesique(imaReconsH,8).imaunsignedchar().sauve_ImaPGM(p[1].p_fname.fname);*/
    }
    return true;
  }
  static bool op_pb4(int nbp, Param_Type* p)
  {
    eltstruct S3(3,3); 
    fichimage_entree ific1(p[0].p_fname.fname), ific2(p[1].p_fname.fname);
    imadata<float> imadon=(imadata<float>)ific1.LoadPGM();
    imabin imaROI((imadata<float>)ific2.LoadPGM(),1);
    int nblig=imaROI.nlig(), nbcol=imaROI.ncol(), nbpix=(int)imaROI.norm();
// Erosion fonctionnnelle (parcellez16_007.pgm, 3x3, 8, parcellez16_007erod.pgm) ;
    imadata<float> imadonE=imadon.erode(S3,p[2].p_sel.sel_cur+1);
// Reconstruction géodésique fonctionnnelle (parcellez16_007.pgm, parcellez16_007erod.pgm, 8, parcellez16_007erodRG.pgm) ;
    imadata<float> imadonRG=imadon.reconst_geod(imadonE,S3);
// + ou - entre images niveaux de gris (parcellez16_007.pgm, parcellez16_007erodRG.pgm, -, sillonsz16_007.pgm) ;
    imadata<float> imadat=imadon-imadonRG;
// Percentile supérieur (sillonsz16_007.pgm, masqparz16_007.pgm, 0.25, si_pct25z16_007.pgm) ;
    float s=imadat.seuil_percentile(imaROI,1.f-(p[3].p_float.val)/100.f);
    imabin imab(imadat,s);
    imadata<BYTE>(imab).sauve_ImaPGM(p[5].p_fname.fname);
// Transformée de Hough (si_pct25z16_007.pgm, Hgsi_pct25z16_007.pgm) ;
    imadata<float> imhough=imab.transformee_Hough();
//	imhough.statbasic(1); imhough.imaunsignedchar().sauve_ImaPGM(p[6].p_fname.fname);
// Transformée de Hough (masqparz16_007.pgm, Hgmsqparz16_007.pgm) ;
    imadata<float> imhoughP=imaROI.transformee_Hough();
    imhoughP.statbasic(); 
    float x_scal=255.f/(float)imhoughP.maxI(); x_scal=mini(1.f,x_scal);
    cout<<" parametre x_scal = "<<x_scal<<"\n";
// Supérieur au seuil (Hgmsqparz16_007.pgm, 20, MqHgparz16_007.pgm) ;
//  imabin binMhough(imhoughP,(p[6].p_float.val)/x_scal);
    s=0.5f*(float)nbpix/pow((float)nblig*nblig+nbcol*nbcol,0.5f); cout<<" seuil propose pour la TH du masque de la parcelle = "<<s<<"\n";
    imabin binMhough(imhoughP,s);
// Application d'un masque (Hgsi_pct25z16_007.pgm, MqHgparz16_007.pgm, Hg2si_pct25z16_007.pgm) ;
    int i,j,nbligH=mini(imhough.nlig(),binMhough.nlig()),nbcolH=mini(imhough.ncol(),binMhough.ncol()); 
    for (i=0; i<nbligH; i++)
      for (j=0; j<nbcolH; j++) 
        if (binMhough(i,j)==0) imhough(i,j)=0;
//+ ou * avec 1 scalaire (Hg2si_pct25z16_007.pgm, 255 / 507, *, Hg2si_pct25z16_007.pgm) ;
    if (x_scal!=1.f) imhough=imhough*x_scal; imhough.statbasic(1);
// + ou - entre images niveaux de gris (Hg2si_pct25z16_007.pgm, Hgmsqparz16_007.pgm, /, HgNsi_pct25z16_007.pgm) ;
    imadata<float> imhoughN=(imhough/imhoughP)*100.f; imhoughN.statbasic(1);
// Percentile supérieur (HgNsi_pct25z16_007.pgm, MqHgparz16_007.pgm, 0.25, HgBsi_pct25z16_007.pgm) ;
    s=imhoughN.seuil_percentile(binMhough,1.f-(p[3].p_float.val)/100.f);
    cout<<" seuil sur la transformee de Hough normalisee = "<<s<<"\n";
    imabin imabH(imhoughN,s);
    imadata<BYTE>(imabH).sauve_ImaPGM(p[6].p_fname.fname);
// Profil en ligne ou en colonne (HgBsi_pct25z16_007.pgm, lig, profil_007.pgm) ;
    int n=16;
    imadata<float> profil(n,nbcolH); profil.mise_a_zero();
    for (j=0; j<nbcolH; j++) {
      for (i=0; i<nbligH; i++) profil(0,j)+=imabH(i,j);
      for (i=1; i<n; i++) profil(i,j)=profil(0,j);
    }
    cout<<" profil :\n"; profil.statbasic(1); //profil.sauve_ImaPGM(p[6].p_fname.fname); 
//Erosion fonctionnnelle (profil_007.pgm, 11x11, 8, profil_007erod.pgm) ;
    int m=5;
    if (p[4].p_sel.sel_cur==1) m=2;
    if (p[4].p_sel.sel_cur==2) m=10;
    S3(0,0)=S3(0,1)=S3(0,2)=S3(2,0)=S3(2,1)=S3(2,2)=0; S3.affiche(); cout<<" a appliquer "<<m<<" fois\n";
    imadata<float> profilE=profil.erode(S3,m); cout<<" profil erode :\n"; profilE.statbasic(1);
//Reconstruction géodésique fonctionnnelle (profil_007.pgm, profil_007erod.pgm, 8, profil_007erodRG.pgm) ;
    imadata<float> profilRG=profil.reconst_geod(profilE,S3); cout<<" profil reconstruit :\n"; profilRG.statbasic(1);
//+ ou - entre images niveaux de gris (profil_007.pgm, profil_007erodRG.pgm, -, sharpprofil_007.pgm) ;
    imadata<float> sharp=profil-profilRG; sharp.statbasic(1);
    sharp.imaunsignedchar(1).sauve_ImaPGM(p[7].p_fname.fname);
    return true;
  }
  static bool op_pb6(int nbp, Param_Type* p)
	{ 
    eltstruct S3(3,3); 
    fichimage_entree ific1(p[0].p_fname.fname)/*, ific2(p[1].p_fname.fname)*/;
		int nblig=p[1].p_int.val, nbcol=p[2].p_int.val;
    imadata<float> imadon=(imadata<float>)ific1.LoadPGM();
		imadon.sousIma(nblig+p[3].p_int.val,nbcol+p[4].p_int.val,p[3].p_int.val,p[4].p_int.val,1,1);
//    imabin imaROI((imadata<float>)ific2.LoadPGM(),1);
//    int nblig=imaROI.nlig(), nbcol=imaROI.ncol(), nbpix=(int)imaROI.norm();
    imadon.imaunsignedchar(1).sauve_ImaPGM("sousImage_LIDAR.pgm");
		imadon.statbasic(1);
		float xminus=(float)pow(imadon.varI(),0.5);
//		imadata<float> imadon2=imadon-xminus;
		imadata<float> imadonRG=imadon.reconst_geod(imadon+(-xminus),S3);
    imadonRG.imaunsignedchar(1).sauve_ImaPGM("sousImage_LIDAR_RG.pgm");
    imadonRG.imaunsignedchar(1).sauve_ImaPGM(p[5].p_fname.fname);
    return true;
	}

  static bool op_pbEES4A(int nbp, Param_Type* p)
  {
    fichimage_entree ifich(p[0].p_fname.fname);
		imadata<float> ima=(imadata<float>)ifich.LoadPGM();

 //   imaobj.detection_disques(p[1].p_float.val,p[2].p_float.val).imaunsignedchar().sauve_ImaPGM(p[1].p_fname.fname);
    return true;
  }

  static Func_Type func[];
  static DWORD GetNbFunc();
  static const char* CatNames[];
  static DWORD GetNbCat();
  static int GetRealIdx(int idx_cat, int idx_fct);
};

#endif // _PARAMS_H