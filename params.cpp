#include "params.h"
#include "params2.h"

Param_Type Param_Type::para[]={
  Param_Type(),
  Param_Type("Image en entrée",true),                                        									// 01
  Param_Type("Image en sortie",false),															               						// 02
  Param_Type("Seuil bas (0..255)",0.0f,255.0f,1.0f),											        						// 03
  Param_Type("Seuil haut (0..255)",0.0f,255.0f,1.0f),														  						// 04
  Param_Type("Ecart-type",0.000001f,100.0f,1.6f),                            									// 05
  Param_Type("Pourcentage",0.0f,1.0f,0.1f),                                 									// 06
  Param_Type("Nombre de classes",1,50,5),																					   					// 07
  Param_Type("Connexité",0,"4|8|"),																														// 08
  Param_Type("Voisinage",0,"3x3|5x5|7x7|9x9|11x11|13x13|15x15|"),															// 09
  Param_Type("Nombre d'iterations",1,10,1),																			  						// 0A
  Param_Type("Noyau transf. distance",0,"4-connexité|8-connexité|5x5|7x7|"),									// 0B
  Param_Type("Element structurant",0,"3x3|5x5|7x7|9x9|11x11|13x13|15x15|17x17|19x19|21x21"),	// 0C
  Param_Type("Image des marqueurs",true),							                               					// 0D
  Param_Type("Nb. de voisins",0,"1|3|5|7|9|"),								                      					// 0E
  Param_Type("Poids beta du terme Markovien",0.0f,99000.0f,0.0f),					           					// 0F 
  Param_Type("Algorithme",0,"ICM|SAG|SAM|"),																			  					// 10
  Param_Type("Noyau gradient",0,"Prewitt|Sobel|MDIF|"),                     									// 11
  Param_Type("Image de la norme du gradient",true),											             					// 12
  Param_Type("Image de la norme du gradient",false),														    					// 13
  Param_Type("Image de la direction du gradient",false),                    									// 14
  Param_Type("Long. prolongation en pixels",1,30,4),				                      						// 15
  Param_Type("Noyau du Laplacien",0,"4|8"),													                					// 16
  Param_Type("Algorithme",0,"Deriche|Shen-Castan|"),												        					// 17
  Param_Type("Paramètre",0.0f,10.0f,1.0f),																										// 18
  Param_Type("Image de contours",true),																      									// 19
  Param_Type("Ima. contours/points d'intérêt",false),															  					// 1A
  Param_Type("Version",0,"Normale|Inversée|"),											                					// 1B
  Param_Type("Enchaînement",0,"Ouv-Ferm|Ferm-Ouv|"),												        					// 1C
  Param_Type("Ordre",1,15,1),																																	// 1D
  Param_Type("Demi-longueur du noyau",1,20,1),                              									// 1E
  Param_Type("Seuil homogeneite",0.0f,10000.0f,0.0f),															  					// 1F 
  Param_Type("Selection germes",0,"histogramme|aleatoire|"),				                					// 20
  Param_Type("Niveau de depart (0=sommet pyramide)",0,10,0),										    					// 21
  Param_Type("Nombre de regions",0,255,8),																										// 22
  Param_Type("Nombre minimal de pixels par region",0,999,5),																	// 23
  Param_Type("Opération",0,"+|-|L1(-)|/|"),																										// 24
  Param_Type("Opération",0,"+|*|"),																														// 25
  Param_Type("Nombre scalaire",-10000.0f,10000.0f,0.0f),																			// 26
  Param_Type("Nom de la 1ère image",true),																										// 27
  Param_Type("Nbre d'images de la séquence",0,999,2),																					// 28
  Param_Type("Espace",0,"IST|LAB|I1I2I3|L1|"),                              									// 29
  Param_Type("Taille de fenêtre",5,31,7),																											// 2A
  Param_Type("N bins",32,256,64),																															// 2B
  Param_Type("Ordre k",1,4,2),																																// 2C
  Param_Type("Orientation",0,"0|1|2|3|"),																											// 2D
  Param_Type("Distance",1,15,1), 																															// 2E
  Param_Type("Paramètre de texture",1,"contraste|entropie|énergie|"),       									// 2F
  Param_Type("Valeur extraite",0,"entropie|contraste|dissimilarité|homogénéité|"),						// 30
  Param_Type("seuil en % pixels",0.0f,100.0f,5.0f),																						// 31
  Param_Type("Image bin. des regions d'intérêt",true),																				// 32
  Param_Type("Direction du profil",0,"lig.|col.|"),																						// 33
  Param_Type("Precision angulaire du profil",0,"+-5°|+-2°|+-10°"),														// 34
  Param_Type("Alpha",0.0f,1.0f,0.5f),																													// 35
  Param_Type("Gamma",1.41f,2.0f,2.0f),																												// 36
  Param_Type("Type de flot optique",0,"Horn&Schunck|.|.|"),																		// 37
	Param_Type("algorithme",0,"rapide| MM|"),																			    					// 38
	Param_Type("Seuil détecteur de Harris",0.f,20000.f,2000.f),																	// 39
	Param_Type("Coordonnées x ou y de la translation",-10000.f,10000.f,0.f),										// 3A
	Param_Type("Angle en degré de la rotation",-180.f,180.f,0.f),																// 3B
	Param_Type("Rapport de l'homothétie",0.1f,10.f,1.f),																				// 3C
	Param_Type("Dim. x ou y de la sous-image",0,10000,0),																				// 3D
	Param_Type("Coord. UL x ou y de la sous-image",0,1000,0),																		// 3E
	Param_Type("Sous-échantillonnage x ou y",1,4,1),																						// 3F
	Param_Type("Amplitude minimale des extrema",1.f,255.f,1.f),																	// 40
  Param_Type("Nombre de pixels minimum",1,500,20),		          															// 41
	Param_Type("Nombre de superpixels",10,10000,200),																						// 42
	Param_Type("Coef. norm. radiometrie vs distance spatiale",10.f,100.f,20.f),								  // 43
  Param_Type("Sens des extrema",0,"max|min|"),																								// 44
  Param_Type("detail ratio (paramètre edge drawing)",1,20,4),																	// 45
};


const char* Func_Type::CatNames[]={
  "00-opérations de base",
  "01-filtrage passe-bas",
  "02-MM binaire",
  "03-MM fonctionnelle",
  "04-classification",
  "05-contours",
  "06-segmentation",
  "07-objets",
  "08-couleur et texture",
  "09-mouvement",
  "10-reconnaissance"
};

// dans l'ordre: les entrées puis les sorties (codes hexa)
  Func_Type Func_Type::func[]={
	Func_Type(0,"+ ou - entre images niveaux de gris","Addition ou soustraction ou valeur absolue de la différence ou division entre valeurs des pixels de 2 images à niveaux de gris.","\x01\x01\x24\x02",op_operation_ima),
  Func_Type(0,"+ ou * avec 1 scalaire","Addition ou multiplication des valeurs des pixels par un scalaire.","\x01\x26\x25\x02",op_operation_scal),
  Func_Type(0,"Inversion de niveaux de gris","Inversion de l'échelle des niveaux de gris.","\x01\x02",op_inv_ngr),
  Func_Type(0,"Etirement de la dynamique","Création de l'image de dynamique étirée entre 0 et 255.","\x01\x02",op_etir_dyn),
  Func_Type(0,"Egalisation d'histogramme","Création de l'image d'histogramme égalisé avec ré-étirement de la dynamique.","\x01\x02",op_egal_histo),
  Func_Type(0,"Distance entre histogrammes","Distance entre histogrammes de 2 image en entrée.","\x01\x01",op_distHist),
  Func_Type(0,"Application d'un masque","Mise à 0 des pixels de l'image à niv. de gris (1ère im. en entrée) ayant la valeur 0 sur l'image binaire 'masque' (2ème im. en entrée).","\x01\x01\x02",op_masq_ima),
  Func_Type(0,"Supérieur au seuil","Création d'1 image binaire des valeurs supérieures au seuil.","\x01\x03\x02",op_bin_sup),
  Func_Type(0,"Inférieur au seuil","Création d'1 image binaire des valeurs inférieures au seuil.","\x01\x04\x02",op_bin_inf),
  Func_Type(0,"Seuillage à hystérésis","Création d'1 image bin. des valeurs sup. au seuil haut ou sup. au seuil bas ET dont la composante connexe a au moins 1 pixel de valeur sup. au seuil haut.",
            "\x01\x03\x04\x02",op_seuil_hysteresis),
	Func_Type(0,"Seuillage entre deux valeurs","Création d'1 image bin. des valeurs comprises entre le seuil haut et le seuil bas.",
            "\x01\x04\x03\x02",op_bin_inter),
  Func_Type(0,"Percentile supérieur","Création d'1 image binaire des p% pixels de valeurs les plus élevées.","\x01\x32\x06\x02",op_bin_pct_sup),
  Func_Type(0,"Percentile inférieur","Création d'1 image binaire des p% pixels de valeurs les plus faibles.","\x01\x32\x06\x02",op_bin_pct_inf),
  Func_Type(0,"Seuillage automatique","Création d'1 image binaire des valeurs supérieures au seuil automatique d'Otsu.","\x01\x02",op_bin_otsu),
  Func_Type(0,"Maxima image","Création d'1 image binaire des valeurs maximales d'1 image à niveaux de gris.","\x01\x02",op_bin_maxAbs),
  Func_Type(0,"Maxima réginaux","Création d'1 image binaire des maxima régionaux d'1 image à niveaux de gris.","\x01\x02",op_bin_maxReg),
  Func_Type(0,"Ajout bruit gaussien","Ajout d'1 bruit gaussien centré sur les valeurs des pixels.","\x01\x05\x02",op_gauss_noise),
  Func_Type(0,"Ajout bruit impulsif","Ajout d'1 bruit impulsif poivre et sel sur les valeurs des pixels.","\x01\x06\x02",op_impul_noise),
  Func_Type(1,"Filtre Gaussien","Filtrage linéaire à noyau gaussien.","\x01\x05\x02",op_filtr_Gauss),
  Func_Type(1,"Filtre moyenne","Filtrage linéaire à noyau identité.","\x01\x09\x02",op_filtr_moyen),
  Func_Type(1,"Filtre médian","Filtrage de rang correspondant à la valeur médiane.","\x01\x09\x02",op_filtr_median),
  Func_Type(1,"Filtre Nagao couplé avec moyenneur","Filtrage multi-noyaux de Nagao.","\x01\x09\x02",op_filtr_Nagao),
  Func_Type(1,"Filtre Nagao couplé avec médian","Filtrage multi-noyaux de Nagao.","\x01\x09\x02",op_filtr_NagaoMed),
  Func_Type(1,"Filtre SNN","Filtrage par filtre 'Symetric Nearest Neighbour'.","\x01\x09\x02",op_filtr_SNN),
  Func_Type(1,"Filtre FAS","Filtrage par Filtre Alterné Séquentiel (en morphologie mathématique fonctionnelle).","\x01\x1C\x08\x1D\x02",op_filtr_FAS),
  Func_Type(1,"Pic Signal Noise Ratio","Mesure de rapport signal à bruit. La première image en entrée est l'image reconstruite et la seconde l'image 'vérité'.","\x01\x01",op_psnr),
  Func_Type(2,"AND","ET entre valeurs des pixels de 2 images binaires.","\x01\x01\x02",op_and_bin),
  Func_Type(2,"XOR","OU exclusif entre valeurs des pixels de 2 images binaires.","\x01\x01\x02",op_xor_bin),
  Func_Type(2,"OR","OU non exclusif entre valeurs des pixels de 2 images binaires.","\x01\x01\x02",op_or_bin),
  Func_Type(2,"NOT","Complémentaire des valeurs d'1 image binaire.","\x01\x02",op_not_bin),
  Func_Type(2,"Transformée en distance","Création d'une image des niveaux de gris des distances du fond aux objets.","\x01\x0B\x02",op_tr_dist),
  Func_Type(2,"Transformée en distance géodésique","Image des distances géodésiques (masque donné par la deuxième image en entrée) du fond aux objets.","\x01\x01\x0B\x02",op_tr_dist_geo),
  Func_Type(2,"Erosion binaire","Application de l'opérateur d'érosion morphologique à 1 image binaire.","\x01\x0C\x08\x02",op_erode_bin),
  Func_Type(2,"Dilatation binaire","Application de l'opérateur de dilatation morphologique à 1 image binaire.","\x01\x0C\x08\x02",op_dilate_bin),
  Func_Type(2,"Ouverture binaire","Application de l'opérateur d'ouverture morphologique à 1 image binaire.","\x01\x0C\x08\x02",op_ouvre_bin),
  Func_Type(2,"Fermeture binaire","Application de l'opérateur de fermeture morphologique à 1 image binaire.","\x01\x0C\x08\x02",op_ferme_bin),
  Func_Type(2,"Tophat binaire","Application de l'opérateur dit 'chapeau haute forme' de morphologie mathématique à 1 image binaire.","\x01\x0C\x08\x1B\x02",op_tophat_bin),
  Func_Type(2,"Reconstruction géodésique","Reconstruction géodésique d'1 image binaire (1ère image en entrée) à partir d'1 image de marqueurs (2nde image en entrée).","\x01\x01\x08\x02",op_reconst_geod_bin),
//  Func_Type(2,"Etiquetage en composantes connexes","Etiquetage en composantes connexes d'1 image binaire.","\x01\x08\x02",op_etiquette_cc),
	Func_Type(2,"Etiquetage en composantes connexes","Etiquetage en composantes connexes d'1 image binaire.","\x01\x08\x38\x02",op_etiquette_cc),
  Func_Type(2,"Elimination des 'petites' composantes","Elimination des composantes connexes de nombre de pixels inférieur à un seuil sur 1 image binaire.","\x01\x41\x02",op_elimin_cc),
	Func_Type(2,"Décomposition en rectangles","Decomposition en rectangles d'1 image binaire.","\x01\x02",op_etiquette_boxes),
  Func_Type(2,"Erosion ultime","Création de l'image des érodés ultimes des objets d'1 image binaire.","\x01\x08\x02",op_erode_ultime),
  Func_Type(2,"Détection de coins","Détection des coins carrés UL (Upper left), UR (Upper right), LL (Lower left) et LR (Lower right) des objets d'1 image binaire.","\x01\x1E\x02",op_detect_coin),
  Func_Type(2,"Enveloppe convexe","Création de l'image des enveloppes convexes des objets d'une image binaire.","\x01\x02",op_env_convexe),
  Func_Type(2,"Squelette","Création de l'image des squelettes des objets d'une image binaire","\x01\x08\x02",op_squelette),
  Func_Type(2,"Elagage/ébardage","Application de la transformation en tout ou rien réalisant l'élagage ou ébardage notamment de squelettes.","\x01\x0A\x02",op_elagage),
  Func_Type(2,"Zones d'influence géodésique","Création de l'image binaire des zones d'influence géodésique des objets d'1 image de marqueurs (1ère image en entrée) dans 1 image binaire (2nde image en entrée).","\x01\x01\x08\x02",op_zon_infl_geod),
  Func_Type(2,"Elimination objets du bord","Elimination des objets touchant le bord de l'image binaire.","\x01\x08\x02",op_elim_objet_bord),
  Func_Type(2,"Bouchage de trous","Bouchage des trous sur l'image binaire.","\x01\x08\x02",op_bouche_trou),
  Func_Type(3,"Erosion fonctionnelle","Application de l'opérateur d'érosion morphologique à 1 image de niveaux de gris.","\x01\x0C\x08\x02",op_erode_fct),
  Func_Type(3,"Dilatation fonctionnelle","Application de l'opérateur de dilatation morphologique à 1 image de niveaux de gris.","\x01\x0C\x08\x02",op_dilate_fct),
  Func_Type(3,"Rehaussement de contraste","Opérateur de rehaussement de contraste de morphologie mathématique de paramètre 'alpha'.","\x01\x08\x0C\x35\x02",op_rehaus_contra),
  Func_Type(3,"Ouverture fonctionnelle","Application de l'opérateur d'ouverture morphologique à 1 image de niveaux de gris.","\x01\x0C\x08\x02",op_ouvre_fct),
  Func_Type(3,"Fermeture fonctionnelle","Application de l'opérateur de fermeture morphologique à 1 image de niveaux de gris.","\x01\x0C\x08\x02",op_ferme_fct),
  Func_Type(3,"Top hat et top hat inversé","Application de l'opérateur de top hat à 1 image de niveaux de gris.","\x01\x0C\x08\x1B\x02",op_top_hat),
  Func_Type(3,"Top hat par érosion reconstruction","Variante de l'opérateur de top hat pour 1 image de niveaux de gris.","\x01\x0C\x08\x1B\x02",op_top_hat_2),
  Func_Type(3,"Reconstruction géodésique fonctionnelle","Reconstruction géodésique d'1 image de niveaux de gris (1ère image en entrée) à partir d'1 image de marqueurs (2nde image en entrée).","\x01\x01\x08\x02",op_reconst_geod_fct),
//  Func_Type(3,"Transformation h-max","Elimination des maxima régionaux d'amplitude inférieure à h.","\x01\x40\x08\x02",op_h_max),
//  Func_Type(3,"Transformation h-min","Elimination des minima régionaux d'amplitude inférieure à h.","\x01\x40\x08\x02",op_h_min),
  Func_Type(3,"Transformation h-max ou h-min","Elimination des extrema (maxima ou minima) régionaux d'amplitude inférieure à h.","\x01\x44\x40\x08\x02",op_h_max_min),
  Func_Type(3,"Maxima régionaux","Detection des maxima régionaux d'amplitude supérieure à h.","\x01\x40\x08\x02",op_max_reg),
  Func_Type(4,"Classification aveugle aux plus proches voisins","Classification supervisée non paramétrique, dite 'k-ppv'. Le 1er fichier en entrée est l'image à classifier, le 2nd est le fichier texte des échantillons d'apprentissage (format : lig., col., label).",
			"\x01\x01\x0E\x02",op_class_ppv),
  Func_Type(4,"Classification aveugle des c-moyennes","Classification non supervisée avec en entrée l'image des données, en sortie l'image des labels. Les caractéristiques des classes sont affichées à la console.",
			"\x01\x07\x02",op_class_cmeans),
  Func_Type(4,"Classification des c-moyennes sur ROI","En entrée : l'image des données et l'image binaire des pixels à considérer, en sortie : l'image des labels et les caractéristiques des classes à la console.",
			"\x01\x32\x07\x02",op_class_cmeans_mask),
  Func_Type(4,"Classification supervisée MRF","Classification basée sur un modèle markovien de l'image des labels. Le 1er fichier en entrée est l'image à classifier, le 2nd est le fichier texte des classes.","\x01\x01\x0F\x10\x02",op_class_mrf),
  Func_Type(4,"Classification supervisée MRF avec processus lignes","Classification basée sur un modèle markovien de l'image des labels. Le 1er fichier en entrée est l'image à classifier, le 2nd est le fichier texte des classes.","\x01\x01\x0F\x10\x02",op_class_mrf_lines),
//	Func_Type(4,"Classification EM Gibssien","Description...","\x01\x02",op_class_emgibbs),
  Func_Type(4,"Longueur des frontières","Image des frontières d'image des labels en entrée.","\x01\x08\x02",op_bord_class),
//  Func_Type(5,"Gradient de l'image","Description...","\x01\x11\x13\x14",op_image_gradient),
  Func_Type(5,"Gradient de l'image","Création de l'image de la norme du gradient par filtrage linéaire passe-haut (Prewitt, Sobel ou MDIF).","\x01\x11\x13",op_image_gradient),
  Func_Type(5,"Détection de contours par gradient","Création de l'image des contours par seuillage de la norme du gradient (obtenu par filtrage linéaire Prewitt, Sobel ou MDIF).","\x01\x11\x15\x1A",op_cont_gradient),
  Func_Type(5,"Détection de contours par laplacien","Création de l'image des contours par détection des passages par zéro du Laplacien (obtenu par filtrage linéaire) et seuillage de la norme du gradient.","\x01\x11\x16\x1A",op_cont_laplacien),
  Func_Type(5,"Détection de contours par filtrage optimal","Création de l'image des contours par seuillage de optimal selon soit le filtre de Canny-Deriche soit celui de Shen-Castan.","\x01\x17\x18\x1A",op_cont_optimal),
  Func_Type(5,"Gradient morphologique","Application de l'opérateur de gradient morphologique (différence entre les résultats de la dilatation et de l'érosion fonctionnelles) en 4 ou 8 connexité.","\x01\x08\x02",op_grad_mm),
  Func_Type(5,"Gradient morphologique multi-échelles","Application de l'opérateur de gradient morphologique multi-échelles (érodé du top-hat du gradient à l'ordre(échelle) choisi(e)) en 4 ou 8 connexité.","\x01\x0C\x08\x02",op_grad_multiechmm),
	Func_Type(5,"Laplacien morphologique","Application de l'opérateur de laplacien morphologique en 4 ou 8 connexité.","\x01\x08\x02",op_lapl_mm),
  Func_Type(5,"Affinement de contours","Sélection des maxima locaux dans la direction du gradient. La 1ère image en entrée est image de contours et la 2nde l'image de données initiale.","\x19\x01\x1A",op_affine_contours),
//  Func_Type(5,"Prolongation de contours","Prolongation des contours sur une distance de N pixels. La 1ère image en entrée est 1 image de contours et la 2nde est 1 image de la norme du gradient.","\x19\x12\x15\x03\x1A",op_prolonge_contours),
  Func_Type(5,"Prolongation de contours","Prolongation des contours sur une distance de N pixels. La 1ère image en entrée est 1 image de contours et la 2nde est 1 image de la norme du gradient.","\x19\x12\x15\x03\x45\x1A",op_prolonge_contours),
  Func_Type(5,"Amélioration de contours","Méthode Edge Drawing sur une distance de N pixels. La 1ère image en entrée est 1 image de contours et la 2nde est l'image des données initiale.","\x19\x01\x15\x45\x1A",op_edge_drawing),
  Func_Type(5,"Transformée de Hough","Création de l'image de la transformée de Hough pour la recherche de droites dans l'espace des coordonnées polaires (rho,theta).","\x01\x02",op_hough),
  Func_Type(5,"Reconstruction à partir de la transformée de Hough","Création de l'image reconstruite à partir des maxima de la transformée de Hough. La deuxième image en entrée est l'image qui sert à donner les dimensions de l'image reconstruite.","\x01\x01\x02",op_rec_hough),
  Func_Type(5,"Profil en ligne ou en colonne","Intégrale sur la direction orthogonale à celle du profil (ligne ou colonne) d'une image en entrée binaire.","\x01\x33\x02",op_profil_bin),  
//  Func_Type("Analyse Hough","Description...","\x01\x02",op_analyse_hough),
  Func_Type(5,"Transformée de Hough sur les cercles","Création de l'image de la transformée de Hough pour la recherche de cercles.","\x01\x02",op_hough_cercles),
//  Func_Type(5,"Détection de points d'intérêt","Création de l'image des points d'intérêt par détecteur de Harris.","\x01\x36\x1A",op_points_interet),
  Func_Type(5,"Détection de points d'intérêt","Création de l'image des points d'intérêt par détecteur de Harris.","\x01\x05\x39\x1A",op_points_interet),
  Func_Type(5,"Mise en correspondance de points d'intérêt","Les deux images en entrée sont des images binaires des points d'intérêt.","\x01\x01\x1A",op_corresp_points_interet),
  Func_Type(6,"Segmentation à partir d'1 classification","Création de l'image des régions correspondant à l'image des labels en entrée.","\x01\x02",op_segm_classif),
  Func_Type(6,"Segmentation par croisement de 2 segmentations","Création de l'image des régions construite à partir des 2 segmentations en entrée.","\x01\x01\x02",op_croise_2segment),
  Func_Type(6,"Segmentation par croissance de régions","Création de l'image des régions par croissance de régions à partir de germes (aléatoiremt ou déduits de l'histo.), jusqu'à invalidation du critère d'homogénéité.","\x01\x1F\x20\x02",op_reg_growing),
  Func_Type(6,"Segmentation contrainte quadtree","Création de l'image des régions correspondant à une représentation de type quadtree.","\x01\x1F\x21\x02",op_reg_quadtree),
	Func_Type(6,"Segmentation par fusion de régions","Segmentation par fusion de régions dans un graphe avec sélection des arêtes de moindre coût par accord mutuel.","\x01\x1F\x02",op_fusion_reg),
  Func_Type(6,"Ligne de partage des eaux","Segmentation par l'algo. de la ligne de partage des eaux de morpho. math. fonctionnelle, ayant en entrée une image de la norme du gradient.","\x01\x23\x02",op_reg_lpe),
	Func_Type(6,"Segmentation selon graphe","Segmentation par fusion de régions dans un graphe avec sélection des arêtes de moindre coût par accord mutuel.","\x01\x22\x02",op_reg_graphe),
//  Func_Type(6,"Fusion de régions selon Mumford&Shah","Fusion de régions dans un graphe selon la fonctionnelle de Mumford et Shah. La 1ère image en entrée est 1 image de données et la 2nde est 1 image de régions.","\x01\x01\x22\x02",op_fus_reg_mumford),
	Func_Type(6,"Segmentation selon Mumford&Shah","Segmentation selon la résolution hiérarchique de Koepfler de la fonctionnelle de Mumford et Shah.","\x01\x22\x02",op_reg_mumford),
	Func_Type(6,"Superpixels SLIC","Super pixels selon l'algorithme Simple Linear Iterative Clustering.","\x01\x42\x43\x02",op_suppix_SLIC),
	Func_Type(6,"Waterpixels","Super pixels selon l'algorithme de Ligne de Partage des Eaux.","\x01\x42\x43\x08\x02",op_waterpix),
	Func_Type(7,"Paramètres de base des objet(s)","Barycentre, surface et boite englobante des objets. L'image des objets est en entrée, les résultats sont affichés à la console.","\x01",op_param_objets),
  Func_Type(7,"Classification selon taille","Classification des objets en fonction de leur taille.","\x01\x02",op_class_tailleCC),
	Func_Type(7,"Détection de disques","Détection des disques dans 1 image de régions (d'objets) avec en entrée les paramètres seuils sur les attributs de forme.","\x01\x18\x18\x02",op_detect_disques),
  Func_Type(8,"Conversion couleur","Conversion d'1 image RGB sous format ppm dans 1 espace couleur IST, LAB, I1I2I3 ou rgb (norme L1).","\x01\x29\x02",op_conv_couleur),
  Func_Type(8,"Composantes","Enregistre dans 3 images pgm chaque plan image couleur (e.g. pour des traitements marginaux).","\x01\x02\x02\x02",op_split_components),
  Func_Type(8,"Superposition","Superposition de 3 images pgm dans chaque plan image couleur.","\x01\x01\x01\x02",op_superpos_can),
  Func_Type(8,"Conversion en fausses couleurs","Conversion d'1 image niv. gris sous format ppm en fausses couleurs.","\x01\x02",op_conv_false_color),
//  Func_Type(8,"Image de Contraste","En entrée, une image ppm. ","\x01\x2A\x2B\x02",op_contraste),
//  Func_Type(8,"Image d'Entropie","En entrée, une image ppm. ","\x01\x2A\x2B\x02",op_entropy),
//  Func_Type(8,"Image d'Energie","En entrée, une image ppm. ","\x01\x2A\x2B\x02",op_energy),
//  Func_Type(8,"Image de moment centré d'ordre k","En entrée, une image ppm.","\x01\x2A\x2B\x2C\x02",op_moment),
//  Func_Type(8,"Image de paramètre de coocurrence","En entrée, une image binaire. Orientation: 0(-), 1(\\), 2(|), 3(/).","\x01\x2A\x2D\x2E\x02",op_cooc_binaire),  
//  Func_Type(8,"Image de Contraste","Création d'une image du paramètre de texture 'contraste' à partir d'une image à niveaux de gris.","\x01\x2A\x02",op_contraste),
//  Func_Type(8,"Image d'Entropie","En entrée, une image ppm. ","\x01\x2A\x02",op_entropy),
//  Func_Type(8,"Image d'Energie","En entrée, une image ppm. ","\x01\x2A\x02",op_energy),
  Func_Type(8,"Image de texture ordre 1","Création d'1 image du paramètre de texture sélectionné (contraste, énergie, entropie) à partir d'1 image à niveaux de gris.","\x01\x2A\x2F\x02",op_texture1),
  Func_Type(8,"Image de moment ordre k","Création d'une image du moment centré d'ordre k à partir d'1 image à niveaux de gris.","\x01\x2A\x2C\x02",op_moment),
  Func_Type(8,"Image de coocurrence (sur image niv. gris)","Création d'1 ima. à niv. de gris du param. calculé sur la mat. de cooccurrence (orient. 0(-), 1(\\), 2(|), ou 3(/) et dist. en entrée).","\x01\x2A\x2D\x2E\x30\x02",op_cooc_gris),
  Func_Type(8,"Image de coocurrence (sur image niv. gris)","Création d'1 ima. à niv. de gris du param. calculé sur la mat. de cooccurrence (largeur de fenêtre en entrée).","\x01\x2A\x30\x02",op_cooc_all_gris),
  Func_Type(8,"Matrice de coocurrence (sur image niv. gris)","Création d'1 ima. à niv. de gris de la mat. de cooccurrence (décalage max. en entrée).","\x01\x2A\x02",op_coocurence),
  Func_Type(8,"Image de coocurrence (sur image binaire)","A partir d'1 ima. bin. en entrée, création d'1 ima. à niv. de gris du param. calculé sur la mat. de cooccurrence (2x2, orient. 0(-), 1(\\), 2(|), ou 3(/) et dist. en entrée).","\x01\x2A\x2D\x2E\x30\x02",op_cooc_binaire),
	Func_Type(9,"Sous-image d'1 image","Découpe une sous-image dans l'image en entrée.","\x01\x3D\x3D\x3E\x3E\x3F\x3F\x02",op_sous_ima),
	Func_Type(9,"Projection à l'ordre 1 d'1 image","Projection de l'image à partir des paramètres de translation, rotation et homothétie en entrée.","\x01\x3A\x3A\x3B\x3C\x02",op_transfo_geom),
	Func_Type(9,"Transformation globale à l'ordre 0 entre 2 images","Estimation d'1 translation entre 2 images par correlation de phase et projection de l'image 2 dans géométrie de l'image 1.","\x01\x01\x02",op_rectif_corrphase),
	Func_Type(9,"Transformation globale à l'ordre 1 entre 2 images","Estimation d'1 transformation globale entre 2 images par Fourrier-Mellin et projection de l'image 2 dans géométrie de l'image 1.","\x01\x01\x02",op_rectif_FourrierMellin),
//  Func_Type(9,"Rectification d'1 image","Rectification pour qu'elle soit carree d'1 image couleur avec fond vert.","\x01\x02",op_rectif_carre_imacol),
  Func_Type(9,"Image de flot optique","A partir de 2 images à niveaux de gris en entrée, création d'1 image 2 canaux du flot optique (composantes H et V).","\x01\x01\x0F\x37\x02",op_flot_optique),  
//  Func_Type(9,"Fonction-Test à 2 entrées et 2 sorties","Description...","\x01\x01\x02\x02",op_test_2_2),
//  Func_Type(10,"Pb 1 : Détection route et transformation géométrique","En entrée, on a 2 images de données, et en sortie une image d'objets et 1 affichage à la console de la transformation.","\x01\x01\x02",op_pb1),
//  Func_Type(10,"Pb 1 : Détection de plot pour e-puck","En entrée : image de données, paramètres de seuillage et élément structurant et en sortie : affichage à la console des caractéristiques géom. du plot.","\x01\x0C\x02",op_pb1),
  Func_Type(10,"Pb 1 : Détection de plot pour e-puck","En entrée : image de données & fichier .txt des paramètres, et en sortie : détection du plot et commande.","\x01\x01\x02\x02",op_pb1),  
  Func_Type(10,"Pb 1b: Détection de plot pour e-puck","En entrée : fichier .txt des paramètres, et en sortie : détection du plot.","\x01\x02",op_pb1bis),  
  Func_Type(10,"Pb 2 : Détection de panneaux routiers","En entrée : image de données, paramètres de seuillage et élément structurant et en sortie : affichage à la console des caractéristiques géom. du plot.","\x01\x04\x03\x04\x0C\x02",op_pb2),
  Func_Type(10,"Pb 3 : Détection et suivi des objets","En entrée, on a le nom générique des images de données et le nombre d'images de la séquence. Les positions successives des objets sont affichées à la console et forment 1 image segmentee.","\x27\x28\x02",op_pb3),
  Func_Type(10,"Pb 4 : Détection des sillons","En entrée : (1) im. données, (2) im. masque de la parcelle. En sortie : (1) im. bin. des sillons, (2) im. bin. Hough normalisée et seuillée, (3) im. étroitesse profil.","\x01\x01\x0C\x31\x34\x02\x02\02",op_pb4),
  Func_Type(10,"Pb 5 : Détection de lignes blanches","En entrée, l'image de données","\x01\x02",op_pb5),
  Func_Type(10,"Pb 6 : Segmentation image LIDAR","En entrée, l'image de données","\x01\x3D\x3D\x3E\x3E\x02",op_pb6),
//	Func_Type(10,"Pb X : Gestalt cluster","En entrée, l'image de données","\x01\x02",detecte_cluster_2),
  Func_Type(10,"Pb miniprojet","????????","\x01\x02",op_pbEES4A),
};

DWORD Func_Type::GetNbFunc() { return sizeof(func)/sizeof(func[0]);}
DWORD Func_Type::GetNbCat() { return sizeof(CatNames)/sizeof(char*);}
int Func_Type::GetRealIdx(int idx_cat, int idx_fct)
{
  if(idx_cat>=(int)GetNbCat()) return idx_fct;
  for(int i=0, loc_idx=0; i<(int)GetNbFunc(); i++)
  {
    if(func[i].cat==idx_cat)
    {
      if(loc_idx==idx_fct) return i;
      else loc_idx++;
    }
  }
  return -1; // cela ne doit jamais arriver !
}

//HANDLE MyStr::hcons=0;
//bool MyStr::done=MyStr::ChangeCout();

void main()
{
	WinMain(::GetModuleHandle(0),0,"",SW_SHOWNORMAL);
}
