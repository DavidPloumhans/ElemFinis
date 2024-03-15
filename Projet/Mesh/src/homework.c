#include "fem.h"


double dist(double x, double y, double x0, double y0) {  // renvoie la distance entre deux points
    return sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0));
}

double geoSize(double x, double y){

    femGeo* theGeometry = geoGetGeometry();
    // faut faire un truc un peu opti
    double h1 = theGeometry->h1;
    double h2 = theGeometry->h2;
    double l1 = theGeometry->l1;
    double l2 = theGeometry->l2;
    double r = theGeometry->r;

    double disklx = l1 - r;
    double diskly = h2 - r;
    double diskrx = l1 + l2 + r;
    double diskry = h2 - r;
    
    double refDist = 3 * r;  // distance de référence par rapport aux centres des cercles
    double hstar = 0.15; // valeur par défault pour le maillage
    double toReturn = hstar; 
    double hd = hstar / 8.0;  // h disk
    
    // distances aux cercles
    double distl = dist(x, y, disklx, diskly);  // distance au centre du cercle de gauche
    double distr = dist(x, y, diskrx, diskry);  // distance au centre du cercle de droite
    // Triangulation de Delaunay
    if (distl < refDist) {
        double a = 2.0 * (hstar - hd) / (refDist * refDist * refDist);
        double b = 3.0 * (hd - hstar) / (refDist * refDist);
        double d = hd;
        toReturn = a * distl * distl * distl + b * distl * distl + d;
    }
    else if (distr < refDist) {
        double a = 2.0 * (hstar - hd) / (refDist * refDist * refDist);
        double b = 3.0 * (hd - hstar) / (refDist * refDist);
        double d = hd;
        toReturn = a * distr * distr * distr + b * 
        distr * distr + d;
    }
    if (y < h2) {
        toReturn = fmin(toReturn, hd * y * y + 0.03);
    }
    // printf("toReturn = %f\n", toReturn);
    // toReturn = hstar;
    return fabs(toReturn);
}


#define ___ 0

void geoMeshGenerate() {

    femGeo* theGeometry = geoGetGeometry();
    // récupérer les valeurs de theGeometry pour les utiliser dans la fonction de création de Mesh
    double h1 = theGeometry->h1;
    double h2 = theGeometry->h2;
    double l1 = theGeometry->l1;
    double l2 = theGeometry->l2;
    double r = theGeometry->r;

 //
//  -1- Construction de la g�om�trie avec OpenCascade
//
 
    int ierr;
    int idPlate = gmshModelOccAddRectangle(0, 0, 0, l2 + 2 * l1, h1 + h2, -1, 0, &ierr);   // c'est x, y et z lower left corner puis dx, dy, tag et rounded_radius
    ErrorGmsh(ierr);
    int idDiskl = gmshModelOccAddDisk(l1 - r, h2 - r, 0, r, r, -1, NULL,0,NULL,0,&ierr);  // c'est xc, yc, zc, rx, ry, tag
    ErrorGmsh(ierr);
    int idDiskr = gmshModelOccAddDisk(l1 + l2 + r, h2 - r, 0, r, r, -1, NULL,0,NULL,0,&ierr);  // c'est xc, yc, zc, rx, ry, tag
    ErrorGmsh(ierr);
    int idRectl1 = gmshModelOccAddRectangle(0, 0, 0, l1, h2-r, -1, 0, &ierr);   // c'est x, y et z lower left corner puis dx, dy, tag et rounded_radius
    ErrorGmsh(ierr);
    int idRectl2 = gmshModelOccAddRectangle(0, h2 - r, 0, l1-r, r, -1, 0, &ierr);   // c'est x, y et z lower left corner puis dx, dy, tag et rounded_radius
    ErrorGmsh(ierr);
    int idRectr1 = gmshModelOccAddRectangle(l1 + l2, 0, 0, l1, h2-r, -1, 0, &ierr);   // c'est x, y et z lower left corner puis dx, dy, tag et rounded_radius
    ErrorGmsh(ierr);
    int idRectr2 = gmshModelOccAddRectangle(l1 + l2 + r, h2 - r, 0, r, r, -1, 0, &ierr);   // c'est x, y et z lower left corner puis dx, dy, tag et rounded_radius
    ErrorGmsh(ierr);

    
    int plate[] = {2, idPlate};
    int Diskl[] = {2, idDiskl};
    int Diskr[] = {2, idDiskr};
    int Rectl1[] = {2, idRectl1};
    int Rectl2[] = {2, idRectl2};
    int Rectr1[] = {2, idRectr1};
    int Rectr2[] = {2, idRectr2};

    gmshModelOccCut(plate, 2, Diskl, 2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);   // c'est objectDimTags, objectDimTags_n, toolDimtags, toolDimtags_n
    ErrorGmsh(ierr);
    gmshModelOccCut(plate, 2, Diskr, 2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);   
    ErrorGmsh(ierr);
    gmshModelOccCut(plate, 2, Rectl1, 2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);  
    ErrorGmsh(ierr);
    gmshModelOccCut(plate, 2, Rectl2, 2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr); 
    ErrorGmsh(ierr);
    gmshModelOccCut(plate, 2, Rectr1, 2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);  
    ErrorGmsh(ierr);
    gmshModelOccCut(plate, 2, Rectr2, 2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);  
    ErrorGmsh(ierr);
    
 
//
//  -2- D�finition de la fonction callback pour la taille de r�f�rence
//      Synchronisation de OpenCascade avec gmsh
//      G�n�ration du maillage (avec l'option Mesh.SaveAll :-)
                  
   
    geoSetSizeCallback(geoSize);
                                  
    gmshModelOccSynchronize(&ierr);       
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr);  
       
//
//  Generation de quads :-)
//
//    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
//    gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
//    gmshOptionSetNumber("Mesh.Algorithm", 8, &ierr);  chk(ierr);
//    gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr);  chk(ierr);
//    gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  chk(ierr);
//    gmshModelMeshGenerate(2, &ierr);  
   
 
//
//  Plot of Fltk
//
//   gmshFltkInitialize(&ierr);
//   gmshFltkRun(&ierr);  chk(ierr);
//
    
}