#include "fem.h"


double geoSize(double x, double y){

    femGeo* theGeometry = geoGetGeometry();
    
    double h = theGeometry->h;
    double toReturn = h;
    double x0 = theGeometry->xNotch;
    double y0 = theGeometry->yNotch;
    double r0 = theGeometry->rNotch;
    double h0 = theGeometry->hNotch;
    double d0 = theGeometry->dNotch;
  
    
    double x1 = theGeometry->xHole;
    double y1 = theGeometry->yHole;
    double r1 = theGeometry->rHole;
    double h1 = theGeometry->hHole;
    double d1 = theGeometry->dHole;

    double dist_notch = sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0)) - r0;
    double dist_hole = sqrt((x-x1)*(x-x1) + (y-y1)*(y-y1)) - r1;

    if (dist_notch < d0 && dist_hole > d1) {
        double a = 2.0 * (h0 - h) / (d0*d0*d0);
        double b = 3.0 * (h - h0) / (d0*d0);
        double d = h0;
        toReturn = (a*dist_notch*dist_notch*dist_notch + b*dist_notch*dist_notch + d) * 5;
    }
    else if (dist_hole < d1 && dist_notch > d0) {
        double a = 2.0 * (h1 - h) / (d1*d1*d1);
        double b = 3.0 * (h - h1) / (d1*d1);
        double d = h1;
        toReturn = (a*dist_hole*dist_hole*dist_hole + b*dist_hole*dist_hole + d)*2;
    }
    else if (dist_hole < d1 && dist_notch < d0) {
        double a = 2.0 * (h0-h) / (d0*d0*d0);
        double b = 3.0 * (h - h0) / (d0*d0);
        double d = h0;
        double toReturn1 = a*dist_notch*dist_notch*dist_notch + b*dist_notch*dist_notch + d;

        a = 2.0 * (h1-h) / (d1*d1*d1);
        b = 3.0 * (h - h1) / (d1*d1);
        d = h1;
        double toReturn2 = a*dist_hole*dist_hole*dist_hole + b*dist_hole*dist_hole + d;
        if (toReturn1 < toReturn2) {
            toReturn = toReturn1;
        } else {
            toReturn = toReturn2;
        }
    }
    else {
        toReturn = h * 5;
    }

    return toReturn;
    
//   
// Your contribution ends here :-)
//

}


#define ___ 0

void geoMeshGenerate() {

    femGeo* theGeometry = geoGetGeometry();

    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;
     
    double x0 = theGeometry->xNotch;
    double y0 = theGeometry->yNotch;
    double r0 = theGeometry->rNotch;
    
    
    double x1 = theGeometry->xHole;
    double y1 = theGeometry->yHole;
    double r1 = theGeometry->rHole;
 
//
//  -1- Construction de la g�om�trie avec OpenCascade
//      On cr�e le rectangle
//      On cr�e les deux cercles
//      On soustrait les cercles du rectangle :-)
//
 
    int ierr;
    int idPlate = gmshModelOccAddRectangle(-w/2, -h/2, 0, w, h, -1, 0,&ierr);   // c'est x, y et z lower left corner puis dx, dy, tag et rounded_radius
    ErrorGmsh(ierr);
    int idNotch = gmshModelOccAddDisk(x0, y0, 0, r0, r0, -1, NULL,0,NULL,0,&ierr);  // c'est xc, yc, zc, rx, ry, tag
    ErrorGmsh(ierr);
    int idHole  = gmshModelOccAddDisk(x1, y1, 0, r1, r1, -1, NULL,0,NULL,0,&ierr);    
    ErrorGmsh(ierr);
    
    int plate[] = {2, idPlate};
    int notch[] = {2, idNotch};
    int hole[]  = {2, idHole};
    gmshModelOccCut(plate, 2, notch, 2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);   // c'est objectDimTags, objectDimTags_n, toolDimtags, toolDimtags_n
    ErrorGmsh(ierr);
    gmshModelOccCut(plate, 2, hole, 2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);   // pareil
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