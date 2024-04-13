#include "fem.h"


void geoMeshGenerate() {  // il ne faut pas y toucher

    femGeo* theGeometry = geoGetGeometry();

    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;
    
    int ierr;
    double r = w/4;
    int idRect = gmshModelOccAddRectangle(0.0,0.0,0.0,w,h,-1,0.0,&ierr); 
    int idDisk = gmshModelOccAddDisk(w/2.0,h/2.0,0.0,r,r,-1,NULL,0,NULL,0,&ierr); 
    int idSlit = gmshModelOccAddRectangle(w/2.0,h/2.0-r,0.0,w,2.0*r,-1,0.0,&ierr); 
    int rect[] = {2,idRect};
    int disk[] = {2,idDisk};
    int slit[] = {2,idSlit};

    gmshModelOccCut(rect,2,disk,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr); 
    gmshModelOccCut(rect,2,slit,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr); 
    gmshModelOccSynchronize(&ierr); 

    if (theGeometry->elementType == FEM_QUAD) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshOptionSetNumber("Mesh.RecombineAll",1,&ierr);
        gmshOptionSetNumber("Mesh.Algorithm",11,&ierr);  
        gmshOptionSetNumber("Mesh.SmoothRatio", 21.5, &ierr);  
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm",1.0,&ierr); 
        gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  
        gmshModelMeshGenerate(2,&ierr);  }
  
    if (theGeometry->elementType == FEM_TRIANGLE) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshModelMeshGenerate(2,&ierr);  }
    
 
    return;
}


double *femElasticitySolve(femProblem *theProblem)
{

    femFullSystem  *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete    *theSpace = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;  // les noeuds
    femMesh        *theMesh = theGeometry->theElements;  // les éléments
    
    
    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    int iElem,iInteg,iEdge,i,j,d,map[4],mapX[4],mapY[4], k, l;
    
    int nLocal = theMesh->nLocalNode;

    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;      
    double rho = theProblem->rho;
    double g   = theProblem->g;
    double **A = theSystem->A;
    double *B  = theSystem->B;
    
    
    //
    //  A faire :-)
    //                
    // Il faut juste construire le système linéaire à résoudre en ayant d'abod calculé toutes les intégrales sur chaque élément
    // Faudra prendre en compte le fait que ça puisse être un élément triangle ou quadrilatère (se sait avec nLocalNode)
    // Faudra aussi prendre en compte qu'on peut être en contrainte planes ou en déformations planes (voir theProblem->planarStrainStress)
    // A priori : a, b et c sont donnés donc faut pas prendre en compte le fait qu'on est en contrainte planes ou en déformations planes
    // même la règle d'intégration et tout est donnée

    // printf("Nombre de noeuds dans theNode : %d\n", theNodes->nNodes);
    // printf("Nombre de noeuds dans theMesh->nodes : %d\n", theMesh->nodes->nNodes);
    // les deux font référence au même pointeur
    // printf("Nombre d'éléments : %d\n", theMesh->nElem);

    // boucle sur les éléments
    for(iElem = 0; iElem < theMesh->nElem; iElem++) {
        // va falloir calculer les intégrales pour l'élément
        // donc pour chaque noeud puis assembler la matrice A avec la numération des noeuds, pareil pour le vecteur B
        // boucle sur les noeuds de l'élément, pour créer le mapping et récupérer les coordonnées x et y
        for(j = 0; j < nLocal; j++) {
            map[j] = theMesh->elem[iElem * nLocal + j];  // juste le numéro global du jème noeud local
            mapX[j] = 2 * map[j];  // Si tu écris la matrice 8x8, tu verras que t'as du xx en haut à gauche et yy en bas à droite donc cette numérotation est logique
            mapY[j] = 2 * map[j] + 1;  // faut pas oublier que A est de taille 2*nbNoeuds x 2*nbNoeuds
            // on récupère les coordonnées des noeuds de l'élément
            x[j] = theNodes->X[map[j]];  // elem contient les numéros des noeuds de l'élément mais de manière fourbe.
            y[j] = theNodes->Y[map[j]];  // le i-ème élément comprend les noeuds dont les numéros vont de elem[i*nLocal] à elem[i*nLocal+nLocal-1]
        }
        // boucle sur les noeuds de l'élément, qui sont aussi les points d'intégration
        for(iInteg = 0; iInteg < theRule->n; iInteg++) {
            // Toutes ces valeurs ne sont valables qu'au point iInteg  !!
            // Obtention des xsi, eta, des poids et des fonctions de forme ainsi que de leurs dérivées sur l'élément parent     
            theSpace->phi2(theRule->xsi[iInteg], theRule->eta[iInteg], phi);  // j'obtiens mes 3 ou 4 fonctions de forme
            theSpace->dphi2dx(theRule->xsi[iInteg], theRule->eta[iInteg], dphidxsi, dphideta);  // j'obtiens les dérivées de ces fonctions de forme
            
            // On calcule le jacobien
            // Tout ça vient de la formule du gradient de la transformation que j'ai faite sur papier
            double dxdxsi = 0.0;
            double dxdeta = 0.0;
            double dydxsi = 0.0;
            double dydeta = 0.0;
            for(i = 0; i < nLocal; i++) {
                dxdxsi += x[i] * dphidxsi[i];
                dxdeta += x[i] * dphideta[i];
                dydxsi += y[i] * dphidxsi[i];
                dydeta += y[i] * dphideta[i];
            }
            double jacobien = fabs(dxdxsi * dydeta - dxdeta * dydxsi);  // fabs nécessaire

            // Calculons A et B
            // On calcule d'abord les dérivées des fonctions de forme par rapport à x et y (4 dphidx et 4 dphidy)
            // d'abord récupéer dxsidx et compagnie par inversion des dxdxsi et compagnie
            // à nouveau, développements faits sur papier
            double dxsidx = 1.0 / jacobien * dydeta;
            double dxsidy = -1.0 / jacobien * dxdeta;
            double detadx = -1.0 / jacobien * dydxsi;
            double detady = 1.0 / jacobien * dxdxsi;
            for(i = 0; i < nLocal; i++) {
                dphidx[i] = dphidxsi[i] * dxsidx + dphideta[i] * detadx;
                dphidy[i] = dphidxsi[i] * dxsidy + dphideta[i] * detady;
            }
            // On peut maintenant calculer les intégrales
            // Calcule et place A
            for (i = 0; i < theSpace->n; i++) {
                for (j = 0; j < theSpace->n; j++) {
                    A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] +
                        dphidy[i] * c * dphidy[j]) * jacobien * theRule->weight[iInteg];
                    A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] +
                        dphidy[i] * c * dphidx[j]) * jacobien * theRule->weight[iInteg];
                    A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] +
                        dphidx[i] * c * dphidy[j]) * jacobien * theRule->weight[iInteg];
                    A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] +
                        dphidx[i] * c * dphidx[j]) * jacobien * theRule->weight[iInteg];
                }
            }

            // vecteur B
            for (i = 0; i < theSpace->n; i++) {
                B[mapY[i]] -= phi[i] * g * rho * jacobien * theRule->weight[iInteg];  // fonction de volume, les autres sont nulles
            }
        }
        
    }
    
  
    int *theConstrainedNodes = theProblem->constrainedNodes;     
    for (int i=0; i < theSystem->size; i++) {
        if (theConstrainedNodes[i] != -1) {
            double value = theProblem->conditions[theConstrainedNodes[i]]->value;
            femFullSystemConstrain(theSystem,i,value); }}
                            
    return femFullSystemEliminate(theSystem);
}
