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

    // Obtention des xsi, eta, des poids et des fonctions de forme ainsi que de leurs dérivées sur l'élément parent
    for (k = 0; k < nLocal; k++) {  // itère sur les points de la règle d'intégration
        theSpace->phi2(theRule->xsi[k], theRule->eta[k], phi);  // j'obtiens mes 3 ou 4 fonctions de forme
        theSpace->dphi2dx(theRule->xsi[k], theRule->eta[k], dphidxsi, dphideta);  // j'obtiens les dérivées de ces fonctions de forme
    }

    // boucle sur les éléments
    for(i = 0; i < theMesh->nElem; i++) {
        // va falloir calculer les intégrales pour l'élément
        // donc pour chaque noeud puis assembler la matrice A avec la numération des noeuds, pareil pour le vecteur B
        // boucle sur les noeuds de l'élément
        for(j = 0; j < nLocal; j++) {
            // on récupère les coordonnées du noeud
            x[j] = theNodes->X[theMesh->elem[i*nLocal+j]];  // elem contient les numéros des noeuds de l'élément mais de manière fourbe.
            y[j] = theNodes->Y[theMesh->elem[i*nLocal+j]];  // le i-ème élément comprend les noeuds dont les numéros vont de elem[i*nLocal] à elem[i*nLocal+nLocal-1]
            // On va devoir calculer nos intégrales
            // Il nous faut de quoi passer sur l'élément parent pour construire la matrice A de raideur et le vecteur B pour chaque noeud
            // Ensuite, il faudra les placer adéquatement dans la grande matrice A et le grand vecteur B
            // Il va y avoir 6 éléments à calculer (4 dans A et 2 dans B) avec 12 intégrales

            double A11;
            double A12;
            double A21;
            double A22;
            double B1;
            double B2;
            // Pour faire l'intégrale, on va calculer le jacobien pour chaque point de la règle d'intégration
            // itération sur les points de la règle d'intégration
            for(k = 0; k < nLocal; k++) {
                // On calcule le jacobien
                // Tout ça vient de la formule du gradient de la transformation que j'ai faite sur papier
                double dxdxsi = 0;
                double dxdeta = 0;
                double dydxsi = 0;
                double dydeta = 0;
                for(l = 0; l < nLocal; l++) {
                    dxdxsi += x[l] * dphidxsi[l];
                    dxdeta += x[l] * dphideta[l];
                    dydxsi += y[l] * dphidxsi[l];
                    dydeta += y[l] * dphideta[l];
                }
                double jacobien = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
                // Calculons A et B
                // On calcule d'abord les dérivées des fonctions de forme par rapport à x et y (dphidx et dphidy)
                // d'abord récupéer detadx et compagnie
                // à nouveau, développement faits sur papier
                double dxsidx = 1.0 / jacobien * dydeta;
                double dxsidy = -1.0 / jacobien * dxdeta;
                double detadx = -1.0 / jacobien * dydxsi;
                double detady = 1.0 / jacobien * dxdxsi;
                for(l = 0; l < nLocal; l++) {
                    dphidx[l] = dphidxsi[l] * detadx + dphideta[l] * detadx;
                    dphidy[l] = dphidxsi[l] * dxsidy + dphideta[l] * detady;
                }
                // On peut maintenant calculer les intégrales des produits des dphidx et dphidy
                double I_dphidx_dphidx = 0;
                double I_dphidy_dphidy = 0;
                double I_dphidx_dphidy = 0;  // je pense que c'est la même que I_dphidy_dphidx
                for(k = 0; k < nLocal; k++) {
                    I_dphidx_dphidx += theRule->weight[k] * dphidx[k] * dphidx[k];
                    I_dphidy_dphidy += theRule->weight[k] * dphidy[k] * dphidy[k];
                    I_dphidx_dphidy += theRule->weight[k] * dphidx[k] * dphidy[k];  // inverser les deux derniers facteurs ne devrait rien changer
                }
                // On peut maintenant calculer les éléments de la matrice A
                A11 = a * I_dphidx_dphidx + c * I_dphidy_dphidy;
                A12 = b * I_dphidx_dphidy + c * I_dphidx_dphidy;
                A21 = b * I_dphidx_dphidy + c * I_dphidx_dphidy;
                A22 = a * I_dphidy_dphidy + c * I_dphidx_dphidx;

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
