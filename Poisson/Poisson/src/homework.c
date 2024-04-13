#include "fem.h"

# ifndef NOPOISSONCREATE

femPoissonProblem *femPoissonCreate(const char *filename)
{
    femGeo* theGeometry = geoMeshCreate(filename);
    femPoissonProblem *theProblem = malloc(sizeof(femPoissonProblem));
    theProblem->geo  = theGeometry;
    femMesh *theMesh = theGeometry->theElements;
    if (theMesh->nLocalNode == 4) {
        theProblem->space = femDiscreteCreate(4,FEM_QUAD);
        theProblem->rule = femIntegrationCreate(4,FEM_QUAD); }
    else if (theMesh->nLocalNode == 3) {
        theProblem->space = femDiscreteCreate(3,FEM_TRIANGLE);
        theProblem->rule = femIntegrationCreate(3,FEM_TRIANGLE); }
    theProblem->system = femFullSystemCreate(theMesh->nodes->nNodes);
    return theProblem;
}

# endif
# ifndef NOPOISSONBOUNDARY

void femPoissonFindBoundaryNodes(femPoissonProblem *theProblem)
{
    femGeo* theGeometry = theProblem->geo;  
    femMesh* theEdges = theGeometry->theEdges; 
    int nBoundary = theEdges->nElem;
    int boundary [nBoundary];
    for (int i = 0; i < nBoundary; i++) {
        int found_1 = 1;
        int found_2 = 1; 
        for (int j = 0; j < nBoundary; j++){
            if (theEdges->elem[2*i] == boundary[j]){
                found_1 = 0;
            } 
            if (theEdges->elem[2*i+1] == boundary[j]){
                found_2 = 0;
            }
            else if(j == nBoundary - 1){
                if(found_1 == 1){
                    boundary[i] = theEdges->elem[2*i];;
                }
                if (found_2 == 1){
                    boundary[i] = theEdges->elem[2*i+1];
                }
            }
        }
    }
    
    /*
    for (int i = 0; i < nBoundary; i++){
        printf("%d\n",boundary[i]);
    }
    */
    
    femDomain *theBoundary = malloc(sizeof(femDomain));
    theGeometry->nDomains++;
    theGeometry->theDomains = realloc(theGeometry->theDomains,theGeometry->nDomains*sizeof(femDomain*));
    theGeometry->theDomains[theGeometry->nDomains-1] = theBoundary;
    theBoundary->nElem = nBoundary;
    theBoundary->elem = malloc(nBoundary*sizeof(int));
    memcpy(theBoundary->elem,boundary,nBoundary*sizeof(int));
    theBoundary->mesh = NULL;
    sprintf(theBoundary->name,"Boundary");

    for (int sheng = 0; sheng < nBoundary;sheng++){
            printf("nBoundary = %d\n",theBoundary->elem[sheng]);
        }
}
    
# endif
# ifndef NOPOISSONFREE

void femPoissonFree(femPoissonProblem *theProblem)
{
    geoMeshFree(theProblem->geo);
    femDiscreteFree(theProblem->space);
    femIntegrationFree(theProblem->rule);
    femFullSystemFree(theProblem->system);
    free(theProblem);
}
    
# endif
# ifndef NOPOISSONLOCAL

void femPoissonLocal(femPoissonProblem *theProblem, const int iElem, int *map, double *x, double *y)
{
    femMesh *theMesh = theProblem->geo->theElements;
    int number_of_nodes = theMesh->nLocalNode;
    for (int i = 0; i < number_of_nodes; i++){
        // On stocke dans map les indices des noeuds du maillage
        map[i] = theMesh->elem[number_of_nodes*iElem+i];
        // On stocke dans x et y les coordonnées des noeuds du maillage
        x[i] = theMesh->nodes->X[map[i]];
        y[i] = theMesh->nodes->Y[map[i]];
    } 

}

# endif
# ifndef NOPOISSONSOLVE

void femPoissonSolve(femPoissonProblem *theProblem)
{

    femMesh *theMesh = theProblem->geo->theElements;
    femDomain *theBoundary = geoGetDomain(theProblem->geo,"Boundary");
    femFullSystem *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete *theSpace = theProblem->space;
 
    if (theSpace->n > 4) Error("Unexpected discrete space size !");  
    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    int iElem,iInteg,iEdge,i,j,map[4];
    int nLocal = theMesh->nLocalNode;
    // printf("Rentre dans 1er boucle\n");
    for (int iElem =0 ; iElem < theMesh->nElem ; iElem++){
        // printf("rentre dans 2e boucle\n");
        //femFullSystemPrint(theSystem);
        // On prend à chaque fois les données qui concerne un seul élément/maille 
        femPoissonLocal(theProblem,iElem,map,x,y); 
        for (int i = 0; i < theMesh->nLocalNode; i++){
            // printf("Rentre dans la 3e boucle\n");
            for (int j = 0; j < theMesh->nLocalNode; j++){
                double integral = 0.0;
                double integral_2 = 0.0;
                // printf("Rentre dans la 4 boucle\n");
                for (int h = 0; h <theRule->n; h++){
                    double xsi_loc = theRule->xsi[h];
                    double eta_loc = theRule->eta[h];
                    double weight_loc = theRule->weight[h];
                    theSpace->phi2(xsi_loc,eta_loc,phi);
                    theSpace->dphi2dx(xsi_loc,eta_loc,dphidxsi,dphideta);
                    double dxdxsi = 0.0;
                    double dxdeta = 0.0;
                    double dydxsi = 0.0;
                    double dydeta = 0.0;
                    // printf("Rentre dans la 5e boucle, %d, %d, %d\n", i, j, h);
                    for (int k = 0; k < nLocal; k++){
                        dxdxsi += dphidxsi[k]*x[k];
                        dxdeta += dphideta[k]*x[k];
                        dydxsi += dphidxsi[k]*y[k];
                        dydeta += dphideta[k]*y[k];
                    }
                    double jacobian = fabs(dxdxsi*dydeta - dxdeta*dydxsi);
                    double dphi_dx = (1/jacobian)*(-dydxsi*dphidxsi[i] + dydeta*dphidxsi[i]);
                    double dphi_dy = (1/jacobian)*(dxdeta*dphidxsi[i] - dxdxsi*dphidxsi[i]);
                    double dphi_dx2 = (1/jacobian)*(-dydxsi*dphidxsi[j] + dydeta*dphidxsi[j]);
                    double dphi_dy2 = (1/jacobian)*(dxdeta*dphidxsi[j] - dxdxsi*dphidxsi[j]);
                    if(i == j){
                        integral_2 += phi[i]*phi[j]*jacobian*weight_loc;
                    }
                    integral += (dphi_dx*dphi_dx2+dphi_dy*dphi_dy2)*jacobian*weight_loc;
                    // printf("%d, %d, %d\n", map[i], map[j], iElem);
                    // printf("sort dans l'integrale\n");
                    
                }
                theSystem->A[map[i]][map[j]] += integral;
                if (i == j){
                theSystem->B[map[i]] += integral_2;
                }
            }
        }
    }
    for (int sheng = 0; sheng < theBoundary->nElem;sheng++){
            femFullSystemConstrain(theSystem,theBoundary->elem[sheng],0.0);
        }
    femFullSystemEliminate(theSystem);
    return; 
}

# endif