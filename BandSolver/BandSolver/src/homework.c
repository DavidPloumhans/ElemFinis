

#include"fem.h"


#ifndef NORENUMBER 

void femMeshRenumber(femMesh *theMesh, femRenumType renumType)
{
    // seul truc louche c'est que j'ai la même pour X et Y au nv du résultat
    int i;
    int j;
    int n = theMesh->nodes->nNodes;
    double orderedX[n];
    int renumbered_alreadyX[n];
    double orderedY[n];
    int renumbered_alreadyY[n];
    // faudra être sûr qu'il soit content avec la déclaration de tout ça avant son "à modifier"
    switch (renumType) {
        case FEM_NO :
            for (i = 0; i < theMesh->nodes->nNodes; i++) 
                theMesh->nodes->number[i] = i;
            break;
// 
// A modifier :-)
// debut
//
        
        case FEM_XNUM:
            // va falloir numéroter les noeuds par ordre croissant de x
            // je copie l'array des x
            for(i = 0; i < n; i++) {
                orderedX[i] = theMesh->nodes->X[i];
            }
            // je print l'array avant le tri
            /*printf("Array avant le tri : ");
            for(i = 0; i < n; i++) {
                printf("%f ", orderedXX[i]);
            }
            */
            // je trie l'array
            for(i = 0; i < n - 1; i++) {
                for(j = i + 1; j < n; j++) {
                    if(orderedX[i] > orderedX[j]) {
                        double temp = orderedX[i];  // zebi t'avais pas mis un double mias un int et ça niquait tout
                        // printf("temp : %f\n", temp);
                        orderedX[i] = orderedX[j];
                        orderedX[j] = temp;
                    }
                }
            }
            // je print l'array après le tri
            /*
            printf("Array après le tri : ");
            for(i = 0; i < n; i++) {
                printf("%f ", orderedX[i]);
            }
            */
            // je renumérote les noeuds, attention aux doublons
            int k = 0;  // indice de renumbered_already
            for(i = 0; i < n; i++) {
                for(j = 0; j < n; j++) {
                    int pass = 1;
                    for (int p = 0; p <k; p++) {
                        if(j == renumbered_alreadyX[p]) {
                            pass = 0;  // le j-ème noeud a déjà été renuméroté avec une autre valeur de i
                            break;
                        }
                    }
                    if(theMesh->nodes->X[j] == orderedX[i] && pass == 1) {
                        theMesh->nodes->number[j] = i;
                        renumbered_alreadyX[k] = j;
                        k++;
                        break;
                    }
                }
            }
            /*
            for(i = 0; i < n; i++) {
                printf("Noeud %d : %d\n", i, theMesh->nodes->number[i]);
            }
            */
            break;  // fallait pas oublier le break !!!!!, sinon ça passait à Y

        case FEM_YNUM:
            // va falloir numéroter les noeuds par ordre croissant de y
            // même chose que pour FEM_XNUM
            for(i = 0; i < n; i++) {
                orderedY[i] = theMesh->nodes->Y[i];
            }
            for(i = 0; i < n - 1; i++) {
                for(j = i + 1; j < n; j++) {
                    if(orderedY[i] > orderedY[j]) {
                        double temp = orderedY[i];
                        orderedY[i] = orderedY[j];
                        orderedY[j] = temp;
                    }
                }
            }
            k = 0;
            for(i = 0; i < n; i++) {
                for(j = 0; j < n; j++) {
                    int pass = 1;
                    for (int p = 0; p <k; p++) {
                        if(j == renumbered_alreadyY[p]) {
                            pass = 0;
                            break;
                        }
                    }
                    if(theMesh->nodes->Y[j] == orderedY[i] && pass == 1) {
                        theMesh->nodes->number[j] = i;
                        renumbered_alreadyY[k] = j;
                        k++;
                        break;
                    }
                }
            }
            break;
        
// 
// end
//

        default : Error("Unexpected renumbering option"); }
    // Print les (X, Y) dans l'ordre de leur numéro
        /*
        for(i = 0; i < n; i++) {
            for(j = 0; j < n; j++) {
                if(theMesh->nodes->number[j] == i) {
                    printf("Noeud %d : (%f, %f)\n", i, theMesh->nodes->X[j], theMesh->nodes->Y[j]);
                    break;
                }
            }
        }
        */
}

#endif
#ifndef NOBAND 

int femMeshComputeBand(femMesh *theMesh)
{
    // je veux regarder les 3 noeuds (à priori 3) noeuds de chaque élément et en calculer la distance que ça ferait à la diagonale de la matrice
    int localNumberNodes = theMesh->nLocalNode;
    // en fait le tableau de int est un tableau 1D mais représentant un tableau de 2D. Il est de longueur nElem * nLocalNode et chaque valeur correspond à un noeud
    int maxDist = 0;
    for(int i = 0; i < theMesh->nElem; i++) {  // pour chaque élément
        int elemNodesNumber[localNumberNodes];  // j'y stocke les numéros des noeuds de l'élément
        for(int j = 0; j < localNumberNodes; j++) {  // pour chaque noeud de l'élément
            elemNodesNumber[j] = theMesh->nodes->number[theMesh->elem[i*localNumberNodes+j]];  // je récupère le numéro du j-ème noeud de l'élément i
        }
        // je veux calculer la plus grand distance entre chaque paire de noeuds
        for(int j = 0; j < localNumberNodes; j++) {
            for(int k = 0; k < localNumberNodes; k++) {
                int temp = abs(elemNodesNumber[j] - elemNodesNumber[k]);
                if(temp > maxDist) {
                    maxDist = temp;
                }
            }
        }
    }
    // printf("maxDist : %d\n", maxDist);  // à priori c'est bon
    return maxDist + 1;  // c'est dans la formule
}


#endif
#ifndef NOBANDASSEMBLE


void femBandSystemAssemble(femBandSystem* myBandSystem, double *Aloc, double *Bloc, int *map, int nLoc)
{
    //printf("%f ", myBandSystem->B[0]);
    // j'ai l'impression que size a déjà été fait et donc que la mémoire a à priori déjà été allouée
    // map va jusqu'à 100 donc on devrait avoir une matrice de taille 101
    // nloc vaut 3 (ou 4 si quadrilatères) donc on devrait avoir une matrice de taille 3x3 à rajouter partout le long de la diagonale
    // les Aloc sont bien de longueur 3*3 (en 1D ça donne 9) donc ça devrait être bon
    for (int i = 0; i < nLoc; i++) {
        myBandSystem->B[map[i]] += Bloc[i];  // met le += et pas = car on peut avoir plusieurs éléments locaux qui contribuent à un même élément global
    }
    for (int i = 0; i < nLoc; i++) {
        for (int j = 0; j < nLoc; j++) {
            if(map[i]<=map[j]) {  // map[i] <= map[j] pour exploiter la symétrie de la matrice
                myBandSystem->A[map[i]][map[j]] += Aloc[i*nLoc+j];  // même logique en n'oubliant pas que Aloc est un tableau 1D de taille 9 qui représente une matrice 3 par 3
            }
        }
    }
    // ça m'a l'air bon
}


#endif
#ifndef NOBANDELIMINATE


double *femBandSystemEliminate(femBandSystem *myBand)
{
    /*  printf("Matrix A:\n");  // juste pour voir la dégaine de la matrice A une fois complétée après plein d'appels à femBandSystemAssemble
    for (int i = 0; i < myBand->size; i++) {
        for (int j = 0; j < myBand->size; j++) {
            printf("%f ", myBand->A[i][j]);
        }
        printf("\n");
    }
    */
    // Version d'Eduardo
    double  **A, *B, factor;
    int     i, j, k, jend, size, band;
    A    = myBand->A;
    B    = myBand->B;
    size = myBand->size;
    band = myBand->band;

    for (int i = 0 ; i < size ; i++){
        if (A[i][i] != 0){
            jend = i+band < size ? i+band : size;
                for (int j = i+1 ; j < jend; j++){
                    //On exploite le fait que notre matrice soit symétrique
                    factor = A[i][j]/A[i][i];
                    for (int k = j ; k < jend ; k++){
                        A[j][k] -= factor*A[i][k];
                    }
                    B[j] = B[j] - factor*B[i];
            }
        }
    }
    // On doit ensuite effectuer une backward substitution
    for (i = size-1; i >= 0; i--) {
        factor = 0;
        jend = i+band < size ? i+band : size;
        for (j = i+1; j < jend; j++)
            factor += A[i][j]*B[j];
        B[i] = (B[i]-factor)/A[i][i];
    }
    return(myBand->B);
}


#endif

