#include <stdio.h>
#include <math.h>
#include "glfem.h"


double integrate(double x[3], double y[3], double (*f) (double, double))
{
    double I = 0;
    double xLoc[3];
    double yLoc[3];
    const double xi[3] = {1.0/6, 1.0/6, 2.0/3};  // Pas oublier les .0 pour que ce soit des doubles !!!!
    const double eta[3] = {1.0/6, 2.0/3, 1.0/6};
    double jacobien = fabs((x[1] - x[0]) * (y[2] - y[0]) - ((x[2] - x[0]) * (y[1] - y[0])));

    for (int i = 0; i < 3; i++)
    {
        xLoc[i] = x[0] * (1 - xi[i] - eta[i]) + x[1] * xi[i] + x[2] * eta[i];
        yLoc[i] = y[0] * (1 - xi[i] - eta[i]) + y[1] * xi[i] + y[2] * eta[i];
        I += f(xLoc[i], yLoc[i]) * (1.0/6.0);  // fallait pas oublier de mettre 1.0/6.0 avec les .0 pour que ca soit un double
        // printf("xLoc[%d] = %f, yLoc[%d] = %f, f(xLoc[%d], yLoc[%d]) = %f\n", i, xLoc[i], i, yLoc[i], i, i, f(xLoc[i], yLoc[i]));
    }

    
    
    // Pour dessiner l'element, les sommets du triangle :-)
    // Decommenter la ligne pour dessiner aussi les points d'integration
    //

    glfemSetColor(GLFEM_BLACK); glfemDrawElement(x,y,3);
    glfemSetColor(GLFEM_BLUE);  glfemDrawNodes(x,y,3);
    glfemSetColor(GLFEM_RED);   glfemDrawNodes(xLoc,yLoc,3);


    return I * jacobien;
}

double integrateRecursive(double x[3], double y[3], double (*f)(double,double), int n)
{

//
// ... A modifier :-)
// y-compris la ligne juste en dessous :-)
//
    double I = 0;
    if (n <= 0) {
        I = integrate(x,y,f);
        return I;
    }
    else {
        // On crée 4 triangles dont l'union est le triangle initial
        // Les midpoints :
        I = 0;
        double M0[2] = {0.5*(x[0]+x[1]), 0.5*(y[0]+y[1])};
        double M1[2] = {0.5*(x[1]+x[2]), 0.5*(y[1]+y[2])};
        double M2[2] = {0.5*(x[2]+x[0]), 0.5*(y[2]+y[0])};
        // crée les 4 triangles :
        double x0[3] = {x[0], M0[0], M2[0]};
        double y0[3] = {y[0], M0[1], M2[1]};
        double x1[3] = {M0[0], x[1], M1[0]};
        double y1[3] = {M0[1], y[1], M1[1]};
        double x2[3] = {M2[0], M1[0], x[2]};
        double y2[3] = {M2[1], M1[1], y[2]};
        double x3[3] = {M0[0], M1[0], M2[0]};
        double y3[3] = {M0[1], M1[1], M2[1]};
        // On somme les intégrales des 4 triangles
        I = integrateRecursive(x0, y0, f, n-1) + integrateRecursive(x1, y1, f, n-1) + integrateRecursive(x2, y2, f, n-1) + integrateRecursive(x3, y3, f, n-1);
        
    }
    return I;
}
