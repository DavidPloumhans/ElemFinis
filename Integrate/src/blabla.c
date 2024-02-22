#include <stdio.h>
#include <math.h>
#include "glfem.h"


double integrate(double x[3], double y[3], double (*f) (double, double))
{
    double I = 0;
    double xLoc[3];
    double yLoc[3];
    double theta[3] = {0.166666666667, 0.166666666667, 0.666666666667}; 
    double etha[3] = {0.166666666667, 0.666666666667, 0.166666666667}; 
    double Jacobian = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]); 
for (int i = 0; i<3 ; i++){
    xLoc[i] = x[0]*(1-theta[i]-etha[i]) + x[1]*theta[i] + x[2]*etha[i];
    yLoc[i] = y[0]*(1-theta[i]-etha[i]) + y[1]*theta[i] + y[2]*etha[i];
}
 
for (int i = 0; i<3 ; i++){
    I += f(xLoc[i], yLoc[i])*(0.166666666667)*Jacobian;
}

//
//
// Pour dessiner l'element, les sommets du triangle :-)
// Decommenter la ligne pour dessiner aussi les points d'integration
//

  glfemSetColor(GLFEM_BLACK); glfemDrawElement(x,y,3);
  glfemSetColor(GLFEM_BLUE);  glfemDrawNodes(x,y,3);
  glfemSetColor(GLFEM_RED);   glfemDrawNodes(xLoc,yLoc,3);
  return I;
}

double integrateRecursive(double x[3], double y[3], double (*f)(double,double), int n)
{
    if (n<=0)return integrate(x, y, f);

    double abscisse [6] = {0.0, 1.0, 0.0, 0.5, 0.5, 0.0}; 
    double ordonnee [6] = {0.0, 0.0, 1.0, 0.0, 0.5, 0.5}; 

    int triangles [4][3] = {{0, 3, 5},{3, 4, 1},{5, 4, 2},{3, 4, 5}}; 

    double I = 0; 

    for (int i = 0; i <4; ++i){
        double xLoc[3]; 
        double yLoc[3]; 
        for (int j = 0; j<3; ++j){
            int points = triangles[i][j]; 
            double thetaLoc  = abscisse [points]; 
            double ethaLoc = ordonnee [points]; 
            xLoc[j] = x[0]*(1-thetaLoc-ethaLoc) + x[1]*thetaLoc + x[2]*ethaLoc;
            yLoc[j] = y[0]*(1-thetaLoc-ethaLoc) + y[1]*thetaLoc + y[2]*ethaLoc;
        }
        I += integrateRecursive(xLoc, yLoc,f, n-1); 
    }
}
