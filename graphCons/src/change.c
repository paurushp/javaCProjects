/*
 *  change.c
 *  MutualInfo
 *
 *  Created by  on 2/17/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functionsMCMC.h"
#define CONVERGENCE 0.050
#define ACCEPT 0.25
#define SAMPLE 10000

int i, j, k, nsgenes;


void sample(int iter1, double** newNet, double*** lList2){
while(iter1 <= SAMPLE)
{
double** temp;
    **newNet = alterNet(temp);
    addToList(newNet, iter1, nsgenes, lList2);
    iter1++;
    temp = newNet;
}
}

void copyNet(int nsgenes, double** net, double** netCopy)
{


for (i = 0; i < nsgenes; i++){
for (j = 0; j < nsgenes; j++){
netCopy[i][j] = net[i][j];
}
}
}

void addToList(double **array, int pos, int size, double ***list)
{


for (i = pos; i <= size; i++){
        for(j = 0; j <= nsgenes; j++){
            for(k = 0; k <= nsgenes; k++){
list[i][j][k]= array[j][k];
}
}
}
}

int covergenceTest(int **arr1, int **arr2)
{
int converged;
double mInfo;
mInfo = mutInfo(arr1, arr2);
    if(mInfo <= CONVERGENCE)
converged = 1;
    else
converged = 0;    
    return(converged);


}
 */
#######################

/*
 *  alterNet.c
 *  MutualInfo
 *
 *  Created by  on 2/17/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
//#include <stdio.h>
//#include <stdlib.h>


#include "functionsMCMC.h"


double alterNet(double** net, int nsgenes, int T, double** newNet){
//int step = 0;
//Memory allocation
double** temp1 = (double**) calloc(nsgenes, sizeof(double*));
double** temp2 = (double**) calloc(nsgenes, double*);
int i, j, k;
for (i = 0; i < nsgenes; i++) {
temp1[i]= (double*) calloc(nsgenes, double);
temp2[i]= (double*) calloc(nsgenes, double);
}


copyNet(nsgenes, net, temp1);//copy the network for modification

 
}
*/ 