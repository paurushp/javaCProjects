/*
 *  update.c
 *  MutualInfo
 *
 *  Created by  on 2/7/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functionsMCMC.h"
#include "netlearn.h"

//Function to get hastings update


double updateFactor (double** netOld, double** netNew, int nsgenes, int T, double*** D, double** egene_prior, int type, int nrep, double alpha, double beta, double hRatio){
double neighOld = 2* counterNonZero(netOld, nsgenes)+ counterNonZeroMax(netOld, nsgenes, T );  //Neighbors of the old network
double neighNew = 2* counterNonZero(netNew, nsgenes)+ counterNonZeroMax(netNew, nsgenes, T );  //Neighbors of the new network
double likLogOld = network_likelihood(netOld, nsgenes, nsgenes, T, D, egene_prior, type, nrep, alpha, beta);     //Likelihood for the old network
double likLogNew = network_likelihood(netNew, nsgenes, nsgenes, T, D, egene_prior, type, nrep, alpha, beta);     //Likelihood for the new network
double ratio = (likLogNew/likLogOld)+log(neighOld/neighNew);   // ****hastings ratio check formula
hRatio = minim(ratio);  // minimum of 1 and hastings ratio
return(hRatio);
}
//****merge to one function counter check the counter for reversible edge
double counterNonZero(double **adjMat, int nsgenes, int count1)// Gives reversible and decreased neighbors
{
int i, j; 
count1 = 0;
for (i=0; i<=nsgenes; i++)
    {
for (j = 0; j <=nsgenes; j++) 
{
        if(adjMat[i][j] > 0)
count1++;
}
    }
return(count1);
}

double counterNonZeroMax(double **adjMat, int nsgenes, int T, int count2)// Gives Increasible neighbors
{
int i, j; 
count2 = 0;
for (i=0; i<=nsgenes; i++)
    {
for (j = 0; j <=nsgenes; j++) 
{
        if((adjMat[i][j] > 0) & (adjMat[i][j] < T))
            count2++;
    }
}
return(count2);
}

//################################# log likelihood calculation (Plugged in from dynoNEM) ############# 

double** getPerturbProb(double** Psi, int T, int nsgenes, int k){
double** perturb_prob = (double**) calloc(nsgenes, sizeof(double*));
int parent_perturb_prob;
int s, t, p;
for(s = 0; s < nsgenes; s++){
perturb_prob[s] = (double*) calloc(T, sizeof(double));
}
for(t = 0; t < T; t++){
for(s = 0; s < nsgenes; s++){
perturb_prob[s][t] = 0; // default: s is unperturbed
perturb_prob[k][t] = 1; // the perturbed gene is always inactive
if(s != k){
for(p = 0; p < nsgenes; p++){
if(Psi[p][s] != 0 && abs(Psi[p][s]) <= t){ // p is a parent
if(t > 0)
parent_perturb_prob = perturb_prob[p][t-1];
else
parent_perturb_prob = (p == k);
if(parent_perturb_prob){
perturb_prob[s][t] = 1;
break;
}
}
}
}
}
}


return(perturb_prob);
}

double network_likelihood(double** Psi, int nsgenes, int negenes, int T, double*** D, double** egene_prior, int type, int nrep, double alpha, double beta){
double*** perturb_prob = (double***) calloc(nsgenes, sizeof(double**));
int s, k, t, i;
for(k = 0; k < nsgenes; k++)
perturb_prob[k] = (double**) getPerturbProb(Psi, T, nsgenes, k);


double tmp;
double loglik0;
double loglik=0;
double loglik_tmp;
for (i=0; i<negenes; i++) {
loglik_tmp=0;
for (s=0; s<nsgenes; s++) {
tmp=0;
for (k=0; k<nsgenes; k++) {
for (t=0; t<T; t++) {
if(egene_prior[i][s] > 0)
if(type == PVAL_DENS) // for p-value densities:
tmp += log((D[t][i][k]*perturb_prob[k][s][t] + (1-perturb_prob[k][s][t])*1) * egene_prior[i][s]);
else if(type == EFFECT_PROB) // for effect probabilities
tmp += log((D[t][i][k]*perturb_prob[k][s][t] + (1-perturb_prob[k][s][t])*(1 - D[t][i][k])) * egene_prior[i][s]);
else // for count data:
tmp += log((pow(1-beta, D[t][i][k]*perturb_prob[k][s][t])*pow(beta, (nrep-D[t][i][k])*perturb_prob[k][s][t]) +
pow(alpha, D[t][i][k]*(1-perturb_prob[k][s][t]))*pow(1-alpha, (nrep-D[t][i][k])*(1-perturb_prob[k][s][t]))) * egene_prior[i][s]);
}
}
if (s==0) {
loglik0 = tmp;
}
else {
loglik_tmp += exp(tmp - loglik0);
}


}
loglik += log(1 + loglik_tmp) + loglik0;
}
for(k = 0; k < nsgenes; k++){
for(s = 0; s < nsgenes; s++)
free(perturb_prob[k][s]);
free(perturb_prob[k]);
}
free(perturb_prob);
return(loglik);
}


