  /*
  *  mi.c
  *
  *  Created by  on 2/6/11.
  *  Copyright 2011 __MyCompanyName__. All rights reserved.
  *
  */
  /*
  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include "functionsMCMC.h"
  */
  /*
  VARIABLES:

  
  v1 and v2 are the vector for which the MI has to be calculated.
  These vectors represent the adjacency matrix for the network reduced to
  a one dimension vector.

  
  ns is the size of vectors that has to be equal for both the vectors
  and MI is the mutual information for the two vectors. In principle ns will
  be the number of genes under consideration for the network.

  
  FUNCTION: 

  adapted from rWMBAT R-package

  

  
  */
/*
  int main()
  {


  int arr1[] = { 1 ,0, 1, 2, 1, 0, 0, 2};
  int arr2[] = { 2, 1, 1, 1, 0, 0, 0, 6};

  int size=sizeof(arr1)/sizeof(int);
  int size2 = sizeof(arr2)/sizeof(int);
  if (size != size2)
    {
    printf("Error: Vectors should be of same size\nThe program will exit now\n");
    exit(0);
    }
  else
    mutInfo(arr1, arr2, size);
  }
*/


  double mutInfo(int *v1, int *v2, int *ns, double *MI)  //function declaration to calculate mutual information
  {  
  int *Pa, *Pb, *Pab, n; // pointer declaration to probability of 'a', 'b' and their joint probability repectively
 
 // Memory allocation to the pointers

  Pa = (int *) calloc(3, sizeof(int));
  Pb = (int *) calloc(3, sizeof(int));
  Pab = (int *) calloc(9, sizeof(int));
  //Assigning size of the vector
  n = (int*)ns ; // cast warning***
  //printf("%d\n\n",n);
  int i,j,k;


  for (k=0; k<9; k++) 
    Pab[k] = 0;
  for (k = 0; k < 3; k++) {
  Pb[k] = 0;
  for (j = 0; j < n; j++)  
  if (v2[j] == k+1) 
  Pb[k]++;
  }
  for (k = 0; k < 3; k++) {
  Pa[k] = 0;
  for (j = 0; j < n; j++) {
  if (v1[j] == k+1) {
  Pa[k]++;
  for (i = 0; i < 3; i++) 
  if (v2[j] == i+1) 
  Pab[i*2 +k]++;
  }
  }
  }



  *MI = n * log(n);
  for (k = 0; k < 9; k++) 
  if (Pab[k]>0) 
  *MI = *MI + (Pab[k] * log(Pab[k])) ;


  for (k = 0; k < 3; k++) {
  if (Pb[k]>0) 
  *MI = *MI - (Pb[k] *log(Pb[k])) ;
  if (Pa[k]>0) 
  *MI = *MI - (Pa[k] *log(Pa[k])) ;
  }


  *MI = *MI / (double) n;
  *MI = *MI/ log(2);
  return (*MI);
  //printf ("%f\n\n",*MI);
  free (Pa);
  free (Pb);
  free (Pab);


  }

