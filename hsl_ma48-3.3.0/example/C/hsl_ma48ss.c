/* Simple example of use of HSL_MA48 */

#include <stdio.h>
#include <stdlib.h>
#include "hsl_ma48s.h"

int main(int argc, char **argv) {
   struct ma48_control control;
   struct ma48_ainfo ainfo;
   struct ma48_finfo finfo;
   struct ma48_sinfo sinfo;
   void *factors;

   float *val, *b, *x, res[2], err;
   int i, m, n, *row, *col;
   long ne;

   /* Read matrix order and number of entries */
   scanf("%d %d %ld\n", &m, &n, &ne);

   /* Allocate arrays of appropriate sizes */
   row = (int *) malloc(ne*sizeof(int));
   col = (int *) malloc(ne*sizeof(int));
   val = (float *) malloc(ne*sizeof(float));
   b = (float *) malloc(n*sizeof(float));
   x = (float *) malloc(n*sizeof(float));

   /* Read matrix and right-hand side */
   for(i=0; i<ne; i++) scanf("%d %d %f\n", &row[i], &col[i], &val[i]);
   for(i=0; i<n; i++) scanf("%f\n", &b[i]);

   /* Initialize the factors and control*/
   ma48_initialize(&factors);
   ma48_default_control(&control);

   /* Analyse and factorize */
   ma48_analyse(m,n,ne,row,col,val,factors,&control,&ainfo,&finfo,NULL,NULL);
   if(ainfo.flag != 0) {
      printf("Failure of ma48_analyse with ainfo.flag = %d\n", ainfo.flag);
      return 1;
   }

   /* Solve without iterative refinement */
   ma48_solve(m,n,ne,row,col,val,factors,b,x,&control,&sinfo,0,NULL,NULL);
   if(sinfo.flag == 0) {
      printf("Solution of first set of equations without refinement is:\n");
      for(i=0; i<n; i++) printf("%10.3lf ", x[i]);
      printf("\n\n");
   }

   /* read new matrix and right-hand side */
   for(i=0; i<ne; i++) scanf("%f\n", &val[i]);
   for(i=0; i<n; i++) scanf("%f\n", &b[i]);

   /* fast factorize */
   ma48_factorize(m,n,ne,row,col,val,factors,&control,&finfo,1,0);
   if(finfo.flag != 0) {
      printf("Failure of ma48_factorize with finfo.flag=%d\n", finfo.flag);
      return 1;
   }

   /* solve with iterative refinement */
   ma48_solve(m,n,ne,row,col,val,factors,b,x,&control,&sinfo,0,res,&err);
   if(sinfo.flag == 0) {
      printf("Solution of second system with refinement is:\n");
      for(i=0; i<n; i++) printf("%10.3f ", x[i]);
      printf("\nScaled residual is %10.3e %10.3e\n", res[0], res[1]);
      printf("Estimated error is %10.3e\n", err);
   }

   /* clean up */
   ma48_finalize(&factors, &control);
   free(row); free(col); free(val);
   free(b); free(x);

   return 0;
}
