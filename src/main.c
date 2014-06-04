/*
 * main.c
 *
 *  Created on: Jun 1, 2014
 *      Author: mahoni
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "sparse.h"
#include "mmio/mmio.h"

#define USE_LAPACK_FORMAT

extern int mm_read(const char *fname, double **__val, int **__I, int **__J, int *_nz, int *_dim1, int *_dim2);

int main(int argc, char *argv[]){

	/* Part 1 - MPI initialization */

	int i,j,k,size,index;
	int index1,index2;
	int mynode, totalnodes;
	double alpha,gamma;
	MPI_Status status;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	size = (int) pow(2,log2(totalnodes+1)+1)-1;
	double **A;

    if (argc < 2) {
		printf("Usage: %s [martix-market-filename]\n", argv[0]);
		exit(1);
	}

	int numactivep = totalnodes;
	int * activep = (int *)malloc(totalnodes * sizeof(int));

#if 1
    int bw = 3;
	int kl, ku;
	kl = ku = (bw-2); /* assume symmetric matrix */

	sparsematrix_t smat;

	if(mynode == 0) {

		/* Read matrix */
		if (mm_read(argv[1], &smat.val, &smat.colidx, &smat.rowidx, &smat.nnz, &smat.ncol, &smat.nrow)) {
			printf("MM read error !!!\n");
			exit(1);
		}

		/* Memory allocation */
		A = (double **)malloc(smat.nrow * sizeof(double *));
		if(A == NULL) {
			fprintf(stderr, "out of memory\n");
			return -1;
		}

		for(i = 0; i < smat.nrow; i++) {
			A[i] = (double *)calloc(bw + 1, sizeof(double));
			A[i][bw] = 1.0;	/* XXX: right size is all one vector */
			if(A[i] == NULL) {
				fprintf(stderr, "out of memory\n");
				return -1;
			}
		}


		/* Convert to LAPACK Band format */
		for (i = 0; i < smat.nnz; i++) {
			A[smat.rowidx[i]][kl + smat.colidx[i] - smat.rowidx[i]] = smat.val[i];
		}

	}

	int numrows = smat.nrow;

#else
	const int numrows = 5;

	A = (double **)malloc(numrows * sizeof(double *));

	for(i=0;i<numrows;i++){
		A[i] = (double *)malloc((size+1) * sizeof(double));
		for(j=0;j<size+1;j++)
			A[i][j] = 0.0;
	}

	if(mynode==0) {
		A[0][0] = -2.0; A[0][1] = 1.0;
		A[1][0] = 1.0; A[1][1] = -2.0; A[1][2] = 1.0;
		A[2][1] = 1.0; A[2][2] = -2.0; A[2][3] = 1.0;
	} else if(mynode==(totalnodes-1)) {
		index = 2*mynode;
		A[0][index-1] = 1.0; A[0][index] = -2.0;
		A[0][index+1] = 1.0;
		index = 2*mynode+1;
		A[1][index-1] = 1.0; A[1][index] = -2.0;
		A[1][index+1] = 1.0;
		A[2][size-2] = 1.0; A[2][size-1] = -2.0;
	} else {
		for(i=0;i<3;i++){
			index = i + 2*mynode;
			A[i][index-1] = 1.0;
			A[i][index] = -2.0;
			A[i][index+1] = 1.0;
		}
	}

	for(i=0;i<3;i++)
		A[i][size] = 2*mynode+i;

	for(j=0;j<size+1;j++){
		A[3][j] = A[0][j];
		A[4][j] = A[2][j];
	}
#endif

	for(j=0;j<numactivep;j++)
		activep[j] = j;

	/*
	Remark 1: Just as in the parallel Gaussian elimination code, we augment the matrix A
	with the right-hand-side (appending A with an extra column). This helps to minimize the
	communication by allowing us to communicate both the row and right-hand-side information
	simultaneously.
	*/

	/* Part 2 - Cyclic reduction */
#ifdef USE_LAPACK_FORMAT
#	define CONV_COLIDX(row_idx, col_idx)	(kl+col_idx-row_idx)
#else
#	define CONV_COLIDX(row_idx, col_idx)	(col_idx)
#endif

	for(i=0;i<log2(size+1)-1;i++){
		for(j=0;j<numactivep;j++){
			if(mynode==activep[j]){
				index1 = 2*mynode + 1 - pow(2,i);
				index2 = 2*mynode + 1 + pow(2,i);
				alpha = A[1][CONV_COLIDX(1,index1)]/A[3][CONV_COLIDX(3,index1)];
				gamma = A[1][CONV_COLIDX(1,index2)]/A[4][CONV_COLIDX(4,index2)];

				for(k=0;k<size+1;k++)
					A[1][CONV_COLIDX(1,k)] -= (alpha*A[3][CONV_COLIDX(3,k)] + gamma*A[4][CONV_COLIDX(4,k)]);

				if(numactivep>1){
					if(j==0) {
						MPI_Send(A[1],size+1,MPI_DOUBLE,activep[1],0, MPI_COMM_WORLD);
					} else if(j==numactivep-1) {
						MPI_Send(A[1],size+1,MPI_DOUBLE,activep[numactivep-2], 1,MPI_COMM_WORLD);
					} else if(j%2==0) {
						MPI_Send(A[1],size+1,MPI_DOUBLE,activep[j-1], 1,MPI_COMM_WORLD);
						MPI_Send(A[1],size+1,MPI_DOUBLE,activep[j+1], 0,MPI_COMM_WORLD);
					} else {
						MPI_Recv(A[3],size+1,MPI_DOUBLE,activep[j-1],0, MPI_COMM_WORLD,&status);
						MPI_Recv(A[4],size+1,MPI_DOUBLE,activep[j+1],1, MPI_COMM_WORLD,&status);
					}
				}
			}
		}
		numactivep = 0;
		for(j=activep[1];j<totalnodes;j=j+pow(2,i+1)){
			activep[numactivep++]=j;
		}
	}

	/*
	Remark 2: The communication is accomplished through a series of MP I Send and MP I Recv
	calls. Each processor is communicating (either sending or receiving) from at most two other
	processors. To keep track of whom is to be sending/receiving, we maintain an active pro-
	cessor list within the integer array activep. A communication schematic for cyclic reduction
	using seven processors is given in figure 9.11.
	*/

	/* Part 3 - Back substitution */

	double * x = (double *)malloc(totalnodes);
	for(j=0;j<totalnodes;j++)
	x[j] = 0.0;
	if(mynode==activep[0]){
		x[mynode] = A[1][CONV_COLIDX(1,size)]/A[1][CONV_COLIDX(1,(size-1)/2)];
	}
	double tmp;
	for(i=log2(size+1)-3;i>=0;i--){
		tmp = x[mynode];
		MPI_Allgather(&tmp,1,MPI_DOUBLE,x,1,MPI_DOUBLE, MPI_COMM_WORLD);
		numactivep = 0;
		for(j=activep[0]-pow(2,i);j<totalnodes;j=j+pow(2,i+1)){
			activep[numactivep++]=j;
		}
		for(j=0;j<numactivep;j++){
			if(mynode == activep[j]){
				x[mynode] = A[1][CONV_COLIDX(1,size)];
				for(k=0;k<totalnodes;k++){
					if(k!=mynode)
						x[mynode] -= A[1][CONV_COLIDX(1,2*k+1)]*x[k];
				}
				x[mynode] = x[mynode]/A[1][CONV_COLIDX(1,2*mynode+1)];
			}
		}
	}
	tmp = x[mynode];
	MPI_Allgather(&tmp,1,MPI_DOUBLE,x,1,MPI_DOUBLE, MPI_COMM_WORLD);

	/*
	Remark 3: A schematic for the backward solve communication is given in figure 9.12. Notice
	that is varies slightly from that of the forward part of the reduction. After a processor has
	found the solution for its row and has communicated that information to the appropriate
	processors, it is no longer active. This can be observed in figure 9.12 - observe that processor
	P 3 no longer has things to compute in the second and third levels.
	/*
	Remark 4: We use the MPI_Allgather command so that at any given level all the processors have
	the available solution up to that point. This all inclusive communication could be
	replaced by MPI_Send/MPI_Recv pairs where only those processors requiring particular
	information would be updated.
	*/

	/* Part 4 - Solving for odd rows */
	for(k=0;k<totalnodes;k++){
		A[0][CONV_COLIDX(0,size)] -= A[0][CONV_COLIDX(0,2*k+1)]*x[k];
		A[2][CONV_COLIDX(2,size)] -= A[2][CONV_COLIDX(2,2*k+1)]*x[k];
	}
	A[0][CONV_COLIDX(0,size)] = A[0][CONV_COLIDX(0,size)]/A[0][CONV_COLIDX(0,2*mynode)];
	A[1][CONV_COLIDX(1,size)] = x[mynode];
	A[2][CONV_COLIDX(2,size)] = A[2][CONV_COLIDX(2,size)]/A[2][CONV_COLIDX(2,2*mynode+2)];
	free(activep);

	for(i=0;i<numrows;i++)
		free(A[i]);

	free(A);
	free(x);
	MPI_Finalize();

	return 0;
}

