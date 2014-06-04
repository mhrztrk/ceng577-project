/*
 * sparse.c
 *
 *  Created on: Jun 2, 2014
 *      Author: mhrztrk
 */
#include "stdio.h"
#include "stdlib.h"
#include "sparse.h"


void spInit(sparsematrix_t *_smat, int nrow, int ncol, int nnz) {
}

void spDestroy(sparsematrix_t *_smat) {

}

void ConvertCSRFormat(sparsematrix_t *_smat, csr_t *_cmat) {

	_cmat->nvtxs = _smat->ncol;
	_cmat->xadj = (int *)malloc((_cmat->nvtxs + 1) * sizeof(int));
	_cmat->adjncy = (int *)malloc(_smat->nnz * sizeof(int));

	int i, idxptr;
	idxptr = 0;
	int k = 0;
	for (i = 0; i < _smat->nrow; i++) {
		_cmat->xadj[i] = k;
		while(_smat->rowidx[k] == i) {
			_cmat->adjncy[idxptr++] = _smat->colidx[k];
			k++;
			if(k >= _smat->nnz) {
				break;
			}
		}
	}
	_cmat->xadj[i] = k;

}


void ConvertLAPACKBandFormat(sparsematrix_t *_smat, double *Mat) {

	int i;

	for (i = 0; i < _smat->nnz; i++);
	/*
	int **array;
	array = malloc(nrows * sizeof(int *));
	if(array == NULL)
		{
		fprintf(stderr, "out of memory\n");
		exit or return
		}
	for(i = 0; i < nrows; i++)
		{
		array[i] = malloc(ncolumns * sizeof(int));
		if(array[i] == NULL)
			{
			fprintf(stderr, "out of memory\n");
			exit or return
			}
		}
	 */
	/*
	int i, j;
	for(i = 0; i < nrows; i++)
		{
		for(j = 0; j < ncolumns; j++)
			array[i][j] = 0;
		}
	*/

}


