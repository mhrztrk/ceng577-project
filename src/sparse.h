/*
 * sparse.h
 *
 *  Created on: Jun 2, 2014
 *      Author: mhrztrk
 */

#ifndef SPARSE_H_
#define SPARSE_H_

typedef struct sparsematrix_s {
	double *val;	/* nonzero entries */
	int *rowidx;	/* nonzero row indices */
	int *colidx;	/* nonzero column indices */
	int nnz;	/* # of nonzero entries */
	int nrow;	/* row length */
	int ncol;	/* column length */
}sparsematrix_t;


/* Compressed Storage Format */
typedef struct csr_s {
	int nvtxs;	/* # of vertices in the graph */
	int *xadj;
	int *adjncy;
} csr_t;

void ConvertCSRFormat(sparsematrix_t *_smat, csr_t *_cmat);

#endif /* SPARSE_H_ */
