/*
 * matrix.h
 *
 *  Created on: Jun 4, 2014
 *      Author: mhrztrk
 */

#ifndef MATRIX_H_
#define MATRIX_H_

/* Band Storage:
 *
 * aij is stored in ab(ku+i-j,j) for max(1,j-ku)<= i <= min(n,j+kl)
 *
 * tridiagonal matrix (kl=ku=1):
 *
 * aij -> ab(1+i-j,j) for max(1,j-1) <= j <= min(n,j+1)
 *
 */


#endif /* MATRIX_H_ */
