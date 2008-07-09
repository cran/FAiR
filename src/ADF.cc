/*  This file is part of FAiR, a program to conduct Factor Analysis in R
    Copyright 2008 Benjamin King Goodrich

    FAiR is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version. 

    FAiR is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with FAiR.  If not, see <http://www.gnu.org/licenses/>.
*/
using namespace std;

extern "C" {

	void
	FAiR_QD_sum(double *sum, const double *x, const double *middle, 
		    const unsigned int *length) {
/*
		Arguments to Quadratic Distance discrepancy function:
			sum is the return value (initialized to zero)
			x is a vector equal to the difference between s and c
			middle is the vecs of the transpose of the Cholesky factor of W

		Note:   (s - c)'W(s - c) = x'[chol(W)' %*% chol(W)]x = z'z
*/
		int count = 0;
		double temp;
		for(unsigned int i = 0; i < *length; ++i) {
			temp = 0;
			for(unsigned int j = i; j < *length; ++j) {
/*
				matrix multiply x by a column of chol(W)'; since chol(W)'
				is lower-triangular, we do fewer operations each pass
*/
				temp += (x[j] * middle[count]);
				++count;
			}
			sum[0] += temp * temp;
		}
	}

	void
	FAiR_QD_grad(double *grad, const double *x, const double *middle, 
			const unsigned int *length, double *holder) {
/*
		Arguments to the gradient of the Quadratic distance discrepany function:
			grad is the return value (wrt to the *covariances*)
			x is a vector equal to the difference between s and c
			middle is the vecs of the transpose of the Cholesky factor of W
			length is the length of x
			holder is just a temporary storage vector (initialized to zero)

		Note:   the gradient of (s - c)'W(s - c) is -2W(s - c), which is equal to
				  -2[chol(W)' %*% chol(W)]x
*/

		int count = 0;
		int mark;
		double temp;
		for(unsigned int i = 0; i < *length; ++i) {
			temp = 0;
			for(unsigned int j = i; j < *length; ++j) {
/*
				matrix multiply a row of chol(W) by x; since chol(W) is
				upper-triangular, we do fewer operations each pass
*/
				temp += (middle[count] * x[j]);
				++count;
			}
			holder[i] = temp;

			temp = 0;
			mark = i;
			for(unsigned int j = 0; j <= i; ++j) {
/*
				matrix multiply a row of chol(W)' by holder; since
				chol(W)' is lower-triangular we only have to go part way
*/
				temp += (middle[mark] * holder[j]);
				mark += (*length - j - 1);
			}
			grad[i] = -2 * temp;
		}
	}
}
