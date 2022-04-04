#include "linear.h"


/*
 * Evaluate the symplectic inner product of two length-2n
 * binary vectors, i.e., as Pauli operators, do they commute?
 *
 * Assume the input is valid: each vector has the same even
 * length and elements that are zero or one.
 */
int symplectic_inner_product(std::vector<int> symplectic_a,
                             std::vector<int> symplectic_b)
{
  int n = symplectic_a.size() / 2;
  int product = 0;
  for(int i=0; i<n; i++) {
    product += symplectic_a[i] * symplectic_b[i+n];
    product += symplectic_a[i+n] * symplectic_b[i];
  }
  return product % 2;
}

/*
 * Decide if a list of length-2n binary vectors span an
 * isotropic subspace W \subseteq W^\perp, i.e., as Pauli operators,
 * do they commute?
 */
bool is_isotropic(std::vector<std::vector<int> > &symplectic_vectors)
{
  if(symplectic_vectors.size() < 2) return true;
  for(int i=0; i<symplectic_vectors.size(); i++)
    for(int j=i+1; j<symplectic_vectors.size(); j++)
      if(symplectic_inner_product(symplectic_vectors[i],
                                  symplectic_vectors[j]) == 1)
        return false;
  return true;
}

/*
 * Given a list of length-2n binary vectors that span an isotropic
 * subspace W, decide if another length-2n binary vector x belongs
 * to W^\perp, i.e., as a Pauli operator, does x commute with the
 * stabilizer group?
 * Assume the input is valid.
 */
bool is_orthogonal(std::vector<int> x,
                   std::vector<std::vector<int> > &symplectic_vectors)
{
  for(int i=0; i<symplectic_vectors.size(); i++)
    if(symplectic_inner_product(symplectic_vectors[i], x) == 1)
      return false;
  return true;
}

/*
 * Compute the dimension of the space spanned by a set of binary vectors.
 * Assume each vector has the same length and elements that are zero or one.
 */
int rank(std::vector<std::vector<int> > vectors)
{
  int rows = vectors.size();
  int cols = vectors[0].size();
  int i = 0;
  for(int k=0; k<cols; k++) {
    int j = 0;
    for(j=i; j<rows; j++)  // find non-zero element in this column
      if(vectors[j][k] == 1) break;
    if(j < rows) {  // if a non-zero element exists
      // swap rows i and j
      std::vector<int> temp = vectors[i];
      vectors[i] = vectors[j];
      vectors[j] = temp;
      // zero the rest of the column
      for(int l=i+1; l<rows; l++)
        if(vectors[l][k] == 1) {
          for(int m=0; m<cols; m++)  // add row i to row l
            vectors[l][m] = (vectors[l][m] + vectors[i][m]) % 2;
        }
      i++;
    }    
  }
  return i;
}

