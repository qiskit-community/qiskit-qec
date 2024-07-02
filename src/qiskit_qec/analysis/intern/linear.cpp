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

bool contradiction(const std::vector<std::vector<bool>> &a, const std::vector<bool> &b) {
  for (int row = 0; row < a.size(); ++row) {
    bool we_good = false;
    for (int col = 0; col < a[0].size(); ++col) {
      if (a[row][col]) {
        we_good = true;
        break;
      }
    }
    if (!we_good && b[row]) {
      return true;
    }
  }
  return false;
}

void debugToFile(const std::string& filename, const std::string& output) {
    std::ofstream debugFile; // Create an ofstream object for file operations
    debugFile.open(filename, std::ios::out | std::ios::app); // Open file in append mode
    
    // Check if file is open and ready for writing
    if (debugFile.is_open()) {
        debugFile << output << "\n";
        debugFile.close(); // Always close the file when you're done writing
    } else {
        // If the file couldn't be opened, you might want to handle the error
        //std::cerr << "Unable to open file for writing.\n";
    }
}

void debugToFile(const std::string& filename, const std::vector<int>& output) {
    std::ofstream debugFile; // Create an ofstream object for file operations
    debugFile.open(filename, std::ios::out | std::ios::app); // Open file in append mode
    
    // Check if file is open and ready for writing
    if (debugFile.is_open()) {
        debugFile << "[";
        // // Example: writing a vector's content to the file
        for (int value : output) {
            debugFile << "" << value << ",";
        }
        
        debugFile << "]";
        debugFile.close(); // Always close the file when you're done writing
    } else {
        // If the file couldn't be opened, you might want to handle the error
        //std::cerr << "Unable to open file for writing.\n";
    }
}

void debugToFile(const std::string& filename, int output) {
    std::ofstream debugFile; // Create an ofstream object for file operations
    debugFile.open(filename, std::ios::out | std::ios::app); // Open file in append mode
    
    // Check if file is open and ready for writing
    if (debugFile.is_open()) {
        debugFile << output << "\n";
        debugFile.close(); // Always close the file when you're done writing
    } else {
        // If the file couldn't be opened, you might want to handle the error
        //std::cerr << "Unable to open file for writing.\n";
    }
}

std::tuple<bool, std::vector<std::vector<bool>>, std::vector<bool>, std::vector<bool>> solve(std::vector<std::vector<bool>> &a, std::vector<bool> &b)
{
  if (contradiction(a, b)) { // 0=1 for at least one row, abort as no solution exists
    return std::make_tuple(false, a, b, b);
  } 
  int g = 0; // how many times gaussian elimination was actually done
  for (int col = 0; col < a[0].size(); ++col) { // go through all columns of A with increment variable i
    for (int row = g; row < a.size(); ++row) {
      if (a[row][col]) { // found first row with 1 in respective column
        // found one such row
        for (int target_row = 0; target_row < a.size(); ++target_row) { // Perform Gaussian elimination with the row found above targeting all other rows that have a 1 at column i
          if (target_row == row) {
            continue;
          }
          if (a[target_row][col]) {
            // actual gaussian elimination step between two rows
            for (int target_col = col; target_col < a[0].size(); ++target_col) {
              a[target_row][target_col] = a[target_row][target_col] ^ a[row][target_col];
            }
            b[target_row] = b[target_row] ^ b[row];
          }
        }

        // after do contradiciton step again
        if (contradiction(a, b)) { // 0=1 for at least one row, abort as no solution exists
          return std::make_tuple(false, a, b, b);
        } 

        // Swap the row that was used for elimination with the one at that has index at the current step 
        std::vector<bool> tmp_a = a[g];
        a[g] = a[row];
        a[row] = tmp_a;
        
        bool tmp_b = b[g];
        b[g] = b[row];
        b[row] = tmp_b;

        ++g; // increment g as gaussian elimination step is done

        break; //break out of search for row and go to next column
      }
    }
  }
  // now a contains the row reduced echolon form of the original a (with remaining all zeros rows)
  // and b the correspondingly transformed target vector
  // TODO: backsubsitution step to populate x still missing
  return std::make_tuple(true, a, b, b);
}

// what to implement this but for sparse matrix representation
// def ker2(a):
//     """ calculates the kernel of the binary matrix 'a' over the field GF(2). Adapted from code from S. Bravyi.
//     Returns a basis for the ker2(a) as rows of a 2d numpy.ndarray. """
//     m,n = a.shape
//     ker = np.identity(n,dtype=int)

//     for i in range(m):
//         y = np.dot(a[i], ker) % 2 # multiplication of current row with all columns of ker
//         good = ker[:,y==0] # columns of ker that are in the kernel of a[i,:] (and thus in the kernel of a[:i+1,:])
//         bad = ker[:, y==1] # colums of ker that are in kernel of a[:i,:] but not in kernel of a[i,:]
//         if bad.shape[1]>0: # in case there are enough columns not in the kernel
//             new_good = (bad[:,:-1] + bad[:,1:]) % 2 # by construction all of these will be in kernel of a[i,:], independent and complete
//             ker = np.concatenate((good, new_good), axis=1) # new basis for kernel of a[:i+1,:]
//     # now columns of ker span the binary null-space of a
//     return np.transpose(ker)