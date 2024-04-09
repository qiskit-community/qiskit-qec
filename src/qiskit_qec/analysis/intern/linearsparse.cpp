#include "linearsparse.h"

/* debuggin stuff */
void debugToFile_B(const std::string& filename, const std::string& output) {
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

void debugToFile_B(const std::string& filename, const std::vector<int>& output) {
    std::ofstream debugFile; // Create an ofstream object for file operations
    debugFile.open(filename, std::ios::out | std::ios::app); // Open file in append mode
    
    // Check if file is open and ready for writing
    if (debugFile.is_open()) {
        debugFile << "[";
        // // Example: writing a vector's content to the file
        for (int value : output) {
            debugFile << "" << value << ",";
        }
        
        debugFile << "]\n";
        debugFile.close(); // Always close the file when you're done writing
    } else {
        // If the file couldn't be opened, you might want to handle the error
        //std::cerr << "Unable to open file for writing.\n";
    }
}

void debugToFile_B(const std::string& filename, int output) {
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

/* conversion operations */
std::vector<int> dense2sparse(const std::vector<bool> a) {
  //debugToFile_B("debug_backst.txt", "in make sparse");
  std::vector<int> res;
  for (int i = 0; i < a.size(); i++) {
    if (a[i]) {
      res.push_back(i);
    }
  }
  return res;
}

std::vector<std::vector<int>> dense2sparse(const std::vector<std::vector<bool>> a) {
  std::vector<std::vector<int>> res;
  for (int i = 0; i < a.size(); i++) {
    res.push_back(dense2sparse(a[i]));
  }
  return res;
}

/* logical operations */

std::vector<int> lor_sparse(const std::vector<int> &a1, const std::vector<int> &a2) {
  // assumes both inputs are sorted. computes set intersection a1 U a2
  std::vector<int> res;
  int idx1 = 0;
  int idx2 = 0;
  while(idx1 < a1.size() && idx2 < a2.size()) {
    if (a1[idx1] < a2[idx2]) {
      res.push_back(a1[idx1++]);
    }
    else if (a1[idx1] > a2[idx2]) {
      res.push_back(a2[idx2++]);
    }
    else { // a1[idx1] == a2[idx2]
      res.push_back(a1[idx1++]);
      idx2++;
    }
  }
  return res;
}

std::vector<int> land_sparse(const std::vector<int> &a1, const std::vector<int> &a2) {
  // assumes both inputs are sorted. computes set intersection a1 U a2
  std::vector<int> res;
  int idx1 = 0;
  int idx2 = 0;
  while(idx1 < a1.size() && idx2 < a2.size()) {
    if (a1[idx1] < a2[idx2]) {
      idx1++;
    }
    else if (a1[idx1] > a2[idx2]) {
      idx2++;
    }
    else { // a1[idx1] == a2[idx2]
      res.push_back(a1[idx1++]);
      idx2++;
    }
  }
  return res;
}

std::vector<int> lxor_sparse(const std::vector<int> &a1, const std::vector<int> &a2) {
  // assumes both inputs are sorted. computes set difference {x| x \in a1 and x \notin a2 or x \notin a1 and x \in a2}
  std::vector<int> res;
  int idx1 = 0;
  int idx2 = 0;
  while(idx1 < a1.size() && idx2 < a2.size()) {
    if (a1[idx1] < a2[idx2]) {
      res.push_back(a1[idx1++]);
    }
    else if (a1[idx1] > a2[idx2]) {
      res.push_back(a2[idx2++]);
    }
    else { // a1[idx1] == a2[idx2]
      idx1++;
      idx2++;
    }
  }
  for (; idx1 < a1.size(); idx1++) {
      res.push_back(a1[idx1]);
  }
  for (; idx2 < a2.size(); idx2++) {
    res.push_back(a2[idx2]);
  }
  return res;
}

std::vector<int> lxor_sparse_a(std::vector<std::vector<int>> &a) {
  /* Computes this for an arbitray number of arrays
  */
  std::vector<int> res;
  std::vector<int> idxs = std::vector<int>(a.size());

  while (idxs.size() > 0) {
    // find minimum value
    int min_val = a[0][idxs[0]];
    for (int i = 1; i < idxs.size(); i++) {
      if (a[i][idxs[i]] < min_val) {
        min_val = a[i][idxs[i]];
      }
    }
    // count parity of occurres and increment counters (could also just keep track and iterate once, but this is less overhead I feel)
    bool parity = 0;
    for (int i = 0; i < idxs.size(); i++) {
      if (a[i][idxs[i]] == min_val) {
        idxs[i]++;
        parity ^= 1;
      }
    }
    // use parity to change result
    if (parity == 1) {
      res.push_back(min_val);
    }
    // remove rows where we reached the end
    for (int i = idxs.size()-1; i >=0; i--) {
      if (idxs[i] == a[i].size()) {
        idxs.erase(idxs.begin() + i);
        a.erase(a.begin() + i);
      }
    }
  }
  return res;
}

int lxor_weight_sparse(const std::vector<int> &a1, const std::vector<int> &a2) {
  // assumes both inputs are sorted. computes |{x| x \in a1 and x \notin a2 or x \notin a1 and x \in a2}|
  // directly without creating the set, for faster computation
  int res = 0;
  int idx1 = 0;
  int idx2 = 0;
  while(idx1 < a1.size() && idx2 < a2.size()) {
    if (a1[idx1] < a2[idx2]) {
      res++;
      idx1++;
    }
    else if (a1[idx1] > a2[idx2]) {
      res++;
      idx2++;
    }
    else { // a1[idx1] == a2[idx2]
      idx1++;
      idx2++;
    }
  }
  if (idx1 < a1.size()) {
    res += (a1.size() - idx1 - 1);
  }
  if (idx2 < a2.size()) {
    res += (a2.size() - idx2 - 1);
  }
  return res;
}

int lxor_weight_sparse_a(std::vector<std::vector<int>> &a) {
  /* Computes the number of integers that are in an odd number of arrays(std::vector<int>), for an arbitrary 
  number of such arrays contained within a.
  assumes the arrays to be sorted.
  */
  int res = 0;
  std::vector<int> idxs = std::vector<int>(a.size());

  while (idxs.size() > 0) {
    // find minimum value
    int min_val = a[0][idxs[0]];
    for (int i = 1; i < idxs.size(); i++) {
      if (a[i][idxs[i]] < min_val) {
        min_val = a[i][idxs[i]];
      }
    }
    // count parity of occurres and increment counters (could also just keep track and iterate once, but this is less overhead I feel)
    bool parity = 0;
    for (int i = 0; i < idxs.size(); i++) {
      if (a[i][idxs[i]] == min_val) {
        idxs[i]++;
        parity ^= 1;
      }
    }
    // add parity to result
    res += parity;
    // remove rows where we reached the end
    for (int i = idxs.size()-1; i >=0; i--) {
      if (idxs[i] == a[i].size()) {
        idxs.erase(idxs.begin() + i);
        a.erase(a.begin() + i);
      }
    }
  }
  return res;
}

std::vector<int> ldif_sparse(const std::vector<int> &a1, const std::vector<int> &a2) {
  // assumes both inputs are sorted. computes set difference a1\a2
  std::vector<int> res;
  int idx1 = 0;
  int idx2 = 0;
  while(idx1 < a1.size() && idx2 < a2.size()) {
    if (a1[idx1] < a2[idx2]) {
      res.push_back(a1[idx1++]);
    }
    else if (a1[idx1] > a2[idx2]) {
      idx2++;
    }
    else { // a1[idx1] == a2[idx2]
      idx1++;
      idx2++;
    }
  }
  for (; idx1 < a1.size(); idx1++) {
      res.push_back(a1[idx1]);
  }
  return res;

}

bool is_orthogonal_sparse(const std::vector<int> &a1, const std::vector<int> &a2) {
  return land_sparse(a1, a2).size() % 2 == 0;
}

bool in_support_sparse(const std::vector<int> &a, int x) {
  for (int i = 0; i < a.size(); ++i) {
    if (a[i] == x) {
      return true;
    }
  }
  return false;
}

/* linalg operations */

std::vector<std::vector<int>> nullspace_sparse(const std::vector<std::vector<int>> &a, int n) {
  // computes the kernel/nullspace of the matrix a over GF(2)
  
  // initialize kernel
  std::vector<std::vector<int>> ker;
  for (int i = 0; i < n; i++) {
    std::vector<int> row;
    row.push_back(i);
    ker.push_back(row);
  }

  // compute kernel
  for (int i = 0; i < a.size(); i++) {
    std::vector<std::vector<int>> new_ker;
    std::vector<std::vector<int>> probs;
    for (int j = 0; j < ker.size(); j++) {
      if (is_orthogonal_sparse(a[i], ker[j])) {
        new_ker.push_back(ker[j]);
      } else {
        probs.push_back(ker[j]);
      }
    }
    for (int j = 0; j+1 < probs.size(); j++) {
      new_ker.push_back(lxor_sparse(probs[j], probs[j+1]));
    }
    ker = new_ker;
  }
  return ker;
}

int nullity_sparse(const std::vector<std::vector<int>> &a, int n) {
  /* computes the nullity (dimension of nullspace) of a*/
  return nullspace_sparse(a, n).size();
}

int rank_sparse(const std::vector<std::vector<int>> &a, int n) {
  /* computes the rank (dimension of rowspan) of a*/
  return n - nullity_sparse(a, n);
}

bool lse_invalid_sparse(const std::vector<std::vector<int>> &a, const std::vector<bool> &b) {
  for (int row = 0; row < a.size(); ++row) {
    if (a[row].size()==0 && b[row]) {
      return true;
    }
  }
  return false;
}

std::tuple<bool, std::vector<std::vector<int>>, std::vector<bool>> gaussian_elimination_sparse(std::vector<std::vector<int>> &a, int n, std::vector<bool> &b)
{ // a is in sparse representation, b is not
  if (lse_invalid_sparse(a, b)) { // 0=1 for at least one row, abort as no solution exists
    return std::make_tuple(false, a, b);
  } 
  int g = 0; // how many times gaussian elimination was actually done
  for (int col = 0; col < n; ++col) { // go through all columns of A with increment variable i
    for (int row = g; row < a.size(); ++row) {
      if (in_support_sparse(a[row], col)) { // found first row with above g that contains a 1 at column
        for (int target_row = 0; target_row < a.size(); ++target_row) { // Perform Gaussian elimination with the row found above targeting all other rows that have a 1 at column i
          if (target_row == row) {
            continue;
          }
          if (in_support_sparse(a[target_row], col)) {
            // actual gaussian elimination step between two rows
            a[target_row] = lxor_sparse(a[row], a[target_row]);
            b[target_row] = b[target_row] ^ b[row];
          }
        }
        // after do contradiciton step again
        if (lse_invalid_sparse(a, b)) { // 0=1 for at least one row, abort as no solution exists
          return std::make_tuple(false, a, b);
        } 

        // Swap the row that was used for elimination with the one at that has index at the current step 
        std::vector<int> tmp_a = a[g];
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
  return std::make_tuple(true, a, b);
}

std::vector<int> back_substitution_sparse(std::vector<std::vector<int>> &a, int n, std::vector<bool> &b) {
  // assumes that a, b are correct output of gaussian elimination
  //debugToFile_B("debug_backst.txt", "at beginning, m,r=");
  int m = a.size();
  //debugToFile_B("debug_backst.txt", m);
  int r = std::min(m, n);
  //debugToFile_B("debug_backst.txt", r);
  std::vector<bool> x = std::vector<bool>(n);
  for (int i = 0; i < r; i++) {
    std::vector<int> row = a[r-1-i];
    //debugToFile_B("debug_backst.txt", "row iteration (i, row)=");
    //debugToFile_B("debug_backst.txt", i);
    //debugToFile_B("debug_backst.txt", row);
    //debugToFile_B("debug_backst.txt", "b[r-1-i]=");
    //debugToFile_B("debug_backst.txt", b[r-1-i]);
    if (b[r-1-i]) {
      //debugToFile_B("debug_backst.txt", "in if statement");
      //debugToFile_B("debug_backst.txt", "row.size()-1=");
      //debugToFile_B("debug_backst.txt", row.size()-1);
      //debugToFile_B("debug_backst.txt", "row[row.size()-1]=");
      //debugToFile_B("debug_backst.txt", row[row.size()-1]);
      x[row[row.size()-1]] = true; // this always exists, otherwise no solution, but we assume solution
      //debugToFile_B("debug_backst.txt", "x set, now changing b vector");
      for (int j = 0; j < r - 1 - i; j++) {
        // now change b vector
        //debugToFile_B("debug_backst.txt", "j,a[j]");
        //debugToFile_B("debug_backst.txt", j);
        //debugToFile_B("debug_backst.txt", a[j]);
        //debugToFile_B("debug_backst.txt", "in_support_sparse(a[j], row[row.size()-1])");
        //debugToFile_B("debug_backst.txt", in_support_sparse(a[j], row[row.size()-1]));
        if (in_support_sparse(a[j], row[row.size()-1])) {
          b[j] = !b[j];
        }
      }
    }
    //debugToFile_B("debug_backst.txt", "// set looked up columns to 0 in all rows above");
    // set looked up columns to 0 in all rows above
    for (int j = 0; j < r - 1 - i; j++) {
      //debugToFile_B("debug_backst.txt", "j,a[j]");
      //debugToFile_B("debug_backst.txt", j);
      //debugToFile_B("debug_backst.txt", a[j]);
      a[j] = ldif_sparse(a[j], row);
      //debugToFile_B("debug_backst.txt", a[j]);
    }
  }
  //debugToFile_B("debug_backst.txt", "got x=");
  for (bool eee : x) {
    //debugToFile_B("debug_backst.txt", int(eee));
  }
  //debugToFile_B("debug_backst.txt", "tryig to make sparse");
  //debugToFile_B("debug_backst.txt", dense2sparse(x));
  return dense2sparse(x);
}

std::tuple<bool, std::vector<int>> solve_sparse(std::vector<std::vector<int>> &a, int n, std::vector<bool> &b)
{ 
  /* a is in sparse representation, b is not. Solves the system of equations Ax=b if it is solvable. */

  // perform gaussian elimimation
  std::tuple<bool, std::vector<std::vector<int>>, std::vector<bool>> res = gaussian_elimination_sparse(a, n, b);

  bool solvable = std::get<0>(res);
  std::vector<std::vector<int>> a_rref = std::get<1>(res);
  std::vector<bool> b_p = std::get<2>(res);

  // perform back substitution
  std::vector<int> x;
  if (solvable) {
    x = back_substitution_sparse(a_rref, n, b_p);
  }
  return std::make_tuple(solvable, x);
}

std::tuple<std::vector<int>, int, bool> minimize_weight_sparse(const std::vector<int> &x_part, const std::vector<std::vector<int>> &nullspace, int max_dof) {
  if (nullspace.size() > max_dof) {
    return std::make_tuple(x_part, nullspace.size(), false);
  }
  std::vector<int> x;
  std::vector<int> x_star = x_part;
  int min_weight = x_part.size();
  for (int i = 1; i < 1 << nullspace.size(); i++) {
    std::vector<std::vector<int>> supported = std::vector<std::vector<int>>();
    supported.push_back(x_part);
    for (int j = 0; j < nullspace.size(); j++) {
      if ((i>>j) & 1) {
        supported.push_back(nullspace[j]);
      }
    }
    int weight = lxor_weight_sparse_a(supported);
    if (weight < min_weight) {
      min_weight = weight;
      x_star = lxor_sparse_a(supported);
    }
  }
  return std::make_tuple(x_star, nullspace.size(), true);
}

std::tuple<bool, std::vector<int>, int, bool> solve_optimal_sparse(std::vector<std::vector<int>> &a, int n, std::vector<bool> &b, int max_dof)
{ 
  /* a is in sparse representation, b is not. Solves the system of equations Ax=b if it is solvable. */

  // perform gaussian elimimation
  std::tuple<bool, std::vector<std::vector<int>>, std::vector<bool>> res = gaussian_elimination_sparse(a, n, b);

  bool solvable = std::get<0>(res);
  std::vector<std::vector<int>> a_rref = std::get<1>(res);
  std::vector<bool> b_p = std::get<2>(res);

  // perform back substitution
  std::vector<int> x;
  if (solvable) {
    // we can use original a or a_rref (which is actually the same, as gaussian elimination operates in place)
    // but it only uses adding and exchanging rows, so nullspace is conserved
    // we have to do it before back_substitution, because that also operates in place, and for efficiency
    // does not conserve the nullspace
    std::vector<std::vector<int>> nullspace = nullspace_sparse(a_rref, n);
    x = back_substitution_sparse(a_rref, n, b_p);
    std::tuple<std::vector<int>, int, bool> min_res = minimize_weight_sparse(x, nullspace, max_dof);

    x = std::get<0>(min_res);
    int dof = std::get<1>(min_res);
    bool optimized = std::get<2>(min_res);
    return std::make_tuple(true, x, dof, optimized);
  } 
  return std::make_tuple(false, std::vector<int>(0), 0, false);
}
