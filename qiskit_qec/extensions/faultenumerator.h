#ifndef __FaultEnumerator__
#define __FaultEnumerator__

#include <tuple>
#include <vector>
#include <string>
#include <algorithm>
#include <exception>

#include "combinations.h"
#include "productiterator.h"
#include "errorpropagator.h"

/*
 * Data for a fault path is stored as a tuple.
 *
 * index = integer index of fault path
 * labels = vector of labels for each failed operation
 * errors = vector of Pauli errors for each failed operation
 * outcome = measurement outcomes for this fault path
 */
typedef std::tuple<long,
                   std::vector<std::string>,
                   std::vector<std::string>,
                   std::vector<int> > FaultPath;

class FaultEnumerator
{
  private:
    // number of faults to insert for each fault path
    int _order;
    // size of quantum register for loaded circuit
    int _qreg_size;
    // size of classical register for loaded circuit
    int _creg_size;
    // loaded circuit, vector of (opcode, arg, ...)
    std::vector<std::vector<int> > _encoded_circ;
    // Operations with potential to fail (i.e., fault operations).
    // Elements of vector are indices into _encoded_circ.
    std::vector<int> _faulty_ops_indices;
    // Vector of labels for faulty operations, in order
    // of elements of faulty_ops_indices.
    std::vector<std::string> _faulty_ops_labels;
    // Vector of vectors of Pauli errors for faulty
    //operations, i.e., "ix", in order of elements
    // of faulty_ops_indices.
    std::vector<std::vector<std::string> > _faulty_ops_pauli_errors;;
    // Pointer to our error propagator
    ErrorPropagator *_ep;
    // Current fault path index
    long _index;
    // Current enumerator state
    std::vector<int> _state;
    // Enumerator completion flag
    bool _done;

  public:
    FaultEnumerator(int order, int qreg_size, int creg_size, std::vector<std::vector<int> > &encoded_circ,
                    std::vector<int> &faulty_ops_indices, std::vector<std::string> &faulty_ops_labels,
                    std::vector<std::vector<std::string> > &faulty_ops_pauli_errors);
    ~FaultEnumerator(void);
    std::vector<FaultPath> enumerate(int blocksize);
    void reset(void);
    long get_index(void) const;
    std::vector<int> get_state(void) const;
    bool done(void) const;
};
#endif

