#include "faultenumerator.h"

/*
 * Construct a fault enumerator.
 */
FaultEnumerator::FaultEnumerator(int order, int qreg_size, int creg_size,
                 std::vector<std::vector<int> > & encoded_circ,
                 std::vector<int> &faulty_ops_indices,
                 std::vector<std::string> &faulty_ops_labels,
                 std::vector<std::vector<std::string> > &faulty_ops_pauli_errors)
{
  _order = order;
  _qreg_size = qreg_size;
  _creg_size = creg_size;
  _encoded_circ = encoded_circ;
  _faulty_ops_indices = faulty_ops_indices;
  _faulty_ops_labels = faulty_ops_labels;
  _faulty_ops_pauli_errors = faulty_ops_pauli_errors;
  _ep = new ErrorPropagator(_qreg_size, _creg_size);
  _ep->load_circuit(_qreg_size, _creg_size, _encoded_circ);
  _index = 0;
  _state.clear();
  _done = false;
}

/*
 * Destroy fault enumerator.
 */
FaultEnumerator::~FaultEnumerator(void)
{
  delete _ep;
}

/*
 * Reset the enumerator state.
 */
void FaultEnumerator::reset(void)
{
  _index = 0;
  _state.clear();
  _done = false;
}

/*
 * Get the current index.
 */
long FaultEnumerator::get_index(void) const {
  return _index;
}

/*
 * Get the current state.
 */
std::vector<int> FaultEnumerator::get_state(void) const {
  return _state;
}

/*
 * Return the "done" flag.
 */
bool FaultEnumerator::done(void) const {
  return _done;
}

/*
 * Enumerate fault paths in blocks.
 *
 * blocksize = minimum number of fault paths to return
 * (may return somewhat more)
 */
std::vector<FaultPath> FaultEnumerator::enumerate(int blocksize)
{
  std::vector<FaultPath> block;
  if(_done) return block;
  // Sanity check on the state
  if(_state.size() != 0 && _state.size() != _order)
    throw std::logic_error("enumerator state must be empty or size _order");
  // Reset if the state is empty
  if(_state.size() == 0) { _index = 0; _done = false; }
  Combinations c(_faulty_ops_indices.size(), _order);
  long initial_index = _index;
  while(c.next_combination(_state))
  {
    std::vector<int> indices;
    std::vector<std::string> labels;
    std::vector<std::vector<std::string> > errors;
    std::vector<int> ni;
    for(int i=0; i<_order; i++) {
      indices.push_back(_faulty_ops_indices[_state[i]]);
      labels.push_back(_faulty_ops_labels[_state[i]]);
      errors.push_back(_faulty_ops_pauli_errors[_state[i]]);
      ni.push_back(_faulty_ops_pauli_errors[_state[i]].size());
    }
    ProductIterator p(ni);
    std::vector<int> pstate;
    while(p.next_element(pstate))
    {
      std::vector<std::string> error;
      for(int j=0; j<_order; j++)
        error.push_back(errors[j][pstate[j]]); 
      FaultPath result = std::make_tuple(_index, labels, error,
                                         _ep->propagate(indices, error));
      block.push_back(result);
      _index++;
    }
    if(_index - initial_index >= blocksize)
      return block;
  }
  _done = true;
  return block;
}

