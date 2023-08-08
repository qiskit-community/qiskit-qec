#include "errorpropagator.h"
#include <iostream>

/*
 * Construct a new error propagator object.
 *
 * qreg_size: number of qubits
 * creg_size: number of classical bits
 */
ErrorPropagator::ErrorPropagator(int qreg_size, int creg_size)
{
  _qreg_size = qreg_size;
  _creg_size = creg_size;
  qubit_array = std::vector<int>(2 * qreg_size);
  clbit_array = std::vector<int>(creg_size);
}

/*
 * Return the clbit_array.
 */
std::vector<int> ErrorPropagator::get_cbits(void) const
{
  return clbit_array;
}

/*
 * Return the qubit_array.
 */
std::vector<int> ErrorPropagator::get_qubit_array(void) const
{
  return qubit_array;
  ;
}

/*
 * Return the qubit_array as a Pauli string.
 */
std::string ErrorPropagator::get_qubits(void) const
{
  std::string pauli;
  for (int i = 0; i < _qreg_size; i++)
  {
    if (qubit_array[i] == 0 &&
        qubit_array[i + _qreg_size] == 0)
      pauli.append("i");
    if (qubit_array[i] == 0 &&
        qubit_array[i + _qreg_size] == 1)
      pauli.append("z");
    if (qubit_array[i] == 1 &&
        qubit_array[i + _qreg_size] == 0)
      pauli.append("x");
    if (qubit_array[i] == 1 &&
        qubit_array[i + _qreg_size] == 1)
      pauli.append("y");
  }
  return pauli;
}

/*
 * Return size of quantum register.
 */
int ErrorPropagator::get_qreg_size(void) const
{
  return _qreg_size;
}

/*
 * Return size of classical register.
 */
int ErrorPropagator::get_creg_size(void) const
{
  return _creg_size;
}

/*
 * Return size of loaded circuit (number of operations).
 */
size_t ErrorPropagator::get_circuit_size(void) const
{
  return _encoded_circ.size();
}

/*
 * Test if q is a valid qubit index.
 *
 * q: qubit index
 *
 * Throw an exception if the index is out of range.
 */
void ErrorPropagator::_range_check(int q)
{
  if (q < 0 || q >= _qreg_size)
    throw std::logic_error("qubit index out of range");
}

/*
 * Load an encoded stabilizer circuit.
 *
 * Each element of the vector encoded_circ is one of:
 * (op_h, q), (op_s, q), (op_x, q), (op_y, q), (op_z, q), (op_id, q),
 * (op_cx, qc, qt), (op_reset, q), (op_barrier), (op_measure, q, c)
 */
void ErrorPropagator::load_circuit(int qreg_size, int creg_size,
                                   std::vector<std::vector<int>> encoded_circ)
{
  _qreg_size = qreg_size;
  _creg_size = creg_size;
  qubit_array = std::vector<int>(2 * qreg_size);
  clbit_array = std::vector<int>(creg_size);
  _encoded_circ = encoded_circ;
  // Sanity check on the encoded circuit data
  for (auto i = encoded_circ.begin(); i != encoded_circ.end(); i++)
  {
    std::vector<int> op = (*i);
    switch (static_cast<opcode>(op[0]))
    {
    case op_h:
    case op_s:
    case op_x:
    case op_y:
    case op_z:
    case op_id:
    case op_reset:
      if (op.size() != 2)
        throw std::logic_error("bad 1-qubit instruction");
      _range_check(op[1]);
      break;
    case op_measure:
      if (op.size() != 3)
        throw std::logic_error("bad measure instruction");
      _range_check(op[1]);
      if (op[2] < 0 || op[2] >= _creg_size)
        throw std::logic_error("bit index out of range");
      break;
    case op_cx:
      if (op.size() != 3)
        throw std::logic_error("bad cx instruction");
      _range_check(op[1]);
      _range_check(op[2]);
      break;
    case op_barrier:
      if (op.size() != 1)
        throw std::logic_error("bad barrier instruction");
      break;
    default:
      throw std::logic_error("bad opcode");
    }
  }
}

/*
 * Apply a Pauli error to a set of qubits.
 *
 * qubits: qubits to apply the error to
 * paulis: Pauli errors to apply, string of 'x', 'y', and 'z' (lowercase)
 */
void ErrorPropagator::apply_error(std::vector<int> qubits, std::string paulis)
{
  if (qubits.size() != paulis.length())
    throw std::logic_error("qubits and pauli length disagree");
  for (int i = 0; i < qubits.size(); i++)
  {
    int q = qubits[i];
    if (paulis[i] == 'x' || paulis[i] == 'y')
      qubit_array[q] ^= 1;
    if (paulis[i] == 'y' || paulis[i] == 'z')
      qubit_array[q + _qreg_size] ^= 1;
  }
}

/*
 * Apply a single-qubit Pauli error to a qubit.
 *
 * q: qubit to apply the error to
 * pauli: Pauli error to apply, one of 'x', 'y', or 'z' (lowercase)
 *
 * Don't do any input checking.
 */
void ErrorPropagator::_apply_error_1(int q, std::string pauli)
{
  if (pauli[0] == 'x' || pauli[0] == 'y')
    qubit_array[q] ^= 1;
  if (pauli[0] == 'y' || pauli[0] == 'z')
    qubit_array[q + _qreg_size] ^= 1;
}

/*
 * Apply a two-qubit Pauli error to a pair of qubits.
 *
 * q: first qubit to apply the error to
 * r: second qubit to apply the error to
 * pauli: two-qubit Pauli error to apply, string containing 'x', 'y', and 'z'
 *
 * Don't do any input checking.
 */
void ErrorPropagator::_apply_error_2(int q, int r, std::string paulis)
{
  if (paulis[0] == 'x' || paulis[0] == 'y')
    qubit_array[q] ^= 1;
  if (paulis[0] == 'y' || paulis[0] == 'z')
    qubit_array[q + _qreg_size] ^= 1;
  if (paulis[1] == 'x' || paulis[1] == 'y')
    qubit_array[r] ^= 1;
  if (paulis[1] == 'y' || paulis[1] == 'z')
    qubit_array[r + _qreg_size] ^= 1;
}

/*
 * Apply a controlled-NOT gate.
 */
void ErrorPropagator::cx(int qc, int qt)
{
  qubit_array[qt] ^= qubit_array[qc];
  qubit_array[qc + _qreg_size] ^= qubit_array[qt + _qreg_size];
}

/*
 * Apply a Hadamard gate.
 */
void ErrorPropagator::h(int q)
{
  int z = qubit_array[q + _qreg_size];
  qubit_array[q + _qreg_size] = qubit_array[q];
  qubit_array[q] = z;
}

/*
 * Apply a Phase gate.
 */
void ErrorPropagator::s(int q)
{
  if (qubit_array[q] == 1)
    qubit_array[q + _qreg_size] ^= 1;
}

/*
 * Apply a Reset operation.
 */
void ErrorPropagator::reset(int q)
{
  qubit_array[q] = 0;
  qubit_array[q + _qreg_size] = 0;
}

/*
 * Apply a measurement.
 *
 * Return the measurement outcome.
 */
int ErrorPropagator::measure(int q, int c)
{
  clbit_array[c] = qubit_array[q];
  return static_cast<int>(clbit_array[c]);
}

/*
 * Propagate a set of Pauli faults through the loaded circuit.
 */
std::vector<int> ErrorPropagator::propagate(std::vector<int> icomb,
                                            std::vector<std::string> error)
{
  qubit_array = std::vector<int>(2 * _qreg_size, 0);
  clbit_array = std::vector<int>(_creg_size, 0);
  // Record which operations fail and the corresponding index into icomb
  std::vector<int> fail(_encoded_circ.size(), -1);
  for (int j = 0; j < icomb.size(); j++)
    fail[icomb[j]] = j;
  // Iterate over the encoded circuit
  for (int i = 0; i < _encoded_circ.size(); i++)
  {
    std::vector<int> op = _encoded_circ[i];
    switch (static_cast<opcode>(op[0]))
    {
    case op_h:
      h(op[1]);
      if (fail[i] >= 0)
        _apply_error_1(op[1], error[fail[i]]);
      break;
    case op_s:
      s(op[1]);
      if (fail[i] >= 0)
        _apply_error_1(op[1], error[fail[i]]);
      break;
    case op_x:
    case op_y:
    case op_z:
    case op_id:
      if (fail[i] >= 0)
        _apply_error_1(op[1], error[fail[i]]);
      break;
    case op_reset:
      reset(op[1]);
      if (fail[i] >= 0)
        _apply_error_1(op[1], error[fail[i]]);
      break;
    case op_measure:
      if (fail[i] >= 0)
        _apply_error_1(op[1], error[fail[i]]);
      measure(op[1], op[2]);
      break;
    case op_cx:
      cx(op[1], op[2]);
      if (fail[i] >= 0)
        _apply_error_2(op[1], op[2], error[fail[i]]);
      break;
    case op_barrier:
      break;
    default:
      throw std::logic_error("bad opcode");
    }
  }
  return clbit_array;
}
