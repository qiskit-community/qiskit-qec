#ifndef __ErrorPropagator__
#define __ErrorPropagator__

#include <vector>
#include <string>
#include <algorithm>
#include <exception>

class ErrorPropagator
{
  private:
    int _qreg_size;
    int _creg_size;
    std::vector<int> qubit_array;
    std::vector<int> clbit_array;
    std::vector<std::vector<int> > _encoded_circ;
    enum opcode { op_h, op_s, op_x, op_y, op_z, op_cx,
                  op_id, op_reset, op_measure, op_barrier };

    void _range_check(int q);
    void _apply_error_1(int q, std::string pauli);
    void _apply_error_2(int q, int r, std::string paulis);

  public:
    ErrorPropagator(int qreg_size, int creg_size);
    void load_circuit(int qreg_size, int creg_size,
                      std::vector<std::vector<int> > encoded_circ);
    void apply_error(std::vector<int> qubits, std::string paulis);
    void cx(int qc, int qt);
    void h(int q);
    void s(int q);
    void reset(int q);
    int measure(int q, int c);
    std::vector<int> propagate(std::vector<int> icomb,
                               std::vector<std::string> error);
    std::vector<int> get_cbits(void) const;
    std::string get_qubits(void) const;
    std::vector<int> get_qubit_array(void) const;
    int get_qreg_size(void) const;
    int get_creg_size(void) const;
    size_t get_circuit_size(void) const;
};
#endif

