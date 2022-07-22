from qiskit_qec.operators.base_pauli import BasePauli
from qiskit_qec.operators.pauli import Pauli
from qiskit_qec.utils import pauli_rep


BasePauli.set_syntax(pauli_rep.LATEX_SYNTAX)
x_op = Pauli("X")
print(x_op)
