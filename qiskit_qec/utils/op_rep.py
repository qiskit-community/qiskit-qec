# This code is part of Qiskit.
#
# (C) Copyright IBM 2017, 2020
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.
"""
N-qubit Pauli Operator/Instruction Conversion Module
"""

# pylint: disable=invalid-name,anomalous-backslash-in-string
# pylint: disable=bad-docstring-quotes  # for deprecate_function decorator

from typing import Tuple, Union, Any

import numpy as np
from qiskit.exceptions import QiskitError
from qiskit.circuit.barrier import Barrier
from qiskit.circuit.delay import Delay
from qiskit.circuit import Instruction
from qiskit.circuit.library.generalized_gates import PauliGate
from qiskit.circuit.library.standard_gates import IGate, XGate, YGate, ZGate
from qiskit.quantum_info.operators.scalar_op import ScalarOp
from qiskit.circuit import Gate
from qiskit.circuit import QuantumCircuit
from qiskit_qec.utils import pauli_rep


# ----------------------------------------------------------------------


def scalar_op2symplectic(
    op: ScalarOp, output_encoding: str = pauli_rep.DEFAULT_EXTERNAL_PHASE_ENCODING
) -> Tuple[np.ndarray, Union[np.array, Any]]:
    """Convert a ScalarOp to symplectic representation with phase.

    TODO: Allow this to work on arrays of ScalarOps

    Args:
        op: Input scalarOp
        output_encoding: Phase encoding to use to encode phase from ScalarOp.
            Default is INTERNAL_PHASE_ENCODING='-i'

    Raises:
        QiskitError: Operator is not an N-qubit identity

    Returns:
        matrix, phase_exponent: GF(2) symplectic matrix and phase_exponent
        representing ScalarOp
    """
    if op.num_qubits is None:
        raise QiskitError(f"{op} is not an N-qubit identity")
    matrix = np.zeros(shape=(1, 2 * op.num_qubits), dtype=np.bool_)
    phase_exp = pauli_rep.cpx2exp(op.coeff, output_encoding=output_encoding)
    return matrix, phase_exp


# ----------------------------------------------------------------------


def gate2symplectic(
    gate: Gate, encoding: str = pauli_rep.INTERNAL_PAULI_ENCODING
) -> Tuple[np.ndarray, Union[np.array, Any]]:
    """Converts a Pauli gate to a symplectic matrix with phase

    Args:
        gate: Gate
        encoding (optional): Pauli encoding to encode symplectic matrix with phase.
            Defaults to DEFAULT_EXTERNAL_PAULI_ENCODING='-iYZX';

    Raises:
        QiskitError: Invalid Pauli instruction

    Returns:
        matrix, phase_exp: phase exponent and symplectic matrix
    """
    if isinstance(gate, PauliGate):
        return pauli_rep.str2symplectic(gate.params[0], output_encoding=encoding)
    if isinstance(gate, IGate):
        return pauli_rep.str2symplectic("I", output_encoding=encoding)
    if isinstance(gate, XGate):
        return pauli_rep.str2symplectic("X", output_encoding=encoding)
    if isinstance(gate, YGate):
        return pauli_rep.str2symplectic("Y", output_encoding=encoding)
    if isinstance(gate, ZGate):
        return pauli_rep.str2symplectic("Z", output_encoding=encoding)
    raise QiskitError("Invalid Pauli instruction.")


# ----------------------------------------------------------------------


def instrs2symplectic(instr: Union[Instruction, QuantumCircuit]):
    """Convert a Pauli circuit to BasePauli data."""
    # Try and convert single instruction
    if isinstance(instr, (PauliGate, IGate, XGate, YGate, ZGate)):
        return gate2symplectic(instr)
    if isinstance(instr, Instruction):
        # Convert other instructions to circuit definition
        if instr.definition is None:
            raise QiskitError(f"Cannot apply Instruction: {instr.name}")
        # Convert to circuit
        instr = instr.definition

    from qiskit_qec.operators.base_pauli import BasePauli
    from qiskit_qec.operators.pauli import Pauli

    # Initialize identity Pauli
    ret = Pauli(np.zeros((1, 2 * instr.num_qubits), dtype=bool), phase_exp=0)
    # Add circuit global phase if specified
    if instr.global_phase:
        ret._phase_exp = pauli_rep.cpx2exp(
            np.exp(1j * float(instr.global_phase)),
            output_encoding=pauli_rep.INTERNAL_PHASE_ENCODING,
        )
    # Recursively apply instructions
    for dinstr, qregs, cregs in instr.data:
        if cregs:
            raise QiskitError(f"Cannot apply instruction with classical registers: {dinstr.name}")
        if not isinstance(dinstr, (Barrier, Delay)):
            next_instr = BasePauli(*instrs2symplectic(dinstr))
            if next_instr is not None:
                qargs = [tup.index for tup in qregs]
                ret = ret.compose(next_instr, qargs=qargs)
    return ret.matrix, ret._phase_exp
