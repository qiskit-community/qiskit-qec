from plerco.pec_python.qec_code.rssc import RSSC
from plerco.pec_python.qec_circuit.rssc_circuit import RSSCCircuit

import time
import cProfile

# pauli_fault_events_error_propagation_older:3:   3.5869 seconds
# pauli_fault_events_error_propagation_old:3:     0.8795 seconds
# pauli_fault_events_error_propagation:3:         0.6654 seconds
# pauli_fault_events_error_propagation_cpp:3:     0.5256 seconds
# pauli_fault_events_error_propagation_older:5:   74.7575 seconds
# pauli_fault_events_error_propagation_old:5:     21.6513 seconds
# pauli_fault_events_error_propagation:5:         16.0609 seconds
# pauli_fault_events_error_propagation_cpp:5:     11.3993 seconds
# pauli_fault_events_error_propagation_old:7:     185.2620 seconds
# pauli_fault_events_error_propagation:7:         133.2774 seconds
# pauli_fault_events_error_propagation_cpp:7:      93.0225 seconds
# pauli_fault_events_error_propagation:9:         630.1768 seconds
# pauli_fault_events_error_propagation_cpp:9:     440.5906 seconds

d = 7
code = RSSC(d)
gen = RSSCCircuit(code)
circ = gen.syndrome_measurement(d, "zx", "z")
tic = time.perf_counter()
cProfile.run("ret = [x for x in pauli_fault_events_error_propagation_cpp(circ, 1)]")
toc = time.perf_counter()
print(f"{toc - tic:0.4f} seconds")
