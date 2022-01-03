from plerco.pec_python.qec_code.rssc import RSSC
from plerco.pec_python.qec_circuit.rssc_circuit import RSSCCircuit
from plerco.pec_python.qec_decoder.decoder_utils.faultenumerator import FaultEnumerator

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
fe = FaultEnumerator(circ, method="propagator")
tic = time.perf_counter()
cProfile.run("ret = [x for x in fe.generate()]")
toc = time.perf_counter()
print(f"{toc - tic:0.4f} seconds")

# A little less than twice as slow as the cpp function.
"""
         352233496 function calls in 167.194 seconds

   Ordered by: standard name

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        1    0.012    0.012  167.194  167.194 <string>:1(<listcomp>)
        1    0.000    0.000  167.194  167.194 <string>:1(<module>)
     5846    0.005    0.000    0.005    0.000 dagnode.py:50(op)
 47500320   43.590    0.000   58.959    0.000 faultpath.py:115(cx)
 16625112   10.633    0.000   13.374    0.000 faultpath.py:122(h)
 18519470    8.589    0.000   11.659    0.000 faultpath.py:135(reset)
 18519470    9.694    0.000   12.670    0.000 faultpath.py:141(measure)
 47500320   15.964    0.000   75.411    0.000 faultpath.py:151(_faulty_cx)
 16625112    5.461    0.000   18.850    0.000 faultpath.py:164(_faulty_h)
 18519470    6.239    0.000   17.902    0.000 faultpath.py:178(_faulty_reset)
 18519470    6.883    0.000   19.559    0.000 faultpath.py:185(_faulty_measure)
   933042    0.091    0.000    0.091    0.000 faultpath.py:192(_faulty_barrier)
    28274   35.275    0.001  167.096    0.006 faultpath.py:196(propagate_faults)
     3578    0.010    0.000    0.015    0.000 faultpath.py:259(_node_name_label)
    28275    0.063    0.000  167.183    0.006 faultpath.py:333(generate)
     3611    0.000    0.000    0.000    0.000 faultpath.py:349(<lambda>)
     3578    0.002    0.000    0.002    0.000 faultpath.py:351(<listcomp>)
     3578    0.003    0.000    0.018    0.000 faultpath.py:352(<listcomp>)
     3578    0.001    0.000    0.001    0.000 faultpath.py:353(<listcomp>)
     3578    0.002    0.000    0.002    0.000 faultpath.py:354(<listcomp>)
148718166   24.166    0.000   24.166    0.000 faultpath.py:40(_range_check)
    28274    0.477    0.000    0.496    0.000 faultpath.py:45(apply_error)
        1    0.000    0.000  167.194  167.194 {built-in method builtins.exec}
   113096    0.017    0.000    0.017    0.000 {built-in method builtins.len}
        1    0.000    0.000    0.000    0.000 {method 'disable' of '_lsprof.Profiler' objects}
    28274    0.019    0.000    0.019    0.000 {method 'index' of 'tuple' objects}


167.2057 seconds
"""
