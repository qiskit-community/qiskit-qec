from qiskit_qec.codes.subsystemcodes import SubSystemCode
from qiskit_qec.structures.gauge import GaugeGroup
from qiskit_qec.operators.pauli_list import PauliList

if __name__ == "__main__":

    pl1 = PauliList("xiiiiz".upper())
    pl2 = PauliList("ixiiiz".upper())
    pl3 = PauliList("zzzzzx".upper())
    pl4_css = PauliList("xxxxxx".upper())

    ssc = SubSystemCode(GaugeGroup(pl1))
    print(pl1)
    print(pl1 == ssc.gauge_group.generators)
