from qiskit.quantum_info.operators import PauliList

pt = PauliList(['X', 'Y', '-X', 'I', 'I', 'Z', 'X', 'iZ'])
unique = pt.unique()
print(unique)