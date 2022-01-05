import numpy as np

from qiskit_qec.geometry.model.qubit_count import QubitCount
from qiskit_qec.geometry.model.qubit_data import QubitData
from qiskit_qec.geometry.model.vertex import Vertex
from qiskit_qec.geometry.model.edge import Edge
from qiskit_qec.geometry.model.wireframe import WireFrame
from qiskit_qec.geometry.model.face import Face
from qiskit_qec.geometry.model.shell import Shell

X_OPERATOR="X"
Z_OPERATOR="Z"
Y_OPERATOR="Y"

qubit_count = QubitCount()
qubit_data = QubitData()

v1 = Vertex(np.array([0,1]))
qubit_data.qubit[v1.id] = qubit_count.new_qubit()
qubit_data.operator[v1.id] = X_OPERATOR

v2=Vertex(np.array([1,0]))
qubit_data.qubit[v2.id] = qubit_count.new_qubit()
qubit_data.operator[v2.id] = X_OPERATOR

v3=Vertex(np.array([1,1]))
qubit_data.qubit[v3.id] = qubit_count.new_qubit()
qubit_data.operator[v3.id] = X_OPERATOR

v4=Vertex(np.array([0,0]))
qubit_data.qubit[v4.id] = qubit_count.new_qubit()
qubit_data.operator[v4.id] = X_OPERATOR


e1=Edge([v4,v2])
e2=Edge([v2,v3])
e3=Edge([v3,v1])
e4=Edge([v1,v4])

wf1=WireFrame([e1,e2,e3,e4])
f1=Face(wf1)

v5=Vertex(np.array([1,0]))
qubit_data.qubit[v5.id] = qubit_count.new_qubit()
qubit_data.operator[v5.id] = Z_OPERATOR

v6=Vertex(np.array([1,1]))
qubit_data.qubit[v6.id] = qubit_count.new_qubit()
qubit_data.operator[v6.id] = Z_OPERATOR

v7=Vertex(np.array([2,0]))
qubit_data.qubit[v7.id] = qubit_count.new_qubit()
qubit_data.operator[v7.id] = Z_OPERATOR

v8=Vertex(np.array([2,1]))
qubit_data.qubit[v8.id] = qubit_count.new_qubit()
qubit_data.operator[v8.id] = Z_OPERATOR


e5=Edge([v5,v6])
e6=Edge([v6,v8])
e7=Edge([v8,v7])
e8=Edge([v7,v5])

wf2=WireFrame([e5,e5,e6,e7])
f2=Face(wf2)

s=Shell([f1,f2])

