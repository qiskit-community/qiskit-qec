Beginners guide
===============

## Content
1. [About project](#about-project)
2. [Installation](#installing-qiskit-qec)
3. [Usage](#running-examples)

----

### About project

**Qiskit Framework for Quantum Error Correction** is an open-source framework for developers, experimentalist and theorists of Quantum Error Correction (QEC).

---

### Installing Qiskit QEC

Refer to [installation guide](./installation.md).

----

### Running examples

```python
from qiskit_qec.codes.codebuilders.triangular_color_code_builder import TriangularColorCodeBuilder

code = TriangularColorCodeBuilder(d=3).build()
code.draw(face_colors=True, show_index=True)
```
