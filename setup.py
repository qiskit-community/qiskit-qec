# This code is part of Qiskit.
#
# (C) Copyright IBM 2021.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

# currently used for backwards compatibility
# and to allow for editable code for development mode
import os
import re
from glob import glob
from setuptools import find_packages, setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

ext_modules = [
    Pybind11Extension(
        "qiskit_qec.extensions.compiledextension",
        sorted(glob("qiskit_qec/extensions/*.cpp")),
    ),
]

with open("requirements.txt") as f:
    REQUIREMENTS = f.read().splitlines()

# Read long description from README.
README_PATH = os.path.join(os.path.abspath(os.path.dirname(__file__)), "README.md")
with open(README_PATH) as readme_file:
    README = re.sub(
        "<!--- long-description-skip-begin -->.*<!--- long-description-skip-end -->",
        "",
        readme_file.read(),
        flags=re.S | re.M,
    )

setup(
    name="qiskit_qec",
    version="v0.0.0-alpha",
    description="Quantum Error Correcting Package",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/Qiskit/qiskit-qec",
    author="Qiskit Development Team",
    author_email="hello@qiskit.org",
    license="Apache 2.0",
    python_requires=">=3.6",
    include_package_data=True,
    install_requires=(REQUIREMENTS,),
    packages=find_packages(exclude=["test*", "qiskit_qec/extensions"]),
    cmdclass={"build_ext": build_ext},
    ext_modules=ext_modules,
)
