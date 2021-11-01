"""Setup file for Qiskit QEC."""
from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

with open('requirements.txt') as fp:
    install_requires = fp.read()

setup(name='qiskit_qec',
      version='0.1',
      description='cool lines',
      long_description=long_description,
      url='https://github.com/Qiskit/qiskit-qec',
      author='nerds',
      author_email='grace.harper@ibm.com',
      license='MIT',
      packages=find_packages(),
      install_requires=install_requires,
      zip_safe=False)
