# Qiskit QEC



<!--- long-description-skip-begin -->[QEC](https://github.com/Qiskit/qiskit-qec)<!--- long-description-skip-end -->

[![License](https://img.shields.io/github/license/Qiskit/qiskit-terra.svg?style=popout-square)](https://opensource.org/licenses/Apache-2.0)

**Qiskit Framework for Quantum Error Correction** is an open-source framework for developers, experimentalist and theorists of Quantum Error Correction (QEC). The framework is in the development stage.

## Fast Access Links

- [Requirements](https://github.com/Qiskit/qiskit-qec/blob/develop/docs/Requirements.md)
- [Discussions of Requirements and Code Ideas](https://github.com/Qiskit/qiskit-qec/discussions)
- [Introduction to QEC papers](#introduction-to-qec-papers)
- List of known [QEC software packages](#publicly-accessible-qec-software-packages) (not Qiskit)

## Documentation

**Note**: Viewing LaTeX rendered on the pages below is possible if your use the Chrome browser and have the following MathJax enabled extension enabled (it is coded to work on github.com) : https://chrome.google.com/webstore/detail/mathjax-plugin-for-github/ioemnmodlmafdkllaclgeombjnmnbima/related


## Contributing 

Contributors to the framework are welcome. The following describes how the framework is being designed and built. The current approach is designed to enable a wide range of disciplines to contribute earlier and to enable better computer engineering earlier on. This approach has a development process and a scope control. The Development process is based on increasing details of requirements and code details. The basic process is outlined in the following figure:

![Development Process](./docs/_static/images/DevelopmentProcessv1.png?raw=true "Development Process")


QEC requirements, both general and specific, are first written down. Each requirement must then me given a detailed precise mathematical description. These descriptions are then use to develop various solution approaches and associated pseudo code. From here prototype code is developed. Once the code has been tested and fits within the agreed upon architectural framework it becomes framework code. A similar process occurs for the architectural side of the project.

Scope control is used initially to establish a good foundation.  The basic process is outlined in the following figure:

![Scope Control](./docs/_static/images/ScopeControlv1.png?raw=true "Scope Control")

Scope control is done by prioritizing which requirements get worked on first. The priorities are determined by selected core QEC or architecure requirements or via requirements required by a few selected projects.


## Ways to Contribute

1. Development QEC or Architectural requirements (general or detailed)
  - Requirements and details are curently stored [here](https://github.com/Qiskit/qiskit-qec/discussions/categories/requirement)
3. Developing approaches and pseudocode for requirements
4. Writing prototype code for requirements
5. Creating and testing unit tests
6. Review requirements, appraoches, pseudo code, prototype code, ...
7. Review and provide architectural recommendations
8. Development detailed documentation (User, Code, Detailed)

## Contributing Flow
1. Create/Choose an issue on github
2. Checkout a new branch named in the format: <issue-number>_<short-description>
3. Make all changes desired on that branch
4. Push that branch to github
5. Make a PR from your branch into **develop** branch. Assign it to yourself and request reviews from Drew Vandeth or Grace Harper

## Code Location

Not all code for qiskit-qec is currently located in this repository since Qiskit already has certain modules that contain some of the code necessary for the QEC framework. This included the quantum_info module in qiskit-terra for example. Below is list of PR's etc for some code under development. This will eventually moved to a better location:

### qiskit-terra.quantum_info

This extends the new PauliList class to have more than one external PauliRepresentation.
https://github.com/dsvandet/qiskit-terra/tree/qi/pauli-list


## Major Projects

There are currently two different major projects to develop the QEC framework. The first is the core CLI/API framework and the second is the visual code GUI.


## CLI/API Framework

## Visual Code GUI

## Introduction to QEC papers

The following provides a list of papers that are good introductions to QEC or various aspects of QEC
  
## Publicly Accessible QEC Software packages
  
* **stabnf**: This source code is the C implementation with a simple text-based user interface of an algorithm that writes stabilizer quantum circuits under normal form by Marc Bataille. Details of the work can be found in this article: [arxov.org:2107.00885](https://arxiv.org/abs/2012.09224). It is available from [github](https://github.com/marcbataille/stabilizer-circuits-normal-forms). Language: C

* **qecsim**: is a Python 3 package for simulating quantum error correction using stabilizer codes by David Tuckett. It provides access to all features via a command-line interface. It can also be used as a library via the fully-documented API. It includes many common codes, error models and decoders, and can be extended with additional components. It is available from [github repo](https://github.com/qecsim/qecsim). Language: Python and C. Status: Active

* **PyMatchingis**  a  fast  Python/C++  library  for  decoding  quantum  error  correcting  codes(QECC) using the Minimum Weight Perfect Matching (MWPM) decoder by Oscar Higgot.It  can  be  downloaded  from [github repo](https://github.com/oscarhiggott/PyMatching).   Language:Python and C++.  Status:  Active

* **qsurface**: is a simulation package for the surface code, and is designed to modularize 3 aspects of a surface code simulation developed by Mark Shui Hu at Delt. Available from [github repo](https://github.com/watermarkhu/qsurface). The documentation is minimal. Language Python. Status: Active

* **Sweep-Decoder-Boundaries**: is a C++ implementation of the sweep decoder for topological quantum codes with and without boundaries by _Michael Vasmer, Dan Browne, and Aleksander Kubica_. Available from [github repo](https://github.com/MikeVasmer/Sweep-Decoder-Boundaries). Language: C++. Status: Last update 6 months ago.
    
*  **Single_shot_3D_HGP**: is a single-shot decoder for three-dimensional homological product codes, using belief propagation with ordered statistics post-processing by _Armanda Quintavalle, Michael Vasmer, Joschka Roffe, and Earl Campbell_. Available from [github repo](https://github.com/MikeVasmer/single_shot_3D_HGP), Langauge: C++, Status: Active.
    
* **Bosonic Codes**: is code for the paper Quantum computing with rotation-symmetric bosonic codes by  _Arne Grimsmo, Joshua Combes, and Ben Q. Baragiola_. Available from [github repo](https://github.com/arnelg/arXiv-1901.08071). Language : Python and MATLAB, Status: No longer maintained.
    
* **DeepNeuralDecoder**: is a package that provides a general framework for training and using deep learning for fault tolerant protocols that use CSS codes by _Christopher Chamberland and Pooya Ronagh_. Available from [github repo](https://github.com/pooya-git/DeepNeuralDecoder). Language: Python and MatLab, Status: No longer maintained.
    
* **DeepQ Decoding**: is a set of reinforcement learning decoders for fault-tolerant quantum computation by _Ryan Sweke, Markus Kesselring, Evert van Nieuwenburg, and Jens Eisert_. Available from [github repo](https://github.com/R-Sweke/DeepQ-Decoding). Language: Juypter Notebook, Python, Status: No longer maintained.
 
* **Autotune**: is surface code software (or at least a stripped down version) by _Andrew Fowler_ and others. Autotune can still be download at _Adam Whiteside's_ Github repository [github repo](https://github.com/adamcw/autotune). The documentation is limited and not designed as a general toolkit for other researchers to use and expand. Language: C with some Python. Status: No longer maintained

* **OpenSurgery**: is for topological assemblies by _Alexandru Paler and Austin Fowler_ is a scalable tool for the preparation of circuits protected by the surface code operated through lattic surgery. Available from [github repo](https://github.com/alexandrupaler/opensurgery). Language: Python and Javascript. Status: Active

* **Quantomatic**: for diagrammatic proof assistant, meaning it provides machine-support for reasoning with diagrammatic languages from the Universities of Oxford and Edinburgh. Available from [github repo](http://quantomatic.github.io/). Language: ML, Java, Scala, ... Status: No longer maintained

* **QTop**: is an open-source python module for simulation and visualization of topological quantum codes by _Jocob Marks_. Available from [github repo](https://github.com/jacobmarks/QTop). Language: Python Status: No longer maintained
    
* **qasm-tools**: is open-source software for studying fault-tolerant quantum circuits by _Andrew Cross_. Available from [github repo](https://www.media.mit.edu/quanta/quanta-web/projects/qasm-tools/). Language: C Status: No longer maintained
