# Qiskit QEC
[![License](https://img.shields.io/github/license/Qiskit/qiskit-terra.svg?style=popout-square)](https://opensource.org/licenses/Apache-2.0)

**Qiskit Framework for Qauntum Error Correction** is an open-source framework for developers, experimentalist and theorists of Quantum Error Correction (QEC). The framework is in the development stage.

## Fast Access Links

- [Requirements](https://github.com/Qiskit/qiskit-qec/blob/develop/docs/Requirements.md)
- [Discussions of Requirements and Code Ideas](https://github.com/Qiskit/qiskit-qec/discussions)
- [Introduction to QEC papers](#introduction-to-qec-papers)

## Contributing 

Contributors to the framework are welcome. The following describes how the framework is being designed and built. The current approach is designed to enable a wide range of disciplines to contribute earlier and to enable better computer engineering earlier on. This approach has a development process and a scope control. The Development process is based on increasing details of requirements and code details. The basic process is outlined in the following figure:

![Development Process](https://github.com/Qiskit/qiskit-qec/blob/develop/docs/images/DevelopmentProcessv1.png?raw=true "Development Process")


QEC requirements, both general and specific, are first written down. Each requirement must then me given a detailed precise mathematical description. These descriptions are then use to develop various solution approaches and associated pseudo code. From here prototype code is developed. Once the code has been tested and fits within the agreed upon architectural framework it becomes framework code. A similar process occurs for the architectural side of the project.

Scope control is used initially to establish a good foundation.  The basic process is outlined in the following figure:

![Scope Control](https://github.com/Qiskit/qiskit-qec/blob/develop/docs/images/ScopeControlv1.png?raw=true "Scope Control")

Scope control is done by prioritizing which requirements get worked on first. The priorities are determined by selected core QEC or architecure requirements or via requirements required by a few selected projects.



## Ways to Contribute

1. Development QEC or Architectural requirements (general or detailed)
  - Requirements and details are curently stored [here](https://github.com/Qiskit/qiskit-qec/blob/develop/docs/Requirements.md)
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


