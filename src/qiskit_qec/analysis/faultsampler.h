#ifndef __FaultSampler__
#define __FaultSampler__

#include <tuple>
#include <utility>
#include <vector>
#include <map>
#include <string>
#include <random>

#include "errorpropagator.h"

/*
 * Data for a fault path is stored as a tuple.
 *
 * index = integer index of fault path
 * labels = vector of labels for each failed operation
 * errors = vector of Pauli errors for each failed operation
 * outcome = measurement outcomes for this fault path
 */
typedef std::tuple<long,
                   std::vector<std::string>,
                   std::vector<std::string>,
                   std::vector<int> > FaultPath;

typedef std::vector<std::pair<std::string, double> > PauliWeight;
typedef std::map<std::string, PauliWeight> LabelToPauliWeightMap;


class FaultSampler
{
  private:
    // size of quantum register for loaded circuit
    int _qreg_size;
    // size of classical register for loaded circuit
    int _creg_size;
    // loaded circuit, vector of (opcode, arg, ...)
    std::vector<std::vector<int> > _encoded_circ;
    // Operations with potential to fail (i.e., fault operations).
    // Elements of vector are indices into _encoded_circ.
    std::vector<int> _faulty_ops_indices;
    // Vector of labels for faulty operations, in order
    // of elements of faulty_ops_indices.
    std::vector<std::string> _faulty_ops_labels;
    // Map from labels to vectors of pairs of Pauli errors and their
    // corresponding weights, i.e., <"ix", 0.5>
    LabelToPauliWeightMap _label_to_pauli_weight;
    // Map from labels to error probabilities
    std::map<std::string, double> _label_to_error_probability;
    // Pointer to our error propagator
    ErrorPropagator *_ep;
    // Pointer to random number generator and uniform distribution
    std::mt19937 *_generator;
    std::uniform_real_distribution<double> *_unif;
    // Map from labels to discrete distributions
    std::map<std::string, std::discrete_distribution<int> > _label_to_dist;

  public:
    FaultSampler(int qreg_size, int creg_size, std::vector<std::vector<int> > &encoded_circ,
                 std::vector<int> &faulty_ops_indices, std::vector<std::string> &faulty_ops_labels,
                 LabelToPauliWeightMap &label_to_pauli_weight,
                 std::map<std::string, double> &label_to_error_probability,
                 unsigned int seed);
    ~FaultSampler(void);
    std::vector<FaultPath> sample(int blocksize);
};
#endif

