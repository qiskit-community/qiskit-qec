#include "faultsampler.h"

/*
 * Construct a fault sampler.
 */
FaultSampler::FaultSampler(int qreg_size, int creg_size, std::vector<std::vector<int> > &encoded_circ,
  std::vector<int> &faulty_ops_indices, std::vector<std::string> &faulty_ops_labels,
  LabelToPauliWeightMap &label_to_pauli_weight,
  std::map<std::string, double> &label_to_error_probability,
  unsigned int seed)
{
  _qreg_size = qreg_size;
  _creg_size = creg_size;
  _encoded_circ = encoded_circ;
  _faulty_ops_indices = faulty_ops_indices;
  _faulty_ops_labels = faulty_ops_labels;
  _label_to_pauli_weight = label_to_pauli_weight;
  _label_to_error_probability = label_to_error_probability;
  _ep = new ErrorPropagator(_qreg_size, _creg_size);
  _ep->load_circuit(_qreg_size, _creg_size, _encoded_circ);
  _generator = new std::mt19937(seed);
  _unif = new std::uniform_real_distribution<double>(0.0, 1.0);
  for(LabelToPauliWeightMap::iterator i=_label_to_pauli_weight.begin();
      i != _label_to_pauli_weight.end(); i++) {
    std::vector<double> probs;
    for(PauliWeight::iterator j=i->second.begin(); j!=i->second.end(); j++) {
      probs.push_back(j->second);
    }
    _label_to_dist[i->first] =
      std::discrete_distribution<int>(probs.begin(), probs.end());
  }
}

/*
 * Destroy fault sampler.
 */
FaultSampler::~FaultSampler(void)
{
  delete _ep;
  delete _generator;
  delete _unif;
}

/*
 * Sample random fault paths in blocks.
 *
 * blocksize = number of fault paths to return
 */
std::vector<FaultPath> FaultSampler::sample(int blocksize)
{
  std::vector<FaultPath> block;
  for(int i=0; i<blocksize; i++) {
    std::vector<int> locations;
    // Sample faulty operations
    for(int j=0; j<_faulty_ops_labels.size(); j++) {
      double r = (*_unif)(*_generator);
      std::string label = _faulty_ops_labels[j];
      if (r < _label_to_error_probability[label])
        locations.push_back(j);
    }
    std::vector<int> indices;
    std::vector<std::string> labels;
    std::vector<std::string> error;
    for(int j=0; j<locations.size(); j++) {
      int index = _faulty_ops_indices[locations[j]];
      std::string label = _faulty_ops_labels[locations[j]];
      indices.push_back(index);
      labels.push_back(label);
      int p = _label_to_dist[label](*_generator);
      error.push_back(_label_to_pauli_weight[label][p].first);
    }
    FaultPath result = std::make_tuple(-1, labels, error,
                                       _ep->propagate(indices, error));
    block.push_back(result);
  }
  return block;
}
