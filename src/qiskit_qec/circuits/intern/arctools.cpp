#include "arctools.h"
#include <iostream>

bool is_cluster_neutral(
    std::vector<std::tuple<int, int, int, bool>> nodes, bool ignore_extras, bool minimal,
    std::map<std::tuple<int, int>, std::set<int>> cycle_dict, std::vector<std::tuple<int, int>> link_graph, std::map<int, std::vector<int>> link_neighbors, std::map<int, int> extras,
    bool linear
    ) {
    if (linear) {
        return nodes.size()%2==0;
    } else {
        std::vector<int> output = check_nodes(nodes, ignore_extras, minimal, cycle_dict, link_graph, link_neighbors, extras);
        bool no_boundary = true;
        for (int j = 2; j < output.size(); j++) {
            no_boundary = no_boundary and (extras[output[j]]<2);
        };
        return (output[0]==1) and no_boundary;
    }
};

std::vector<int> check_nodes(
    std::vector<std::tuple<int, int, int, bool>> nodes, bool ignore_extras, bool minimal,
    std::map<std::tuple<int, int>, std::set<int>> cycle_dict, std::vector<std::tuple<int, int>> link_graph, std::map<int, std::vector<int>> link_neighbors, std::map<int, int> extras
    ) {
    
    // output[0] is neutral (as int), output[1] is num_errors, rest is list of given extras
    std::vector<int> output;

    // we convert to flat nodes, which are a std::tuple with (q0, q1, boundary/logical)
    // here boundary/logical is 0 if the node is neither boundary or logical, 1 for logical, 2 for bounary and 3 for both
    // if we have an even number of corresponding nodes, they cancel
    std::map<std::tuple<int, int, int>, int> node_counts;
    for (auto & node : nodes) {
        node_counts[std::make_tuple(std::get<0>(node), std::get<1>(node), std::get<3>(node))] = 0;
    }
    for (auto & node : nodes) {
        node_counts[std::make_tuple(std::get<0>(node), std::get<1>(node), std::get<3>(node))] ++;
    }
    // make a std::vector of the net flat nodes
    std::vector<std::tuple<int, int, bool>> flat_nodes;
    for (auto & node_count : node_counts) {
        if (node_count.second % 2 == 1) {
            flat_nodes.push_back(node_count.first);
        }
    }
    // see what extras and bulk nodes are given
    std::set<int> given_extras;
    std::set<std::tuple<int, int, bool>> bulk_nodes;
    for (auto & node : flat_nodes) {
        if (std::get<2>(node) > 0) {
            given_extras.insert(std::get<0>(node));
        } else {
            bulk_nodes.insert(node);
        }
    }

    if (bulk_nodes.size()==0){
        // without bulk nodes, neutral only if no logical nodes are given (unless this is ignored)
        int neutral = (ignore_extras || given_extras.size()==0);
        int num_errors = 0;
        // compile the output
        output.push_back(neutral);
        output.push_back(num_errors);
        // no flipped extras need to be added
    } else {
        std::map<int, int> parities;
        // check how many times the bulk nodes turn up in each cycle
        for (int c = 0; c < cycle_dict.size(); c++){
            parities[c] = 0;
        }
        for (auto & node: bulk_nodes) {
            for (auto & c: cycle_dict[(std::make_tuple(std::get<0>(node), std::get<1>(node)))]){
                parities[c]++;
            }
        }
        // we get frustration if any of these is odd
        bool frust = false;
        for (auto & parity: parities){
            frust = frust || (parity.second % 2 == 1);
        }
        if (frust) {
            // if it's frustrated, it's not neutral
            output.push_back(0);
            // number of errors not counted
            output.push_back(-1);
            // no flipped extras need to be added
        } else {
            // now we must bicolor the qubits of the link graph, such that node edges connect unlike edges
            
            // first make a list of the qubits that definitely need to be covered
            // (those in the bulk nodes) and see how often each comes up in the nodes
            std::set<int> node_qubits;
            std::map<int, int> nq_nums;
            for (auto & node: bulk_nodes){
                std::vector<int> qs = {
                    std::get<0>(node), 
                    std::get<1>(node)
                };
                for (auto & q: qs) {
                    if (node_qubits.insert(q).second){
                        nq_nums[q] = 1;
                    } else {
                        nq_nums[q] += 1;
                    }
                }
            }
            // find the most mentioned qubit
            int root;
            int max_num = 0;
            for (auto & nq_num: nq_nums){
                if (nq_num.second > max_num){
                    root = nq_num.first;
                    max_num = nq_num.second;
                }
            }
            // start colouring with the most mentioned qubit
            std::map<int, bool> color;
            color[root] = 0;
            std::vector<int> newly_colored = {root};
            std::set<int> colored = {root};
            // stop once all node qubits are coloured and one color has stopped growing
            bool converged = false;
            node_qubits.erase(root);
            std::map<bool, int> num_nodes = {
                {0, 1},
                {1, 0}
            };
            std::map<bool, int> last_num_nodes = num_nodes;
            while (not converged){
                // for each newly coloured qubit
                std::vector<int> very_newly_colored;
                for (auto & n: newly_colored){
                    // loop through all the neighbours
                    for (auto & nn: link_neighbors[n]){
                        // if they haven't yet been coloured
                        if (colored.find(nn) == colored.end()){
                            // if this pair don't correspond to a bulk node, the new one is the same colour
                            if ((bulk_nodes.find(std::make_tuple(n,nn,false)) == bulk_nodes.end()) and (bulk_nodes.find(std::make_tuple(nn,n,false)) == bulk_nodes.end())){
                                color[nn] = color[n];
                            // otherwise, it's the opposite color
                            } else {
                                color[nn] = not color[n];
                            }
                            very_newly_colored.push_back(nn);
                            colored.insert(nn);
                            num_nodes[color[nn]]++;
                            node_qubits.erase(nn);
                        }
                    }
                }
                converged = (node_qubits.size() == 0) and ((num_nodes[0] == last_num_nodes[0]) or (num_nodes[1] == last_num_nodes[1]));
                newly_colored = very_newly_colored;
                if (not converged){
                    last_num_nodes = num_nodes;
                }
            }
            // see which colour has converged
            bool conv_color = (num_nodes[1] == last_num_nodes[1]);
            // calculate the number of nodes for the other
            num_nodes[not conv_color] = link_neighbors.size() - num_nodes[conv_color];
            // see which colour has the fewer qubits
            int min_color = (num_nodes[1] <= num_nodes[0]);
            // list the colours with the max error one first
            // (unless we do min only)
            std::vector<int> cs;
            cs.push_back(min_color);
            if (not minimal){
                cs.push_back((min_color+1)%2);
            }
            // determine which flipped extras correspond to which colour
            std::vector<std::set<int>> color_extras = {{}, {}};
            for (auto & qe: extras){
                if (color.find(qe.first) == color.end()){
                    color[qe.first] = (conv_color+1)%2;
                }
                color_extras[color[qe.first]].insert(qe.first);
            }
            // see if we can find a color for which we have no extra extras
            // and see what additional extras are required
            std::set<int> flipped_extras;
            std::set<int> flipped_ng_extras;
            std::vector<int> extra_extras;
            bool done = false;
            int j = 0;
            while (not done){
                flipped_extras = {};
                flipped_ng_extras = {};
                // see which extras for this colour have been flipped
                for (auto & q: color_extras[cs[j]]){
                    flipped_extras.insert(q);
                    // and which of those were not given
                    if (given_extras.find(q) == given_extras.end()) {
                        flipped_ng_extras.insert(q);
                    }
                }
                // see which extra extras are given
                extra_extras = {};
                if (not ignore_extras) {
                    for (auto & q: given_extras){
                        if ((flipped_extras.find(q) == flipped_extras.end())) {
                            extra_extras.push_back(q);
                        }
                    }
                }
                // if we have no extra extras, or we've run out of colours, we move on
                // otherwise we try the next colour
                done = (extra_extras.size()==0) || (j+1)==cs.size();
                if (not done){
                    j++;
                }
            }

            // construct output
            output.push_back(extra_extras.size()==0); // neutral
            output.push_back(num_nodes[cs[j]]); // num_errors
            for (auto & q: flipped_ng_extras){
                output.push_back(q);
            }

        }

    }

    return output;
};