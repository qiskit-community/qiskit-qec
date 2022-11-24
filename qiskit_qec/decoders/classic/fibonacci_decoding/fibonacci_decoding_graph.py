import pymatching as pm
import rustworkx as rx


class FibonacciDecoderGraph:
    def __init__(
        self, matching_graph, hori_probe_fault_id, verti_probe_fault_id, stab2node
    ) -> None:

        self.matching_decoder = pm.Matching(matching_graph)
        self.matching_graph = matching_graph
        self.hori_probe_fault_id = hori_probe_fault_id
        self.verti_probe_fault_id = verti_probe_fault_id
        self.stab2node = stab2node

    def graph_label_attr_fn(self, fund_stab_labels=True):
        def node_attr(node):
            return {"label": str(node["element"])}

        def edge_attr(edge):
            return {"label ": str(list(edge["fault_ids"])[0])}

        return node_attr, edge_attr

    def decode_prob(self, syndrome):
        """Returns whether the horizontal probe edge (aka the bottom middle bit of the triangle) was in an even (0 aka no flip) or odd (1 aka flip) number of times in the matching graph && the vertical probe"""
        res = self.matching_decoder.decode(syndrome)
        return res[self.hori_probe_fault_id], res[self.verti_probe_fault_id], res
