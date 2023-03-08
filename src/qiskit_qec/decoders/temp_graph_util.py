"""Temporary module with methods for graphs."""
import json
import networkx as nx
import rustworkx as rx
import os


def ret2net(graph: rx.PyGraph):
    """Convert rustworkx graph to equivalent networkx graph."""
    nx_graph = nx.Graph()
    for j, node in enumerate(graph.nodes()):
        nx_graph.add_node(j)
        for k, v in node:
            if isinstance(v, dict):
                for _k, _v in v.items():
                    nx.set_node_attributes(nx_graph, {j: _v}, _k)
            else:
                nx.set_node_attributes(nx_graph, {j: v}, k)
    for j, (n0, n1) in enumerate(graph.edge_list()):
        nx_graph.add_edge(n0, n1)
        for k, v in graph.edges()[j]:
            if isinstance(v, dict):
                for _k, _v in v.items():
                    nx.set_edge_attributes(nx_graph, {(n0, n1): _v}, _k)
            else:
                nx.set_edge_attributes(nx_graph, {(n0, n1): v}, k)
    return nx_graph


def write_graph_to_json(graph: rx.PyGraph, filename: str):
    """Export a JSON formatted file with the graph data."""
    with open(filename, "w", encoding="utf-8") as fp:
        # Any fields that are not serializable are converted to strings
        from_ret = ret2net(graph)
        json.dump(nx.node_link_data(from_ret), fp, indent=4, default=str)
    fp.close()


def get_cached_graph(file):
    if os.path.isfile(file) and not os.stat(file) == 0:
        with open(file, "r+") as f:
            json_data = json.loads(f.read())
            net_graph = nx.node_link_graph(json_data)
        ret_graph = rx.networkx_converter(net_graph, keep_attributes=True)
        for node in ret_graph.nodes():
            del node["__networkx_node__"]
        return ret_graph
    return None


def cache_graph(graph, file):
    net_graph = ret2net(graph)
    with open(file, "w+") as f:
        json.dump(nx.node_link_data(net_graph), f)
