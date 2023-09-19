"""Temporary module with methods for graphs."""
import json
import os
import networkx as nx
import rustworkx as rx

from qiskit_qec.utils.decoding_graph_attributes import DecodingGraphEdge, DecodingGraphNode


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


def get_cached_decoding_graph(path):
    """
    Returns graph cached in file at path "file" using cache_graph method.
    """
    if os.path.isfile(path) and not os.stat(path) == 0:
        with open(path, "r+", encoding="utf-8") as file:
            json_data = json.loads(file.read())
            net_graph = nx.node_link_graph(json_data)
        ret_graph = rx.networkx_converter(net_graph, keep_attributes=True)
        for node_index, node in zip(ret_graph.node_indices(), ret_graph.nodes()):
            del node["__networkx_node__"]
            qubits = node.pop("qubits")
            time = node.pop("time")
            index = node.pop("index")
            is_boundary = node.pop("is_boundary")
            properties = node.copy()
            node = DecodingGraphNode(is_boundary=is_boundary, time=time, index=index, qubits=qubits)
            node.properties = properties
            ret_graph[node_index] = node
        for edge_index, edge in zip(ret_graph.edge_indices(), ret_graph.edges()):
            weight = edge.pop("weight")
            qubits = edge.pop("qubits")
            properties = edge.copy()
            edge = DecodingGraphEdge(weight=weight, qubits=qubits)
            edge.properties = properties
            ret_graph.update_edge_by_index(edge_index, edge)
        return ret_graph
    return None


def cache_decoding_graph(graph, path):
    """
    Cache rustworkx PyGraph to file at path.
    """
    net_graph = ret2net(graph)
    with open(path, "w+", encoding="utf-8") as file:
        json.dump(nx.node_link_data(net_graph), file)
