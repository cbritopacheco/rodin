import os
import re
import subprocess
import numpy as np
import argparse
import networkx as nx
from pathlib import Path
from shutil import which
from termcolor import colored
from networkx.drawing.nx_agraph import write_dot

PT_TO_INCH = (1.0 / 72)

MIN_NODE_SIZE = 2.0
MAX_NODE_SIZE = 7.0

MIN_FONT_SIZE = 6
MAX_FONT_SIZE = 42


COMMON_CPP_EXTENSIONS = ['.h', '.hpp', '.hxx']

COLORSCHEME = 'spectral9'
NUMBER_OF_COLORS = 9

RODIN_TOP_LEVEL_INCLUDES = {
    'Rodin.h',
    'Rodin/Alert.h',
    'Rodin/Variational.h',
    'Rodin/Geometry.h',
    'RodinExternal.h',
    'RodinExternal/MMG.h',
}

RODIN_TOP_LEVEL_INCLUDES_COLOR_MAP = {
    'Alert.h': 'red',
    'Variational.h': 'green',
    'Geometry.h': 'blue'
}


class GraphvizNotFoundError(Exception):
    pass


def merge_matching_paths(path1, path2):
    parts1 = Path(path1).parts
    parts2 = Path(path2).parts
    if parts2[:len(parts1)] == parts1:
        merged_parts = parts1 + parts2[len(parts1):]
        merged_path = Path(*merged_parts)
        return str(merged_path)
    else:
        joined_path = Path(path1) / path2
        return str(joined_path)


def extract_includes(file_path, context_path):
    includes = []
    with open(file_path, 'r') as file:
        for line in file:
            # system_include_match = re.match(r'#include\s*[<](.*)[>]', line)
            user_include_match = re.match(r'#include\s*["](.*)["]', line)
            if user_include_match:
                include = user_include_match.group(1)
                if include == 'ForwardDecls.h':
                    continue
                include_path = os.path.dirname(include)
                if include_path:
                    if not include_path.startswith('Rodin'):
                        include = merge_matching_paths(context_path, include)
                else:
                    include = os.path.join(context_path, include)
                includes.append(include)
    return includes


def generate_dependency_graph(root_path, extensions):
    G = nx.DiGraph()
    for root, _, files in os.walk(root_path):
        for file in files:
            if any(file.endswith(ext) for ext in extensions):
                source_path = os.path.join(root, file)
                relative_path = os.path.relpath(source_path, root_path)
                context_path = os.path.dirname(relative_path)
                includes = extract_includes(source_path, context_path)
                G.add_node(relative_path)
                for include in includes:
                    G.add_edge(relative_path, include)
    return G

def get_predecessors_map(graph):
    sinks = set()
    node_to_set = {}
    for node in graph.nodes():
        node_to_set[node] = set()
        if graph.out_degree(node) == 0:
            sinks.add(node)
    for sink in sinks:
        visited = set()
        q = { sink }
        while len(q) > 0:
            node = q.pop()
            visited.add(node)
            for n in graph.predecessors(node):
                node_to_set[n].add(node)
                node_to_set[n] |= node_to_set[node]
                if n not in visited:
                    q.add(n)
    return node_to_set

def get_accumulated_map(graph):
    predecessors_map = get_predecessors_map(graph)
    node_to_accumulated = {}
    for node, predecessors in predecessors_map.items():
        node_to_accumulated[node] = len(predecessors)
    return node_to_accumulated

def get_flow_map(graph):
    node_to_flow = {}
    for node in graph.nodes():
        out_degree = graph.out_degree(node)
        in_degree = graph.in_degree(node)
        node_to_flow[node] = in_degree - out_degree
    return node_to_flow

def get_laplacian_map(graph):
    laplacian_matrix = nx.directed_laplacian_matrix(graph)
    node_to_laplacian = {}
    for node_idx, node in enumerate(graph.nodes()):
        laplacian_value = laplacian_matrix[node_idx, node_idx]
        node_to_laplacian[node] = laplacian_value
    return node_to_laplacian

def get_pagerank_map(graph):
    pagerank_scores = nx.pagerank(graph, alpha=0.85)
    return pagerank_scores

def get_color_map(m):
    max_v = max(m.values())
    min_v = min(m.values())
    color_map = {}
    for node, v in m.items():
        std = float(v - min_v) / float(max_v - min_v)
        color_map[node] = int(NUMBER_OF_COLORS - std * (NUMBER_OF_COLORS - 1))
    return color_map

def get_nodesize_map(m):
    max_v = max(m.values())
    min_v = min(m.values())
    size_map = {}
    for node, v in m.items():
        std = float(v - min_v) / float(max_v - min_v)
        size_map[node] = int((MAX_NODE_SIZE - MIN_NODE_SIZE) * std + MIN_NODE_SIZE)
    return size_map

def get_fontsize_map(m):
    max_v = max(m.values())
    min_v = min(m.values())
    fontsize_map = {}
    for node, v in m.items():
        std = float(v - min_v) / float(max_v - min_v)
        fontsize_map[node] = int((MAX_FONT_SIZE - MIN_FONT_SIZE) * std + MIN_FONT_SIZE)
    return fontsize_map


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Generates #include dependency graph.')
    parser.add_argument('root_path', help='Root path of the codebase.')
    parser.add_argument('--layout_engine',
                        default='dot',
                        help='Layout engine.')
    parser.add_argument('--dot',
                        help='Output DOT file path.')
    parser.add_argument('--png',
                        help='Output PNG file path.')
    parser.add_argument('--svg',
                        help='Output SVG file path.')
    parser.add_argument('--extensions', nargs='+',
                        default=COMMON_CPP_EXTENSIONS,
                        help='List of file extensions to consider.')
    args = parser.parse_args()

    root_path = args.root_path
    folder_name = os.path.basename(os.path.normpath(root_path))

    dot_filename = args.dot if args.dot else f'{folder_name}_include_dependency_graph.dot'

    root_path = args.root_path
    dependency_graph = generate_dependency_graph(root_path, args.extensions)
    weight_map = get_accumulated_map(dependency_graph)
    color_map = get_color_map(weight_map)
    fontsize_map = get_fontsize_map(weight_map)
    nodesize_map = get_nodesize_map(weight_map)

    nx.set_node_attributes(dependency_graph, 'filled', name='style')
    nx.set_node_attributes(dependency_graph, 'circle', name='shape')
    nx.set_node_attributes(dependency_graph, 'true', name='fixedsize')
    nx.set_node_attributes(dependency_graph, COLORSCHEME, name='colorscheme')
    nx.set_node_attributes(dependency_graph, color_map, name='fillcolor')
    nx.set_node_attributes(dependency_graph, weight_map, name='weight')
    nx.set_node_attributes(dependency_graph, nodesize_map, name='width')
    nx.set_node_attributes(dependency_graph, nodesize_map, name='height')
    nx.set_node_attributes(dependency_graph, fontsize_map, name='fontsize')

    write_dot(dependency_graph, dot_filename)
    print(colored(f'DOT output generated: {dot_filename}', 'green'))

    layout_engine = args.layout_engine
    print(colored(f'Using layout engine {layout_engine}', 'blue'))
    dot_command = which(layout_engine)
    if not dot_command:
        raise GraphvizNotFoundError(
            "Graphviz %s command not found. Please install Graphviz to generate graphical output." % command)

    command_line = [f'{dot_command}', dot_filename]
    command_line += ['-Goverlap=prism', '-Gbeautify=true', '-Gsmoothing=true ']
    if args.png:
        png_filename = args.png
        command_line += ["-Tpng -o %s " % png_filename]

    if args.svg:
        svg_filename = args.svg
        command_line += ["-Tsvg -o %s " % svg_filename]


    result = subprocess.run(command_line, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    if result.returncode == 0:
        if args.png:
            print(colored(f'PNG output generated: {png_filename}', 'green'))
        if args.svg:
            print(colored(f'SVG output generated: {svg_filename}', 'green'))
    else:
        raise result.stderr

