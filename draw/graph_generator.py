"""
    Waxman model - graph generator for TOC paper

"""
import sys
import networkx as nx
import random
import numpy
import math
from math import ceil

RANGE = 300
G = nx.Graph()
vis = {}

def dfs(node):
    global vis
    global G
    if vis[node] == True:
        return
    vis[node] = True
    for v in G.neighbors(node):
        if vis[v] == False:
            dfs(v)
    return

def is_connect():
    global vis
    global G
    # print(type(G.nodes()))
    vis = {node : False for node in G.nodes()}
    for node in G.nodes():
        dfs(node)
        break
    for node in G.nodes():
        if(vis[node] == False):
            return False
    return True

def dist(p1, p2):
    (x1, y1) = p1
    (x2, y2) = p2
    return ((x1 - x2) ** 2 + (y1 - y2) ** 2) ** (1 / 2)


num_of_node = 100

while True:
    G = nx.waxman_graph(num_of_node, beta=0.85, alpha=0.03, domain=(0, 0, 0.5, 1))
    # G = nx.waxman_graph(num_of_node, beta=0.85, alpha=10, domain=(0, 0, 0.5, 1))
    positions = nx.get_node_attributes(G, 'pos')
    add_edge = []
    # connect the u's nearest vertex with u
    for u in range(G.order()-1):
        mi_dist = dist(positions[u], positions[G.order()-1])
        mi_idx = G.order()-1
        for v in range(u+1, G.order()):
            if(G.has_edge(u, v)):
                continue
            if(mi_dist > dist(positions[u], positions[v])):
                mi_dist = dist(positions[u], positions[v])
                mi_idx = v
        if(G.has_edge(u, mi_idx)):
            continue
        add_edge.append((u, mi_idx))

    for e in add_edge:
        (u, v) = e
        G.add_edge(u, v)

    if is_connect():
        break 
    else:
        print("topo is not connected", file=sys.stderr)

positions = nx.get_node_attributes(G, 'pos')
# write node
print(num_of_node)
for n in G.nodes():
    (x, y) = positions[n]
    pos_x = str(x*RANGE)
    pos_y = str(y*RANGE)
    num_of_memory = random.randint(-1, 1)
    print(x, y)

# write edge
num_of_edge = 0
for e in G.edges():
    if e[0] != e[1]:
        num_of_edge += 1
print(num_of_edge)
avg_l = 0
for e in G.edges():
    if e[0] != e[1]:
        print(e[0], e[1])