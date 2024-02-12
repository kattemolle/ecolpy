#!/usr/bin/env python3

import ast
from itertools import product
import json
import os
from time import time
import numpy as np
import networkx as nx
import z3


class LatticeGraph:
    """
    Defines a lattice graph.

    Parameters:
    -----------
    basis_graph : nx.Graph or list of edges.
        Graph with edges of the form ((0,0,i),(dx,dy,j)), with dx,dy the position of the other cell relative to the base unit cell. It can also be a coloring basis graph, with edges of the form ((x,y,i),(x+dx,y+dy,j)) with additionally an added "color" edge attribute that specifies the color as an int.

    patch size : tuple
        Patch size (n,m).

    name : str

    coloring_class : int
        Set to 1 if the lattice graph induced by the basis graph is to be class I colored, to 2 if it is to be class II colored, and to 3 if a mere proper edge coloring is sought (latter not tested). The int coloring_class is called $t$ in the paper.

    Attributes:
    -----------
    basis_graph
    patch_size
    name
    coloring_class

    Methods:
    --------
    get_max_x_y
    get_patch
    edge_color
    from_edge_colored_json
    from_json

    Example:
    --------
    import ecolpy as ep
    edges = [
        ((0, 0, 0), (1, -1, 0)),
        ((0, 0, 0), (1, 0, 0)),
        ((0, 0, 0), (1, 1, 0)),
        ((0, 0, 0), (0, 1, 0)),
    ]
    lat = ep.LatticeGraph(edges, name="diagonal")


    # Find a class 1  edge coloring of the lattice graph and store it as `clat` (an edge colored periodic lattice graph).
    clat = lat.edge_color(1)

    # Inspect the result.
    print(clat.basis_graph.edges(data=True))

    # Save plot of the lattice in the workding directory.
    ep.plot_edge_coloring(lat, clat, ".")
    """

    def __init__(self, basis_graph, patch_size=(1, 1), name=None, coloring_class=None):
        assert type(basis_graph) == nx.Graph or type(basis_graph) == list
        self.basis_graph = nx.Graph(basis_graph)
        assert (
            nx.number_of_selfloops(self.basis_graph) == 0
        )  # A nx.Graph is undirected and has no multiedges. If it has also no self loops it is simple.
        self.patch_size = patch_size
        self.name = name
        self.coloring_class = coloring_class

    def __str__(self):
        s = "Lattice graph with\nname: {}\npatch size: {}\ncoloring class: {}\nbasis graph: {}".format(
            self.name, self.patch_size, self.coloring_class, self.basis_graph.edges()
        )
        return s

    def get_max_x_y(self):
        """
        Return the maximum |x| and maximum |y| among the vertices (x,y,s) of self.basis_graph.
        """
        v = self.basis_graph.nodes()
        v = list(v)
        v.sort(key=lambda x: abs(x[0]))
        maxx = abs(v[-1][0])
        v.sort(key=lambda x: abs(x[1]))
        maxy = abs(v[-1][1])

        return maxx, maxy

    def get_patch(self, ax, bx, ay, by, boundaries="open"):
        """
        Return nx.Graph or nx.MultiGraph (which can have, but need not have, multiedges or self loops) of bx-ax by by-ay basis graphs or coloring basis graphs (which are themselves patches), with open or periodic boundary conditions. If boundaries=="open", no periodic boundary conditions are imposed. Then, returned graph is a nx.Graph.

        If boundaries=="periodic", perdiodic boundary conditions are applied and old edge data is stored as an edge attribute (wrapping). The returned graph is a nx.Multigraph.
        """

        def pbc(edge):
            """
            Given an edge as a tuple, impose periodic boundary conditions. Discards any edge attributes.
            """
            a = edge[0]
            b = edge[1]
            wa = (
                a[0] % (bx * self.patch_size[0]) + ax,
                a[1] % (by * self.patch_size[1]) + ay,
                a[2],
            )
            wb = (
                b[0] % (bx * self.patch_size[0]) + ax,
                b[1] % (by * self.patch_size[1]) + ay,
                b[2],
            )
            return (wa, wb)

        def wrap(graph):
            ng = nx.MultiGraph()
            for edge in graph.edges(data=True):
                if edge == pbc(edge):
                    ng.add_edges_from([edge])
                else:
                    new_edge = pbc(edge)
                    new_data = edge[2].copy()  # inherit attributes from old edge
                    new_data["old_edge"] = edge
                    new_edge = (*new_edge, new_data)
                    ng.add_edges_from([new_edge])
            return ng

        def lattice_cell(x, y):  # Return basis graph translated by (x,y).
            mapping = {
                node: (
                    node[0] + x * (self.patch_size[0]),
                    node[1] + y * (self.patch_size[1]),
                    node[2],
                )
                for node in self.basis_graph.nodes()
            }
            cell = nx.relabel_nodes(self.basis_graph, mapping)

            return cell

        g = nx.Graph()
        for x in range(ax, bx):
            for y in range(ay, by):
                cell = lattice_cell(x, y)
                g = nx.compose(g, cell)

        if boundaries == "periodic":
            g = wrap(g)

        return g

    def edge_color(self, cl):
        """
        Find class I coloring of lattice lat if cl==1 (and such a coloring exist). Find class 2 coloring if cl==2. Return edge colored LatticeGraph with LatticeGraph.coloring_class set to cl. Return None if no class cl edge coloring is found. The edge colored LatticeGraph inherits the name attribute from self.
        """
        maxn = 4
        maxm = 4
        pss = list(
            product(range(1, maxn), range(1, maxm))
        )  # Try various patch sizes. This can be replaced by the sequence of nondecreasing patch sizes ((1,1),(1,2),(2,1),(1,3),(3,1),(2,2),...).
        g = None
        i = 0
        while g == None and i < len(
            pss
        ):  # Edge colors super untit cells untill a class cl coloring is found or untill we exhaust the list of candidate coloring patchs.
            if (
                self.name == "t5145" or self.name == "t5282" or self.name == "t6369"
            ):  # This is a hack. For an unknown reason, z3 hangs for the these graphs from the galebach collection when n,m=1.
                n, m = 1, 2
            elif self.name == "NNN-square":
                n, m = 4, 4
            else:
                n, m = pss[i]
            print("      (n,m) =", (n, m))
            cng = self.get_patch(0, n, 0, m, boundaries="periodic")  # Candidate g.
            if (
                nx.number_of_selfloops(cng) == 0
            ):  # A nececerry and sufficient condition for a proper edge coloring to be extendable to the entire lattice graph and remain proper.
                g = edge_color_graph(cng, cl)
            i += 1

        if g == None:
            return None

        g = unwrap(g)
        clat = LatticeGraph(g, patch_size=(n, m), name=self.name, coloring_class=cl)
        return clat

    @classmethod
    def from_edge_colored_json(self, fname, clatname):
        """
        Make self the edge colored lattice object from edge colored json.
        """
        with open(fname, "r") as f:
            cdb = json.load(f)
        clatbg = nx.node_link_graph(cdb[clatname]["clatbg"])
        clatbg = string_nodes_to_arrays(clatbg)
        patch_size = cdb[clatname]["patch_size"]
        coloring_class = cdb[clatname]["coloring_class"]
        return self(clatbg, patch_size, clatname, coloring_class)

    @classmethod
    def from_json(self, fname, latname):
        """
        Make self the LatticeGraph object from the uncolored json at fname.
        """
        with open(fname, "r") as f:
            db = json.load(f)

        if (
            "T1" in db[latname] and "T2" in db[latname]
        ):  # If the db entry is a tiling symbol.
            edges = tiling_symbol_to_basis_graph(db[latname])
        else:  # otherwise, assume it is a graph.
            edges = nx.node_link_graph(db[latname])
            edges = string_nodes_to_arrays(edges)

        return self(edges, name=latname)


def edge_color_graph(g, coloring_class):
    """
    Return a class-`coloring_class` edge coloring of a nx.Graph or nx.MultiGraph g. If no such coloring exists, return None. If not None, for each edge of g, the 'color' attribute is set to the color (an int). The returned graph is always a nx.Multigraph. In the paper, coloring_class is called t.

    Parameters
    ----------
    g : nx.Graph or nx.MultiGraph or list of edges.

    coloring_class : int

    Returns
    -------
    cg : nx.MultiGraph or None
        Colored graph. Colors are stored in the "color" multiedge attribute.
    """

    def model_to_colored_graph(g, m):  # z3 model on graph g to colored graph cg
        cg = nx.MultiGraph()
        for e in g.edges(data=True, keys=True):
            ne = list(e)
            color = m[e[3]["color"]].as_long()
            ne[3]["color"] = color
            cg.add_edges_from([ne])
        return cg

    g = nx.MultiGraph(g)
    deg = max_degree(g)
    edges = list(g.edges(data=True))
    cs = z3.IntVector("c", len(edges))

    ng = nx.MultiGraph()
    for i, edge in enumerate(edges):
        new_edge = list(edge)
        new_edge[2]["color"] = cs[i]
        ng.add_edges_from([new_edge])
    g = ng

    s = z3.Solver()

    # Two colors cannot touch
    edges = list(g.edges(data=True, keys=True))
    for i in range(g.number_of_edges()):
        for j in range(i, g.number_of_edges()):
            e = edges[i]
            ep = edges[j]
            if e != ep:
                s.add(
                    z3.Implies(
                        e[3]["color"] == ep[3]["color"],
                        e[0] != ep[0]
                        and e[0] != ep[1]
                        and e[1] != ep[0]
                        and e[1] != ep[1],
                    )
                )

    cg = None
    if coloring_class == 1:
        s.add([z3.And(0 <= c, c < deg) for c in cs])
        if s.check() == z3.sat:
            m = s.model()
            cg = model_to_colored_graph(g, m)
    elif coloring_class == 2:
        s.add([z3.And(0 <= c, c < deg + 1) for c in cs])
        if s.check() == z3.sat:
            m = s.model()
            cg = model_to_colored_graph(g, m)

    return cg


def unwrap(g):
    """
    Given nx.MultiGraph g created with periodic boundary conditions, undo any wrapping of the edges, retaining any coloring of these edges. It is asserted that the unwrapped graph has no multiple edges and no self loops, and the graph is converted to a normal nx.Graph. Return the unwrapped nx.Graph.
    """
    assert type(g) == nx.MultiGraph, "Wrapped patch must be a networkx.MultiGraph"

    def has_multiple_edges(ug):
        for edge in ug.edges(keys=True):
            if edge[2] == 1:
                return True
        return False

    ug = nx.MultiGraph()
    for e in g.edges(data=True):
        if not "old_edge" in e[2]:  # If the edge was not wrapped
            ug.add_edges_from([e])
        else:
            oe = e[2]["old_edge"]
            oe[2]["color"] = e[2]["color"]
            ug.add_edges_from([oe])

    assert nx.number_of_selfloops(ug) == 0
    assert not has_multiple_edges(ug)
    ug = nx.Graph(ug)
    return ug


def max_degree(g):
    """
    Return max degree of the nx.Graph or nx.Multigraph g.
    """
    assert nx.number_of_selfloops(g) == 0
    degs = list(g.degree())
    degs.sort(key=lambda x: x[1])
    maxdeg = degs[-1][1]
    return maxdeg


def update_json(path, lat):
    """
    Update the json at path with the the LatticeGraph lat.

    Example:
    --------
    Create the diagonal lattice, store it to some new json file or update existing json.
    ```
    import networkx as nx
    import ecolpy as ep

    intra = []
    inter = [
        ((0, 0, 0), (1, -1, 0)),
        ((0, 0, 0), (1, 0, 0)),
        ((0, 0, 0), (1, 1, 0)),
        ((0, 0, 0), (0, 1, 0)),
    ]
    edges = intra + inter
    basis_graph = nx.Graph(edges)
    lat = ep.LatticeGraph(basis_graph, name="diagonal")
    ep.update_json("colorings/misc/misc.json", lat)
    ```
    """

    def nodes_to_string(g):
        """
        Returns graph with all node labels mapped to strings.
        Needed when reloading a graphs from a JSON. JSON edxport casts tuples to arrays and arrays cannot be used as nodes.
        """
        mapping = {node: str(node) for node in g}
        g = nx.relabel_nodes(g, mapping)
        return g

    if os.path.isfile(path):
        with open(path, "r") as f:
            db = json.load(f)
    else:
        db = {}

    assert lat.name is not None, "Specify a lattice name (as name attribute)"
    lat.basis_graph = nodes_to_string(lat.basis_graph)
    db[lat.name] = nx.node_link_data(lat.basis_graph)

    with open(path, "w") as f:
        json.dump(db, f)


def verify_coloring(clat):
    """
    Assert that the edge coloring of the LatticeGraph clat is proper and uses a number of colors consistent with the coloring class clat.cl. It is assumed 'color' attibutes are defined for the edges of the coloring basis graph of clat and that these colors are integers.
    """
    assert (
        clat.coloring_class != None
    ), "A coloring class must be provided as clat.coloring_class."
    cl = clat.coloring_class

    def is_proper(g):
        for n in g:
            colors = [g[n][nbr]["color"] for nbr in g[n]]
            if len(colors) != len(set(colors)):
                return False
        return True

    maxx, maxy = clat.get_max_x_y()
    n = 2 * maxx + 1
    m = 2 * maxy + 1
    g = clat.get_patch(0, n, 0, m)
    assert type(g) == nx.Graph
    colors = [edge[2]["color"] for edge in g.edges(data=True)]
    nc = len(set(colors))
    if cl == 1:
        return nc == max_degree(g) and is_proper(g)
    if cl == 2:
        return nc == max_degree(g) + 1 and is_proper(g)
    if cl == 3:
        return is_proper(g)


def tiling_symbol_to_basis_graph(ts):
    """
    Input a definition using a tiling symbol `ts`, a dict of e.g. the form
    {
    "T1":[2,0,-1,0],
    "T2":[-1,0,2,0],
    "Seed":[
    [0,0,0,0],
    [0,0,1,0],
    ],
    }

    See the folowing resources for how this defines a lattice.
    https://chequesoto.info/thesis.html,
    https://observablehq.com/@esperanc/synthesizing-periodic-tilings-of-the-plane,
    https://w3.impa.br/~cheque/tiling/,
    https://ieeexplore.ieee.org/document/8614306


    Return a list of edges of the form [(x,y,i),..] as used to define the basis graph in `LatticeGraph`.
    """
    t1 = np.array(ts["T1"])
    t2 = np.array(ts["T2"])
    seeds = np.array(ts["Seed"])
    dis = [
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1],
        [-1, 0, 1, 0],
        [0, -1, 0, 1],
        [-1, 0, 0, 0],
        [0, -1, 0, 0],
        [0, 0, -1, 0],
        [0, 0, 0, -1],
        [1, 0, -1, 0],
        [0, 1, 0, -1],
    ]  # Directions

    # Create a patch of 3x3 basis graphs.
    # Keep the x,y,seed data as entries of xyseed. That is, xyseed[i] is the (x,y,seed) of the point patch[i].

    patch = []
    xyseed = []
    for x in [0, -1, 1]:
        for y in [0, -1, 1]:
            for seed in seeds:
                point = x * t1 + y * t2 + seed
                point = tuple(point)
                patch += [point]
                xyseed += [(x, y, tuple(seed.tolist()))]

    # For every seed, check the neighbors and add edges using the old xyseed data.

    def no_equivalent_edge(edge, edges):
        # Translate the edge to the surroundings. If already present, return false.
        for x, y in product([0, -1, 1], [0, -1, 1]):
            an1 = (edge[0][0] + x, edge[0][1] + y, edge[0][2])
            an2 = (edge[1][0] + x, edge[1][1] + y, edge[1][2])
            ae = (an1, an2)
            if ae in edges or ae[::-1] in edges:
                return False
        return True

    edges = []
    for seed in seeds:
        for di in dis:
            ca = tuple(seed + np.array(di))  # candidate point
            if ca in patch:
                ind = patch.index(ca)
                edge = ((0, 0, tuple(seed.tolist())), xyseed[ind])
                if no_equivalent_edge(edge, edges):
                    edges += [edge]

    return edges


def plot_edge_coloring(lat, clat, folder, symbol=None):
    """
    Use Graphviz to make an image of the colored lattice graph clat. Edges of clat.basis_graph should have attribute "color" which is an int. Saves the plot under folder/<clat.name>.pdf.

    Parameters:
    -----------
    lat : LatticeGraph
        Lattice graph with primal basis graph. Need not be colored. Any colors are ignored.

    clat : LatticeGraph
        Lattice graph with a coloring basis graph.

    folder : str
        Output is saved under folder/<clat.name>.pdf.

    symbol: None or dict
        If a lattice symbol is given, use the t1 and t2 vectors in this symbol (in seed notation) to reconstruct the geometric location of all points to be plotted.
    """
    pyplotc = [
        "#1f77b4",
        "#ff7f0e",
        "#2ca02c",
        "#d62728",
        "#9467bd",
        "#8c564b",
        "#e377c2",
        "#7f7f7f",
        "#bcbd22",
        "#17becf",
    ]
    otherc = ["#2f17a9", "#99f14a"]
    colors = pyplotc + otherc

    def xyseed_to_cartesian(xyseed, t1, t2):
        om = complex(np.sqrt(3) / 2 + 1j / 2)

        def seed_to_compl(seed):
            a = [seed[i] * om**i for i in range(4)]
            return sum(a)

        x, y, seed = xyseed
        ct1 = seed_to_compl(t1)
        ct2 = seed_to_compl(t2)
        cseed = seed_to_compl(seed)
        cpos = x * ct1 + y * ct2 + cseed
        return cpos.real, cpos.imag

    def cartesian_to_string(c):
        spos = str(c[0]) + "," + str(c[1])
        return spos

    print("    Creating plot at", folder + "/" + clat.name + ".pdf")
    penwidth = 5
    wide_penwidth = 16
    width = 0.1
    headclip = "false"
    tailclip = "false"
    dashed_style = "dashed"

    lat.basis_graph = string_nodes_to_arrays(lat.basis_graph)
    clat.basis_graph = string_nodes_to_arrays(clat.basis_graph)
    g = clat.get_patch(-1, 2, -1, 2, boundaries="open")
    nx.set_node_attributes(g, "point", "shape")
    nx.set_node_attributes(g, width, "width")
    nx.set_edge_attributes(g, headclip, "headclip")
    nx.set_edge_attributes(g, tailclip, "tailclip")
    nx.set_edge_attributes(g, penwidth, "penwidth")
    args = "-Goutputorder=edgesfirst"

    if clat.name == "NNN-square":
        penwidth = 5
        wide_penwidth = 5
        width = 0.1
        g = clat.get_patch(0, 1, 0, 1, boundaries="open")
        nx.set_node_attributes(g, "point", "shape")
        nx.set_node_attributes(g, width, "width")
        dashed_style = "solid"
        args = "-Goutputorder=nodesfirst -Gsplines=true"

    for edge in g.edges(data=True):
        color = colors[edge[2]["color"]]
        g[edge[0]][edge[1]]["color"] = color
        if clat.basis_graph.has_edge(edge[0], edge[1]):
            g[edge[0]][edge[1]]["penwidth"] = wide_penwidth
            g[edge[0]][edge[1]]["style"] = dashed_style

        if lat.basis_graph.has_edge(edge[0], edge[1]):
            g[edge[0]][edge[1]]["style"] = "solid"  # Overwrites earlier dashed setting.

    if symbol == None:
        assert all(
            type(node[0]) == type(node[1]) == int for node in clat.basis_graph
        ), "First two entries of a node must be ints, representing the x,y coordinate of the cell the node is in."
        mapping = {node: {"pos": cartesian_to_string(node[:2])} for node in g}
        if clat.name in ["square", "diagonal", "NNN-square"]:
            nx.set_node_attributes(g, "true", "pin")
        else:
            nx.set_node_attributes(g, "false", "pin")

    else:
        for node in g.nodes():
            assert len(node) == 3
            assert len(node[2]) == 4
        t1 = symbol["T1"]
        t2 = symbol["T2"]
        mapping = {
            xyseed: {"pos": cartesian_to_string(xyseed_to_cartesian(xyseed, t1, t2))}
            for xyseed in g
        }

    nx.set_node_attributes(g, mapping)
    largest_cc = max(nx.connected_components(g), key=len)
    if clat.name != "NNN-square":
        g = g.subgraph(largest_cc)
    g = nx.nx_agraph.to_agraph(g)
    g.draw(
        folder + "/" + clat.name + ".pdf",
        prog="neato",
        args=args,
    )


def string_nodes_to_arrays(g):
    """
    If the nodes were arrays converted to strings (because of JSON export), unconvert to get arrays again. Note that *every* node that is of string type will be interpreted as a literal python expression.
    """
    mapping = {node: ast.literal_eval(node) for node in g if type(node) == str}
    g = nx.relabel_nodes(g, mapping)
    return g


def edge_color_json(db, cl, f=lambda x: True):
    """
    Find and store a class `cl` edge coloring of every lattice graph defined in the json database db, applying filter f. Filter f is a function that takes a lattice name (str) and returns a Bool. It can be used to edge color parts of the database.

    Returns a json containing all coloring basis graphs. Results are keyed by the lattice name in the original json and contain the basis graph under "latbg", the coloring basis graph under "clatbg", and the patch size of the coloring basis graph under "patch_size" and the coloring class under "coloing_class". The graph data is in the form of the output of `networkx.node_link_data`.
    """

    def nodes_to_string(g):
        """
        Returns graph with all node labels mapped to strings.
        Needed when reloading a graphs from a JSON. JSON export casts tuples to arrays and arrays cannot be used as nodes.
        """
        mapping = {node: str(node) for node in g}
        g = nx.relabel_nodes(g, mapping)
        return g

    cdb = {}
    for name in db:
        if f(name) == True:
            print("    Coloring", name, flush=True)
            if (
                "T1" in db[name] and "T2" in db[name]
            ):  # If the db entry is a tiling symbol
                edges = tiling_symbol_to_basis_graph(db[name])
            else:  # otherwise, assume it is a graph.
                edges = nx.node_link_graph(db[name])
                edges = string_nodes_to_arrays(edges)

            lat = LatticeGraph(edges, name=name)
            start = time()
            clat = lat.edge_color(cl)
            assert clat != None, (
                "No class "
                + cl
                + " coloring found! No such coloring exists or maxn or maxm was set too low."
            )
            end = time()
            print("      Found coloring with patch size", clat.patch_size)

            lat.basis_graph = nodes_to_string(lat.basis_graph)
            clat.basis_graph = nodes_to_string(clat.basis_graph)

            cdb[name] = {}
            cdb[name]["latbg"] = nx.node_link_data(lat.basis_graph)
            cdb[name]["clatbg"] = nx.node_link_data(clat.basis_graph)
            cdb[name]["patch_size"] = clat.patch_size
            cdb[name]["coloring_class"] = clat.coloring_class
            cdb[name]["wall_clock"] = end - start

    return cdb


def plot_edge_colored_json(cdb, folder, db=None):
    """
    Plot all graphs from a json dict `cdb` containing edge coloring basis graphs. Put the plots in folder/name.pdf, where name is derived from clat.name.
    If a json dict `db` is given, use this dict to get symbols for geometric plotting.
    """
    for name in cdb:
        latbg = nx.node_link_graph(cdb[name]["latbg"])
        lat = LatticeGraph(latbg, name)
        clatbg = nx.node_link_graph(cdb[name]["clatbg"])
        patch_size = cdb[name]["patch_size"]
        clat = LatticeGraph(clatbg, patch_size, name)
        if db != None:
            symbol = db[name]
        else:
            symbol = None
        plot_edge_coloring(lat, clat, folder, symbol=symbol)
