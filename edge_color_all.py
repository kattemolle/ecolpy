#!/usr/bin/env python3
"""
This script edge colors, plots and verifies all edge colored lattices contained in the galebach, saesa, laves, misc. and classII databases.

For the production run, the timing data was
Timing data:
  Timing data for misc
    4 edge colored lattice graphs found in database.
    Total wall clock time 25.614869832992554 s
    Mean wall clock time 6.403717458248138 s
    Maximum wall clock time occurs for NNN-square and is 25.52237105369568 s
  Timing data for classII
    1 edge colored lattice graphs found in database.
    Total wall clock time 0.025203943252563477 s
    Mean wall clock time 0.025203943252563477 s
    Maximum wall clock time occurs for wheel_decorated_honeycomb and is 0.025203943252563477 s
  Timing data for laves
    8 edge colored lattice graphs found in database.
    Total wall clock time 0.3970973491668701 s
    Mean wall clock time 0.049637168645858765 s
    Maximum wall clock time occurs for dualstar and is 0.16474699974060059 s
  Timing data for clean_saesa
    57 edge colored lattice graphs found in database.
    Total wall clock time 24.067044973373413 s
    Mean wall clock time 0.4222288591819897 s
    Maximum wall clock time occurs for NQTUVW and is 2.4447808265686035 s
  Timing data for galebach
    1248 edge colored lattice graphs found in database.
    Total wall clock time 608.2434813976288 s
    Mean wall clock time 0.4873745844532282 s
    Maximum wall clock time occurs for t6665 and is 74.43041372299194 s
  Timing data for Archimedean lattices (subset of galebach)
    11 Archimedean lattice graphs found in galebach.
    Total wall clock time 0.3268449306488037 s
    Mean wall clock time 0.02971317551352761 s
    Maximum wall clock time occurs for t1007 and is 0.05106782913208008 s
"""
from importlib import reload
import json
from time import time
import networkx as nx
import ecolpy as ep
import os


def edge_color_databases(dbnames):
    """
    Minimal edge color the graphs in the json databases with a name in dbnames. E.g. dbnames=["galebach", "clean_saesa","laves","misc","classII"] if you want to edge-color all databases.
    """
    for name in dbnames:
        fname = "colorings/{name}/{name}.json".format(name=name)
        with open(fname, "r") as f:
            db = json.load(f)

        print("  Edge coloring", fname)
        cl = 1
        if (
            name == "classII"
        ):  # The lattices in the classII database get a special treatment: these lattices are known to be class II and a type II coloring is sought for, which is hence minimal.
            cl = 2
        if (
            name == "galebach"
        ):  # The galebach collection actually also contains, seperately, the x-uniform-archimedean lattices. Exclude those superflous entries.
            f = lambda x: "u" not in x
        else:
            f = lambda x: True
        start = time()
        cdb = ep.edge_color_json(db, cl, f)
        end = time()
        t = end - start
        print("  Edge-colored {} in {} s".format(name, t))

        ecfname = "colorings/{name}/edge_colored_{name}.json".format(name=name)
        with open(ecfname, "w") as f:
            json.dump(cdb, f)

        print("  Results stored in", ecfname)


def verify_edge_colored_databases(dbnames):
    """
    Verify that the edge coloring of all lattice graph in the databases in dbnames agrees with the coloring_class attribute of those lattices.
    """
    for name in dbnames:
        fname = "colorings/{name}/edge_colored_{name}.json".format(name=name)
        print("  Verifying edge colorings in", name)
        with open(fname, "r") as f:
            cdb = json.load(f)

        for clatname in cdb:
            print("    Verifying ", clatname)
            clatbg = ep.string_nodes_to_arrays(
                nx.node_link_graph(cdb[clatname]["clatbg"])
            )
            ps = cdb[clatname]["patch_size"]
            cl = cdb[clatname]["coloring_class"]
            clat = ep.LatticeGraph(clatbg, ps, clatname, cl)
            assert ep.verify_coloring(
                clat
            ), "Error in verification of the latice graph {} from database {}".format(
                clatname, name
            )
    print("  All edge colorings are correct.")


def plot_edge_colored_databases(dbnames):
    """
    Plot the edge colored graphs in the edge colored json databases with name in dbnames. E.g. dbnames=["misc", "laves", "clean_saesa", "galebach", "classII"] if you want to plot the edge-colored lattices in all databases.
    """
    for name in dbnames:
        fname = "colorings/{name}/edge_colored_{name}.json".format(name=name)
        with open(fname, "r") as f:
            cdb = json.load(f)

        print("  Plotting edge colorings from", fname)
        start = time()
        ep.plot_edge_colored_json(cdb, "colorings/{name}".format(name=name))
        end = time()
        t = end - start
        print("  Plotted edge-colored {name} in {t} s".format(name=name, t=t))


def timing_data(dbnames):
    """
    Print wall clock times used to compute the lattice graphs in databases. E.g. dbnames=["misc", "laves", "clean_saesa", "galebach", "classII"].
    """
    data = {}
    for name in dbnames:
        fname = "colorings/{name}/edge_colored_{name}.json".format(name=name)
        with open(fname, "r") as f:
            cdb = json.load(f)

        wcs = {}
        for clatname in cdb:
            wc = cdb[clatname]["wall_clock"]
            wcs[clatname] = wc

        data[name] = wcs

    for name in data:
        print("  Timing data for", name)
        ts = [data[name][clatname] for clatname in data[name]]
        print("   ", len(ts), "edge colored lattice graphs found in database.")
        tnwc = sum(ts)
        print("    Total wall clock time", tnwc, "s")
        mnwc = tnwc / len(data[name])
        print("    Mean wall clock time", mnwc, "s")
        srt = sorted(data[name].items(), key=lambda item: item[1])
        print(
            "    Maximum wall clock time occurs for",
            srt[-1][0],
            "and is",
            srt[-1][1],
            "s",
        )


def timing_data_archimedean():
    """
    The Archimedean lattices are a subset of the k uniform lattices so they need seperate timing data.
    """
    data = {}
    if os.path.isfile("colorings/galebach/edge_colored_galebach.json"):
        with open("colorings/galebach/edge_colored_galebach.json", "r") as f:
            cdb = json.load(f)

        wcs = {}
        for clatname in cdb:
            if clatname[1] == "1":
                wc = cdb[clatname]["wall_clock"]
                wcs[clatname] = wc

        print("  Timing data for Archimedean lattices (subset of galebach)")
        ts = [wcs[clatname] for clatname in wcs]
        print("   ", len(ts), "Archimedean lattice graphs found in galebach.")
        tnwc = sum(ts)
        print("    Total wall clock time", tnwc, "s")
        mnwc = tnwc / len(wcs)
        print("    Mean wall clock time", mnwc, "s")
        srt = sorted(wcs.items(), key=lambda item: item[1])
        print(
            "    Maximum wall clock time occurs for",
            srt[-1][0],
            "and is",
            srt[-1][1],
            "s",
        )


if __name__ == "__main__":
    dbnames = ["misc", "classII", "laves", "clean_saesa", "galebach"]
    print("Edge coloring databases", dbnames)
    edge_color_databases(dbnames)

    print("Plotting databases", dbnames)
    plot_edge_colored_databases(dbnames)

    print("Verifying edge colored databases", dbnames)
    verify_edge_colored_databases(dbnames)

    print("")
    print("Timing data:")
    timing_data(dbnames)
    timing_data_archimedean()

    print("Done")
