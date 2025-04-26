#!/usr/bin/env python

import tskit, pyslim
import xml.dom.minidom # to fix tskit#3144

ts = tskit.load("out.trees")

focal_pos = 0
focal_site = ts.site(position=focal_pos)
tsk_to_slim = {-1 : ''}
slim_to_type = {'' : -1}
site_muts = focal_site.mutations
for m in site_muts:
    ml = m.metadata['mutation_list']
    assert len(ml) == 1
    slim_to_type[m.derived_state] = ml[0]['mutation_type']
    tsk_to_slim[m.id] = m.derived_state

type_labels = {1: "A", 2: "B", -1: "anc"}
type_colors = {"A" : "blue", "B": "red", "anc": "black"}
mut_labels = {m.id : type_labels[slim_to_type[tsk_to_slim[m.id]]] for m in site_muts}
nonmuts = [m.id for m in site_muts if m.parent >= 0 and mut_labels[m.id] == mut_labels[m.parent]]

t = ts.at(focal_pos)

styles = []
# get mutations in order by time so that descendant styling takes precedence
mids = list(set(m.id for m in site_muts) - set(nonmuts))
mids.sort(key=lambda m: -ts.mutation(m).time)
for mid in mids:
    col = type_colors[mut_labels[mid]]
    s = (
        f".m{mid} .node .edge, "
        f".mut.m{mid} line, "
        f".mut.m{mid} .sym "
        f"{{stroke: {col}; stroke-width: 2px}}"
        f".mut.m{mid} .lab "
        f"{{stroke: {col}}}"
        f".m{mid} .sym "
        f"{{stroke: {col}; fill: {col}}}"    )
    styles.append(s)
for mid in nonmuts + [m.id for m in ts.mutations() if m.site != focal_site.id]:
    s = f".mut.m{mid} .sym {{display: none}}"
    styles.append(s)
css_string = " ".join([f"#myUID {s}" for s in styles])

svg_size = (600, 500)
svg_string = t.draw_svg(path="tree.svg",
    size=svg_size,
    y_axis=True, y_label=" ", y_ticks=[1, 10, 100],
    x_axis=False,
    node_labels={},
    time_scale="log_time",
    mutation_labels={k: v for k,v in mut_labels.items() if k not in nonmuts},
    root_svg_attributes={'id': "myUID"},
    style=css_string,
)

