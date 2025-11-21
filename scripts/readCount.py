import gffutils
import pysam
import numpy as np
import plotly.graph_objs as go
from plotly.subplots import make_subplots
from plotly.colors import DEFAULT_PLOTLY_COLORS
import plotly.express as px
import pandas as pd
import re
from typing import Dict

#from GenBro.gff_helper import read_gff3, filter_gff_by_interval, find_gene_coordinates, coords_from_gff_row
#from GenBro.process_coverage import precompute_from_design, filter_strand

# Regex patterns to match different CSS color formats
hex_pattern = re.compile(r'^#([0-9a-fA-F]{3,8})$')
rgb_pattern = re.compile(r'^rgb\((\d{1,3}),\s*(\d{1,3}),\s*(\d{1,3})\)$')
rgba_pattern = re.compile(r'^rgba\((\d{1,3}),\s*(\d{1,3}),\s*(\d{1,3}),\s*(\d+(\.\d+)?)\)$')
hsl_pattern = re.compile(r'^hsl\((\d{1,3}),\s*(\d{1,3})%,\s*(\d{1,3})%\)$')
hsla_pattern = re.compile(r'^hsla\((\d{1,3}),\s*(\d{1,3})%,\s*(\d{1,3})%,\s*(\d+(\.\d+)?)\)$')


GFF3_TYPE_TO_COLOR = {
    "gene": px.colors.qualitative.Light24[0],
    "mRNA": px.colors.qualitative.Light24[1],
    "exon": px.colors.qualitative.Light24[2],
    "CDS": px.colors.qualitative.Light24[3],
    "start_codon": px.colors.qualitative.Light24[4],
    "stop_codon": px.colors.qualitative.Light24[5],
    "five_prime_UTR": px.colors.qualitative.Light24[6],
    "5UTR": px.colors.qualitative.Light24[6],
    "5'UTR": px.colors.qualitative.Light24[6],
    "three_prime_UTR": px.colors.qualitative.Light24[7],
    "3UTR": px.colors.qualitative.Light24[7],
    "3'UTR": px.colors.qualitative.Light24[7],
    "ncRNA": px.colors.qualitative.Light24[8],
    "rRNA": px.colors.qualitative.Light24[9],
    "tRNA": px.colors.qualitative.Light24[10],
    "repeat_region": px.colors.qualitative.Light24[11],
    "default": "#7CB9E8"
}
DO_NOT_ANNOTATE = [
   key for key in GFF3_TYPE_TO_COLOR.keys() if "UTR" in key or "codon" in key
]


def empty_figure(annotation: str = None):
    fig = make_subplots(rows=3, cols=1)
    fig.update_yaxes(showticklabels=False, showgrid=False)
    fig.update_xaxes(showgrid=False, showticklabels=False)
    fig.update_layout(
        margin={"t": 0, "b": 0, "r": 50},
        font=dict(
            size=16,
        ),
        yaxis=dict(zeroline=False),
        xaxis=dict(zeroline=False),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",

    )
    if annotation is not None:
        fig.add_annotation(
            xref="paper",
            yref="y2 domain",
            xanchor="center",
            yanchor="middle",
            x=0.5,
            y=0.5,
            text=annotation,
            showarrow=False,
            font=(dict(size=28))
        )
    fig.layout.template = "plotly_white"
    for row in range(1, 4):
        fig.add_trace(go.Scatter(

            x=[None],
            y=[None],
            showlegend=False,

        ),
            row=row,
            col=1
        )
    fig.update_layout(
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),

    )
    fig.update_xaxes(
        showline=True,
        mirror=True,

    )
    fig.update_yaxes(
        showline=True,
        mirror=True,
        row=2

    )
    return fig

def hex_to_rgb(hex_color):
    hex_color = hex_color.lstrip('#')
    length = len(hex_color)
    if length == 3:
        return tuple(int(hex_color[i] * 2, 16) for i in range(3))
    elif length == 4:
        return tuple(int(hex_color[i] * 2, 16) for i in range(3)) + (int(hex_color[3] * 2, 16) / 255,)
    elif length == 6:
        return tuple(int(hex_color[i:i + 2], 16) for i in (0, 2, 4))
    elif length == 8:
        return tuple(int(hex_color[i:i + 2], 16) for i in (0, 2, 4)) + (int(hex_color[6:8], 16) / 255,)


def hsl_to_rgb(h, s, l):
    s /= 100
    l /= 100
    c = (1 - abs(2 * l - 1)) * s
    x = c * (1 - abs((h / 60) % 2 - 1))
    m = l - c / 2
    if 0 <= h < 60:
        r, g, b = c, x, 0
    elif 60 <= h < 120:
        r, g, b = x, c, 0
    elif 120 <= h < 180:
        r, g, b = 0, c, x
    elif 180 <= h < 240:
        r, g, b = 0, x, c
    elif 240 <= h < 300:
        r, g, b = x, 0, c
    elif 300 <= h < 360:
        r, g, b = c, 0, x
    else:
        raise ValueError
    return (int((r + m) * 255), int((g + m) * 255), int((b + m) * 255))


def css_color_to_rgba(color, opacity):
    """
    Convert any acceptable CSS color to the same color with the specified opacity.

    :param color: A string representing the CSS color.
    :param opacity: A float representing the opacity value (between 0 and 1).
    :return: A string representing the color in RGBA format.
    """
    # Ensure opacity is a valid float between 0 and 1
    if not (0 <= opacity <= 1):
        raise ValueError("Opacity must be a float between 0 and 1.")


    match = hex_pattern.match(color)
    if match:
        rgba = hex_to_rgb(match.group(1))
        if len(rgba) == 4:
            return f'rgba({rgba[0]}, {rgba[1]}, {rgba[2]}, {opacity})'
        else:
            return f'rgba({rgba[0]}, {rgba[1]}, {rgba[2]}, {opacity})'

    match = rgb_pattern.match(color)
    if match:
        return f'rgba({match.group(1)}, {match.group(2)}, {match.group(3)}, {opacity})'

    match = rgba_pattern.match(color)
    if match:
        return f'rgba({match.group(1)}, {match.group(2)}, {match.group(3)}, {opacity})'

    match = hsl_pattern.match(color)
    if match:
        rgb = hsl_to_rgb(int(match.group(1)), int(match.group(2)), int(match.group(3)))
        return f'rgba({rgb[0]}, {rgb[1]}, {rgb[2]}, {opacity})'

    match = hsla_pattern.match(color)
    if match:
        rgb = hsl_to_rgb(int(match.group(1)), int(match.group(2)), int(match.group(3)))
        return f'rgba({rgb[0]}, {rgb[1]}, {rgb[2]}, {opacity})'

    raise ValueError("Invalid CSS color format.")


def get_trace(file, name, color, strand: str = None, size_factor: float = None,  **kwargs):
    samfile = pysam.AlignmentFile(file, "rb")
    kwargs["start"] = max(kwargs["start"], 0)
    y = samfile.count_coverage(read_callback=filter_strand(strand), **kwargs)
    y = np.asarray(y).sum(axis=0)
    if size_factor:
        y = y * size_factor
    s = kwargs["start"] if "start" in kwargs else 0
    trace = go.Scatter(
        x=np.arange(s, s+len(y)),
        y=y,
        name=name,
        line=dict(color=color),
        mode="lines"
    )
    return trace


def _sum_trace(y, color, name, start, step: int = 1):
    x = np.arange(start, start + y.shape[-1])
    x = x[::step]
    qx = np.concatenate((x, np.flip(x)), axis=0)
    kernel = np.ones(step) / step
    y = np.apply_along_axis(lambda m: np.convolve(m, kernel, mode="same"), axis=-1, arr=y)
    y = y[:, ::step]
    mean_y = np.mean(y, axis=0)
    median_y = np.median(y, axis=0)
    q25 = np.quantile(y, q=0.25, axis=0)
    q75 = np.quantile(y, q=0.75, axis=0)
    quantiles = np.concatenate((q75, np.flip(q25)), axis=0)
    scatter = go.Scatter if step < 5 else go.Scattergl

    traces = [
        scatter(
            x=x,
            y=median_y,
            name=f"Median",
            line=dict(dash="dot", color=color),
            mode="lines",
            legendgroup=name,
            legendgrouptitle=dict(text=name),

        ),
        scatter(
            x=x,
            y=mean_y,
            name=f"Mean",
            line=dict(color=color),
            mode="lines",
            legendgroup=name

        )
    ]
    a_color = css_color_to_rgba(color, 0.4)

    errors = [
        scatter(
            x=qx,
            y=quantiles,
            marker=dict(color=color),
            name=f"Q.25-Q.75",
            fill="toself",
            fillcolor=a_color,
            line=dict(color='rgba(255,255,255,0)'),
            legendgroup=name

        )
    ]
    return traces, errors


def get_summary_trace(files, name, color, strand, size_factors=None, step=1, **kwargs):
    t = []
    start = kwargs["start"] if "start" in kwargs else 0
    kwargs["start"] = max(kwargs["start"], 0)
    for idx, file in enumerate(files):
        samfile = pysam.AlignmentFile(file, "rb")
        y = samfile.count_coverage(read_callback=filter_strand(strand), **kwargs)
        y = np.asarray(y).sum(axis=0)
        if size_factors:
            y = y * size_factors[idx]

        t.append(y)
    y = np.asarray(t)
    traces, errors = _sum_trace(y, color, name, start, step=step)

    return traces, errors

def _compute_lines(df):
    lines = []
    df = df.sort_values(by=["strand", 'start'])

    # List to hold line assignments
    df['line'] = np.nan

    for idx, row in df.iterrows():
        start, end = row['start'], row['end']
        placed = False

        # Try to place the span in an existing line
        for i, line_end in enumerate(lines):
            if start > line_end:
                df.at[idx, 'line'] = i
                lines[i] = end
                placed = True
                break

        # If no suitable line was found, create a new line
        if not placed:
            df.at[idx, 'line'] = len(lines)
            lines.append(end)
    df["line"] = df["line"].astype(int)
    return df

def _collect_children(gff, parent_id, idx):
    children = [r for r in gff.children(parent_id) if r.id != parent_id]  # Get the direct children of the parent
    result = [(idx, child) for child in children]  # Store the children with their idx
    for child in children:
        # Recursively collect the children of each child
        result += _collect_children(gff, child.id, idx)

    return result

def deduplicate_preserve_order(items):
    return list(dict.fromkeys(items))


def _add_gffutils_entries(gff, gff_name, wstart, wend, contig, type_colors, arrow_size=None, line_mapping=None):
    if type_colors is None:
        type_colors = GFF3_TYPE_TO_COLOR
    else:
        type_colors = GFF3_TYPE_TO_COLOR | type_colors
    v = 0.4
    if arrow_size is None:
        arrow_size = 0.01 * (wend - wstart)

    data = list(gff.region(start=wstart, end=wend, seqid=contig, featuretype="gene"))
    traces = []
    annotations = []
    genes = []
    remainder = []
    idx = 0
    forbidden = set()
    m_idx = 0
    for row in data:
        if line_mapping:
            idx = line_mapping[row.id]
            if idx > m_idx:
                m_idx = idx
        genes.append((idx, row))
        annotations.append(
            dict(
                text=f'{",".join(row.attributes[gff_name]) if gff_name in row.attributes else row.id}',
                x=(row.end - row.start) / 2 + row.start,
                y=idx,
                showarrow=False,
            )

        )
        children = deduplicate_preserve_order(_collect_children(gff, row.id, idx))
        #children = list(gff.children(row.id))
        genes += children
        forbidden = forbidden | set([child.id for _, child in children] + [row.id])
        idx += 1
    data = list(gff.region(start=wstart, end=wend, seqid=contig))
    for row in data:
        if row.id not in forbidden:
            if line_mapping:
                if row.id not in line_mapping:

                    print(f'KeyError(f"ID {row.id} not found in the line mapping")')
                    continue
                idx = line_mapping[row.id]
                if idx > m_idx:
                    m_idx = idx
            genes.append((idx, row))
            idx += 1
    idx = 0
    for idx, row in genes:
        trace = _trace_from_row(row, idx, v, arrow_size, type_colors, gff_name)
        traces.append(
            trace
        )
    max_idx = m_idx if line_mapping else idx
    return traces, annotations, max_idx



def _add_gff_entries(gff, gff_name, wstart, wend, type_colors, arrow_size=None):
    if type_colors is None:
        type_colors = GFF3_TYPE_TO_COLOR
    else:
        type_colors = GFF3_TYPE_TO_COLOR | type_colors
    gff = gff[~gff[gff_name].isna()]
    gff = _compute_lines(gff)
    traces, annotations = [], []
    v = 0.4
    if arrow_size is None:
        arrow_size = 0.01 * (wend - wstart)

    for (_, row) in gff.iterrows():
        idx = row.line
        trace = _trace_from_row(row, idx, v, arrow_size, type_colors, gff_name)
        traces.append(
            trace
        )

        if row["featuretype"] not in DO_NOT_ANNOTATE:
            annotations.append(
                dict(
                    text=f'{row[gff_name]}',
                    x=(row["end"] - row["start"]) / 2 + row["start"],
                    y=idx,
                    showarrow=False,
                )

            )
    return traces, annotations, gff["line"].max()


def _trace_from_row(row, idx, v, arrow_size, type_colors, gff_name):
    if row.strand == "+":
        y = [idx - v, idx - v, idx, idx + v, idx + v, idx - v]
        pe = max(row.end - arrow_size, row.start)
        x = [row.start, pe, row.end, pe, row.start, row.start]
    elif row.strand == "-":
        y = [idx, idx - v, idx - v, idx + v, idx + v, idx]
        ps = min(row.start + arrow_size, row.end)
        x = [row.start, ps, row.end, row.end, ps, row.start]
    else:
        y = [idx - v, idx - v, idx + v, idx + v, idx - v]
        x = [row.start, row.end, row.end, row.start, row.start]
    color = type_colors.get(row.featuretype, type_colors["default"])
    if isinstance(row, gffutils.Feature):
        name = row.attributes[gff_name][0] if gff_name in row.attributes else row.id
        hovertext = "<br>".join([f"{key}: {','.join(value)}" for key, value in row.attributes.items()])
    else:
        name = row[gff_name]
        hovertext=None
    trace = go.Scatter(
        y=y,
        x=x,
        mode="lines",
        fillcolor=color,
        line=dict(color="black", width=1),
        fill="toself",
        name=f'{row.featuretype}-{name}',
        text=hovertext,
        legendgroup="Annotations",
        legendgrouptitle=dict(text="Annotations")

    )
    return trace


def traces_from_precomputed_coverage(design, contig, coverages: Dict[str, np.ndarray], wstart, wend, strand, step, colors=None):
    tmp = design.groupby(["Treatment"], as_index=False, observed=False).agg({"File": list, "index": list})
    if colors is None:
        colors = {trace: px.colors.qualitative.Light24[idx] for idx, trace in enumerate(design.Treatment)}
    traces = []
    errors = []
    for idx, row in tmp.iterrows():
        name = row["Treatment"]
        color = colors[name]
        y = coverages[contig][strand][row["index"], wstart:wend]
        t, e = _sum_trace(y, color, name, wstart, step)
        traces += t
        errors += e
    return traces, errors


def traces_from_design(design, contig, wstart, wend, strand, colors, step):
    tmp = design.groupby(["Treatment"], as_index=False, observed=False).agg({"File": list, "Size factors": list}
                                                                                ).dropna().reset_index()
    t = []
    e = []
    for idx, row in tmp.iterrows():
        name = row["Treatment"]
        color = colors[row["Treatment"]]
        sfs = row["Size factors"]

        traces, errors = get_summary_trace(row["File"], name=name, contig=contig, start=wstart, stop=wend,
                                           strand=strand, color=color, size_factors=sfs, step=step)
        t += traces
        e += errors
    return t, e


def filter_strand(strand):
    if strand == "+":
        def filter_read(read):
         return not read.is_reverse
    elif strand == "-":
        def filter_read(read):
            return read.is_reverse
    else:
        return "all"
    return filter_read



def plot_precomputed_coverage(
        design,
        contig,
        wstart,
        wend,
        gff = None,
        gff_name = None,
        coverages: Dict[str, np.ndarray] = None,
        colors=None,
        type_colors=None,
        arrow_size=None,
        step: int = 1,
        show_annotations: bool = True,
        line_mapping: Dict = None,
        force: bool =False,
        strand: str = "both",
        **kwargs):
    wstart = max(int(wstart), 0)
    rows = 2 if strand in ["+", "-"] else 3
    fig = make_subplots(rows=rows, shared_xaxes=True, **kwargs)
    traces = []
    errors = []
    if strand != "-":
        if coverages is not None :
            wend = min(int(wend), coverages[contig]["+"].shape[-1])
            traces, errors = traces_from_precomputed_coverage(design, contig, coverages, wstart, wend, "+", step, colors)

        else:
            if force or (wend - wstart) <= 4001:
                traces, errors = traces_from_design(design, contig, wstart, wend, "+", colors, step)
            else:
                traces = [go.Scatter(x=[None], y=[None], showlegend=False)]
                errors = []
                fig.add_annotation(
                    xref="paper",
                    yref="y domain",
                    xanchor="center",
                    yanchor="middle",
                    x=0.5,
                    y=0.5,
                    text="Decrease window size or use precomputed coverage for plotting",
                    showarrow=False,
                    font=(dict(size=28))
                )

        fig.add_traces(traces, rows=1, cols=1)
        fig.add_traces(errors, rows=1, cols=1)
    if strand != "+":
        if coverages is not None:
            traces, errors = traces_from_precomputed_coverage(design, contig, coverages, wstart, wend, "-", step, colors)
        else:
            if force or (wend - wstart) <= 4001:
                traces, errors = traces_from_design(design, contig, wstart, wend, "-", colors, step)
            else:
                traces = [go.Scatter(x=[None], y=[None], showlegend=False)]
                errors = []
        if strand != "-":
            for trace in traces:
                trace.update(showlegend=False)
            for error in errors:
                error.update(showlegend=False)
        row = 1 if strand == "-" else 3
        fig.add_traces(traces, rows=row, cols=1)
        fig.add_traces(errors, rows=row, cols=1)
    if gff is not None:
        if isinstance(gff, pd.DataFrame):
            gff = filter_gff_by_interval(gff, contig, wstart, wend)
            gff_traces, annotations, max_idx = _add_gff_entries(gff, gff_name, wstart, wend, type_colors, arrow_size)
        elif isinstance(gff, gffutils.FeatureDB):
            gff_traces, annotations, max_idx = _add_gffutils_entries(gff, gff_name, wstart, wend, contig, type_colors, arrow_size, line_mapping)
        else:
            raise NotImplementedError("Non supported feature database")
        fig.add_traces(
            gff_traces,
            rows=2,
            cols=1
        )
        if show_annotations:
            for anno in annotations:
                fig.add_annotation(
                    anno,
                    row=2, col=1
                )
    else:
        max_idx = np.nan
    fig.update_xaxes(
        range=[wstart, wend],
        minallowed=[wstart, wend],
        #maxallowed=coverages[contig]["+"].shape[-1]
    )
    if np.isnan(max_idx):
        fig.add_annotation(
            text="", showarrow=False,
            row=2, col=1
        )
        max_idx = 2
    fig.update_yaxes(
        minallowed=-0.5-1,
        maxallowed=max_idx+0.5+1,
        showticklabels=False,
        row=2,
    )
    fig.update_yaxes(
        title=dict(text=f"Average Read count<br>[kernel size {step}]", standoff=50),
        row=2
    )
    fig.update_layout(legend=dict(groupclick="toggleitem"))
    return fig



def plot_coverage(design, contig, start, stop, extend: int = 50, strand: str = None, colors = None, mode: str = "replicate", gff: pd.DataFrame = None, gff_types = None, gff_name: str = "locus_tag", type_colors = None, **kwargs):
    if colors is None:
        colors = DEFAULT_PLOTLY_COLORS
    sf = True if "Size factors" in design else False

    fig = make_subplots(rows=2, cols=1, shared_xaxes=True, **kwargs)
    wstart = start - extend
    wend = stop + extend
    if gff_types:
        gff = gff[gff["type"].isin(gff_types)]
    if strand:
        gff = gff[gff["strand"] == strand]

    color_idx = {treat: idx for idx, treat in enumerate(design["Treatment"].unique())}
    if mode == "replicate":
        for idx, row in design.iterrows():
            name = f"{row['Treatment']} - {row['Replicate']}"
            color = colors[color_idx[row["Treatment"]]]
            size_factor = row["Size factor"] if sf else None
            trace = get_trace(row["File"], name=name, contig=contig, start=wstart, stop=wend, strand=strand, color=color, size_factor=size_factor)
            fig.add_trace(trace)
    else:
        if sf:
            tmp = design.groupby(["Treatment"], as_index=False, observed=False).agg({"File": list, "Size factors": list}
                ).dropna().reset_index()
        else:
            tmp = design.groupby(["Treatment"], as_index=False, observed=False)["File"].agg(
                list).dropna().reset_index()

        t = []
        e = []
        for idx, row in tmp.iterrows():
            name = row["Treatment"]
            color = colors[color_idx[row["Treatment"]]]
            sfs = row["Size factors"] if sf else None

            traces, errors = get_summary_trace(row["File"], name=name, contig=contig,  start=wstart, stop=wend, strand=strand, color=color, size_factors=sfs)
            t += traces
            e += errors
        fig.add_traces(
            e+t,
            rows=1,
            cols=1
        )
    if gff:
        gff = filter_gff_by_interval(gff, contig, wstart, wend)
        gff_traces, annotations, indices = _add_gff_entries(gff, gff_name, wstart, wend, type_colors)
        fig.add_traces(
            gff_traces,
            rows=2,
            cols=1
        )
        for anno in annotations:
            fig.add_annotation(
                anno,
                row=2, col=1
            )
        fig.update_yaxes(
            tickvals=list(range(len(indices))),
            ticktext=list(indices.keys()),
            tickmode="array",
            row=2
        )

    fig.update_xaxes(range=[wstart, wend])

    return fig


if __name__ == '__main__':

    chr = "BA000022.2"
    start = 7179
    stop = 8311
    import os

    files = [os.path.abspath(os.path.join("../testData/files/", file)) for file in os.listdir("../testData/files/") if file.endswith(".bam")]

    design = {
        "File": [],
        "Replicate": [],
        "Treatment": [],

    }
    coverage = {
    }
    tmp = []
    for file in files:
        s = os.path.basename(file)
        s = s.split("_")
        if len(s) == 8:
            rep = s[4]
            treat = f"+p {s[3]}"
            if treat == "+p TC":
                continue
        elif len(s) == 7:
            rep = s[3]
            treat = f"-p {s[2]}"
        else:
            rep = s[-1].split(".")[0]
            treat = s[0] + "_" + s[1]
            continue
            # if "totalRNA" in treat:
            #     continue
        design["Replicate"].append(rep)
        design["Treatment"].append(treat)
        design["File"].append(file)


    design = pd.DataFrame(design)
    design["index"] = list(range(len(design)))

    import pickle
    file = "pickled_cov.pckl"
    if not os.path.exists(file):
        coverage = precompute_from_design(design)

        with open(file, "wb") as handle:
            pickle.dump(coverage, handle)
    else:
        with open(file ,"rb") as handle:
            coverage = pickle.load(handle)
    #design["Size factors"] = list(range(len(design)))
    gff = "../testData/20210217_syne_onlyUnique_withFeat.gff3"
    lt2name = "../testData/TableWithSequences.tsv"
    lt2name = pd.read_csv(lt2name, sep="\t")
    gff = read_gff3(gff)
    gff = gff.merge(lt2name, on="locus_tag", how="left")
    gff.loc[gff.gene_name.isna(), "gene_name"] = gff.loc[gff.gene_name.isna()].locus_tag
    sl = "psba2"
    d = find_gene_coordinates(gff, gene=sl, anno_type="gene_name")
    #d = d.iloc[1]
    #d = d[d["type"] == "CDS"]
    for _, row in d.iterrows():
        _, start, stop, strand = coords_from_gff_row(row)
        plot_precomputed_coverage(design, contig=chr, coverages=coverage, wstart=7000, wend=9000,gff=gff, gff_name="gene_name")
        exit()
        fig = plot_coverage(design, chr, start, stop, extend=3000, strand=strand, mode="f", gff=gff, gff_types=None,
                            gff_name="gene_name",
                            colors=None,  # ("rgb(0, 142, 151)", "rgb(252, 76, 2)"),
                            type_colors={"CDS": "rgb(252, 76, 2)"},
                            vertical_spacing=0)
        fig.show()
        break
