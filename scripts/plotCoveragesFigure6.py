import pandas as pd
from scripts.readCount import plot_precomputed_coverage
from scripts.gff_helper import read_gff3, read_gff_via_gffutils, optimize_gff_lines
import plotly.graph_objs as go

def main(design, gff, db_name, output_path, cfg):
    gff = read_gff_via_gffutils(gff, dbname=db_name, fix=False)
    lmap, attributes = optimize_gff_lines(gff)
    psba_start = 7229
    psba_end = 8311

    rpl10_12start = 922844
    rpl10_12end = 923839



    seqid = "BA000022.2"

    lmap["5'UTR_5"] = 0
    lmap["asRNA-slr1311-as-gene"] = 1

    fig = plot_precomputed_coverage(
        design,
        contig=seqid,
        wstart=psba_start - 100,
        wend=psba_end + 50,
        strand="+",
        gff=gff,
        gff_name="gene",
        colors=cfg["TRACE_COLORS"],
        type_colors=cfg["RNA_TYPE_COLORS"],
        force=True,
        line_mapping=lmap,
        vertical_spacing=0,
        row_heights=[0.9, 0.1],
        cols=2,
    )
    fig.update_layout(cfg["LayoutDefault"])
    already = set()
    for trace in fig.data:
        if trace.legendgroup == "Annotations":
            type = trace.name.split("-")[0]
            if type not in already and type != "gene" and type != "asRNA":
                trace.update(name=type)
                already.add(type)
            else:
                trace.update(showlegend=False)
        else:
            trace.update(line=dict(width=1), showlegend=False)
    fig.layout.annotations = [anno for anno in fig.layout.annotations if not anno.text.startswith("asRNA")]

    fig.update_yaxes(row=2, showgrid=False)
    fig2 = plot_precomputed_coverage(
        design,
        contig=seqid,
        wstart=rpl10_12start - 50,
        wend=rpl10_12end + 50,
        strand="-",
        gff=gff,
        gff_name="gene",
        colors=cfg["TRACE_COLORS"],
        type_colors=cfg["RNA_TYPE_COLORS"],
        force=True,
        line_mapping=lmap,
        vertical_spacing=0,
        row_heights=[0.9, 0.1],
        cols=1
    )
    for trace in fig2.data:
        if trace.xaxis == "x":
            trace.update(line=dict(width=1), showlegend=False)
            #trace.update(xaxis="x2", yaxis="y2", name=trace.name + "foo", showlegend=False)
            fig.add_trace(trace, col=2, row=1)
        else:
            trace.update(showlegend=False)
            fig.add_trace(trace, col=2, row=2)
    for anno in fig2.layout.annotations:
        anno.update(xref="x4", yref="y4")
        fig.add_annotation(anno)
    fig.update_yaxes(range=[-.5, .5], row=2)
    fig.update_yaxes(title=None)

    fig.update_yaxes(title="normalized read count", row=1, col=1)
    fig.update_yaxes(title=None, col=2)
    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="markers",
            name=f"Membrane",
            marker=dict(color=cfg["TRACE_COLORS"]["-p M"], symbol="square", size=10),
            legendgroup="untreated",
            legendgrouptitle=dict(text="untreated"),

        )
    )
    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="markers",
            name=f"Cytosol",
            marker=dict(color=cfg["TRACE_COLORS"]["-p C"], symbol="square", size=10),
            legendgroup="untreated",
            legendrank=1

        )
    )
    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="markers",
            name=f"Membrane",
            marker=dict(color=cfg["TRACE_COLORS"]["+p M"], symbol="square", size=10),
            legendgroup="+ Puromycin",
            legendrank=2

        )
    )
    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="markers",
            name=f"Cytosol",
            marker=dict(color=cfg["TRACE_COLORS"]["+p C"], symbol="square", size=10),
            legendgroup="+ Puromycin",
            legendgrouptitle=dict(text="+ Puromycin"),

            legendrank=1

        )
    )
    lgroup = "Lines"
    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            legendgroup="Lines",
            legendgrouptitle=dict(text="Lines"),
            mode="lines",
            name=f"Mean",
            line=dict(color='black'),

        )
    )
    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="lines",
            name=f"Median",
            line=dict(dash="dot", color="black"),
            legendgroup=lgroup,

        )
    )
    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="lines",
            name=f"Q.25-Q.75",
            fill="toself",
            fillcolor='rgba(0,0,0,0.4)',
            line=dict(color='rgba(255,255,255,0)'),
            legendgroup=lgroup,

        )
    )

    fig.update_xaxes(range=[rpl10_12start - 50, rpl10_12end + 50],  col=2, )
    fig.update_layout(cfg["F6Coverage"]["layout"])
    fig.write_image(output_path)

    p = 0

if __name__ == '__main__':
    file = snakemake.input.design
    gff = snakemake.input.gff
    cfg = snakemake.config
    sf = snakemake.input.sf

    design = pd.read_csv(file, sep="\t")
    sf = pd.read_csv(sf, sep="\t")
    sf = sf.rename({"x": "Size factors"}, axis=1)
    design = design[design["Replicate"].isin(sf.index) & design["Replicate"] != "plusP_M_R3"]
    design = design.merge(sf, left_on="Replicate", right_index=True)
    main(design, gff, output_path=snakemake.output.svg, cfg=cfg, db_name=snakemake.output.db)
