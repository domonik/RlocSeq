import pandas as pd
from scripts.readCount import plot_precomputed_coverage
from scripts.gff_helper import read_gff3, read_gff_via_gffutils, optimize_gff_lines
import plotly.graph_objs as go


def main(design, gff, db_name, output_path, cfg):
    gff = read_gff_via_gffutils(gff, dbname="foofile.db", fix=False)
    lmap, attributes = optimize_gff_lines(gff)
    psaa_start = 941686
    psab_end = 946382

    cupA_start = 935409
    cupA_end = 936704

    rnpb_start = 153165

    rnpb_end = 153603

    t13_start = 1329163

    t13_end = 1329236

    #design["Size factors"] = 1


    seqid = "BA000022.2"

    fig = plot_precomputed_coverage(
        design,
        contig=seqid,
        wstart=psaa_start - 150,
        wend=psab_end + 100,
        strand="+",
        gff=gff,
        gff_name="gene",
        colors=cfg["TRACE_COLORS"],
        type_colors=cfg["RNA_TYPE_COLORS"],
        force=True,
        line_mapping=lmap,
        vertical_spacing=0,
        row_heights=[0.9, 0.1],
        cols=3,

    )
    fig.update_layout(cfg["LayoutDefault"])


    f2 = plot_precomputed_coverage(
        design,
        contig=seqid,
        wstart=cupA_start - 200,
        wend=cupA_end + 100,
        strand="-",
        gff=gff,
        gff_name="gene",
        colors=cfg["TRACE_COLORS"],
        type_colors=cfg["RNA_TYPE_COLORS"],
        force=True,
        line_mapping=lmap,
        vertical_spacing=0,
        row_heights=[0.9, 0.1]
    )
    f3 = plot_precomputed_coverage(
        design,
        contig=seqid,
        wstart=t13_start - 50,
        wend=t13_end + 50,
        strand="-",
        gff=gff,
        gff_name="gene",
        colors=cfg["TRACE_COLORS"],
        type_colors=cfg["RNA_TYPE_COLORS"],
        force=True,
        line_mapping=lmap,
        vertical_spacing=0,
        row_heights=[0.9, 0.1]
    )


    fig.update_layout(cfg["LayoutDefault"])
    fig.data = [d for d in fig.data if d.name not in ("5'UTR-5'UTR_274", "5'UTR-5'UTR_275") ]
    already = set()
    for trace in fig.data:
        if trace.legendgroup == "Annotations":
            type = trace.name.split("-")[0]
            if type not in already and type != "gene":
                trace.update(name=type)
                already.add(type)
            else:
                trace.update(showlegend=False)
        else:
            trace.update(line=dict(width=1), showlegend=False)
    fig.layout.annotations = [anno for anno in fig.layout.annotations if not anno.text.startswith("ncRNA")]
    fig.update_yaxes(title=None)
    fig.update_yaxes(title="normalized read count", row=1, col=1)
    fig.update_yaxes(row=2, showgrid=False)
    for trace in f2.data:
        if trace.xaxis == "x":
            trace.update(line=dict(width=1), showlegend=False)
            fig.add_trace(trace, col=2, row=1)
        else:
            trace.update(showlegend=False)
            fig.add_trace(trace, col=2, row=2)
    for anno in f2.layout.annotations:
        anno.update(xref="x5", yref="y5")
        fig.add_annotation(anno)
    fig.update_xaxes(range=[cupA_start - 50, cupA_end + 50],  col=2, )
    for trace in f3.data:
        if trace.xaxis == "x":
            trace.update(line=dict(width=1), showlegend=False)
            fig.add_trace(trace, col=3, row=1)
        else:
            trace.update(showlegend=False)
            fig.add_trace(trace, col=3, row=2)
    for anno in f3.layout.annotations:
        anno.update(xref="x6", yref="y6", text="tRNA13")
        fig.add_annotation(anno)
    fig.update_xaxes(range=[t13_start - 25, t13_end + 25],  col=3, )
    fig.update_yaxes(range=[-.5, .5], row=2)
    fig.update_yaxes(title=None, col=2)
    fig.update_yaxes(title=None, col=3)

    lgroup = "Read Counts"
    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="markers",
            name=f"Membrane",
            marker=dict(color=cfg["membraneColor"], symbol="square", size=10),
            legendgroup=lgroup,

        )
    )
    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="markers",
            name=f"Cytosol",
            marker=dict(color=cfg["cytosolColor"], symbol="square", size=10),
            legendgroup=lgroup,
            legendrank=1

        )
    )
    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            legendgroup=lgroup,
            legendgrouptitle=dict(text=lgroup),
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
    fig.update_layout(cfg["LayoutDefault"])

    fig.update_layout(cfg["MinusPuroCoverage"]["layout"])

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

    main(design, gff, output_path=snakemake.output.svg, cfg=cfg, db_name=snakemake.input.db)
