import pandas as pd
import plotly.graph_objects as go
import numpy as np
from scripts.plotly_template import register_template

register_template()

def _color_to_calpha(color: str, alpha: float = 0.2):
    color = color.split("(")[-1].split(")")[0]
    return f"rgba({color}, {alpha})"

def main(file, probes, cfg, layout, subplot_params):
    df = pd.read_csv(file, sep="\t")
    id_vars = ["Mask", "Image", "Probe", "Figure", "Treatment"]
    mode = "mean"
    cols = [col for col in df.columns if mode in col and "bg_mean" not in col]
    d = {col: col.split("_")[1] for col in cols}
    df = df.rename(d, axis=1)
    cols = [p for p in d.values() if int(p) <=20]

def main(tsv, cfg, outfile):
    df = pd.read_csv(tsv, sep="\t")
    fig = go.Figure()
    mp_m = df[df["minusPMembran"]]
    mp_c = df[df["minusPCytosol"]]
    others = df[~(df["minusPMembran"] | df["minusPCytosol"])]

    colors = [cfg["seahawksGrey"], cfg["membraneColor"], cfg["cytosolColor"]]
    names = ["Not sig.", "Membrane", "Cytosol"]
    lg  = "untreated"
    for idx, subdf in enumerate([others, mp_m, mp_c, ]):
        color = colors[idx]
        fig.add_trace(
            go.Scatter(
                x=subdf["log2FoldChange_plusPuromycin"],
                y=subdf["log2FoldChange_tc"],
                mode="markers",
                marker = dict(
                    color=np.abs(subdf["log2FoldChange_minusPuromycin"]),
                    colorscale=[_color_to_calpha(color, 0.2), _color_to_calpha(color, 1)],
                    cmin=0,
                    cmax=np.abs(subdf["log2FoldChange_minusPuromycin"]).max(),
                    showscale=False  # <-- do NOT show a colorbar
                ),
                name=names[idx],
                legendgroup=lg,
                legendgrouptitle=dict(text=lg),
                hovertext=subdf["gene_name"],
                customdata=np.stack((
                    subdf["gene_name"],
                    subdf["log2FoldChange_minusPuromycin"],
                    subdf["padj_plusPuromycin"],
                    subdf["padj_tc"],
                ), axis=-1),

                # Custom hover template
                hovertemplate=(
                    "<b>%{customdata[0]}</b><br>"
                    "log2FC +Puromycin: %{x:.2f}<br>"
                    "log2FC TC: %{y:.2f}<br>"
                    "log2FC –Puromycin: %{customdata[1]:.2f}<br>"
                    "padj +Puromycin: %{customdata[2]:.2f}<br>"
                    "padj TC: %{customdata[3]:.2f}"
                    "<extra></extra>"
                ),

            )
        )
    fig.update_xaxes(title="Puromycin log<sub>2</sub>FoldChange")
    fig.update_yaxes(title="Total Cell log<sub>2</sub>FoldChange")
    color = "rgb(0, 0,0)"
    cmin = 0
    cmax = np.abs(subdf["log2FoldChange_minusPuromycin"]).max()
    fig.add_trace(
        go.Scatter(
            x=[None], y=[None],
            mode="markers",
            marker=dict(
                color=np.linspace(cmin, cmax, 10),
                colorscale="Greys",  # colorbar gradient
                cmin=0,
                cmax=cmax,
                showscale=True,
                colorbar=dict(
                    title="abs(log2FC untreated)",
                    len=0.4
                )
            ),
            showlegend=False
        )
    )

    fig.update_layout(
        coloraxis=dict(
            colorscale=[_color_to_calpha(color, 0.2), _color_to_calpha(color, 1)],
            cmin=np.abs(subdf["log2FoldChange_minusPuromycin"]).min(),
            cmax=np.abs(subdf["log2FoldChange_minusPuromycin"]).max(),
            colorbar=dict(
                title="abs(log2FC untreated)",
                thickness=15,
                len=0.8,
                x=1.05,
                y=0,
                yref="paper",
                yanchor="bottom",
            )
        ),
        xaxis_title="+ Puromycin log<sub>2</sub>FoldChange",
        yaxis_title="Total Cell log<sub>2</sub>FoldChange"
    )
    fig.update_layout(cfg["LayoutDefault"])
    fig.update_layout(cfg["SupplementScatter"]["layout"])
    fig.write_image(outfile)
    p = 0




if __name__ == '__main__':
    if "snakemake" in globals():
        file = snakemake.input.tsv
        cfg = snakemake.config
        outfile = snakemake.output.svg

    else:
        import yaml
        file = "Pipeline/Summary/joinedDESeq_Table.tsv"
        intensities = "Pipeline/FISH/single_spots.tsv"
        outfile = "foo.svg"
        with open("configFile.yml", 'r') as f:
            cfg = yaml.safe_load(f)

    main(file, cfg, outfile)

