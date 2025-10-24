import os

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import tifffile
import plotly.express as px
import numpy as np
import plotly.io as pio
from scipy import stats
from scripts.plotFISHIntensity import calc_summary_and_pval

pio.renderers.default = "browser"


def _color_to_calpha(color: str, alpha: float = 0.2):
    color = color.split("(")[-1].split(")")[0]
    return f"rgba({color}, {alpha})"

def main(file, probes, cfg, layout, subplot_params, intensities=None):
    df = pd.read_csv(file, sep="\t")
    if intensities is not None:
        intensities = pd.read_csv(intensities, sep="\t")
        probe_names = [probe[0] for probe in probes]
        intensities, results_df = calc_summary_and_pval(intensities, probe_names)
        specs = [
            [{"type": "image"} for _ in range(len(probes))],
            [{"type": "box", "colspan": 4}, None, None, None] ,# Row 2: one boxplot spanning both cols
            # Row 1: two image plots
            [{"type": "xy"} for _ in range(len(probes))],  # Row 1: two image plots
        ]
        row = 3
    else:
        specs = None
        row = 2
    id_vars = ["Mask", "Image", "Probe", "Figure", "Treatment"]
    mode = "mean"
    cols = [col for col in df.columns if mode in col and "bg_mean" not in col]
    d = {col: col.split("_")[1] for col in cols}
    df = df.rename(d, axis=1)
    cols = [p for p in d.values() if int(p) <=20]
    df = df[cols + id_vars]

    df[cols] = (df[cols].values -  df[cols].values.min(axis=1)[:, None])/ (df[cols].values.max(axis=1)[:, None] -  df[cols].values.min(axis=1)[:, None])


    df["max"] = df[cols].values.argmax(axis=-1) * 26

    sdf = df.melt(id_vars=id_vars)

    summary = sdf.groupby(['variable', 'Mask', "Probe", "Figure", "Treatment"])['value'].agg(
        mean='mean',
        q25=lambda x: np.nanquantile(x, 0.25),
        q75=lambda x: np.nanquantile(x, 0.75),
        count=lambda x: x.notna().sum(),
        std=lambda x: np.nanstd(x, ddof=1)
    ).reset_index()
    summary['sem'] = summary['std'] / np.sqrt(summary['count'])
    summary['ci95'] = summary['sem'] * stats.t.ppf(0.975, df=summary['count'] - 1)
    summary = summary[summary["variable"] != "max"]
    summary = summary[summary["variable"].astype(int) <= 20]
    col_titles = [p[0] for p in probes]
    fig = make_subplots(
        rows=2 if intensities is None else 3,
        cols=len(probes),
        row_heights=[0.6, 0.4]if intensities is None else [0.4, 0.3, 0.3],
        x_title="distance to probe [nm]",
        specs=specs,
        **subplot_params


    )
    fig.update_layout(cfg["LayoutDefault"])

    color_map = {"Replicate 1": "rgb(0, 106, 56)", "Replicate 2": "rgb(0, 48, 135)", "Replicate 3": "rgb( 111, 160, 130)"}
    for i, (probe, treatment, replicate) in enumerate(probes):
        traces = []
        errors = []
        image_dir = f"initialData/AllFISHData/{probe}/{replicate}"
        files = [f for f in os.listdir(image_dir) if f.endswith(f"{treatment}_green.tif") and "Image" in f]
        print(image_dir)
        image_file = os.path.join(image_dir, files[0])
        img = tifffile.imread(image_file)
        tmp = px.imshow(img)
        fig.add_trace(tmp.data[0], row=1, col=i + 1)
        amaxes = []
        text = "<b>n: </b>"
        subdf = summary[summary["Probe"] == probe]
        subdf = subdf[subdf["Treatment"] == treatment]
        for j, f in enumerate(sorted(subdf["Figure"].unique())):
            color = color_map[f]
            sdf = subdf[subdf["Figure"] == f]
            sdf = sdf[sdf["Mask"] == "Measured"]
            print(sdf)
            if len(sdf) == 0:
                continue
            if j > 0:
                text += ", "
            c = sdf["count"].iloc[0]
            text += f"<span style='color:{color}'>{c}</span>"

            sdf["variable"] = sdf.variable.astype(int)
            sdf = sdf.sort_values("variable")
            x = sdf['variable'] * 26
            p = go.Scatter(
                x=x,
                y=sdf['mean'],
                mode='lines',
                name=f'{probe} Mean',
                line=dict(width=2, color=color),
                showlegend=False,
            )
            amaxes.append(x.iloc[int(sdf["mean"].argmax())])
            traces.append(p)
            iqr = pd.concat([sdf['q75'], sdf['q25'][::-1]])
            ci_mean = pd.concat([sdf['mean'] - sdf["ci95"], (sdf['mean'] + sdf["ci95"])[::-1]])
            errors.append(go.Scatter(
                x=pd.concat([sdf['variable'], sdf['variable'][::-1]]) * 26,
                y=ci_mean,
                fill='toself',
                fillcolor=_color_to_calpha(color, 0.2),  # You can make this color dynamic per mask
                line=dict(color='rgba(255,255,255,0)'),
                hoverinfo="skip",
                name=f'{probe} 95% CI',
                showlegend=False,
            ))
        fig.add_traces(
            traces + errors, cols=i+1, rows=2 if intensities is None else 3
        )
        if intensities is not None:
            group = intensities[(intensities.Treatment == treatment) & (intensities.Probe == probe)]
            values = group["ProbeIntensity"].tolist()
            mean_val = np.mean(values)
            std_val = np.std(values)
            fig.add_trace(
                go.Bar(
                    x=[f"{treatment} – {probe}"],
                    y=[mean_val],
                    name=f"{treatment} – {probe}",
                    marker=dict(color="#69BE28"),
                    showlegend=False,

                    error_y=dict(
                        type="data",
                        array=[std_val],
                        visible=True,
                        color="black",  # <-- sets error bar line color
                        thickness=1.5,  # optional: makes the line a bit thicker
                        width=4
                    )
                ),
                row=2,
                col=1
            )
        vint = float(np.median(np.asarray(amaxes)))
        didx = i+1+len(probes) if intensities is None else i+1 + len(probes) + 1
        domain = f"{didx} domain"
        fig.add_annotation(
            x=1, y=1,
            yref=f"y{domain}",
            xref=f"x{domain}",
            yanchor="top",
            xanchor="right",
            text=text,
            showarrow=False,
            bgcolor="white",
            opacity=1,

        )
        fig.add_vline(vint, col=i+1, row=row, line=dict(dash="dash", color="black", width=1), layer="above")
    colors = list(color_map.values())
    fig.add_traces(
        [
            go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                name=f"{i+1}",
                marker=dict(color=colors[i], symbol="square", size=10),
                legendgroup="rep",
                legendgrouptitle=dict(text="Replicate", font=dict(size=10)),

                legendrank=2

            ) for i in range(len(colors))
        ]
    )
    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            legendgroup="Data",
            legendgrouptitle=dict(text="Data", font=dict(size=10)),
            mode="lines",
            name=f"Mean",
            line=dict(color='black'),
            legendrank=1

        )
    )
    if intensities is not None:
        fig.add_trace(
            go.Bar(
                x=[None],
                y=[None],
                name=f"Intensity",
                marker=dict(color="#69BE28"),
                showlegend=True,
                legendgroup="Intensities",
                legendgrouptitle=dict(text="Intensities"),

            ),
            row=2,
            col=1
        )
        fig.update_xaxes(
            showticklabels=False, row=2
        )
        fig.update_yaxes(
            title=dict(text="mean<br>Intensity"), row=2
        )


    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="lines",
            name=f"95% CI",
            fill="toself",
            fillcolor='rgba(0,0,0,0.4)',
            line=dict(color='rgba(255,255,255,0)'),
            legendgroup="Data",
            legendrank=2

        )
    )

    fig.update_xaxes(
        showgrid=False, showticklabels=False, row=1
    )

    fig.update_yaxes(
        showgrid=False, showticklabels=False, row=1
    )
    fig.update_yaxes(
        range=[0, 1], row=row, showticklabels=False,
    )

    fig.update_layout(
        width=621,
        height=400,
    )

    fig.update_yaxes(
        title=dict(text="normalized<br>intensity"),
        row=row,
        col=1,
        showticklabels=True,
    )
    fig.update_layout(**layout)
    for shape in fig.layout.shapes:
        shape['layer'] = 'above'
    for annotation in fig.layout["annotations"]:
        annotation.update(
            layout["annotations"]
        )
    return fig








if __name__ == '__main__':
    if "snakemake" in globals():
        file = snakemake.input.tsv
        probes = snakemake.params.probes
        layout = snakemake.params.layout
        cfg = snakemake.config
        outfile = snakemake.output.svg
        subplot_params = snakemake.params.subplot
    else:
        import yaml
        file = "Pipeline/FISH/all_distributions.tsv"
        probes = [("psaA", "untreated", "Replicate 1"), ("cupA", "untreated", "Replicate 1"), ("rnpB", "untreated", "Replicate 1"), ("fusB", "untreated", "Replicate 1") ]
        outfile = "foo.svg"
        subplot_params = dict(vertical_spacing=0.05, column_titles=[p[0] for p in probes],)
        with open("configFile.yml", 'r') as f:
            cfg = yaml.safe_load(f)
        layout = cfg["F2FISH"]["layout"]
    fig = main(file, probes, cfg, layout, subplot_params)
    if "snakemake" not in globals():
        fig.show()
    fig.write_image(outfile)