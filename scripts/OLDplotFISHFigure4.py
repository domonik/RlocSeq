import numpy as np
from scipy import stats
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import pandas as pd
import os
import tifffile
import plotly.express as px
import plotly.io as pio

pio.renderers.default = "browser"

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

    fig = make_subplots(
        rows=2,
        cols=2,
        row_heights=[0.5, 0.5],
        subplot_titles=["untreated", "Puromycin", None, None],
        specs=[
            [{"type": "xy"}] * 2,  # row 1
            [{"type": "xy"}, {"type": "table"}]
        ],
        **subplot_params


    )
    fig.update_layout(cfg["LayoutDefault"])
    treatments = ["untreated", "Puromycin"]

    color_map = {"untreated": "rgb(0, 106, 56)", "Puromycin": "rgb(0, 48, 135)"}

    traces = []
    errors = []
    probe="16S rRNA"
    for j, treatment in enumerate(treatments):

        image_dir = f"initialData/AllFISHData/{probe}/Replicate 2"
        files = [f for f in os.listdir(image_dir) if f.endswith(f"{treatment}_green.tif") and "Image" in f]
        print(image_dir)
        image_file = os.path.join(image_dir, files[0])
        img = tifffile.imread(image_file)
        tmp = px.imshow(img)
        fig.add_trace(tmp.data[0], row=1, col=j+1)
        amaxes = []
        sdf = summary[summary["Probe"] == probe]
        sdf = sdf[sdf["Mask"] == "Measured"]
        sdf = sdf[sdf["Treatment"] == treatment]
        sdf["variable"] = sdf.variable.astype(int)

        sdf = sdf.sort_values("variable")
        sdf = sdf.groupby("variable")["mean"].agg(
            mean='mean',
            q25=lambda x: np.nanquantile(x, 0.25),
            q75=lambda x: np.nanquantile(x, 0.75),
            count=lambda x: x.notna().sum(),
            std=lambda x: np.nanstd(x, ddof=1)
        ).reset_index()
        color = color_map[treatment]

        if len(sdf) == 0:
            continue
        c = sdf["count"].iloc[0]
        x = sdf['variable'] * 26
        p = go.Scatter(
            x=x,
            y=sdf['mean'],
            mode='lines',
            name=f'{probe} Mean',
            line=dict(width=2, color=color),
            showlegend=False,
        )
        amax = x.iloc[int(sdf["mean"].argmax())]
        traces.append(p)
        iqr = pd.concat([sdf['q75'], sdf['q25'][::-1]])
        #ci_mean = pd.concat([sdf['mean'] - sdf["ci95"], (sdf['mean'] + sdf["ci95"])[::-1]])
        #fig.add_vline(amax, col=i + 1, row=4, line=dict(dash="dash", color=color, width=1), layer="above", )

        errors.append(
            go.Scatter(
            x=pd.concat([sdf['variable'], sdf['variable'][::-1]]) * 26,
            y=iqr,
            fill='toself',
            fillcolor=_color_to_calpha(color, 0.2),  # You can make this color dynamic per mask
            line=dict(color='rgba(255,255,255,0)'),
            hoverinfo="skip",
            name=f'{probe} IQR',
            showlegend=False,
        )
        )
    p = summary[(summary["Probe"] == probe)]
    p = p.groupby(["Treatment", "Figure"])["count"].min().reset_index()
    p = p.pivot(index="Figure", columns="Treatment", values="count")
    fig.add_trace(
            go.Table(
                header=dict(
                    values=[t[0].upper() for t in treatments],
                    fill_color=[color_map[t] for t in treatments],
                    align="center",
                    font=dict(
                        size=8,
                        color=[
                            "white", "white", "black"
                        ]
                    ),
                    height=20,
                ),
                cells=dict(
                    values=[[("-" if pd.isna(val) else val) for val in p[col]] for col in treatments],                        align="center",
                    font=dict(size=8)
                )
            ),
        row=2, col=2
    )
    fig.add_traces(
        traces + errors, cols=1, rows=2
    )





    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="lines",
            name=f"Mean",
            line=dict(color='black'),
            legendrank=1

        )
    )
    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="lines",
            name=f"IQR",
            fill="toself",
            fillcolor='rgba(0,0,0,0.4)',
            line=dict(color='rgba(255,255,255,0)'),
            legendrank=2

        )
    )


    fig.update_yaxes(
        range=[0, 1], row=2, col=1, showticklabels=False,
    )
    for row in range(1, 2):
        fig.update_xaxes(
            showgrid=False, showticklabels=False, row=row
        )
        fig.update_yaxes(
            showgrid=False, showticklabels=False, row=row
        )

    fig.update_yaxes(
        title=dict(text="normalized intensity", font=dict(size=8)),
        row=2,
        col=1,
        showticklabels=True,
        title_standoff=1
    )
    fig.update_xaxes(
        title=dict(text="distance to probe [nm]", font=dict(size=8)),
        row=2,
        col=1,
        showticklabels=True,
        title_standoff=1
    )
    fig.update_layout(**layout)
    for shape in fig.layout.shapes:
        shape['layer'] = 'above'
    for annotation in fig.layout["annotations"]:
        annotation.update(
            layout["annotations"]
        )
    fig.add_annotation(text="# spots<br>per replicate (n)", showarrow=False, xref="paper", yref="paper", x=0.75, y=0,
                       yanchor="top", xanchor="center", font=dict(size=8))

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
        outfile = "foo.svg"
        with open("configFile.yml", 'r') as f:
            cfg = yaml.safe_load(f)
        probes = ["16S rRNA", "cupA", "psaA", "psbA2", "fusB"]

        layout = cfg["F5FISH"]["layout"]
        layout["height"] = 550
        subplot_params = dict(vertical_spacing=0.01, horizontal_spacing=0.02,
                              column_titles=probes,)

    hs = subplot_params["horizontal_spacing"]
    l = (1- (len(probes)-1)*hs) / len(probes)
    fig = main(file, probes, cfg, layout, subplot_params)

    if "snakemake" not in globals():
        fig.show()
    fig.write_image(outfile)
