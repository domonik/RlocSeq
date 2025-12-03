import os
import pytesseract
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import tifffile
import plotly.express as px
import numpy as np
from scripts.plotly_template import register_template, PT8, PT7
from scipy import stats
import plotly.io as pio


register_template()
pio.renderers.default = "browser"


def _color_to_calpha(color: str, alpha: float = 0.2):
    color = color.split("(")[-1].split(")")[0]
    return f"rgba({color}, {alpha})"

def main(file, probes, cfg, layout, subplot_params, intensities=None, test_output=None):
    df = pd.read_csv(file, sep="\t")
    if intensities is not None:
        intensities = pd.read_csv(intensities, sep="\t")
        probe_names = [probe[0] for probe in probes]
        intensities, results_df = calc_summary_and_pval(intensities, probe_names)

        lr = 0.04 / len(probes)
        r = 0.15 / len(probes)
        specs = [
            [{"type": "image", "l": lr, "r": lr}, {"type": "image", "l": lr, "r": lr}, {"type": "xy", "b": 0.01, "t": 0.05, "r": r} ] * (len(probes) // 2),
            [{"type": "box", "colspan": 3, "r": r+ lr, "b": 0.1}, None, None] * (len(probes) // 2),
            [{"type": "xy", "colspan": 1, "r": r+ lr}, {"type": "xy", "colspan": 1, "r": r+ lr}, None] * (len(probes) //2),  # Row 1: two image plots
        ]
        base = 1 / (len(probes) // 2)

        row_heights =  [0.4, 0.3, 0.3]
        column_widths = [base * .4, base*.4, base * .2] * (len(probes) // 2)

        col = int(len(probes) * 1.5)
        row = 3
    else:
        lr = 0.05 / len(probes)
        r = 0.15 / len(probes)
        col = len(probes)*2
        specs = [
            [{"type": "image", "l": lr, "r": lr}, {"type": "xy", "b": 0.12, "t": 0.2, "r": r} ] * len(probes),
            [{"type": "xy", "colspan": 2, "r": r+ lr}, None] * len(probes),  # Row 1: two image plots
        ]
        row_heights = [0.7, 0.3]
        base = 1 / len(probes)

        column_widths = [base * .8, base * .2] * len(probes)
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
    col_titles = []
    for title in subplot_params["column_titles"]:
        col_titles.append(title)
        col_titles.append(None)
    del subplot_params["column_titles"]
    fig = make_subplots(
        rows=row,
        cols=col,
        row_heights=row_heights,
        x_title="Distance to probe [nm]",
        specs=specs,
        column_widths=column_widths,
        column_titles = col_titles,
        **subplot_params


    )
    if intensities is None:
        for annotation in fig.layout.annotations:
            annotation.update(font=dict(size=PT7), yanchor="top")
    else:
        for annotation in fig.layout.annotations:
            annotation.update(font=dict(size=PT7))

    fig.update_layout(cfg["LayoutDefault"])

    color_map = {"Replicate 1": "rgb(0, 106, 56)", "Replicate 2": "rgb(0, 48, 135)", "Replicate 3": "rgb( 111, 160, 130)"}
    for i, (probe, treatment, replicate) in enumerate(probes):
        traces = []
        errors = []
        image_dir = f"initialData/AllFISHData/{probe}/{replicate}"
        treatl = "untreated" if treatment == "Untreated" else treatment
        files = [f for f in os.listdir(image_dir) if f.endswith(f"{treatl}_green.tif") and "Image" in f]
        cbar_files = [f for f in os.listdir(image_dir) if "CBar" in f]
        print(image_dir)
        if len(cbar_files) > 1:
            cbarfile = [f for f in cbar_files if treatl in f][0]
        else:
            cbarfile = cbar_files[0]

    
        image_file = os.path.join(image_dir, files[0])
        img = tifffile.imread(image_file)
        tmp = px.imshow(img)
        fig.add_trace(tmp.data[0], row=1, col=(i * 2) + 1)

        image_file = os.path.join(image_dir, cbarfile)
        img = tifffile.imread(image_file)
        max_v = int(pytesseract.image_to_string(img).split("\n")[0])

        img = np.rot90(img)
        mask = (img[..., 0] != 255) &  (img[..., 2] != 255)
        max_color = img[mask, 1].max()

        fake_heat = np.linspace(0, max_v, 100).reshape(100, 1)
        didx = i+1+(len(probes) * 1) if intensities is None else i+1 + (len(probes) * 2) // 2
        domain = f"{didx} domain"
        if treatment == "Puromycin" or intensities is None:
            fig.add_trace(
                go.Heatmap(
                    z=fake_heat,
                    colorscale=["rgb(0, 0, 0)", f"rgb(0, {max_color}, 0)"],
                    showscale=False,
                    # Hide the actual heatmap by making it tiny and offset
                    y=fake_heat.flatten(),  # off the visible range
                    x=[-1],
                    colorbar=dict(
                        len=0.2,
                        x=(1 / len(probes)) * i,
                        y=1 - row_heights[0],
                        orientation="h",
                        yref=f"paper",
                        xref=f"paper",
                        xanchor="left",
                        yanchor="middle"
                    )
                    # Configure colorbar
                ),
                row=1, col=i * 2 + 2
            )
            fig.update_yaxes(
                tickvals=[0, max_v // 2, max_v], row=1, col=i*2+2,
                range=[0, max_v],  # ensure all tickvals are within the axis
            )
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
            traces + errors, cols=(i * 2)+1, rows=2 if intensities is None else 3
        )
        if intensities is not None:
            group = intensities[(intensities.Treatment == treatment) & (intensities.Probe == probe)]
            values = group["subtractedIntensity"].tolist()
            mean_val = np.mean(values)
            std_val = np.std(values)
            fig.add_trace(
                go.Bar(
                    x=[f"{treatment} – {probe}"],
                    y=[mean_val],
                    name=f"{treatment} – <i>{probe}</i>",
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
                col=1 if probe == probes[0][0] else 5
            )
        vint = float(np.median(np.asarray(amaxes)))
        didx = i+1+(len(probes) * 2) if intensities is None else i+1 + (len(probes) * 2) + len(probes) // 2
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
        fig.add_vline(vint, col=i*2+1, row=row, line=dict(dash="dash", color="black", width=1), layer="above")
        fig.add_annotation(
            x=vint,
            yref=f"y{domain}",
            xref=f"x{didx}",
            y=0,  # position near top of the subplot (in data coords)
            text=f"{vint}",  # or any custom label
            showarrow=False,
            xanchor="left",  # position text to the right of the line
            yanchor="bottom",
            font=dict(color="black", size=7),
        )
    colors = list(color_map.values())
    # fig.add_traces(
    #     [
    #         go.Scatter(
    #             x=[None],
    #             y=[None],
    #             mode="markers",
    #             name=f"{i+1}",
    #             marker=dict(color=colors[i], symbol="square", size=10),
    #             legendgroup="rep",
    #             legendgrouptitle=dict(text="Replicate", font=dict(size=10)),
    #
    #             legendrank=2
    #
    #         ) for i in range(len(colors))
    #     ]
    # )
    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            legendgroup="Data",
            legendgrouptitle=dict(text=""),
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
        fig.update_yaxes(
            title=dict(text="Mean Fluorescence<br>Intensity [A.U.]"), row=3, col=1
        )
        results_df.to_csv(test_output, sep="\t", index=False)


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
        showgrid=False, showticklabels=False, row=1, ticks="", showline=False
    )

    fig.update_yaxes(
        showgrid=False, showticklabels=False, row=1, ticks="", showline=False
    )
    for col in range(2, len(probes)* 2 + 1, 2):
        fig.update_yaxes(
            showgrid=False, row=1, col=col, ticks="outside", ticklen=3, showline=False, showticklabels=True, side="right"
        )
        fig.update_xaxes(
            showgrid=False, showticklabels=False, row=1, col=col, ticks="", showline=False
        )
    fig.update_yaxes(
        range=[0, 1], row=row, showticklabels=False, title_standoff=30
    )

    fig.update_layout(
        width=621,
        height=400,
    )

    fig.update_yaxes(
        title=dict(text="Norm. intensity"),
        row=row,
        col=1,
        showticklabels=True,
    )
    fig.update_layout(bargap=0.7)
    fig.update_layout(**layout)
    for shape in fig.layout.shapes:
        shape['layer'] = 'above'
    df = df.rename({"Figure": "Replicate"}, axis=1)
    return fig, df








if __name__ == '__main__':
    if "snakemake" in globals():
        file = snakemake.input.tsv
        probes = snakemake.params.probes
        layout = snakemake.params.layout
        cfg = snakemake.config
        outfile = snakemake.output.svg
        subplot_params = snakemake.params.subplot
        distribution_out = snakemake.output.distribution_out
    else:
        import yaml
        file = "Pipeline/FISH/all_distributions.tsv"
        probes = [("psaA", "untreated", "Replicate 1"), ("cupA", "untreated", "Replicate 1"), ("rnpB", "untreated", "Replicate 1"), ("fusB", "untreated", "Replicate 1") ]
        outfile = "foo.svg"
        subplot_params = dict(vertical_spacing=0.0, column_titles=[p[0] for p in probes], horizontal_spacing=0.)
        with open("configFile.yml", 'r') as f:
            cfg = yaml.safe_load(f)
        layout = cfg["F2FISH"]["layout"]
    fig, df = main(file, probes, cfg, layout, subplot_params)
    df.to_csv(distribution_out, sep="\t", index=False)
    for annotation in fig.layout.annotations:
        if annotation.text.startswith("Distance to"):
            annotation.update(yanchor="bottom")
    if "snakemake" not in globals():
        fig.show()
    fig.write_image(outfile)
