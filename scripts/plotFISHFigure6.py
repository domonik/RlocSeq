import os
import pytesseract
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import tifffile
import plotly.express as px
import numpy as np
from scripts.plotly_template import register_template, PT7
from scipy import stats
import plotly.io as pio
from scripts.plotlyTrash import replace_square_images


register_template()
pio.renderers.default = "browser"


def _color_to_calpha(color: str, alpha: float = 0.2):
    color = color.split("(")[-1].split(")")[0]
    return f"rgba({color}, {alpha})"

def main(file, probes, cfg, layout, subplot_params, intensities, test_output=None):
    df = pd.read_csv(file, sep="\t")
    probe_names = [probe[0] for probe in probes]

    lr = 0.04 / len(probes)
    r = 0.12 / (len(probes) // 2)
    specs = [
        [{"type": "image", "l": lr, "r": lr}, {"type": "image", "l": lr, "r": lr}, {"type": "xy", "b": 0.05, "t": 0.1, "r": r} ] * (len(probes) // 2),
        [{"type": "xy", "colspan": 1, "r": lr}, {"type": "xy", "colspan": 1, "l": lr}, None] * (len(probes) //2),  # Row 1: two image plots
    ]
    base = 1 / len(probes)

    row_heights =  [0.6, 0.4]
    column_widths = [0.43 * base, 0.43 * base, 0.14 *base] * (len(probes) // 2)

    col = int(len(probes) * 1.5)
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
    col_titles = ["untreated", "Puromycin", None] * (len(probes) // 2)
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
            annotation.update(font=dict(size=PT7))
    else:
        for annotation in fig.layout.annotations:
            annotation.update(font=dict(size=PT7))

    fig.update_layout(cfg["LayoutDefault"])
    images_to_replace = []
    color_map = {"Replicate 1": "rgb(0, 106, 56)", "Replicate 2": "rgb(0, 48, 135)", "Replicate 3": "rgb( 111, 160, 130)"}
    for i, (probe, treatment, replicate) in enumerate(probes):
        traces = []
        errors = []
        image_dir = f"initialData/AllFISHData/{probe}/{replicate}"
        treatl = "untreated" if treatment == "Untreated" else treatment

        files = [f for f in os.listdir(image_dir) if f.endswith(f"{treatl}_green.tif") and "Image" in f]
        cbar_files = [f for f in os.listdir(image_dir) if "CBar" in f]
        if len(cbar_files) > 1:
            cbarfile = [f for f in cbar_files if treatment in f][0]
        else:
            cbarfile = cbar_files[0]

        print(image_dir)
        image_file = os.path.join(image_dir, files[0])
        img = tifffile.imread(image_file)
        fig.add_trace(
            go.Image(z=img),
            row=1,
            col=int(i*1.5) + 1
        )
        images_to_replace.append(image_file)


        image_file = os.path.join(image_dir, cbarfile)
        img = tifffile.imread(image_file)
        max_v = int(pytesseract.image_to_string(img).split("\n")[0])

        img = np.rot90(img)
        mask = (img[..., 0] != 255) &  (img[..., 2] != 255)
        max_color = img[mask, 1].max()

        fake_heat = np.linspace(0, max_v, 100).reshape(100, 1)
        didx = i+1+(len(probes) * 1)
        domain = f"{didx} domain"
        if treatment == "Puromycin":
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
                row=1, col=int(i * 1.5) + 2
            )
            fig.update_yaxes(
                tickvals=[0, max_v // 2, max_v], row=1, col=int(i*1.5)+2,
                range=[0, max_v + 5],  # ensure all tickvals are within the axis
            )
        amaxes = []

        text = "<b>n: </b>"
        subdf = summary[summary["Probe"] == probe]
        subdf = subdf[subdf["Treatment"] == treatment]
        for j, f in enumerate(sorted(subdf["Figure"].unique())):
            color = color_map[f]
            sdf = subdf[subdf["Figure"] == f]
            sdf = sdf[sdf["Mask"] == "Measured"]
            if len(sdf) == 0:
                continue
            if j > 0:
                text += ","
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
                line=dict(width=1, color=color),
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
        print(i)
        fig.add_traces(
            traces + errors, cols=int(i * 1.5)+1, rows=2
        )


        vint = float(np.median(np.asarray(amaxes)))
        didx = i+1 + int(len(probes) * 1.5)
        domain = f"{didx} domain"
        fig.add_annotation(
            x=1, y=1,
            yref=f"y{domain}",
            xref=f"x{domain}",
            yanchor="top",
            xanchor="right",
            text=text,
            showarrow=False,
            borderpad= 0.2,
            bgcolor="white",
            opacity=1,

        )
        fig.add_vline(vint, col=int(i*1.5)+1, row=row, line=dict(dash="dash", color="black", width=1), layer="above")
        fig.add_annotation(
            x=vint,
            yref=f"y{domain}",
            xref=f"x{didx}",
            y=0,  # position near top of the subplot (in data coords)
            text=f"{vint}",  # or any custom label
            showarrow=False,
            xanchor="left",  # position text to the right of the line
            yanchor="bottom",
            font=dict(color="black"),
        )
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
    for col in range(3, int(len(probes)* 1.5) + 1, 3):
        fig.update_yaxes(
            showgrid=False, row=1, col=col, ticks="outside", ticklen=2, showline=False, showticklabels=True, side="right", linewidth=0
        )
        fig.update_xaxes(
            showgrid=False, showticklabels=False, row=1, col=col, ticks="", showline=False
        )
    fig.update_yaxes(
        range=[0, 1], row=2, showticklabels=False, title_standoff=30
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
    return fig, df, images_to_replace







if __name__ == '__main__':
    if "snakemake" in globals():
        file = snakemake.input.tsv
        probes = snakemake.params.probes
        layout = snakemake.params.layout
        cfg = snakemake.config
        outfile = snakemake.output.svg
        subplot_params = snakemake.params.subplot
        html_outfile = snakemake.output.html
        distribution_out = snakemake.output.distribution_out

    else:
        import yaml
        file = "Pipeline/FISH/all_distributions.tsv"
        outfile = "foo.svg"
        with open("configFile.yml", 'r') as f:
            cfg = yaml.safe_load(f)
        probes =[("psbA2","untreated","Replicate 1"), ("psbA2","Puromycin","Replicate 1"), ("rpl10_12", "untreated", "Replicate 1"), ("rpl10_12", "Puromycin", "Replicate 1")]

        layout = cfg["F6FISH"]["layout"]
        subplot_params = dict(vertical_spacing=0., horizontal_spacing=0.,
                              column_titles=["untreated", "Puromycin", None, "untreated", "Puromycin", None])

    hs = subplot_params["horizontal_spacing"]
    l = (1- (len(probes)-1)*hs) / len(probes)
    fig, df, images_to_replace = main(file, probes, cfg, layout, subplot_params, intensities=None)
    df.to_csv(distribution_out, sep="\t", index=False)
    y = 1.1
    fig.add_shape(
        type="line",
        xref="paper",
        yref="paper",
        x0=fig.layout.annotations[0].x,
        x1=fig.layout.annotations[1].x,
        y0=y,
        y1=y,
        line=dict(color="black", width=1),
        opacity=1
    )
    fig.add_shape(
        type="line",
        xref="paper",
        yref="paper",
        x0=fig.layout.annotations[2].x,
        x1=fig.layout.annotations[3].x,
        y0=y,
        y1=y,
        line=dict(color="black", width=1),
        opacity=1
    )
    fig.add_annotation(
        xref="paper",
        yref="paper",
        x=(fig.layout.annotations[1].x + fig.layout.annotations[0].x) /2,
        xanchor="center",
        y=y,
        yanchor="bottom",
        text=f"<i>{probes[0][0] if probes[0][0] != "psbA2" else "psbA"}</i>",
        showarrow=False,
    )
    fig.add_annotation(
        xref="paper",
        yref="paper",
        x=(fig.layout.annotations[3].x + fig.layout.annotations[2].x) /2,
        xanchor="center",
        y=y,
        yanchor="bottom",
        text="<i>rpl10, rpl12</i>",
        showarrow=False,
    )
    if "snakemake" not in globals():
        fig.show()

    for annotation in fig.layout.annotations:
        if annotation.text.startswith("Distance to"):
            annotation.update(yanchor="bottom")

    tmp_out = outfile + ".tmp"
    fig.write_image(tmp_out, format="svg")
    replace_square_images(tmp_out, images_to_replace, outfile)
    fig.write_html(html_outfile)

