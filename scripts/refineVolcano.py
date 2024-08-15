


def main(json, output, params, config):
    lout = config[params["figureConfig"]]
    from plotly.io import read_json
    fig = read_json(json)
    fig.update_layout(config["LayoutDefault"])
    fig.update_layout(lout["layout"])
    for idx, shape in enumerate(fig.layout.shapes):
        if idx == 0:
            shape.update(fillcolor=config["membraneColor"], opacity=0.2)
        elif idx == 1:
            shape.update(fillcolor=config["cytosolColor"], opacity=0.2)
    fig.layout.annotations = None

    fig.update_traces(
        marker=dict(size=3),
    )
    fig.add_annotation(
        text=params["baseline"],
        font=dict(color=config["cytosolColor"]),
        showarrow=False,
        x=0,
        y=1,
        xref="x domain",
        yref="y domain",
        xanchor="left",
        yanchor="top"
    )
    fig.add_annotation(
        text=params["condition"],
        font=dict(color=config["membraneColor"]),
        showarrow=False,
        x=1,
        y=1,
        xref="x domain",
        yref="y domain",
        xanchor="right",
        yanchor="top"
    )
    fig.data[0].update(name="mRNA", showlegend=True)
    fig.update_xaxes(dtick=2)
    fig.write_image(output.svg)

if __name__ == '__main__':
    main(
        snakemake.input.json,
        snakemake.output,
        snakemake.params,
        snakemake.config
    )