import plotly.graph_objects as go
import yaml

def main(cfg, outfile):
    x = ["Control", "Puromycin"]
    y = [156.575, 169.039]
    std = [6.52287, 7.00263]
    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=x,
            y=y,
            error_y=dict(
                type='data',
                visible=True,
                array=std
            ),
            marker=dict(line=dict(color="black", width=1), color=cfg["cytosolColor"])
        )
    )
    fig.update_yaxes(title=dict(text="2<sup>–ΔΔCt</sup>"))
    fig.update_layout(cfg["LayoutDefault"])
    fig.update_layout(cfg["F3QTPCR"]["layout"])

    fig.write_image(outfile)



if __name__ == '__main__':
    with open("configFile.yml", 'r') as f:
        cfg = yaml.safe_load(f)
    cfg = snakemake.config
    outfile = snakemake.output.svg
    main(cfg, outfile)