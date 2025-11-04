from importlib_metadata import distribution

from scripts.plotFISH import main





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
        subplot_params = dict(vertical_spacing=0.05, horizontal_spacing=0.05,
                              column_titles=["untreated", "Puromycin", "untreated", "Puromycin"])

    hs = subplot_params["horizontal_spacing"]
    l = (1- (len(probes)-1)*hs) / len(probes)
    fig, df = main(file, probes, cfg, layout, subplot_params, intensities=None)
    df.to_csv(distribution_out, sep="\t", index=False)
    y = 1.1
    fig.add_shape(
        type="line",
        xref="paper",
        yref="paper",
        x0=l/4,
        x1=l*1.75+hs,
        y0=y,
        y1=y,
        line=dict(color="black", width=2),
        opacity=1
    )
    fig.add_shape(
        type="line",
        xref="paper",
        yref="paper",
        x0=2.25*l+2*hs,
        x1=3.75*l+3*hs,
        y0=y,
        y1=y,
        line=dict(color="black", width=2),
        opacity=1
    )
    fig.add_annotation(
        xref="paper",
        yref="paper",
        x=1/4,
        xanchor="center",
        y=y,
        yanchor="bottom",
        text=probes[0][0],
        showarrow=False,
    )
    fig.add_annotation(
        xref="paper",
        yref="paper",
        x=3 / 4,
        xanchor="center",
        y=y,
        yanchor="bottom",
        text=probes[2][0],
        showarrow=False,
    )
    if "snakemake" not in globals():
        fig.show()
    fig.write_image(outfile)
    fig.write_html(html_outfile)

