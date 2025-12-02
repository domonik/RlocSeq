from scripts.plotFISHWithIntensities import main





if __name__ == '__main__':
    if "snakemake" in globals():
        file = snakemake.input.tsv
        probes = snakemake.params.probes
        layout = snakemake.params.layout
        cfg = snakemake.config
        outfile = snakemake.output.svg
        subplot_params = snakemake.params.subplot
        intensities = snakemake.input.intensities
        res_tsv = snakemake.output.res_tsv
        html_outfile = snakemake.output.html
        distribution_out = snakemake.output.distribution_out


    else:
        import yaml
        file = "Pipeline/FISH/all_distributions.tsv"
        intensities = "Pipeline/FISH/single_spots.tsv"
        outfile = "foo.svg"
        distribution_out = "foo.tsv"
        with open("configFile.yml", 'r') as f:
            cfg = yaml.safe_load(f)
        probes =[("fusB","untreated","Replicate 1"), ("fusB","Puromycin","Replicate 1"), ("psaA", "untreated", "Replicate 1"), ("psaA", "Puromycin", "Replicate 1")]

        layout = cfg["F5FISH"]["layout"]
        subplot_params = dict(vertical_spacing=0., horizontal_spacing=0.,
                              column_titles=["untreated", "Puromycin", "untreated", "Puromycin"])
        res_tsv = None

    hs = subplot_params["horizontal_spacing"]
    l = (1- (len(probes)-1)*hs) / len(probes)
    fig, df = main(file, probes, cfg, layout, subplot_params, intensities=intensities, test_output=res_tsv)
    df.to_csv(distribution_out, sep="\t", index=False)
    y = 1.075
    fig.add_shape(
        type="line",
        xref="x7",
        yref="paper",
        x0=0,
        x1=1,
        y0=y,
        y1=y,
        line=dict(color="black", width=1),
        opacity=1
    )
    fig.add_shape(
        type="line",
        xref="x8",
        yref="paper",
        x0=0,
        x1=1,
        y0=y,
        y1=y,
        line=dict(color="black", width=1),
        opacity=1
    )
    fig.add_annotation(
        xref="x7",
        yref="paper",
        x=0.5,
        xanchor="center",
        y=y,
        yanchor="bottom",
        text=probes[0][0],
        showarrow=False,
    )
    fig.add_annotation(
        xref="x8",
        yref="paper",
        x=0.5,
        xanchor="center",
        y=y,
        yanchor="bottom",
        text=probes[2][0],
        showarrow=False,
    )
    if "snakemake" not in globals():
        fig.show()
    for annotation in fig.layout.annotations:
        if annotation.text.startswith("Distance to"):
            annotation.update(yanchor="bottom")
    fig.write_image(outfile)
    fig.write_html(html_outfile)
