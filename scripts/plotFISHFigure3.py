from scripts.plotFISH import main





if __name__ == '__main__':
    if "snakemake" in globals():
        file = snakemake.input.tsv
        probes = snakemake.params.probes
        layout = snakemake.params.layout
        cfg = snakemake.config
        outfile = snakemake.output.svg
        html_outfile = snakemake.output.html
        subplot_params = snakemake.params.subplot
        distribution_out = snakemake.output.distribution_out
    else:
        import yaml
        file = "Pipeline/FISH/all_distributions.tsv"
        outfile = "foo.svg"
        raise ValueError
        with open("configFile.yml", 'r') as f:
            cfg = yaml.safe_load(f)
        probes = [("16S rRNA", "untreated", "Replicate 1"), ("16S rRNA", "Puromycin", "Replicate 1"),]

        layout = cfg["F5FISH"]["layout"]
        subplot_params = dict(vertical_spacing=0.05, horizontal_spacing=0.05,
                              column_titles=["untreated", "Puromycin"])

    hs = subplot_params["horizontal_spacing"]
    l = (1- (len(probes)-1)*hs) / len(probes)
    fig, df = main(file, probes, cfg, layout, subplot_params)
    df.to_csv(distribution_out, sep="\t", index=False)
    for annotation in fig.layout.annotations:
        if annotation.text.startswith("distance to"):
            annotation.update(yanchor="bottom", text="")
    fig.update_xaxes(title=dict(text="distance to probe [nm]"), row=2, col=1)
    if "snakemake" not in globals():
        print("showing")
        fig.show()
    fig.write_image(outfile)
    fig.write_html(html_outfile)
