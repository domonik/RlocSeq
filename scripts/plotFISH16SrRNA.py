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
        distribution_out = snakemake.output.distribution_out
    else:
        import yaml
        file = "Pipeline/FISH/all_distributions.tsv"
        intensities = "Pipeline/FISH/single_spots.tsv"

        outfile = "foo.svg"
        with open("configFile.yml", 'r') as f:
            cfg = yaml.safe_load(f)
        probes =[("16S rRNA","untreated","Replicate 1"), ("16S rRNA","Puromycin","Replicate 1")]

        layout = cfg["Supplement16SFISH"]["layout"]
        subplot_params = dict(vertical_spacing=0.05, horizontal_spacing=0.05,
                              column_titles=["- Puromycin", "+ Puromycin", "- Puromycin", "+ Puromycin"])

    hs = subplot_params["horizontal_spacing"]
    l = (1- (len(probes)-1)*hs) / len(probes)
    fig, df = main(file, probes, cfg, layout, subplot_params, intensities=intensities)
    df.to_csv(distribution_out, sep="\t", index=False)


    if "snakemake" not in globals():
        fig.show()
    fig.write_image(outfile)
