from svgutils.compose import Figure, Panel, SVG, Text


def main(input, config, output):
    b_y = config["RBP3Venn"]["layout"]["height"]
    ges_y = config["RBP3Venn"]["layout"]["height"] + config["RBP3GSEA"]["layout"]["height"]
    f = Figure("624px", f"{ges_y}px",
               Panel(
                   SVG(input.f1),
                   Text("A", 2, 15, size=config["fonts"]["multipanel"]["size"], weight='bold', font="Arial")

               ),
               Panel(
                   SVG(input.f2),
                   Text("B", 2, 2, size=config["fonts"]["multipanel"]["size"], weight='bold', font="Arial")
               ).move(0, b_y),
               )
    svg_string = f.tostr()
    svg_string = svg_string.decode().replace("encoding='ASCII'", "encoding='utf-8'")
    with open(output.svg, "w") as handle:
        handle.write(svg_string)


if __name__ == '__main__':
    main(
        snakemake.input,
        snakemake.config,
        snakemake.output,

    )