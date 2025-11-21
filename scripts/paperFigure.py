from svgutils.compose import Figure, Panel, SVG, Text


def main(params, output, config):
    ges_y = params.ges_y
    figures = params.figures
    panels = []
    letters = []
    for letter, value in figures.items():
        fconfig = config[value]
        x = fconfig["x"]
        y = fconfig["y"]
        p = Panel(
            SVG(fconfig["file"]),
        ).move(x, y)
        letters.append(Text(letter, 2, 15, size=config["fonts"]["multipanel"]["size"], weight='bold', font="Arial").move(x, y))
        panels.append(p)

    panels = panels + letters
    f = Figure("680px", f"{ges_y}px", *panels)
    svg_string = f.tostr()
    svg_string = svg_string.decode().replace("encoding='ASCII'", "encoding='utf-8'")
    with open(output.svg, "w") as handle:
        handle.write(svg_string)


if __name__ == '__main__':
    main(
        snakemake.params,
        snakemake.output,
        snakemake.config

    )