import plotly.graph_objects as go
import plotly.io as pio


def register_template():
    base = pio.templates["plotly_white"]

    pio.templates["simple_white_custom"] = go.layout.Template(
        layout={
            **base.layout.to_plotly_json(),
            "xaxis": {
                "title_font": {"size": 8},
                "tickfont": {"size": 7},
                "color": "black",
                "title_standoff": 5,
                "showgrid": False,
                "zeroline": False,
                "showline": True,
                "linecolor": 'black',
                "linewidth": 1,
                "ticks": "outside",
                "tickwidth": 1

            },
            "yaxis": {
                "title_font": {"size": 8},
                "tickfont": {"size": 7},
                "color": "black",
                "title_standoff": 5,
                "showgrid": False,
                "zeroline": False,
                "showline": True,
                "linecolor": 'black',
                "linewidth": 1,
                "ticks": "outside",
                "tickwidth": 1

            },
            "legend": {
                "font": {"size": 7},
                "title": {"font": {"size": 7}},
            },
            "legend2": {
                "font": {"size": 7},
                "title": {"font": {"size": 7}},
            },
            "font": {"family": "Arial", "size": 8, "color": "black"},
            "margin": {"b": 50, "t": 10, "l":30 , "r": 30   },
            "annotationdefaults": {"font": {"family": "Arial", "size": 8},},
            "coloraxis": {
                "colorbar": {
                    "title": {"font": {"size": 8}},
                    "tickfont": {"size": 7},
                }

            },
        }

    )