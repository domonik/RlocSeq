import plotly.graph_objects as go
import plotly.io as pio

PT7 = 9.33333333333
PT8 = 10.666666666

def register_template():
    base = pio.templates["plotly_white"]

    pio.templates["simple_white_custom"] = go.layout.Template(
        layout={
            **base.layout.to_plotly_json(),
            "xaxis": {
                "title_font": {"size": PT8},
                "tickfont": {"size": PT7},
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
                "title_font": {"size": PT8},
                "tickfont": {"size": PT7},
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
                "font": {"size": PT7},
                "title": {"font": {"size": PT7}},
            },
            "legend2": {
                "font": {"size": PT7},
                "title": {"font": {"size": PT7}},
            },
            "font": {"family": "Arial", "size": PT7, "color": "black"},
            "margin": {"b": 50, "t": 10, "l":30 , "r": 30   },
            "annotationdefaults": {"font": {"family": "Arial", "size": PT7},},
            "coloraxis": {
                "colorbar": {
                    "title": {"font": {"size": PT7}},
                    "tickfont": {"size": PT7},
                }

            },
        }

    )