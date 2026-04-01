from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd
import plotly.express as px


def plotly_spatial_cells(
    adata,
    image_id: str = "patient_exp",
    x: str = "Cell Center X",
    y: str = "Cell Center Y",
    phenotype: str = "cell_type",
    size: float = 2,
    color_map: dict[str, str] | None = None,
    save_html: str | Path | None = None,
    show_browser: bool = True,
    **kwargs: Any,
):
    """
    Plot spatial cells from an AnnData object using Plotly.

    Parameters
    ----------
    adata
        AnnData object containing metadata in `adata.obs`.
    image_id
        Column in `adata.obs` identifying the sample/patient/image.
    x
        Column in `adata.obs` containing x coordinates.
    y
        Column in `adata.obs` containing y coordinates.
    phenotype
        Column in `adata.obs` containing cell labels/classes.
    size
        Marker size.
    color_map
        Optional mapping from phenotype labels to colors.
    save_html
        Optional path to save the interactive figure as HTML.
    show
        Whether to display the figure.
    **kwargs
        Additional keyword arguments passed to `plotly.express.scatter`.

    Returns
    -------
    plotly.graph_objects.Figure
        The generated Plotly figure.
    """
    required_cols = [image_id, x, y, phenotype]
    missing = [col for col in required_cols if col not in adata.obs.columns]
    if missing:
        raise ValueError(f"Missing columns in adata.obs: {missing}")

    patient_values = adata.obs[image_id].dropna().unique()
    patient_exp = patient_values[0] if len(patient_values) > 0 else "unknown"

    data = pd.DataFrame(
        {
            "x": adata.obs[x],
            "y": adata.obs[y],
            "col": adata.obs[phenotype],
        }
    ).sort_values(by="col")

    scatter_kwargs = dict(x="x", y="y", color="col")
    if color_map is not None:
        scatter_kwargs["color_discrete_map"] = color_map

    fig = px.scatter(data, **scatter_kwargs, **kwargs)

    fig.update_traces(
        marker=dict(size=size),
        selector=dict(mode="markers"),
        hoverlabel=dict(namelength=-1),
    )

    fig.update_yaxes(autorange="reversed", tickformat="g")
    fig.update_xaxes(tickformat="g")
    fig.update_layout(
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        title_text=f"Patient: {patient_exp}",
        title_x=0.5,
        title_font_size=20,
    )

    if save_html is not None:
        fig.write_html(str(save_html), auto_open=show)

    if show_browser:
        fig.show(renderer="browser")

    return fig