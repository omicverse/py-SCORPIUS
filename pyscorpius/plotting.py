"""Visualisation helpers — 1:1 ports of SCORPIUS::draw_trajectory_plot
and SCORPIUS::draw_trajectory_heatmap using ggplot2-python.

Render with ``ggsave(path, plot=p, width=..., height=..., dpi=...)``.
"""

from __future__ import annotations

from typing import Iterable

import numpy as np
import pandas as pd

from ggplot2_py import (
    aes,
    coord_equal,
    geom_density_2d,
    geom_path,
    geom_point,
    geom_polygon,
    ggplot,
    labs,
    scale_color_gradientn,
    scale_color_manual,
    scale_fill_gradientn,
    scale_fill_manual,
    stat_density_2d,
    theme_classic,
    xlim,
    ylim,
)


_SET1_BREWER = [
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
    "#FFFF33", "#A65628", "#F781BF", "#999999",
]

# Reverse of RColorBrewer::brewer.pal(9, "RdYlBu")
_RDYLBU_REV = [
    "#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#FFFFBF",
    "#FEE090", "#FDAE61", "#F46D43", "#D73027",
]


def _gg_color_hue(n: int) -> list[str]:
    """ggplot2's default discrete palette — equally-spaced HCL(h, c=100, l=65)."""
    hues = np.linspace(15, 375, n + 1)[:-1]
    # HCL → CIELAB → linear sRGB
    L, C = 65, 100
    out = []
    for h in hues:
        h_rad = np.deg2rad(h)
        a, b = C * np.cos(h_rad), C * np.sin(h_rad)
        # CIELAB → XYZ (D65)
        Y = (L + 16) / 116
        X = a / 500 + Y
        Z = Y - b / 200
        def f_inv(t):
            return t**3 if t**3 > 216/24389 else (116*t - 16) / (24389/27)
        X_xyz = f_inv(X) * 0.95047
        Y_xyz = f_inv(Y) * 1.00000
        Z_xyz = f_inv(Z) * 1.08883
        # XYZ → sRGB
        r =  3.2404542 * X_xyz - 1.5371385 * Y_xyz - 0.4985314 * Z_xyz
        g = -0.9692660 * X_xyz + 1.8760108 * Y_xyz + 0.0415560 * Z_xyz
        b_ = 0.0556434 * X_xyz - 0.2040259 * Y_xyz + 1.0572252 * Z_xyz
        # Gamma-correct + clamp
        def gamma(v):
            v = max(0.0, min(1.0, v))
            return v * 12.92 if v <= 0.0031308 else 1.055 * v**(1/2.4) - 0.055
        rgb = [gamma(r), gamma(g), gamma(b_)]
        out.append("#{:02x}{:02x}{:02x}".format(
            int(round(rgb[0]*255)), int(round(rgb[1]*255)), int(round(rgb[2]*255))))
    return out


def _default_discrete_palette(groups: list) -> list[str]:
    n = len(set(groups))
    if n <= 9:
        return _SET1_BREWER[:n]
    return _gg_color_hue(n)


def _default_continuous_palette() -> list[str]:
    return _RDYLBU_REV


def draw_trajectory_plot(
    space: np.ndarray | pd.DataFrame,
    progression_group=None,
    path: np.ndarray | pd.DataFrame | None = None,
    contour: bool = False,
    progression_group_palette: dict | list | None = None,
    point_size: float = 2.0,
    point_alpha: float = 1.0,
    path_size: float = 0.5,
    path_alpha: float = 1.0,
    contour_alpha: float = 0.2,
):
    """1:1 port of SCORPIUS::draw_trajectory_plot.

    Args:
        space: (n × 2) reduced coordinates.
        progression_group: categorical labels (list/array) OR numeric values.
        path: (m × 2) trajectory points.
        contour: overlay 2-D KDE contour per group.
        progression_group_palette: dict name→colour for discrete groups, or
            list/iterable of colours for numeric.
        ...sizes & alphas as in R.
    """
    arr = np.asarray(space, dtype=np.float64)[:, :2]
    df = pd.DataFrame(arr, columns=["Comp1", "Comp2"])

    is_numeric_group = (
        progression_group is not None
        and pd.api.types.is_numeric_dtype(pd.Series(progression_group))
    )
    if progression_group is not None:
        df["progression_group"] = list(progression_group)

    lo, hi = float(arr.min()), float(arr.max())
    diff = (hi - lo) / 2.0
    if contour:
        lim = (lo - 0.1 * diff, hi + 0.1 * diff)
    else:
        lim = (lo, hi)

    g = (
        ggplot(df)
        + theme_classic()
        + labs(x="Component 1", y="Component 2", colour="Group", fill="Group")
        + xlim(lo - diff, hi + diff)
        + ylim(lo - diff, hi + diff)
        + coord_equal(xlim=lim, ylim=lim)
    )

    if contour and progression_group is not None and not is_numeric_group:
        g = g + stat_density_2d(
            aes(x="Comp1", y="Comp2", fill="progression_group"),
            geom="polygon",
            data=df,
            alpha=contour_alpha,
            contour=True,
        )
    elif contour:
        g = g + stat_density_2d(
            aes(x="Comp1", y="Comp2"),
            data=df,
            alpha=contour_alpha,
            contour=True,
        )

    # Points
    if progression_group is not None:
        g = g + geom_point(
            aes(x="Comp1", y="Comp2", colour="progression_group"),
            data=df,
            size=point_size,
            alpha=point_alpha,
        )
    else:
        g = g + geom_point(
            aes(x="Comp1", y="Comp2"),
            data=df,
            size=point_size,
            alpha=point_alpha,
        )

    # Path
    if path is not None:
        path_arr = np.asarray(path, dtype=np.float64)[:, :2]
        path_df = pd.DataFrame(path_arr, columns=["Comp1", "Comp2"])
        g = g + geom_path(
            aes(x="Comp1", y="Comp2"),
            data=path_df,
            size=path_size,
            alpha=path_alpha,
        )

    # Scale
    if progression_group is not None:
        if is_numeric_group:
            palette = (
                list(progression_group_palette)
                if progression_group_palette is not None
                else _default_continuous_palette()
            )
            g = g + scale_color_gradientn(colours=palette)
        else:
            unique = sorted(set(progression_group), key=str)
            if isinstance(progression_group_palette, dict):
                palette = {k: progression_group_palette[k] for k in unique}
            elif progression_group_palette is not None:
                palette = dict(zip(unique, progression_group_palette))
            else:
                palette = dict(zip(unique, _default_discrete_palette(progression_group)))
            g = g + scale_color_manual(values=palette)
            if contour:
                g = g + scale_fill_manual(values=palette)

    return g


def draw_trajectory_heatmap(
    x: np.ndarray | pd.DataFrame,
    time: np.ndarray,
    progression_group=None,
    modules: pd.DataFrame | None = None,
    show_labels_row: bool = False,
    show_labels_col: bool = False,
    scale_features: bool = True,
    progression_group_palette: dict | None = None,
    title: str | None = None,
    **pheatmap_kwargs,
):
    """1:1 port of SCORPIUS::draw_trajectory_heatmap.

    Uses bio-babel/pheatmap-python (R pheatmap port). Cells are ordered along
    inferred time; rows are features (optionally grouped by modules).

    Returns the pheatmap object — call .save() to write.
    """
    from pheatmap_python import pheatmap

    expr = pd.DataFrame(np.asarray(x, dtype=np.float64))
    if hasattr(x, "columns"):
        expr.columns = list(x.columns)
    else:
        expr.columns = [f"feature_{i}" for i in range(expr.shape[1])]
    if hasattr(x, "index"):
        expr.index = list(x.index)

    if scale_features:
        col_std = expr.std(axis=0).replace(0, 1)
        col_mean = expr.mean(axis=0)
        expr_scaled = (expr - col_mean) / col_std
    else:
        expr_scaled = expr

    # Order cells by time
    order = np.argsort(np.asarray(time))
    mat = expr_scaled.iloc[order].T  # features × cells

    # Build annotation_col
    ann_col = None
    if progression_group is not None:
        ann_col = pd.DataFrame(
            {"Group": np.asarray(progression_group)[order]},
            index=mat.columns,
        )

    # Build annotation_row from modules
    ann_row = None
    if modules is not None:
        if "module" in modules.columns:
            ann_row = modules[["module"]].copy()
        # Re-order rows by module
        if "feature" in modules.columns:
            mat = mat.reindex(modules["feature"].values)

    return pheatmap(
        mat,
        cluster_rows=False,
        cluster_cols=False,
        show_rownames=show_labels_row,
        show_colnames=show_labels_col,
        annotation_col=ann_col,
        annotation_row=ann_row,
        annotation_colors={"Group": progression_group_palette} if progression_group_palette else None,
        main=title,
        **pheatmap_kwargs,
    )
