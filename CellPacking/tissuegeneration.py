import numpy as np
import pandas as pd

from tyssue import Sheet
# from tyssue import PlanarGeometry as geom
from .dynamics import ShearPlanarGeometry as geom
from tyssue.generation import config, AnnularSheet


def sheet_init(nx, ny, gamma_0=0.5, phi=np.pi / 3, noise=0.2):
    sheet = Sheet.planar_sheet_2d(
        'sheet', nx=nx, ny=ny, distx=1, disty=1, noise=noise)

    # Cut at the border to avoid strange cells
    to_cut = sheet.cut_out([(0.1, nx), (0.1, ny)])

    sheet.update_specs({"edge": {"gamma_0": gamma_0,
                                 "phi0": phi},
                        })

    geom.update_all(sheet)
    sheet.remove(to_cut, trim_borders=True)
    sheet.sanitize(trim_borders=True)
    geom.update_all(sheet)

    # Center sheet to (0,0)
    geom.center(sheet)

    # Add specs for periodic boundary condition
    # sheet.update_specs(
    #     {'settings': {
    #         'boundaries': {'x': [-nx / 2 - 1, nx / 2 + 1],
    #                        'y': [-ny / 2 - 1, ny / 2 + 1]}
    #     }})
    geom.update_all(sheet)
    sheet.edge_df["opposite"] = sheet.get_opposite()
    return sheet, geom


def symetric_circular(radius, gamma_0=0.5, phi_apical=np.pi / 2, phi_basal=0, noise=0.0):
    ncells = radius * 3
    sheet = Sheet.planar_sheet_2d(
        "planar", nx=ncells, ny=ncells, distx=1, disty=1, noise=noise,
    )

    sheet.specs['settings']['dt'] = 0.01
    sheet.edge_df['density'] = 0.
    sheet.edge_df['weight'] = 1

    sheet.edge_df['line_tension'] = sheet.edge_df['weight'] * 1
    sheet.edge_df['is_active'] = 1
    sheet.face_df['prefered_area'] = 1
    sheet.face_df['prefered_perimeter'] = 3. * np.sqrt(sheet.face_df['prefered_area'])
    sheet.face_df['perimeter_elasticity'] = 1.
    sheet.face_df['area_elasticity'] = 1.
    sheet.vert_df['viscosity'] = 1.0
    sheet.update_specs({"vert": {"compression": 0.0}})  # 0.01
    sheet.update_specs({"edge": {"gamma_0": gamma_0,
                                 "phi0_apical": phi_apical,
                                 "phi0_basal": phi_basal},
                        })
    sheet.edge_df['strain'] = 0.

    sheet.get_opposite()
    sheet.edge_df.loc[sheet.edge_df.opposite == -1, 'line_tension'] = 5.0

    geom.update_all(sheet)

    geom.scale(sheet, sheet.face_df.area.median() ** (-0.5), sheet.coords)
    geom.center(sheet)
    geom.update_all(sheet)
    # Put the center most cell at the origin
    central_cell = (sheet.face_df.x ** 2 + sheet.face_df.y ** 2).idxmin()
    sheet.vert_df[["x", "y"]] -= sheet.face_df.loc[central_cell, ["x", "y"]]
    geom.update_all(sheet)

    out = sheet.edge_df[
        (sheet.edge_df.fx ** 2 + sheet.edge_df.fy ** 2) > radius ** 2
        ].index
    sheet.remove(out)
    geom.update_all(sheet)
    geom.center(sheet)
    sheet.sanitize(trim_borders=True)
    geom.update_all(sheet)

    sheet.edge_df["opposite"] = sheet.get_opposite()
    border_edges = sheet.edge_df[sheet.edge_df["opposite"] == -1].index
    return sheet, border_edges


def generate_ellipsis(Nf, a, b, cell_height, R_vit=None, apical="in"):
    """Generates a 2D tyssue object aranged in a ring of Nf tetragonal cells
    with inner diameter R_in and outer diameter R_out

    Parameters
    ----------
    Nf : int
        The number of cells in the tissue
    a : float
        semi major axis
    b : float
        semi_minor axis
    cell_height : float
        height of cells
    R_vit : float
        The vitelline membrane diameter
        (a non strechable membrane around the annulus)
    apical : str {'in' | 'out'}
        The side of the apical surface
        if "in", the apical surface is inside the annulus, facing
        the lumen as in an organoid; if 'out': the apical side is facing the
        exterior of the tissue, as in an embryo
    Returns
    -------
    eptm : :class:`AnnularSheet`
        2D annular tissue. The `R_in` and `R_out` parameters
        are stored in the class `settings` attribute.
    """
    specs = config.geometry.planar_spec()
    specs["settings"] = specs.get("settings", {})
    specs["settings"]["a"] = a
    specs["settings"]["b"] = b
    specs["settings"]["cell_height"] = cell_height

    Ne = Nf * 4
    Nv = Nf * 2
    vert_df = pd.DataFrame(
        index=pd.Index(range(Nv), name="vert"),
        columns=specs["vert"].keys(),
        dtype=float,
    )
    edge_df = pd.DataFrame(
        index=pd.Index(range(Ne), name="edge"),
        columns=specs["edge"].keys(),
        dtype=float,
    )
    face_df = pd.DataFrame(
        index=pd.Index(range(Nf), name="face"),
        columns=specs["face"].keys(),
        dtype=float,
    )

    inner_edges = np.array(
        [
            [f0, v0, v1]
            for f0, v0, v1 in zip(range(Nf), range(Nf), np.roll(range(Nf), -1))
        ]
    )

    outer_edges = np.zeros_like(inner_edges)
    outer_edges[:, 0] = inner_edges[:, 0]
    outer_edges[:, 1] = inner_edges[:, 2] + Nf
    outer_edges[:, 2] = inner_edges[:, 1] + Nf

    left_spokes = np.zeros_like(inner_edges)
    left_spokes[:, 0] = inner_edges[:, 0]
    left_spokes[:, 1] = outer_edges[:, 2]
    left_spokes[:, 2] = inner_edges[:, 1]

    right_spokes = np.zeros_like(inner_edges)
    right_spokes[:, 0] = inner_edges[:, 0]
    right_spokes[:, 1] = inner_edges[:, 2]
    right_spokes[:, 2] = outer_edges[:, 1]

    edges = np.concatenate([inner_edges, outer_edges, left_spokes, right_spokes])

    edge_df[["face", "srce", "trgt"]] = edges
    edge_df[["face", "srce", "trgt"]] = edge_df[["face", "srce", "trgt"]].astype(int)

    thetas = np.linspace(0, 2 * np.pi, Nf, endpoint=False)
    thetas += thetas[1] / 2

    thetas = thetas[::-1]
    # Setting vertices position (turning clockwise for correct orientation)
    vert_df.loc[range(Nf), "x"] = a * np.cos(thetas)
    vert_df.loc[range(Nf), "y"] = b * np.sin(thetas)
    vert_df.loc[range(Nf, 2 * Nf), "x"] = (a + cell_height) * np.cos(thetas)
    vert_df.loc[range(Nf, 2 * Nf), "y"] = (b + cell_height) * np.sin(thetas)

    vert_df["segment"] = "basal"
    edge_df["segment"] = "basal"
    if apical == "out":
        edge_df.loc[range(Nf, 2 * Nf), "segment"] = "apical"
        vert_df.loc[range(Nf, 2 * Nf), "segment"] = "apical"
    elif apical == "in":
        edge_df.loc[range(Nf), "segment"] = "apical"
        vert_df.loc[range(Nf), "segment"] = "apical"
    else:
        raise ValueError(
            "apical argument not understood, "
            f"should be either 'in' or 'out', got {apical}"
        )
    edge_df.loc[range(2 * Nf, 4 * Nf), "segment"] = "lateral"

    datasets = {"vert": vert_df, "edge": edge_df, "face": face_df}
    ring = AnnularSheet("Ellipsis", datasets, specs, coords=["x", "y"])

    ring.reset_topo()
    return ring







