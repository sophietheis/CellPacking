import numpy as np
from tyssue import Sheet
# from tyssue import PlanarGeometry as geom
from .dynamics import ShearPlanarGeometry as geom


def sheet_init(nx, ny, gamma_0=0.5, phi=np.pi / 3, noise=0.2):
    sheet = Sheet.planar_sheet_2d(
        'sheet', nx=nx, ny=ny, distx=1, disty=1, noise=noise)

    # Cut at the border to avoid strange cells
    to_cut = sheet.cut_out([(0.1, nx), (0.1, ny)])

    sheet.update_specs({"edge": {"gamma_0": gamma_0,
                                 "phi": phi},
                        })

    geom.update_all(sheet)
    sheet.remove(to_cut, trim_borders=True)
    sheet.sanitize(trim_borders=True)
    geom.update_all(sheet)

    # Center sheet to (0,0)
    geom.center(sheet)

    # Add specs for periodic boundary condition
    sheet.update_specs(
        {'settings': {
            'boundaries': {'x': [-nx / 2 - 1, nx / 2 + 1],
                           'y': [-ny / 2 - 1, ny / 2 + 1]}
        }})
    geom.update_all(sheet)
    sheet.edge_df["opposite"] = sheet.get_opposite()
    return sheet, geom
