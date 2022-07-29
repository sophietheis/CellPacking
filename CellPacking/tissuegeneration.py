from tyssue import Sheet
from tyssue import SheetGeometry


def sheet_init (nx, ny):
    sheet = Sheet.planar_sheet_3d(
                'sheet', nx=nx, ny=ny, distx=1, disty=1, noise=0.2)

    # Cut at the border to avoid strange cells
    to_cut = sheet.cut_out([(0.1, nx), (0.1, ny)])
    SheetGeometry.update_all(sheet)
    sheet.remove(to_cut, trim_borders=True)
    sheet.sanitize(trim_borders=True)

    # Center sheet to (0,0)
    SheetGeometry.center(sheet)

    # Add specs for periodic boundary condition
    sheet.update_specs(
                    {'settings':{
                        'boundaries':{'x' : [-nx/2-1, nx/2+1], 
                                      'y': [-ny/2-1, ny/2+1]}
                    }})
    SheetGeometry.update_all(sheet)
    return sheet, SheetGeometry