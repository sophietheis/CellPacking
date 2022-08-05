import pandas as pd
import numpy as np
import plotly.express as px

from tyssue.draw.plt_draw import _get_lines
from tyssue.draw import highlight_cells
from tyssue.draw.ipv_draw import sheet_view as sheet_view_3d


def sheet_view(sheet, name_sheet=None):
    if name_sheet is None:
        name_sheet = ["sheet"]

    a_lines_x, a_lines_y = _get_lines(sheet, list('xy'))

    df = pd.DataFrame([a_lines_x,
                       a_lines_y,
                       np.repeat(name_sheet[0], len(a_lines_x))])
    df = df.T
    df.columns = ['lines_x', 'lines_y', 'name']

    fig = px.line(df, x='lines_x', y='lines_y', color='name')

    fig.update_layout(
        autosize=True,
        width=1000,
        height=1000,
    )
    return fig


def superimpose_sheet_view(sheet1, sheet2, name_sheet=None):
    if name_sheet is None:
        name_sheet = ["sheet1", "sheet2"]

    a_lines_x, a_lines_y = _get_lines(sheet1, list('xy'))
    b_lines_x, b_lines_y = _get_lines(sheet2, list('xy'))

    df = pd.DataFrame([np.concatenate((a_lines_x, b_lines_x)),
                       np.concatenate((a_lines_y, b_lines_y)),
                       np.concatenate(
                           (np.repeat(name_sheet[0], len(a_lines_x)), np.repeat(name_sheet[1], len(b_lines_x))))])
    df = df.T
    df.columns = ['lines_x', 'lines_y', 'name']

    fig = px.line(df, x='lines_x', y='lines_y', color='name')

    fig.update_layout(
        autosize=True,
        width=1000,
        height=1000,
    )
    return fig


def view3d(mono, color='darkturquoise'):
    ipv.clear()
    draw_spec = config.draw.sheet_spec()
    draw_spec['face']['visible'] = True
    draw_spec['face']['color'] = color

    fig, meshes = sheet_view_3d(mono, **draw_spec)
    fig = ipv.gcf()

    fig.anglex = 1.0
    fig.angley = 0.2
    fig.anglez = 0.1

    ipv.show()
