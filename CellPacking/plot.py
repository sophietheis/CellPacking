import pandas as pd
import numpy as np
import plotly.express as px

from tyssue.draw.plt_draw import _get_lines


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
