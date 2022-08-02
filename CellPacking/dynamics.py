import numpy as np
from tyssue.utils.utils import to_nd
from tyssue.dynamics import effectors
from tyssue.geometry.planar_geometry import PlanarGeometry


class Compression(effectors.AbstractEffector):

    @staticmethod
    def energy(sheet):
        return sheet.vert_df.eval('compression * x**2')

    @staticmethod
    def gradient(sheet):
        grad = sheet.vert_df[sheet.coords].copy()
        grad.columns = ['gx', 'gy']
        grad['gy'] = 0
        grad['gx'] = sheet.vert_df.eval("compression * x")
        return grad, None


class Shear(effectors.AbstractEffector):
    @staticmethod
    def energy(sheet):
        return sheet.edge_df.eval('gamma * length /2 * is_active')  # accounts for half edges

    @staticmethod
    def gradient(sheet):
        grad_srce = -sheet.edge_df[sheet.ucoords] * to_nd(
            sheet.edge_df.eval("gamma * is_active / 2"), len(sheet.coords)
        )
        grad_srce.columns = ["g" + u for u in sheet.coords]
        grad_trgt = -grad_srce
        return grad_srce, grad_trgt


class ShearPlanarGeometry(PlanarGeometry):
    @classmethod
    def update_all(cls, sheet):
        PlanarGeometry.update_all(sheet)
        cls.update_gamma(cls, sheet)

    @staticmethod
    def update_gamma(cls, sheet):
        gamma_0 = sheet.specs['edge']['gamma_0']
        phi = sheet.specs['edge']['phi']

        e_angle = cls.get_phis(sheet)
        sheet.edge_df['angle'] = e_angle
        # sheet.edge_df["angle"] = [np.pi + a if a < 0 else a for a in e_angle]
        sheet.edge_df["gamma"] = gamma_0 * np.cos(2 * (e_angle * phi))

    @classmethod
    def get_phis(cls, sheet):
        if "dx" not in sheet.edge_df:
            cls.update_dcoords(sheet)
            cls.update_centroid(sheet)

        return np.arctan2(sheet.edge_df["dy"], sheet.edge_df["dx"])
