import numpy as np
import pandas as pd
from tyssue.utils.utils import to_nd
from tyssue.dynamics import units
from tyssue.dynamics import effectors
from tyssue.geometry.planar_geometry import PlanarGeometry
from tyssue import MonolayerGeometry


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


class AnisotropicLineTension(effectors.AbstractEffector):
    label = "Anisotropic Line Tension"

    @staticmethod
    def energy(sheet):
        return sheet.edge_df.eval('gamma * length * is_active')  # accounts for half edges

    @staticmethod
    def gradient(sheet):
        grad_srce = -sheet.edge_df[sheet.ucoords] * to_nd(
            sheet.edge_df.eval("gamma * is_active"), len(sheet.coords)
        )
        grad_srce.columns = ["g" + u for u in sheet.coords]
        grad_trgt = -grad_srce
        return grad_srce, grad_trgt


from tyssue.dynamics.effectors import elastic_force, elastic_energy


class PlaneBarrierElasticity(effectors.AbstractEffector):
    """
    Barrier use to maintain the tissue integrity.
    """

    dimensions = units.line_elasticity
    magnitude = "plane barrier_elasticity"
    label = "Plane barrier elasticity"
    element = "vert"
    specs = {
        "vert": {"barrier_elasticity": 1.0, "is_active": 1, "z": 0.0, "z_barrier": 1}
    }  # distance to a barrier membrane

    @staticmethod
    def energy(eptm):
        return eptm.vert_df.eval("0.5*(z_distance)**2 * barrier_elasticity")

    @staticmethod
    def gradient(eptm):
        kl_l0 = elastic_force(eptm.vert_df, "z_distance", "barrier_elasticity", 0)
        grad = eptm.vert_df[eptm.coords] * to_nd(kl_l0, eptm.dim)
        grad.columns = ["g" + u for u in eptm.coords]
        grad["gx"] = 0
        grad["gy"] = 0
        # if "z" in eptm.coords:
        #     grad["gz"] = 0
        return grad, grad


class BarrierElasticity(effectors.AbstractEffector):
    """
    Barrier use to maintain the tissue integrity.
    """

    dimensions = units.line_elasticity
    magnitude = "barrier_elasticity"
    label = "Barrier elasticity"
    element = "vert"
    specs = {
        "vert": {"barrier_elasticity": 1.0, "is_active": 1, "delta_rho": 0.0}
    }  # distance to a barrier membrane

    @staticmethod
    def energy(eptm):
        return eptm.vert_df.eval("0.5*delta_rho**2 * barrier_elasticity")

    @staticmethod
    def gradient(eptm):
        kl_l0 = elastic_force(eptm.vert_df, "delta_rho", "0", "0")
        grad = eptm.vert_df[eptm.coords] * to_nd(kl_l0, eptm.dim)
        grad.columns = ["g" + u for u in eptm.coords]
        grad["gy"] = 0
        if "z" in eptm.coords:
            grad["gz"] = 0
        return grad, grad


class ShearPlanarGeometry(PlanarGeometry):
    @classmethod
    def update_all(cls, sheet):
        PlanarGeometry.update_all(sheet)
        cls.update_gamma(cls, sheet)

    @staticmethod
    def update_gamma(cls, sheet):
        gamma_0 = sheet.edge_df['gamma_0']
        phi = sheet.specs['edge']['phi0_apical']

        e_angle = cls.get_phis(sheet)
        sheet.edge_df['angle'] = e_angle
        # sheet.edge_df["angle"] = [np.pi + a if a < 0 else a for a in e_angle]
        sheet.edge_df["gamma"] = gamma_0 * np.cos(2 * (e_angle - phi))

    @classmethod
    def get_phis(cls, sheet):
        if "dx" not in sheet.edge_df:
            cls.update_dcoords(sheet)
            cls.update_centroid(sheet)

        return np.arctan2(sheet.edge_df["dy"], sheet.edge_df["dx"])


class ShearMonolayerGeometry(MonolayerGeometry):
    @classmethod
    def update_all(cls, sheet):
        MonolayerGeometry.update_all(sheet)
        cls.update_gamma(cls, sheet)
        cls.update_zdistance(cls, sheet)

    @staticmethod
    def update_zdistance(cls, sheet):
        z_barrier = sheet.specs['cell']['z_barrier']
        sheet.vert_df.loc[sheet.vert_df['segment'] == 'apical', 'z_barrier'] = z_barrier
        sheet.vert_df.loc[sheet.vert_df['segment'] == 'basal', 'z_barrier'] = -z_barrier
        sheet.vert_df.loc[sheet.vert_df['segment'] == 'lateral', 'z_barrier'] = z_barrier * 10

        sheet.vert_df['z_distance'] = (np.abs(sheet.vert_df["z"])
                                       - np.abs(sheet.vert_df["z_barrier"])
                                       ).clip(lower=0)

    @staticmethod
    def update_gamma(cls, sheet):
        gamma_0 = sheet.edge_df['gamma_0']
        # phi = sheet.specs['edge']['phi0']

        phi = pd.to_numeric(sheet.edge_df['gamma_0']) * 0
        phi[sheet.edge_df[sheet.edge_df['segment'] == 'apical'].index] = sheet.specs['edge']['phi0_apical']
        phi[sheet.edge_df[sheet.edge_df['segment'] == 'basal'].index] = sheet.specs['edge']['phi0_basal']
        phi = phi.to_numpy()

        e_angle = cls.get_phis(sheet)
        sheet.edge_df['angle'] = e_angle
        # sheet.edge_df["angle"] = [np.pi + a if a < 0 else a for a in e_angle]
        sheet.edge_df["gamma"] = gamma_0 * np.cos(2 * (e_angle - phi))

    @classmethod
    def get_phis(cls, sheet):
        if "dx" not in sheet.edge_df:
            cls.update_dcoords(sheet)
            cls.update_centroid(sheet)

        return np.arctan2(sheet.edge_df["dy"], sheet.edge_df["dx"])


from tyssue import PlanarGeometry


class EllipsisGeometry(PlanarGeometry):
    """ """

    @classmethod
    def update_all(cls, eptm):
        PlanarGeometry.update_all(eptm)
        cls.update_lumen_volume(eptm)
        cls.update_height(eptm)
        cls.update_tilt(eptm)

    @staticmethod
    def update_lumen_volume(eptm):
        if eptm.settings['inside'] == 'apical':
            inside_edge = eptm.apical_edges
        else:
            inside_edge = eptm.basal_edges
        srce_pos = eptm.upcast_srce(eptm.vert_df[["x", "y"]]).loc[inside_edge]
        trgt_pos = eptm.upcast_trgt(eptm.vert_df[["x", "y"]]).loc[inside_edge]
        apical_edge_pos = (srce_pos + trgt_pos) / 2
        apical_edge_coords = eptm.edge_df.loc[inside_edge, ["dx", "dy"]]
        eptm.settings["lumen_volume"] = (
                -apical_edge_pos["x"] * apical_edge_coords["dy"]
                + apical_edge_pos["y"] * apical_edge_coords["dx"]
        ).values.sum()
        eptm.settings["lumen_vol"] = eptm.settings["lumen_volume"]

    @staticmethod
    def update_height(eptm):
        a = eptm.settings["a"]
        b = eptm.settings["b"]
        h = eptm.settings["barrier_height"]
        eptm.vert_df["theta"] = np.arctan((eptm.vert_df.y / eptm.vert_df.x).clip(-1, 1))
        eptm.vert_df["barrier_rho"] = np.sqrt(
            ((a + h) * np.cos(eptm.vert_df["theta"])) ** 2 + ((b + h) * np.sin(eptm.vert_df["theta"])) ** 2)
        # eptm.vert_df["basal_shift"] = (
        #         eptm.vert_df["barrier_rho"] - eptm.specs["vert"]["basal_shift"]
        # )
        eptm.vert_df["delta_rho"] = (
                np.linalg.norm(eptm.vert_df[["x", "y"]], axis=1)
                - eptm.vert_df["barrier_rho"]
        ).clip(lower=0)

    # def update_tilt(eptm):
