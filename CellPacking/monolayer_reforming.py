from tyssue.config.geometry import bulk_spec
from tyssue.topology.base_topology import collapse_edge
from tyssue.topology.bulk_topology import split_vert
from tyssue.core.history import HistoryHdf5

import pandas as pd
import numpy as np


def history_monolayer_from_sheets(apical_history, basal_history, hf5file=None):
    # 1
    init_monolayer = monolayer_from_sheets(apical_history.retrieve(0).datasets, basal_history.retrieve(0).datasets,
                                           distance=1)
    monolayer = init_monolayer.copy(deep_copy=True)
    # 2
    if hf5file is None:
        monolayer_history = HistoryHdf5(init_monolayer)
    else:
        monolayer_history = HistoryHdf5(init_monolayer, hf5file=hf5file)

    # 3
    for i in apical_history.time_stamps:
        apical_sheet = apical_history.retrieve(i)
        basal_sheet = basal_history.retrieve(i)

        # 3a
        if apical_history.trackevent[i]['remove_edge'][0] != -1:
            for e in apical_history.trackevent[i]['remove_edge']:
                remove_e = monolayer_history.sheet.edge_df[
                    (monolayer_history.sheet.edge_df['segment'] == 'apical') & (
                                monolayer_history.sheet.edge_df['id_sheet'] == e)].index
                collapse_edge(monolayer_history.sheet, remove_e, allow_two_sided=False)
                monolayer_history.sheet.edge_df.loc[
                    monolayer_history.sheet.edge_df['segment'] == 'apical', 'id_sheet'] = apical_sheet.edge_df.index
                monolayer_history.sheet.vert_df.loc[
                    monolayer_history.sheet.vert_df['segment'] == 'apical', 'id_sheet'] = apical_sheet.vert_df.index

        monolayer_history.record(time_stamp=i)
    return monolayer_history


def monolayer_from_sheets(apical_datasets, basal_datasets, distance=1):
    coords = list("xyz")
    datasets = {}

    apical_vert = apical_datasets["vert"].copy()
    apical_face = apical_datasets["face"].copy()
    apical_edge = apical_datasets["edge"].copy()

    apical_vert["segment"] = "apical"
    apical_face["segment"] = "apical"
    apical_edge["segment"] = "apical"

    apical_vert["id_sheet"] = apical_vert.index
    apical_face["id_sheet"] = apical_face.index
    apical_edge["id_sheet"] = apical_edge.index

    apical_vert['z'] = distance / 2
    apical_face['z'] = distance / 2
    apical_edge['z'] = distance / 2

    Nv = apical_vert.index.max() + 1
    Ne = apical_edge.index.max() + 1
    Nf = apical_face.index.max() + 1

    cell_df = apical_face[coords].copy()
    cell_df.index.name = "cell"
    cell_df["is_alive"] = 1

    basal_vert = basal_datasets["vert"].copy()
    basal_face = basal_datasets["face"].copy()
    basal_edge = basal_datasets["edge"].copy()

    basal_vert["segment"] = "basal"
    basal_face["segment"] = "basal"
    basal_edge["segment"] = "basal"

    basal_vert["id_sheet"] = basal_vert.index
    basal_face["id_sheet"] = basal_face.index
    basal_edge["id_sheet"] = basal_edge.index

    basal_vert['z'] = -distance / 2
    basal_face['z'] = -distance / 2
    basal_edge['z'] = -distance / 2

    basal_vert.index = basal_vert.index + Nv
    basal_face.index = basal_face.index + Nf

    apical_edge["cell"] = apical_edge["face"]
    basal_edge["cell"] = basal_edge["face"]
    # ## Flip edge so that normals are outward
    basal_edge[["srce", "trgt"]] = basal_edge[["trgt", "srce"]] + Nv
    basal_edge["face"] = basal_edge["face"] + Nf
    basal_edge.index = basal_edge.index + Ne

    lateral_face = pd.DataFrame(
        index=apical_edge.index + 2 * Nf, columns=apical_face.columns
    )
    lateral_face["segment"] = "lateral"
    lateral_face["is_alive"] = 1

    lateral_edge = pd.DataFrame(
        index=np.arange(2 * Ne, 6 * Ne), columns=apical_edge.columns
    )

    lateral_edge["cell"] = np.repeat(apical_edge["cell"].values, 4)
    lateral_edge["face"] = np.repeat(lateral_face.index.values, 4)
    lateral_edge["segment"] = "lateral"

    lateral_edge.loc[np.arange(2 * Ne, 6 * Ne, 4), ["srce", "trgt"]] = apical_edge[
        ["trgt", "srce"]
    ].values

    lateral_edge.loc[np.arange(2 * Ne + 1, 6 * Ne, 4), "srce"] = apical_edge[
        "srce"
    ].values
    lateral_edge.loc[np.arange(2 * Ne + 1, 6 * Ne, 4), "trgt"] = basal_edge[
        "trgt"
    ].values

    lateral_edge.loc[np.arange(2 * Ne + 2, 6 * Ne, 4), ["srce", "trgt"]] = basal_edge[
        ["trgt", "srce"]
    ].values

    lateral_edge.loc[np.arange(2 * Ne + 3, 6 * Ne, 4), "srce"] = basal_edge[
        "srce"
    ].values
    lateral_edge.loc[np.arange(2 * Ne + 3, 6 * Ne, 4), "trgt"] = apical_edge[
        "trgt"
    ].values

    datasets["cell"] = cell_df
    datasets["vert"] = pd.concat([apical_vert, basal_vert])
    datasets["vert"]["is_active"] = 1
    datasets["edge"] = pd.concat([apical_edge, basal_edge, lateral_edge])
    datasets["face"] = pd.concat([apical_face, basal_face, lateral_face])
    datasets["edge"]["is_active"] = 1
    specs = bulk_spec()

    for elem in ["vert", "edge", "face", "cell"]:
        datasets[elem].index.name = elem
        for col, value in specs[elem].items():
            if col not in datasets[elem]:
                datasets[elem][col] = value

    return datasets


def IO_transition(eptm, edge, recenter=False):
    """
    I → H transition as defined in Okuda et al. 2013
    (DOI 10.1007/s10237-012-0430-7).
    See tyssue/doc/illus/IH_transition.png for the algorithm
    """
    srce, trgt, face, cell = eptm.edge_df.loc[edge, ["srce", "trgt", "face", "cell"]]
    vert = min(srce, trgt)
    collapse_edge(eptm, edge)

    return vert, face


def OH_transition(eptm, vert, face, recenter=False):
    """
    I → H transition as defined in Okuda et al. 2013
    (DOI 10.1007/s10237-012-0430-7).
    See tyssue/doc/illus/IH_transition.png for the algorithm
    """
    split_vert(eptm, vert, face, recenter=recenter)
    return 0
