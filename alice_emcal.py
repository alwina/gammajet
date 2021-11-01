import math
import numpy as np


def getSuperModule(cellId):
    if cellId < 11520:
        return math.floor(cellId / 1152)
    elif cellId < 12288:
        return 10 + math.floor((cellId - 11520) / 384)
    elif cellId < 16896:
        return 12 + math.floor((cellId - 12288) / 768)
    else:
        return 18 + (cellId - 16896) / 384


def getNphi(sm):
    if sm < 10:
        return 24
    elif sm < 12:
        return 8
    elif sm < 18:
        return 24
    else:
        return 8


def getIetaIphi(cellId, sm, nphi):
    if sm < 10:
        n0 = sm * 1152
    elif sm < 12:
        n0 = 11520 + (sm - 10) * 384
    elif sm < 18:
        n0 = 12288 + (sm - 12) * 768
    else:
        n0 = 16896 + (sm - 18) * 384

    n1 = cellId - n0

    ieta = 2 * math.floor(n1 / (2 * nphi)) + 1 - (n1 % 2)
    iphi = math.floor(n1 / 2) % nphi

    return ieta, iphi


def get5x5Cells(cellId):
    sm = getSuperModule(cellId)
    nphi = getNphi(sm)
    cells = []

    for dn in range(0, 10, 2):
        cells.append(cellId - 2 * nphi - 4 + dn)

        if cellId % 2 == 0:
            cells.append(cellId - 3 + dn)
        else:
            cells.append(cellId - 2 * nphi - 5 + dn)

        cells.append(cellId - 4 + dn)

        if cellId % 2 == 0:
            cells.append(cellId + 2 * nphi - 3 + dn)
        else:
            cells.append(cellId - 5 + dn)

        cells.append(cellId + 2 * nphi - 4 + dn)

    return cells


def getCrossCells(cellId):
    sm = getSuperModule(cellId)
    nphi = getNphi(sm)
    cells = []

    if cellId % 2 == 0:
        cells.append(cellId + 1)
    else:
        cells.append(cellId - 2 * nphi - 1)

    cells.append(cellId - 2)
    cells.append(cellId + 2)

    if cellId % 2 == 0:
        cells.append(cellId + 2 * nphi + 1)
    else:
        cells.append(cellId - 1)

    return cells


# port from http://alidoc.cern.ch/AliRoot/v5-09-11/_ali_e_m_c_a_l_rec_point_8cxx_source.html#l01002
def calculateShowerShapeFromCells(cellIds, cell_e, cluster_e):
    wtot = 0.0
    x = 0.0
    z = 0.0
    dxx = 0.0
    dzz = 0.0
    dxz = 0.0

    for cellId in cellIds:
        sm = getSuperModule(cellId)
        nphi = getNphi(sm)
        ieta, iphi = getIetaIphi(cellId, sm, nphi)

        # this handles shared clusters but also doesn't do anything for
        # clusters in a single supermodule, so don't bother checking
        # specifically for shared clusters
        if (sm % 2):
            ieta = ieta + 48

        w = max(0, 4.5 + math.log(cell_e[cellId] / cluster_e))
        dxx = dxx + w * ieta * ieta
        x = x + w * ieta
        dzz = dzz + w * iphi * iphi
        z = z + w * iphi
        dxz = dxz + w * ieta * iphi
        wtot = wtot + w

    if wtot > 0:
        x = x / wtot
        z = z / wtot
        dxx = dxx / wtot - x * x
        dzz = dzz / wtot - z * z
        dxz = dxz / wtot - x * z

    return 0.5 * (dxx + dzz) + math.sqrt(0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz)


def calculateShowerShapes5x5(icluster, cluster_e, cellMaxId, cell_e, cell_cluster_index):
    # first, get the cell Ids of the surrounding 5x5 grid
    cells5x5 = get5x5Cells(cellMaxId)

    # this is not supposed to happen, so return a nonsensical number to see how there are
    # and also to not count them
    if cell_e[cellMaxId] < 0.1:
        return {
            '5x5contiguouscluster': -1.0,
            '5x5contiguous': -1.0,
            '5x5cluster': -1.0,
            '5x5all': -1.0
        }

    cells_ss5x5contiguouscluster = set()
    cells_ss5x5contiguous = set()
    cells_ss5x5cluster = set()
    cells_ss5x5all = set()

    # shower shape that requires contiguous cells in the 5x5 in the cluster:
    # recursively, starting with the max energy cell (cellMaxId):
    # 1. look for the crossing cells
    # 2. for each cell, check that it is not already in the set of cells to be used in the calculation
    # 3. if so, determine if it is in the 5x5
    # 4. if so, determine if it is in the cluster
    # 5. if so, add to the set of cells to be used in the calculation
    # 6. recurse
    cells_ss5x5contiguouscluster.add(cellMaxId)
    checkedCells = set()
    checkedCells.add(cellMaxId)

    def getContiguousClusterCells(cellId):
        crossCells = getCrossCells(cellId)
        for crossCellId in crossCells:
            if crossCellId in checkedCells:
                continue
            checkedCells.add(crossCellId)
            if crossCellId > 17663 or crossCellId < 0:
                continue
            if cell_e[crossCellId] < 0.1 or np.isnan(cell_e[crossCellId]):
                continue
            if crossCellId not in cells5x5:
                continue
            if cell_cluster_index[crossCellId] != icluster:
                continue
            cells_ss5x5contiguouscluster.add(crossCellId)
            getContiguousClusterCells(crossCellId)

    getContiguousClusterCells(cellMaxId)

    # shower shape that requires contiguous cells in the 5x5, not necessarily in the cluster:
    # recursively, starting with the max energy cell (cellMaxId):
    # 1. look for the crossing cells
    # 2. for each cell, check that it is not already in the set of cells to be used in the calculation
    # 3. if so, determine if it is in the 5x5
    # 4. if so, add it to the set of cells to be used in the calculation
    # 5. recurse
    # calculate
    cells_ss5x5contiguous.add(cellMaxId)
    checkedCells = set()
    checkedCells.add(cellMaxId)

    def getContiguousCells(cellId):
        crossCells = getCrossCells(cellId)
        for crossCellId in crossCells:
            if crossCellId in checkedCells:
                continue
            checkedCells.add(crossCellId)
            if crossCellId > 17663 or crossCellId < 0:
                continue
            if cell_e[crossCellId] < 0.1 or np.isnan(cell_e[crossCellId]):
                continue
            if crossCellId not in cells5x5:
                continue
            cells_ss5x5contiguous.add(crossCellId)
            getContiguousCells(crossCellId)

    getContiguousCells(cellMaxId)

    # shower shape that takes all cells in the 5x5 in the cluster, not necessarily contiguous
    # for each cell in the 5x5, add it to the set of cells to be used in the calculation if it is in the cluster

    # shower shape that takes all cells in the 5x5, not necessarily contiguous or in the cluster
    # for each cell in the 5x5, add it to the set of cells to be used in the calculation

    for cellId in cells5x5:
        if cellId > 17663 or cellId < 0:
            continue
        if cell_e[cellId] < 0.1 or np.isnan(cell_e[cellId]):
            continue
        cells_ss5x5all.add(cellId)
        if cell_cluster_index[cellId] == icluster:
            cells_ss5x5cluster.add(cellId)

    ss5x5contiguouscluster = calculateShowerShapeFromCells(cells_ss5x5contiguouscluster, cell_e, cluster_e)
    ss5x5contiguous = calculateShowerShapeFromCells(cells_ss5x5contiguous, cell_e, cluster_e)
    ss5x5cluster = calculateShowerShapeFromCells(cells_ss5x5cluster, cell_e, cluster_e)
    ss5x5all = calculateShowerShapeFromCells(cells_ss5x5all, cell_e, cluster_e)

    return {
        '5x5contiguouscluster': ss5x5contiguouscluster,
        '5x5contiguous': ss5x5contiguous,
        '5x5cluster': ss5x5cluster,
        '5x5all': ss5x5all
    }
