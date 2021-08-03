import math


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
