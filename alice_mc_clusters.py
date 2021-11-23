def getNMCPhotons(cluster_nmc_truth, cluster_mc_truth_index, mc_truth_pdg_code):
    photonCount = 0
    for i in range(min(cluster_nmc_truth, len(cluster_mc_truth_index))):
        mcindex = cluster_mc_truth_index[i]
        if mcindex == 65535:
            continue
        if mcindex >= len(mc_truth_pdg_code):
            continue
        if mc_truth_pdg_code[mcindex] == 22:
            photonCount = photonCount + 1

    return photonCount


# see the ntuplizer for how mc_truth_is_prompt_photon is determined
# it is only filled for photons and electrons (and positrons), so that needs to be checked first
# current definition of prompt photon cluster: any cluster that contains
# at least one photon/electron/positron that came from a prompt photon
def getIsPrompt(cluster_nmc_truth, cluster_mc_truth_index, mc_truth_pdg_code, mc_truth_is_prompt_photon):
    for i in range(min(cluster_nmc_truth, len(cluster_mc_truth_index))):
        mcindex = cluster_mc_truth_index[i]
        # if this is the dummy index, ignore it
        if mcindex == 65535:
            continue
        # if this is somehow out of range, ignore it
        if mcindex >= len(mc_truth_pdg_code):
            continue
        if mcindex >= len(mc_truth_is_prompt_photon):
            continue
        # if this is not a photon or electron, ignore it
        if mc_truth_pdg_code[mcindex] not in (11, -11, 22):
            continue
        # if this is from a prompt photon, then the cluster is considered a prompt photon cluster
        if mc_truth_is_prompt_photon[mcindex]:
            return True
    return False


def getTruthPtAndComponents(cluster_nmc_truth, cluster_mc_truth_index, mc_truth_pdg_code, mc_truth_is_prompt_photon, mc_truth_pt):
    truthpt = 0
    truthcomponents = 0
    for i in range(min(cluster_nmc_truth, len(cluster_mc_truth_index))):
        mcindex = cluster_mc_truth_index[i]
        # if this is the dummy index, ignore it
        if mcindex == 65535:
            continue
        # if this is somehow out of range, ignore it
        if mcindex >= len(mc_truth_pdg_code):
            continue
        if mcindex >= len(mc_truth_is_prompt_photon):
            continue
        # if this is not a photon or electron, ignore it
        if mc_truth_pdg_code[mcindex] not in (11, -11, 22):
            continue
        # if it is not from a prompt photon, ignore it
        if not mc_truth_is_prompt_photon[mcindex]:
            continue

        # account for the truth particles in this cluster: 1 = photon, 100 = electron, 10000 = positron
        # surely I could do 1/10/100, but juuuuuust in case
        if mc_truth_pdg_code[mcindex] == 22:
            truthcomponents += 1
        elif mc_truth_pdg_code[mcindex] == 11:
            truthcomponents += 100
        elif mc_truth_pdg_code[mcindex] == -11:
            truthcomponents += 10000

        # and then add the truth pt to the total
        truthpt += mc_truth_pt[mcindex]

    return truthpt, truthcomponents


def getParentPi0Pt(cluster_nmc_truth, cluster_mc_truth_index, mc_truth_first_parent_pdg_code, mc_truth_first_parent_pt):
    parentpi0pt = -1
    for i in range(min(cluster_nmc_truth, len(cluster_mc_truth_index))):
        mcindex = cluster_mc_truth_index[i]
        # if this is the dummy index, ignore it
        if mcindex == 65535:
            continue
        # if this is somehow out of range, ignore it
        if mcindex >= len(mc_truth_first_parent_pdg_code):
            continue
        if mcindex >= len(mc_truth_first_parent_pt):
            continue
        # if the parent is not a pi0, ignore it
        # not sure if pi0 can have negative pdg code, but check just in case
        if abs(mc_truth_first_parent_pdg_code[mcindex]) != 111:
            continue
        # for now, take the min pt of the pi0s we see
        # unclear how often we end up with more than one, but just in case,
        # we'll be consistent about which one we take
        if parentpi0pt == -1:
            parentpi0pt = mc_truth_first_parent_pt[mcindex]
        else:
            parentpi0pt = min(parentpi0pt, mc_truth_first_parent_pt[mcindex])

    return parentpi0pt
