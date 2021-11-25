# -- Cluster MC info -- #
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


# -- MC weights -- #
fixedWeightProductions = ['17g6a3', '17g6a1']


def getFixedWeight(filename):
    if '17g6a3' in filename:
        if 'pthat1' in filename:
            weight = 4.47e-11
        elif 'pthat2' in filename:
            weight = 9.83e-11
        elif 'pthat3' in filename:
            weight = 1.04e-10
        elif 'pthat4' in filename:
            weight = 1.01e-10
        elif 'pthat5' in filename:
            weight = 6.93e-11
        elif 'pthat6' in filename:
            weight = 5.13e-11
        elif 'pthat7' in filename:
            weight = 3.03e-11
        elif 'pthat8' in filename:
            weight = 1.89e-11
    elif '17g6a1' in filename:
        if('pthat1' in filename):
            weight = 1.60e-11
        elif('pthat2' in filename):
            weight = 2.72e-12
        elif('pthat3' in filename):
            weight = 3.96e-13
        elif('pthat4' in filename):
            weight = 6.14e-14
        elif('pthat5' in filename):
            weight = 1.27e-14

    return weight


# -- PDG code info -- #
# for some reason, python can't read the mc_truth_charge branch
# because it's a Char_t and the \0xfd character is not valid in utf-8
# actually figuring this out is too hard, so instead, we use the PDG code
# so far, it seems there's a relatively small number of codes that
# actually exist in the embedded ntuples
# [-3334 -3322 -3312 -3222 -3122 -3112 -2212 -2112  -321  -211   -16   -14
#    -13   -12   -11    11    12    13    14    16    22   130   211   310
#    321  2112  2212  3112  3122  3222  3312  3322  3334]
pdgCodeToCharge = {}
pdgCodeToCharge[-3334] = 1   # Omega-bar
pdgCodeToCharge[-3322] = 0   # Xi0bar
pdgCodeToCharge[-3312] = 1   # Xi-bar
pdgCodeToCharge[-3222] = -1  # Sigma+bar
pdgCodeToCharge[-3122] = 0   # Lambdabar
pdgCodeToCharge[-3112] = 1   # Sigma-bar
pdgCodeToCharge[-2212] = -1  # pbar
pdgCodeToCharge[-2112] = 0   # nbar
pdgCodeToCharge[-321] = -1   # K-
pdgCodeToCharge[-211] = -1   # pi-
pdgCodeToCharge[-16] = 0     # nubar (tau)
pdgCodeToCharge[-14] = 0     # nubar (mu)
pdgCodeToCharge[-13] = 1     # mu+
pdgCodeToCharge[-12] = 0     # nubar (e)
pdgCodeToCharge[-11] = 1     # e+
pdgCodeToCharge[11] = -1     # e-
pdgCodeToCharge[12] = 0      # nu (e)
pdgCodeToCharge[13] = -1     # mu-
pdgCodeToCharge[14] = 0      # nu (mu)
pdgCodeToCharge[16] = 0      # nu (tau)
pdgCodeToCharge[22] = 0      # gamma
pdgCodeToCharge[130] = 0     # K0L
pdgCodeToCharge[211] = 1     # pi+
pdgCodeToCharge[310] = 0     # K0S
pdgCodeToCharge[321] = 1     # K+
pdgCodeToCharge[2112] = 0    # n
pdgCodeToCharge[2212] = 1    # p
pdgCodeToCharge[3112] = -1   # Sigma-
pdgCodeToCharge[3122] = 0    # Lambda
pdgCodeToCharge[3222] = 1    # Sigma+
pdgCodeToCharge[3312] = -1   # Xi-
pdgCodeToCharge[3322] = 0    # Xi0
pdgCodeToCharge[3334] = -1   # Omega-
