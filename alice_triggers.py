import csv


allTriggerClasses = {}
with open('all-trigger-classes.csv') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        # skip the header
        if row[0] == 'Run #':
            continue
        # remove empty strings
        allTriggerClasses[int(row[0])] = [ts for ts in row[1:] if ts]


def getTriggerIds(runNumber, triggerStrings):
    triggerClass = allTriggerClasses[runNumber]
    triggerIDs = [triggerClass.index(ts) for ts in triggerStrings if ts in triggerClass]

    return triggerIDs


def getCentralTriggerIds(runNumber):
    # kCentral,"+CV0H7-[B|U|T]-NOPF-[CENT|CENTNOTRD|CENTNOPMD]"
    triggerStrings = []
    triggerStrings.append("CV0H7-B-NOPF-CENT")
    triggerStrings.append("CV0H7-B-NOPF-CENTNOTRD")
    triggerStrings.append("CV0H7-B-NOPF-CENTNOPMD")
    triggerStrings.append("CV0H7-U-NOPF-CENT")
    triggerStrings.append("CV0H7-U-NOPF-CENTNOTRD")
    triggerStrings.append("CV0H7-U-NOPF-CENTNOPMD")
    triggerStrings.append("CV0H7-T-NOPF-CENT")
    triggerStrings.append("CV0H7-T-NOPF-CENTNOTRD")
    triggerStrings.append("CV0H7-T-NOPF-CENTNOPMD")

    return getTriggerIds(runNumber, triggerStrings)


def getINT7TriggerIds(runNumber):
    # kINT7,"+[CINT7|CINT7ZAC|CV0L7]-[B|U|T]-NOPF-[CENT|CENTNOTRD|CENTNOPMD]"
    triggerStrings = []
    triggerStrings.append("CINT7-B-NOPF-CENT")
    triggerStrings.append("CINT7-B-NOPF-CENTNOTRD")
    triggerStrings.append("CINT7-B-NOPF-CENTNOPMD")
    triggerStrings.append("CINT7-U-NOPF-CENT")
    triggerStrings.append("CINT7-U-NOPF-CENTNOTRD")
    triggerStrings.append("CINT7-U-NOPF-CENTNOPMD")
    triggerStrings.append("CINT7-T-NOPF-CENT")
    triggerStrings.append("CINT7-T-NOPF-CENTNOTRD")
    triggerStrings.append("CINT7-T-NOPF-CENTNOPMD")
    triggerStrings.append("CINT7ZAC-B-NOPF-CENT")
    triggerStrings.append("CINT7ZAC-B-NOPF-CENTNOTRD")
    triggerStrings.append("CINT7ZAC-B-NOPF-CENTNOPMD")
    triggerStrings.append("CINT7ZAC-U-NOPF-CENT")
    triggerStrings.append("CINT7ZAC-U-NOPF-CENTNOTRD")
    triggerStrings.append("CINT7ZAC-U-NOPF-CENTNOPMD")
    triggerStrings.append("CINT7ZAC-T-NOPF-CENT")
    triggerStrings.append("CINT7ZAC-T-NOPF-CENTNOTRD")
    triggerStrings.append("CINT7ZAC-T-NOPF-CENTNOPMD")
    triggerStrings.append("CV0L7-B-NOPF-CENT")
    triggerStrings.append("CV0L7-B-NOPF-CENTNOTRD")
    triggerStrings.append("CV0L7-B-NOPF-CENTNOPMD")
    triggerStrings.append("CV0L7-U-NOPF-CENT")
    triggerStrings.append("CV0L7-U-NOPF-CENTNOTRD")
    triggerStrings.append("CV0L7-U-NOPF-CENTNOPMD")
    triggerStrings.append("CV0L7-T-NOPF-CENT")
    triggerStrings.append("CV0L7-T-NOPF-CENTNOTRD")
    triggerStrings.append("CV0L7-T-NOPF-CENTNOPMD")

    return getTriggerIds(runNumber, triggerStrings)


def getSemiCentralTriggerIds(runNumber):
    # kSemiCentral,"+CMID7-[B|U|T]-NOPF-[CENT|CENTNOTRD|CENTNOPMD]"
    triggerStrings = []
    triggerStrings.append("CMID7-B-NOPF-CENT")
    triggerStrings.append("CMID7-B-NOPF-CENTNOTRD")
    triggerStrings.append("CMID7-B-NOPF-CENTNOPMD")
    triggerStrings.append("CMID7-U-NOPF-CENT")
    triggerStrings.append("CMID7-U-NOPF-CENTNOTRD")
    triggerStrings.append("CMID7-U-NOPF-CENTNOPMD")
    triggerStrings.append("CMID7-T-NOPF-CENT")
    triggerStrings.append("CMID7-T-NOPF-CENTNOTRD")
    triggerStrings.append("CMID7-T-NOPF-CENTNOPMD")

    return getTriggerIds(runNumber, triggerStrings)


def getEMCEGATriggerIds(runNumber):
    # kEMCEGA,"+C[INT|EMC|DMC]7[E|D]G[1|2]-B-NOPF-[CENT|CENTNOTRD|CENTNOPMD],C[EMC|DMC]7[E|D]G2PER-B-NOPF-CENTNOPMD"
    triggerStrings = []
    triggerStrings.append("CINT7EG1-B-NOPF-CENT")
    triggerStrings.append("CINT7EG1-B-NOPF-CENTNOTRD")
    triggerStrings.append("CINT7EG1-B-NOPF-CENTNOPMD")
    triggerStrings.append("CINT7EG2-B-NOPF-CENT")
    triggerStrings.append("CINT7EG2-B-NOPF-CENTNOTRD")
    triggerStrings.append("CINT7EG2-B-NOPF-CENTNOPMD")
    triggerStrings.append("CINT7DG1-B-NOPF-CENT")
    triggerStrings.append("CINT7DG1-B-NOPF-CENTNOTRD")
    triggerStrings.append("CINT7DG1-B-NOPF-CENTNOPMD")
    triggerStrings.append("CINT7DG2-B-NOPF-CENT")
    triggerStrings.append("CINT7DG2-B-NOPF-CENTNOTRD")
    triggerStrings.append("CINT7DG2-B-NOPF-CENTNOPMD")
    triggerStrings.append("CEMC7EG1-B-NOPF-CENT")
    triggerStrings.append("CEMC7EG1-B-NOPF-CENTNOTRD")
    triggerStrings.append("CEMC7EG1-B-NOPF-CENTNOPMD")
    triggerStrings.append("CEMC7EG2-B-NOPF-CENT")
    triggerStrings.append("CEMC7EG2-B-NOPF-CENTNOTRD")
    triggerStrings.append("CEMC7EG2-B-NOPF-CENTNOPMD")
    triggerStrings.append("CEMC7DG1-B-NOPF-CENT")
    triggerStrings.append("CEMC7DG1-B-NOPF-CENTNOTRD")
    triggerStrings.append("CEMC7DG1-B-NOPF-CENTNOPMD")
    triggerStrings.append("CEMC7DG2-B-NOPF-CENT")
    triggerStrings.append("CEMC7DG2-B-NOPF-CENTNOTRD")
    triggerStrings.append("CEMC7DG2-B-NOPF-CENTNOPMD")
    triggerStrings.append("CDMC7EG1-B-NOPF-CENT")
    triggerStrings.append("CDMC7EG1-B-NOPF-CENTNOTRD")
    triggerStrings.append("CDMC7EG1-B-NOPF-CENTNOPMD")
    triggerStrings.append("CDMC7EG2-B-NOPF-CENT")
    triggerStrings.append("CDMC7EG2-B-NOPF-CENTNOTRD")
    triggerStrings.append("CDMC7EG2-B-NOPF-CENTNOPMD")
    triggerStrings.append("CDMC7DG1-B-NOPF-CENT")
    triggerStrings.append("CDMC7DG1-B-NOPF-CENTNOTRD")
    triggerStrings.append("CDMC7DG1-B-NOPF-CENTNOPMD")
    triggerStrings.append("CDMC7DG2-B-NOPF-CENT")
    triggerStrings.append("CDMC7DG2-B-NOPF-CENTNOTRD")
    triggerStrings.append("CDMC7DG2-B-NOPF-CENTNOPMD")
    triggerStrings.append("CEMC7EG2PER-B-NOPF-CENTNOPMD")
    triggerStrings.append("CEMC7DG2PER-B-NOPF-CENTNOPMD")
    triggerStrings.append("CDMC7EG2PER-B-NOPF-CENTNOPMD")
    triggerStrings.append("CDMC7DG2PER-B-NOPF-CENTNOPMD")

    return getTriggerIds(runNumber, triggerStrings)


def isEventSelected(triggerIds, triggerMask):
    # apparently this is the faster way to test what's in the position of an integer
    # https://stackoverflow.com/questions/9298865/get-n-th-bit-precision-of-integer
    # for an event, only one of the triggers has to fire for it to be selected
    for triggerId in triggerIds:
        if not not(int(triggerMask) & (1 << triggerId)):
            return True
    return False
