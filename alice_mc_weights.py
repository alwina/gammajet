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
