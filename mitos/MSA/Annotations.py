from datetime import datetime


def annotate_modification(msa, modification):
    feature = "MIToS_{:}".format(datetime.now().date())
    msa.annotations[feature] = '{:}\n{:}'.format(msa.annotations.get(feature, ''),
                                                 modification)
    return True
