from xgboost import XGBClassifier
from sklearn.datasets import make_classification
from sklearn.externals import joblib
import numpy as np

def get_model(pkl):
    model = joblib.load(pkl)

    def _monkey_patch():
        return model._Booster

    if isinstance(model.booster, basestring):
        model.booster = _monkey_patch
    return model

base = get_model(
    './'
    '/xgb_fulldata_02May2021_lowq2RK1_bothSB_tightPreselectionWODmass_onePerEvent_pauc2_pf_fold1.model.pkl')

print(base)
print ""

booster = base._Booster
print booster

booster.dump_model("xgb_fulldata_02May2021_lowq2RK1_bothSB_tightPreselectionWODmass_onePerEvent_pauc2_pf_fold1.model.txt")
