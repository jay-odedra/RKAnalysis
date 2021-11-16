### source /cvmfs/sft.cern.ch/lcg/views/LCG_96python3/x86_64-centos7-gcc8-opt/setup.sh
from sklearn.externals import joblib
from xgboost import XGBClassifier
import glob
import xgboost as xgb  

if __name__ == "__main__":
  directory = 'models_LowPtPF_v7.2'
  for f in glob.glob('{}/*.pkl'.format(directory)):
    print('dumping BDT output reading from...'.format(f))
    bdt = joblib.load(f)

    print("bdt:")
    print(bdt)

    print("booster:")
    booster = bdt._Booster
    print(booster)

    print("dtest:")
    dtest = xgb.DMatrix('test_otto.txt')  
    print(dtest)

    print("now pred:")
    ypred = booster.predict(dtest) 
    print(ypred)

    print('Done')
