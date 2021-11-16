### source /cvmfs/sft.cern.ch/lcg/views/LCG_96python3/x86_64-centos7-gcc8-opt/setup.sh
from sklearn.externals import joblib
from xgboost import XGBClassifier
import glob

if __name__ == "__main__":
  directory = 'models_LowPtPF_v7.2'
  for f in glob.glob('{}/*.pkl'.format(directory)):
    print('converting model {} to txt file...'.format(f))
    bdt = joblib.load(f)

    print("bdt")
    print(bdt)

    booster = bdt._Booster
    print("booster")
    print(booster)

    print("now dump")
    booster.dump_model(f.replace('.pkl', '.txt'))

    print('Done')
