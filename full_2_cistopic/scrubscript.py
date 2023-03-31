import pycisTopic
pycisTopic.__version__

import warnings
warnings.filterwarnings('ignore')
warnings.simplefilter('ignore')

import pickle
import pandas as pd

import os
wdir = '/dodrio/scratch/projects/starting_2022_023/benchmark/scatac_benchmark/full_2_cistopic/'
os.chdir(wdir)

import glob
from collections import OrderedDict

import scrublet as scr
import pandas as pd
import matplotlib.pyplot as plt

filenames = sorted(glob.glob('cistopic_objects/*__cto.pkl'))
cto_dict = {}
for filename in filenames:
    sample = filename.split('/')[-1].split('__cto.pkl')[0]
    cto_dict[sample] = filename
    print(sample)


# load objects into dict:
for sample in cto_dict.keys():
    with open(cto_dict[sample], 'rb') as f:
        cto = pickle.load(f)
    print(f"Loaded {cto_dict[sample]}")
    scrub = scr.Scrublet(cto.fragment_matrix.T, expected_doublet_rate=0.1)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()
    scrub.plot_histogram()
