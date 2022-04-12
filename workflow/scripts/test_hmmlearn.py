import sys,os
import logging, traceback
logging.basicConfig(filename=snakemake.log[0],
                    level=logging.INFO,
                    format='%(asctime)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    )
def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger.error(''.join(["Uncaught exception: ",
                         *traceback.format_exception(exc_type, exc_value, exc_traceback)
                         ])
                 )
# Install exception handler
sys.excepthook = handle_exception


# Import libraries

import numpy as np
import pandas as pd
from hmmlearn import hmm

# Import variables

## Debug
IN = "/hps/nobackup/birney/users/ian/somites/hmm_in/F2/hdrr/F1_het_min_DP/15/5000.csv"
MOD = "A"
LOW_COV_SAMPLES = [26,89,166,178,189,227,239,470,472,473,490,501,502,503,504,505,506,507,508,509,510,511]

## True
IN = snakemake.input[0]
MOD = snakemake.params.mod
LOW_COV_SAMPLES = snakemake.params.low_cov_samples
OUT = snakemake.output[0]

######################
# Setup
######################

# Read in file

df = pd.read_csv(IN, comment = "#")

# Remove low coverage samples

df = df[~df['SAMPLE'].isin(LOW_COV_SAMPLES)]

# get length (number of rows) for each sample

s_lens = df.groupby('SAMPLE').size().values.tolist()
s_chr_lens = df.groupby(['SAMPLE', 'CHROM']).size().values.tolist()

# Filter DF for first five samples

s_keep = df['SAMPLE'].unique()[0:5].tolist()
d_trim = df[df['SAMPLE'].isin(s_keep)]

s_lens_trim = d_trim.groupby('SAMPLE').size().values.tolist()
s_chr_lens_trim = d_trim.groupby(['SAMPLE', 'CHROM']).size().values.tolist()

# Create HMM input

hmm_in = d_trim.loc[:, ['PROP_CAB']] # Needs to remain as a data frame

######################
# HMM mods
######################

# Set means so that the states are in the right order (0,1,2 for Cab, Het, Kaga)
means = np.array([[0],
                 [0.5],
                 [1]])

## Mod A

### Single sequence

if MOD == "A":
  model = hmm.GaussianHMM(n_components=3,
                          covariance_type="diag",
                          n_iter=100,
                          init_params = 'sct' # 'm' should be excluded from initialisations so that we can set them
                          ) 
  model.means_ = means
  model.fit(hmm_in) 
  states = model.predict(hmm_in)

## Mod B

### Separate chromosomes and samples

if MOD == "B":
  model = hmm.GaussianHMM(n_components=3,
                          covariance_type="diag",
                          n_iter=100,
                          init_params = 'sct'
                          )
  model.means_ = means
  model.fit(hmm_in, lengths = s_chr_lens_trim)
  states = model.predict(hmm_in)

## Mod C

### Transition matrix:

#### There are ~80k SNPs per sample
#### We want 2-4 blocks per chromosome = 1-3 transitions = 80000 / (3 * 24) = up to 1 change every say 1000 SNPs
#### 1 - (1/1000) = 0.999 probability of staying the same
#### Double probability to transition to HET state than HOM

if MOD == "C":
  tmat = np.array([[0.999, 0.00066, 0.00033],
                   [0.0005, 0.999, 0.0005],
                   [0.00033, 0.00066, 0.999]])
  model = hmm.GaussianHMM(n_components=3,
                          covariance_type="diag",
                          n_iter=100,
                          algorithm = 'viterbi',
                          params = 'smc', # 't' excluded so that it doesn't get updated during training
                          init_params = 'sc' # 's' (startprob), 't' (transmat), 'm' (means), or 'c' (covars)
                          )
  model.transmat_ = tmat
  model.means_ = means
  model.fit(hmm_in, lengths = s_chr_lens_trim)
  states = model.predict(hmm_in)
  
## Mod D

#### 'Error' states

if MOD == "D":
  # set means
  means = np.array([[0],
                   [0],
                   [0.5],
                   [0.5],
                   [1],
                   [1]])
  tmat = np.array([[0.999, 0.00066, 0.00033], # 
                   [0.0005, 0.999, 0.0005],
                   [0.00033, 0.00066, 0.999]])
  covars = np.array([[[]]])
  # Note changed to MultinomialHMM
  model = hmm.GaussianHMM(n_components=3,
                          covariance_type="diag",
                          n_iter=100,
                          algorithm = 'viterbi',
                          params = 'smc',
                          init_params = 'sc' # 's' (startprob), 't' (transmat), 'm' (means), or 'c' (covars)
                          )
  model.transmat_ = tmat
  model.fit(hmm_in, lengths = s_chr_lens_trim)
  states = model.predict(hmm_in)





