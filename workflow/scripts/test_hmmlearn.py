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
import pickle

# Import variables

## Debug
IN = "/hps/nobackup/birney/users/ian/somites/hmm_in/F2/hdrr/F1_het_min_DP/15/5000.csv"
MOD = "A"
LOW_COV_SAMPLES = [26,89,166,178,189,227,239,470,472,473,490,501,502,503,504,505,506,507,508,509,510,511]

## True
IN = snakemake.input[0]
MOD = snakemake.params.mod[0]
LOW_COV_SAMPLES = snakemake.params.low_cov_samples
LOW_COV_SAMPLES = [int(i) for i in LOW_COV_SAMPLES] # convert to integers
OUT_CSV = snakemake.output.csv
OUT_PCK = snakemake.output.pck

######################
# Setup
######################

# Read in file

df = pd.read_csv(IN, comment = "#")

# Remove low coverage samples

df = df[~df['SAMPLE'].isin(LOW_COV_SAMPLES)]

# get length (number of rows) for each sample / chr

s_chr_lens = df.groupby(['SAMPLE', 'CHROM']).size().values.tolist()

######################
# Filter DF for first five samples (hash out when running on full sample)
######################

s_keep = df['SAMPLE'].unique()[0:10].tolist()
df = df[df['SAMPLE'].isin(s_keep)]

s_chr_lens = df.groupby(['SAMPLE', 'CHROM']).size().values.tolist()

######################
######################

# Create HMM input

hmm_in = df.loc[:, ['PROP_KAGA']] # Needs to remain as a data frame

######################
# HMM mods
######################

# Set initial means so that the states are in the right order (0,1,2 for Cab, Het, Kaga)
means = np.array([[0],
                 [0.5],
                 [1]])

# Params: GaussianHMM(params='stmc', init_params='stmc')
# From the docs:
# The parameters that get updated during (params) or initialized before (init_params) the training. 
# Can contain any combination of ‘s’ for startprob, ‘t’ for transmat, and other characters for subclass-specific emission parameters. 
# Defaults to all parameters.

# So if you want to initialize the means, set `init_params = "stc"`

## Mod A

### Single sequence

if MOD == "A":
  model = hmm.GaussianHMM(n_components=3,
                          covariance_type="diag",
                          n_iter=100,
                          init_params = "stc"
                          )
  model.means_ = means
  model.fit(hmm_in, lengths = s_chr_lens)
  states = model.predict(hmm_in)


## Mod B

### Transition matrix:

#### There are ~80k SNPs per sample
#### We want 2-4 blocks per chromosome = 1-3 transitions = 80000 / (3 * 24) = up to 1 change every say 1000 bins
#### 1 - (1/1000) = 0.999 probability of staying the same
#### Double probability to transition to HET state than HOM

if MOD == "B":
  model = hmm.GaussianHMM(n_components=3,
                          covariance_type="diag",
                          n_iter=100,
                          algorithm = 'viterbi',
                          init_params = 'stc' # 's' (startprob), 't' (transmat), 'm' (means), or 'c' (covars)
                          )
  model.means_ = means
  model.fit(hmm_in, lengths = s_chr_lens)
  model.transmat_ = np.array([[0.999, 0.00066, 0.00033],
                              [0.0005, 0.999, 0.0005],
                              [0.00033, 0.00066, 0.999]]
                              )
  states = model.predict(hmm_in)
  
## Mod C

### Wide, even covars

if MOD == "C":
  model = hmm.GaussianHMM(n_components=3,
                          covariance_type="diag",
                          n_iter=100,
                          algorithm = 'viterbi',
                          init_params = 'stc' # 's' (startprob), 't' (transmat), 'm' (means), or 'c' (covars)
                          )
  model.means_ = means
  model.fit(hmm_in, lengths = s_chr_lens)
  model.transmat_ = np.array([[0.999, 0.00066, 0.00033],
                              [0.0005, 0.999, 0.0005],
                              [0.00033, 0.00066, 0.999]]
                              )
  # Note: in the models above, these tend to be 
  # [0.001, 0.09, 0.001]
  model.covars_ = np.array([
                     [0.333],   # Cab true 
                     [0.333],   # Het true
                     [0.333],   # Kaga true 
                    ])
  states = model.predict(hmm_in)

## Mod D

### Narrow, even covars

if MOD == "D":
  model = hmm.GaussianHMM(n_components=3,
                          covariance_type="diag",
                          n_iter=100,
                          algorithm = 'viterbi',
                          init_params = 'stc' # 's' (startprob), 't' (transmat), 'm' (means), or 'c' (covars)
                          )
  model.means_ = means
  model.fit(hmm_in, lengths = s_chr_lens)
  model.transmat_ = np.array([[0.999, 0.00066, 0.00033],
                              [0.0005, 0.999, 0.0005],
                              [0.00033, 0.00066, 0.999]]
                              )
  # Note: in the models above, these tend to be 
  # [0.001, 0.09, 0.001]
  model.covars_ = np.array([
                     [0.01],   # Cab true 
                     [0.01],   # Het true
                     [0.01],   # Kaga true 
                    ])
  states = model.predict(hmm_in)
  
#### 'Error' states (no HET bin)

# 1 error state for each true state
# with broader emission properties
# and a cost to go in and quite an aggressive decay (ie exit probability is high)

if MOD == "E":
  
  model = hmm.GaussianHMM(n_components=5,
                          covariance_type="diag",
                          n_iter=100,
                          algorithm = 'viterbi',
                          init_params = 'stc' # 's' (startprob), 't' (transmat), 'm' (means), or 'c' (covars)
                          )

  model.means_ = np.array([
                           [0],   # Cab true 
                           [0.1], # Cab error
                           [0.5], # Het true
                           [0.9], # Kaga error
                           [1]    # Kaga true
                           ])
  
  model.fit(hmm_in, lengths = s_chr_lens)
  
  model.transmat_ = np.array([
                             [0.85, 0.149, 0.00066, 0, 0.00033],  # Cab true 
                             [0.999, 1e-03, 0, 0, 0],  # Cab error
                             [0.0005, 0, 0.999, 0, 0.0005],   # Het true
                             [0, 0, 0, 1e-03, 0.999], # Kaga error
                             [0.00033, 0, 0.00066, 0.149, 0.85]  # Kaga true 
                            ])
                            
  model.covars_ = np.array([
                     [0.2],   # Cab true 
                     [1],     # Cab error
                     [0.2],   # Het true
                     [1],   # Kaga error 
                     [0.2]      # Kaga true                  
                    ])
                    
  states = model.predict(hmm_in)

## Mod F

### Wide, even covars

if MOD == "F":
  model = hmm.GaussianHMM(n_components=3,
                          covariance_type="diag",
                          n_iter=100,
                          algorithm = 'viterbi',
                          init_params = 'stc' # 's' (startprob), 't' (transmat), 'm' (means), or 'c' (covars)
                          )
  model.means_ = means
  model.fit(hmm_in, lengths = s_chr_lens)
  model.transmat_ = np.array([[0.999, 0.00066, 0.00033],
                              [0.0005, 0.999, 0.0005],
                              [0.00033, 0.00066, 0.999]]
                              )
  # Note: in the models above, these tend to be 
  # [0.001, 0.09, 0.001]
  model.covars_ = np.array([
                     [0.8],   # Cab true 
                     [0.8],   # Het true
                     [0.8],   # Kaga true 
                    ])
  states = model.predict(hmm_in)

## Mod G

### Wide, even covars

if MOD == "G":
  model = hmm.GaussianHMM(n_components=3,
                          covariance_type="diag",
                          n_iter=100,
                          algorithm = 'viterbi',
                          init_params = 'stc' # 's' (startprob), 't' (transmat), 'm' (means), or 'c' (covars)
                          )
  model.means_ = means
  model.fit(hmm_in, lengths = s_chr_lens)
  model.transmat_ = np.array([[0.999, 0.00066, 0.00033],
                              [0.0005, 0.999, 0.0005],
                              [0.00033, 0.00066, 0.999]]
                              )
  # Note: in the models above, these tend to be 
  # [0.001, 0.09, 0.001]
  model.covars_ = np.array([
                     [1],   # Cab true 
                     [1],   # Het true
                     [1],   # Kaga true 
                    ])
  states = model.predict(hmm_in)


######################
# Add to df
######################

df = df.assign(STATE = states)

######################
# Write to file
######################

## Model
pickle.dump(model, open(OUT_PCK, 'wb'))
## csv
df.to_csv(OUT_CSV, index = False)

