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
IN = "/hps/nobackup/birney/users/ian/somites/processed_recomb/F2/all_sites/20000.csv"
REPORTER_LOC = "16:28706898-28708417"
REP_PHENO = 

# Extract reporter chr, start, and end

REP_CHR = int(REPORTER_LOC.split(':')[0])
REP_STA = int(REPORTER_LOC.split(':')[1].split('-')[0])
REP_END = int(REPORTER_LOC.split(':')[1].split('-')[1])

# Read in file

df = pd.read_csv(IN, comment = "#")

# Filter for REP_CHR chromosome

df_chr = df.loc[df['CHROM'] == REP_CHR, ]
df_chr = df_chr.reset_index()

# Remove NAs
df_chr.dropna(inplace = True)

# Create column with proportion of Cab reads

df_chr['PROP_CAB'] = df_chr.apply(
    lambda row: row.CAB_READS /  (row.CAB_READS + row.KAGA_READS), axis=1)
    
# Create input

hmm_in = df_chr['PROP_CAB'].values.reshape(-1, 1)

# get length (number of rows) for each sample

sample_lengths = df_chr.groupby('SAMPLE').size().values

# Set up model

model = hmm.GaussianHMM(n_components=3, covariance_type="diag", n_iter=100)

# Train

model.fit(hmm_in, lengths = sample_lengths)

# Predict

df_chr['STATE'] = model.predict(hmm_in)



