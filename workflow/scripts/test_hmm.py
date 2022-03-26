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
REP_PHENO = "data/20220321_phenotypes.xlsx"

# Extract reporter chr, start, and end

REP_CHR = int(REPORTER_LOC.split(':')[0])
REP_STA = int(REPORTER_LOC.split(':')[1].split('-')[0])
REP_END = int(REPORTER_LOC.split(':')[1].split('-')[1])

# Read in files

df = pd.read_csv(IN, comment = "#")
rep_pheno = pd.read_excel(REP_PHENO)

# Filter for REP_CHR chromosome

df_chr = df.loc[df['CHROM'] == REP_CHR, ]
df_chr = df_chr.reset_index()

# Remove NAs
df_chr.dropna(inplace = True)

# get length (number of rows) for each sample

sample_lengths = df_chr.groupby('SAMPLE').size().values

# Pull out rows with the reporter

key_cols = ['CHROM', 'BIN_START', 'BIN_END', 'SAMPLE', 'STATE_IMP']
rep_df = df_chr.loc[(df_chr['BIN_START'] < REP_STA) & (df_chr['BIN_END'] > REP_STA), key_cols]

# Recode `STATE_IMP`
recode_dict = {0:-1, 1:0, 2:1}
rep_df = rep_df.replace(dict(STATE_IMP=recode_dict))

# Join with the reporter phenotype

rep_df = rep_df.merge(rep_pheno[['SAMPLE', 'reporter_pheno']], how = 'left')

# Concordance between `STATE_IMP` and `reporter_pheno`

sum(rep_df['STATE_IMP'] == rep_df['reporter_pheno']) / len(rep_df)
## 82.57%

###################
# Proportion of Cab
###################

# Create column with proportion of Cab reads

df_chr['PROP_CAB'] = df_chr.apply(
    lambda row: row.CAB_READS /  (row.CAB_READS + row.KAGA_READS), axis=1
    )
    
# Create input

hmm_propc = df_chr['PROP_CAB'].values.reshape(-1, 1)

# Set up model

model_propc = hmm.GaussianHMM(n_components=3, covariance_type="diag", n_iter=100)

# Train

model_propc.fit(hmm_propc, lengths = sample_lengths)

# Predict

df_chr['STATE_PROPC'] = model_propc.predict(hmm_propc)



