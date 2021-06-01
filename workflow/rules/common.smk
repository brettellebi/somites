import pandas as pd

###### Config file and sample sheets #####
configfile: "config/config.yaml"


samples = pd.read_table(config["samples"]).set_index("sample", drop=False)