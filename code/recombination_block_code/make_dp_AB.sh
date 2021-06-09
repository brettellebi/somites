SITES_FILE=sites.HdrR_Ho5.filtered
dp_dir=../dp4_all
dpAB_dir=../dpAB_all

FILES=($( ls $dp_dir  ))
f=${FILES[$LSB_JOBINDEX-1]}
fbname=$(basename "$f" .txt)

perl map_parent_genotypes.pl $dp_dir/$f $SITES_FILE > $dpAB_dir/$fbname.txt

