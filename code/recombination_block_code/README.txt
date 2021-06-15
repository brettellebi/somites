Step1.

Script: extract_dp4.sh
Inputs: path to all sorted BAMs, the reference FASTA file and a SITES_FILE (the sites file is the selected set of SNP i.e. well behaved homozygous divergent calls in F0) and output dir.
Format of sites file is this:
1	21200	21200	T	C	0/0	1/1
1	21497	21497	A	G	0/0	1/1
1	22027	22027	T	G	0/0	1/1
1	22345	22345	A	T	0/0	1/1
1	22662	22662	T	C	0/0	1/1
1	23344	23344	A	C	0/0	1/1
1	23345	23345	T	C	0/0	1/1
1	24358	24358	T	C	0/0	1/1
1	24884	24884	A	G	0/0	1/1

i.e. chr,pos,pos,ref,alt,A_genotype(F0_1), B_genotype(F0_2)

—> outputs “dp4” files for all input BAM.
Nb. currently uses jobs array and the LSB_JOBINDEX variable


Step2.

Script: make_dp_AB.sh
Inputs: “dp4_dir” (location of dp4 files), “dpAB_dir” (output location), SITES_FILE.
Calls the following perl script: map_parent_genotypes.pl

—> output dp_AB files for all dp4 files (ie. maps the allele to A and B based on genotype in the SITES file).
Nb. currently uses jobs array and the LSB_JOBINDEX variable

Step3.

Script: run_rc_block.R
Inputs: “dpAB_dir”
Sources: rc_block_hmm.R

—> a file containing all called recombination block across all samples.
NB. this currently joint calls all blocks (this is most reliable) but it would be possible to call in batch to save mem / speed.
Try joint call everything first - i think it should be fine.

Nb. after this i have some more steps that we could automate - i.e. generate “segmented recombination blocks” and realised relationship matrix based on segmented blocks.
Plus a load of plotting functionally and prepare the association tests.


But i think start with these main processing steps will be a good starting point :)

Let me know any questions - happy to help out on any of the items - just let me know.

Hopefully once you get into the scripts the flow should be quite clear and basically just need to change the hard coded inputs and outputs into arguments and hoping this can all to stringed together into a pipeline process.

Cheers,

Tom

p.s. if its useful i can provide some test datasets to you - although this should work off any set of sorted BAMs.

Oh and finally there is a samtools index command in the first script - you may want to remove this if its already done prior to these steps.