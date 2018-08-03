# cacao_gbs

Please verify accuracy of inputs and outputs. These scripts have been tested on multiple sequencing datasets. 

Perform in silico digestion of Theobroma cacao cv. Matina  

Before running the following scripts its necessary to perform the digestion of the reference genome with the enzyme of interest. To do this:

1. Install rebaseextract (http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/rebaseextract.html#input.1) 

2. export PATH "/usr/local/emboss/bin"

3. Download latest versions of withrefm (ftp://ftp.neb.com/pub/rebase/withrefm.txt) and proto (ftp://ftp.neb.com/pub/rebase/proto.txt) files. 

4. Move withrefm and proto file to PATH=/usr/local/emboss/share/EMBOSS/data

5. Process the REBASE database for use by restriction enzyme applications
REBASE database withrefm file: /usr/local/emboss/share/EMBOSS/data/withrefm
REBASE database proto file: /usr/local/emboss/share/EMBOSS/data/proto

6. Perform restrict command line: restrict -sequence reference_genome.fa -enzyme 'BsaXI' -outfile cacao.restrict

7. split_InSilico.sh

This script is meant to split the resulting output file (cacao.restrict) into a folder containing digestion files per chromosome.  

8. Summarize_InSilico_Digest_clean2.R

This script takes the folder digestion files per chromosome and performs the summarize of in silico digestions per chromosome and plot the distribution of fragments after restrict and the distributions of enzymes cuts throughout the chromosomes.

cacao_GBS_pipeline.sh

This script takes performs check quality, trimming, aligement and SNP calling process. 
