INTRODUCTION
------------





REQUIREMENTS
------------
Python 3.0 and newer versions
 




INSTALLATION
------------
Download the archive from : (github page)

unzip the files : gunzip archive.gz



USAGE
------------
Before launching any analysis make sure all the required packages are install
try : pip install <module>


The program is a python3 script and can be launched in the shell.

e.g.: python3 <scriptName> 



OPTIONS
-----------

List of possible options;

-f / --file : File name and its extension. The script requires an alignement file in fasta format.

-s / --size : Minimum size of the stretch of N. Default value is 5.

-w / --window : Size of the windows surroundings the stretches of N in which the program will search for sequencing errors. Default value is 10. (see figure 1)

-r / --reference : File name of the reference genome in fasta format. Default file is Reference_SARS-CoV-2.fasta which is the reference from <....>.

-p / --protected : File name of the positions in the genome to protect in bed format. If not specify, all positions are subject to changes.

-l / --look : File name of the positions in the genome to look at for possible sequencing errors in bed format. The 4th column can specify an allele. If not specify, all positions are subject to changes.

-o / --output : Prefix of the outputs files. Default value is <output>.



OUTPUTS
-----------

<output>.txt : File that has additionnal informations on every Flagged site in every sequence. The columns are in order: SequenceID, Position of the Flag site, Distance from the closest stretch of N, An attribute for protected attribute (0 if not protected, 1 if it correspond to the protected alternative allele (if specified) or to any alternative allele (if not specified), 2 if it is protected (allele specified) and doesn't correspond to the alternative allele and reference), Allele at this Flag site.  

<output>_nb_flag_per_seq.txt : Statistic file of the number of Flag sites in every sequence from the input file.

<output>.fasta : Corrected alignement file in fasta format. This files correspond to the alignement input file with the Flagged SNP changed for a "N" (unknown nucleotide).



EXAMPLES
-----------

> python3 flag_SNP_FINAL_v5_cmd_line.py -f <filename>.fasta -s 5 -w 10 -r <Reference_SARS-CoV-2>.fasta -p <SNP_to_protect>.bed -o <output>
Will look at all sites in every sequence in the <filename> file.


> python3 flag_SNP_FINAL_v5_cmd_line.py -f <filename>.fasta -s 5 -w 10 -r <Reference_SARS-CoV-2>.fasta -p <SNP_to_protect>.bed -o <output>
Will look at all sites in every sequence in the <filename> file but won't modified the Flag sites if they are in the <SNP_to_protect> file.


> python3 flag_SNP_FINAL_v5_cmd_line.py -f <filename>.fasta -s 5 -w 10 -r <Reference_SARS-CoV-2>.fasta -l <Sites_to_look_at>.bed -o <output>
Will only look at the sites specified in the <Sites_to_look_at> file in every sequence in the <filename>.


> python3 flag_SNP_FINAL_v5_cmd_line.py -f <filename>.fasta -s 5 -w 10 -r <Reference_SARS-CoV-2>.fasta -p <SNP_to_protect>.bed -l <Sites_to_look_at>.bed -o <output>
Will only look at the sites specified in the <Sites_to_look_at> file in every sequence in the <filename> but won't modified the Flag sites if they are in the <SNP_to_protect> file.



AUTHOR
-----------
PELLETIER Justin (https://mhi-omics.org/people/justin-pelletier/)
email: justin.pelletier@umontreal.ca
Montreal Heart Institute
Julie Hussin's computational biomedecine lab (https://mhi-omics.org/people/julie-hussin/)

