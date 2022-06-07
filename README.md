# Genetic diversity and connectivity of chemosynthetic cold seep mussels from the U.S. Atlantic margin
### Supplementary code for DeLeo et al. 2022 (https://doi.org/10.1186/s12862-022-02027-4)

### Provides supplementary information on the code used to process and analyze the RADseq data use in this study

Analyses were run locally and/or on the Smithsonian's High-Performance Cluster (HPC) Hydra

----

## 1) Processing RADseq data 
##### QC raw data
`nohup fastqc L*_plate*.fastq.gz`

##### Create file that matches sample IDs to RAD sequencing barcodes

`mkdir barcodes`

`cd barcodes`

`nano plate3_sample_barcodes.txt`

TCCGGAGCGC	RB-19-114
CTAACACGGC	RB-19-115
AGCTTCGATT	RB-19-116
TCGCCGCAAT	RB-19-117
TCAGTTCCGG	RB-19-118
CGGAAGTGAG	RB-19-119
GTTGCTAGAC	RB-19-120
AATAGATTCA	RB-19-121
AGCTGATACA	RB-19-122
ATCAGTAGAA	RB-19-123
TCGTCTTAGT	RB-19-124
GCTCAGCCAG	RB-19-125
CGGCTACTTC	RB-19-126
CAAGCCGGTT	RB-19-127
TTGCGCAAGC	RB-19-128
TACGATGGAG	RB-19-129
GCAATATACA	RB-19-130
AAGAATTCGG	RB-19-131
TCGGCAGTCG	RB-19-132
AGTTCCATTG	RB-19-133
TTCTTGCGCT	RB-19-134
AGCAATCTAA	RB-19-135
GAATTGTCGC	RB-19-136
CTTCGACATA	RB-19-137
GAGATATGGT	RB-19-138
CTCCTTGGAG	RB-19-139
GTGTCTCTTG	RB-19-140
TGCAGTTATC	RB-19-141
TTCTGGAATA	RB-19-142
ACGCAACACA	RB-19-143
ACTGCCTCAA	RB-19-144
ACATCAATAT	RB-19-145
CCTCTTATCA	RB-19-150
TATCGTTAGT	RB-19-151
TAGTGCGGTC	RB-19-152
GGCCGGTAAC	RB-19-153
AGGAACCTCG	RB-19-154
TTATCCGTAG	RB-19-155
CGCTATACGG	RB-19-156
CACGCAACGA	RB-19-157
TGTCCTAGGA	RB-19-158
ATCCGTCTAC	RB-19-159
GGACTCACGG	RB-19-160
GCGTCCTGCC	RB-19-161
ACTTGACCGG	RB-19-162
AATGGTGACT	RB-19-164
CTAACAGTAT	RB-19-165
TCATAGGCTA	RB-19-166
GCTGCACGGT	RB-19-167
GCCGCAATGC	RB-19-168
CGCTTCTCTG	RB-19-169
CTCATTAACC	RB-19-170
TTGATGGTGC	RB-19-171
CAACATGAAG	RB-19-172
ATGAAGGCAG	RB-19-173
GATGGACTAA	RB-19-174
GCATGGAGGT	RB-19-175
GTATATCCAC	RB-19-176
CGTACCTTGC	RB-19-177
TATCGCGGAG	RB-19-178
CATGCATACT	CM-00128
TCACTGAGAA	CM-00129
ATCCATAAGA	CM-00130
CTGTTAGATT	CM-00131
TAACTGGTAC	CM-00132
GCCGGTGATT	CM-00133
AGACGAATAG	CM-00134
GGTCATTGTA	CM-00135
CGGTCGTTAC	CM-00136
GACGGACAGG	CM-00137
ATTGCCACCG	CM-00138
CAGAACCAGC	CM-00139
GCTGTGCAGA	CM-00140
AAGACCAATC	CM-00141
CGCGCGGCTG	CM-00142
GCATGAGGCG	CM-00143
TGACGGTGAT	CM-00144
CTTACCGGAG	CM-00146
ACACACATCA	CM-00148
TGTTGTCCGC	CM-00149
ATCCGCGACG	CM-00150
AAGTCGAGTA	CM-00151
CCAATCAAGA	CM-00153
TGAGCCAGCT	CM-00154
TCTAGAGAAG	CM-00155
GCCGAGGTGA	CM-00157
GGTGAGTCGG	CM-00160
CTTATTCTAC	CM-00161
CAAGAGACGT	CM-00163
GCACGTCTCC	CM-00165
GTGCTCTCTA	CM-00166
GCGTCAGATG	CM-00167
AATGAATCAG	CM-00168
ATTAGGAGGC	CM-00169
TTCTTCAGAC	CM-00158
CGCACTTGAT	FGXCONTROL

`nano plate4_sample_barcodes.txt`

TCCGGAGCGC	MASm32
CTAACACGGC	MASm36
AGCTTCGATT	HRS-1704-CM-004
TCGCCGCAAT	HRS-1704-CM-005
TCAGTTCCGG	HRS-1704-CM-007
CGGAAGTGAG	HRS-1704-CM-009
GTTGCTAGAC	HRS-1704-CM-011
AATAGATTCA	HRS-1704-CM-015
AGCTGATACA	HRS-1704-CM-017
ATCAGTAGAA	HRS-1704-CM-019
TCGTCTTAGT	HRS-1704-CM-021
GCTCAGCCAG	HRS-1704-CM-023
CGGCTACTTC	HRS-1704-CM-025
CAAGCCGGTT	HRS-1704-CM-028
TTGCGCAAGC	HRS-1704-CM-029
TACGATGGAG	HRS-1704-CM-031
GCAATATACA	HRS-1704-CM-033
AAGAATTCGG	HRS-1704-CM-037
TCGGCAGTCG	HRS-1704-CM-039
AGTTCCATTG	HRS-1704-CM-041
TTCTTGCGCT	HRS-1704-CM-043
AGCAATCTAA	HRS-1704-CM-045
GAATTGTCGC	HRS-1704-CM-047
CTTCGACATA	HRS-1704-CM-049
GAGATATGGT	HRS-1704-CM-051
CTCCTTGGAG	HRS-1704-CM-055
GTGTCTCTTG	HRS-1704-CM-058
TGCAGTTATC	HRS-1704-CM-061
TTCTGGAATA	HRS-1704-CM-063
ACGCAACACA	HRS-1704-CM-067
ACTGCCTCAA	HRS-1704-CM-069
ACATCAATAT	MAS283
CCTCTTATCA	MAS284
TATCGTTAGT	MAS285
TAGTGCGGTC	MAS286
GGCCGGTAAC	MAS288
AGGAACCTCG	MAS289
TTATCCGTAG	MAS290
CGCTATACGG	MAS291
CACGCAACGA	MAS292
TGTCCTAGGA	MAS293
ATCCGTCTAC	MAS295
GGACTCACGG	MAS297
GCGTCCTGCC	MAS298
ACTTGACCGG	MAS299
AATGGTGACT	MAS306
CTAACAGTAT	MAS310
TCATAGGCTA	MAS311
GCTGCACGGT	MAS313
GCCGCAATGC	MAS314
CGCTTCTCTG	MAS320
CTCATTAACC	MAS321
TTGATGGTGC	MAS322
CAACATGAAG	MAS323
ATGAAGGCAG	MAS326
GATGGACTAA	MAS327
GCATGGAGGT	MAS338
GTATATCCAC	MAS339
CGTACCTTGC	MAS340
TATCGCGGAG	MAS341
CATGCATACT	MAS343
TCACTGAGAA	MAS537
ATCCATAAGA	MAS539
CTGTTAGATT	MAS540
TAACTGGTAC	MAS541
GCCGGTGATT	MAS542
AGACGAATAG	MAS543
GGTCATTGTA	MAS544
CGGTCGTTAC	MAS545
GACGGACAGG	MAS547
ATTGCCACCG	MAS548
CAGAACCAGC	MAS549
GCTGTGCAGA	MAS550
AAGACCAATC	MAS551
CGCGCGGCTG	MAS552
GCATGAGGCG	MAS553
TGACGGTGAT	MAS554
CTTACCGGAG	MAS555
ACACACATCA	MAS556
TGTTGTCCGC	MAS557
ATCCGCGACG	MAS561
AAGTCGAGTA	HRS-1704-CM-35
CCAATCAAGA	CM-00201
TGAGCCAGCT	CM-00204
TCTAGAGAAG	CM-00205
GCCGAGGTGA	JSL-05-4890-02
GGTGAGTCGG	JSL-05-4890-07
CTTATTCTAC	JSL-05-4891-02
CAAGAGACGT	JSL-05-4891-06
GCACGTCTCC	JSL-05-4892-15
GTGCTCTCTA	JSL-05-4892-17
GCGTCAGATG	JSL-05-4893-02
AATGAATCAG	JSL-05-4894-04
ATTAGGAGGC	JSL-05-4894-18
TTCTTCAGAC	JSL-05-4894-23
CGCACTTGAT	FGXCONTROL

`cd ..`

`mkdir samples_plate3 samples_plate4`


##### Combine lane data for each plate
`cat L1_plate3.fastq.gz L2_plate3.fastq.gz L3_plate3.fastq.gz L4_plate3.fastq.gz > Plate3.fastq.gz &` 

`cat L1_plate4.fastq.gz L2_plate4.fastq.gz L3_plate4.fastq.gz L4_plate4.fastq.gz > Plate4.fastq.gz &`

##### Run STACKS process_radtags
> This program examines raw reads from an Illumina sequencing run and first, checks that the barcode and the RAD cutsite are intact, and demultiplexes the data. 
> If there are errors in the barcode or the RAD site within a certain allowance process_radtags can correct them.
> READ MORE: https://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php

Ran on HYDRA in interactive mode (qrsh = enter interactive mode, navigate to working directory)
	  
`qrsh` 

`module load gcc/7.3.0`

`module load bioinformatics/stacks` #load STACKS

`process_radtags -f /pool/genomics/deleod/bathymodiolus/Plate3.fastq.gz -o /pool/genomics/deleod/bathymodiolus/samples_plate3/ -b /pool/genomics/deleod/bathymodiolus/barcodes/plate3_sample_barcodes.txt -e pstI --inline_null -r -c -q -i gzfastq > process_radtags_plate3.out &`
> Example of output:
> 
> 383208399 total sequences                                                                                                      
> 12661721 barcode not found drops (3.3%)                                                                                       
> 273000 low quality read drops (0.1%)                                                                                        
> 1759364 RAD cutsite not found drops (0.5%)                                                                                   
> 368514314 retained reads (96.2%)   
		                                                                                            
``process_radtags -f /pool/genomics/deleod/bathymodiolus/Plate4.fastq.gz -o /pool/genomics/deleod/bathymodiolus/samples_plate4/ -b /pool/genomics/deleod/bathymodiolus/barcodes/plate4_sample_barcodes.txt -e pstI --inline_null -r -c -q -i gzfastq > process_radtags_plate4.out &
``
> output:
> 			  
> 381701674 total sequences                                                                                                      
> 12117674 barcode not found drops (3.2%)                                                                                       
> 71800 low quality read drops (0.1%)                                                                                        
> 1794663 RAD cutsite not found drops (0.5%)                                                                                   
> 367517537 retained reads (96.3%)                                                                                               
				  
##### Breakdown of parameters

> -f /PATH/TO/INPUT/FILE
> -o path to output the processed files  
> -b path to a file containing barcodes for this run
> -c clean data, remove any read with an uncalled base.
> -q discard reads with low quality scores.
> -r rescue barcodes and RAD-Tags.           
> -e provide the restriction enzyme used (cut site occurs on single-end read)
> -i input file type: 'gzfastq' (gzipped fastq)
> --inline_null: barcode is inline with sequence, occurs only on single-end read (default)
				  

##### Sort files into approriate folders
For this study, only interested in B.childressi (= G.childressi) and b. heckerae samples

`mkdir samples`

`mv samples_plate3/*.gz ./samples`

`mv samples_plate4/*.gz ./samples`

`cd samples`

`mkdir Lophelia`

`mv JSL-05-489* Lophelia/`

`mkdir desmophylum`

`mv CM-002* desmophylum/`

`mkdir bheckerae`

`mv HRS-1704-CM-35* bheckerae/`

`mv RB-19-1* bheckerae/`

`mv CM-001* bheckerae/`

`mkdir bchildressi`

`mv HRS-1704-CM-0* bchildressi/`

`mv MAS* bchildressi/`

----
## 2a) Assembing RADseq data with iPYRAD locally with jupyter notebooks

> Optional: create new environment and download packages there
> environment location: /Users/Odontodactylus/miniconda2/envs/myenv

``conda create --name myenv
conda activate myenv``

##### Install required packages 
conda install ipyrad -c conda-forge -c bioconda
conda install notebook -c conda-forge
conda install toytree -c conda-forge
conda install sra-tools -c bioconda
conda install entrez-direct -c bioconda
conda install mpi4py -c conda-forge

##### launch jupyter notebook#
jupyter notebook 

> Optional: In a seperate terminal window run ipcluster
> This will essentially start a number of independent Python processes (kernels) which we can then send bits of work to do. 

``ipcluster start --n=20``
> n = max number of cores you want to use
> 
##### Create an Assembly object
> This object stores the parameters of the assembly and the organization of data files.

``data = ip.Assembly("bathymodiolus")
data.set_params("project_dir", "/Users/Odontodactylus/Smithsonian/bathymodiolus")
data.set_params("sorted_fastq_path", "./fastq/*fq.gz")
data.set_params("clust_threshold", "0.95")
data.set_params("min_samples_locus", "51")
data.set_params("max_alleles_consens", 2)
data.set_params("max_shared_Hs_locus", "0.25")
data.set_params("filter_adapters", "2")
data.set_params("output_formats", ("p","s","v","n","k", "u"))``

``data.get_params()``

> example of output:
> 
> 	0   assembly_name               bathymodiolus                                
> 	1   project_dir                 ~/Smithsonian/bathymodiolus                  
> 	2   raw_fastq_path                                                           
> 	3   barcodes_path                                                            
> 	4   sorted_fastq_path           ./fastq/*fq.gz                               
> 	5   assembly_method             reference                                       
> 	6   reference_sequence          /pool/genomics/deleod/bathymodiolus/GCA_002080005.1_Bpl_v1.0_genomic.fna                                             
> 	7   datatype                    rad                                          
> 	8   restriction_overhang        ('TGCAG', '')                                
> 	9   max_low_qual_bases          5                                            
> 	10  phred_Qscore_offset         33                                           
> 	11  mindepth_statistical        6                                            
> 	12  mindepth_majrule            6                                            
> 	13  maxdepth                    10000                                        
> 	14  clust_threshold             0.85                                         
> 	15  max_barcode_mismatch        0                                            
> 	16  filter_adapters             2                                            
> 	17  filter_min_trim_len         35                                           
> 	18  max_alleles_consens         2                                            
> 	19  max_Ns_consens              0.05                                         
> 	20  max_Hs_consens              0.05                                         
> 	21  min_samples_locus           51                                           
> 	22  max_SNPs_locus              0.2                                          
> 	23  max_Indels_locus            8                                            
> 	24  max_shared_Hs_locus         0.50                                         
> 	25  trim_reads                  (0, 0, 0, 0)                                 
> 	26  trim_loci                   (0, 0, 0, 0)                                 
> 	27  output_formats              ['p', 's', 'v', 'n', 'k', 'u']               
> 	28  pop_assign_file                                                          
> 	29  reference_as_filter

###### Run iPYRAD
``data.run("1234567")``

> Alternatively, run steps separately, to check steps

##### run steps 1 & 2 of the assembly
``data.run("12")``

##### check the stats of the assembly (so far) from the .stats attribute
``data.stats``

##### run steps 3-6 of the assembly
> Note: Step 3 & 6 are the “clustering” steps, most data intensive part - will take a while

``data.run("3456")``

## 2b) Assembing RADseq data with iPYRAD on the HPC

> Files REQUIRED: 
> 
> 1) params.txt
> Edit ipyrad parameters file in working directory 
> update all paths to fastq files and working directories 
> 
> 2) ipyrad.job
> again, update all paths 
> 
> See example files

##### Run iPYRAD 
``qsub ipyrad.job``

##### Check output/progress log
``less ipyrad.log``

##### OPTIONAL: To re-run removing specific samples & low read samples
> create a new branch assembly called "subdata" including only the samples listed in the file samples_to_keep.txt

``ipyrad -p params.txt -b subdata samples_to_keep.txt``

>  loading Assembly: bchildressi
>  from saved path: /pool/genomics/deleod/bathymodiolus/samples/bchildressi/bchildressi.json
>  creating a new branch called 'subdata' with 70 Samples
>  writing new params file to params-subdata.txt
  
##### Then edit the desired parameters in the newly generated parameters file
``nano params-subdata.txt ``

##### Then resubmit job, running only the steps needed (this will depend on the parameters you change)
https://ipyrad.readthedocs.io/en/latest/6-params.html

``ipyrad -p params-subdata.txt -d -s 34567 -f -c $NSLOTS``
