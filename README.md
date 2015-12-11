# samfire
SAMFIRE carries out multi-locus variant calling for time-resolved sequence data.  Instructions for using the code follow.

Compiling the code

Samfire is code for calling multi-locus variants from time-resolved sequence data.  It requires prior installation of the GSL library, and is compiled with the command line

make samfire

Running the code

Samfire processes data via a series of subprograms:

./samfire filter carries out QC on a set of data, returning a list of filtered reads
./samfire sl_traj calls single locus trajectories from the filtered data
./samfire sl_noise estimates the extent of noise in the single-locus trajectory data, in addition to that caused by the finite sampling of data
./samfire sl_neutrality identifies trajectories that deviate significantly from a pattern of evolution under neutrality
./samfire call_ml calls multi-locus trajectories
./samfire ml_noise estimates the extent of noise in the multi-locus trajectory data, in addition to that caused by the finite sampling of data

Details of the individual programs follows:
  

./samfire filter --ref <ref> [options]: carries out QC on a set of data, returning a list of filtered reads.  Here <ref> is the name of a reference file in fasta format.

Required inputs:

Input_files.in: A list of files in .sam format, ordered by the time at which they were collected.  In each case the full path to the file should be specified

e.g. 

Data/Data_t0.sam
Data/Data_t1.sam
Data/Data_t2.sam

A reference sequence, in fasta format.

Outputs:

Aligned?.out is a list of aligned single reads in the format:

<Position of first base> <Sequence>

Files are output for all inputs specified in Input_files.in.  Reads in each file report a null “-“ allele for all positions that fail the applied tests of quality.

Joined?.out is a list of joined and aligned sequences.  Sequences are reported as in Aligned?.out, albeit that paired-end reads are joined by ‘-‘ nucleotides.  Files are produced for all inputs specified in Input_files.in.  If paired end reads are not specified, these files will be identical to the Aligned?.out files.


Code options:

--qlib n	[default 1]  : specifies the format used in each .sam file for base quality: 1 indicates Illumina 1.3+ PHRED values; 2 indicates Illumina 1.8+ PHRED values

--pairs p	[default 0]  : Flag to indicate that the data contains paired-end read information

Explanation:

Where paired-end reads cover distinct regions of the genome, Samfire exploits this information to create joined reads, which span an increased proportion of the genome; such information is often valuable in providing information about partial haplotype composition.  Setting p to 1 instructs Samfire to search for, and join, paired reads.

--pl p		[default 3]  : Position of first element in the .sam file QNAME (comma delimited) containing paired-end information (counting from zero).
	
Explanation:

Information in the QNAME string is stored in a colon-separated manner, however the position of the paired-end information is not always consistent.  The --pl flag allows the position from which data is retrieved to be specified.

Example 1:
M00196:19:000000000-A0YC6:1:4:8245:10740_1:N:0:37
M00196:19:000000000-A0YC6:1:4:8245:10740_2:N:0:37

Here the paired-end read data are given as 8245:10740_1 and 8245:10740_2.  Setting --pl 5 recovers these data.

Example 2:
N73019:2:13:10640:17123#ACAGTG

Here the values 10640:17123 specify paired-end reads.  The default value --pl 3 recovers these.

--min_qual q	[default 30] : Minimum required PHRED score to include a nucleotide

Explanation: Base quality is used in two distinct ways in Samfire.  Firstly, reads are processed according to their median read quality.  Where the median read quality is lower than the minimum quality, the read is trimmed from either end until either the median quality is sufficiently high, or the sequence is discarded.  The end from which the least nucleotides need to be trimmed to achieve this is chosen for trimming.  Secondly, individual nucleotides in a read are processed according to the minimum quality.  Nucleotides not reaching the required quality are discarded.

--ali_qual q	[default 30] : Minimum required alignment quality.  Details of alignment quality are read from the .sam file

--ali_inc a		[default 0]  : Flag to include all sequences with unknown alignment quality.  Necessary where no alignment quality is recorded.

--min_rlen l	[default 30] : Minimum number of alleles that have to be reported by a short read for the read to be included in the output

Explanation: Once all other filtering steps are completed, the number of alleles reported by a read is measured against a minimum value, removing reads that report very few high quality nucleotides.

--ddup d	[default -1] : Removes duplicate sequences, namely reads spanning identical loci (in both paired reads where appropriate) with no more than d differences in sequence

--alm a	[default 1] : Setting a=2 ignores the alignment contained within the .sam file and realigns sequences using a simple comparison of nucleotide strings

--in i	[default Input_files.in] : File name from which the input .sam file names are read

--sorted s	[default 0] : Setting s=1 flags that the .sam files have been sorted in order of QNAME.  This is exploited to give a more rapid identification of paired reads.  Setting s=0 guarantees that paired reads will be identified, initiating a comprehensive, albeit slower, search.

--verb v	[default 0]  : Flag to produce more verbose output



./samfire sl_traj --ref <ref> [options]: calls single locus trajectories from previously filtered data

Required inputs:

Files Joined?.out generated by the previous program.

Times.in: A file containing time points (in integer format) at which samples were collected.  The default unit of time is days.

e.g. 

1
3
4

Outputs:

Single_locus_trajectories.out is a list of identified polymorphisms recorded over time.  The format of this is:

Locus Consensus Variant #Time-points [#Time #A #C #G #T N] where the time of sampling, the numbers of each identifiable nucleotide (only A, C, G, and T are considered), and the total number of samples are recorded for each time point for each trajectory.

For example:

312 A T 3 0 312 0 1 23 326 2 273 1 0 56 330 4 513 1 2 452 968

would indicate an A to T mutation at locus 312, that data was observed at three time points (i.e. 0, 2, 4), and gives details of the numbers of A,C,G, and T alleles observed at each time point: (312, 0, 1, 23), (273, 1 0 56), and (513, 1, 2, 452) respectively.  The total number of reads at each time point is specified as 326, 330, and 968 respectively.

By default, trajectories are shown for all loci at which a polymorphism was identified for at least one time point.

Variants?.out is a list of variants identified at every locus for each time point.  The format is:

Locus #A #C #G #T

Variants are shown in this file irrespective of whether they are formally called as polymorphisms.

Code options:

--q_cut q [default 0.01] : minimum frequency at which a polymorphism is identified.

Explanation:

Each polymorphism has to reach the frequency q at a given time point in order to be called.

--min_qual q [default 30] : Minimum required PHRED score to include a nucleotide

Explanation: 

The minimum PHRED score is utilised within the calling of SNPs (see next option).  If a non-default value is specified in the filtering, the same value should be applied in the calling of single locus trajectories.

--qp_cut p [default 0.001] : minimum probability at which a polymorphism is identified

Explanation: 

A probabilistic measure is applied in order to filter polymorphisms.  A base call error probability pr_e is calculated as pr_e = 10^-q, where q is the minimum base quality score used in filtering the data.  Now, supposing that at a given time point n out of a total of N bases are called to have a specific nucleotide, the probability that this occurred as a result of errors in base calling can be modelled as pr_E = \sum_{i=n}^{N} Binom(N,n,pr_e) where Binom is the binomial distribution.  A polymorphism is only called if pr_E < p.

--nmin n [default 10] : minimum number of variant counts required to call a polymorphism

Explanation: 

In a given time point, a polymorphism can only be called if at least n copies of a given variant are found in the filtered data

--repq r [default 1] : Number of time-point at which a polymorphism has to be called in order to call a trajectory.

Explanation: 

By default, identifying a polymorphism at a single time-point is sufficient to record a trajectory, which reports the number of variant alleles found at each time point.   This can be altered, reporting only trajectories that are found to show polymorphism at multiple time points.

--in i	[default Input_files.in] : File name from which the input .sam file names are read

--out o 	[default Single_locus_trajectories.out] : File name to which identified trajectories are saved

--verb v	[default 0]  : Flag to produce more verbose output



./samfire sl_noise : estimates the extent of noise in a dataset of single locus trajectories, in addition to that caused by the finite sampling of data

Noise is estimated via a two-step process.  Firstly, trajectories are filtered, removing trajectories which exhibit large changes in frequency and specific observations of polymorphisms that do not meet a frequency threshold.  Secondly, under the assumption that the remaining trajectories change in frequency only as the result of noise, a Dirichlet multinomial distribution is fitted to the data, the distribution parameter C describing the extent of overdispersion of the distribution.  Here, a value of C=0 is equivalent to a uniform distribution; as C tends to infinity, the distribution converges to a standard multinomial distribution.

Required inputs:

File Single_locus_trajectories.out generated by the previous program

Outputs:

Csl.out: A file containing the value of the overdispersion parameter C.

Code options:

--q_cut q 		[default 0.01] : minimum frequency at which a polymorphism is identified.

Here this parameter is used to filter out low-frequency polymorphisms.  

--dq_cut q	[default 0.05] : mean change in allele frequency per day above which a trajectory is excluded from noise calculations

This cutoff is conservative in nature, designed for application to data from within-host influenza infection, where the time-scale of infection is only a few days.

--seed s 	[default machine time] : integer seed to initialise the random number generator

--in i		[default Single_locus_trajectories.out] : File name from which trajectories are read

--verb v 	[default 0] : Flag to produce more verbose output.


./samfire sl_neutrality: Identifies loci at which trajectories evolve in a potentially non-neutral manner

Identifies trajectories that deviate significantly from a deterministic model of neutral evolution.  Models of evolution under neutrality, constant selection, and time-dependent selection are fitted to the data using a likelihood model informed by the characterisation of noise in observed allele frequencies.  The Bayesian Information Criterion is then applied to identify the best fitting model for each trajectory.  Trajectories for which a non-neutral model produces the best fit to the data are then reported.

Required inputs:

File Single_locus_trajectories.out, describing reads at identified polymorphic loci over time.

File Csl.out, containing an overdispersion parameter C, derived above.  If this file is not provided, the parameter C can be specified in the options to the code.

Outputs:

Single_locus_trajectories_sl.out, a filtered list of trajectories at which a non-neutral model produced the best fit to the observed single-locus data

Code options:

--specify_csl [By default this is read from Csl.out] : Overdispersion parameter C

--det [default 0] : Specifies the model used for fitting allele frequencies.  0 fits a model designed for within-host influenza evolution, with two rounds of mutation and selection per day.  1 fits a generic model of selection for a specific nucleotide with no selection; this is theoretically faster to calculate given trajectories of several days in length.

--in i	[default Single_locus_trajectories.out] : File name from which trajectories are read

--out o 	[default Single_locus_trajectories_sl.out] : File name to which identified trajectories are saved

--verb v [default 0] : Flag to produce more verbose output.



./samfire call_ml: Calls multi-locus trajectories from short data and single-locus trajectories

In this step, the content of each short read is considered across a set of loci at which polymorphism has been identified.  Short reads are processed into sets of partial haplotype information, according to which of these loci they span.

Required inputs:

Joined?.out, a list of joined and aligned sequences created above

Single_locus_trajectories_sl.out is used by default to identify a list of loci at which to call multi-locus haplotypes.

Times.in a list of times at which samples were collected

Outputs:

Output is given in one of two sets of formats.  By default a single output file is generated:

Multi_locus_trajectories.out is a list of multi-locus variants, recorded over time.  The format is:

#Loci [Set of loci] Variant #times [Time #observations]

where the set of loci specifies at which loci information was available from the dataset.  Within this file, longer partial haplotypes take priority when reporting variants; for example reports of the variant AG at loci 1 and 2 are not included in the reports of A at position 1 or of G at position 2.

Options may be specified to output a set of files expressing inferred full haplotypes, with partial haplotype data loci spanned separated, as follows:

Files Hap_data?.dat contain time-resolved, partial haplotype information.  The format is

P1 	n11 	n12 	n13
P2 	n21	n22	n23

Where P? is a partial haplotype, namely a string of alleles at a particular subset, and n?? are the number of each observed partial haplotype at each point in time.

Files Loci?.dat specify which loci are spanned by partial haplotypes in each set Hap_data?.dat.

Where full haplotypes are constructed, or imported, and partial haplotypes called against these, further files may be output:

Files Contribs?.dat specify which of the full haplotypes each of the partial haplotypes contributes to.  For example, if the derived full haplotypes are {ACG, ACT, AGG, TAA}, the partial haplotype AC- spanning the first two loci would contribute to the first two haplotypes.  For the purpose of an inference calculation performed at the level of full haplotypes, the frequency of a given partial haplotype is equal to the sum of the frequencies of the full haplotypes that contribute to that partial haplotype.

File Haps1.dat contains a list of inferred full haplotypes, spanning all of the loci identified by the code.

Code options:

--in i	[default Single_locus_trajectories_sl.out] : File name from which trajectories are read.  Note: Only the first column of this file is used by this part of the code.

--hap_q_cut [default 0.01] : Specifies the minimum frequency within a subset at which a multi-locus variant is called.

Explanation: Partial haplotypes are considered as subsets, according to which loci they span.  Within each subset, partial haplotype variants are filtered, removing observations that fall below a given frequency threshold within a set, and those with a low absolute number of observations.

--hap_n_min [default 10] : Specifies the minimum number of times a partial haplotype must be observed for it to be included in the set of called partial haplotypes.

--full_rep f [default 1] : Flag to specify output type.  By default, the file Multi_locus_trajectories.out is produced.  Setting f to 0 divides the output into multiple files according to which loci are spanned by each set of partial haplotypes, and whether full haplotypes are called or read in. 

--full_haps f [default 0] : Flag to call full haplotypes.  Setting f to 1 reports multi-locus variants by reference to a set of full haplotypes

Explanation: Potential haplotypes may be called against a set of full haplotypes.  Where partial haplotypes are called against full haplotypes, unfiltered partial haplotypes are reported if they could be emitted from a sequence contained in the full haplotype list.

--readhap [default 0] : Relevant when full haplotypes are specified.  Flag to read in full haplotypes from an external file.  If r is equal to 0, full haplotypes are constructed from the partial haplotype data.  If r equals 1, partial haplotypes are called against a list of haplotypes contained in an external file.

Explanation: Full haplotypes may be constructed from the partial haplotypes, as follows:

1.  All observed, filtered partial haplotypes are collected for each time point, being divided into sets according to the loci spanned by each partial haplotype.

2.  Partial haplotypes are combined to form a set of full haplotypes from which the filtered partial haplotypes could potentially have been emitted.  A set of rules is applied to do this:
	i) 	Loci at which only a single nucleotide is observed are identified
	ii) 	Partial haplotypes which are contained within other partial haplotypes are removed.  For example, -A-- is contained in CAT-.
	iii) 	Overlapping partial haplotypes are combined where they are in agreement.  For example CA-- and -AT- are combined to form CAT-.
	iv)	Combinations of partial haplotypes are created where they are not in agreement.  For example CA-- and -GT- are combined to give CAT- and CGT-.
	v)	Sets of possible spanning haplotypes across partial haplotypes are created.  For example, CA--, CG--, --TC and --CA give CATC, CGTC, CACA and CGCA.

3. Multi-locus variants are called against the inferred full haplotypes.  Observed partial haplotypes that do not match to full haplotypes are collected under the nominal partial haplotype X within each set.

--hap_file [default ../Haps1.dat] : Specifies the name of the file containing full haplotypes, if partial haplotypes are called against an externally provided list.

Full haplotypes are specified as a list of nucleotides at the given set of loci e.g.

ACTA
AGTA
AGCA
TGCA
etc.



./samfire ml_noise: Estimates the extent of noise in the multi-locus data, in addition to that caused by the finite sampling of data

The extent of noise in the multi-locus partial haplotype frequency data may not equal that in the single-locus allele frequency data.  A similar process to that used to estimate noise in the single-locus trajectories is applied to the multi-locus data.  Before calculating an estimate of noise, time-points at which there is no significant polymorphism, and partial haplotypes which do not exhibit significant polymorphism, are filtered from the dataset.  Both a maximum likelihood estimate of the noise parameter C, and a conservative estimate of this parameter, are reported.

Required inputs:

Partial_hap_list.dat: A list of files in the format of the Hap_data?.dat files, describing multi-locus trajectories.

Partial_hap_times.dat: A list of files containing the time points at which data was recorded for the partial haplotype sets described in the files listed in Multi_locus_trajectories.out.  Note that calculations can be performed across trajectories from sets of data that were collected at different time points.  The format of each file is the same as that of the file Times.in, used in in the sl_traj code above.  One time file should be specified for each multi-locus trajectory file.

Outputs:

Cml.out : Optimal value of C_ml

Conservative value of C_ml : Upon identifying the maximum likelihood value of C for the multi-locus dataset, a more conservative (i.e. lower) value of C is optionally calculated by finding that value of C corresponding to a likelihood 5 log units below the maximum inferred value.  This is recorded in Cml_cons.out.

Code options:

--full_haps f [default 0] : Flag to use data resulting from calling full haplotypes

--hap_q_cut [default 0.01] : Specifies the minimum frequency within a subset at which a multi-locus variant is included in the noise calculation.  This should be consistent with values for this parameter used earlier in the inference process.

--dq_cut q [default 0.05] : mean change in partial haplotype frequency per day above which a trajectory is excluded from the noise calculations

--conservative c [default 0] : Specifying c=1 returns a more conservative estimate of C

--verb v [default 0] : Flag to produce more verbose output.



./samfire ld_calc: Estimates linkage disequilibrium between alleles in a dataset

Required inputs:

Single_locus_trajectories_sl.out : A list of single-locus trajectories, the loci of which will be used to calculate linkage disequilibrium.

Cml.out : Optimal value of C_ml

Outputs:

LD_info.out: Information is reported in the format:
Time point	Locus i	Locus j	D_ij	D’
where D_ij is the linkage disequilibrium between alleles, and D’ is the normalised linkage disequilibrium.

Code options:

--verb v [default 0] : Flag to produce more verbose output.

