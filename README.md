# samfire
SAMFIRE carries out multi-locus variant calling for time-resolved sequence data.  Instructions for using the code follow.

Compiling the code

Samfire is code for calling multi-locus variants from time-resolved sequence data.  It requires prior installation of the GSL library, and is compiled with the command line

make samfire

Running the code

Samfire processes data via a series of subprograms:

./samfire filter carries out QC on a set of data, returning a list of filtered readsm
./samfire sl_traj calls single locus trajectories from the filtered data
./samfire sl_noise estimates the extent of noise in the single-locus trajectory data, in addition to that caused by the finite sampling of data
./samfire sl_neutrality identifies trajectories that deviate significantly from a pattern of evolution under neutrality
./samfire call_ml calls multi-locus trajectories
./samfire ml_noise estimates the extent of noise in the multi-locus trajectory data, in addition to that caused by the finite sampling of data

Details of the individual programs follows:
  

./samfire filter —ref <ref> [options]: carries out QC on a set of data, returning a list of filtered reads.  Here <ref> is the name of a reference file in fasta format.

Required inputs:

Input_files.in: A list of files in .sam format, ordered by the time at which they were collected.  In each case the full path to the file (relative to the directory in which samfire is executed) should be specified one file per line as follows:

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

--qlib n	[default 2]  : specifies the format used in each .sam file for base quality: 1 indicates Illumina 1.3+ PHRED values; 2 indicates Illumina 1.8+ PHRED values; 3 indicates Oxford Nanopore data

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

--trim t	[default 0]: Trims the ends of each short read by an amount t nucleotides.  Can be used in cases where it is suspected that the ends of reads contain spurious variation.

--sorted s	[default 0] : Setting s=1 flags that the .sam files have been sorted in order of QNAME.  This is exploited to give a more rapid identification of paired reads.  Setting s=0 guarantees that paired reads will be identified, initiating a comprehensive, albeit slower, search.

--verb v	[default 0]  : Flag to produce more verbose output

--printali a	[default 0] : Setting a=1 prints out the files Aligned?.out.  By default just the Joined?.out files are printed.

--printqual q [default 0] : Setting a=1 prints out Quality?.out files, giving the string of read quality values from the .sam file.  These match the Joined?.out data; i.e. these are post-processing of CIGAR data.  Setting a=2 outputs the same data, but processed into space-separated numerical PHRED scores.


./samfire sl_traj --ref <ref> [options]: calls single locus trajectories from previously filtered data  

Required inputs:

Files Joined?.out generated by the previous program.

A reference sequence against which to call polymorphisms, in .fasta format.

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

--repq r [default 1] : Number of time-points at which a polymorphism has to be called in order to call a trajectory.

Explanation: 

By default, identifying a polymorphism at a single time-point is sufficient to record a trajectory, which reports the number of variant alleles found at each time point.   This can be altered, reporting only trajectories that are found to show polymorphism at multiple time points.

—first f [default 0] : Reports only polymorphisms which are found to be polymorphic in the data from the first sample.

--vs_ref v [default 0] : Flag to call polymorphisms with respect to the reference sequence.  The initial allele is defined according to the specified consensus sequence, rather than the majority allele in the first sample.

--gmaf g [default 0] : Flag to report fixation events relative to the consensus sequence as polymorphisms.  Used with vs_ref 1.

--dep_cut q	[default 1] : minimum read depth for a variant to be called.

--in i	[default Input_files.in] : File name from which the input .sam file names are read

--nosam n [default 0] : Run the code with no .sam file list.  The value of nosam designates the number of Joined?.out files to read.

--out o 	[default Single_locus_trajectories.out] : File name to which identified trajectories are saved

--verb v	[default 0]  : Flag to produce more verbose output

--uniq v	[default 0]  : Only output one trajectory per locus (in the case where there are multiple alles at a single locus)



./samfire consensus : Calculates the consensus sequence at each time point and across time points

The Variants?.out files created by the sl_traj code are used to calculate a sequence consensus.  The across-time consensus for any locus is calculated as the nucleotide reported at that locus by the greatest total number of reads.

A consensus sequence is reported for each Variants?.out file; and is given in the file Consensus?.fa.  A total consensus is given in the file Consensus_all.fa.  Output files are in fasta format.

Code options:

--translate t	[default 0] : Setting the flag to 1 translates the consensus nucleotide sequences into amino acid sequences.  Translation does not require a start codon to commence but ceases one a stop codon is reached.  Amino acid consensus sequences are reported in the files Consensus?_AA.fa.

--trans_start s 	[default 1]: Sets the point at which SAMFIRE begins translating the sequence.  The default value of 1 begins translation at the start of the nucleotide sequence.

--repair_consensus r [default 1]: Flag to repair nucleotide sequences before translation to an amino acid sequence.  Ambiguous nucleotides are reported as ’N’ in the consensus and cause problems for translation to protein sequence.  Setting r=1 uses the global nucleotide consensus sequence to repair such ambiguous nucleotides.  Setting r=2 attempts to use previous consensus sequences to repair these nucleotides, resorting to the global consensus if no prior consensus nucleotide can be found.
	
--get_variants g [default 0]: Flag to generate Variant_mask files.  The amino acid consensus is used to characterise the type of change to the consensus made by each nucleotide at each position.  These are reported as 0 (no change), 1 (synonymous change), 2 (non-synonymous change), 3 (nonsense change) in the files Variant_mask?.out.

Note: The perl script combine_masks.pl can be used to combine Variant_mask files contained within two subdirectories: Use e.g. perl combine_masks.pl -a dir1 -b dir2.  This is designed for genomic data with multiple or overlapping read frames.


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

--dep_cut q	[default 1] : minimum read depth for a variant to be considered in the noise calculation

--seed s 	[default machine time] : integer seed to initialise the random number generator

--in i		[default Single_locus_trajectories.out] : File name from which trajectories are read

--verb v 	[default 0] : Flag to produce more verbose output.


./samfire ef_depth : Calculates the effective depth of sequencing at each locus

Takes as input the file Csl.out and Variants?.out files.  Produces files Depths?.out


./samfire sl_neutrality: Identifies loci at which trajectories evolve in a potentially non-neutral manner

Identifies trajectories that deviate significantly from a deterministic model of neutral evolution.  Models of evolution under neutrality, constant selection, and time-dependent selection are fitted to the data using a likelihood model informed by the characterisation of noise in observed allele frequencies.  The Bayesian Information Criterion is then applied to identify the best fitting model for each trajectory.  Trajectories for which a non-neutral model produces the best fit to the data are then reported.

Note: As with any optimisation algorithm, this approach is not guaranteed to find a global maximum in every case.  Some tweaking of the code may be required for any specific set of data.  A good place to start with this is the SetInitialSelection routine, part of optimisation.cpp.  Altering the range in which initial estimates are generated can sometimes improve the inference.

Required inputs:

File Single_locus_trajectories.out, describing reads at identified polymorphic loci over time.

File Csl.out, containing an overdispersion parameter C, derived above.  If this file is not provided, the parameter C can be specified in the options to the code.

Outputs:

Single_locus_trajectories_sl.out, a filtered list of trajectories at which a non-neutral model produced the best fit to the observed single-locus data

Code options:

--specify_csl [By default this is read from Csl.out] : Overdispersion parameter C

--detprop [default 0] : Specifies the model used for fitting allele frequencies.  0 fits a model designed for within-host influenza evolution, with two rounds of mutation and selection per day.  1 fits a generic model of selection for a specific nucleotide with no selection; this is theoretically faster to calculate given trajectories of several days in length.

--in i	[default Single_locus_trajectories.out] : File name from which trajectories are read

--out o 	[default Single_locus_trajectories_sl.out] : File name to which identified trajectories are saved

--verb v [default 0] : Flag to produce more verbose output.

--printnl p [default 0] : Flag to output trajectories which do not appear to be under selection.  Trajectories are printed in Single_locus_trajectories_nl.out.


./samfire contamination: Checks for contamination between two samples

This routine compares short read data from datasets to distinct reference sequences, looking for potential mixing of reads between samples.

 Required inputs:

Two files of the form Joined?.out.  These must be specified by the flags --jn1 and --jn2.

Two reference sequences in fasta file format.  These must be specified by the flags --ref1 and --ref2.  The file given by ref1 should be the consensus sequence of the data from the jn1 flag, while the file given by ref2 should be the consensus sequence of the data from the jn2 flag.  See elsewhere in the instructions for getting consensus sequences from data in SAMFIRE format.

By default the two consensus sequences are compared, identifying loci at which they possess different alleles.  Individual reads are then compared to the consensus sequences at points where they overlap these loci.

The manner in which the routine works is to scan individual sequence reads, as processed by the filter routine above.  For reads which span at least two variant sites, the number of sequence mismatches between each read and and each of the two reference sequences is calculated.  This is conducted both for the variant sites and across all sites.  Reads which are more similar to the ‘wrong’ reference, measured by an increased total number of nucleotide differences, are reported.

Outputs:

Two files are outputted, Flagged_reads1.out and Flagged_reads2.out.  Each file contains the index number of each potentially mismatched read, along with the number of reverse variant sites spanned, the total number of sites spanned (and for which the read reports a variant), and the number of differences with each reference at each of the two sets of sites.

As an option, variants of the Joined files can be outputted, from which flagged reads have been removed.

Code options:

--in i	: File name from which trajectories are read.  By default the code compares the two 

--ref1: File containing first reference sequence

--ref2: File containing second reference sequence

--jn1: File containing first reference sequence

--jn2: File containing second reference sequence

--dec: Output Joined files

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

--hap_n_min [default 0] : Specifies the minimum number of times a partial haplotype must be observed for it to be included in the set of called partial haplotypes.

--multi_gap m [default 0] : By default, multi-locus variants are called as reads containing information for sets of continuous variant loci.  This option reverts to only calling full haplotypes, including indels.  It is suitable for data where the sequences span the entirety of a given region.

--maxgap m [default 1] : Works with multi_gap.  Allows for haplotypes with up to m indels.  Currently implemented for m<=2.

--full_rep f [default 1] : Flag to specify output type.  By default, the file Multi_locus_trajectories.out is produced.  Setting f to 0 divides the output into multiple files according to which loci are spanned by each set of partial haplotypes, and whether full haplotypes are called or read in. 

--printx p [default 0] : Flag to print X haplotype data in Multi_locus_trajectories.out.  These are reads which do not match any of the identified haplotypes.  Sometimes these lines of the file were annoying so I got rid of them by default.

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


./samfire stats: Calculate some statistics from the data

This routine requires at least one more input to do something

Code options:

--calc_lengths [default 0]: Flag to calculate the lengths of joined read sequences.  Data are output in histogram format in files Seq_lengths?.out

--calc_qual [default 0]: Flag to calculate the distribution of PHRED scores.  Requires the prior calculations of Quality?.out files, done using the filter option above.  Data are output in the Base_quality?.out files.

--calc_pi [default 0]: Flag to calculate the pi statistic of population diversity.  Produces output from each sample in the file Pi_diversity.dat.  Bootstrap replicate samples are given in the file Pi_bootstraps.dat if these are calculated.

	--sns[default 0]: Flag to recalculate pi over four-fold synonymous and non-synonymous sites.  Requires input Variant_types.dat from variant_types option.

--calc_var_comp [default 0]: Flag to calculate the proportion of variants in each of the twelve mutational classes A to C, A to G, A to T, etc.  By default the routine calculates this statistic for low-frequency variation, that is, under the cutoff specified by the option q_cut.  The fraction of variation in each class is reported in the file Variant_composition.dat.  Bootstrap replicate samples are given in the file Variant_composition_bootstrap.dat

Produces output from each sample in the file Pi_diversity.dat.  Bootstrap replicate samples are given in the file Pi_bootstraps.dat if these are calculated.

--calc_ham_cons [default 0]: Flag to calculate Hamming distances between the consensus sequences at each timepoint.  Requires consensus sequences to have been found using the consensus command.  Output is to the file Consensus_Hamming_Distances.dat

--calc_ham_var [default 0]: Flag to calculate Hamming distances between the variant allele frequency data at each timepoint.  Requires consensus sequences to have been found using the consensus command.  Output is to the file Variant_Hamming_Distances.dat

--calc_tot_var [default 0]: Flag to calculate sums of variant frequencies across the whole genome.  Output is to the file Total_Frequency_Stats.dat.  Allows for bootstrapping.

	--sns[default 0]: Flag (if set to 1) to recalculate total variant frequencies over four-fold synonymous and non-synonymous sites.  Requires input Variant_types.dat from variant_types option.  Setting this parameter to 2 calculates total variant frequencies over all S and NS sites in the genome.  Outputs to Total_Frequency_Stats_S.dat and Total_Frequency_Stats_NS.dat

--bootstrap [default 0]: Flag to calculate bootstrap version of the Variants?.out files via multinomial sampling of the original with the observed allele frequencies.  Applies to the calculations of diversity and variant composition.  100 bootstrap samples are calculated; statistics are then calculated for each bootstrap population.

./samfire reading_frames: Identify reading frames within the consensus sequence

This reads through an input consensus sequence (by default that at the first time point), and identifies reading frames within the data.  For translated protein sequences longer than a given threshold the code outputs the nucleotide positions of the reading frame, the nucleotides involved, and the translated sequence.  Data is output to the file Reading_Frames.dat.  The file ORFs.dat contains simply the positions of the reading frames.  Note: This currently assumes that proteins are translated from continuous stretches of nucleotide sequence.

Code options:

--len [default 30] : Sets the minimum number of amino acids for a potential reading frame to be reported.

--in [default Consensus0.fa] : Sets the file which is read in as input to generate the translations.


./samfire variant_types: Identifies whether variants are synonymous, non-synonymous, or fall in regions between reading frames.  This requires a corrected version of ORFs.dat; this can be generated, for example, by searching for translated proteins using BLAST.

This reads through the sequence and identifies whether substitutions are each position are consensus (0), synonymous (1), non-synonymous (2), or in a non-coding region (3).  These data are output in the file Variant_types.dat.  The code outputs to command line whether variants in single locus trajectories are S, NS, or NS as above.

--in [default Single_locus_trajectories.out] : Sets the file which is read in as input for the single-locus variants.



./samfire calc_distances: Calculate distances between Variants files

This routine calculates the distance between populations at different times using information from the Variants?.out files.  Given the nucleotide counts {n_A,n_C,n_G,n_T} at a given locus, these are converted to frequencies.  Comparing two loci can then be achieved using the metric:

d = [ abs (q^1_A - q^2_A) + abs (q^1_C - q^2_C) + abs (q^1_G - q^2_G) + abs (q^1_T - q^2_T) ] / 2

where division by two allows for the differences being counted twice in this sum.  Summing up these values across all nucleotide sites gives the total distance between populations; this is a variant of the Hamming distance which accounts for low frequency polymorphism in a population.

Distances are reported in the matrix Distances.out.  The number of nucleotide positions used in the calculation is given for each point of the matrix in the file Counts.out; positions in the genome for which no nucleotides were observed are excluded from the calculation.

Code options:

--sns_distances s [default 0] : Flag to calculate separate synonymous and non-synonymous distances.  Setting s=1 activates the flag.  This option requires the previous generation of Variant_mask?.out files, achieved under the consensus routine, detailed above.  Here, calculations are performed identifying which nucleotides represent synonymous and non-synonymous changes to the population, using the consensus sequence of the first sequence in the calculation to define this (i.e. the distance function is not necessarily commutative).  Outputs are given in the matrices Distances_n.out and Distances_s.out.  The number of nucleotide positions in each calculation are given in the files Counts_n.out and Counts_s.out.  Here the metric for d is not divided by two; the sum is calculated only over nucleotides which represent a change in the consensus.

