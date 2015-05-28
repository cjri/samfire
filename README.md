# samfire
SAMFIRE (SAm FIle REalignment) processes and realigns data from a .sam file.  The code is at an early stage of development.

The code requires installation of the MUSCLE sequence alignment package.  MUSCLE is currently freely available from http://www.drive5.com/muscle/

Overview of the method: Indidual short reads from the .sam file are processed as follows:
i) Assessment of sequence read quality.
    The median read quality of each read is calculated.  If this value is below a default minimum of 30, the read is trimmed from the end until either the median read quality reaches the required minimum, or the read length falls below 30 nucleotides.  In the latter case, the read is discarded.
    Surviving reads are aligned to the reference file.  Nucleotides extending beyond the end of the reference are discarded.  The aligned sequence is retained if the sequence identity between the read and the reference is greater or equal to 95%, and for which no insertion or deletion events were inferred.
    Surviving reads are processed, replacing all nucleotides with base quality score less than 30 with blank '-' nucleotides.  Sequences with at least 30 non-blank nucleotides are retained, and reported.
    
- Compiling the code:
Compilation is achived with the command 'make samfire'.

- Running the code:
The executable line is of the form ./samfire [reference] [.sam file] 

where [reference] is a reference sequence for the alignment data in .fasta format, and [.sam file] comprises sequence data in .sam format.

- Results produced:
Three files are produced by the code.

Allele_frequencies.out: Describes the number of each nucleotide found at each locus in the reference sequence, in the format:
	Locus	#nucleotides	#A	#C	#G	#T

Included.out: Index of short reads that were not discarded in the sequence quality checking process.

Sequences.out: Processed read data, in the format:
	Initial locus	Sequence

Within the reported sequences, blank ‘-‘ nucleotides denote positions in each sequence for which the base quality was insufficient to unambiguously call a nucleotide.
 
