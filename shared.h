//Shared information for linked optimisation
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <vector>
#include <string>
#include <algorithm>
#include <map>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multimin.h>

using namespace std;

struct run_params {
	string ref; //Reference sequence
	string ref1;  //References for contamination check
	string ref2;
	string jn1;  //Joined.out file names for contamination check
	string jn2;
	int plines;	  //Number of elements in string defining paired-end data in specific .sam file
    string delimit;  //Delimiter used for paired-end data.
	int almethod; //Method of gathering sequence alignment data
	int min_qual; //Minimum sequence quality
	int max_qual; //Maximum sequence quality
	int min_rlen; //Minimum number of alleles that must be contained in a single read
	int ali_qual; //Minimum alignment quality
    int print_qual; //Print quality data
    int print_ali; //Prints unjoined Aligned?.out files
	int ali_inc; //Flag to include reads with uncertain alignment quality
	int trim; //Amount to trim from each end of reads
	int sorted; //Flag to indicate that .sam files have been sorted by their first column
	int qlib;  //Base quality indicators
	int pairs; //Include paired end reads
	int ddup; //Remove duplicate reads
	double specify_csl; //Read value of noise parameter C_sl;
	double q_cut; //Cutoff to identify polymorphism
	double qp_cut; //Probability required to identify polymorphism
	int n_min; //Minimum number of variant alleles to call polymorphism
	int rep_q; //Number of time points at which a SNP needs to be observed in order for it to be recorded
	int first; //Flag to require first time-point to be polymorphic when calling SNPs
	double dq_cut; //Cutoff change in allele frequency per day for potential neutrality
	int dep_cut; //Cutoff overall read depth to call polymorphism
	double seed; //Random seed
	int skip;
	int printnl; //In sl_neutrality, print trajectories which are not inferred to be under selection
	int no_sam; //Skip reading in list of .sam files in calling single-locus variants
	int pos; //Default position for mutational scan
	int det; //Flag to use deterministic model of single-locus evolution; no mutation or drift
	double mu; //Mutation rate for single-locus evolution model
	double hap_q_cut;  //Minimum freuqency within a dataset at which to include a partial haplotype
	int hap_n_min;  //Minimum number of calls within a dataset required to include a partial haplotype
	int hap_index; //Flag to use the original haplotype numbers in numbering Hap_data?.dat, Contribs?.dat
	int conservative; //Flag to calculate additional conservative estimate of noise
	int readhap; //Flag to call multi-locus polymorphisms against specific haplotypes, contained in a file
	string hap_file; //File name for haplotype data
	string in_file; //File name for output
	int get_in; //Flag to use read input file
	string out_file; //File name for output
	int get_out; //Flag to use read output file
	int full_haps; //Flag to generate full haplotypes
	int printx; //Flag to print X haplotypes in Multi_locus_haplotypes
	int full_rep; //Flag to give Multi_locus_haplotypes output
	int vs_ref; //Call against the given reference sequence, rather than against the consensus in the first time point at each position
	int gmaf; //Call polymorphisms that are fixed against the reference sequence
	int uniq; //Only print one trajectory per locus
	int multi_gap; //Option in multilocus calling to call full haplotypes with gaps, rather than consecutive loci
	int maxgap; //Maximum number of gaps to allow in multi_gap code
	int translate; //Translate consensus nucleotide sequences to amino acids.  May require start codon to be specified
	int trans_start; //Start position for translation of consensus sequence.
	int get_variants; //Option for consensus sequence - generates mask of mutational types for each variant file conditional on the start position for translation
	int repair_consensus; //Option to repair consensus sequences which contain Ns - ambiguous nucleotides
	int sns_distances; //Option to calculate distances between variant files in terms of S and NS variants
	
	//Contamination
	int decon; //Option to decontaminate Joined files in contamination code
	
    //LD_simple
    int noopt; //Flag allows the printing of variant statistics (with --verb 1) while not doing the likelihood calculation
    
	int calc_len; //Print stats for read lengths
	int calc_qual; //Print stats for base qualities
	int calc_pi; //Print diversity statistics for datasets
	int sns; //Calculate diversity for four-fold synonymous and non-synonymous sites
	int calc_var_comp; //Calcuation composition of variation
	int calc_ham_cons; //Calculation of Hamming distances between consensus sequences
	int calc_ham_var; //Calculation of Hamming distances between variant data
	int calc_sum_var; //Calculation of the sum of variant frequencies
	int bootstrap; //Calculate bootstrap samples of Variant files for statistics
	int len; //Minimum length of protein to report in reading frames calculation
	int verb; //Verbose output
};

struct rseq {
	string name;
	string seq;
	int size;
};

struct rd {
	int alpos;
	int alq;
	int inc;
	int flag;
	vector<int> flagv;
	int rev;
	int refno;
	int pairno;
	int del;
	string ref;
	string cigar;
	string cig_string; //Cigar data by nucleotide
	string seq;
	string revseq;
	string qual;
	string revqual;
	string paircode;
	vector<int> rmqual;
};

struct alldat {
	int s_length;
	vector<rd> data;
};

struct poly {
	int locus;
	int first;
	char nuc;
	char cons;
};

struct mpoly {
	vector<char> vars;
	int count;
};

struct npoly {
	int cons;
	int locus;
	int cod;
	vector<int> count;
	vector<int> ccount;
	vector<int> tot;
};

struct str {  //Single-locus trajectory through time
	int locus;
	int inc;
	char cons;
	char nuc;
	vector<int> times;
	vector<int> nA;
	vector<int> nC;
	vector<int> nG;
	vector<int> nT;
	vector<int> nN;
	vector<double> qA;
	vector<double> qC;
	vector<double> qG;
	vector<double> qT;
	double mA;	//Mean frequencies - sum of observations over number of observations
	double mC;
	double mG;
	double mT;
	vector<double> logL;
	vector<double> BIC;
};

struct mtr { //Multi-locus trajectory through time
	vector<char> seq;
	vector<int> n; //Number of observations
	double m; //Mean frequency
	vector<double> logL;
};

struct joined {
	int alpos;
	string seq;
};

struct ql {
        int alpos;
        vector<int> seq;
};

struct pr {
	string s1;
	string s2;
};

struct par {
	int i1;
	int i2;
};

struct ld_info {
	int i;
	int j;
	vector<int> n_i1;
	vector<int> n_i0;
	vector<int> n_j1;
	vector<int> n_j0;
	vector<int> n_11;
	vector<int> n_10;
	vector<int> n_01;
	vector<int> n_00;
};

struct det {
	int start;
	int length;
	int used;
	int sam_ref;
	string seq;
};

struct nuc {
	vector<int> nA;
	vector<int> nC;
	vector<int> nG;
	vector<int> nT;
	vector<int> nN;
};



double DirichletMultiCalc (int N, double c, vector<int> obs, vector<double> inf, vector<double>& fact_store);


