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
	string ref;
	int plines;	  //Number of elements in string defining paired-end data in specific .sam file
	int almethod; //Method of gathering sequence alignment data
	int min_qual; //Minimum sequence quality
	int max_qual; //Maximum sequence quality
	int min_rlen; //Minimum number of alleles that must be contained in a single read
	int ali_qual; //Minimum alignment quality
	int ali_inc; //Flag to include reads with uncertain alignment quality
	int sorted; //Flag to indicate that .sam files have been sorted by their first column
	int qlib;  //Base quality indicators
	int pairs; //Include paired end reads
	int ddup; //Remove duplicate reads
	double specify_csl; //Read value of noise parameter C_sl;
	double q_cut; //Cutoff to identify polymorphism
	double qp_cut; //Probability required to identify polymorphism
	int n_min; //Minimum number of variant alleles to call polymorphism
	int rep_q; //Number of time points at which a SNP needs to be observed in order for it to be recorded
	double dq_cut; //Cutoff change in allele frequency per day for potential neutrality
	double seed; //Random seed
	int skip;
	int no_sam; //Skip reading in list of .sam files in calling single-locus variants
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
	int full_rep; //Flag to give Multi_locus_haplotypes output
	int vs_ref; //Call against the given reference sequence, rather than against the consensus in the first time point at each position
	int gmaf; //Call polymorphisms that are fixed against the reference sequence
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
	
};

struct alldat {
	int s_length;
	vector<rd> data;
};

struct poly {
	int locus;
	char nuc;
	char cons;
};

struct mpoly {
	vector<char> vars;
	int count;
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


