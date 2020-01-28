//Shared information for linked optimisation
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <vector>
#include <string>
#include <algorithm>
#include <map>

using namespace std;

struct cutoff {
	int time;
	vector<double> freq;
	vector<double> top;
	vector<double> bottom;
};

//Call SNPs
void CountNucleotides (run_params p, rseq refseq, vector<string> sam_files, vector< vector<joined> > t_reads, vector<nuc>& ref_counts);
void SetupRefCounts (rseq refseq, vector<string> sam_files, vector<nuc>& ref_counts);
void PerSequenceSNPs (run_params p, rseq refseq, vector<string> sam_files);
void CountNucs (run_params p, rseq refseq, vector<joined> t_read, nuc& r_count);
void SetupRCount (rseq refseq, nuc& r_count);

void CallPolymorphisms (run_params p, rseq refseq, vector<nuc> ref_counts, vector<poly>& polys);
void ProcessP (run_params p, int j, int t, vector<double> freqs, vector<char> nuc, char cons, vector<int> counts, int N, vector<poly>& polys);
void SpreadFirst(vector<poly>& polys);
void ConstructSLTrajs (run_params p, vector<poly> polys, vector<nuc> ref_counts, vector<str>& sltrajs);
void SLTFreqs (run_params p, vector<str>& sltrajs);
void SLTMeanFreqs (run_params p, vector<str>& sltrajs);
void FilterSLTrajs (run_params p, vector<str>& sltrajs);
void FilterSLTrajs2 (run_params p, vector<str>& sltrajs);
void PotentialNonNeutrality (run_params p, double Csl_opt, vector<str>& sltrajs, vector<double> fact_store, gsl_rng *rgen);
void DeleteDuplicatePolymporphisms (run_params p, vector<poly>& poly);
void CallPolymorphismsVsRef (run_params p, rseq refseq, vector<nuc> ref_counts, vector<poly>& polys);
void ProcessPVsRef (run_params p, int j, vector<double> freqs, vector<char> nuc, char cons, vector<int> counts, int N, vector<poly>& polys);
