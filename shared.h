#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <vector>
#include <string>
#include <algorithm>
#include <map>

using namespace std;

struct refseq {
	string name;
	string seq;
	int size;
};

struct rd {
	int alpos;
	int inc;
	int rev;
	int refno;
	string ref;
	string seq;
	string revseq;
	string qual;
};

struct pr {
	string s1;
	string s2;
};

struct par {
	int i1;
	int i2;
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

//utilities
void GetRefSeqs (ifstream& ref_file, vector<refseq>& refs);
void ReadSamFile (ifstream& sam_file, int& s_length, vector<rd>& data, vector<refseq> refs);
void ReadSamFile2 (ifstream& sam_file, int& s_length, vector<pr>& data);
void DoForwardAlignment (int i, int min_qual, int max_qual, int& firstpos, int&isok, vector<char> qual, vector<refseq> refs, vector<rd>& data);
void DoReverseAlignment (int i, int min_qual, int max_qual, int& firstpos, int&isok, vector<char> qual, vector<refseq> refs, vector<rd>& data);
void MeasureInsertions (int& ins, string aliseq);
void OverlapStart (int& refin_s, string refseq);
void OverlapEnd (int& refin_e, string refseq);
void CheckOK (int ins, int& isok, string star, string seq);
void FindFirstPos (int& firstpos, int refin_s, string refseq, string aliseq);
void AlignSequences (int s_length, int min_qual, int scope, vector<rd>& data, vector<refseq> refs);
void makequal (vector<char>& qual);
int findqual2 (int i, int min_qual, int max_qual, vector<char> qual, vector<rd>& data);
int ScoreSim (string a, string b, int p, int rsize);
string RevTr (string a);
int GetMedian (int a, int b, vector<int> qvec);
int ScoreSim2 (string a, string b, int p, int rsize);
