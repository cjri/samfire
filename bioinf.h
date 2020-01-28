//Shared information for linked optimisation
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <map>

using namespace std;

char TranslateCodon (vector<char> c, vector< vector<char> > cod, vector<char> aa);
char FastTranslateCodon (vector<char> c, vector<char> aa);
int FastTranslateCodonN (char c1, char c2, char c3);
char TranslateNAA (int n,vector<char> aa);
void RepairConsensusSequence (run_params p, int i, vector< vector<char> >& all_cons, vector<char>& global_cons);
void GetVarTypes (vector<char> c, char s, vector< vector<char> >& cod, vector<char>& aa, vector< vector<int> >& s_types);
void CalculateVariantFrequencies (vector< vector< vector<int> > >& all_nucs, vector< vector< vector<double> > >& all_freqs);
void CalculateVariantDistances (vector< vector< vector<double> > >& all_freqs, vector< vector<double> >& distances, vector< vector<double> >& all_counts);
void CalculateVariantDistancesSNS (vector< vector< vector<double> > >& all_freqs, vector< vector< vector<int> > >& all_types, vector< vector<double> >& distances_s, vector< vector<double> >& distances_n, vector< vector<double> >& all_counts_s, vector< vector<double> >& all_counts_n);
void Find_Reading_Frames (run_params p);
void Find_Variant_Types (run_params p);
void PrintType(int locus, int i, vector< vector<int> > types);
void process_nuc(int pos, char orig_aa, vector<char> codon, vector<char> nuc, vector<char>& aa, vector<int>& replace);
void add_to_replace(int pos, char orig_aa, char v, vector<char> codon, vector<char> nuc, vector<char>& aa, vector<int>& replace);
void aa_compare(char new_aa, char orig_aa, vector<int>& replace);




