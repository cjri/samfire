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

//Stats
void Calc_Seq_Lengths (vector<joined> t_read, vector<int>& l_dist);
void Calc_Base_Quality (vector<ql> q_read, vector<long>& q_dist);
void MakeBootstraps (vector< vector< vector<int> > >& all_nucs, vector< vector< vector< vector<int> > > >& bootstrap, gsl_rng *rgen);
void Calc_Pi_Diversity (run_params p, vector< vector< vector<int> > >& all_nucs, vector< vector< vector< vector<int> > > >& bootstrap);
void PiCalculation(vector< vector< vector<int> > >& all_nucs, vector<double>& pi);
void Construct_FourFold_Variants (vector< vector<int> >& var, vector< vector< vector<int> > >& all_nucs, vector< vector< vector<int> > >& all_nucs_s, vector< vector< vector<int> > >& all_nucs_ns);
void Construct_AllFreqsSNS_Variants (vector< vector<int> >& var, vector< vector< vector<int> > >& all_nucs, vector< vector< vector<int> > >& all_nucs_s, vector< vector< vector<int> > >& all_nucs_ns);
void Calc_Variant_Composition (run_params p, vector< vector< vector<int> > >& all_nucs, vector< vector<char> >& all_cons, vector< vector< vector< vector<int> > > >& bootstrap);
void CVC (run_params p, vector< vector< vector<int> > >& all_nucs, vector< vector<char> >& all_cons, ofstream& cv_file);
void Calculate_Consensus_Hamming_Distances (vector<string>& sam_files);
void Calculate_Variant_Hamming_Distances (vector< vector< vector<int> > >& all_nucs);
void Calc_Tot_Freqs (run_params p, vector< vector< vector<int> > >& all_nucs, vector< vector< vector< vector<int> > > >& bootstrap);
void TF_Calculation (vector< vector< vector<int> > >& all_nucs, vector< vector<double> > total_freqs, ofstream& tot_file);




