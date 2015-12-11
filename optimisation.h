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

void GetFactVectorSL (vector<str> sltrajs, vector<double>& fact_store);
void GetFactVectorML (vector< vector<mtr> > mltrajs, vector<double>& fact_store);
void FindLogFact(vector<double>& fact_store, int N);
void OptimiseCsl (run_params p, double& Csl_opt, vector<str> sltrajs, vector<double> fact_store, gsl_rng *rgen);
void OptimiseCml (run_params p, double& Cml_opt, vector<int> include_ml, vector< vector<mtr> > mltrajs, vector<double> fact_store, gsl_rng *rgen);
void ConservativeCml (run_params p, double& Cml_opt, vector<int> include_ml, vector< vector<mtr> > mltrajs, vector<double> fact_store, gsl_rng *rgen);
void GetInfMeanFreq(int traj, vector<str> sltrajs, vector<double>&init_freqs);
void GetInfStartFreq(int traj, vector<str> sltrajs, vector<double>&init_freqs);
double OptimiseSLTraj (int verb, run_params p, int maxrep, int traj, char sel_model, int tds, int n_times, vector<int> times, vector<int> N, double Csl_opt, vector<str> sltrajs, vector<double> fact_store, gsl_rng *rgen);
void GetObsVector(int traj, vector<str> sltrajs, vector< vector<int> >& obs);
void GetMutationMatrix (double mu, vector <vector <double> >& mut);
void SquareMutationMatrix (vector<vector<double> > &m );
void SetInitialFreqs (vector<double>& init_freqs, int n_haps, gsl_rng *rgen);
void NormaliseFreqs (vector<double>& init_freqs, int n_haps);
void SetInitialSelection (int tds, char sel_model, int n_times, vector<double>& sigs, gsl_rng *rgen);
void SigsHapSigs (int tds, char sel_model, int n_times, vector<double> sigs, vector< vector<double> >& hapsigs);
void Propagation(int verb, vector<int> times, vector<double> init_freqs, vector< vector<double> > mut, vector< vector<double> > hapsigs, vector< vector<double> >& inf);
void DeterministicPropagation(int verb, vector<int> times, vector<double> init_freqs, vector< vector<double> > hapsigs, vector< vector<double> >& inf);
void GetSI (int index, vector< vector<double> > hapsigs, vector<double>& si);
void GrowGeneration(vector<vector<double> >& m, vector<double>& q, vector<double>& sigma);
void MatMult (vector<vector<double> >& m, vector<double>& v, vector<double>& mv);
void GrowSig (vector<double>& q, vector<double>& sigma);
void OptimiseLD (run_params p, double C, ld_info ld, int t, vector<double>& fact_store, ofstream& ld_file, gsl_rng *rgen);

