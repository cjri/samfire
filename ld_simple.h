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

void ConstructNSites (vector <vector<char> >& m_vars, vector< vector<int> >& n_sites);
int CountSNPsPerRead (vector< vector<int> >& n_sites);
void ConstructLDData (vector< vector<int> >& n_sites, vector<ld_info>& ld_data);
void MakePairgrid (vector<int> polys, vector< vector<int> >& n_sites, vector< vector<int> >& pairgrid);
void ConstructLDPairsData (vector< vector<int> >& pairgrid, vector<ld_info>& ld_data);
void CollateLDData (int t, vector<str>& sltrajs, vector <vector<char> >& m_vars, vector< vector<int> >& n_sites, vector<ld_info>& ld_data);
void ConvertLocusNumbers (vector<str>& sltrajs, vector<ld_info>& ld_data);
void OptimiseLDSetup (vector<double>& fact_store, vector<ld_info>& ld_data);
void RunLDOptimisation (run_params p, int t, vector<double>& fact_store, vector<ld_info>& ld_data, gsl_rng *rgen);


