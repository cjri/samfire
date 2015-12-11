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

//Call MNPs
void MaxReadLength (run_params p, int& max_r_size, vector< vector<joined> >& t_reads);
void MaxRdLength (run_params p, int& max_r_size,  vector<joined>& t_read);
void MaxReadSpan(run_params p, int max_r_size, int& max_l, vector<int> polys);
void MapSequenceToSNPLoci (run_params p, int n_times, vector<int> polys, vector< vector<joined> > t_reads, vector <vector <vector<char> > >& t_m_vars);

void MapSeqToSNPLoci (run_params p, int n_times, vector<int> polys, vector<joined> t_read, vector <vector<char> >& m_vars);
void CallReadTypesConsecutive (run_params p, int max_l, int max_r_size, vector<int> polys, vector< vector<int> >& l_combs);
void CallMNPs2 (run_params p, vector< vector<int> >& l_combs, vector <vector<char> > m_vars, vector< vector<mpoly> >& m_pol);
void FindMaxLoc(int& max_loc, vector<vector< vector<int> > > all_l_combs);
void ConstructLocVec (run_params p, int max_loc, int max_l_store, vector<vector< vector<int> > > all_l_combs, vector< vector<int> >& l_vec);
void ConstructMPoly (vector< vector<int> >& l_vec, vector<vector< vector<int> > > all_l_combs, vector< vector< vector<mpoly> > >& m_polys);
void MakeAssignedList (vector <vector <vector<char> > > t_m_vars, vector< vector<int> >& done);
void FilterMNPs (run_params p, vector< vector< vector<mpoly> > > m_polys, vector< vector< vector<mpoly> > >& m_polys_f);
void CombineMNPsTime (run_params p, vector< vector< vector<mpoly> > > m_polys, vector< vector<mpoly> >& c_m_polys);
void MLPMeanFreqs (run_params p, vector< vector<mtr> >& mltrajs);
void GetIncludeML (run_params p, vector<int>& include_ml, vector< vector<mtr> > mltrajs);
void FilterMLFreq (run_params p, vector< vector<int> >& times, vector< vector<mtr> >& mltrajs);
void FilterMLFreq2 (run_params p, vector< vector<int> >& times, vector< vector<mtr> >& mltrajs);
void FilterMLFreq3 (run_params p, vector<int>& include_ml, vector< vector<int> >& times, vector< vector<mtr> >& mltrajs);
void GetDepths (int i, vector<int>& depths, vector< vector<mtr> >& mltrajs);
void CompleteMLCalls (vector< vector<char> > haps, vector< vector<mpoly> >& c_m_polys);
