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

void MakeInitialHaps (run_params p, vector< vector<mpoly> > c_m_polys, vector< vector<char> >& haps);
void PrintHaps (vector< vector<char> > haps);
void ConstructFullHaps2 (run_params p, vector< vector<char> >& haps);
void OverlapStepShared1 (run_params p, vector< vector<char> >& haps);
void OverlapStepShared2 (run_params p, vector< vector<char> >& haps);
void OverlapStepNoShare (run_params p, vector< vector<char> >& haps);
void GetStartFinish(vector<par>& sf, vector< vector<char> >& haps);
void BuildPaths(vector<par>& sf, vector< vector<char> > haps, vector< vector<char> >& new_haps);
void AddPath (int prev, int start, vector<char> v, vector<par>& sf, vector< vector<char> >& haps, vector< vector<char> >& new_haps);
void DuplicateStep (run_params p, vector< vector<char> >& haps);
void ResetHaps (vector<int> incl, vector< vector<char> >& haps);




