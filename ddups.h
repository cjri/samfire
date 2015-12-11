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

//Delete duplicate reads
void DelDupSequences (run_params p, int& found_pairs, vector<char> qual, alldat& a);
void DelDups (run_params p, vector<char> qual, vector<rd>& data);
void ProcessMatch (int i, int j, int pair, int& n_del, run_params p, vector<char> qual, vector<rd>& data);
int GetQual(run_params p, int pair, vector<char> qual, int i, vector<rd>& data);
void SQual (run_params p, string q, string s, vector<char> qual, vector<int>& qvec);

