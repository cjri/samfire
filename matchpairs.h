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

//Match paired end reads
void MatchPairs (run_params p, int ireq, vector<rd>& data);
void CheckOldPairs (vector<rd>& data);
void JoinPairs (int i, run_params p, vector<joined>& jseqs, vector<rd> data);
void JoinWithPairs (vector<joined>& jseqs, vector<rd> data, ofstream& jp_file);
void JoinNoPairs (vector<joined>& jseqs, vector<rd> data, ofstream& jp_file);

