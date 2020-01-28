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
void JoinPairs (int i, run_params p, vector<joined>& jseqs, vector<char>& qual, vector<rd> data);
void JoinWithPairs (run_params p, vector<joined>& jseqs, vector<char>& qual, vector<rd> data, ofstream& jp_file, ofstream& qp_file);
void JoinNoPairs (run_params p, vector<joined>& jseqs, vector<char>& qual, vector<rd> data, ofstream& jp_file, ofstream& qp_file);
void TranslateQualityVector(string s, vector<char>& qual, ofstream& qp_file);

