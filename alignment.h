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

void AlignSequencesSam (run_params& p, int s_length, vector<char> qual, rseq refseq, vector<rd>& data);
void MinBaseQual (vector<char> qual, string& q0);
void ReadCigar (int i, vector<rd>& data);
void FilterCigar (int i, vector<rd>& data);
int findqual (run_params p, int i, int min_qual, int max_qual, vector<char> qual, vector<rd>& data);
int GetMedian (int a, int b, vector<int> qvec);
void RemoveInitialSoftClipping (int i, vector<rd>& data);
void FixDeletions (int i, string q0, vector<rd>& data);
void FixInsertions (int i, vector<rd>& data);
void RemoveSoftClipping (int i, vector<rd>& data);
void ProcessReadQual (int i, run_params p, vector<char> qual, vector<rd>& data);
void AlignSequencesSlide (int s_length, run_params& p, vector<char> qual, rseq refseq, vector<rd>& data);
void DoForwardSlideAlignment (run_params p, int i, int& firstpos, int&isok, vector<char> qual, rseq refseq, vector<rd>& data);
void DoReverseSlideAlignment (run_params p, int i, int& firstpos, int&isok, vector<char> qual, rseq refseq, vector<rd>& data);
void SlideAlign (string read, string ref, double id, int& isok, int i, vector<rd>& data, int dir);

