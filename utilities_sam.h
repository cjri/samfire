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

void GetOptions (run_params& p, int argc, const char **argv);
void makequal (run_params& p, vector<char>& qual);
int ScoreSim (string a, string b, int p, int rsize);
int ScoreSim2 (string a, string b, int p, int rsize);
string RevTr (string a);
string RevString (string a);
void PrintInstructions (int a);
bool compareJoined (joined& a, joined& b);

