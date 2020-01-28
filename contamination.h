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
void FindRefLen (int& len, string ref1, string ref2);
void FindPolymorphisms (run_params p, int len, string ref1, string ref2, int& check, int& start_c, vector<int>& polys);
void FindDifferences (run_params p, int len, int start_c, string ref1, string ref2, string c1, vector<int>& polys, vector<joined>& t_read, vector<int>& filter);
int CountPolys (int i, vector<int>& polys, vector<joined>& t_read);

