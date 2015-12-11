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

//Likelihood
double DirichletMultiCalc (int N, double c, vector<int> obs, vector<double> inf, vector<double>& fact_store);

