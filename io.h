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

void GetRefSeq (ifstream& ref_file, rseq& ref);
void ImportSamFileNames (run_params p, vector<string>& sam_files);
void ImportSamFile (run_params p, int i, rseq refseq, vector<string> sam_files, alldat& a);
void ReadSamFile (run_params p, ifstream& sam_file, int& s_length, vector<rd>& data, rseq refseq);
void BinDecomp (int i, vector<int>& v);
void OutputAlign (int i, vector<rd> data);
void OutAl (vector<rd> data, ofstream& al_file);
void InputVariantData (vector<string> sam_files, vector< vector< vector<int> > >& all_nucs);
void OutputParticularConsensus (string name2, string name3, vector<char> cons);
void OutputGlobalConsensus (vector< vector< vector<int> > >& all_nucs);
void GetVariantTotals (vector<string> sam_files, vector< vector<int> >& all_tots);
void InputJoinedData (int n_times, vector< vector<joined> >& t_reads);
void InputJnData (int i, vector<joined>& t_read);
void OutputVarFile (int i, rseq refseq, nuc r_count);
void OutputSLTData (run_params p, const char* filename, vector<str> sltrajs);
void ImportSLTData (run_params p, vector<str>& sltrajs);
void ImportSLTLoci (const char* filename, int& check, vector<int>& loci);
void ImportTimeData (vector<int>& times);
void InputHaplotypeData (run_params p, vector<int> polys, vector< vector<char> >& haps);
void GenerateOutputFiles (run_params p, vector<int> polys, vector< vector<char> > haps, vector< vector< vector<mpoly> > > m_polys, vector< vector<mpoly> > c_m_polys, vector< vector<mtr> >& mltrajs);
void GenerateOutputFilesNoHaps (run_params p, vector<int> polys, vector< vector< vector<mpoly> > >& m_polys, vector< vector<mpoly> > c_m_polys);
void ReadMultiMLTrajs (run_params p, vector< vector<int> >& times, vector< vector<mtr> >& mltrajs);
void ReadMultiLocTrajs (run_params p, vector< vector<int> >& times, vector< vector<mtr> >& mltrajs);
void ReadMultiLocTrajs2 (run_params p, vector<str> sltrajs, vector<ld_info>& ld_data);

