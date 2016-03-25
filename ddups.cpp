#include "shared.h"
#include "matchpairs.h"
#include "alignment.h"
#include "ddups.h"
#include <iostream>
#include <string>
#include <sstream>


void DelDupSequences (run_params p, int& found_pairs, vector<char> qual, alldat& a) {
	//Find pairs
	if (p.pairs==1) {
		MatchPairs(p,0,a.data);
	}
	found_pairs=1;
	DelDups(p,qual,a.data);
}


void DelDups (run_params p, vector<char> qual, vector<rd>& data) {
	int n_del=0;
	cout << "Removing duplicate sequences...\n";
	for (int i=0;i<data.size();i++) {
//		cout << i << " " << data[i].pairno << "\n";
		if (data[i].del==0) {
			for (int j=i+1;j<data.size();j++) {
				if (data[j].del==0) {
					if (data[j].alpos==data[i].alpos) {
						if (p.ddup==0) {
							//Sequence equality measure
							if (data[j].seq==data[i].seq) {
								//Check paired reads
								if (data[i].pairno>=0&&data[j].pairno>=0) { //Both reads are paired; check for equality across paired reads
									if (data[data[i].pairno].alpos==data[data[j].pairno].alpos&&data[data[i].pairno].seq==data[data[j].pairno].seq) {
										//Compare base quality
										ProcessMatch (i,j,1,n_del,p,qual,data);
									}
								} else if (data[i].pairno==-1&&data[j].pairno==-1){ //Neither read has a pair
									ProcessMatch (i,j,0,n_del,p,qual,data);
								}
							}
						} else {
							//Sequence similarity measure
							int n_diff=0;
							//Check paired reads
							if (data[i].pairno>=0&&data[j].pairno>=0) {
								if (data[data[i].pairno].alpos==data[data[j].pairno].alpos) {
									//Count differences in first reads.  Require reads to be same length; cover same bases
									if (data[i].seq.length()==data[j].seq.length()) {
										for (int k=0;k<data[i].seq.length();k++) {
											if (data[i].seq[k]!=data[j].seq[k]) {
												n_diff++;
											}
										}
									} else {
										n_diff=n_diff+1000;
									}
									
									//Count differences in second reads
									if (data[data[i].pairno].seq.length()==data[data[j].pairno].seq.length()) {
										for (int k=0;k<data[data[i].pairno].seq.length();k++) {
											if (data[data[i].pairno].seq[k]!=data[data[j].pairno].seq[k]) {
												n_diff++;
											}
										}
									} else {
										n_diff=n_diff+1000;
									}
									
									if (n_diff<=p.ddup) {
										ProcessMatch (i,j,1,n_del,p,qual,data);
									}
								}
							} else if (data[i].pairno==-1&&data[j].pairno==-1) { //No paired reads; compare at single read level
								//Count differences in first reads
								if (data[i].seq.length()==data[j].seq.length()) {
									for (int k=0;k<data[i].seq.length();k++) {
										if (data[i].seq[k]!=data[j].seq[k]) {
											n_diff++;
										}
									}
								} else {
									n_diff=n_diff+1000;
								}
								if (n_diff<=p.ddup) {
									ProcessMatch (i,j,0,n_del,p,qual,data);
								}
							}
						}
					}
				}
			}
		}
	}
	cout << "Number of deletions " << n_del << "\n";
}

void ProcessMatch (int i, int j, int pair, int& n_del, run_params p, vector<char> qual, vector<rd>& data) {
	int median_i=GetQual(p,pair,qual,i,data);
	int median_j=GetQual(p,pair,qual,j,data);
	if (median_j>median_i) { //Delete i
		data[i].del=1;
		if (pair==1) {
			data[data[i].pairno].del=1;
		}
	} else {  //Delete j
		data[j].del=1;
		if (pair==1) {
			data[data[j].pairno].del=1;
		}
	}
	n_del++;
}

int GetQual(run_params p, int pair, vector<char> qual, int i, vector<rd>& data) {
	vector<int> qvec;
	string q=data[i].qual;
	string s=data[i].seq;
	SQual(p,q,s,qual,qvec);
	if (pair==1) {
		q=data[data[i].pairno].qual;
		s=data[data[i].pairno].seq;
		SQual(p,q,s,qual,qvec);
	}
	int m = GetMedian(0,0,qvec);
	return m;
}

void SQual (run_params p, string q, string s, vector<char> qual, vector<int>& qvec) {
	for (int k=0;k<s.size();k++) {
		for (int l=0;l<=p.max_qual;l++) {
			if (q[k]==qual[l]) {
				qvec.push_back(l);
				break;
			}
		}
	}
}