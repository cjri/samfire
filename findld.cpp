#include "shared.h"
#include "optimisation.h"
#include <iostream>
#include <string>

void CollectLDInfo (vector<str> sltrajs, vector<ld_info>& ld_data) {
	for (int i=0;i<sltrajs.size();i++) {
		for (int j=0;j<sltrajs.size();j++) {
			if (i!=j) {
				ld_info ld;
				ld.i=sltrajs[i].locus;
				ld.j=sltrajs[j].locus;
				ld_data.push_back(ld);
			}
		}
	}
}

void ProcessLDInfo (vector<double>& fact_store, vector<ld_info>& ld_data) {
	vector<ld_info> ld_data2;
	int max_tot=0;
	for (int i=0;i<ld_data.size();i++) {
		int tot=0;
		for (int j=0;j<ld_data[i].n_11.size();j++) {
			tot=tot+ld_data[i].n_11[j];
		}
		for (int j=0;j<ld_data[i].n_10.size();j++) {
			tot=tot+ld_data[i].n_10[j];
		}
		for (int j=0;j<ld_data[i].n_01.size();j++) {
			tot=tot+ld_data[i].n_01[j];
		}
		for (int j=0;j<ld_data[i].n_00.size();j++) {
			tot=tot+ld_data[i].n_00[j];
		}
		if (tot>100&&ld_data[i].i!=ld_data[i].j) {
			ld_data2.push_back(ld_data[i]);
			//cout << ld_data[i].i << " " << ld_data[i].j << "\n";
		}
		for (int j=0;j<ld_data[i].n_i1.size();j++) {
			tot=tot+ld_data[i].n_i1[j]+ld_data[i].n_i0[j]+ld_data[i].n_j1[j]+ld_data[i].n_j0[j];
		}
		if (tot>max_tot) {
			max_tot=tot;
		}
	}
	ld_data=ld_data2;
	FindLogFact(fact_store,max_tot);
}







