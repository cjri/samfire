#include "shared.h"
#include "optimisation.h"
#include <iostream>
#include <string>
#include <sstream>

void ConstructNSites (vector <vector<char> >& m_vars, vector< vector<int> >& n_sites) {
    for (int i=0;i<m_vars.size();i++) {//Count polymorphisms covered in each read
        vector<int> sites;
        for (int j=0;j<m_vars[i].size();j++){
            if (m_vars[i][j]!='-') {
                sites.push_back(j);
            }
        }
        n_sites.push_back(sites);
    }
}

int CountSNPsPerRead (vector< vector<int> >& n_sites) {
    int maxs=0;
    for (int i=0;i<n_sites.size();i++) {
       // n_polys.push_back(n_sites[i].size());
        if (n_sites[i].size()>maxs) {
            maxs=n_sites[i].size();
        }
    }
    return maxs;
}

void ConstructLDData (vector< vector<int> >& n_sites, vector<ld_info>& ld_data) {
    //Create structure to contain LD information.
    //Here we work out which pairs of sites are to be represented in the structure
    
    //Question around here : What is in n_sites?  Better way to find unique entries and list in ld_info?  At the moment this is taking a long time for long lists...
    
    
    for (int i=0;i<n_sites.size();i++) {
        if (i%10000==0) {
            cout << i << " " << n_sites.size() << " " << ld_data.size() << "\n";
        }
        if (n_sites[i].size()>1) {
            for (int j=0;j<n_sites[i].size();j++) {
                for (int k=j+1;k<n_sites[i].size();k++) {
                    ld_info l;
                    l.i=n_sites[i][j];
                    l.j=n_sites[i][k];
                    //Check uniqueness
                    int uniq=1;
                    for (int ll=0;ll<ld_data.size();ll++) {
                        if (ld_data[ll].i==n_sites[i][j]&&ld_data[ll].j==n_sites[i][k]) {
                            uniq=0;
                            break;
                        }
                    }
                    if (uniq==1) {
                        ld_data.push_back(l);
                    }
                }
            }
        }
    }
    cout << "LD data size " << ld_data.size() << "\n";
    for (int i=0;i<ld_data.size();i++) {
        ld_data[i].n_i1.push_back(0);
        ld_data[i].n_i0.push_back(0);
        ld_data[i].n_j1.push_back(0);
        ld_data[i].n_j0.push_back(0);
        ld_data[i].n_11.push_back(0);
        ld_data[i].n_10.push_back(0);
        ld_data[i].n_01.push_back(0);
        ld_data[i].n_00.push_back(0);
    }
}

void MakePairgrid (vector<int> polys, vector< vector<int> >& n_sites, vector< vector<int> >& pairgrid) {
    //Find pairs of loci for which there is potentially LD information
    vector<int> pp(polys.size(),0);
    for (int i=0;i<polys.size();i++) {
        pairgrid.push_back(pp);
    }
    for (int i=0;i<n_sites.size();i++) {
        if (n_sites[i].size()>1) {
            for (int j=0;j<n_sites[i].size();j++) {
                for (int k=j+1;k<n_sites[i].size();k++) {
                    pairgrid[n_sites[i][j]][n_sites[i][k]]=1;
                }
            }
        }
    }
}

void ConstructLDPairsData (vector< vector<int> >& pairgrid, vector<ld_info>& ld_data) {
    //Construct dataset of LD information
    for (int i=0;i<pairgrid.size();i++) {
        for (int j=0;j<pairgrid[i].size();j++) {
            if (pairgrid[i][j]==1) {
                ld_info l;
                l.i=i;
                l.j=j;
                ld_data.push_back(l);
            }
        }
    }
    cout << "LD data size " << ld_data.size() << "\n";
    for (int i=0;i<ld_data.size();i++) {
        ld_data[i].n_i1.push_back(0);
        ld_data[i].n_i0.push_back(0);
        ld_data[i].n_j1.push_back(0);
        ld_data[i].n_j0.push_back(0);
        ld_data[i].n_11.push_back(0);
        ld_data[i].n_10.push_back(0);
        ld_data[i].n_01.push_back(0);
        ld_data[i].n_00.push_back(0);
    }
}


void CollateLDData (int t, vector<str>& sltrajs, vector <vector<char> >& m_vars, vector< vector<int> >& n_sites, vector<ld_info>& ld_data) {
    //Fill in the structure with counts of numbers of variants
    for (int i=0;i<n_sites.size();i++) {
        for (int l1=0;l1<ld_data.size();l1++) {
            for (int j=0;j<n_sites[i].size();j++) {
                //One-locus stats
                if (n_sites[i][j]==ld_data[l1].i) {
                    //Check data in m_vars
                    if (m_vars[i][n_sites[i][j]]==sltrajs[ld_data[l1].i].cons) {
                        ld_data[l1].n_i0[t]++;
                    } else if (m_vars[i][n_sites[i][j]]==sltrajs[ld_data[l1].i].nuc) {
                        ld_data[l1].n_i1[t]++;
                    }
                }
                if (n_sites[i][j]==ld_data[l1].j) {
                    //Check data in m_vars
                    if (m_vars[i][n_sites[i][j]]==sltrajs[ld_data[l1].j].cons) {
                        ld_data[l1].n_j0[t]++;
                    } else if (m_vars[i][n_sites[i][j]]==sltrajs[ld_data[l1].j].nuc) {
                        ld_data[l1].n_j1[t]++;
                    }
                }
                //Two-locus stats
                for (int k=j+1;k<n_sites[i].size();k++) {
                    if (n_sites[i][j]==ld_data[l1].i&&n_sites[i][k]==ld_data[l1].j) {
                        if (m_vars[i][n_sites[i][j]]==sltrajs[ld_data[l1].i].cons&&m_vars[i][n_sites[i][k]]==sltrajs[ld_data[l1].j].cons) {
                            ld_data[l1].n_00[t]++;
                            ld_data[l1].n_i0[t]--;
                            ld_data[l1].n_j0[t]--;
                        }
                        if (m_vars[i][n_sites[i][j]]==sltrajs[ld_data[l1].i].cons&&m_vars[i][n_sites[i][k]]==sltrajs[ld_data[l1].j].nuc) {
                            ld_data[l1].n_01[t]++;
                            ld_data[l1].n_i0[t]--;
                            ld_data[l1].n_j1[t]--;
                        }
                        if (m_vars[i][n_sites[i][j]]==sltrajs[ld_data[l1].i].nuc&&m_vars[i][n_sites[i][k]]==sltrajs[ld_data[l1].j].cons) {
                            ld_data[l1].n_10[t]++;
                            ld_data[l1].n_i1[t]--;
                            ld_data[l1].n_j0[t]--;
                        }
                        if (m_vars[i][n_sites[i][j]]==sltrajs[ld_data[l1].i].nuc&&m_vars[i][n_sites[i][k]]==sltrajs[ld_data[l1].j].nuc) {
                            ld_data[l1].n_11[t]++;
                            ld_data[l1].n_i1[t]--;
                            ld_data[l1].n_j1[t]--;
                        }
                    }
                }
            }
        }
    }
}

void ConvertLocusNumbers (vector<str>& sltrajs, vector<ld_info>& ld_data) {
    //Change variant numbers from 0, 1, 2, 3 etc to positions in genome
    for (int i=0;i<ld_data.size();i++) {
        ld_data[i].i=sltrajs[ld_data[i].i].locus;
        ld_data[i].j=sltrajs[ld_data[i].j].locus;
    }
}

void OptimiseLDSetup (vector<double>& fact_store, vector<ld_info>& ld_data) {
    int n=0;
    for(int i=0;i<ld_data.size();i++) {
        int c=ld_data[i].n_i1[0]+ld_data[i].n_i0[0]+ld_data[i].n_j1[0]+ld_data[i].n_j0[0]+ld_data[i].n_11[0]+ld_data[i].n_10[0]+ld_data[i].n_01[0]+ld_data[i].n_00[0];
        if (c>n) {
            n=c;
            
        }
    }
    FindLogFact(fact_store,n);
}

void RunLDOptimisation (run_params p, int t, vector<double>& fact_store, vector<ld_info>& ld_data, gsl_rng *rgen) {
    ofstream ld_file;
    ostringstream convert;
    convert << t;
    string temp=convert.str();
    string name = "LD_info"+temp+".out";
    ld_file.open(name.c_str());
    for (int i=0;i<ld_data.size();i++) {
        int ndouble=ld_data[i].n_11[0]+ld_data[i].n_10[0]+ld_data[i].n_01[0]+ld_data[i].n_00[0];
        if (ndouble>=100) {
            OptimiseLD(p,200,ld_data[i],0,fact_store,ld_file,rgen);
        }
    }

}
