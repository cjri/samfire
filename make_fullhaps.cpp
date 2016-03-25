#include "shared.h"
#include "make_fullhaps.h"
#include <iostream>
#include <string>

void MakeInitialHaps (run_params p, vector< vector<mpoly> > c_m_polys, vector< vector<char> >& haps) {
	for (int i=0;i<c_m_polys.size();i++) {
		for (int j=0;j<c_m_polys[i].size();j++) {
			haps.push_back(c_m_polys[i][j].vars);
		}
	}
	if (p.verb==1) {
		PrintHaps(haps);
	}
}

void PrintHaps (vector< vector<char> > haps) {
	cout << "Number of haps " << haps.size() << "\n";
	for (int i=0;i<haps.size();i++) {
		for (int j=0;j<haps[i].size();j++) {
			cout << haps[i][j];
		}
		cout << "\n";
	}
}

void ConstructFullHaps2 (run_params p, vector< vector<char> >& haps) {
	//Cycle through overlap steps
	vector<vector<char> > haps2=haps;
	
	//Matching overlap joins
	for (int i=0;i<20;i++) {
		haps2=haps;
		OverlapStepShared1(p,haps);
		OverlapStepShared2(p,haps);
		if (haps==haps2) {
			break;
		}
	}
	
	//Non-matching overlap joins
	for (int i=0;i<20;i++) {
		haps2=haps;
		OverlapStepNoShare(p,haps);
		OverlapStepShared1(p,haps);
		OverlapStepShared2(p,haps);
		if (haps==haps2) {
			break;
		}
	}

	//Non-overlap joins
	vector<par> sf;
	GetStartFinish(sf,haps);
	vector< vector<char> > new_haps;
	BuildPaths(sf,haps,new_haps);
	haps=new_haps;
	DuplicateStep(p,haps);
	
	if (p.verb==1) {
		PrintHaps(haps);
	}

	ofstream haps_file;
	haps_file.open("Haps1.dat");
	for (int i=0;i<haps.size();i++) {
		for (int j=0;j<haps[i].size();j++) {
			haps_file << haps[i][j];
		}
		haps_file << "\n";
	}
	haps_file.close();

}



void OverlapStepShared1 (run_params p, vector< vector<char> >& haps) {
	//One contained within another - try doing these first before other overlaps...
	if (p.verb==1) {
		cout << "Overlap step\n";
	}
	vector<int> incl (haps.size(),1);
	int size=haps.size();
	for (int i=0;i<size;i++) {
		for (int j=0;j<size;j++) {
			if (i!=j) {
				int match=1;   //Match at shared positions
				int n_match=0;  //Number of shared positions
				for	(int k=0;k<haps[i].size();k++) {
					//Check whether j is contained in i.
					if (haps[j][k]!='-') {
						if (haps[i][k]!='-') {
							n_match++;
						}
						if (haps[j][k]!=haps[i][k]) {
							match=0;
						}
					}
				}
				if (match==1&&n_match>0) {
					if (p.verb==1) {
						for	(int k=0;k<haps[i].size();k++) {
							cout << haps[i][k];
						}
						cout << " contains ";
						for	(int k=0;k<haps[j].size();k++) {
							cout << haps[j][k];
						}
						cout << "\n";
					}
					//Remove haps[j]
					incl[j]=0;
				}
			}
		}
	}
	ResetHaps(incl,haps);
	DuplicateStep(p,haps);
	if (p.verb==1) {
		PrintHaps(haps);
	}
}


void OverlapStepShared2 (run_params p, vector< vector<char> >& haps) {
	if (p.verb==1) {
		cout << "Overlap step\n";
	}
	vector<int> incl (haps.size(),1);
	int size=haps.size();
	for (int i=0;i<size;i++) {
		for (int j=i+1;j<size;j++) {
			int match=1;   //Match at shared positions
			int n_match=0;  //Number of shared positions
			for	(int k=0;k<haps[i].size();k++) {
				//Check whether i and j have the same nucleotides wherever they overlap.
				if (haps[j][k]!='-'&&haps[i][k]!='-') { //Number of overlaps
					n_match++;
					if (haps[j][k]!=haps[i][k]) {
						match=0;
					}
				}
			}
			if (match==1&&n_match>0) {
				for	(int k=0;k<haps[i].size();k++) {
					cout << haps[i][k];
				}
				cout << " overlaps ";
				for	(int k=0;k<haps[j].size();k++) {
					cout << haps[j][k];
				}
				cout << "\n";
				vector<char> newhap;
				for	(int k=0;k<haps[i].size();k++) {
					if (haps[j][k]!='-') {
						newhap.push_back(haps[j][k]);
					} else if (haps[i][k]!='-') {
						newhap.push_back(haps[i][k]);
					} else {
						newhap.push_back('-');
					}
				}
				haps.push_back(newhap);
				incl.push_back(1);
				incl[i]=0;
				incl[j]=0;
			}
		}
	}
	ResetHaps(incl,haps);
	DuplicateStep(p,haps);
	if (p.verb==1) {
		PrintHaps(haps);
	}
}


void OverlapStepNoShare (run_params p, vector< vector<char> >& haps) {
	if (p.verb==1) {
		cout << "Overlap step 2\n";
	}
	vector<int> incl (haps.size(),1);
	int size=haps.size();
	for (int i=0;i<size;i++) {
		for (int j=i+1;j<size;j++) {
			int match=1;   //Match at shared positions
			int n_match=0;  //Number of shared positions
			for	(int k=0;k<haps[i].size();k++) {
				//Check whether i and j have the same nucleotides wherever they overlap.  N.B. Must have at least one locus not in common
				if (haps[j][k]!='-'&&haps[i][k]!='-') { //Number of overlaps
					n_match++;
					if (haps[j][k]!=haps[i][k]) {
						match=0;
					}
				}
			}
			int n_over=0;
			if (match==0&&n_match<haps[i].size()&&n_match>0) {
				for	(int k=0;k<haps[i].size();k++) {
					if (haps[j][k]!='-'&&haps[i][k]=='-') {
						n_over=1;
					}
					if (haps[i][k]!='-'&&haps[j][k]=='-') {
						n_over=1;
					}
				}
			}
			
			if (n_over==1) {
				if (p.verb==1) {
					for	(int k=0;k<haps[i].size();k++) {
						cout << haps[i][k];
					}
					cout << " misoverlaps ";
					for	(int k=0;k<haps[j].size();k++) {
						cout << haps[j][k];
					}
					cout << "\n";
				}
				
				vector<char> nh_start;
				//Make beginning
				int ov=0;
				for	(int k=0;k<haps[i].size();k++) {
					if (haps[i][k]!='-'&&haps[j][k]!='-') {
						ov=1;
					}
					if (ov==0) {
						if (haps[j][k]!='-') {
							nh_start.push_back(haps[j][k]);
						} else if (haps[i][k]!='-') {
							nh_start.push_back(haps[i][k]);
						} else {
							nh_start.push_back('-');
						}
					} else {
						nh_start.push_back('-');
					}
				}
				
				//Make end
				vector<char> nh_end;
				ov=0;
				for	(int k=0;k<haps[i].size();k++) {
					if (haps[i][k]!='-'&&haps[j][k]!='-') {
						ov=1;
					}
					if (ov==1) {
						if (haps[j][k]!='-'&&haps[i][k]=='-') {
							nh_end.push_back(haps[j][k]);
						} else if (haps[i][k]!='-'&&haps[j][k]=='-') {
							nh_end.push_back(haps[i][k]);
						} else {
							nh_end.push_back('-');
						}
					} else {
						nh_end.push_back('-');
					}
				}

				//Make c1 and c2
				vector<char> nh_c1;
				vector<char> nh_c2;
				for	(int k=0;k<haps[i].size();k++) {
					if (nh_start[k]=='-'&&nh_end[k]=='-') {
						nh_c1.push_back(haps[i][k]);
						nh_c2.push_back(haps[j][k]);
					} else {
						nh_c1.push_back('-');
						nh_c2.push_back('-');
					}
				}

				//Construct start + c1 + end
				vector<char> newhap1;
				vector<char> newhap2;
				for (int k=0;k<haps[i].size();k++){
					if (nh_start[k]!='-') {
						newhap1.push_back(nh_start[k]);
						newhap2.push_back(nh_start[k]);
					} else if (nh_c1[k]!='-') {
						newhap1.push_back(nh_c1[k]);
						newhap2.push_back(nh_c2[k]);
					} else if (nh_end[k]!='-') {
						newhap1.push_back(nh_end[k]);
						newhap2.push_back(nh_end[k]);
					} else {
						newhap1.push_back('-');
						newhap2.push_back('-');
					}
				}
				
				haps.push_back(newhap1);
				haps.push_back(newhap2);
				incl.push_back(1);
				incl.push_back(1);
				cout << newhap1.size() << " " << newhap2.size() << "\n";
				incl[i]=0;
				incl[j]=0;
			}
		}
	}
	ResetHaps(incl,haps);
	DuplicateStep(p,haps);
	if (p.verb==1) {
		PrintHaps(haps);
	}
}

void GetStartFinish(vector<par>& sf, vector< vector<char> >& haps) {
	for (int i=0;i<haps.size();i++) {
		par p;
		p.i1=-1;
		p.i2=-1;
		for (int j=0;j<haps[i].size();j++) {
			//cout << "i " << i << " j " << j << " " << haps[i][j] << "\n";
			if (haps[i][j]!='-'&&p.i1<0) {
				p.i1=j;
				//	cout << "Found start\n";
			}
			if (haps[i][j]=='-'&&p.i1>=0&&p.i2<0) {
				p.i2=j-1;
				//	cout << "Found end\n";
			}
		}
		if (p.i2<0) {
			p.i2=haps[i].size();
		}
		sf.push_back(p);
	}
}

void BuildPaths(vector<par>& sf, vector< vector<char> > haps, vector< vector<char> >& new_haps) {
	vector<char> v;
	for (int i=0;i<sf.size();i++) {
		if (sf[i].i1==0) {
			//cout << "Detected " << i << "\n";
			AddPath(i,sf[i].i2+1,v,sf,haps,new_haps);
		}
	}
}

void AddPath (int prev, int start, vector<char> v, vector<par>& sf, vector< vector<char> >& haps, vector< vector<char> >& new_haps) {
	//Add last haplotype to vector
	vector<char> v2=v;
	//cout << "Previous is " << prev << " " << sf[prev].i1 << " " << sf[prev].i2 << "\n";
	for (int j=sf[prev].i1;j<=sf[prev].i2;j++) {
		v2.push_back(haps[prev][j]);
	}
	//cout << "Vector v is: ";
	//for (int i=0;i<v2.size();i++) {
	//	cout << v2[i];
	//}
	//cout << "\n";
	
	//Check if this is the end of the haplotype and if so push to new_haps
	if (sf[prev].i2>=haps[prev].size()) {
		//	cout << "Doing push_back\n";
		vector<char> v3;
		for (int i=0;i<haps[prev].size();i++) {//Trim length
			v3.push_back(v2[i]);
		}
		new_haps.push_back(v3);
	} else {
		//If not search for new haplotype and recall
		for (int i=0;i<sf.size();i++) {
			if (sf[i].i1==start) {
				AddPath(i,sf[i].i2+1,v2,sf,haps,new_haps);
			}
		}
	}
}


void DuplicateStep (run_params p, vector< vector<char> >& haps) {
	if (p.verb==1) {
		cout << "Duplicate step\n";
	}
	vector< vector<char> > new_haps;
	vector<int> uniq;
	for (int i=0;i<haps.size();i++) {
		uniq.push_back(1);
	}
	cout << "Haplotypes " << haps.size() << "\n";
	for (int i=0;i<haps.size();i++) {
		for (int j=i;j<haps.size();j++) {
			if (i!=j&&haps[i]==haps[j]) {
				uniq[j]=0;
			}
		}
	}
	for (int i=0;i<haps.size();i++) {
		if (uniq[i]==1) {
			new_haps.push_back(haps[i]);
		}
	}
	cout << "Unique haplotypes " << new_haps.size() << "\n";
	haps=new_haps;
}

void ResetHaps (vector<int> incl, vector< vector<char> >& haps) {
	vector< vector<char> > new_haps;
	for (int i=0;i<haps.size();i++) {
		if (incl[i]==1) {
			new_haps.push_back(haps[i]);
		}
	}
	haps=new_haps;
}







/*

void ConstructFullHaps (run_params p, vector< vector<char> >& haps) {
	cout << "Inferring full haplotypes...\n";
	
	
	//Step zero: Identify loci with only one nucleotide
	UniqueStep(p,haps);

	//Cycle over first two steps
	for (int c=0;c<=20;c++) {
		vector< vector<char> > haps2=haps;
		//Step one: redundancy
		RedundancyStep (p,haps);

		//Step two: overlap
		OverlapStep (p,haps);
		
	//	DuplicateStep(p,haps);

		if (haps2==haps) {
			cout << "Done\n";
			break;
		}
	}
	//Step 3: Combination 1: Make full haplotypes from disjoint sequences
	CombinationOne(p,haps);
	RedundancyStep (p,haps);
	//DuplicateStep(p,haps);
	
	cout << "Combination Two\n";
	CombinationTwo(p,haps);
	
	RedundancyStep(p,haps);
	
	vector<par> sf;
	GetStartFinish(sf,haps);
//	cout << "Start-finish\n";
//	for (int i=0;i<sf.size();i++) {
//		cout << i << " " << sf[i].i1 << " " << sf[i].i2 << "\n";
//	}

	vector< vector<char> > new_haps;
	BuildPaths(sf,haps,new_haps);
	haps=new_haps;
	if (p.verb==1) {
		PrintHaps(haps);
	}

//	DuplicateStep(p,haps);

	//Step 4: Combination 2: Merge non-full haplotypes with full haplotypes
//	CombinationTwo(p,haps);
//	RedundancyStep(p,haps);
	
	ofstream haps_file;
	haps_file.open("Haps1.dat");
	for (int i=0;i<haps.size();i++) {
		for (int j=0;j<haps[i].size();j++) {
			haps_file << haps[i][j];
		}
		haps_file << "\n";
	}
	haps_file.close();
}

void UniqueStep (run_params p, vector< vector<char> >& haps) {
	for (int j=0;j<haps[0].size();j++) {
		int uniq=1;
		char c='-';
		for (int i=0;i<haps.size();i++) {
			if (c=='-'&&haps[i][j]!='-') {
				c=haps[i][j];
			}
			if (c!='-'&&haps[i][j]!='-'&&haps[i][j]!=c) {
				uniq=0;
			}
		}
		cout << "Unique " << j << " " << uniq << "\n";
		if (uniq==1) {
			for (int i=0;i<haps.size();i++) {
				if (haps[i][j]=='-') {
					if (j==0&&haps[i][1]!='-') {
						haps[i][j]=c;
					} else if (j==haps[0].size()&&haps[i][haps[0].size()-1]!='-') {
						haps[i][j]=c;
					} else if (haps[i][j-1]!='-'||haps[i][j+1]!='-') {
						haps[i][j]=c;
					}
				}
			}
		}
	}
	if (p.verb==1) {
		PrintHaps(haps);
	}
}


void RedundancyStep (run_params p, vector< vector<char> >& haps) {
	if (p.verb==1) {
		cout << "Redundancy step\n";
	}
	vector<int> incl (haps.size(),1);
	for (int i=0;i<haps.size();i++) {
		if (incl[i]==1) {
			for (int j=0;j<haps.size();j++) {
				if (incl[j]==1&&i!=j) {
					int match=1;
					for	(int k=0;k<haps[i].size();k++) {
						//Check whether j has same nucleotides as i whenever j has a nucleotide
						if (haps[j][k]!='-'&&haps[j][k]!=haps[i][k]) {
							match=0;
						}
					}
					if (match==1) {
						incl[j]=0;
					}
				}
			}
		}
	}
	ResetHaps(incl,haps);
	if (p.verb==1) {
		PrintHaps(haps);
	}
}


void OverlapStep (run_params p, vector< vector<char> >& haps) {
	if (p.verb==1) {
		cout << "Overlap step\n";
	}
	vector<int> incl (haps.size(),1);
	int size=haps.size();
	for (int i=0;i<size;i++) {
		if (incl[i]==1) {
			for (int j=0;j<size;j++) {
				if (incl[j]==1&&i!=j) {
					int match=1;
					int n_match=0;
					for	(int k=0;k<haps[i].size();k++) {
						//Check whether i and j have the same nucleotides wherever they overlap.  If so, create a new haplotype
						if (haps[j][k]!='-'&&haps[i][k]!='-') { //Number of overlaps
							n_match++;
						}
						
						
						if (haps[j][k]!='-'&&haps[i][k]!='-'&&haps[j][k]!=haps[i][k]) {
							match=0;
						}
					}
					if (match==1&&n_match>0) {
						vector<char> newhap;
						for	(int k=0;k<haps[i].size();k++) {
							if (haps[j][k]!='-') {
								newhap.push_back(haps[j][k]);
							} else if (haps[i][k]!='-') {
								newhap.push_back(haps[i][k]);
							} else {
								newhap.push_back('-');
							}
						}
						haps.push_back(newhap);
						incl.push_back(1);
					}
				}
			}
		}
	}
	ResetHaps(incl,haps);
	if (p.verb==1) {
		PrintHaps(haps);
	}
}

void CombinationOne (run_params p, vector< vector<char> >& haps) {
	//Merge sets of reads that collectively span all loci
	vector<int> incl (haps.size(),1);
	int size=haps.size();
	for (int i=0;i<size;i++) {
		if (incl[i]==1) {
			for (int j=0;j<size;j++) {
				if (incl[j]==1&&i!=j) {
					int match=1;
					int n_i=0;
					int n_j=0;
					for	(int k=0;k<haps[i].size();k++) {
						//Check whether i and j collectively span all loci.  If so, create new haplotype
						if (haps[i][k]!='-') {
							n_i++;
						}
						if (haps[j][k]!='-') {
							n_j++;
						}
						
						if (haps[j][k]=='-'&&haps[i][k]=='-') {
							match=0;
						}
					}
					if (match==1&&n_i<haps[i].size()&&n_j<haps[i].size()) {
						vector<char> newhap;
						for	(int k=0;k<haps[i].size();k++) {
							if (haps[j][k]!='-') {
								newhap.push_back(haps[j][k]);
							} else if (haps[i][k]!='-') {
								newhap.push_back(haps[i][k]);
							} else {
								newhap.push_back('-');
							}
						}
						haps.push_back(newhap);
						incl.push_back(1);
					}
				}
			}
		}
	}
	ResetHaps(incl,haps);
	if (p.verb==1) {
		PrintHaps(haps);
	}
}




void CombinationTwo (run_params p, vector< vector<char> >& haps) {
	vector<int> incl (haps.size(),1);
	int size=haps.size();
	for (int i=0;i<size;i++) {
		if (incl[i]==1) {
			for (int j=0;j<size;j++) {
				if (incl[j]==1&&i!=j) {
					int n_i=0;
					int n_j=0;
					for	(int k=0;k<haps[i].size();k++) {
						//Check whether i and j collectively span all loci.  If so, create new haplotype
						if (haps[i][k]!='-') {
							n_i++;
						}
						if (haps[j][k]!='-') {
							n_j++;
						}
					}
					if (n_i==haps[i].size()&&n_j<haps[i].size()) {
						//Add in haplotype i with alleles replaced by those of j
						vector<char> newhap=haps[i];
						for	(int k=0;k<haps[i].size();k++) {
							if (haps[j][k]!='-') {
								newhap[k]=haps[j][k];
							}
						}
						haps.push_back(newhap);
						incl.push_back(1);
					}
				}
			}
		}
	}
	ResetHaps(incl,haps);
	if (p.verb==1) {
		PrintHaps(haps);
	}
}
*/