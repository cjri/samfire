#include "shared.h"
#include "data.h"
#include "bioinf.h"
#include <iostream>
#include <sstream>
#include <string>

char TranslateCodon (vector<char> c, vector< vector<char> > cod, vector<char> aa) {
	//Routine to translate a codon (slow but effective)
	char s='X';
	for (int k=0;k<64;k++) {
		if (c[0]==cod[k][0]&&c[1]==cod[k][1]&&c[2]==cod[k][2]) {
			s=aa[k];
			break;
		}
	}
	return s;
}

char FastTranslateCodon (vector<char> c, vector<char> aa) {
	char s='X';
	int k=0;
	if (c[0]=='C') {
		k=k+16;
	}
	if (c[0]=='G') {
		k=k+32;
	}
	if (c[0]=='T') {
		k=k+48;
	}
	if (c[1]=='C') {
		k=k+4;
	}
	if (c[1]=='G') {
		k=k+8;
	}
	if (c[1]=='T') {
		k=k+12;
	}
	if (c[2]=='C') {
		k=k+1;
	}
	if (c[2]=='G') {
		k=k+2;
	}
	if (c[2]=='T') {
		k=k+3;
	}
	s=aa[k];
	return s;
}

int FastTranslateCodonN (char c1, char c2, char c3) {
	int k=0;
	int chk=1;
	if (c1=='-') {
		chk=0;
	} else {
		if (c1=='C') {
			k=k+16;
		}
		if (c1=='G') {
			k=k+32;
		}
		if (c1=='T') {
			k=k+48;
		}
	}
	if (c2=='-') {
		chk=0;
	} else {
		if (c2=='C') {
			k=k+4;
		}
		if (c2=='G') {
			k=k+8;
		}
		if (c2=='T') {
			k=k+12;
		}
	}
	if (c3=='-') {
		chk=0;
	} else {
		if (c3=='C') {
			k=k+1;
		}
		if (c3=='G') {
			k=k+2;
		}
		if (c3=='T') {
			k=k+3;
		}
	}
	if (chk==0) {
		k=-1;
	}
	return k;
}

char TranslateNAA (int n,vector<char> aa) {
	vector<char> c (3,' ');
	if (n>48) {
		c[0]='T';
		n=n-48;
	} else if (n>32) {
		c[0]='G';
		n=n-32;
	} else if (n>16) {
		c[0]='C';
		n=n-16;
	} else {
		c[0]='A';
	}
	if (n>12) {
		c[1]='T';
		n=n-12;
	} else if (n>8) {
		c[1]='G';
		n=n-8;
	} else if (n>4) {
		c[1]='C';
		n=n-4;
	} else {
		c[1]='A';
	}
	if (n==3) {
		c[2]='T';
	} else if (n==2) {
		c[2]='G';
	} else if (n==1) {
		c[2]='C';
	} else {
		c[2]='A';
	}
	char k=FastTranslateCodon (c,aa);
	return k;
}

void RepairConsensusSequence (run_params p, int i, vector< vector<char> >& all_cons, vector<char>& global_cons) {
	int repaired=1;
	if (p.repair_consensus==1) {
		if (i==0) {
			cout << "Repairing gaps in the nucleotide consensus for translation using global consensus sequence\n";
		}
		//Repair ambiguous nucleotides in the aaseq consensus using the global consensus sequence
		for (int j=p.trans_start-1;j<all_cons[i].size();j++) {
			repaired=1;
			if (all_cons[i][j]=='N') {
				repaired=0;
				//cout << "Position " << j << " is " << all_cons[i][j] << " with global " << global_cons[j] << "\n";
				if (global_cons[j]!='N') {
					all_cons[i][j]=global_cons[j];
					repaired=1;
				}
			}
			if (repaired==0) {
				cout << "Error in repair process\n";
				cout << "Nucleotide consensus sequence " << i << " contains an ambiguous nucleotide at position " << j << "\n";
			}
		}
	}
	if (p.repair_consensus==2) {
		if (i==0) {
			cout << "Repairing gaps in the nucleotide consensus for translation using previous consensus sequence\n";
		}
		//Repair ambiguous nucleotides in the aaseq consensus using the global consensus sequence
		for (int j=p.trans_start-1;j<all_cons[i].size();j++) {
			if (all_cons[i][j]=='N') {
				repaired=0;
				for (int k=0;k<i;k++) {
					if (all_cons[k][j]!='N') {
						all_cons[i][j]=all_cons[k][j];
						repaired=1;
						break;
					}
					if (global_cons[j]!='N') {
						all_cons[i][j]=global_cons[j];
						repaired=1;
					}
				}
				if (repaired==0) {
					cout << "Error in repair process\n";
					cout << "Nucleotide consensus sequence " << i << " contains an ambiguous nucleotide\n";
				}
			}
		}
	}
}

void GetVarTypes (vector<char> c, char s, vector< vector<char> >& cod, vector<char>& aa, vector< vector<int> >& s_types) {
	vector <char> c_orig=c;
	vector <char> nucs;
	nucs.push_back('A');
	nucs.push_back('C');
	nucs.push_back('G');
	nucs.push_back('T');
	for (int pos=0;pos<3;pos++) { //Codon position
		vector<int> types;
		for (int n=0;n<4;n++) { //Nucleotide substitution
			int type=0;
			c=c_orig;
			if (c[pos]!=nucs[n]) {
				c[pos]=nucs[n];
				char s_new=FastTranslateCodon(c,aa);
				if (s_new==s) {
					type=1; //Synonymous mutation
				} else {
					type=2;  //Non-synonymous mutation
					if (s_new=='X') {
						type=3; //Nonsense mutation
					}
				}
			}
			types.push_back(type);
		}
		s_types.push_back(types);
	}
}

void CalculateVariantFrequencies (vector< vector< vector<int> > >& all_nucs, vector< vector< vector<double> > >& all_freqs) {
	for (int i=0;i<all_nucs.size();i++) {
		vector< vector<double> > all_f;
		for (int j=0;j<all_nucs[i].size();j++) {
			vector<double> f;
			double p;
			for (int k=0;k<4;k++) {
				if (all_nucs[i][j][4]==0) {
					p=-1;
				} else {
					p=(all_nucs[i][j][k]+0.)/(all_nucs[i][j][4]+0.);
				}
				f.push_back(p);
			}
			all_f.push_back(f);
		}
		all_freqs.push_back(all_f);
	}
}

void CalculateVariantDistances (vector< vector< vector<double> > >& all_freqs, vector< vector<double> >& distances, vector< vector<double> >& all_counts) {
	for (int i=0;i<all_freqs.size();i++) {
		vector<double> d;
		vector<double> c;
		for (int j=0;j<all_freqs.size();j++) {
			double dist=0;
			int counts=0;
			for (int k=0;k<all_freqs[i].size();k++) {
				double d=0;
				int add_count=0;
				for (int l=0;l<4;l++) {
					if (all_freqs[i][k][l]>=0&&all_freqs[j][k][l]>=0) {
						d=d+abs(all_freqs[i][k][l]-all_freqs[j][k][l]);
						add_count++;
					}
				}
				counts=counts+add_count;
				dist=dist+(d/2);	//Note here divide by two - count across all nucleotides so get distances twice
			}
			d.push_back(dist);
			c.push_back(counts);
			//cout << i << " " << j << " " << dist << " " << counts << "\n";
		}
		distances.push_back(d);
		all_counts.push_back(c);
	}
}

void CalculateVariantDistancesSNS (vector< vector< vector<double> > >& all_freqs, vector< vector< vector<int> > >& all_types, vector< vector<double> >& distances_s, vector< vector<double> >& distances_n, vector< vector<double> >& all_counts_s, vector< vector<double> >& all_counts_n) {
	for (int i=0;i<all_freqs.size();i++) {
		vector<double> d_n;
		vector<double> d_s;
		vector<double> c_n;
		vector<double> c_s;
		for (int j=0;j<all_freqs.size();j++) {
			double dist_n=0;
			double dist_s=0;
			int counts_n=0;
			int counts_s=0;
			for (int k=0;k<all_freqs[i].size();k++) {
				double d_n=0;
				double d_s=0;
				int add_count_n=0;
				int add_count_s=0;
				for (int l=0;l<4;l++) {
					if (all_freqs[i][k][l]>=0&&all_freqs[j][k][l]>=0) {
						if (all_types[i][k][l]==1) {
							d_s=d_s+abs(all_freqs[i][k][l]-all_freqs[j][k][l]);
							add_count_s++;
						} else if (all_types[i][k][l]>1) {
							d_n=d_n+abs(all_freqs[i][k][l]-all_freqs[j][k][l]);
							add_count_n++;
						}
					}
				}
				counts_n=counts_n+add_count_n;
				counts_s=counts_s+add_count_s;
				dist_n=dist_n+d_n;	//Note here - don't divide by two as we have ignored the consensus differences
				dist_s=dist_s+d_s;
			}
			d_n.push_back(dist_n);
			d_s.push_back(dist_s);
			c_n.push_back(counts_n);
			c_s.push_back(counts_s);
			//cout << i << " " << j << " " << dist_n << " " << counts_n << " " << dist_s << " " << counts_s << "\n";
		}
		distances_n.push_back(d_n);
		distances_s.push_back(d_s);
		all_counts_n.push_back(c_n);
		all_counts_s.push_back(c_s);
	}
}

void Find_Reading_Frames (run_params p) {
	//Get main consensus file
	ifstream cons_file;
	cons_file.open(p.in_file);
	string cons;
	cons_file >> cons;
	cons_file >> cons;
	
	//Set up arrays
	vector< vector<char> > cod;
	vector<char> aa;
	SetupAATrans(cod,aa);

	//Look for reading frames
	int start=0;
	const char *nucs = cons.c_str();
	//Make mask vector - don't consider things in frame already evaluated...
	vector<int> mask;
	for (int i=0;i<cons.length();i++) {
		mask.push_back(0);
	}

	//Scan nucleotide sequence
	ofstream rfram_file;
	rfram_file.open("Reading_Frames.dat");
	ofstream orf_file;
	orf_file.open("ORFs.dat");
	while (start<cons.length()-2) {
		vector<char> nuc;
		if (mask[start]==0) {
			nuc.push_back(nucs[start]);
			nuc.push_back(nucs[start+1]);
			nuc.push_back(nucs[start+2]);
			char a=FastTranslateCodon(nuc,aa);
			if (a=='M') {
				int pos=start;
				vector<char> trans;
				vector<char> rframe;
                int finished=0; //Issue here if a reading frame goes over the end of the gene
				while (pos<cons.length()-2) {
					vector<char> nuc;
					nuc.push_back(nucs[pos]);
					nuc.push_back(nucs[pos+1]);
					nuc.push_back(nucs[pos+2]);
					char a=FastTranslateCodon(nuc,aa);
					trans.push_back(a);
					rframe.push_back(nucs[pos]);
					rframe.push_back(nucs[pos+1]);
					rframe.push_back(nucs[pos+2]);
					if (a=='X') {
                        finished=1;
						if (trans.size()>p.len) {
							rfram_file << start << " " << pos+2 << "\n";
							orf_file << start << " " << pos+2 << "\n";
							for (int i=0;i<rframe.size();i++) {
								rfram_file << rframe[i];
							}
							rfram_file << "\n";
							for (int i=0;i<trans.size()-1;i++) {
								rfram_file << trans[i];
							}
							rfram_file << "\n";
						}
						break;
					} else {
						mask[pos]=1;
					}
					pos=pos+3;
				}
                if (finished==0) {
                    if (trans.size()>p.len) {
                        rfram_file << start << " " << pos+2 << "\n";
                        orf_file << start << " " << pos+2 << "\n";
                        for (int i=0;i<rframe.size();i++) {
                            rfram_file << rframe[i];
                        }
                        rfram_file << "\n";
                        for (int i=0;i<trans.size()-1;i++) {
                            rfram_file << trans[i];
                        }
                        rfram_file << "\n";
                    }
                }
			}
		}
		start++;
	}
}

void Find_Variant_Types_Multi (run_params p, vector<string>& sam_files) { //Code for all consensus sequences
	vector< vector<int> > types;

	//Read ORF windows
	ifstream orf_file;
	orf_file.open("ORFs.dat");
	vector< vector<int> > orf;
	int i;
	int j;
	for (int k=0;k<1000;k++) {
		vector<int> o;
		if (!(orf_file >> i)) break;
		if (!(orf_file >> j)) break;
		o.push_back(i);
		o.push_back(j);
		orf.push_back(o);
	}

	//Read total consensus sequence.  Use as backup
	ifstream ca_file;
	ca_file.open("Consensus_all.fa");
	string ca;
	ca_file >> ca;
	ca_file >> ca;
	const char *con_all = ca.c_str();
	
	//Read Consensus sequence from time zero.
	for (int t=0;t<sam_files.size();t++) {
		types.clear();
		ifstream cons_file;
		ostringstream convert;
		convert << t;
		string temp=convert.str();
		string name = "Consensus"+temp+".fa";
		cons_file.open(name.c_str());
		string cons;
		cons_file >> cons;
		cons_file >> cons;
		const char *con_seq = cons.c_str();
		//cout << "Length " << cons.length() << "\n";
		//cout << con_seq << "\n";

		//Set up types vector
		for (int i=0;i<cons.length();i++) {
			vector<int> tt;
			tt.push_back(i+1);
			for (int j=0;j<4;j++) {
				tt.push_back(3);
			}
			types.push_back(tt);
		}

		//Set up arrays
		vector< vector<char> > cod;
		vector<char> aa;
		SetupAATrans(cod,aa);

		//Add in S and NS changes.  Consensus is 0, s is 1, ns is 2.
		for (int i=0;i<orf.size();i++) {
			for (int j=orf[i][0];j<=orf[i][1];j=j+3) {
				//cout << "j= " << j << "\n";
				vector<char> codon;
				vector<char> nuc;
				codon.push_back(con_seq[j]);
				codon.push_back(con_seq[j+1]);
				codon.push_back(con_seq[j+2]);
				if (codon[0]!='A'&&codon[0]!='C'&&codon[0]!='G'&&codon[0]!='T') {
					codon[0]=con_all[j];
				}
				if (codon[1]!='A'&&codon[1]!='C'&&codon[1]!='G'&&codon[1]!='T') {
					codon[1]=con_all[j+1];
				}
				if (codon[2]!='A'&&codon[2]!='C'&&codon[2]!='G'&&codon[2]!='T') {
					codon[2]=con_all[j+2];
				}

				char orig_aa=FastTranslateCodon(codon,aa);
				//if (j<=756&&j>=754) {
					//cout << codon[0] << codon[1] << codon[2] << " " << orig_aa << "\n";
				//}
				vector<int> replace;
				replace.push_back(j+1);
				process_nuc(0,orig_aa,codon,nuc,aa,replace);
				types[j]=replace;
				replace.clear();
				replace.push_back(j+2);
				process_nuc(1,orig_aa,codon,nuc,aa,replace);
				types[j+1]=replace;
				replace.clear();
				replace.push_back(j+3);
				process_nuc(2,orig_aa,codon,nuc,aa,replace);
				types[j+2]=replace;
			}
		}
		
		ofstream type_file;
		name = "Variant_types"+temp+".dat";
		type_file.open(name.c_str());
		for (int i=0;i<types.size();i++) {
			for (int j=0;j<types[i].size();j++) {
				type_file << types[i][j] << " ";
			}
			type_file << "\n";
		}

		type_file.close();
		cons_file.close();

		if (t==0) {
		
			//Next step - Calculate Single locus variants
			if (p.get_in==0) {
				p.in_file="Single_locus_trajectories.out";
			}

			ifstream traj_file;
			traj_file.open(p.in_file);
			int locus;
			char from;
			char to;
			int n;
			int temp;
			string line;
			for (int l=0;l<100000;l++) {
				if (!(traj_file >> locus)) break;
				if (!(traj_file >> from)) break;
				if (!(traj_file >> to)) break;
				if (!(traj_file >> n)) break;
				for (int i=0;i<n;i++) {
					for (int j=0;j<6;j++) {
						if (!(traj_file >> temp)) break;
					}
				}
				cout << locus << " " << from << " " << to << " ";
				//for (int k=0;k<types[locus-1].size();k++) {
				//cout << types[locus-1][k] << " ";
				//}

				if (to=='A') {
					PrintType(locus-1,1,types);
				} else if (to=='C') {
					PrintType(locus-1,2,types);
				} else if (to=='G') {
					PrintType(locus-1,3,types);
				} else if (to=='T') {
					PrintType(locus-1,4,types);
				}
			}
		}
	}
}

void Find_Variant_Types (run_params p) {
	vector< vector<int> > types;

	//Read ORF windows
	ifstream orf_file;
	orf_file.open("ORFs.dat");
	vector< vector<int> > orf;
	int i;
	int j;
	for (int k=0;k<1000;k++) {
		vector<int> o;
		if (!(orf_file >> i)) break;
		if (!(orf_file >> j)) break;
		o.push_back(i);
		o.push_back(j);
		orf.push_back(o);
	}

	//Read Consensus sequence from time zero.
	
	ifstream cons_file;
	cons_file.open("Consensus0.fa");
	string cons;
	cons_file >> cons;
	cons_file >> cons;
	const char *con_seq = cons.c_str();
	//cout << "Length " << cons.length() << "\n";
	//cout << con_seq << "\n";

	//Set up types vector
	for (int i=0;i<cons.length();i++) {
		vector<int> t;
		t.push_back(i+1);
		for (int j=0;j<4;j++) {
			t.push_back(3);
		}
		types.push_back(t);
	}

	//Set up arrays
	vector< vector<char> > cod;
	vector<char> aa;
	SetupAATrans(cod,aa);

	//Add in S and NS changes.  Consensus is 0, s is 1, ns is 2.
	for (int i=0;i<orf.size();i++) {
		//cout << "i= " << i << "\n";
		for (int j=orf[i][0];j<=orf[i][1];j=j+3) {
			//cout << "j= " << j << "\n";
			vector<char> codon;
			vector<char> nuc;
			codon.push_back(con_seq[j]);
			codon.push_back(con_seq[j+1]);
			codon.push_back(con_seq[j+2]);
			
			
			char orig_aa=FastTranslateCodon(codon,aa);
			//if (j<=756&&j>=754) {
				//cout << codon[0] << codon[1] << codon[2] << " " << orig_aa << "\n";
			//}
			vector<int> replace;
			replace.push_back(j);
			process_nuc(0,orig_aa,codon,nuc,aa,replace);
			types[j-1]=replace;
			replace.clear();
			replace.push_back(j+1);
			process_nuc(1,orig_aa,codon,nuc,aa,replace);
			types[j]=replace;
			replace.clear();
			replace.push_back(j+2);
			process_nuc(2,orig_aa,codon,nuc,aa,replace);
			types[j+1]=replace;
		}
	}
	ofstream type_file;
	type_file.open("Variant_types.dat");
	for (int i=0;i<types.size();i++) {
		for (int j=0;j<types[i].size();j++) {
			type_file << types[i][j] << " ";
		}
		type_file << "\n";
	}

	//Next step - Calculate Single locus variants
	if (p.get_in==0) {
		p.in_file="Single_locus_trajectories.out";
	}

	ifstream traj_file;
	traj_file.open(p.in_file);
	int locus;
	char from;
	char to;
	int n;
	int temp;
	string line;
	for (int l=0;l<100000;l++) {
		if (!(traj_file >> locus)) break;
		if (!(traj_file >> from)) break;
		if (!(traj_file >> to)) break;
		if (!(traj_file >> n)) break;
		for (int i=0;i<n;i++) {
			for (int j=0;j<6;j++) {
				if (!(traj_file >> temp)) break;
			}
		}
		cout << locus << " " << from << " " << to << " ";
		//for (int k=0;k<types[locus-1].size();k++) {
			//cout << types[locus-1][k] << " ";
		//}

		if (to=='A') {
			PrintType(locus-1,1,types);
		} else if (to=='C') {
			PrintType(locus-1,2,types);
		} else if (to=='G') {
			PrintType(locus-1,3,types);
		} else if (to=='T') {
			PrintType(locus-1,4,types);
		}
	}
}


void PrintType(int locus, int i, vector< vector<int> > types) {
	if (types[locus-1][i]==3) {
		cout << "NC\n";
	} else if (types[locus-1][i]==2) {
		cout << "NS\n";
	} if (types[locus-1][i]==1) {
		cout << "S\n";
	}
}
void process_nuc(int pos, char orig_aa, vector<char> codon, vector<char> nuc, vector<char>& aa, vector<int>& replace) {
	//cout << pos << " " << replace.size() << " ";
	if (codon[pos]=='A') {
		replace.push_back(0);
		add_to_replace(pos,orig_aa,'C',codon,nuc,aa,replace);
		add_to_replace(pos,orig_aa,'G',codon,nuc,aa,replace);
		add_to_replace(pos,orig_aa,'T',codon,nuc,aa,replace);
		//cout << "A " << replace.size() << "\n";
	} else if (codon[pos]=='C') {
		add_to_replace(pos,orig_aa,'A',codon,nuc,aa,replace);
		replace.push_back(0);
		add_to_replace(pos,orig_aa,'G',codon,nuc,aa,replace);
		add_to_replace(pos,orig_aa,'T',codon,nuc,aa,replace);
		//cout << "C " << replace.size() << "\n";
	} else if (codon[pos]=='G') {
		add_to_replace(pos,orig_aa,'A',codon,nuc,aa,replace);
		add_to_replace(pos,orig_aa,'C',codon,nuc,aa,replace);
		replace.push_back(0);
		add_to_replace(pos,orig_aa,'T',codon,nuc,aa,replace);
		//cout << "G " << replace.size() << "\n";
	} else if (codon[pos]=='T') {
		add_to_replace(pos,orig_aa,'A',codon,nuc,aa,replace);
		add_to_replace(pos,orig_aa,'C',codon,nuc,aa,replace);
		add_to_replace(pos,orig_aa,'G',codon,nuc,aa,replace);
		replace.push_back(0);
		//cout << "T " << replace.size() << "\n";
	}
	
}

void add_to_replace(int pos, char orig_aa, char v, vector<char> codon, vector<char> nuc, vector<char>& aa, vector<int>& replace) {
	nuc=codon;
	nuc[pos]=v;
	char new_aa=FastTranslateCodon(nuc,aa);
	aa_compare(new_aa,orig_aa,replace);
}

void aa_compare(char new_aa, char orig_aa, vector<int>& replace) {
	if (new_aa==orig_aa) {
		replace.push_back(1);
	} else {
		replace.push_back(2);
	}
}
