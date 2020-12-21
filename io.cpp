#include "shared.h"
#include "io.h"
#include <iostream>
#include <string>
#include <sstream>

void GetRefSeq (ifstream& ref_file, rseq& refseq) {
	//Import the reference sequence, in fasta format
	cout << "Getting reference sequence...\n";
	string rs;
	string s;
	int rsize;
	ref_file >> s;
	if (s[0]=='>') {
		string rname=s.substr(1,s.length()-1);
		cout << "Reference name " << rname << "\n";
		refseq.name=rname;
		//cout << rname << "\n";
		getline(ref_file,s);
	}
	ref_file >> rs;
	refseq.seq=rs;
	//cout << rs << "\n";
	//cout << ref << "\n";
	rsize=rs.size();
	//cout << rsize << "\n";
	refseq.size=rsize;
	getline(ref_file,s);
}

void ImportSamFileNames (run_params p, vector<string>& sam_files) {
	ifstream in_file;
	if (p.get_in==0) {
		p.in_file="Input_files.in";
	}
	in_file.open(p.in_file.c_str());
	for (int i=0;i<1000;i++) {
		string name;
		if (!(in_file >> name)) break;
		sam_files.push_back(name);
		cout << name << "\n";
	}
}

void ImportSamFileNamesFix (vector<string>& sam_files) {
    ifstream in_file;
    in_file.open("Input_files.in");
    for (int i=0;i<1000;i++) {
        string name;
        if (!(in_file >> name)) break;
        sam_files.push_back(name);
        cout << name << "\n";
    }
}


void ImportSamFile (run_params p, int i, rseq refseq, vector<string> sam_files, alldat& a) {
	ifstream sam_file;
	int s_length=0;
	sam_file.open(sam_files[i].c_str());
	vector<rd> data;
	s_length=0;
	ReadSamFile (p,sam_file,s_length,data,refseq);
	a.s_length=s_length-1;
	a.data=data;
}

void ReadSamFile (run_params p, ifstream& sam_file, int& s_length, vector<rd>& data, rseq refseq) {
	cout << "Reading .sam file... \n";
	string line;
	const string temp="@SQ";
	string lstore;
	int k=0;
	while (!sam_file.eof()) {
		k++;
		rd r;
		r.pairno=-1;
		for (int i=1;i<=11;i++) { //Columns of .sam file
			sam_file >> line;
			if (i==1) {
				lstore=line;
				string del=p.delimit;
				string tok=lstore.substr(0,lstore.find(del));
				if (lstore.compare(0,1,"@")!=0) {
					for (int j=1;j<=p.plines;j++) {  //Number of lines here depends on specific .sam file formatting
						lstore.erase(0,lstore.find(del)+1);
						tok=lstore.substr(0,lstore.find(del));
					}
					int k=0;
					while (lstore.compare(k,1,"0")==0||lstore.compare(k,1,"1")==0||lstore.compare(k,1,"2")==0||lstore.compare(k,1,"3")==0||lstore.compare(k,1,"4")==0||lstore.compare(k,1,"5")==0||lstore.compare(k,1,"6")==0||lstore.compare(k,1,"7")==0||lstore.compare(k,1,"8")==0||lstore.compare(k,1,"9")==0||lstore.compare(k,1,":")==0) {
						k++;
					}
					lstore=lstore.substr(0,k);
					r.paircode=lstore;
					//cout << r.paircode << "\n";
				}
			}
			if (i==2) {
				r.flag=atoi(line.c_str());
				BinDecomp(r.flag,r.flagv);
			}
			if (i==3) {
				r.ref=line;
				r.refno=1;
			}
			if (i==4) {
				r.alpos=atoi(line.c_str());
			//					cout << r.alpos << " ";
			}
			if (i==5) {
				r.alq=atoi(line.c_str());
			//					cout << r.alq << " ";
			}
			if (i==6) {
				r.cigar=line;
			//					cout << r.cigar << " ";
			}

			if (i==10) {
				r.seq=line;
			//					cout << r.seq << "\n";
			}
			if (i==11) {
				r.qual=line;
			//					cout << r.qual << "\n";
			}
		}
		r.rev=0;
		r.inc=0;
		r.pairno=-1;
		r.del=0; //Flag for deletion
		if (lstore.compare(0,1,"@")!=0) {
			s_length++;
			data.push_back(r);
		}
		getline(sam_file,line);
	}
	cout << "Number of reads " << s_length << "\n";
}

void BinDecomp (int i, vector<int>& v) {
	vector<int> vv (12,0);
	for (int j=11;j>=0;j--) {
		int p=pow(2,j);
		if (i>=p) {
			vv[j]=1;
			i=i-p;
		}
	}
	v=vv;
}


void OutputAlign (int i, vector<rd> data) {
	cout << "Output aligned reads...\n";
	ofstream al_file;
	ostringstream convert;
	convert << i;
	string temp=convert.str();
	string name = "Aligned"+temp+".out";
	cout << name << "\n";
	al_file.open(name.c_str());
	OutAl(data,al_file);
}

void OutAl (vector<rd> data, ofstream& al_file) {
	for (unsigned int i=0;i<data.size();i++) {
		if (data[i].del==0) {
			if (data[i].inc==1) {
				al_file << data[i].alpos << " " << data[i].seq << "\n";
			} else if (data[i].inc==-1) {
				al_file << data[i].alpos << " " << data[i].revseq << "\n";
			}
		}
	}
}

void InputConsensusSequences (vector<string> sam_files, vector<string>& consensus) {
	for (int i=0;i<sam_files.size();i++) {
		ifstream var_file;
		ostringstream convert;
		convert << i;
		string temp=convert.str();
		string name = "Consensus"+temp+".fa";
		var_file.open(name.c_str());
		string c1;
		string c2;
		var_file >> c1;
		var_file >> c2;
		consensus.push_back(c2);
	}
}


void InputVariantData (int output_consensus, vector<string> sam_files, vector< vector< vector<int> > >& all_nucs, vector< vector<char> >& all_cons) {
	cout << "Reading in variant files...\n";
	for (int i=0;i<sam_files.size();i++) {
		ifstream var_file;
		ostringstream convert;
		convert << i;
		string temp=convert.str();
		string name = "Variants"+temp+".out";
		string name2 = "Consensus"+temp+".fa";
		string name3 = ">Consensus_"+temp;
		//cout << name << "\n";
		var_file.open(name.c_str());
	
		vector< vector<int> > nucs;
		vector<char> cons;
		for (int j=0;j<1000000000;j++) {
			vector<int> ns;
			int n;
			char max='N';
			if (!(var_file >> n)) break;
			if (!(var_file >> n)) break;
			ns.push_back(n);
			if (n>0) {
				max='A';
			}
			if (!(var_file >> n)) break;
			ns.push_back(n);
			if (n>ns[0]) {max='C';}
			if (!(var_file >> n)) break;
			ns.push_back(n);
			if (n>ns[0]&&n>ns[1]) {
				max='G';
			}
			if (!(var_file >> n)) break;
			ns.push_back(n);
			if (n>ns[0]&&n>ns[1]&&n>ns[2]) {
				max='T';
			}
			if (!(var_file >> n)) break;
			ns.push_back(n);
			nucs.push_back(ns);
			cons.push_back(max);
		}
		if (output_consensus==1) {
			OutputParticularConsensus(name2,name3,cons);
		}
		all_nucs.push_back(nucs);
		all_cons.push_back(cons);
	}
}

void InputVariantTypes (vector< vector<int> >& var) {
	ifstream type_file;
	type_file.open("Variant_types.dat");
	int n;
	for (int i=0;i<1000000;i++) {
		vector<int> v;
		if (!(type_file >> n)) break;
		if (!(type_file >> n)) break;
		v.push_back(n);
		if (!(type_file >> n)) break;
		v.push_back(n);
		if (!(type_file >> n)) break;
		v.push_back(n);
		if (!(type_file >> n)) break;
		v.push_back(n);
		var.push_back(v);
	}
}

void OutputParticularConsensus (string name2, string name3, vector<char> cons) {
	ofstream cons_file;
	cons_file.open(name2.c_str());
	cons_file << name3 << "\n";
	for (int j=0;j<cons.size();j++) {
		cons_file << cons[j];
	}
	cons_file << "\n";
}

void OutputCodingConsensii (run_params p, vector<string>& sam_files, vector< vector <char> >& all_cons) {
	for (int i=0;i<sam_files.size();i++) {
		ifstream var_file;
		ostringstream convert;
		convert << i;
		string temp=convert.str();
		string name2 = "Consensus"+temp+"_cod.fa";
		string name3 = ">Consensus_"+temp;
		ofstream cons_file;
		cons_file.open(name2.c_str());
		cons_file << name3 << "\n";
		for (int j=p.trans_start;j<all_cons[i].size();j++) {
			cons_file << all_cons[i][j];
		}
		cons_file << "\n";
		cons_file.close();
	}
}

void OutputParticularAAConsensus (string name2, string name3, vector<char> aaseq) {
	ofstream cons_file;
	cons_file.open(name2.c_str());
	cons_file << name3 << "\n";
	for (int j=0;j<aaseq.size();j++) {
		cons_file << aaseq[j];
	}
	cons_file << "\n";
}

void OutputVarMask (string name4, int pre, int total, vector< vector<int> > s_types) {
	ofstream mask_file;
	mask_file.open(name4.c_str());
	//Variants before the start of the coding region - give these -1.
	pre=pre-1;
	for (int i=0;i<pre;i++) {
		mask_file << "-1 -1 -1 -1\n";
	}
	for (int i=0;i<s_types.size();i++) {
		for (int j=0;j<s_types[i].size();j++) {
			mask_file << s_types[i][j] << "  ";
		}
		mask_file << "\n";
	}
	int remain=total-pre-s_types.size();
	for (int i=0;i<remain;i++) {
		mask_file << "-1 -1 -1 -1\n";
	}
	mask_file.close();
}


void OutputGlobalConsensus (vector< vector< vector<int> > >& all_nucs, vector<char>& global_cons) {
	cout << "Output Global Consensus\n";
	cout << all_nucs[0].size() << "\n";
	ofstream cons_file;
	cons_file.open("Consensus_all.fa");
	cons_file << ">Consensus_all\n";
	vector<char> total_cons;
	for (int i=0;i<all_nucs[0].size();i++) {
		vector<int> ns;
		for (int j=0;j<4;j++) {
			int tot=0;
			for (int k=0;k<all_nucs.size();k++) {
				tot=tot+all_nucs[k][i][j];
			}
			ns.push_back(tot);
		}
		char max='N';
		if (ns[0]>0) {
			max='A';
		}
		if (ns[1]>ns[0]) {
			max='C';
		}
		if (ns[2]>ns[0]&&ns[2]>ns[0]) {
			max='G';
		}
		if (ns[3]>ns[2]&&ns[3]>ns[1]&&ns[3]>ns[0]) {
			max='T';
		}
		cons_file << max;
		global_cons.push_back(max);
	}
	cons_file << "\n";
}

void GetVariantTotals (vector<string> sam_files, vector< vector<int> >& all_tots) {
	for (int i=0;i<sam_files.size();i++) {
		ifstream var_file;
		ostringstream convert;
		convert << i;
		string temp=convert.str();
		string name = "Variants"+temp+".out";
		cout << name << "\n";
		var_file.open(name.c_str());
		vector<int> ns;
		for (int j=0;j<1000000000;j++) {
			int n;
			if (!(var_file >> n)) break;
			if (!(var_file >> n)) break;
			if (!(var_file >> n)) break;
			if (!(var_file >> n)) break;
			if (!(var_file >> n)) break;
			if (!(var_file >> n)) break;
			ns.push_back(n);
		}
		all_tots.push_back(ns);
	}
}

void InputJoinedData (int n_times, vector< vector<joined> >& t_reads) {
	for (int i=0;i<n_times;i++) {
		ifstream jp_file;
		ostringstream convert;
		convert << i;
		string temp=convert.str();
		string name = "Joined"+temp+".out";
		jp_file.open(name.c_str());
		int x;
		string s;
		vector<joined> reads;
		for (int j=0;j<=100000000;j++) {
			joined jn;
			if (!(jp_file >> x)) break;
			if (!(jp_file >> s)) break;
			jn.alpos=x;
			jn.seq=s;
			reads.push_back(jn);
		}
		t_reads.push_back(reads);
	}
}

void InputJnData (int i, vector<joined>& t_read) {
	ifstream jp_file;
	ostringstream convert;
	convert << i;
	string temp=convert.str();
	string name = "Joined"+temp+".out";
	cout << "Read data from " << name << "\n";
	jp_file.open(name.c_str());
	int x;
	string s;
	for (int j=0;j<=100000000;j++) {
		joined jn;
		if (!(jp_file >> x)) break;
		if (!(jp_file >> s)) break;
		jn.alpos=x;
		jn.seq=s;
		t_read.push_back(jn);
	}
}

void InputJnDataGeneral (string name, vector<joined>& t_read) {
	ifstream jp_file;
	jp_file.open(name.c_str());
	int x;
	string s;
	for (int j=0;j<=100000000;j++) {
		joined jn;
		if (!(jp_file >> x)) break;
		if (!(jp_file >> s)) break;
		jn.alpos=x;
		jn.seq=s;
		t_read.push_back(jn);
	}
}

void InputQualData (run_params p, int i, vector<ql>& q_read, vector<char> qual) {
	ifstream jp_file;
	ostringstream convert;
	convert << i;
	string temp=convert.str();
	string name = "Quality"+temp+".out";
	cout << "Read data from " << name << "\n";
	jp_file.open(name.c_str());
	int x;
	string s;
	
	for (int j=0;j<=100000000;j++) {
		ql q;
		if (!(jp_file >> x)) break;
		if (!(jp_file >> s)) break;
		q.alpos=x;
		for (int k=0;k<s.size();k++) {
			for (int l=0;l<qual.size();l++) {
				if (s[k]==qual[l]) {
					q.seq.push_back(l);
					break;
				}
			}
		}
		q_read.push_back(q);
	}
	//int max=0;
	if (p.qlib==3) {
		for (int i=0;i<q_read.size();i++) {
			for (int j=0;j<q_read[i].seq.size();j++) {
				if (q_read[i].seq[j]>0&&q_read[i].seq[j]<93) {
					q_read[i].seq[j]=floor(10.*log10(q_read[i].seq[j]));
					/*if (q_read[i].seq[j]>max) {
					 max=q_read[i].seq[j];
					 cout << "Max qual " << max << "\n";
					 }*/
				} else {
					q_read[i].seq[j]=0;
				}
			}
		}
	}
}



void OutputVarFile (int i, rseq refseq, nuc r_count) {
	ostringstream convert;
	convert << i;
	string temp=convert.str();
	string name = "Variants"+temp+".out";
	ofstream sl_file;
	sl_file.open(name.c_str());
	for (int j=1;j<=refseq.size;j++) {
		sl_file << j << " " << r_count.nA[j] << " " << r_count.nC[j] << " " << r_count.nG[j] << " " << r_count.nT[j] << " " << r_count.nN[j] << "\n";
	}
	sl_file.close();
}


void OutputSLTData (run_params p, int sel, const char* filename, vector<str> sltrajs) {
	if (sel!=0&&sel!=1) {
		cout << "Error in call of OutputSLTData: Value must be 0 or 1\n";
	}
	ofstream slt_file;
	slt_file.open(filename);
	int min=1e9;
	int m_index=-1;
	int tot=0;
	for (unsigned int i=0;i<sltrajs.size();i++) {
		if (sltrajs[i].inc==sel) {
			tot++;
		}
	}
	cout << "Total " << tot << "\n";
	
	vector<int> seen;
	
	for (int i=0;i<tot;i++) {
		//cout << "T " << i << "\n";
		min=1e9;
		for (unsigned int j=0;j<sltrajs.size();j++) {
			if (sltrajs[j].locus<min&&sltrajs[j].inc==sel) {
				min=sltrajs[j].locus;
				m_index=j;
			}
		}

		sltrajs[m_index].inc=1-sel;
		char cons=sltrajs[m_index].cons;
		
		if (cons!=sltrajs[m_index].nuc) {
			int s=0;
			if (p.uniq==1) {
				for (int j=0;j<seen.size();j++) {
					if (seen[j]==sltrajs[m_index].locus) {
						s=1;
					}
				}
			}
			if (s==0) {
				slt_file << sltrajs[m_index].locus << " " << cons << " " << sltrajs[m_index].nuc << " ";
				slt_file <<sltrajs[m_index].times.size() << " ";
				for (unsigned int j=0;j<sltrajs[m_index].times.size();j++) {
					slt_file << sltrajs[m_index].times[j] << " ";
					slt_file << sltrajs[m_index].nA[j] << " ";
					slt_file << sltrajs[m_index].nC[j] << " ";
					slt_file << sltrajs[m_index].nG[j] << " ";
					slt_file << sltrajs[m_index].nT[j] << " ";
					slt_file << sltrajs[m_index].nN[j] << " ";
				}
				slt_file << "\n";
				seen.push_back(sltrajs[m_index].locus);
			}
		}
	}
}

void ImportSLTData (run_params p, vector<str>& sltrajs) {
	ifstream slt_file;
	slt_file.open(p.in_file.c_str());
	cout << "Opened file " << p.in_file << "\n";
	char c;
	int n;
	int t;
	for (int i=0;i<100000;i++) {
		str s;
		if (!(slt_file >> n)) break;
		s.locus=n;
		slt_file >> c;
		s.cons=c;
		slt_file >> c;
		s.nuc=c;
		slt_file >> t;
		for (int j=0;j<t;j++) {
			slt_file >> n;
			s.times.push_back(n);
			slt_file >> n;
			s.nA.push_back(n);
			slt_file >> n;
			s.nC.push_back(n);
			slt_file >> n;
			s.nG.push_back(n);
			slt_file >> n;
			s.nT.push_back(n);
			slt_file >> n;
			s.nN.push_back(n);
		}
		s.inc=1;
		sltrajs.push_back(s);
		if (p.verb==1) {
			cout << s.locus << " " << s.cons << " " << s.nuc << " " << s.times.size() << " ";
			for (int j=0;j<s.times.size();j++) {
				cout << s.times[j] << " " << s.nA[j] << " " << s.nC[j] << " " << s.nG[j] << " " << s.nT[j] << " " << s.nN[j] << " ";
			}
			cout << "\n";
		}
	}
}

void ImportSLTLoci (const char* filename, int& check, vector<int>& loci) {
	//Read in a list of loci from a trajectory file
	cout << "Import \n";
	ifstream slt_file;
	slt_file.open(filename);
	cout << "Opened file " << filename << "\n";
	string ss;
	while (getline(slt_file,ss)) {
		int loc;
		istringstream in(ss);
		in >> loc;
		loci.push_back(loc);
	}
	if (loci.size()==0) {
		cout << "Error: No single-locus trajectories in input file\n";
		check=1;
	}
}

void ImportTimeData (vector<int>& times) {
	ifstream t_file;
	t_file.open("Times.in");
	int t;
	for (int i=0;i<100000;i++) {
		if (!(t_file >> t)) break;
		times.push_back(t);
	}
}



void InputHaplotypeData (run_params p, vector<int> polys, vector< vector<char> >& haps) {
	ifstream hap_file;
	cout << "Read haplotypes from " << p.hap_file.c_str() << "\n";
	hap_file.open(p.hap_file.c_str());
	for (int i=0;i<100000;i++) {
		char c;
		vector<char> h;
		if (!(hap_file >> c)) break;
		for (unsigned int j=0;j<polys.size()-1;j++) {
			h.push_back(c);
			hap_file >> c;
		}
		h.push_back(c);
		haps.push_back(h);
	}
	if (p.verb==1) {
		cout << "Imported external haplotypes\n";
		for (unsigned int i=0;i<haps.size();i++) {
			for (unsigned int j=0;j<haps[i].size();j++) {
				cout << haps[i][j];
			}
			cout << "\n";
		}
	}
}



void GenerateOutputFilesNoHaps (run_params p, vector<int> polys, vector< vector< vector<mpoly> > >& m_polys, vector< vector<mpoly> > c_m_polys) {
	//Don't produce Contribs data, or X-described calls
	int index=0;
	cout << "Building output files...\n";
	ofstream ah_file;
	vector<int> times;
	if (p.full_rep==1) {
		ah_file.open("Multi_locus_trajectories.out");
		ImportTimeData(times);
	}
	
	//Elements of m_polys are ordered by time, then subset of read types
	//Elements of c_m_polys are ordered by subset of read types
	
	for (unsigned int i=0;i<c_m_polys.size();i++) {
		//cout << "i= " << i << " " << c_m_polys[i].size() << "\n";
		//	//Go through partial haplotypes in this set
		vector<mtr> mhap;
		if (c_m_polys[i].size()>0) {
			//cout << "Check types\n";
			ofstream haps_file;
			ofstream locus_file;
			ostringstream convert;
			string name_h;
			string name_l;
			if (p.hap_index==1) {
				convert << i;				//Include all classes in the labelling
				string temp=convert.str();
				name_h = "Hap_data"+temp+".dat";
				name_l = "Loci"+temp+".dat";
			} else {
				convert << index;				//Include only observed classes in the labelling
				string temp=convert.str();
				name_h = "Hap_data"+temp+".dat";
				name_l = "Loci"+temp+".dat";
			}
			index++;
			if (p.full_rep==0) {
				haps_file.open(name_h.c_str());
				locus_file.open(name_l.c_str());
			}
			
			int done_loc=0;
			vector<int> x_vect (m_polys.size(),0);  //Collect numbers of non-included haplotypes
			for (unsigned int j=0;j<c_m_polys[i].size();j++) {
				//Output partial haplotype name and relevant loci
				if (p.full_rep==1) {
					vector<int> loci;
					for (unsigned int l=0;l<c_m_polys[i][j].vars.size();l++) {
						if (c_m_polys[i][j].vars[l]!='-') {
							loci.push_back(l);
							//ah_file << l << " ";
						}
					}
					ah_file << loci.size() << " ";
					for (int q=0;q<loci.size();q++) {
						ah_file << polys[loci[q]] << " ";
					}
				}
				
				for (unsigned int l=0;l<c_m_polys[i][j].vars.size();l++) {
					if (c_m_polys[i][j].vars[l]!='-') {
						if (done_loc<10) {
							if (p.full_rep==0) {
								locus_file << l << " ";
							}
							done_loc=1;
						}
						if (p.full_rep==0) {
							haps_file << c_m_polys[i][j].vars[l];
						}
						if (p.full_rep==1) {
							ah_file << c_m_polys[i][j].vars[l];
						}
					}
				}
				if (done_loc==1) {
					locus_file << "\n";
					done_loc=10;
				}
				if (p.full_rep==0) {
					haps_file << " ";
				}
				if (p.full_rep==1) {
					ah_file << " " << times.size() << " ";
				}
				//Find number of matching haplotypes with time
				int t_index=0;
				for (unsigned int ii=0;ii<m_polys.size();ii++) { //Time point
					//cout << "Time " << ii << "\n";
					//Known class is i
					int found=0;
					for (unsigned int jj=0;jj<m_polys[ii][i].size();jj++) { //Partial haplotypes at time ii
						int match=1;
						for (unsigned int kk=0;kk<m_polys[ii][i][jj].vars.size();kk++) {
							if (c_m_polys[i][j].vars[kk]!=m_polys[ii][i][jj].vars[kk]) {
								match=0;
							}
						}
						if (match==1) {
							if (p.full_rep==0) {
								haps_file << m_polys[ii][i][jj].count << " ";
							}
							if (p.full_rep==1) {
								ah_file << times[t_index] << " ";
								ah_file << m_polys[ii][i][jj].count << " ";
								t_index++;
							}
							found=1;
						}
					}
					if (found==0) { //No identical haplotype identified at this time-point
						if (p.full_rep==0) {
							haps_file << "0 ";
						}
						if (p.full_rep==1) {
							ah_file << times[t_index] << " ";
							ah_file << "0 ";
							t_index++;
						}
					}
				}
				if (p.full_rep==0) {
					haps_file << "\n";
				}
				if (p.full_rep==1) {
					ah_file << "\n";
				}
			}
				
		}
	}
}


void GenerateOutputFiles (run_params p, vector<int> polys, vector< vector<char> > haps, vector< vector< vector<mpoly> > > m_polys, vector< vector<mpoly> > c_m_polys, vector< vector<mtr> >& mltrajs) {
	int index=0;
	cout << "Building output files...\n";
	ofstream ah_file;
	vector<int> times;
	if (p.full_rep==1) {
		ah_file.open("Multi_locus_trajectories.out");
		ImportTimeData(times);
	}
	for (unsigned int i=0;i<c_m_polys.size();i++) {
		//	//Go through partial haplotypes in this set
		if (c_m_polys[i].size()>0) {
			ofstream haps_file;
			ofstream contribs_file;
			ofstream locus_file;
			ostringstream convert;
			string name_h;
			string name_c;
			string name_l;
			if (p.hap_index==1) {
				convert << i;				//Include all classes in the labelling
				string temp=convert.str();
				name_h = "Hap_data"+temp+".dat";
				name_c = "Contribs"+temp+".dat";
				name_l = "Loci"+temp+".dat";
			} else {
				convert << index;				//Include all classes in the labelling
				string temp=convert.str();
				name_h = "Hap_data"+temp+".dat";
				name_c = "Contribs"+temp+".dat";
				name_l = "Loci"+temp+".dat";
			}
			index++;
			int done_loc=0;
			if (p.full_rep==0) {
				haps_file.open(name_h.c_str());
				locus_file.open(name_l.c_str());
				contribs_file.open(name_c.c_str());
			}
			vector<int> x_vect (m_polys.size(),0);  //Collect numbers of non-included haplotypes
			int pr=0;
			
			for (unsigned int j=0;j<c_m_polys[i].size();j++) {
				
				if (p.full_rep==1) {
					vector<int> loci;
					for (unsigned int l=0;l<c_m_polys[i][j].vars.size();l++) {
						if (c_m_polys[i][j].vars[l]!='-') {
							loci.push_back(l);
							//ah_file << l << " ";
						}
					}
					ah_file << loci.size() << " ";
					for (int q=0;q<loci.size();q++) {
						ah_file << polys[loci[q]] << " ";
					}
					
				}
				
				vector<int> matches;
				pr=0;
				for	(unsigned int k=0;k<haps.size();k++) {
					int mmatch=1;
					for (int l=0;l<c_m_polys[i][j].vars.size();l++) {
						if (c_m_polys[i][j].vars[l]!='-'&&c_m_polys[i][j].vars[l]!=haps[k][l]) {
							mmatch=0;
						}
					}
					if (mmatch==1) {//match haplotype
						pr=1;
						matches.push_back(k);
					}
				}
				if (matches.size()>0) { //Matches some haplotype
					//Output time-dependent frequencies to Hap_data?
					for (unsigned int l=0;l<c_m_polys[i][j].vars.size();l++) { //Partial haplotype name
						if (c_m_polys[i][j].vars[l]!='-') {
							if (done_loc<10) {
								if (p.full_rep==0) {
									locus_file << l << " ";
								}
								done_loc=1;
							}
							if (p.full_rep==0) {
								haps_file << c_m_polys[i][j].vars[l];
							}
							if (p.full_rep==1) {
								ah_file << c_m_polys[i][j].vars[l];
							}
						}
					}
					if (done_loc==1) {
						if (p.full_rep==0) {
							locus_file << "\n";
						}
						done_loc=10;
					}
					if (p.full_rep==0) {
						haps_file << " ";
					}
					if (p.full_rep==1) {
						ah_file << " " << times.size() << " ";
					}

					//Find number of matching haplotypes with time
					int t_index=0;
					for (unsigned int ii=0;ii<m_polys.size();ii++) { //Time point
						//Known class is i
						int found=0;
						for (unsigned int jj=0;jj<m_polys[ii][i].size();jj++) { //Partial haplotypes at time ii
							int match=1;
							for (unsigned int kk=0;kk<m_polys[ii][i][jj].vars.size();kk++) {
								if (c_m_polys[i][j].vars[kk]!=m_polys[ii][i][jj].vars[kk]) {
									match=0;
								}
							}
							if (match==1) {
								if (p.full_rep==0) {
									haps_file << m_polys[ii][i][jj].count << " ";
								}
								if (p.full_rep==1) {
									ah_file << times[t_index] << " ";
									ah_file << m_polys[ii][i][jj].count << " ";
									t_index++;
								}
								found=1;
							}
						}
						if (found==0) { //No identical haplotype identified at this time-point
							if (p.full_rep==0) {
								haps_file << "0 ";
							}
							if (p.full_rep==1) {
								ah_file << times[t_index] << " ";
								ah_file << "0 ";
								t_index++;
							}
						}
					}
					if (p.full_rep==0) {
						haps_file << "\n";
					}
					if (p.full_rep==1) {
						ah_file << "\n";
					}
					//Output hap matches to Contribs
					for (unsigned int k=0;k<matches.size();k++) {
						contribs_file << matches[k] << " ";
					}
					contribs_file << "\n";
				} else { //No matches - count the number of occurrences of this haplotype and add to final x vector
					for (unsigned int ii=0;ii<m_polys.size();ii++) { //Time point
						//Known class is i
						for (unsigned int jj=0;jj<m_polys[ii][i].size();jj++) { //Partial haplotypes at time ii
							int match=1;
							for (unsigned int kk=0;kk<m_polys[ii][i][jj].vars.size();kk++) {
								if (c_m_polys[i][j].vars[kk]!=m_polys[ii][i][jj].vars[kk]) {
									match=0;
								}
							}
							if (match==1) {
								x_vect[ii]=x_vect[ii]+m_polys[ii][i][jj].count;
							}
						}
					}
				}
			}
			if (p.full_rep==0) {
				haps_file << "X ";
			}
			if (p.full_rep==1) {
				if (p.printx==1) {
					vector<int> loci;
					for (unsigned int l=0;l<c_m_polys[i][0].vars.size();l++) {
						if (c_m_polys[i][0].vars[l]!='-') {
							loci.push_back(l);
						}
					}
					ah_file << loci.size() << " ";
					for (int q=0;q<loci.size();q++) {
						ah_file << polys[loci[q]] << " ";
					}
					ah_file << "X ";
				}
			}
			int t_index=0;
			for (unsigned int k=0;k<x_vect.size();k++) {
				if (p.full_rep==0) {
					haps_file << x_vect[k] << " ";
				}
				if (p.full_rep==1) {
					if (p.printx==1) {
						ah_file << times[t_index] << " ";
						ah_file << x_vect[k] << " ";
						t_index++;
					}
				}
			}
			if (p.full_rep==0) {
				haps_file << "\n";
			}
			if (p.full_rep==1) {
				if (p.printx==1) {
					ah_file << "\n";
				}
			}
		}
	}
}

void ReadMultiMLTrajs (run_params p, vector< vector<int> >& times, vector< vector<mtr> >& mltrajs) {
	ifstream ml_file;
	ifstream mlt_file;
	ifstream dat_file;
	ifstream tim_file;
	ml_file.open("Partial_hap_list.dat");
	mlt_file.open("Partial_hap_times.dat");
	for (int i=0;i<=10000;i++) {
		//Haplotype data
		vector<mtr> mt;
		string s;
		if (!(ml_file >> s)) break;
		dat_file.open(s.c_str());
		string ss;
		while (getline(dat_file,ss)) {
			mtr m;
			istringstream in(ss);
			string sss;
			in >> sss;
			int n;
			int tot=0;
			while (in >> n) {
				m.n.push_back(n);
				tot=tot+n;
			}
			if (tot>p.hap_n_min) {  //Exclude types with too few observations across time; may occur if calling is done to external haplotypes
				mt.push_back(m);
			}
		}
		
		//Time data
		mlt_file >> s;
		tim_file.open(s.c_str());
		int n;
		vector<int> t;
		for (int j=0;j<100000;j++) {
			if (!(tim_file >> n)) break;
			t.push_back(n);
		}
		dat_file.close();
		tim_file.close();
		mltrajs.push_back(mt);
		times.push_back(t);
	}
}

void ReadMultiLocTrajs (run_params p, vector< vector<int> >& times, vector< vector<mtr> >& mltrajs) {
	if (p.get_in==0) {
		p.in_file="Multi_locus_trajectories.out";
	}
	ifstream ml_file;
	ml_file.open(p.in_file.c_str());
	vector<int> loc_ref;
	vector<mtr> mt;
	vector<int> t;
	int g_time=0;
	for (int i=0;i<=1000000;i++) {
		mtr m;
		//Haplotype data
		int n_loc;
		if (!(ml_file >> n_loc)) break;
		int k;
		vector<int> loci;
		for (int j=0;j<n_loc;j++) {
			ml_file >> k;
			loci.push_back(k);
		}
		if (loci!=loc_ref) {
			if (mt.size()>1) {
				mltrajs.push_back(mt);
				times.push_back(t);
			}
			loc_ref=loci;
			g_time=0;
			t.clear();
			mt.clear();
		}
		string s;
		ml_file >> s;
		int n_times;
		ml_file >> n_times;
		int tot=0;
		for (int j=0;j<n_times;j++) {
			int n;
			ml_file >> n;
			if (g_time==0) {
				t.push_back(n);
			}
			ml_file >> n;
			m.n.push_back(n);
			tot=tot+n;
		}
		g_time=1;
		mt.push_back(m);
	}
	if (p.verb==1) {
		cout << "Number of trajectories " << mltrajs.size() << "\n";
	}
}

void ReadMultiLocTrajs2 (run_params p, vector<str> sltrajs, vector<ld_info>& ld_data) {
	ifstream ml_file;
	ml_file.open("Multi_locus_trajectories.out");
	for (int i=0;i<=1000000;i++) {
		vector<int> loci;
		int n_loc;
		int k;
		if (!(ml_file >> n_loc)) break;
		for (int j=0;j<n_loc;j++) {
			ml_file >> k;
			loci.push_back(k);
		}
		string s;
		ml_file >> s;
		int n_times;
		ml_file >> n_times;
		vector<int> samples;
		for (int j=0;j<n_times;j++) {
			int n;
			ml_file >> n;
			ml_file >> n;
			samples.push_back(n);
		}
		
		if (i==0) {
			for (int j=0;j<ld_data.size();j++) {
				for (int k=0;k<n_times;k++) {
					ld_data[j].n_i1.push_back(0);
					ld_data[j].n_i0.push_back(0);
					ld_data[j].n_j1.push_back(0);
					ld_data[j].n_j0.push_back(0);
					ld_data[j].n_11.push_back(0);
					ld_data[j].n_10.push_back(0);
					ld_data[j].n_01.push_back(0);
					ld_data[j].n_00.push_back(0);

				}
			}
		}
		
	//	cout << "# loci " << n_loc << "\n" << "Loci: ";
		//Two locus measurements
		for (int j=0;j<n_loc;j++) {
			for (int k=0;k<n_loc;k++) {
				if (j!=k) {
				//	cout << "j " << loci[j] << " k " << loci[k] << "\n";
					//Find equivalent marker
					int mark=-1;
					for (int m=0;m<ld_data.size();m++) {
						if (ld_data[m].i==loci[j]&&ld_data[m].j==loci[k]) {
							mark=m;
							break;
						}
					}
				//	cout << "Mark " << mark << "\n";
					int j_mark=0;
					int k_mark=0;
					for (int m=0;m<sltrajs.size();m++) {
						if (sltrajs[m].locus==loci[j]) {
							if (s[j]==sltrajs[m].nuc) {
								j_mark=1;
							}
						}
					}
					for (int m=0;m<sltrajs.size();m++) {
						if (sltrajs[m].locus==loci[k]) {
							if (s[k]==sltrajs[m].nuc) {
								k_mark=1;
							}
						}
					}
					for (int l=0;l<samples.size();l++) {
						if (j_mark==1) {
							if (k_mark==1) {
								ld_data[mark].n_11[l]=ld_data[mark].n_11[l]+samples[l];
							} else {
								ld_data[mark].n_10[l]=ld_data[mark].n_10[l]+samples[l];
							}
						} else {
							if (k_mark==1) {
								ld_data[mark].n_01[l]=ld_data[mark].n_01[l]+samples[l];
							} else {
								ld_data[mark].n_00[l]=ld_data[mark].n_00[l]+samples[l];
							}
						}
					}
				}
			}
		}
		
		//Loci that have not been seen
		vector<int> not_seen;
		for (int j=0;j<sltrajs.size();j++) {
			int seen=0;
			for (int k=0;k<loci.size();k++) {
				if (sltrajs[j].locus==loci[k]) {
					seen=1;
				}
			}
			if (seen==0) {
				not_seen.push_back(sltrajs[j].locus);
			}
		}
		
		for (int j=0;j<n_loc;j++) {
			for (int k=0;k<not_seen.size();k++) {
				//Find equivalent marker
			//	cout << "2 j " << loci[j] << " k " << not_seen[k] << "\n";

				int mark=-1;
				for (int m=0;m<ld_data.size();m++) {
					if (ld_data[m].i==loci[j]&&ld_data[m].j==not_seen[k]) {
						mark=m;
						break;
					}
				}
				int j_mark=0;
				for (int m=0;m<sltrajs.size();m++) {
					if (sltrajs[m].locus==loci[j]) {
						if (s[j]==sltrajs[m].nuc) {
							j_mark=1;
						}
					}
				}
				for (int l=0;l<samples.size();l++) {
					if (j_mark==1) {
						ld_data[mark].n_i1[l]=ld_data[mark].n_i1[l]+samples[l];
					} else {
						ld_data[mark].n_i0[l]=ld_data[mark].n_i0[l]+samples[l];
					}
				}
			}
		}
		
		for (int j=0;j<not_seen.size();j++) {
			for (int k=0;k<loci.size();k++) {
				int mark=-1;
				for (int m=0;m<ld_data.size();m++) {
					if (ld_data[m].i==not_seen[j]&&ld_data[m].j==loci[k]) {
						mark=m;
						break;
					}
				}
				int k_mark=0;
				for (int m=0;m<sltrajs.size();m++) {
					if (sltrajs[m].locus==loci[k]) {
						if (s[k]==sltrajs[m].nuc) {
							k_mark=1;
						}
					}
				}

				for (int l=0;l<samples.size();l++) {
					if (k_mark==1) {
						ld_data[mark].n_j1[l]=ld_data[mark].n_j1[l]+samples[l];
					} else {
						ld_data[mark].n_j0[l]=ld_data[mark].n_j0[l]+samples[l];
					}
				}
			}
		}
	}
}

void PrintSeqLens (int i, vector<int> l_dist) {
	ofstream stat_file;
	ostringstream convert;
	convert << i;
	string temp=convert.str();
	string name = "Seq_lengths"+temp+".out";
	stat_file.open(name.c_str());
	for (int k=0;k<l_dist.size();k++) {
		stat_file << k << " " << l_dist[k] << "\n";
	}
	stat_file.close();
}

void PrintSeqQual (int i, vector<long> q_dist) {
	ofstream stat_file;
	ostringstream convert;
	convert << i;
	string temp=convert.str();
	string name = "Base_quality"+temp+".out";
	stat_file.open(name.c_str());
	int max=0;
	for (int j=0;j<q_dist.size();j++) {
		if (q_dist[j]>0) {
			max=j;
		}
	}
	for (int j=0;j<=max;j++) {
		stat_file << j << " " << q_dist[j] << "\n";
	}
	stat_file.close();
}


void OutputDistancesCounts (vector< vector<double> >& distances, vector< vector<double> >& all_counts) {
	ofstream dist_file;
	ofstream counts_file;
	dist_file.open("Distances.out");
	counts_file.open("Counts.out");
	for (int i=0;i<distances.size();i++) {
		for (int j=0;j<distances[i].size();j++) {
			dist_file << distances[i][j] << " ";
		}
		dist_file << "\n";
	}
	for (int i=0;i<all_counts.size();i++) {
		for (int j=0;j<all_counts[i].size();j++) {
			counts_file << all_counts[i][j] << " ";
		}
		counts_file << "\n";
	}
}

void OutputDistancesCountsSNS (vector< vector<double> >& distances_n, vector< vector<double> >& distances_s, vector< vector<double> >& all_counts_n, vector< vector<double> >& all_counts_s) {
	ofstream dist_n_file;
	ofstream dist_s_file;
	ofstream counts_n_file;
	ofstream counts_s_file;
	dist_n_file.open("Distances_n.out");
	dist_s_file.open("Distances_s.out");
	counts_n_file.open("Counts_n.out");
	counts_s_file.open("Counts_s.out");
	for (int i=0;i<distances_n.size();i++) {
		for (int j=0;j<distances_n[i].size();j++) {
			dist_n_file << distances_n[i][j] << " ";
		}
		dist_n_file << "\n";
	}
	for (int i=0;i<distances_s.size();i++) {
		for (int j=0;j<distances_s[i].size();j++) {
			dist_s_file << distances_s[i][j] << " ";
		}
		dist_s_file << "\n";
	}
	for (int i=0;i<all_counts_n.size();i++) {
		for (int j=0;j<all_counts_n[i].size();j++) {
			counts_n_file << all_counts_n[i][j] << " ";
		}
		counts_n_file << "\n";
	}
	for (int i=0;i<all_counts_s.size();i++) {
		for (int j=0;j<all_counts_s[i].size();j++) {
			counts_s_file << all_counts_s[i][j] << " ";
		}
		counts_s_file << "\n";
	}
}

void ImportVariantMasks (vector<string>& sam_files, vector< vector< vector<int> > >& all_types) {
	for (int i=0;i<sam_files.size();i++) {
		ifstream var_file;
		ostringstream convert;
		convert << i;
		string temp=convert.str();
		string name = "Variant_mask"+temp+".out";
		var_file.open(name.c_str());
		vector< vector<int> > all_t;
		for (int j=0;j<1000000000;j++) {
			vector<int> type;
			int t;
			if (!(var_file >> t)) break;
			type.push_back(t);
			if (!(var_file >> t)) break;
			type.push_back(t);
			if (!(var_file >> t)) break;
			type.push_back(t);
			if (!(var_file >> t)) break;
			type.push_back(t);
			all_t.push_back(type);
		}
		all_types.push_back(all_t);
	}
}

void OutputDeconJoined (run_params p, vector<joined> t_read1, vector<joined> t_read2, vector<int> filter1, vector<int> filter2) {
	ofstream filt_file1;
	ofstream filt_file2;
	string ff1=p.jn1+".dec";
	string ff2=p.jn2+".dec";
	filt_file1.open(ff1.c_str());
	filt_file2.open(ff2.c_str());
	int index=0;
	for (int i=0;i<t_read1.size();i++) {
		if (filter1[index]==i) {
			index++;
		} else {
			filt_file1 << t_read1[i].alpos << " " << t_read1[i].seq << "\n";
		}
	}
	index=0;
	for (int i=0;i<t_read2.size();i++) {
		if (filter2[index]==i) {
			index++;
		} else {
			filt_file2 << t_read2[i].alpos << " " << t_read2[i].seq << "\n";
		}
	}
}

void OutputLDRawInformation (int t, vector<ld_info>& ld_data) {
    ofstream ldf_file;
    ostringstream convert;
    convert << t;
    string temp=convert.str();
    string name = "LD_freqs"+temp+".out";
    ldf_file.open(name.c_str());
    for (int i=0;i<ld_data.size();i++) {
        ldf_file << ld_data[i].i  << " " << ld_data[i].j << " " << ld_data[i].n_i1[0] << " " << ld_data[i].n_i0[0] << " " << ld_data[i].n_j1[0] << " " << ld_data[i].n_j0[0] << " " << ld_data[i].n_11[0] << " " << ld_data[i].n_10[0] << " " << ld_data[i].n_01[0] << " " << ld_data[i].n_00[0] << "\n";
    }
}
