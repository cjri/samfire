#include "shared2.h"
#include <iostream>
#include <string>

void GetOptions (run_params& p, int argc, const char **argv) {
	string p_switch;
	p.read=0;
	p.unobs=0;
	p.seed=0;
	p.timefile=0;
	p.mcmc=0;
	p.mu=0.333333333333e-5;
	p.c=60.1155;
	p.seed=atoi(argv[1]);
	p.data=argv[2];
	int x=3;
	while (x < argc && (argv[x][0]=='-')) {
		p_switch=argv[x];
		if (p_switch.compare("--mu")==0) {
			x++;
			p.mu=atof(argv[x]);
		} else if (p_switch.compare("--c")==0) {
			x++;
			p.c=atof(argv[x]);
		} else if (p_switch.compare("--read")==0) {
			x++;
			p.read=atoi(argv[x]);
		} else if (p_switch.compare("--mcmc")==0) {
			x++;
			p.mcmc=atoi(argv[x]);
		} else if (p_switch.compare("--unobs")==0) {
			x++;
			p.unobs=atoi(argv[x]);
		} else if (p_switch.compare("--readtimes")==0) {
			x++;
			p.timefile=atoi(argv[x]);
		}
		else {
			cout << "Incorrect usage\n ";
			exit(1);
		}
		p_switch.clear();
		x++;
	}
}

void GetOptionsDual (run_params& p, int argc, const char **argv) {
	string p_switch;
	p.read=0;
	p.unobs=0;
	p.seed=0;
	p.timefile=0;
	p.mu=0.333333333333e-5;
	p.c=60.1155;
	p.seed=atoi(argv[1]);
	p.data1=argv[2];
	p.data2=argv[3];
	int x=4;
	while (x < argc && (argv[x][0]=='-')) {
		p_switch=argv[x];
		if (p_switch.compare("--mu")==0) {
			x++;
			p.mu=atof(argv[x]);
		} else if (p_switch.compare("--c")==0) {
			x++;
			p.c=atof(argv[x]);
		} else if (p_switch.compare("--read")==0) {
			x++;
			p.read=atoi(argv[x]);
		} else if (p_switch.compare("--unobs")==0) {
			x++;
			p.unobs=atoi(argv[x]);
		} else if (p_switch.compare("--readtimes")==0) {
			x++;
			p.timefile=atoi(argv[x]);
		}
		else {
			cout << "Incorrect usage\n ";
			exit(1);
		}
		p_switch.clear();
		x++;
	}
}

void GetOptionsThree (run_params& p, int argc, const char **argv) {
	string p_switch;
	p.read=0;
	p.unobs=0;
	p.seed=0;
	p.timefile=0;
	p.mu=0.333333333333e-5;
	p.c=60.1155;
	p.seed=atoi(argv[1]);
	p.data1=argv[2];
	p.data2=argv[3];
	p.data3=argv[4];
	int x=5;
	while (x < argc && (argv[x][0]=='-')) {
		p_switch=argv[x];
		if (p_switch.compare("--mu")==0) {
			x++;
			p.mu=atof(argv[x]);
		} else if (p_switch.compare("--c")==0) {
			x++;
			p.c=atof(argv[x]);
		} else if (p_switch.compare("--read")==0) {
			x++;
			p.read=atoi(argv[x]);
		} else if (p_switch.compare("--unobs")==0) {
			x++;
			p.unobs=atoi(argv[x]);
		} else if (p_switch.compare("--readtimes")==0) {
			x++;
			p.timefile=atoi(argv[x]);
		}
		else {
			cout << "Incorrect usage\n ";
			exit(1);
		}
		p_switch.clear();
		x++;
	}
}


void GetOptionsQuad (run_params& p, int argc, const char **argv) {
	string p_switch;
	p.read=0;
	p.unobs=0;
	p.seed=0;
	p.timefile=0;
	p.mu=0.333333333333e-5;
	p.c=60.1155;
	p.seed=atoi(argv[1]);
	p.data1=argv[2];
	p.data2=argv[3];
	p.data3=argv[4];
	p.data4=argv[5];
	int x=6;
	while (x < argc && (argv[x][0]=='-')) {
		p_switch=argv[x];
		if (p_switch.compare("--mu")==0) {
			x++;
			p.mu=atof(argv[x]);
		} else if (p_switch.compare("--c")==0) {
			x++;
			p.c=atof(argv[x]);
		} else if (p_switch.compare("--read")==0) {
			x++;
			p.read=atoi(argv[x]);
		} else if (p_switch.compare("--unobs")==0) {
			x++;
			p.unobs=atoi(argv[x]);
		} else if (p_switch.compare("--readtimes")==0) {
			x++;
			p.timefile=atoi(argv[x]);
		}
		else {
			cout << "Incorrect usage\n ";
			exit(1);
		}
		p_switch.clear();
		x++;
	}
}

void GetOptionsSix (run_params& p, int argc, const char **argv) {
	string p_switch;
	p.read=0;
	p.unobs=0;
	p.seed=0;
	p.timefile=0;
	p.mu=0.333333333333e-5;
	p.c=60.1155;
	p.seed=atoi(argv[1]);
	p.data1=argv[2];
	p.data2=argv[3];
	p.data3=argv[4];
	p.data4=argv[5];
	p.data5=argv[6];
	p.data6=argv[7];
	int x=8;
	while (x < argc && (argv[x][0]=='-')) {
		p_switch=argv[x];
		if (p_switch.compare("--mu")==0) {
			x++;
			p.mu=atof(argv[x]);
		} else if (p_switch.compare("--c")==0) {
			x++;
			p.c=atof(argv[x]);
		} else if (p_switch.compare("--read")==0) {
			x++;
			p.read=atoi(argv[x]);
		} else if (p_switch.compare("--unobs")==0) {
			x++;
			p.unobs=atoi(argv[x]);
		} else if (p_switch.compare("--readtimes")==0) {
			x++;
			p.timefile=atoi(argv[x]);
		}
		else {
			cout << "Incorrect usage\n ";
			exit(1);
		}
		p_switch.clear();
		x++;
	}
}

void GetOptionsEight (run_params& p, int argc, const char **argv) {
	string p_switch;
	p.read=0;
	p.unobs=0;
	p.seed=0;
	p.timefile=0;
	p.mu=0.333333333333e-5;
	p.c=60.1155;
	p.seed=atoi(argv[1]);
	p.data1=argv[2];
	p.data2=argv[3];
	p.data3=argv[4];
	p.data4=argv[5];
	p.data5=argv[6];
	p.data6=argv[7];
	p.data7=argv[8];
	p.data8=argv[9];
	int x=10;
	while (x < argc && (argv[x][0]=='-')) {
		p_switch=argv[x];
		if (p_switch.compare("--mu")==0) {
			x++;
			p.mu=atof(argv[x]);
		} else if (p_switch.compare("--c")==0) {
			x++;
			p.c=atof(argv[x]);
		} else if (p_switch.compare("--read")==0) {
			x++;
			p.read=atoi(argv[x]);
		} else if (p_switch.compare("--unobs")==0) {
			x++;
			p.unobs=atoi(argv[x]);
		} else if (p_switch.compare("--readtimes")==0) {
			x++;
			p.timefile=atoi(argv[x]);
		}
		else {
			cout << "Incorrect usage\n ";
			exit(1);
		}
		p_switch.clear();
		x++;
	}
}

void ObsTimes (run_params p, vector<int>& times) {
	if (p.timefile==1) {
		ifstream t_file;
		t_file.open("../Times.in");
		int s;
		for (int t=0;t<1000;t++) {
			if (!(t_file >> s)) break;
			times.push_back(s);
			cout << s << "\n";
		}
	} else {
		times.push_back(0);
		times.push_back(1);
		times.push_back(3);
		times.push_back(5);
	}

}

void InputDataHap (vector<int> times, ifstream& in_file, vector<rec>& dat) {  //Should change to flag up non-zero frequency observations.  Wish to convert these to a list of initial frequencies that can be set to zero and excluded as parameters.  Have to figure out map to marginals...
	cout << "Reading allele frequency data...\n";
	string s;
	do {
		rec r;
		if (!(in_file >> s)) break;
		cout << s << "\n";
		r.st=s;
		int k;
		for (int i=0;i<times.size();i++) {
			in_file >> k;
			r.obs.push_back(k);
		}
		for (int i=0;i<s.size();i++) {
			r.seq.push_back(s[i]);
		}
		dat.push_back(r);
	} while (1==1);
}

void FindLogFact(vector<double>& fact_store,int N){
	double logN=0;
	fact_store.push_back(0);
	for (int i=1;i<=N;i++) {
		logN=logN+log(i);
		fact_store.push_back(logN);
		//cout << "fact_store "<<i<<" "<<gsl_vector_get(fact_store,i)<<"\n";
	}
}

void GetObsStats (vector<int> times, vector<rec> dat1, vector<int>& N, vector< vector<int> >& obs) {
	for (int j=0;j<times.size();j++) {
		vector<int> o;
		for (int i=0;i<dat1.size();i++) {
			o.push_back(dat1[i].obs[j]);
		}
		obs.push_back(o);
	}
	for (int i=0;i<obs.size();i++) {
		int n=0;
		for (int j=0;j<obs[i].size();j++) {
			n=n+obs[i][j];
		}
		N.push_back(n);
		cout << n << "\n";
	}
}

void GetConsensusBase (vector<int>& consensus, vector<int>& bases, vector<rec>& dat1) {
	vector<int> count;
	for (int i=0;i<dat1[0].seq.size();i++) {
		count.clear();
		for (int k=0;k<4;k++) {
			count.push_back(0);
		}
		for (int j=0;j<dat1.size();j++) {
			if (dat1[j].seq[i]=='A') {
				count[0]=count[0]+dat1[j].obs[0];
			}
			if (dat1[j].seq[i]=='C') {
				count[1]=count[1]+dat1[j].obs[0];
			}
			if (dat1[j].seq[i]=='G') {
				count[2]=count[2]+dat1[j].obs[0];
			}
			if (dat1[j].seq[i]=='T') {
				count[3]=count[3]+dat1[j].obs[0];
			}
		}
		if ((count[0]>=count[1])&&(count[0]>=count[2])&&(count[0]>=count[3])) {
			consensus.push_back(0);
		} else if ((count[1]>=count[0])&&(count[1]>=count[2])&&(count[1]>=count[3])) {
			consensus.push_back(1);
		} else if ((count[2]>=count[0])&&(count[2]>=count[1])&&(count[2]>=count[3])) {
			consensus.push_back(2);
		} else {
			consensus.push_back(3);
		}
//		cout << consensus[i] << "\n";
	}
	
//	cout << "Bases\n";
	
	//Get minority bases
	for (int i=0;i<dat1[0].seq.size();i++) {
		count.clear();
		for (int k=0;k<4;k++) {
			count.push_back(0);
		}
		for (int j=0;j<dat1.size();j++) {
			//		cout <<dat1[j].seq[i] << " " << consensus[i] << "\n";
			if (dat1[j].seq[i]=='A'&&consensus[i]!=0) {
				for (int k=0;k<dat1[k].obs.size();k++) {
					count[0]=count[0]+dat1[j].obs[k];
				}
			}
			if (dat1[j].seq[i]=='C'&&consensus[i]!=1) {
				for (int k=0;k<dat1[k].obs.size();k++) {
					count[1]=count[1]+dat1[j].obs[k];
				}
			}
			if (dat1[j].seq[i]=='G'&&consensus[i]!=2) {
				for (int k=0;k<dat1[k].obs.size();k++) {
					count[2]=count[2]+dat1[j].obs[k];
				}
			}
			if (dat1[j].seq[i]=='T'&&consensus[i]!=3) {
				for (int k=0;k<dat1[k].obs.size();k++) {
					count[3]=count[3]+dat1[j].obs[k];
				}
			}
		}
		//	cout << count[0] << " " << count[1] << " " << count[2] << " " << count[3] << "\n";
		if ((count[0]==0)&&(count[1]==0)&&(count[2]==0)&&(count[3]==0)) {
			//		cout << "No non-consensus allele\n";
			bases.push_back(consensus[i]);
		} else {
			if ((count[0]>=count[1])&&(count[0]>=count[2])&&(count[0]>=count[3])) {
				bases.push_back(0);
			} else if ((count[1]>=count[0])&&(count[1]>=count[2])&&(count[1]>=count[3])) {
				bases.push_back(1);
			} else if ((count[2]>=count[0])&&(count[2]>=count[1])&&(count[2]>=count[3])) {
				bases.push_back(2);
			} else {
				bases.push_back(3);
			}
		}
//		cout << bases[i] << "\n";
	}
	
	//Edit sequences, bases and consensus if there is no variation at a given site.
	vector<int> remove;
	for (int i=0;i<bases.size();i++) {
		if (bases[i]==consensus[i]) {
			cout << "No variation at position " << i << "\n";
			remove.push_back(i);
		}
	}
	
	 /*cout << "Consensus\n";
	 for (int i=0;i<consensus.size();i++) {
	 cout << consensus[i] << " ";
	 }
	 cout << "\n";
	 cout << "Bases\n";
	 for (int i=0;i<bases.size();i++) {
	 cout << bases[i] << " ";
	 }
	 cout << "\n";*/
	
	
	if (remove.size()>0) {
		vector<int> newcons;
		vector<int> newbases;
		for (int i=0;i<consensus.size();i++) {
			int keep=1;
			for (int j=0;j<remove.size();j++) {
				if (remove[j]==i) {
					keep=0;
					//				cout << "Remove " << i << "\n";
				}
			}
			if (keep==1) {
				//			cout << "Push " << i << "\n";
				newcons.push_back(consensus[i]);
				newbases.push_back(bases[i]);
				for (int j=0;j<dat1.size();j++) {
					dat1[j].newseq.push_back(dat1[j].seq[i]);
				}
			}
		}
		consensus=newcons;
		bases=newbases;
		for (int j=0;j<dat1.size();j++) {
			dat1[j].seq=dat1[j].newseq;
		}
		
	}
}

void GetConsensusBaseAllHaps (vector<int>& consensus, vector<int>& bases, vector<rec>& dat1) {
	vector<int> count;
	for (int i=0;i<dat1[0].seq.size();i++) {
		count.clear();
		for (int k=0;k<4;k++) {
			count.push_back(0);
		}
		for (int j=0;j<dat1.size();j++) {
			if (dat1[j].seq[i]=='A') {
				count[0]=count[0]+dat1[j].obs[0];
			}
			if (dat1[j].seq[i]=='C') {
				count[1]=count[1]+dat1[j].obs[0];
			}
			if (dat1[j].seq[i]=='G') {
				count[2]=count[2]+dat1[j].obs[0];
			}
			if (dat1[j].seq[i]=='T') {
				count[3]=count[3]+dat1[j].obs[0];
			}
		}
		if ((count[0]>=count[1])&&(count[0]>=count[2])&&(count[0]>=count[3])) {
			consensus.push_back(0);
		} else if ((count[1]>=count[0])&&(count[1]>=count[2])&&(count[1]>=count[3])) {
			consensus.push_back(1);
		} else if ((count[2]>=count[0])&&(count[2]>=count[1])&&(count[2]>=count[3])) {
			consensus.push_back(2);
		} else {
			consensus.push_back(3);
		}
		//		cout << consensus[i] << "\n";
	}
	
	//	cout << "Bases\n";
	
	//Get minority bases
	for (int i=0;i<dat1[0].seq.size();i++) {
		count.clear();
		for (int k=0;k<4;k++) {
			count.push_back(0);
		}
		for (int j=0;j<dat1.size();j++) {
			//		cout <<dat1[j].seq[i] << " " << consensus[i] << "\n";
			if (dat1[j].seq[i]=='A'&&consensus[i]!=0) {
				for (int k=0;k<dat1[k].obs.size();k++) {
					count[0]=count[0]+dat1[j].obs[k];
				}
			}
			if (dat1[j].seq[i]=='C'&&consensus[i]!=1) {
				for (int k=0;k<dat1[k].obs.size();k++) {
					count[1]=count[1]+dat1[j].obs[k];
				}
			}
			if (dat1[j].seq[i]=='G'&&consensus[i]!=2) {
				for (int k=0;k<dat1[k].obs.size();k++) {
					count[2]=count[2]+dat1[j].obs[k];
				}
			}
			if (dat1[j].seq[i]=='T'&&consensus[i]!=3) {
				for (int k=0;k<dat1[k].obs.size();k++) {
					count[3]=count[3]+dat1[j].obs[k];
				}
			}
		}
		//	cout << count[0] << " " << count[1] << " " << count[2] << " " << count[3] << "\n";
		if ((count[0]==0)&&(count[1]==0)&&(count[2]==0)&&(count[3]==0)) {
			//		cout << "No non-consensus allele\n";
			bases.push_back(consensus[i]);
		} else {
			if ((count[0]>=count[1])&&(count[0]>=count[2])&&(count[0]>=count[3])) {
				bases.push_back(0);
			} else if ((count[1]>=count[0])&&(count[1]>=count[2])&&(count[1]>=count[3])) {
				bases.push_back(1);
			} else if ((count[2]>=count[0])&&(count[2]>=count[1])&&(count[2]>=count[3])) {
				bases.push_back(2);
			} else {
				bases.push_back(3);
			}
		}
		//		cout << bases[i] << "\n";
	}
	
	//Keep all sequences, even where no variation is seen.  Variation should be picked up in the non-observed sites

}

void MakeMutationMatrix (int dim, double mu, vector< vector<double> >& mut, vector<rec> dat1) {
	//Check distances between haplotypes
	vector< vector<int> > m;
	for (int i=0;i<dat1.size();i++) {
		vector<int> row;
		for (int j=0;j<dat1.size();j++) {
			int dist=0;
			for (int k=0;k<dat1[0].seq.size();k++) {
				if (dat1[i].seq[k]!=dat1[j].seq[k]) {
					dist++;
				}
			}
			if (dist==1) {
				if (i<j) {
					cout << i << " " << j << " " << dist << "\n";
				}
				row.push_back(1);
			} else {
				row.push_back(0);
			}
		}
		m.push_back(row);
	}
	
	vector<int> tot;
	for (int i=0;i<m.size();i++) {
		int t=0;
		for (int j=0;j<m[i].size();j++) {
			//                      cout << m[i][j] << " ";
			if (m[i][j]==1) {
				t++;
			}
		}
		cout << t << " " << (3*dim)-t << " ";
		tot.push_back(t);
		cout << "\n";
	}
	
	for (int i=0;i<m.size();i++) {
		vector<double> mrow;
		for (int j=0;j<m[i].size();j++) {
			if (m[i][j]==1) {
				mrow.push_back(mu);
			} else if (i==j) {
				mrow.push_back(pow(1-mu,3*dim));  //Incorporate mutations to other haplotypes
			} else {
				mrow.push_back(0);
			}
		}
		mrow.push_back(0);
		mut.push_back(mrow);
	}
	vector<double> mrow;
	for (int i=0;i<m[0].size();i++) {
		mrow.push_back(mu*((3*dim)-tot[i]));
	}
	mrow.push_back(1);
	mut.push_back(mrow);
	
	cout << "Mutation matrix\n";
	for (int i=0;i<mut.size();i++) {
		for (int j=0;j<mut[i].size();j++) {
			cout << mut[i][j] << " ";
		}
		cout << "\n";
	}
	
}

void SquareMutationMatrix (vector<vector<double> > &m ) {
//	cout << "M size " << m.size() << "\n";
	vector<vector<double> > m2;
	for (unsigned int i=0;i<m.size();i++) {
		vector<double> m2i;
		m2i.assign(m.size(),0);
		m2.push_back(m2i);
	}
//	cout << "Now multiplying...\n";
	for (unsigned int i=0;i<m.size();i++) {
		for (unsigned int j=0;j<m[i].size();j++) {
			for (unsigned int k=0;k<m[i].size();k++) {
				m2[i][j]=m2[i][j]+(m[i][k]*m[k][j]);
			}
		}
	}
	m=m2;
}



void MakeMutationMatrixLinear (int dim, double mu, vector< vector<double> >& mut, vector<seldiv> sd, vector<rec> dat1) {
	//Calculate distances between haplotypes
	vector< vector<int> > m;
	for (int i=0;i<dat1.size();i++) {
		vector<int> row;
		for (int j=0;j<dat1.size();j++) {
			int dist=0;
			for (int k=0;k<dat1[0].seq.size();k++) {
				if (dat1[i].seq[k]!=dat1[j].seq[k]) {
					dist++;
				}
			}
			if (dist==1) {
				if (i<j) {
					cout << i << " " << j << " " << dist << "\n";
				}
				row.push_back(1);
			} else {
				row.push_back(0);
			}
		}
		m.push_back(row);
	}
	
	vector<int> tot;
	for (int i=0;i<m.size();i++) {
		int t=0;
		for (int j=0;j<m[i].size();j++) {
//			cout << m[i][j] << " ";
			if (m[i][j]==1) {
				t++;
			}
		}
		cout << t << " " << (3*dim)-t << " ";
		tot.push_back(t);
		cout << "\n";
	}
	
	
	//Matrix T: to implicit states
	cout << "Matrix T transpose\n";
	vector< vector<double> > t_implic;
	for (int i=0;i<dat1.size();i++) {
		vector<double> mir;
		//Construct single mutants
		vector< vector<char> > smut;
		vector<char> orig=dat1[i].seq;
		vector<char> mut=orig;
	//	cout << "Original ";
	//	for (int k=0;k<orig.size();k++) {
	//		cout << orig[k] << " ";
	//	}
	//	cout << "\n";

		
		for (int j=0;j<dat1[i].seq.size();j++) {
			vector<char> mut=orig;
			if (mut[j]=='A') {
				mut[j]='C';
				smut.push_back(mut);
				mut[j]='G';
				smut.push_back(mut);
				mut[j]='T';
				smut.push_back(mut);
			} else if (mut[j]=='C') {
				mut[j]='A';
				smut.push_back(mut);
				mut[j]='G';
				smut.push_back(mut);
				mut[j]='T';
				smut.push_back(mut);
			} else if (mut[j]=='G') {
				mut[j]='C';
				smut.push_back(mut);
				mut[j]='A';
				smut.push_back(mut);
				mut[j]='T';
				smut.push_back(mut);
			} else if (mut[j]=='T') {
				mut[j]='C';
				smut.push_back(mut);
				mut[j]='G';
				smut.push_back(mut);
				mut[j]='A';
				smut.push_back(mut);
			}
		}
	//	cout << "S mutation size " << smut.size() << "\n";;
		for (int k=0;k<sd.size();k++) {
			sd[k].count=0;
		}
		for (int j=0;j<smut.size();j++) {
			//Is this mutation in the explicit set
			int exp=0;
			for (int k=0;k<dat1.size();k++) {
				int d=0;
				for (int l=0;l<smut[j].size();l++) {
					if (dat1[k].seq[l]!=smut[j][l]) {
						d++;
					}
				}
				if (d==0) {
				//	cout << "Match ";
				//	for (int l=0;l<smut[j].size();l++) {
				//		cout << smut[j][l] << " ";
				//	}
				//	cout << " with ";
				//	for (int l=0;l<dat1[k].seq.size();l++) {
				//		cout << dat1[k].seq[l] << " ";
				//	}
				//	cout << "\n";
					exp=1;
				}
			}
			if (exp==0) { //Not in the set of observed haplotypes
				for (int k=0;k<sd.size();k++) {
					int d=0;
					for (int l=0;l<smut[j].size();l++) {
						if ((sd[k].seq[l]=='A'&&smut[j][l]!='A')||(sd[k].seq[l]=='1'&&smut[j][l]=='A')) {
							d++;
						}
						if ((sd[k].seq[l]=='C'&&smut[j][l]!='C')||(sd[k].seq[l]=='2'&&smut[j][l]=='C')) {
							d++;
						}
						if ((sd[k].seq[l]=='G'&&smut[j][l]!='G')||(sd[k].seq[l]=='3'&&smut[j][l]=='G')) {
							d++;
						}
						if ((sd[k].seq[l]=='T'&&smut[j][l]!='T')||(sd[k].seq[l]=='4'&&smut[j][l]=='T')) {
							d++;
						}
					}
					//cout << "i " << i << " j " << j << " k " << k << " d " << d << "\n";

					if (d==0) {
						sd[k].count++;
		//				cout << "i " << i << " j " << j << " k " << k << " d " << d << " ";
		//				cout << "Match ";
		//				for (int l=0;l<smut[j].size();l++) {
		//					cout << smut[j][l] << " ";
		//				}
		//				cout << " with ";
		//				for (int l=0;l<sd[k].seq.size();l++) {
		//					cout << sd[k].seq[l] << " ";
		//				}
		//				cout << "\n";
					}
				}
			}
		}
		for (int k=0;k<sd.size();k++) {
		//	cout << k << " ";
		//	for (int l=0;l<sd[k].seq.size();l++) {
		//		cout << sd[k].seq[l];
		//	}
			double z=sd[k].count*mu;
			mir.push_back(z);
			cout << z << " ";
		}
		//for (int k=0;k<mir.size();k++) {
		//	cout << mir[k] << "  ";
		//}
		cout << "\n";
		t_implic.push_back(mir);
	}
	
//	cout << "Matrix T\n";
//	for (int i=0;i<t_implic.size();i++) {
//		for (int j=0;j<t_implic[i].size();j++) {
//			cout << t_implic[i][j] << " ";
//		}
//		cout << "\n";
//	}

	
	//Matrix I: between implicit states
	vector< vector<double> > i_implic;
	for (int i=0;i<sd.size();i++) {
		vector<double> irow;
		double tot=0;
		for (int j=0;j<sd.size();j++) {
			int d=0;
			for (int k=0;k<sd[i].seq.size();k++) {
				if (sd[i].seq[k]!=sd[j].seq[k]) {
					d++;
				}
			}
			if (d==1) {
				for (int k=0;k<sd[i].seq.size();k++) {
					if (sd[i].seq[k]!=sd[j].seq[k]) {
						if (sd[i].seq[k]=='1'||sd[i].seq[k]=='2'||sd[i].seq[k]=='3'||sd[i].seq[k]=='4') {
							irow.push_back(3*mu);
							tot=tot+(3*mu);
							//cout << "i " << i << " j " << j << " r " << 3 << "\n";
						} else {
							irow.push_back(mu);
							tot=tot+mu;
							//cout << "i " << i << " j " << j << " r " << 1 << "\n";

						}
					}
				}
			} else if (d>1) {
				irow.push_back(0);
			} else if (d==0) {
				irow.push_back(1);
			}
			//cout << i << " " << irow[i];
		}
		for (int j=0;j<irow.size();j++) {
			if (irow[j]==1){
				irow[j]=1-tot;
			}
		}
		i_implic.push_back(irow);
	}
	cout << "Matrix I\n";
	for (int i=0;i<i_implic.size();i++) {
		for (int j=0;j<i_implic[i].size();j++) {
			cout << i_implic[i][j] << " ";
		}
		cout << "\n";
	}
	
	for (int i=0;i<m.size();i++) {
		vector<double> mrow;
		for (int j=0;j<m[i].size();j++) {
			if (m[i][j]==1) {
				mrow.push_back(mu);
			} else if (i==j) {
				mrow.push_back(pow(1-mu,3*dim));  //Incorporate mutations to other haplotypes
			} else {
				mrow.push_back(0);
			}
		}
//		mrow.push_back(0);
		mut.push_back(mrow);
	}
//	vector<double> mrow;
//	for (int i=0;i<m[0].size();i++) {
//		mrow.push_back(mu*((3*dim)-tot[i]));
//	}
//	mrow.push_back(1);
//	mut.push_back(mrow);
	
	cout << "Matrix E\n";
	
	//vector< vector<int> > m
	
	for (int i=0;i<mut.size();i++) {
		for (int j=0;j<mut[i].size();j++) {
			cout << mut[i][j] << " ";
		}
		cout << "\n";
	}
	
	cout << "Matrix A\n";
	vector< vector<double> > a;
	for (int i=0;i<mut.size();i++) {
		vector<double> row;
		for (int j=0;j<mut[i].size();j++) {
			row.push_back(mut[i][j]);
		}
		for (int j=0;j<t_implic[i].size();j++) {
			row.push_back(0);
		}
		a.push_back(row);
	}
	
//	for (int i=0;i<a.size();i++) {
//		for (int j=0;j<a[i].size();j++) {
//			cout << a[i][j] << " ";
//		}
//		cout << "\n";
//	}

	for (int i=0;i<t_implic[0].size();i++) {
		vector<double> row;
		for (int j=0;j<t_implic.size();j++) {
			row.push_back(t_implic[j][i]);
		}
		for (int j=0;j<i_implic[i].size();j++) {
			row.push_back(i_implic[i][j]);
		}
		a.push_back(row);
	}

	
	mut=a;

	for (int i=0;i<mut.size();i++) {
		for (int j=0;j<mut[i].size();j++) {
			cout << mut[i][j] << " ";
		}
		cout << "\n";
	}
}

void ReadSelectionModel(int dim, vector<char>& sel_model, vector<int>& tds, vector<epi>& epistat2, vector<epi>& epistat3, vector<epi>& epistat4, vector<epi>& epistat5) {
	ifstream sel_file;
	sel_file.open("Selection.in");
	for (int i=0;i<dim;i++) {
		char k='0'; //Default of no selection
		sel_file >> k;
		sel_model.push_back(k);
		cout << k << " ";
	}
	cout << "\n";
	sel_file.close();
	ifstream tds_file;
	tds_file.open("TDS.in");
	for (int i=0;i<dim;i++) {
		int k=0; //Default of no selection
		tds_file >> k;
		tds.push_back(k);
		cout << k << " ";
	}
	cout << "\n";
	tds_file.close();
	
	ifstream epi_file;
	epi_file.open("TwoWay.in");
	for (int i=0;i<1000;i++) {
		epi e;
		int l1=0;
		int l2=0;
		if (!(epi_file >> l1)) break;
		if (!(epi_file >> l2)) break;
		l1--;
		l2--;
		e.loc.push_back(l1);
		e.loc.push_back(l2);
		e.x=0;
		epistat2.push_back(e);
	}
	epi_file.close();
	
	epi_file.open("ThreeWay.in");
	for (int i=0;i<1000;i++) {
		epi e;
		int l1=0;
		int l2=0;
		int l3=0;
		if (!(epi_file >> l1)) break;
		if (!(epi_file >> l2)) break;
		if (!(epi_file >> l3)) break;
		l1--;
		l2--;
		l3--;
		e.loc.push_back(l1);
		e.loc.push_back(l2);
		e.loc.push_back(l3);
		e.x=0;
		epistat3.push_back(e);
	}
	epi_file.close();

	epi_file.open("FourWay.in");
	for (int i=0;i<1000;i++) {
		epi e;
		int l1=0;
		int l2=0;
		int l3=0;
		int l4=0;
		if (!(epi_file >> l1)) break;
		if (!(epi_file >> l2)) break;
		if (!(epi_file >> l3)) break;
		if (!(epi_file >> l4)) break;
		l1--;
		l2--;
		l3--;
		l4--;
		e.loc.push_back(l1);
		e.loc.push_back(l2);
		e.loc.push_back(l3);
		e.loc.push_back(l4);
		e.x=0;
		epistat4.push_back(e);
	}
	epi_file.close();

	epi_file.open("FiveWay.in");
	for (int i=0;i<1000;i++) {
		epi e;
		int l1=0;
		int l2=0;
		int l3=0;
		int l4=0;
		int l5=0;
		if (!(epi_file >> l1)) break;
		if (!(epi_file >> l2)) break;
		if (!(epi_file >> l3)) break;
		if (!(epi_file >> l4)) break;
		if (!(epi_file >> l5)) break;
		l1--;
		l2--;
		l3--;
		l4--;
		l5--;
		e.loc.push_back(l1);
		e.loc.push_back(l2);
		e.loc.push_back(l3);
		e.loc.push_back(l4);
		e.loc.push_back(l5);
		e.x=0;
		epistat5.push_back(e);
	}
	epi_file.close();

}

void SetInitialFreqs (vector<double>& init_freqs, vector< vector<double> > mut, vector<rec> dat1, gsl_rng *rgen) {
	for (int i=0;i<mut.size();i++) {
		double x=1e-20;
		if (i<dat1.size()) {
			x=gsl_rng_uniform(rgen);
		}
		init_freqs.push_back(x);
	//	cout << x << " ";
	}
	NormaliseFreqs(init_freqs,mut);
	//cout << "\n";
}

void SetInitialFreqsLinear (vector<double>& init_freqs, vector< vector<double> > mut, vector<rec> dat1, gsl_rng *rgen) {
	cout << "Initial freqs\n";
	for (int i=0;i<mut.size();i++) {
		double x=1e-20;
		if (i<dat1.size()) {
			x=gsl_rng_uniform(rgen);
		}
		init_freqs.push_back(x);
	//	cout << x << " ";
	}
	NormaliseFreqs(init_freqs,mut);
	//cout << "\n";
}


void SetInitialFreqs2 (vector<double>& init_freqs, vector< vector<double> > mut, vector<int> obs, vector<int> N) {
	for (int i=0;i<mut.size();i++) {
		double x=0;
		if (N[i]>0) {
			x=(obs[i]+0.)/(N[i]+0.);
		} else {
			x=1e-5;
		}
		if (x<1e-5) {
			x=1e-5;
		}
	//	cout << x << "\n";
		init_freqs.push_back(x);
	}
	NormaliseFreqs(init_freqs,mut);
}


void NormaliseFreqs (vector<double>& init_freqs,vector< vector<double> > mut) {
	double tot=0;
	for (int i=0;i<mut.size();i++) {
		tot=tot+init_freqs[i];
	}
	for (int i=0;i<mut.size();i++) {
		init_freqs[i]=init_freqs[i]/tot;
	//	cout << init_freqs[i] << " ";
	}
	//cout << "\n";
}

void SetInitialSelection (int dim, vector<int> tds, vector<char> sel_model, vector<int> times, vector< vector<double> >& sigs, gsl_rng *rgen) {
	for (int i=0;i<dim;i++) {
		vector<double> sig;
		if (sel_model[i]!='0') {
			if (tds[i]==0) {
				double s=gsl_rng_uniform(rgen)-0.5;
				sig.push_back(s);
			} else {
				for (int t=0;t<times.size()-1;t++) {
					double s=gsl_rng_uniform(rgen)-0.5;
					sig.push_back(s);
				}
			}
		}
		sigs.push_back(sig);
	}
	for (int i=0;i<sigs.size();i++) {
		cout << i << " ";
		for (int j=0;j<sigs[i].size();j++) {
			cout << sigs[i][j] << " ";
		}
		cout << "\n";
	}
}

void SetInitialEpistasis (vector<epi>& epistat2, vector<epi>& epistat3, vector<epi>& epistat4, vector<epi>& epistat5, gsl_rng *rgen) {
	for (int i=0;i<epistat2.size();i++) {
		double s=gsl_rng_uniform(rgen)-0.5;
		epistat2[i].x=s;
		cout << "Epi2 " << epistat2[i].loc[0] << " " << epistat2[i].loc[1] << " " << epistat2[i].x << "\n";
	}
	for (int i=0;i<epistat3.size();i++) {
		double s=gsl_rng_uniform(rgen)-0.5;
		epistat3[i].x=s;
		cout << "Epi3 " << epistat3[i].loc[0] << " " << epistat3[i].loc[1] << " " << epistat3[i].loc[2] << " " << epistat3[i].x << "\n";
	}
	for (int i=0;i<epistat4.size();i++) {
		double s=gsl_rng_uniform(rgen)-0.5;
		epistat4[i].x=s;
		cout << "Epi4 " << epistat4[i].loc[0] << " " << epistat4[i].loc[1] << " " << epistat4[i].loc[2] << " " << epistat4[i].loc[3] << " " << epistat4[i].x << "\n";
	}
	for (int i=0;i<epistat5.size();i++) {
		double s=gsl_rng_uniform(rgen)-0.5;
		epistat5[i].x=s;
		cout << "Epi5 " << epistat5[i].loc[0] << " " << epistat5[i].loc[1] << " " << epistat5[i].loc[2] << " " << epistat5[i].loc[3] << " " << epistat5[i].loc[4] << " " << epistat5[i].x << "\n";
	}
}

void CheckEpistasis (int& go, int dim, vector<char> sel_model, vector<rec> dat1, vector<epi>& epistat2, vector<epi>& epistat3, vector<epi>& epistat4, vector<epi>& epistat5) {
	cout << "Check epistasis\n";
	int sugg=0;
	int seen=0;
	for (int i=0;i<epistat2.size();i++) {
		sugg++;
		int f2=0;
		for (int j=0;j<dim;j++) {
			for (int k=0;k<dim;k++) {
				if (epistat2[i].loc[0]==j&&epistat2[i].loc[1]==k) { //Have an epistatic effect of epistat2[i].x on double mutants j,k
					//cout << "j= " << j << " k= " << k << " " << sel_model[j] << " " << sel_model[k] << "\n";
					for (int s=0;s<dat1.size();s++) {
						if (dat1[s].seq[j]==sel_model[j]&&dat1[s].seq[k]==sel_model[k]) {
							int d=0;
							for (int t=0;t<dat1[s].obs.size();t++) {
								if (dat1[s].obs[t]>0) {
									d=1;
								}
							}
							if (f2==0&&d==1) {
								cout << "Have two-way interaction " << i << " " << epistat2[i].loc[0] << " " << epistat2[i].loc[1] << "\n";
								f2=1;
								seen++;
							}
						}
					}
				}
			}
		}
	}
	
	for (int i=0;i<epistat3.size();i++) {
		sugg++;
		int f3=0;
		for (int k1=0;k1<dim;k1++) {
			for (int k2=0;k2<dim;k2++) {
				for (int k3=0;k3<dim;k3++) {
					if (epistat3[i].loc[0]==k1&&epistat3[i].loc[1]==k2&&epistat3[i].loc[2]==k3) { //Have an epistatic effect of epistat3[i].x on triple mutants k1,k2,k3
						for (int s=0;s<dat1.size();s++) {
							if (dat1[s].seq[k1]==sel_model[k1]&&dat1[s].seq[k2]==sel_model[k2]&&dat1[s].seq[k3]==sel_model[k3]) {
								int d=0;
								for (int t=0;t<dat1[s].obs.size();t++) {
									if (dat1[s].obs[t]>0) {
										d=1;
									}
								}
								if (f3==0&&d==1) {
									cout << "Have three-way interaction " << i << " " << epistat3[i].loc[0] << " " << epistat3[i].loc[1] << " " << epistat3[i].loc[2] << "\n";
									f3=1;
									seen++;
								}
							}
						}
					}
				}
			}
		}
	}
	
	for (int i=0;i<epistat4.size();i++) {
		sugg++;
		int f4=0;
		for (int k1=0;k1<dim;k1++) {
			for (int k2=0;k2<dim;k2++) {
				for (int k3=0;k3<dim;k3++) {
					for (int k4=0;k4<dim;k4++) {
						if (epistat4[i].loc[0]==k1&&epistat4[i].loc[1]==k2&&epistat4[i].loc[2]==k3&&epistat4[i].loc[3]==k4) { //Have an epistatic effect of epistat2[i].x on double mutants j,k
							for (int s=0;s<dat1.size();s++) {
								if (dat1[s].seq[k1]==sel_model[k1]&&dat1[s].seq[k2]==sel_model[k2]&&dat1[s].seq[k3]==sel_model[k3]&&dat1[s].seq[k4]==sel_model[k4]) {
									int d=0;
									for (int t=0;t<dat1[s].obs.size();t++) {
										if (dat1[s].obs[t]>0) {
											d=1;
										}
									}
									if (f4==0&&d==1) {
										cout << "Have four-way interaction " << i << " " << epistat4[i].loc[0] << " " << epistat4[i].loc[1] << " " << epistat4[i].loc[2] << " " << epistat4[i].loc[3] << "\n";
										f4=1;
										seen++;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	
	for (int i=0;i<epistat5.size();i++) {
		sugg++;
		int f5=0;
		for (int k1=0;k1<dim;k1++) {
			for (int k2=0;k2<dim;k2++) {
				for (int k3=0;k3<dim;k3++) {
					for (int k4=0;k4<dim;k4++) {
						for (int k5=0;k5<dim;k5++) {
							if (epistat5[i].loc[0]==k1&&epistat5[i].loc[1]==k2&&epistat5[i].loc[2]==k3&&epistat5[i].loc[3]==k4&&epistat5[i].loc[4]==k5) { //Have an epistatic effect of epistat2[i].x on double mutants j,k
								for (int s=0;s<dat1.size();s++) {
									if (dat1[s].seq[k1]==sel_model[k1]&&dat1[s].seq[k2]==sel_model[k2]&&dat1[s].seq[k3]==sel_model[k3]&&dat1[s].seq[k4]==sel_model[k4]&&dat1[s].seq[k5]==sel_model[k5]) {
										int d=0;
										for (int t=0;t<dat1[s].obs.size();t++) {
											if (dat1[s].obs[t]>0) {
												d=1;
											}
										}
										if (f5==0&&d==1) {
											cout << "Have five-way interaction " << i << " " << epistat5[i].loc[0] << " " << epistat5[i].loc[1] << " " << epistat5[i].loc[2] << " " << epistat5[i].loc[3] << " " << epistat5[i].loc[4] << "\n";
											f5=1;
											seen++;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	if (seen<sugg) {
		cout << "Not all epistatic interactions were observed\n";
		go=0;
	}
}


void CheckDatEpi2 (int i, int j, int k, int& f2, int& seen, vector<char> sel_model, vector<epi> epistat2, vector<rec> dat) {
	for (int s=0;s<dat.size();s++) {
		if (dat[s].seq[j]==sel_model[j]&&dat[s].seq[k]==sel_model[k]) {
			int d=0;
			for (int t=0;t<dat[s].obs.size();t++) {
				if (dat[s].obs[t]>0) {
					d=1;
				}
			}
			if (f2==0&&d==1) {
				cout << "Have two-way interaction " << i << " " << epistat2[i].loc[0] << " " << epistat2[i].loc[1] << "\n";
				f2=1;
				seen++;
			}
		}
	}
}

void CheckDatEpi3 (int i, int k1, int k2, int k3, int& f3, int& seen, vector<char> sel_model, vector<epi> epistat3, vector<rec> dat) {
	for (int s=0;s<dat.size();s++) {
		if (dat[s].seq[k1]==sel_model[k1]&&dat[s].seq[k2]==sel_model[k2]&&dat[s].seq[k3]==sel_model[k3]) {
			int d=0;
			for (int t=0;t<dat[s].obs.size();t++) {
				if (dat[s].obs[t]>0) {
					d=1;
				}
			}
			if (f3==0&&d==1) {
				cout << "Have three-way interaction " << i << " " << epistat3[i].loc[0] << " " << epistat3[i].loc[1] << " " << epistat3[i].loc[2] << "\n";
				f3=1;
				seen++;
			}
		}
	}
}

void CheckDatEpi4 (int i, int k1, int k2, int k3,int k4, int& f4, int& seen, vector<char> sel_model, vector<epi> epistat4, vector<rec> dat) {
	for (int s=0;s<dat.size();s++) {
		if (dat[s].seq[k1]==sel_model[k1]&&dat[s].seq[k2]==sel_model[k2]&&dat[s].seq[k3]==sel_model[k3]&&dat[s].seq[k4]==sel_model[k4]) {
			int d=0;
			for (int t=0;t<dat[s].obs.size();t++) {
				if (dat[s].obs[t]>0) {
					d=1;
				}
			}
			if (f4==0&&d==1) {
				cout << "Have four-way interaction " << i << " " << epistat4[i].loc[0] << " " << epistat4[i].loc[1] << " " << epistat4[i].loc[2] << " " << epistat4[i].loc[3] << "\n";
				f4=1;
				seen++;
			}
		}
	}
}

void CheckDatEpi5 (int i, int k1, int k2, int k3,int k4, int k5, int& f5, int& seen, vector<char> sel_model, vector<epi> epistat5, vector<rec> dat) {
	for (int s=0;s<dat.size();s++) {
		if (dat[s].seq[k1]==sel_model[k1]&&dat[s].seq[k2]==sel_model[k2]&&dat[s].seq[k3]==sel_model[k3]&&dat[s].seq[k4]==sel_model[k4]&&dat[s].seq[k5]==sel_model[k5]) {
			int d=0;
			for (int t=0;t<dat[s].obs.size();t++) {
				if (dat[s].obs[t]>0) {
					d=1;
				}
			}
			if (f5==0&&d==1) {
				cout << "Have five-way interaction " << i << " " << epistat5[i].loc[0] << " " << epistat5[i].loc[1] << " " << epistat5[i].loc[2] << " " << epistat5[i].loc[3] << " " << epistat5[i].loc[4] << "\n";
				f5=1;
				seen++;
			}
		}
	}
}



void CheckEpistasisThree (int& go, int dim, vector<char> sel_model, vector<rec> dat1, vector<rec> dat2, vector<rec> dat3, vector<epi>& epistat2, vector<epi>& epistat3, vector<epi>& epistat4, vector<epi>& epistat5) {
	cout << "Check epistasis\n";
	int sugg=0;
	int seen=0;
	for (int i=0;i<epistat2.size();i++) {
		sugg++;
		int f2=0;
		for (int j=0;j<dim;j++) {
			for (int k=0;k<dim;k++) {
				if (epistat2[i].loc[0]==j&&epistat2[i].loc[1]==k) { //Have an epistatic effect of epistat2[i].x on double mutants j,k
					//cout << "j= " << j << " k= " << k << " " << sel_model[j] << " " << sel_model[k] << "\n";
					CheckDatEpi2 (i,j,k,f2,seen,sel_model,epistat2,dat1);
					CheckDatEpi2 (i,j,k,f2,seen,sel_model,epistat2,dat2);
					CheckDatEpi2 (i,j,k,f2,seen,sel_model,epistat2,dat3);
				}
			}
		}
	}
	
	for (int i=0;i<epistat3.size();i++) {
		sugg++;
		int f3=0;
		for (int k1=0;k1<dim;k1++) {
			for (int k2=0;k2<dim;k2++) {
				for (int k3=0;k3<dim;k3++) {
					if (epistat3[i].loc[0]==k1&&epistat3[i].loc[1]==k2&&epistat3[i].loc[2]==k3) { //Have an epistatic effect of epistat3[i].x on triple mutants k1,k2,k3
						CheckDatEpi3 (i,k1,k2,k3,f3,seen,sel_model,epistat3,dat1);
						CheckDatEpi3 (i,k1,k2,k3,f3,seen,sel_model,epistat3,dat2);
						CheckDatEpi3 (i,k1,k2,k3,f3,seen,sel_model,epistat3,dat3);
					}
				}
			}
		}
	}
	
	for (int i=0;i<epistat4.size();i++) {
		sugg++;
		int f4=0;
		for (int k1=0;k1<dim;k1++) {
			for (int k2=0;k2<dim;k2++) {
				for (int k3=0;k3<dim;k3++) {
					for (int k4=0;k4<dim;k4++) {
						if (epistat4[i].loc[0]==k1&&epistat4[i].loc[1]==k2&&epistat4[i].loc[2]==k3&&epistat4[i].loc[3]==k4) { //Have an epistatic effect of epistat2[i].x on double mutants j,k
							CheckDatEpi4 (i,k1,k2,k3,k4,f4,seen,sel_model,epistat4,dat1);
							CheckDatEpi4 (i,k1,k2,k3,k4,f4,seen,sel_model,epistat4,dat2);
							CheckDatEpi4 (i,k1,k2,k3,k4,f4,seen,sel_model,epistat4,dat3);

						}
					}
				}
			}
		}
	}
	
	for (int i=0;i<epistat5.size();i++) {
		sugg++;
		int f5=0;
		for (int k1=0;k1<dim;k1++) {
			for (int k2=0;k2<dim;k2++) {
				for (int k3=0;k3<dim;k3++) {
					for (int k4=0;k4<dim;k4++) {
						for (int k5=0;k5<dim;k5++) {
							if (epistat5[i].loc[0]==k1&&epistat5[i].loc[1]==k2&&epistat5[i].loc[2]==k3&&epistat5[i].loc[3]==k4&&epistat5[i].loc[4]==k5) { //Have an epistatic effect of epistat2[i].x on double mutants j,k
								CheckDatEpi5 (i,k1,k2,k3,k4,k5,f5,seen,sel_model,epistat5,dat1);
								CheckDatEpi5 (i,k1,k2,k3,k4,k5,f5,seen,sel_model,epistat5,dat2);
								CheckDatEpi5 (i,k1,k2,k3,k4,k5,f5,seen,sel_model,epistat5,dat3);

							}
						}
					}
				}
			}
		}
	}
	if (seen<sugg) {
		cout << "Not all epistatic interactions were observed\n";
		go=0;
	}
}


void CheckEpistasisEight (int& go, int dim, vector<char> sel_model, vector<rec> dat1, vector<rec> dat2, vector<rec> dat3, vector<rec> dat4, vector<rec> dat5, vector<rec> dat6, vector<rec> dat7, vector<rec> dat8, vector<epi>& epistat2, vector<epi>& epistat3, vector<epi>& epistat4, vector<epi>& epistat5) {
	cout << "Check epistasis\n";
	int sugg=0;
	int seen=0;
	for (int i=0;i<epistat2.size();i++) {
		sugg++;
		int f2=0;
		for (int j=0;j<dim;j++) {
			for (int k=0;k<dim;k++) {
				if (epistat2[i].loc[0]==j&&epistat2[i].loc[1]==k) { //Have an epistatic effect of epistat2[i].x on double mutants j,k
					//cout << "j= " << j << " k= " << k << " " << sel_model[j] << " " << sel_model[k] << "\n";
					CheckDatEpi2 (i,j,k,f2,seen,sel_model,epistat2,dat1);
					CheckDatEpi2 (i,j,k,f2,seen,sel_model,epistat2,dat2);
					CheckDatEpi2 (i,j,k,f2,seen,sel_model,epistat2,dat3);
					CheckDatEpi2 (i,j,k,f2,seen,sel_model,epistat2,dat4);
					CheckDatEpi2 (i,j,k,f2,seen,sel_model,epistat2,dat5);
					CheckDatEpi2 (i,j,k,f2,seen,sel_model,epistat2,dat6);
					CheckDatEpi2 (i,j,k,f2,seen,sel_model,epistat2,dat7);
					CheckDatEpi2 (i,j,k,f2,seen,sel_model,epistat2,dat8);
				}
			}
		}
	}
	
	for (int i=0;i<epistat3.size();i++) {
		sugg++;
		int f3=0;
		for (int k1=0;k1<dim;k1++) {
			for (int k2=0;k2<dim;k2++) {
				for (int k3=0;k3<dim;k3++) {
					if (epistat3[i].loc[0]==k1&&epistat3[i].loc[1]==k2&&epistat3[i].loc[2]==k3) { //Have an epistatic effect of epistat3[i].x on triple mutants k1,k2,k3
						CheckDatEpi3 (i,k1,k2,k3,f3,seen,sel_model,epistat3,dat1);
						CheckDatEpi3 (i,k1,k2,k3,f3,seen,sel_model,epistat3,dat2);
						CheckDatEpi3 (i,k1,k2,k3,f3,seen,sel_model,epistat3,dat3);
						CheckDatEpi3 (i,k1,k2,k3,f3,seen,sel_model,epistat3,dat4);
						CheckDatEpi3 (i,k1,k2,k3,f3,seen,sel_model,epistat3,dat5);
						CheckDatEpi3 (i,k1,k2,k3,f3,seen,sel_model,epistat3,dat6);
						CheckDatEpi3 (i,k1,k2,k3,f3,seen,sel_model,epistat3,dat7);
						CheckDatEpi3 (i,k1,k2,k3,f3,seen,sel_model,epistat3,dat8);
					}
				}
			}
		}
	}
	
	for (int i=0;i<epistat4.size();i++) {
		sugg++;
		int f4=0;
		for (int k1=0;k1<dim;k1++) {
			for (int k2=0;k2<dim;k2++) {
				for (int k3=0;k3<dim;k3++) {
					for (int k4=0;k4<dim;k4++) {
						if (epistat4[i].loc[0]==k1&&epistat4[i].loc[1]==k2&&epistat4[i].loc[2]==k3&&epistat4[i].loc[3]==k4) { //Have an epistatic effect of epistat2[i].x on double mutants j,k
							CheckDatEpi4 (i,k1,k2,k3,k4,f4,seen,sel_model,epistat4,dat1);
							CheckDatEpi4 (i,k1,k2,k3,k4,f4,seen,sel_model,epistat4,dat2);
							CheckDatEpi4 (i,k1,k2,k3,k4,f4,seen,sel_model,epistat4,dat3);
							CheckDatEpi4 (i,k1,k2,k3,k4,f4,seen,sel_model,epistat4,dat4);
							CheckDatEpi4 (i,k1,k2,k3,k4,f4,seen,sel_model,epistat4,dat5);
							CheckDatEpi4 (i,k1,k2,k3,k4,f4,seen,sel_model,epistat4,dat6);
							CheckDatEpi4 (i,k1,k2,k3,k4,f4,seen,sel_model,epistat4,dat7);
							CheckDatEpi4 (i,k1,k2,k3,k4,f4,seen,sel_model,epistat4,dat8);
							
						}
					}
				}
			}
		}
	}
	
	for (int i=0;i<epistat5.size();i++) {
		sugg++;
		int f5=0;
		for (int k1=0;k1<dim;k1++) {
			for (int k2=0;k2<dim;k2++) {
				for (int k3=0;k3<dim;k3++) {
					for (int k4=0;k4<dim;k4++) {
						for (int k5=0;k5<dim;k5++) {
							if (epistat5[i].loc[0]==k1&&epistat5[i].loc[1]==k2&&epistat5[i].loc[2]==k3&&epistat5[i].loc[3]==k4&&epistat5[i].loc[4]==k5) { //Have an epistatic effect of epistat2[i].x on double mutants j,k
								CheckDatEpi5 (i,k1,k2,k3,k4,k5,f5,seen,sel_model,epistat5,dat1);
								CheckDatEpi5 (i,k1,k2,k3,k4,k5,f5,seen,sel_model,epistat5,dat2);
								CheckDatEpi5 (i,k1,k2,k3,k4,k5,f5,seen,sel_model,epistat5,dat3);
								CheckDatEpi5 (i,k1,k2,k3,k4,k5,f5,seen,sel_model,epistat5,dat4);
								CheckDatEpi5 (i,k1,k2,k3,k4,k5,f5,seen,sel_model,epistat5,dat5);
								CheckDatEpi5 (i,k1,k2,k3,k4,k5,f5,seen,sel_model,epistat5,dat6);
								CheckDatEpi5 (i,k1,k2,k3,k4,k5,f5,seen,sel_model,epistat5,dat7);
								CheckDatEpi5 (i,k1,k2,k3,k4,k5,f5,seen,sel_model,epistat5,dat8);
							}
						}
					}
				}
			}
		}
	}
	if (seen<sugg) {
		cout << "Not all epistatic interactions were observed\n";
		go=0;
	}
}




void SigsHapSigs (int verb, int dim, run_params p, vector<int> times, vector<char> sel_model, vector<rec> dat1, vector< vector<double> > sigs, vector<epi>& epistat2, vector<epi>& epistat3, vector<epi>& epistat4, vector<epi>& epistat5, vector< vector<double> >& hapsigs) {
	for (int i=0;i<dat1.size();i++) {
		vector<double> hs;
		for (int j=0;j<times.size()-1;j++) {
			hs.push_back(0);
		}
		for (int j=0;j<dim;j++) {
			if (sigs[j].size()==1) {
				if (dat1[i].seq[j]==sel_model[j]) {
					for (int k=0;k<times.size()-1;k++) {
						hs[k]=hs[k]+sigs[j][0];
					}
				}
			}
			if (sigs[j].size()>1) {
				if (dat1[i].seq[j]==sel_model[j]) {
					for (int k=0;k<times.size()-1;k++) {
						hs[k]=hs[k]+sigs[j][k];
					}
				}
			}
		}
		hapsigs.push_back(hs);
	}
	
	//Epistasis 2-way
	for (int i=0;i<epistat2.size();i++) {
		for (int j=0;j<dim;j++) {
			for (int k=0;k<dim;k++) {
				if (epistat2[i].loc[0]==j&&epistat2[i].loc[1]==k) { //Have an epistatic effect of epistat2[i].x on double mutants j,k
					//cout << "j= " << j << " k= " << k << " " << sel_model[j] << " " << sel_model[k] << "\n";
					for (int s=0;s<dat1.size();s++) {
						if (dat1[s].seq[j]==sel_model[j]&&dat1[s].seq[k]==sel_model[k]) {
							//cout << s << " " << dat1[s].seq[j] << " " << dat1[s].seq[k] << "\n";
							for (int t=0;t<times.size()-1;t++) {
								hapsigs[s][t]=hapsigs[s][t]+epistat2[i].x;
							}
						}
					}
				}
			}
		}
	}
	
	for (int i=0;i<epistat3.size();i++) {
		for (int k1=0;k1<dim;k1++) {
			for (int k2=0;k2<dim;k2++) {
				for (int k3=0;k3<dim;k3++) {
					if (epistat3[i].loc[0]==k1&&epistat3[i].loc[1]==k2&&epistat3[i].loc[2]==k3) { //Have an epistatic effect of epistat2[i].x on double mutants j,k
						for (int s=0;s<dat1.size();s++) {
							if (dat1[s].seq[k1]==sel_model[k1]&&dat1[s].seq[k2]==sel_model[k2]&&dat1[s].seq[k3]==sel_model[k3]) {
								for (int t=0;t<times.size()-1;t++) {
									hapsigs[s][t]=hapsigs[s][t]+epistat3[i].x;
								}
							}
						}
					}
				}
			}
		}
	}

	for (int i=0;i<epistat4.size();i++) {
		for (int k1=0;k1<dim;k1++) {
			for (int k2=0;k2<dim;k2++) {
				for (int k3=0;k3<dim;k3++) {
					for (int k4=0;k4<dim;k4++) {
						if (epistat4[i].loc[0]==k1&&epistat4[i].loc[1]==k2&&epistat4[i].loc[2]==k3&&epistat4[i].loc[3]==k4) { //Have an epistatic effect of epistat2[i].x on double mutants j,k
							for (int s=0;s<dat1.size();s++) {
								if (dat1[s].seq[k1]==sel_model[k1]&&dat1[s].seq[k2]==sel_model[k2]&&dat1[s].seq[k3]==sel_model[k3]&&dat1[s].seq[k4]==sel_model[k4]) {
									for (int t=0;t<times.size()-1;t++) {
										hapsigs[s][t]=hapsigs[s][t]+epistat4[i].x;
									}
								}
							}
						}
					}
				}
			}
		}
	}

	for (int i=0;i<epistat5.size();i++) {
		for (int k1=0;k1<dim;k1++) {
			for (int k2=0;k2<dim;k2++) {
				for (int k3=0;k3<dim;k3++) {
					for (int k4=0;k4<dim;k4++) {
						for (int k5=0;k5<dim;k5++) {
							if (epistat4[i].loc[0]==k1&&epistat5[i].loc[1]==k2&&epistat5[i].loc[2]==k3&&epistat5[i].loc[3]==k4&&epistat5[i].loc[4]==k5) { //Have an epistatic effect of epistat2[i].x on double mutants j,k
								for (int s=0;s<dat1.size();s++) {
									if (dat1[s].seq[k1]==sel_model[k1]&&dat1[s].seq[k2]==sel_model[k2]&&dat1[s].seq[k3]==sel_model[k3]&&dat1[s].seq[k4]==sel_model[k4]&&dat1[s].seq[k5]==sel_model[k5]) {
										for (int t=0;t<times.size()-1;t++) {
											hapsigs[s][t]=hapsigs[s][t]+epistat5[i].x;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	
	//Non-observed state
	vector<double> hs;
	for (int j=0;j<times.size()-1;j++) {
		if (p.unobs==1) {
			hs.push_back(0);
		} else if (p.unobs==2) {
			hs.push_back(-10);
		}
	}
	hapsigs.push_back(hs);
	
	if (verb==1) {
		cout << "Hapsigs\n";
		for (int i=0;i<hapsigs.size();i++) {
			cout << "i " << i << " ";
			for (int j=0;j<hapsigs[i].size();j++) {
				cout << hapsigs[i][j] << " ";
			}
			cout << "\n";
		}
	}

}

void SigsHapSigsLinear (int verb, int dim, vector<int> times, vector<char> sel_model, vector<rec> dat1, vector< vector<double> > sigs, vector<seldiv> sd, vector<epi>& epistat2, vector<epi>& epistat3, vector<epi>& epistat4, vector<epi>& epistat5, vector< vector<double> >& hapsigs) {
	//Additive selection
	for (int i=0;i<dat1.size();i++) {
		vector<double> hs;
		for (int j=0;j<times.size()-1;j++) {
			hs.push_back(0);
		}
		for (int j=0;j<dim;j++) {
			if (sigs[j].size()==1) {
				if (dat1[i].seq[j]==sel_model[j]) {
					for (int k=0;k<times.size()-1;k++) {
						hs[k]=hs[k]+sigs[j][0];
					}
				}
			}
			if (sigs[j].size()>1) {
				if (dat1[i].seq[j]==sel_model[j]) {
					for (int k=0;k<times.size()-1;k++) {
						hs[k]=hs[k]+sigs[j][k];
					}
				}
			}
		}
		hapsigs.push_back(hs);
	}
	//cout << "Getting epistasis\n";
	
	//Epistasis
	for (int i=0;i<epistat2.size();i++) {
		for (int j=0;j<dim;j++) {
			for (int k=0;k<dim;k++) {
				if (epistat2[i].loc[0]==j&&epistat2[i].loc[1]==k) { //Have an epistatic effect of epistat2[i].x on double mutants j,k
					//cout << "j= " << j << " k= " << k << " " << sel_model[j] << " " << sel_model[k] << "\n";
					for (int s=0;s<dat1.size();s++) {
						if (dat1[s].seq[j]==sel_model[j]&&dat1[s].seq[k]==sel_model[k]) {
							//cout << s << " " << dat1[s].seq[j] << " " << dat1[s].seq[k] << "\n";
							for (int t=0;t<times.size()-1;t++) {
								hapsigs[s][t]=hapsigs[s][t]+epistat2[i].x;
							}
						}
					}
				}
			}
		}
	}
	
	for (int i=0;i<epistat3.size();i++) {
		for (int k1=0;k1<dim;k1++) {
			for (int k2=0;k2<dim;k2++) {
				for (int k3=0;k3<dim;k3++) {
					if (epistat3[i].loc[0]==k1&&epistat3[i].loc[1]==k2&&epistat3[i].loc[2]==k3) { //Have an epistatic effect of epistat3[i].x on triple mutants k1,k2,k3
						for (int s=0;s<dat1.size();s++) {
							if (dat1[s].seq[k1]==sel_model[k1]&&dat1[s].seq[k2]==sel_model[k2]&&dat1[s].seq[k3]==sel_model[k3]) {
								for (int t=0;t<times.size()-1;t++) {
									hapsigs[s][t]=hapsigs[s][t]+epistat3[i].x;
								}
							}
						}
					}
				}
			}
		}
	}
	
	for (int i=0;i<epistat4.size();i++) {
		for (int k1=0;k1<dim;k1++) {
			for (int k2=0;k2<dim;k2++) {
				for (int k3=0;k3<dim;k3++) {
					for (int k4=0;k4<dim;k4++) {
						if (epistat4[i].loc[0]==k1&&epistat4[i].loc[1]==k2&&epistat4[i].loc[2]==k3&&epistat4[i].loc[3]==k4) { //Have an epistatic effect of epistat2[i].x on double mutants j,k
							for (int s=0;s<dat1.size();s++) {
								if (dat1[s].seq[k1]==sel_model[k1]&&dat1[s].seq[k2]==sel_model[k2]&&dat1[s].seq[k3]==sel_model[k3]&&dat1[s].seq[k4]==sel_model[k4]) {
									for (int t=0;t<times.size()-1;t++) {
										hapsigs[s][t]=hapsigs[s][t]+epistat4[i].x;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	
	for (int i=0;i<epistat5.size();i++) {
		for (int k1=0;k1<dim;k1++) {
			for (int k2=0;k2<dim;k2++) {
				for (int k3=0;k3<dim;k3++) {
					for (int k4=0;k4<dim;k4++) {
						for (int k5=0;k5<dim;k5++) {
							if (epistat4[i].loc[0]==k1&&epistat5[i].loc[1]==k2&&epistat5[i].loc[2]==k3&&epistat5[i].loc[3]==k4&&epistat5[i].loc[4]==k5) { //Have an epistatic effect of epistat2[i].x on double mutants j,k
								for (int s=0;s<dat1.size();s++) {
									if (dat1[s].seq[k1]==sel_model[k1]&&dat1[s].seq[k2]==sel_model[k2]&&dat1[s].seq[k3]==sel_model[k3]&&dat1[s].seq[k4]==sel_model[k4]&&dat1[s].seq[k5]==sel_model[k5]) {
										for (int t=0;t<times.size()-1;t++) {
											hapsigs[s][t]=hapsigs[s][t]+epistat5[i].x;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	//cout << "Non-observed add " << dat1.size() << "\n";

	//Non-observed states - additive selection
	for (int i=0;i<sd.size();i++) {
		vector<double> hs;
		for (int j=0;j<times.size()-1;j++) {
			hs.push_back(0);
		}
		for (int j=0;j<dim;j++) {
			if (sigs[j].size()==1) {
				if (sd[i].seq[j]==sel_model[j]) {
					for (int k=0;k<times.size()-1;k++) {
						hs[k]=hs[k]+sigs[j][0];
					}
				}
			}
			if (sigs[j].size()>1) {
				if (sd[i].seq[j]==sel_model[j]) {
					for (int k=0;k<times.size()-1;k++) {
						hs[k]=hs[k]+sigs[j][k];
					}
				}
			}
		}
		hapsigs.push_back(hs);
	}
	//cout << "Non-observed epi\n";

	//Non-observed states - epistasis
	for (int i=0;i<epistat2.size();i++) {
		for (int j=0;j<dim;j++) {
			for (int k=0;k<dim;k++) {
				if (epistat2[i].loc[0]==j&&epistat2[i].loc[1]==k) { //Have an epistatic effect of epistat2[i].x on double mutants j,k
					//cout << "j= " << j << " k= " << k << " " << sel_model[j] << " " << sel_model[k] << "\n";
					for (int s=0;s<sd.size();s++) {
						if (sd[s].seq[j]==sel_model[j]&&sd[s].seq[k]==sel_model[k]) {
							//cout << s << " " << dat1[s].seq[j] << " " << dat1[s].seq[k] << "\n";
							for (int t=0;t<times.size()-1;t++) {
								int ss=dat1.size()+s;
								hapsigs[ss][t]=hapsigs[ss][t]+epistat2[i].x;
							}
						}
					}
				}
			}
		}
	}


	if (verb==1) {
		cout << "Hapsigs\n";
		for (int i=0;i<hapsigs.size();i++) {
			cout << "i " << i << " ";
			for (int j=0;j<hapsigs[i].size();j++) {
				cout << hapsigs[i][j] << " ";
			}
			cout << "\n";
		}
	}
	
}


void Propagation(int verb, vector<int> times, vector<double> init_freqs, vector< vector<double> > mut, vector< vector<double> > hapsigs, vector< vector<double> >& inf) {
	//cout << "Size " << init_freqs.size() << "\n";
	inf.clear();
	int start=times[0];
	int fin=times[times.size()-1];
	int index=0;
	vector<double> si;
	GetSI(index,hapsigs,si);
	vector<double> q=init_freqs;
	for (int j=0;j<=fin;j++) {
		if (j==0&&j==times[index]) {
		//	cout << "Time " << times[index] << "\n";
			if (verb==1) {
				for (int i=0;i<q.size();i++) {
					cout << q[i] << " ";
				}
				cout << "\n";
			}
			inf.push_back(q);
			index++;
		}
		if (j>=start) {
			GrowGeneration(mut,q,si);
		}
		if (j==times[index]) {
		//	cout << "Time " << times[index] << "\n";
			GetSI(index,hapsigs,si);
			index++;
			if (verb==1) {
				for (int i=0;i<q.size();i++) {
					cout << q[i] << " ";
				}
				cout << "\n";
			}
			inf.push_back(q);
		}
	}
}

void PropagationDual(int verb, vector<int> times, vector<double> init_freqs1, vector<double> init_freqs2, vector< vector<double> > mut, vector< vector<double> > hapsigs, vector< vector<double> >& inf1, vector< vector<double> >& inf2) {
	//cout << "Size " << init_freqs.size() << "\n";
	int start=times[0];
	int fin=times[times.size()-1];
	int index=0;
	vector<double> si;
	GetSI(index,hapsigs,si);

	inf1.clear();
	vector<double> q=init_freqs1;
	for (int j=0;j<=fin;j++) {
		if (j==0&&j==times[index]) {
			//	cout << "Time " << times[index] << "\n";
			if (verb==1) {
				for (int i=0;i<q.size();i++) {
					cout << q[i] << " ";
				}
				cout << "\n";
			}
			inf1.push_back(q);
			index++;
		}
		if (j>=start) {
			GrowGeneration(mut,q,si);
		}
		if (j==times[index]) {
			//	cout << "Time " << times[index] << "\n";
			GetSI(index,hapsigs,si);
			index++;
			if (verb==1) {
				for (int i=0;i<q.size();i++) {
					cout << q[i] << " ";
				}
				cout << "\n";
			}
			inf1.push_back(q);
		}
	}
	if (verb==1) {
		cout << "Second system\n";
	}
	inf2.clear();
	q.clear();
	index=0;
	q=init_freqs2;
	for (int j=0;j<=fin;j++) {
		if (j==0&&j==times[index]) {
			cout << "Time " << times[index] << "\n";
			if (verb==1) {
				for (int i=0;i<q.size();i++) {
					cout << q[i] << " ";
				}
				cout << "\n";
			}
			inf2.push_back(q);
			index++;
		}
		if (j>=start) {
			GrowGeneration(mut,q,si);
		}
		if (j==times[index]) {
			cout << "Time " << times[index] << "\n";
			GetSI(index,hapsigs,si);
			index++;
			if (verb==1) {
				for (int i=0;i<q.size();i++) {
					cout << q[i] << " ";
				}
				cout << "\n";
			}
			inf2.push_back(q);
		}
	}

}



void GetSI (int index, vector< vector<double> > hapsigs, vector<double>& si) {
	si.clear();
	for (int i=0;i<hapsigs.size();i++) {
		si.push_back(hapsigs[i][index]);
	}
}


void GrowGeneration(vector<vector<double> >& m, vector<double>& q, vector<double>& sigma) {
	vector<double> temp;
	//Change in frequencies from two rounds of growth in single cells
//	MatMult(m,q,temp);
//	q=temp;
	MatMult(m,q,temp);  //Now using m^2, calculated at the beginning...
	q=temp;
	GrowSig(q,sigma);

//	MatMult(m,q,temp);
//	q=temp;
	MatMult(m,q,temp);
	q=temp;
	GrowSig(q,sigma);
}

void MatMult (vector<vector<double> >& m, vector<double>& v, vector<double>& mv) {
	mv.clear();
	mv.assign(v.size(),0);
	for (unsigned int i=0;i<v.size();i++) {
		for (unsigned int j=0;j<v.size();j++) {
			mv[i]=mv[i]+m[i][j]*v[j];
		}
	}
}

void GrowSig (vector<double>& q, vector<double>& sigma) {
	double tot=0;
	for (unsigned int i=0;i<q.size();i++) {
		q[i]=q[i]*exp(sigma[i]);
		tot=tot+q[i];
	}
	for (unsigned int i=0;i<q.size();i++) {
		q[i]=q[i]/tot;
	}
}

double MultiCalc(int N, vector<int> obs, vector<double> inf, vector<double>& fact_store) {
	double bin=fact_store[N];
	for (unsigned int i=0;i<obs.size();i++) {
		bin=bin-fact_store[obs[i]];
	}
	for (unsigned int i=0;i<obs.size();i++) {
		bin=bin+(obs[i]*log(inf[i]));
	}
	return(bin);
}


double DirichletMultiCalc(int N, int c, vector<int> obs, vector<double> inf, vector<double>& fact_store) {
//	cout << "Obs\n";
	obs.push_back(0);
//	for (unsigned int i=0;i<obs.size();i++) {
//		cout << obs[i] << " ";
//	}
//	cout << "\n";
//	cout << "Inf\n";
//	for (unsigned int i=0;i<inf.size();i++) {
//		cout << inf[i] << " ";
//	}
//	cout << "\n";
	
	double bin=fact_store[N];
	for (unsigned int i=0;i<obs.size();i++) {
		bin=bin-fact_store[obs[i]];
	}
	vector<double> alpha;
	for (unsigned int i=0;i<inf.size();i++) {
		alpha.push_back(c*inf[i]);
	}
	double a=0;
	for (unsigned int i=0;i<alpha.size();i++) {
		a=a+alpha[i];
		bin=bin-gsl_sf_lngamma(alpha[i]);
	}
	bin=bin+gsl_sf_lngamma(a);
	a=0;
	for (unsigned int i=0;i<alpha.size();i++) {
		double b=alpha[i]+obs[i];
		a=a+b;
		bin=bin+gsl_sf_lngamma(b);
	}
	bin=bin-gsl_sf_lngamma(a);
	return(bin);
}

double DirichletMultiCalcLinear(int N, int c, vector<int> obs, vector<double> inf, vector<seldiv> sd, vector<double>& fact_store) {
//	cout << "Obs \n";
	for (int i=0;i<sd.size();i++) {
		obs.push_back(0);
	}
//	for (unsigned int i=0;i<obs.size();i++) {
//		cout << obs[i] << " ";
//	}
//	cout << "\n";
//	cout << obs.size() << "\n";
//	cout << "Inf\n";
//	for (unsigned int i=0;i<inf.size();i++) {
//		cout << inf[i] << " ";
//	}
//	cout << "\n";
//	cout << inf.size() << "\n";
	
	double bin=fact_store[N];
	for (unsigned int i=0;i<obs.size();i++) {
		bin=bin-fact_store[obs[i]];
	}
	vector<double> alpha;
	for (unsigned int i=0;i<inf.size();i++) {
		alpha.push_back(c*inf[i]);
	}
	double a=0;
	for (unsigned int i=0;i<alpha.size();i++) {
		a=a+alpha[i];
		bin=bin-gsl_sf_lngamma(alpha[i]);
	}
	bin=bin+gsl_sf_lngamma(a);
	a=0;
	for (unsigned int i=0;i<alpha.size();i++) {
		double b=alpha[i]+obs[i];
		a=a+b;
		bin=bin+gsl_sf_lngamma(b);
	}
	bin=bin-gsl_sf_lngamma(a);
	return(bin);
}


