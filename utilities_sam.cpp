#include "shared.h"
#include "utilities_sam.h"
#include <iostream>
#include <string>

//Control parameters in the code
void GetOptions (run_params& p, int argc, const char **argv) {
	string p_switch;
	p.ref="";
	p.almethod=1;
	p.min_qual=30;
	p.max_qual=0;
	p.ali_qual=30;
	p.ali_inc=0;
	p.min_rlen=30;
	p.sorted=0;
	p.qlib=1;
	p.q_cut=0.01;
	p.qp_cut=1e-3;
	p.rep_q=1;
	p.n_min=10;
	p.specify_csl=0;
	p.pairs=0;
	p.plines=3;
	p.ddup=-1;
	p.dq_cut=0.05;
	p.hap_q_cut=0.01;
	p.hap_n_min=10;
	p.seed=(int) time(NULL);
	p.skip=0;
	p.det=0;
	p.hap_index=0;
	p.conservative=0;
	p.readhap=0;
	p.verb=0;
	p.in_file="";
	p.out_file="";
	p.get_in=0;
	p.get_out=0;
	p.hap_file="../Haps1.dat";
	p.full_haps=0;
	p.full_rep=1;
	p.mu=0.33333333333e-5;
	int x=2;
	while (x < argc && (argv[x][0]=='-')) {
		p_switch=argv[x];
		if (p_switch.compare("--min_qual")==0) {
			x++;
			p.min_qual=atoi(argv[x]);
		} else if (p_switch.compare("--ref")==0) {
			x++;
			p.ref=argv[x];
		} else if (p_switch.compare("--alm")==0) {
			x++;
			p.almethod=atoi(argv[x]);
		} else if (p_switch.compare("--aliq")==0) {
			x++;
			p.ali_qual=atoi(argv[x]);
		} else if (p_switch.compare("--ali_inc")==0) {
			x++;
			p.ali_inc=atoi(argv[x]);
		} else if (p_switch.compare("--qlib")==0) {
			x++;
			p.qlib=atoi(argv[x]);
		} else if (p_switch.compare("--ddup")==0) {
			x++;
			p.ddup=atoi(argv[x]);
		} else if (p_switch.compare("--readhap")==0) {
			x++;
			p.readhap=atoi(argv[x]);
		} else if (p_switch.compare("--sorted")==0) {
			x++;
			p.sorted=atoi(argv[x]);
		} else if (p_switch.compare("--conservative")==0) {
			x++;
			p.conservative=atoi(argv[x]);
		} else if (p_switch.compare("--repq")==0) {
			x++;
			p.rep_q=atoi(argv[x]);
		} else if (p_switch.compare("--csl")==0) {
			x++;
			p.specify_csl=atof(argv[x]);
		} else if (p_switch.compare("--q_cut")==0) {
			x++;
			p.q_cut=atof(argv[x]);
		} else if (p_switch.compare("--nmin")==0) {
			x++;
			p.n_min=atoi(argv[x]);
		} else if (p_switch.compare("--pl")==0) {
			x++;
			p.plines=atoi(argv[x]);
		} else if (p_switch.compare("--hap_file")==0) {
			x++;
			p.hap_file=argv[x];
		} else if (p_switch.compare("--in")==0) {
			x++;
			p.in_file=argv[x];
			p.get_in=1;
			cout << "Set in file\n";
		} else if (p_switch.compare("--out")==0) {
			x++;
			p.out_file=argv[x];
			p.get_out=1;
		} else if (p_switch.compare("--qp_cut")==0) {
			x++;
			p.qp_cut=atof(argv[x]);
		} else if (p_switch.compare("--seed")==0) {
			x++;
			p.seed=atoi(argv[x]);
		} else if (p_switch.compare("--min_rlen")==0) {
			x++;
			p.min_rlen=atoi(argv[x]);
		} else if (p_switch.compare("--dq_cut")==0) {
			x++;
			p.dq_cut=atof(argv[x]);
		} else if (p_switch.compare("--hap_q_cut")==0) {
			x++;
			p.hap_q_cut=atof(argv[x]);
		} else if (p_switch.compare("--hap_n_min")==0) {
			x++;
			p.hap_n_min=atof(argv[x]);
		} else if (p_switch.compare("--hap_index")==0) {
			x++;
			p.hap_index=atoi(argv[x]);
		} else if (p_switch.compare("--full_haps")==0) {
			x++;
			p.full_haps=atoi(argv[x]);
		} else if (p_switch.compare("--full_rep")==0) {
			x++;
			p.full_rep=atoi(argv[x]);
		} else if (p_switch.compare("--pairs")==0) {
			x++;
			p.pairs=atoi(argv[x]);
		} else if (p_switch.compare("--detprop")==0) {
			x++;
			p.det=atoi(argv[x]);
		} else if (p_switch.compare("--verb")==0) {
			x++;
			p.verb=atoi(argv[x]);
		} else if (p_switch.compare("--minfreq")==0) {
			x++;
			p.q_cut=atof(argv[x]);
		} else if (p_switch.compare("--mu")==0) {
			x++;
			p.mu=atof(argv[x]);
		} else if (p_switch.compare("--skip")==0) {
			x++;
			p.skip=atoi(argv[x]);
		} else {
			cout << "Incorrect usage\n ";
			exit(1);
		}
		p_switch.clear();
		x++;
	}
}

//Translation between base qualities and PHRED scores.
//qlib==1: Illumina 1.3+ PHRED
//qlib==2: Illumina 1.8+ PHRED
void makequal (run_params& p, vector<char>& qual) {
	if (p.qlib==1) {
		qual.push_back('@');
		qual.push_back('A');
		qual.push_back('B');
		qual.push_back('C');
		qual.push_back('D');
		qual.push_back('E');
		qual.push_back('F');
		qual.push_back('G');
		qual.push_back('H');
		qual.push_back('I');
		qual.push_back('J');
		qual.push_back('K');
		qual.push_back('L');
		qual.push_back('M');
		qual.push_back('N');
		qual.push_back('O');
		qual.push_back('P');
		qual.push_back('Q');
		qual.push_back('R');
		qual.push_back('S');
		qual.push_back('T');
		qual.push_back('U');
		qual.push_back('V');
		qual.push_back('W');
		qual.push_back('X');
		qual.push_back('Y');
		qual.push_back('Z');
		qual.push_back('[');
		qual.push_back('\\');
		qual.push_back(']');
		qual.push_back('^');
		qual.push_back('_');
		qual.push_back('`');
		qual.push_back('a');
		qual.push_back('b');
		qual.push_back('c');
		qual.push_back('d');
		qual.push_back('e');
		qual.push_back('f');
		qual.push_back('g');
		qual.push_back('h');
		p.max_qual=40;
	}
	if (p.qlib==2) {
		qual.push_back('!');
		qual.push_back('"');
		qual.push_back('#');
		qual.push_back('$');
		qual.push_back('%');
		qual.push_back('&');
		qual.push_back('\'');
		qual.push_back('(');
		qual.push_back(')');
		qual.push_back('*');
		qual.push_back('+');
		qual.push_back(',');
		qual.push_back('-');
		qual.push_back('.');
		qual.push_back('/');
		qual.push_back('0');
		qual.push_back('1');
		qual.push_back('2');
		qual.push_back('3');
		qual.push_back('4');
		qual.push_back('5');
		qual.push_back('6');
		qual.push_back('7');
		qual.push_back('8');
		qual.push_back('9');
		qual.push_back(':');
		qual.push_back(';');
		qual.push_back('<');
		qual.push_back('=');
		qual.push_back('>');
		qual.push_back('?');
		qual.push_back('@');
		qual.push_back('A');
		qual.push_back('B');
		qual.push_back('C');
		qual.push_back('D');
		qual.push_back('E');
		qual.push_back('F');
		qual.push_back('G');
		qual.push_back('H');
		qual.push_back('I');
		qual.push_back('J');
		qual.push_back('K');
		qual.push_back('L');
		qual.push_back('M');
		qual.push_back('N');
		p.max_qual=45;
	}
}

int ScoreSim (string a, string b, int p, int rsize) {
	int s=a.size();
	int fulls=rsize;
	int sim=0;
	int l=0;

	for (int i=0;i<s;i++) {
		if ((p+i-1>=0)&&(p+i-1<fulls)) {
			l++;
			if (a[i]!=b[p+i-1]) {
				sim++;
			}
		}
	}
	
	if (l<50) {  //Minimum aligned read
		sim=50;
	}
	return sim;
}

int ScoreSim2 (string a, string b, int p, int rsize) {
	int s=a.size();
	int fulls=rsize;
	int sim=0;
	int l=0;
	cout << a << "\n";
	for (int i=0;i<s;i++) {
		if ((p+i-1>=0)&&(p+i-1<fulls)) {
			l++;
			cout << a[i] << " " << b[p+i-1] << "\n";
			if (a[i]!=b[p+i-1]) {
				sim++;
			}
		}
	}
	
	if (l<50) {  //Minimum aligned read
		sim=50;
	}
	cout << "Align2 " << l << " " << sim << "\n";
	return sim;
}

string RevTr (string a) {
	int s=a.size();
	string b=a;
	for (int i=s-1;i>=0;i--) {
		if (a[i]=='A') {
			b[s-i-1]='T';
		}
		if (a[i]=='C') {
			b[s-i-1]='G';
		}
		if (a[i]=='G') {
			b[s-i-1]='C';
		}
		if (a[i]=='T') {
			b[s-i-1]='A';
		}
		if (a[i]=='N') {
			b[s-i-1]='N';
		}
	}
	return b;
}

string RevString (string a) {
	int s=a.size();
	string b=a;
	for (int i=s-1;i>=0;i--) {
		b[s-i-1]=a[i];
	}
	return b;
}

void PrintInstructions (int a) {
	cout << "\n";
	cout << "---SAMFIRE Guide to Use---\n";
	cout << "\n";
	cout << "./process_sam filter        [options] : Performs QC on .sam format data.\n";
	cout << "./process_sam sl_traj       [options] : Identifies single-locus polymorphisms in data.\n";
	cout << "./process_sam sl_noise	     [options] : Calculates noise parameter for identified single-locus polymorphisms.\n";
	cout << "./process_sam sl_neutrality [options] : Identifies single-locus polymorphisms that evolve in a potentially non-neutral manner\n";
	cout << "./process_sam call_ml       [options] : Calls multi-locus polymorphisms\n";
	cout << "./process_sam ml_noise      [options] : Calculats noise parameter for identified multi-locus polymorphisms\n";
	cout << "./process_sam ld_calc       [options] : Calculats linkage disequilibrium statistics from multi-locus data\n";
	cout << "\n";
	cout << "See accompanying instructions for more information.\n";
	cout << "\n";
}


