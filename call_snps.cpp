#include "shared.h"
#include "call_snps.h"
#include "optimisation.h"
#include "io.h"
#include <iostream>
#include <string>
#include <sstream>


void CountNucleotides (run_params p, rseq refseq, vector<string> sam_files, vector< vector<joined> > t_reads, vector<nuc>& ref_counts) {
	if (p.verb==1) {
		cout << "Size of reference sequence " << refseq.seq.size() << "\n";
	}
	SetupRefCounts(refseq,sam_files,ref_counts);
	for (int i=0;i<t_reads.size();i++) {  //Dataset
		cout << t_reads[i].size() << "\n";
		for (int j=0;j<t_reads[i].size();j++) {  //Read
		//	cout << j << " " << t_reads[i][j].alpos << "\n";
			for (int k=0;k<t_reads[i][j].seq.size();k++) {
				if (t_reads[i][j].seq[k]=='A') {
					ref_counts[i].nA[t_reads[i][j].alpos+k]++;
					ref_counts[i].nN[t_reads[i][j].alpos+k]++;
				} else if (t_reads[i][j].seq[k]=='C') {
					ref_counts[i].nC[t_reads[i][j].alpos+k]++;
					ref_counts[i].nN[t_reads[i][j].alpos+k]++;
				} else if (t_reads[i][j].seq[k]=='G') {
					ref_counts[i].nG[t_reads[i][j].alpos+k]++;
					ref_counts[i].nN[t_reads[i][j].alpos+k]++;
				} else if (t_reads[i][j].seq[k]=='T') {
					ref_counts[i].nT[t_reads[i][j].alpos+k]++;
					ref_counts[i].nN[t_reads[i][j].alpos+k]++;
				}
			}
		}
	}
}

void SetupRefCounts (rseq refseq, vector<string> sam_files, vector<nuc>& ref_counts) {
	for (int t=0;t<sam_files.size();t++) {
		nuc counts;
		for (int i=0;i<refseq.size;i++) {
			counts.nN.push_back(0);
			counts.nA.push_back(0);
			counts.nC.push_back(0);
			counts.nG.push_back(0);
			counts.nT.push_back(0);
		}
		ref_counts.push_back(counts);
	}
}

void PerSequenceSNPs (run_params p, rseq refseq, vector<string> sam_files) {
	for (int i=0;i<sam_files.size();i++) {
		vector<joined> t_read;
		InputJnData (i,t_read);
		ofstream mut_file;
		ostringstream convert;
		convert << i;
		string temp=convert.str();
		string name = "Mutations"+temp+".out";
		mut_file.open(name.c_str());
		for (int j=0;j<t_read.size();j++) {
			int index=0;
			int a_index=0;
			int mut=0;
			if (p.pos==-1) {
				for (int k=t_read[j].alpos-1;k<t_read[j].alpos+t_read[j].seq.size()-1;k++) {
					if (t_read[j].seq[index]!='-'&&t_read[j].seq[index]!=refseq.seq[k]) {
						mut++;
					}
					index++;
				}
				mut_file << index << " " << mut << "\n";
			} else {
				if (t_read[j].alpos<=p.pos-1&&t_read[j].alpos+t_read[j].seq.size()>=p.pos-1) {
					char all='-';
					for (int k=t_read[j].alpos-1;k<t_read[j].alpos+t_read[j].seq.size()-1;k++) {
						if (k==p.pos-1) {
							all=t_read[j].seq[index];
						}
						if (t_read[j].seq[index]!='-'&&t_read[j].seq[index]!=refseq.seq[k]) {
							mut++;
						}
						index++;
						if (t_read[j].seq[index]!='-') {
							a_index++;
						}
					}
					mut_file << a_index << " " << all << " " << mut << "\n";
				}
			}
		}
	}
}


void CountNucs (run_params p, rseq refseq, vector<joined> t_read, nuc& r_count) {
	if (p.verb==1) {
		cout << "Size of reference sequence " << refseq.seq.size() << "\n";
	}
	SetupRCount(refseq,r_count);
	cout << t_read.size() << "\n";
	for (int j=0;j<t_read.size();j++) {  //Read
		//	cout << j << " " << t_reads[i][j].alpos << "\n";
		for (int k=0;k<t_read[j].seq.size();k++) {
			if (t_read[j].seq[k]=='A') {
				r_count.nA[t_read[j].alpos+k]++;
				r_count.nN[t_read[j].alpos+k]++;
			} else if (t_read[j].seq[k]=='C') {
				r_count.nC[t_read[j].alpos+k]++;
				r_count.nN[t_read[j].alpos+k]++;
			} else if (t_read[j].seq[k]=='G') {
				r_count.nG[t_read[j].alpos+k]++;
				r_count.nN[t_read[j].alpos+k]++;
			} else if (t_read[j].seq[k]=='T') {
				r_count.nT[t_read[j].alpos+k]++;
				r_count.nN[t_read[j].alpos+k]++;
			}
		}
	}
}

void SetupRCount (rseq refseq, nuc& r_count) {
	for (int i=0;i<refseq.size;i++) {
		r_count.nN.push_back(0);
		r_count.nA.push_back(0);
		r_count.nC.push_back(0);
		r_count.nG.push_back(0);
		r_count.nT.push_back(0);
	}
}

void CallPolymorphismsVsRef (run_params p, rseq refseq, vector<nuc> ref_counts, vector<poly>& polys) {
	if (p.verb==1) {
		cout << "Calling polymorphisms versus reference sequence...\n";
	}
	vector<char> consensus;
	//Define consensus as the reference sequence
	consensus.push_back('X');
	for (int j=0;j<refseq.size;j++) {
		char n=refseq.seq[j];
		cout << n;
		consensus.push_back(n);
	}
	cout << "\n";
	cout << "Consensus size " << consensus.size() << "\n";
	for (int i=0;i<ref_counts.size();i++) { //Time points
		if (p.verb==1) {
			cout << "Time point #" << i << "\n";
		}
		for (int j=0;j<refseq.size;j++) {
			//Here, know the consensus.  For all other nucleotides, check to see whether the SNP criteria are fulfilled
			double fA=(ref_counts[i].nA[j]+0.)/(ref_counts[i].nN[j]+0.);
			double fC=(ref_counts[i].nC[j]+0.)/(ref_counts[i].nN[j]+0.);
			double fG=(ref_counts[i].nG[j]+0.)/(ref_counts[i].nN[j]+0.);
			double fT=(ref_counts[i].nT[j]+0.)/(ref_counts[i].nN[j]+0.);
			char cons=consensus[j];
			if (consensus[j]=='A') {
				vector<double> freqs;
				vector<char> nuc;
				vector<int> counts;
				freqs.push_back(fC);
				freqs.push_back(fG);
				freqs.push_back(fT);
				nuc.push_back('C');
				nuc.push_back('G');
				nuc.push_back('T');
				counts.push_back(ref_counts[i].nC[j]);
				counts.push_back(ref_counts[i].nG[j]);
				counts.push_back(ref_counts[i].nT[j]);
				ProcessPVsRef(p,j,freqs,nuc,cons,counts,ref_counts[i].nN[j],polys);
			}
			if (consensus[j]=='C') {
				vector<double> freqs;
				vector<char> nuc;
				vector<int> counts;
				freqs.push_back(fA);
				freqs.push_back(fG);
				freqs.push_back(fT);
				nuc.push_back('A');
				nuc.push_back('G');
				nuc.push_back('T');
				counts.push_back(ref_counts[i].nA[j]);
				counts.push_back(ref_counts[i].nG[j]);
				counts.push_back(ref_counts[i].nT[j]);
				ProcessPVsRef(p,j,freqs,nuc,cons,counts,ref_counts[i].nN[j],polys);
			}
			if (consensus[j]=='G') {
				vector<double> freqs;
				vector<char> nuc;
				vector<int> counts;
				freqs.push_back(fA);
				freqs.push_back(fC);
				freqs.push_back(fT);
				nuc.push_back('A');
				nuc.push_back('C');
				nuc.push_back('T');
				counts.push_back(ref_counts[i].nA[j]);
				counts.push_back(ref_counts[i].nC[j]);
				counts.push_back(ref_counts[i].nT[j]);
				ProcessPVsRef(p,j,freqs,nuc,cons,counts,ref_counts[i].nN[j],polys);
			}
			if (consensus[j]=='T') {
				vector<double> freqs;
				vector<char> nuc;
				vector<int> counts;
				freqs.push_back(fA);
				freqs.push_back(fC);
				freqs.push_back(fG);
				nuc.push_back('A');
				nuc.push_back('C');
				nuc.push_back('G');
				counts.push_back(ref_counts[i].nA[j]);
				counts.push_back(ref_counts[i].nC[j]);
				counts.push_back(ref_counts[i].nG[j]);
				ProcessPVsRef(p,j,freqs,nuc,cons,counts,ref_counts[i].nN[j],polys);
			}
		}
	}
	//Delete duplicate polymorphisms
	DeleteDuplicatePolymporphisms(p,polys);
}

void DeleteDuplicatePolymporphisms (run_params p, vector<poly>& polys) {
	vector<poly> poly2;
	vector<int> x; //Vectors of previously observed trajectories
	vector<char> c;
	for (int i=0;i<polys.size();i++) {
		int inc=0;
		for (int j=0;j<x.size();j++) { //Previously recorded polymorphisms
			if (polys[i].locus==x[j]&&polys[i].nuc==c[j]) {
				inc++;
			}
		}
		x.push_back(polys[i].locus);
		c.push_back(polys[i].nuc);
		if (inc==p.rep_q-1) { //Push back only once after required number of observations
			poly2.push_back(polys[i]);
		}
	}
	polys=poly2;
}


void ProcessPVsRef (run_params p, int j, vector<double> freqs, vector<char> nuc, char cons, vector<int> counts, int N, vector<poly>& polys) {
	double prob=-(p.min_qual+0.)/10.;
	prob=pow(10,prob);
	for (int i=0;i<3;i++) { //Cycle through non-consensus alleles
		double x=freqs[i];
		int inc=0;
		poly y;
		y.locus=j;
		y.cons=cons;
		y.nuc=nuc[i];
		if (x>p.q_cut) { //Frequency threshold
			if (p.gmaf==1||x<1-p.q_cut) { //gmaf=1 calls fixation events relative to the reference sequence
				inc=1;
				if (counts[i]<p.n_min) { //Absolute nucleotide count threshold
					inc=0;
				} else {
					double pr=0;
					for (int k=counts[i];k<=N;k++) {
						pr=pr+gsl_ran_binomial_pdf(k,prob,N);
					}
					if (pr>p.qp_cut) { //Probability threshold
						inc=0;
					}
				}
			}
		}
		if (inc==1) {
			polys.push_back(y);
			if (p.verb==1) {
				cout << "Polymorphism " << y.locus << " " << y.cons << " " << y.nuc << " " << x << "\n";
			}
		}
	}
}



void CallPolymorphisms (run_params p, rseq refseq, vector<nuc> ref_counts, vector<poly>& polys) {
	if (p.verb==1) {
		cout << "Calling polymorphisms...\n";
	}
	//Define consensus sequence from most common base in first time point
	vector<char> consensus;
	for (int j=0;j<refseq.size;j++) {
		int t0=0;
		while (t0<ref_counts.size()&&ref_counts[t0].nN[j]==0) {
			t0++;
		}
		if (t0<ref_counts.size()) {
			char n='A';
			if (ref_counts[t0].nC[j]>ref_counts[t0].nA[j]) {
				if (ref_counts[t0].nG[j]>ref_counts[t0].nC[j]) {
					if (ref_counts[t0].nT[j]>ref_counts[t0].nG[j]) {
						n='T';
					} else {
						n='G';
					}
				} else if (ref_counts[t0].nT[j]>ref_counts[t0].nC[j]) {
					n='T';
				} else {
					n='C';
				}
			} else if (ref_counts[t0].nG[j]>ref_counts[t0].nA[j]) {
				if (ref_counts[t0].nT[j]>ref_counts[t0].nG[j]) {
					n='T';
				} else {
					n='G';
				}
			} else if (ref_counts[t0].nT[j]>ref_counts[t0].nA[j]) {
				n='T';
			} else {
				n='A';
			}
			consensus.push_back(n);
		} else {
			consensus.push_back('N');
		}
	}
	
	for (int i=0;i<ref_counts.size();i++) {
		if (p.verb==1) {
			cout << "Time point #" << i << "\n";
		}
		for (int j=0;j<refseq.size;j++) {
			//Here, know the consensus.  For all other nucleotides, check to see whether the SNP criteria are fulfilled
			double fA=(ref_counts[i].nA[j]+0.)/(ref_counts[i].nN[j]+0.);
			double fC=(ref_counts[i].nC[j]+0.)/(ref_counts[i].nN[j]+0.);
			double fG=(ref_counts[i].nG[j]+0.)/(ref_counts[i].nN[j]+0.);
			double fT=(ref_counts[i].nT[j]+0.)/(ref_counts[i].nN[j]+0.);
			char cons=consensus[j];
			if (consensus[j]=='A') {
				vector<double> freqs;
				vector<char> nuc;
				vector<int> counts;
				freqs.push_back(fC);
				freqs.push_back(fG);
				freqs.push_back(fT);
				nuc.push_back('C');
				nuc.push_back('G');
				nuc.push_back('T');
				counts.push_back(ref_counts[i].nC[j]);
				counts.push_back(ref_counts[i].nG[j]);
				counts.push_back(ref_counts[i].nT[j]);
				ProcessP(p,j,freqs,nuc,cons,counts,ref_counts[i].nN[j],polys);
			}
			if (consensus[j]=='C') {
				vector<double> freqs;
				vector<char> nuc;
				vector<int> counts;
				freqs.push_back(fA);
				freqs.push_back(fG);
				freqs.push_back(fT);
				nuc.push_back('A');
				nuc.push_back('G');
				nuc.push_back('T');
				counts.push_back(ref_counts[i].nA[j]);
				counts.push_back(ref_counts[i].nG[j]);
				counts.push_back(ref_counts[i].nT[j]);
				ProcessP(p,j,freqs,nuc,cons,counts,ref_counts[i].nN[j],polys);
			}
			if (consensus[j]=='G') {
				vector<double> freqs;
				vector<char> nuc;
				vector<int> counts;
				freqs.push_back(fA);
				freqs.push_back(fC);
				freqs.push_back(fT);
				nuc.push_back('A');
				nuc.push_back('C');
				nuc.push_back('T');
				counts.push_back(ref_counts[i].nA[j]);
				counts.push_back(ref_counts[i].nC[j]);
				counts.push_back(ref_counts[i].nT[j]);
				ProcessP(p,j,freqs,nuc,cons,counts,ref_counts[i].nN[j],polys);
			}
			if (consensus[j]=='T') {
				vector<double> freqs;
				vector<char> nuc;
				vector<int> counts;
				freqs.push_back(fA);
				freqs.push_back(fC);
				freqs.push_back(fG);
				nuc.push_back('A');
				nuc.push_back('C');
				nuc.push_back('G');
				counts.push_back(ref_counts[i].nA[j]);
				counts.push_back(ref_counts[i].nC[j]);
				counts.push_back(ref_counts[i].nG[j]);
				ProcessP(p,j,freqs,nuc,cons,counts,ref_counts[i].nN[j],polys);
			}
		}
	}
	
	//Delete duplicate polymorphisms
	DeleteDuplicatePolymporphisms(p,polys);
}

void ProcessP (run_params p, int j, vector<double> freqs, vector<char> nuc, char cons, vector<int> counts, int N, vector<poly>& polys) {
	double prob=-(p.min_qual+0.)/10.;
	prob=pow(10,prob);
	for (int i=0;i<3;i++) {
		double x=freqs[i];
		int inc=0;
		poly y;
		y.locus=j;
		y.cons=cons;
		y.nuc=nuc[i];
		if (x>p.q_cut&&x<1-p.q_cut) { //Frequency threshold
			inc=1;
			if (counts[i]<p.n_min) { //Absolute nucleotide count threshold
				inc=0;
			} else {
				double pr=0;
				for (int k=counts[i];k<=N;k++) {
					pr=pr+gsl_ran_binomial_pdf(k,prob,N);
				}
				if (pr>p.qp_cut) { //Probability threshold
					inc=0;
				}
			}
		}
		if (inc==1) {
			polys.push_back(y);
			if (p.verb==1) {
				cout << "Polymorphism " << y.locus << " " << y.cons << " " << y.nuc << " " << x << "\n";
			}
		}
	}
}

void ConstructSLTrajs (run_params p, vector<poly> polys, vector<nuc> ref_counts, vector<str>& sltrajs) {
	if (p.verb==1) {
		cout << "Constructing trajectories...\n";
	}
	vector<int> times;
	ImportTimeData(times);
	for (int i=0;i<polys.size();i++) {
		str sltraj;
		sltraj.times=times;
		sltraj.locus=polys[i].locus;
		sltraj.nuc=polys[i].nuc;
		sltraj.cons=polys[i].cons;
		for (int j=0;j<ref_counts.size();j++) {
			if (p.verb==1) {
				cout << polys[i].locus << " " << polys[i].cons << " " << polys[i].nuc << " " << j << " " << ref_counts[j].nA[polys[i].locus] << " " << ref_counts[j].nC[polys[i].locus] << " " << ref_counts[j].nG[polys[i].locus] << " " << ref_counts[j].nT[polys[i].locus] << " " << ref_counts[j].nN[polys[i].locus] << "\n";
			}
			sltraj.nA.push_back(ref_counts[j].nA[polys[i].locus]);
			sltraj.nC.push_back(ref_counts[j].nC[polys[i].locus]);
			sltraj.nG.push_back(ref_counts[j].nG[polys[i].locus]);
			sltraj.nT.push_back(ref_counts[j].nT[polys[i].locus]);
			sltraj.nN.push_back(ref_counts[j].nN[polys[i].locus]);
			double qA=(ref_counts[j].nA[polys[i].locus]+0.)/(ref_counts[j].nN[polys[i].locus]+0.);
			double qC=(ref_counts[j].nC[polys[i].locus]+0.)/(ref_counts[j].nN[polys[i].locus]+0.);
			double qG=(ref_counts[j].nG[polys[i].locus]+0.)/(ref_counts[j].nN[polys[i].locus]+0.);
			double qT=(ref_counts[j].nT[polys[i].locus]+0.)/(ref_counts[j].nN[polys[i].locus]+0.);
			sltraj.qA.push_back(qA);
			sltraj.qC.push_back(qC);
			sltraj.qG.push_back(qG);
			sltraj.qT.push_back(qT);
			sltraj.inc=1;
		}
		sltrajs.push_back(sltraj);
	}
	if (p.verb==1) {
		cout << "Number of trajectories " << sltrajs.size() << "\n";
	}
}

void SLTFreqs (run_params p, vector<str>& sltrajs) {
	for (int i=0;i<sltrajs.size();i++) {
		for (int j=0;j<sltrajs[i].nA.size();j++) {
			double q=(sltrajs[i].nA[j]+0.)/(sltrajs[i].nN[j]+0.);
			sltrajs[i].qA.push_back(q);
			q=(sltrajs[i].nC[j]+0.)/(sltrajs[i].nN[j]+0.);
			sltrajs[i].qC.push_back(q);
			q=(sltrajs[i].nG[j]+0.)/(sltrajs[i].nN[j]+0.);
			sltrajs[i].qG.push_back(q);
			q=(sltrajs[i].nT[j]+0.)/(sltrajs[i].nN[j]+0.);
			sltrajs[i].qT.push_back(q);
			/*if (p.verb==1) {
				cout << sltrajs[i].nA[j] << " " << sltrajs[i].nC[j] << " " << sltrajs[i].nG[j] << " " << sltrajs[i].nT[j] << "\n";
				cout << sltrajs[i].qA[j] << " " << sltrajs[i].qC[j] << " " << sltrajs[i].qG[j] << " " << sltrajs[i].qT[j] << "\n";
			}*/
		}
	}
}

void SLTMeanFreqs (run_params p, vector<str>& sltrajs) {
	if (p.verb==1) {
		cout << "Calculating mean polymorphism frequencies over time...\n";
	}
	//cout << "Cutoff " << p.dep_cut << "\n";
	for (int i=0;i<sltrajs.size();i++) {
		double N=0;
		sltrajs[i].mA=0;
		sltrajs[i].mC=0;
		sltrajs[i].mG=0;
		sltrajs[i].mT=0;
		for (int j=0;j<sltrajs[i].qA.size();j++) {
			if (sltrajs[i].nN[j]>p.dep_cut) {
				sltrajs[i].mA=sltrajs[i].mA+sltrajs[i].nA[j];
				sltrajs[i].mC=sltrajs[i].mC+sltrajs[i].nC[j];
				sltrajs[i].mG=sltrajs[i].mG+sltrajs[i].nG[j];
				sltrajs[i].mT=sltrajs[i].mT+sltrajs[i].nT[j];
				N=N+sltrajs[i].nN[j];
				if (p.verb==1) {
					cout << i << " " << j << " " << sltrajs[i].nA[j] << " " << sltrajs[i].nC[j] << " " << sltrajs[i].nG[j] << " " << sltrajs[i].nT[j] << " " << sltrajs[i].nN[j] << "\n";
				}

			}
		}
		if (N>1) {
			sltrajs[i].mA=sltrajs[i].mA/N;
			sltrajs[i].mC=sltrajs[i].mC/N;
			sltrajs[i].mG=sltrajs[i].mG/N;
			sltrajs[i].mT=sltrajs[i].mT/N;
		}
		if (p.verb==1) {
			cout << i << " " << sltrajs[i].mA << " " << sltrajs[i].mC << " " << sltrajs[i].mG << " " << sltrajs[i].mT << "\n";
		}
	}
}

void FilterSLTrajs (run_params p, vector<str>& sltrajs) {
	if (p.verb==1) {
		cout << "Filter trajectories by mean change in frequency with time...\n";
	}
	for (int i=0;i<sltrajs.size();i++) {
		double dA=0;
		double dC=0;
		double dG=0;
		double dT=0;
		double q=0;
		double deltat=0;
		//Mean change in frequency per day
		for (int j=0;j<sltrajs[i].qA.size()-1;j++) {
			//Absolute frequency changes; occasional compiler issue with abs()
			q=sltrajs[i].qA[j+1]-sltrajs[i].qA[j];
			if (q<0) {q=-q;}
			dA=dA+q;
			q=sltrajs[i].qC[j+1]-sltrajs[i].qC[j];
			if (q<0) {q=-q;}
			dC=dC+q;
			q=sltrajs[i].qG[j+1]-sltrajs[i].qG[j];
			if (q<0) {q=-q;}
			dG=dG+q;
			q=sltrajs[i].qT[j+1]-sltrajs[i].qT[j];
			if (q<0) {q=-q;}
			dT=dT+q;
			deltat=deltat+(sltrajs[i].times[j+1]-sltrajs[i].times[j]);
		}
		//cout << "Trajectory " << i << " " << dA << " " << dC << " " << dG << " " << dT << " " << deltat << "\n";
		dA=dA/deltat;
		dC=dC/deltat;
		dG=dG/deltat;
		dT=dT/deltat;
		double m1=max(dA,dC);
		double m2=max(dG,dT);
		double m=max(m1,m2);
	//	m=m/deltat;
		//cout << " Move " << m << "\n";
		if (m>p.dq_cut) {  //Remove trajectories that move more than the cutoff
			sltrajs[i].inc=0;
			if (p.verb==1) {
				cout << "Remove trajectory " << i << "\n";
			}
		}
	}
}

void FilterSLTrajs2 (run_params p, vector<str>& sltrajs) {
	if (p.verb==1) {
		cout << "Filter trajectories by frequency cutoff...\n";
	}
	vector<str> new_sltrajs;
	for (int i=0;i<sltrajs.size();i++) {
		if (sltrajs[i].inc==1) {
		//	cout << "Check " << i << "\n";
			str s;
			for (int j=0;j<sltrajs[i].qA.size();j++) {
				if (sltrajs[i].nN[j]>p.dep_cut) {
					int inc=0;
					if (sltrajs[i].qA[j]>p.q_cut&&sltrajs[i].qA[j]<(1-p.q_cut)) {
						inc=1;
					}
					if (sltrajs[i].qC[j]>p.q_cut&&sltrajs[i].qC[j]<(1-p.q_cut)) {
						inc=1;
					}
					if (sltrajs[i].qG[j]>p.q_cut&&sltrajs[i].qG[j]<(1-p.q_cut)) {
						inc=1;
					}
					if (sltrajs[i].qT[j]>p.q_cut&&sltrajs[i].qT[j]<(1-p.q_cut)) {
						inc=1;
					}
					//if (p.verb==1) {
					//	cout << sltrajs[i].qA[j] << " " << sltrajs[i].qC[j] << " " << sltrajs[i].qG[j] << " " << sltrajs[i].qT[j] << "\n";
					//}
					if (inc==1) {
						s.qA.push_back(sltrajs[i].qA[j]);
						s.qC.push_back(sltrajs[i].qC[j]);
						s.qG.push_back(sltrajs[i].qG[j]);
						s.qT.push_back(sltrajs[i].qT[j]);
						s.nA.push_back(sltrajs[i].nA[j]);
						s.nC.push_back(sltrajs[i].nC[j]);
						s.nG.push_back(sltrajs[i].nG[j]);
						s.nT.push_back(sltrajs[i].nT[j]);
						s.nN.push_back(sltrajs[i].nN[j]);
						s.inc=1;
					}
				}
			}
		//	cout << "Size " << s.qA.size() << "\n";
			if (s.qA.size()>1) {
				new_sltrajs.push_back(s);
			} else {
				if (p.verb==1) {
					cout << "Remove trajectory " << i << "\n";
				}
			}
		}
	}
	if (p.verb==1) {
		cout << "Number of remaining trajectories " << new_sltrajs.size() << "\n";
	}
	sltrajs=new_sltrajs;
}


void PotentialNonNeutrality (run_params p, double Csl_opt, vector<str>& sltrajs, vector<double> fact_store, gsl_rng *rgen) {
	cout << "Calculating fit under single-locus models...\n";
	cout << sltrajs.size() << "\n";
	for (int traj=0;traj<sltrajs.size();traj++) {
		cout << "Trajectory " << traj << " locus " << sltrajs[traj].locus << " ";
		
		int n_times=sltrajs[traj].times.size();
		vector<int> times = sltrajs[traj].times;
		vector<int> N=sltrajs[traj].nN;
		double Ntot=0;
		for (int i=0;i<N.size();i++) {
			Ntot=Ntot+N[i];
		}
		
		cout << " Neutral model ";
		char sel_model='0';
		int tds=0;
		double L;
		double bic;
		L=OptimiseSLTraj(p.verb,p,1,traj,sel_model,tds,n_times,times,N,Csl_opt,sltrajs,fact_store,rgen);
		sltrajs[traj].logL.push_back(L);
		bic=(-2*L);
		cout << "BIC " << bic << " ";
		sltrajs[traj].BIC.push_back(bic);
		cout << "Constant selection ";
		sel_model=sltrajs[traj].nuc;
		L=OptimiseSLTraj(p.verb,p,5,traj,sel_model,tds,n_times,times,N,Csl_opt,sltrajs,fact_store,rgen);
		sltrajs[traj].logL.push_back(L);
		bic=(-2*L)+log(Ntot);
		cout << " BIC " << bic << " ";
		sltrajs[traj].BIC.push_back(bic);
		
		cout << "Time-dependent selection ";
		tds=1;
		L=OptimiseSLTraj(p.verb,p,5,traj,sel_model,tds,n_times,times,N,Csl_opt,sltrajs,fact_store,rgen);
		sltrajs[traj].logL.push_back(L);
		bic=(-2*L)+(times.size()*log(Ntot));
		cout << "BIC " << bic << "\n";
		sltrajs[traj].BIC.push_back(bic);
		
		if (sltrajs[traj].BIC[1]<sltrajs[traj].BIC[0] || sltrajs[traj].BIC[2]<sltrajs[traj].BIC[0]) {
			if (sltrajs[traj].BIC[2]<sltrajs[traj].BIC[1]) {
				cout << "Model T is best fit\n";
			} else {
				cout << "Model S is best fit\n";
			}
		} else {
			cout << "Model 0 is best fit\n";
			sltrajs[traj].inc=0;
		}
	}
}






