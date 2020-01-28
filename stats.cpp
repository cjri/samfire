#include "shared.h"
#include "stats.h"
#include "io.h"
#include <iostream>
#include <string>

void Calc_Seq_Lengths (vector<joined> t_read, vector<int>& l_dist) {
	vector<int> lengths;
	int max=0;
	for (int i=0;i<t_read.size();i++) {
		lengths.push_back(t_read[i].seq.size());
		if (t_read[i].seq.size()>max) {
			max=t_read[i].seq.size();
		}
	}
	for (int i=0;i<=max;i++) {
		l_dist.push_back(0);
	}
	for (int i=0;i<lengths.size();i++) {
		l_dist[lengths[i]]++;
	}
}

void Calc_Base_Quality (vector<ql> q_read, vector<long>& q_dist) {
	for (int i=0;i<45;i++) {
		q_dist.push_back(0);
	}
	for (int j=0;j<q_read.size();j++) {
		for (int k=0;k<q_read[j].seq.size();k++) {
			if (q_read[j].seq[k]>q_dist.size()) {
				cout << "Base quality error " << q_read[j].seq[k] << "\n";
			} else {
				q_dist[q_read[j].seq[k]]++;
			}
		}
	}
}

void MakeBootstraps (vector< vector< vector<int> > >& all_nucs, vector< vector< vector< vector<int> > > >& bootstrap, gsl_rng *rgen) {
	cout << "Generating bootstraps\n";
	for (int b=0;b<100;b++) {
		vector< vector< vector<int> > > bs3;
		for (int i=0;i<all_nucs.size();i++) {
			vector< vector<int> > bs2;
			for (int j=0;j<all_nucs[i].size();j++) {
				int n=all_nucs[i][j][4];
				double* prob_select_mirror = new double[4];
				unsigned int* draw = new unsigned int[4];
				vector<int> bs1;
				prob_select_mirror[0]=(all_nucs[i][j][0]+0.)/(n+0.);
				prob_select_mirror[1]=(all_nucs[i][j][1]+0.)/(n+0.);
				prob_select_mirror[2]=(all_nucs[i][j][2]+0.)/(n+0.);
				prob_select_mirror[3]=(all_nucs[i][j][3]+0.)/(n+0.);
				gsl_ran_multinomial(rgen,4,n,prob_select_mirror,draw);
				for (int i=0;i<4;i++) {
					bs1.push_back(draw[i]);
				}
				bs1.push_back(n);
				bs2.push_back(bs1);
			}
			bs3.push_back(bs2);
		}
		bootstrap.push_back(bs3);
	}
}

void Calc_Pi_Diversity (run_params p, vector< vector< vector<int> > >& all_nucs, vector< vector< vector< vector<int> > > >& bootstrap) {
	cout << "Calculating diversity statistic pi\n";
	vector<double> pi;
	PiCalculation(all_nucs,pi);
	ofstream pi_file;
	pi_file.open("Pi_diversity.dat");
	for (int i=0;i<pi.size();i++) {
		pi_file << pi[i] << " ";
	}
	pi_file << "\n";
	pi_file.close();
	if (p.bootstrap==1) {
		pi_file.open("Pi_bootstraps.dat");
		for (int i=0;i<bootstrap.size();i++) {
			vector<double> pi;
			PiCalculation(bootstrap[i],pi);
			for (int i=0;i<pi.size();i++) {
				pi_file << pi[i] << " ";
			}
			pi_file << "\n";

		}
	}
	
	if (p.sns==1) {
		cout << "Recalculating pi over four-fold S and NS sites\n";

		//Import SNS types
		vector< vector<int> > var;
		InputVariantTypes(var);
		
		//Make copies of nucleotide count files: Initial setup
		vector< vector< vector<int> > > all_nucs_s;
		vector< vector< vector<int> > > all_nucs_ns;
		for (int i=0;i<all_nucs.size();i++) {
			vector< vector<int> > anucs;
			all_nucs_s.push_back(anucs);
			all_nucs_ns.push_back(anucs);
		}
		
		//Find four-fold variant sites
		Construct_FourFold_Variants(var,all_nucs,all_nucs_s,all_nucs_ns);
		cout << "Found " << all_nucs_s[0].size() << " S sites and " << all_nucs_ns[0].size() << " NS sites\n";

		pi.clear();
		PiCalculation(all_nucs_s,pi);
		pi_file.open("Pi_S_diversity.dat");
		for (int i=0;i<pi.size();i++) {
			pi_file << pi[i] << " ";
		}
		pi_file << "\n";
		pi_file.close();
		pi.clear();
		PiCalculation(all_nucs_ns,pi);
		pi_file.open("Pi_NS_diversity.dat");
		for (int i=0;i<pi.size();i++) {
			pi_file << pi[i] << " ";
		}
		pi_file << "\n";
		pi_file.close();
	}
}

void Construct_FourFold_Variants (vector< vector<int> >& var, vector< vector< vector<int> > >& all_nucs, vector< vector< vector<int> > >& all_nucs_s, vector< vector< vector<int> > >& all_nucs_ns) {
	for (int i=0;i<var.size();i++) {
		//cout << i << " " << var.size() << "\n";
		if (var[i][1]!=3) { //Ignore non-coding positions
			//Count S and NS sites
			int count_s=0;
			int count_ns=0;
			//cout << var[i][1] << " " << var[i][2] << " " << var[i][3] << " " << var[i][4] << "\n";
			for (int j=1;j<=4;j++) {
				if (var[i][j]==1) {
					count_s++;
				}
				if (var[i][j]==2) {
					count_ns++;
				}
			}
			if (count_ns==3) {
				//Push this line back to all_nucs_ns
				for (int j=0;j<all_nucs.size();j++) {
					all_nucs_ns[j].push_back(all_nucs[j][i]);
				}
			}
			if (count_s==3) {
				//Push this line back to all_nucs_ns
				for (int j=0;j<all_nucs.size();j++) {
					all_nucs_s[j].push_back(all_nucs[j][i]);
				}
			}
		}
	}
}

void Construct_AllFreqsSNS_Variants (vector< vector<int> >& var, vector< vector< vector<int> > >& all_nucs, vector< vector< vector<int> > >& all_nucs_s, vector< vector< vector<int> > >& all_nucs_ns) {
	//Designed for count of total frequencies
	//Retains all sites containing an S or NS polymorphism
	//Sets non-compliant variants to zero
	//Retains original read counts N
	for (int i=0;i<var.size();i++) {
		//cout << i << " " << var.size() << "\n";
		if (var[i][1]!=3) { //Ignore non-coding positions
			//Count S and NS sites
			int count_s=0;
			int count_ns=0;
			for (int j=1;j<=4;j++) {
				if (var[i][j]==1) {
					count_s++;
				}
				if (var[i][j]==2) {
					count_ns++;
				}
			}
			if (count_ns>0) {
				for (int j=0;j<all_nucs.size();j++) {
					vector<int> temp=all_nucs[j][i];
					//Remove synonymous counts (not consensus though)
					for (int k=1;k<=4;k++) {
						if (var[i][k]==1) {
							temp[k]=0;
						}
					}
					//Push this line back to all_nucs_ns
					all_nucs_ns[j].push_back(temp);
				}
			}
			if (count_s>0) {
				for (int j=0;j<all_nucs.size();j++) {
					vector<int> temp=all_nucs[j][i];
					//Remove synonymous counts (not consensus though)
					for (int k=1;k<=4;k++) {
						if (var[i][k]==2) {
							temp[k]=0;
						}
					}
					//Push this line back to all_nucs_ns
					all_nucs_s[j].push_back(temp);
				}
			}
		}
	}
}

void PiCalculation(vector< vector< vector<int> > >& all_nucs, vector<double>& pi) {
	for (int i=0;i<all_nucs.size();i++) {
		vector<double> div;
		double p=0;
		double count=0;
		for (int j=0;j<all_nucs[i].size();j++) {
			double n=all_nucs[i][j][4];
			if (n>1) {
				double d=n*(n-1);
				for (int k=0;k<=3;k++) {
					d=d-(all_nucs[i][j][k]*(all_nucs[i][j][k]-1));
				}
				d=d/(n*(n-1));
				div.push_back(d);
				p=p+d;
				count++;
			}
		}
		p=p/all_nucs[i].size();
		pi.push_back(p);
	}
}


void Calc_Variant_Composition (run_params p, vector< vector< vector<int> > >& all_nucs, vector< vector<char> >& all_cons, vector< vector< vector< vector<int> > > >& bootstrap) {
	ofstream cv_file;
	cv_file.open("Variant_composition.dat");
	CVC(p,all_nucs,all_cons,cv_file);
	cv_file.close();
	if (p.bootstrap==1) {
		cv_file.open("Variant_composition_bootstrap.dat");
		for (int i=0;i<bootstrap.size();i++) {
			CVC(p,bootstrap[i],all_cons,cv_file);
			cv_file << "\n";
		}
	}
}


void CVC (run_params p, vector< vector< vector<int> > >& all_nucs, vector< vector<char> >& all_cons, ofstream& cv_file) {
	vector< vector<double> > twelve_class;
	for (int i=0;i<all_nucs.size();i++) {
		vector<double> tc;
		for(int i=0;i<12;i++) {
			tc.push_back(0);
		}
		for (int j=0;j<all_nucs[i].size();j++) {
			double n=all_nucs[i][j][4];
			if (n>0) {
				double qa=(all_nucs[i][j][0]+0.)/(n+0.);
				double qc=(all_nucs[i][j][1]+0.)/(n+0.);
				double qg=(all_nucs[i][j][2]+0.)/(n+0.);
				double qt=(all_nucs[i][j][3]+0.)/(n+0.);
			
				if (qa>p.q_cut) {
					qa=0;
				}
				if (qc>p.q_cut) {
					qc=0;
				}
				if (qg>p.q_cut) {
					qg=0;
				}
				if (qt>p.q_cut) {
					qt=0;
				}
				if (all_cons[i][j]=='A') {
					tc[0]=tc[0]+qc;
					tc[1]=tc[1]+qg;
					tc[2]=tc[2]+qt;
				}
				if (all_cons[i][j]=='C') {
					tc[3]=tc[3]+qa;
					tc[4]=tc[4]+qg;
					tc[5]=tc[5]+qt;
				}
				if (all_cons[i][j]=='G') {
					tc[6]=tc[6]+qa;
					tc[7]=tc[7]+qc;
					tc[8]=tc[8]+qt;
				}
				if (all_cons[i][j]=='T') {
					tc[9]=tc[9]+qa;
					tc[10]=tc[10]+qc;
					tc[11]=tc[11]+qg;
				}
			}
		}
		double tot=0;
		for (int k=0;k<12;k++) {
			tot=tot+tc[k];
		}
		for (int k=0;k<12;k++) {
			tc[k]=tc[k]/tot;
			cv_file << tc[k] << " ";
		}
		
		cv_file << "\n";
		twelve_class.push_back(tc);
	}

}


void Calculate_Consensus_Hamming_Distances (vector<string>& sam_files) {
	vector<string> consensus;
	InputConsensusSequences(sam_files,consensus);
	vector< vector<int> > ham_dist;
	cout << consensus.size() << "\n";
	for (int i=0;i<consensus.size();i++) {
		vector<int> d;
		cout << "i= " << i << "\n";
		for (int j=0;j<consensus.size();j++) {
			cout << "j= "<< j << "\n";
			int h=0;
			for (int k=0;k<consensus[i].length();k++) {
				if (consensus[i].compare(k,1,"N")!=0&&consensus[j].compare(k,1,"N")!=0&&consensus[i].compare(k,1,consensus[j],k,1)!=0) {
					h++;
				}
			}
			d.push_back(h);
		}
		ham_dist.push_back(d);
	}
	ofstream ham_file;
	ham_file.open("Consensus_Hamming_Distances.dat");
	for (int i=0;i<ham_dist.size();i++) {
		for (int j=0;j<ham_dist[i].size();j++){
			ham_file << ham_dist[i][j] << " ";
		}
		ham_file << "\n";
	}
}

void Calculate_Variant_Hamming_Distances (vector< vector< vector<int> > >& all_nucs) {
	vector< vector<double> > ham_dist;
	for (int i=0;i<all_nucs.size();i++) {
		vector<double> vec_d;
		for (int j=0;j<all_nucs.size();j++) {
			double all_d=0;
			//cout << "i= " << i << " j= " << j << "\n";
			for (int k=0;k<all_nucs[i].size();k++) {
				vector<double> freqs_i;
				vector<double> freqs_j;
				for (int l=0;l<4;l++) {
					double q=-1;
					if (all_nucs[i][k][4]>0) {
						q=(all_nucs[i][k][l]+0.)/(all_nucs[i][k][4]+0.);
					}
					freqs_i.push_back(q);
				}
				for (int l=0;l<4;l++) {
					double q=-1;
					if (all_nucs[j][k][4]>0) {
						q=(all_nucs[j][k][l]+0.)/(all_nucs[j][k][4]+0.);
					}
					freqs_j.push_back(q);
				}
				double d=0;
				for (int l=0;l<4;l++) {
					if (freqs_i[l]>=0&&freqs_j[l]>=0) {
						d=d+abs(freqs_i[l]-freqs_j[l]);
					} else {
						d=-10;
					}
				}
				if (d>0) {
					d=d/2;
					all_d=all_d+d;
				}
			}
			//cout << "Sum distance " << all_d << "\n";
			vec_d.push_back(all_d);
		}
		ham_dist.push_back(vec_d);
	}
	ofstream ham_file;
	ham_file.open("Variant_Hamming_Distances.dat");
	for (int i=0;i<ham_dist.size();i++) {
		for (int j=0;j<ham_dist[i].size();j++){
			ham_file << ham_dist[i][j] << " ";
		}
		ham_file << "\n";
	}
}
	
void Calc_Tot_Freqs (run_params p, vector< vector< vector<int> > >& all_nucs, vector< vector< vector< vector<int> > > >& bootstrap) {
	cout << "Calculating total frequency statistics...\n";
	vector< vector<double> > total_freqs;
	ofstream tot_file;
	tot_file.open("Total_frequency_stats.dat");
	TF_Calculation(all_nucs,total_freqs,tot_file);
	tot_file.close();
	
	if (p.bootstrap==1) {
		tot_file.open("Total_frequency_bootstrap.dat");
		for (int i=0;i<bootstrap.size();i++) {
			TF_Calculation(bootstrap[i],total_freqs,tot_file);
			tot_file << "\n";
		}
	}
	
	if (p.sns>0) {
		cout << "Repeat for four-fold S and NS sites\n";

		//Import SNS types
		vector< vector<int> > var;
		InputVariantTypes(var);
		
		//Find four-fold variant sites
		vector< vector< vector<int> > > all_nucs_s;
		vector< vector< vector<int> > > all_nucs_ns;
		for (int i=0;i<all_nucs.size();i++) {
			vector< vector<int> > anucs;
			all_nucs_s.push_back(anucs);
			all_nucs_ns.push_back(anucs);
		}
		vector< vector<double> > total_freqs_s;
		vector< vector<double> > total_freqs_ns;
		if (p.sns==1) {
			Construct_FourFold_Variants(var,all_nucs,all_nucs_s,all_nucs_ns);
		} else {
			Construct_AllFreqsSNS_Variants(var,all_nucs,all_nucs_s,all_nucs_ns);
		}
		cout << "Found " << all_nucs_s[0].size() << " S sites and " << all_nucs_ns[0].size() << " NS sites\n";
		tot_file.open("Total_frequency_stats_S.dat");
		TF_Calculation(all_nucs_s,total_freqs_s,tot_file);
		tot_file.close();
		tot_file.open("Total_frequency_stats_NS.dat");
		TF_Calculation(all_nucs_ns,total_freqs_ns,tot_file);
		tot_file.close();
		
		if (p.bootstrap==1) {
			tot_file.open("Total_frequency_bootstrap_S.dat");
			ofstream tot2_file;
			tot2_file.open("Total_frequency_bootstrap_NS.dat");
			for (int i=0;i<bootstrap.size();i++) {
				all_nucs_s.clear();
				all_nucs_ns.clear();
				total_freqs_s.clear();
				total_freqs_ns.clear();
				for (int i=0;i<all_nucs.size();i++) {
					vector< vector<int> > anucs;
					all_nucs_s.push_back(anucs);
					all_nucs_ns.push_back(anucs);
				}
				if (p.sns==1) {
					Construct_FourFold_Variants(var,bootstrap[i],all_nucs_s,all_nucs_ns);
				} else {
					Construct_AllFreqsSNS_Variants(var,bootstrap[i],all_nucs_s,all_nucs_ns);
				}
				TF_Calculation(all_nucs_s,total_freqs_s,tot_file);
				tot_file << "\n";
				if (p.sns==1) {
					Construct_FourFold_Variants(var,bootstrap[i],all_nucs_s,all_nucs_ns);
				} else {
					Construct_AllFreqsSNS_Variants(var,bootstrap[i],all_nucs_s,all_nucs_ns);
				}
				TF_Calculation(all_nucs_ns,total_freqs_ns,tot2_file);
				tot2_file << "\n";
			}
			
		}
	}
}

void TF_Calculation (vector< vector< vector<int> > >& all_nucs, vector< vector<double> > total_freqs, ofstream& tot_file) {
	for (int i=0;i<all_nucs.size();i++) {
		vector<double> tot_freq (12,0);
		for (int j=0;j<all_nucs[i].size();j++) {
			double n=all_nucs[i][j][4];
			double max=1;
			if (n>0) {
				double qa=(all_nucs[i][j][0]+0.)/(n+0.);
				double qc=(all_nucs[i][j][1]+0.)/(n+0.);
				if (qc>qa) {max=2;}
				double qg=(all_nucs[i][j][2]+0.)/(n+0.);
				if (qg>qc&&qg>qa) {max=3;}
				double qt=(all_nucs[i][j][3]+0.)/(n+0.);
				if (qt>qg&&qt>qc&&qt>qa) {max=4;}
				if (max==1) {
					tot_freq[0]=tot_freq[0]+qc;
					tot_freq[1]=tot_freq[1]+qg;
					tot_freq[2]=tot_freq[2]+qt;
				}
				if (max==2) {
					tot_freq[3]=tot_freq[3]+qa;
					tot_freq[4]=tot_freq[4]+qg;
					tot_freq[5]=tot_freq[5]+qt;
				}
				if (max==3) {
					tot_freq[6]=tot_freq[6]+qa;
					tot_freq[7]=tot_freq[7]+qc;
					tot_freq[8]=tot_freq[8]+qt;
				}
				if (max==4) {
					tot_freq[9]=tot_freq[9]+qa;
					tot_freq[10]=tot_freq[10]+qc;
					tot_freq[11]=tot_freq[11]+qg;
				}
			}
		}
		total_freqs.push_back(tot_freq);
	}
	for (int i=0;i<total_freqs.size();i++) {
		double tot=0;
		for (int j=0;j<total_freqs[i].size();j++) {
			tot_file << total_freqs[i][j] << " ";
			tot=tot+total_freqs[i][j];
		}
		tot_file << tot << "\n";
	}
}
