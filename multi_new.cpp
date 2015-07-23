//Applies linear model to unobserved haplotypes.  Retains loci for which there is no variation.

#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

#include "shared2.h"

vector<rec> dat;  //Here, sequence data consists of multiple independent reads

vector<rec> dat1;  //Sequence data consists of multiple independent reads
vector<rec> dat2;  //Sequence data consists of multiple independent reads
vector<rec> dat3;  //Sequence data consists of multiple independent reads
vector<rec> dat4;  //Sequence data consists of multiple independent reads
vector<rec> dat5;  //Sequence data consists of multiple independent reads

rec data;
int main(int argc, const char **argv) {

	run_params p;
	GetOptions (p,argc,argv);
	if (p.mcmc==1) {
		p.read=1;  //MCMC step does MCMC around previous optimum
	}
	
	//Initialise random number generator
	gsl_rng_env_setup();
	gsl_rng *rgen = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set (rgen, p.seed);
	
	//Observation times
	vector<int> times;
	ObsTimes(p,times);
	
	vector<double> fact_store;
	FindLogFact(fact_store,300000);
	
	ifstream in_file;
	in_file.open(p.data);
	InputDataHap(times,in_file,dat1);

	cout << "Times\n";
	for (int i=0;i<times.size();i++) {
		cout << times[i] << " ";
	}
	cout << "\n";
	
	cout << dat1.size() << "\n";

	vector< vector<int> > obs;
	vector<int> N;
	GetObsStats (times,dat1,N,obs);

	
	cout << "Read information\n";
	for (int i=0;i<dat1.size();i++) {
		cout << i << " ";
		for (int j=0;j<dat1[i].seq.size();j++) {
			cout << dat1[i].seq[j] << " ";
		}
		cout << "\n";
		for (int j=0;j<dat1[i].obs.size();j++) {
			cout << dat1[i].obs[j] << " ";
		}
		cout << "\n";
	}
	
	//Calculate consensus and base types.  Look for loci with no variation
	vector<int> consensus;
	vector<int> bases;
	GetConsensusBaseAllHaps (consensus,bases,dat1);

	cout << "Edited information\n";
	for (int i=0;i<dat1.size();i++) {
		cout << i << " ";
		for (int j=0;j<dat1[i].seq.size();j++) {
			cout << dat1[i].seq[j] << " ";
		}
		cout << "\n";
		for (int j=0;j<dat1[i].obs.size();j++) {
			cout << dat1[i].obs[j] << " ";
		}
		cout << "\n";
	}
	cout << "Consensus\n";
	for (int i=0;i<consensus.size();i++) {
		cout << consensus[i] << " ";
	}
	cout << "\n";
	cout << "Bases\n";
	for (int i=0;i<bases.size();i++) {
		cout << bases[i] << " ";
	}
	cout << "\n";
	
	int dim=consensus.size();
	

	//Square mutation matrix here...
	
	//Idea - selection is nucleotide-specific. Don't need consensus and base in this respect...
	
	vector<char> sel_model; //Code describing selection or non-selection on each locus e.g. A C 0 0.  0: Neutral, A: Constant selection for allele A
	vector<int> tds;
	vector<epi> epistat2;
	vector<epi> epistat3;
	vector<epi> epistat4;
	vector<epi> epistat5;
	cout << "Read selection model\n";
	ReadSelectionModel(dim,sel_model,tds,epistat2,epistat3,epistat4,epistat5);

	//At this point, divide unobserved space into regions
	vector<seldiv> sd;
	
	if (p.unobs==0) {
		seldiv s;
		s.seq.push_back('0');
		sd.push_back(s);
		for (int i=0;i<sel_model.size();i++) {
			if (sel_model[i]!='0') {
				vector<seldiv> newsd;
				for (int j=0;j<sd.size();j++) {
					seldiv s1=sd[j];
					seldiv s2=sd[j];
					s1.seq.push_back(sel_model[i]);
					if (sel_model[i]=='A') {
						s2.seq.push_back('1');
					}
					if (sel_model[i]=='C') {
						s2.seq.push_back('2');
					}
					if (sel_model[i]=='G') {
						s2.seq.push_back('3');
					}
					if (sel_model[i]=='T') {
						s2.seq.push_back('4');
					}
					newsd.push_back(s1);
					newsd.push_back(s2);
				}
				sd=newsd;
			} else {
				for (int j=0;j<sd.size();j++) {
					sd[j].seq.push_back('0');
				}
			}
		}
	
		cout << "Seldiv\n";
		for (int i=0;i<sd.size();i++) {
			sd[i].seq.erase(sd[i].seq.begin());
			for (int j=0;j<sd[i].seq.size();j++) {
				cout << sd[i].seq[j] << " ";
			}
			cout << "\n";
		}
	}
	
	
	//Construct mutational information.  First check distances
	double mu=p.mu;
	vector< vector<double> > mut;
	if (p.unobs==0) {
		MakeMutationMatrixLinear (dim,mu,mut,sd,dat1);
	} else {
		MakeMutationMatrix (dim,mu,mut,dat1);
	}
	
	SquareMutationMatrix(mut);
	
	//Setup initial data for selection model and initial frequency model
	cout << "Initial Frequencies\n";
	vector<double> init_freqs;
	if (p.unobs==0) {
		SetInitialFreqsLinear (init_freqs,mut,dat1,rgen);
	} else {
		SetInitialFreqs (init_freqs,mut,dat1,rgen);
	}
	
	if (p.read==1) {
		ifstream ifreqs_file;
		ifreqs_file.open("Initial_freqs.in");
		for (int i=0;i<100;i++) {
			double x=0;
			if (!(ifreqs_file >> x)) break;
			init_freqs[i]=x;
		}
		ifreqs_file.close();
		cout << "Initial freqs read\n";
		for (int i=0;i<init_freqs.size();i++) {
			cout << init_freqs[i] << " ";
		}
		cout << "\n";
		ifreqs_file.close();
	}
	//SetInitialFreqs2 (init_freqs,mut,obs[0],N);

	cout << "Initial selection coefficients\n";
	vector< vector<double> > sigs;
	SetInitialSelection (dim,tds,sel_model,times,sigs,rgen);
	if (p.read==1) {
		ifstream sc_file;
		sc_file.open("Sel_coeffs.in");
		cout << "Selection coefficients read\n";
	
		for (int i=0;i<sigs.size();i++) {
			double x=0;
			for (int j=0;j<sigs[i].size();j++) {
				if (!(sc_file >> x)) break;
				sigs[i][j]=x;
				cout << x << "\n";
			}
		}
		sc_file.close();
	}

	SetInitialEpistasis (epistat2,epistat3,epistat4,epistat5,rgen);
	if (p.read==1) {
		ifstream ep_file;
		ep_file.open("Epi_stats.in");
		for (int i=0;i<epistat2.size();i++) {
			double x=0;
			if (!(ep_file >> x)) break;
			epistat2[i].x=x;
		}
		ep_file.close();

		ep_file.open("Epi_stats3.in");
		for (int i=0;i<epistat3.size();i++) {
			double x=0;
			if (!(ep_file >> x)) break;
			epistat3[i].x=x;
		}
		ep_file.close();
		
		ep_file.open("Epi_stats4.in");
		for (int i=0;i<epistat4.size();i++) {
			double x=0;
			if (!(ep_file >> x)) break;
			epistat4[i].x=x;
		}
		ep_file.close();

		ep_file.open("Epi_stats5.in");
		for (int i=0;i<epistat5.size();i++) {
			double x=0;
			if (!(ep_file >> x)) break;
			epistat5[i].x=x;
		}
		ep_file.close();
	}
	
	//Check epistasis - are there any cases where this is relevant?
	int go=1;
	CheckEpistasis(go,dim,sel_model,dat1,epistat2,epistat3,epistat4,epistat5);
	if (go==0) {
		return 0;
	}
	vector< vector<double> > hapsigs;
	if (p.unobs==0) {
		SigsHapSigsLinear (1,dim,times,sel_model,dat1,sigs,sd,epistat2,epistat3,epistat4,epistat5,hapsigs);
	} else {
		SigsHapSigs (1,dim,p,times,sel_model,dat1,sigs,epistat2,epistat3,epistat4,epistat5,hapsigs);
	}

	cout << "Here\n";
	
	//Get selection coefficients for unobserved states
	
	//Do evaluation
	int first=1;
	vector< vector<double> > inf;
	double L=-1e-10;
	double bestL=-1e10;
	double optL=-1e10;
	double storeL=-1e10;
	double dL=-1;
	vector<double> init_freqs_best;
	vector< vector<double> > sigs_best;
	vector<epi> epistat2_best;
	vector<epi> epistat3_best;
	vector<epi> epistat4_best;
	vector<epi> epistat5_best;
	vector<double> init_freqs_store;
	vector< vector<double> > sigs_store;
	vector<epi> epistat2_store;
	vector<epi> epistat3_store;
	vector<epi> epistat4_store;
	vector<epi> epistat5_store;
	vector<double> progress;
	double beta=1;
	double changex=0.05;
	double changes=0.1;
	double c=p.c;
	int acceptx=0;
	int accepts=0;
	int trys=0;
	int tryx=0;
	int movex=-1;
	int maxrep=1;
	if (p.read==1) {
		maxrep=5;
	}
	if (p.mcmc==1) {
		maxrep=1;
	}
	double a_rate=0;
	
	//Check selection exists
	vector<int> sel_alls;
	for (int i=0;i<sigs.size();i++) {
		if (sigs[i].size()>0) {
			sel_alls.push_back(i);
			cout << "Selection at " << i << "\n";
		}
	}

	for (int rep=0;rep<maxrep;rep++) {
		changex=0.05;
		changes=0.1;
		acceptx=0;
		tryx=0;
		for (int it=0;it<1000000;it++) {
			//cout << it << "\n";
			if (first==0) {
				//Evaluate likelihood
				if (p.mcmc==1) {
					dL=L-optL;
				} else {
					dL=L-storeL;
				}
				if (p.mcmc==0) {
					if (dL>0) {
						cout << "New freqs\n";
						for (int i=0;i<init_freqs.size();i++) {
							cout << init_freqs[i] << " ";
						}
						cout << "\n";
						cout << "New sigma\n";
						for	(int i=0;i<sigs.size();i++) {
							cout << i << " ";
							for (int j=0;j<sigs[i].size();j++) {
								cout << sigs[i][j] << " ";
							}
							cout << "\n";
						}
						for (int i=0;i<epistat2.size();i++) {
							cout << "Epistasis2 " << epistat2[i].loc[0] << " " << epistat2[i].loc[1] << " " << epistat2[i].x << "\n";
						}
						for (int i=0;i<epistat3.size();i++) {
							cout << "Epistasis3 " << epistat3[i].loc[0] << " " << epistat3[i].loc[1] << " "  << epistat3[i].loc[2] << " " << epistat3[i].x << "\n";
						}
						for (int i=0;i<epistat4.size();i++) {
							cout << "Epistasis4 " << epistat4[i].loc[0] << " " << epistat4[i].loc[1] << " "  << epistat4[i].loc[2] << " " << epistat4[i].loc[3] << " " << epistat4[i].x << "\n";
						}
						for (int i=0;i<epistat5.size();i++) {
							cout << "Epistasis5 " << epistat5[i].loc[0] << " " << epistat5[i].loc[1] << " "  << epistat5[i].loc[2] << " " << epistat5[i].loc[3] << " " << epistat5[i].loc[4] << " " << epistat5[i].x << "\n";
						}

						
						if (movex==1) {
							acceptx++;
						}
						if (movex==0) {
							accepts++;
						}
						//if (p.mcmc==0) {
							Propagation(1,times,init_freqs,mut,hapsigs,inf);
							cout << "Better L = " << L << "\n";
						//}
						
						//if (p.mcmc==0) {
							progress.push_back(L);
							if (progress.size()>200&&progress[progress.size()-21]-progress[progress.size()-1]>-0.00001) {
								break;
							}
						//}
						
						storeL=L;
						init_freqs_store=init_freqs;
						sigs_store=sigs;
						epistat2_store=epistat2;
						epistat3_store=epistat3;
						epistat4_store=epistat4;
						epistat5_store=epistat5;
						//Check for best so far
						if (L>bestL) {
							bestL=L;
							init_freqs_best=init_freqs;
							sigs_best=sigs;
							epistat2_best=epistat2;
							epistat3_best=epistat3;
							epistat4_best=epistat4;
							epistat5_best=epistat5;
						}
						
					} else {
						init_freqs=init_freqs_store;
						sigs=sigs_store;
						epistat2=epistat2_store;
						epistat3=epistat3_store;
						epistat4=epistat4_store;
						epistat5=epistat5_store;
					}
				}
				
				if (p.mcmc==1) {
					
					if (dL<-3) {
						if (movex==1) {
							acceptx++;
						}
						if (movex==0) {
							accepts++;
						}

						storeL=L;
						init_freqs_store=init_freqs;
						sigs_store=sigs;
						epistat2_store=epistat2;
						epistat3_store=epistat3;
						epistat4_store=epistat4;
						epistat5_store=epistat5;
						//Check for best so far
						if (L>bestL) {
							bestL=L;
							init_freqs_best=init_freqs;
							sigs_best=sigs;
							epistat2_best=epistat2;
							epistat3_best=epistat3;
							epistat4_best=epistat4;
							epistat5_best=epistat5;
						}
						if (it%100==0) {
							cout << it << "\n";
							cout << "Freqs ";
							for (int i=0;i<init_freqs_store.size();i++) {
								cout << init_freqs_store[i] << " ";
							}
							cout << "\n";
							cout << "Sigma ";
							for	(int i=0;i<sigs_store.size();i++) {
								for (int j=0;j<sigs_store[i].size();j++) {
									cout << sigs_store[i][j] << " ";
								}
							}
							cout << "\n";
							for (int i=0;i<epistat2.size();i++) {
								cout << "Epistasis2 " << epistat2[i].loc[0] << " " << epistat2[i].loc[1] << " " << epistat2[i].x << "\n";
							}
							for (int i=0;i<epistat3.size();i++) {
								cout << "Epistasis3 " << epistat3[i].loc[0] << " " << epistat3[i].loc[1] << " "  << epistat3[i].loc[2] << " " << epistat3[i].x << "\n";
							}
							for (int i=0;i<epistat4.size();i++) {
								cout << "Epistasis4 " << epistat4[i].loc[0] << " " << epistat4[i].loc[1] << " "  << epistat4[i].loc[2] << " " << epistat4[i].loc[3] << " " << epistat4[i].x << "\n";
							}
							for (int i=0;i<epistat5.size();i++) {
								cout << "Epistasis5 " << epistat5[i].loc[0] << " " << epistat5[i].loc[1] << " "  << epistat5[i].loc[2] << " " << epistat5[i].loc[3] << " " << epistat5[i].loc[4] << " " << epistat5[i].x << "\n";
							}
							cout << "Log L = " << storeL << "\n";
							
						}
						
					} else {
						init_freqs=init_freqs_store;
						sigs=sigs_store;
						epistat2=epistat2_store;
					}
				}
				

			
				//Modify parameter changes
				if (p.mcmc==0) {
					if (it>10000) {
						if (it%1000==0) {
							a_rate=(acceptx+0.)/(tryx+0.);
							if (a_rate<0.05) {
								changex=changex*(0.95+a_rate);
							} else {
								changex=changex*2*(a_rate+0.95);
							}
							if (trys>0) {
								a_rate=(accepts+0.)/(trys+0.);
								if (a_rate<0.05) {
									changes=changes*(0.95+a_rate);
								} else {
									changes=changes*2*(a_rate+0.95);
								}
							}
							acceptx=0;
							tryx=0;
							accepts=0;
							trys=0;
						}
					}
				}
					
				//Change parameters init_freqs and sigs
				if (p.mcmc==0) {
					if (it%2==0) {
						//Frequency
						int j=floor(gsl_rng_uniform(rgen)*(dat1.size()));
						init_freqs[j]=init_freqs[j]+(gsl_rng_uniform(rgen)*changex)-(changex/2);
						if(init_freqs[j]<0) {
							init_freqs[j]=1e-10;
						}
						NormaliseFreqs (init_freqs,mut);
						tryx++;
						movex=1;
					}
				}
			
				if (it%2==1) {
					//Selection
					if (sel_alls.size()>0) {
						int s=sel_alls.size()+epistat2.size()+epistat3.size()+epistat4.size()+epistat5.size();
						int r=floor(gsl_rng_uniform(rgen)*s);
					//	cout << r << " " << sel_alls.size() << " " << epistat2.size() << " "  << epistat3.size() << "\n";
						if (r<sel_alls.size()) {
					//		cout << "Here 1\n";
							int i=sel_alls[r];
							for (int j=0;j<sigs[i].size();j++) {
								sigs[i][j]=sigs[i][j]+(gsl_rng_uniform(rgen)*changes)-(changes/2);
							}
						} else if (r<sel_alls.size()+epistat2.size()) {
					//		cout << "Here 2\n";
							int i=r-sel_alls.size();
							epistat2[i].x=epistat2[i].x+(gsl_rng_uniform(rgen)*changes)-(changes/2);
							
						} else if (r<sel_alls.size()+epistat2.size()+epistat3.size()) {
					//		cout << "Here 3\n";
							int i=r-sel_alls.size()-epistat2.size();
							epistat3[i].x=epistat3[i].x+(gsl_rng_uniform(rgen)*changes)-(changes/2);
							
						} else if (r<sel_alls.size()+epistat2.size()+epistat3.size()+epistat4.size()) {
							int i=r-sel_alls.size()-epistat2.size()-epistat3.size();
							epistat4[i].x=epistat4[i].x+(gsl_rng_uniform(rgen)*changes)-(changes/2);
							
						} else {
							int i=r-sel_alls.size()-epistat2.size()-epistat3.size()-epistat4.size();
							epistat5[i].x=epistat5[i].x+(gsl_rng_uniform(rgen)*changes)-(changes/2);
						}
						trys++;
						movex=0;
					}
				}
			}
			first=0;
			//Convert sigmas to haplotype sigmas
			hapsigs.clear();
			if (p.unobs==0) {
				SigsHapSigsLinear (0,dim,times,sel_model,dat1,sigs,sd,epistat2,epistat3,epistat4,epistat5,hapsigs);
			} else {
				SigsHapSigs (0,dim,p,times,sel_model,dat1,sigs,epistat2,epistat3,epistat4,epistat5,hapsigs);
			}

			//Propagate under mutation/selection
			Propagation(0,times,init_freqs,mut,hapsigs,inf);

			//Evaluate likelihood
			L=0;
			for (int i=0;i<times.size();i++) {
				//double lL=MultiCalc(N[i],obs[i],inf[i],fact_store);
				double lL=DirichletMultiCalcLinear(N[i],c,obs[i],inf[i],sd,fact_store);
				L=L+lL;
//				cout << lL << " " << L << "\n";
			}
			if (p.mcmc==1&&it==1) {
				optL=L;
				cout << "Opt L = " << L << "\n";
			}
		}
		init_freqs=init_freqs_best;
		sigs=sigs_best;
		epistat2=epistat2_best;
		epistat3=epistat3_best;
		epistat4=epistat4_best;
		epistat5=epistat5_best;
	}

	return 0;
}
							
 
						
	
