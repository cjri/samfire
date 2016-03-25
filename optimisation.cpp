#include "shared.h"
#include "optimisation.h"
#include "likelihood.h"
#include <iostream>
#include <string>

void GetFactVectorSL (vector<str> sltrajs, vector<double>& fact_store) {
	int maxN=0;
	for (int i=0;i<sltrajs.size();i++) {
		
		for (int j=0;j<sltrajs[i].qA.size();j++) {
			if (maxN<sltrajs[i].nN[j]) {
				maxN=sltrajs[i].nN[j];
			}
		}
	}
	maxN=maxN*2;
	cout << maxN << "\n";
	FindLogFact(fact_store,maxN);
}

void GetFactVectorML (vector< vector<mtr> > mltrajs, vector<double>& fact_store) {
	int maxN=0;
	for (int i=0;i<mltrajs.size();i++) {
		for (int j=0;j<mltrajs[i].size();j++) {
			int tot=0;
			for (int k=0;k<mltrajs[i][j].n.size();k++) {
				tot=tot+mltrajs[i][j].n[k];
			}
			if (maxN<tot) {
				maxN=tot;
			}
		}
	}
	maxN=maxN+10;
	FindLogFact(fact_store,maxN);
}

void FindLogFact(vector<double>& fact_store,int N){
	double logN=0;
	fact_store.push_back(0);
	for (int i=1;i<=N;i++) {
		logN=logN+log(i);
		fact_store.push_back(logN);
	}
}


void OptimiseCsl (run_params p, double& Csl_opt, vector<str> sltrajs, vector<double> fact_store, gsl_rng *rgen) {
	cout << "Find optimal value of C_sl\n";
	Csl_opt=100;
	double C=100;
	double logL=-1e10;
	double logL_store=-1e10;
	double changec=10;
	int tryc=0;
	int acceptc=0;
	int max_it=100000;
	int first=1;
	vector<double> progress;
	
	for (int it=0;it<max_it;it++) {
		//cout << "Iteration " << it << "\n";
		if (it==max_it-1) {
			cout << "Warning: Possible convergence failure in C\n";
		}
		//cout << "Iteration " << it << "\n";
		if (first==0) {
			if (logL-logL_store>0) {
				Csl_opt=C;
				if (p.verb==1) {
					cout << "C " << C << " logL " << logL << "\n";
				}
				logL_store=logL;
				progress.push_back(logL);
				acceptc++;
			} else {
				C=Csl_opt;
			}
		}
		first=0;
		if (progress.size()>11&&progress[progress.size()-1]-progress[progress.size()-10]<0.001) {
			break;
		}
		
		if (it%100==0&&it>0) {
			double a_rate=(acceptc+0.)/(tryc+0.);
			changec=changec*(0.95+a_rate);
			acceptc=0;
			tryc=0;
		}
		
		C=C+(gsl_rng_uniform(rgen)*changec)-(changec/2);
		tryc++;
		if (it==max_it-1) {
			C=Csl_opt;
		}
		
		//Calculate likelihood
		logL=0;
		for (int i=0;i<sltrajs.size();i++) {
			if (sltrajs[i].inc==1) {
				vector<double> inf;
				vector<int> inc (4,0);
				if (sltrajs[i].mA>0) {
					inf.push_back(sltrajs[i].mA);
					inc[0]=1;
				}
				if (sltrajs[i].mC>0) {
					inf.push_back(sltrajs[i].mC);
					inc[1]=1;
				}
				if (sltrajs[i].mG>0) {
					inf.push_back(sltrajs[i].mG);
					inc[2]=1;
				}
				if (sltrajs[i].mT>0) {
					inf.push_back(sltrajs[i].mT);
					inc[3]=1;
				}
				for (int j=0;j<sltrajs[i].qA.size();j++) {
					vector<int> obs;
					if (inc[0]==1) {
						obs.push_back(sltrajs[i].nA[j]);
					}
					if (inc[1]==1) {
						obs.push_back(sltrajs[i].nC[j]);
					}
					if (inc[2]==1) {
						obs.push_back(sltrajs[i].nG[j]);
					}
					if (inc[3]==1) {
						obs.push_back(sltrajs[i].nT[j]);
					}
					int N=sltrajs[i].nN[j];
					logL=logL+DirichletMultiCalc(N,C,obs,inf,fact_store);
				}
			}
		}
	}
	
	cout << "Optimal C is " << Csl_opt << "\n";

}


void OptimiseCml (run_params p, double& Cml_opt, vector<int> include_ml, vector< vector<mtr> > mltrajs, vector<double> fact_store, gsl_rng *rgen) {
	cout << "Find optimal value of C_ml\n";
	Cml_opt=100;
	double C=100;
	double logL=-1e10;
	double logL_store=-1e10;
	double changec=10;
	int tryc=0;
	int acceptc=0;
	int max_it=100000;
	int first=1;
	vector<double> progress;
	
	for (int it=0;it<max_it;it++) {
		if (it==max_it-1) {
			cout << "Warning: No convergence in C\n";
		}
		if (first==0) {
			if (logL-logL_store>0) {
				Cml_opt=C;
				if (p.verb==1) {
					cout << "C " << C << " logL " << logL << "\n";
				}
				logL_store=logL;
				progress.push_back(logL);
				acceptc++;
			} else {
				C=Cml_opt;
			}
		}
		first=0;
		if (progress.size()>11&&progress[progress.size()-1]-progress[progress.size()-10]<0.001) {
			break;
		}
		
		if (it%100==0&&it>0) {
			double a_rate=(acceptc+0.)/(tryc+0.);
			changec=changec*(0.95+a_rate);
			acceptc=0;
			tryc=0;
		}
		
		C=C+(gsl_rng_uniform(rgen)*changec)-(changec/2);
		tryc++;
		if (it==max_it-1) {
			C=Cml_opt;
		}
		
		//Calculate likelihood
		logL=0;
		for (int i=0;i<mltrajs.size();i++) {  //Each partial haplotype set
			
			if (include_ml[i]==1) {
				
				//Inferred i.e. mean values
				vector<double> inf;
				for (int j=0;j<mltrajs[i].size();j++) { //Each haplotype
					inf.push_back(mltrajs[i][j].m);
				}
				
				//Observed values for each time point
				for (int k=0;k<mltrajs[i][0].n.size();k++) {  //Each time point
					vector<int> obs;
					int N=0;
					for (int j=0;j<mltrajs[i].size();j++) { //Each haplotype
						obs.push_back(mltrajs[i][j].n[k]);
						N=N+mltrajs[i][j].n[k];
					}
					logL=logL+DirichletMultiCalc(N,C,obs,inf,fact_store);
				}
				
			//	cout << i << " " << logL << "\n";

			}
			
		}
		//cout << "Log L " << logL << "\n";
		
	}
	
	cout << "Optimal C is " << Cml_opt << "\n";
	
}




void ConservativeCml (run_params p, double& Cml_opt, vector<int> include_ml, vector< vector<mtr> > mltrajs, vector<double> fact_store, gsl_rng *rgen) {
	cout << "Find conservative value of C_ml:\n";
	double C=Cml_opt;
	int max_it=100000;
	double logL=0;
	double logL_target=0;
	double changec=10;
	int first=1;
	vector<double> logs;
	for (int it=0;it<max_it;it++) {
		//cout << "Iteration " << it << "\n";
		if (first==0) {
			logs.push_back(logL);
			double diff=logL-logL_target;
			if (diff<0) {diff=-diff;}
			//cout << "Diff " << diff << "\n";
			if (diff<0.00001) {
				cout << "Conservative C = " << C << "\n";
				Cml_opt=C;
				break;
			} else if (logL>logL_target) {
				//cout << logL << " " << logs[it-2] << "\n";
				if (it>=2&&logs[it-2]<logL_target) {
					changec=changec*0.9;
				}
				C=C-(gsl_rng_uniform(rgen)*changec);
			} else if (C<0) {
				cout << "Conservative C = " << C << "\n";
				Cml_opt=C;
				break;
			} else {
				if (it>=2&&logs[it-2]>logL_target) {
					changec=changec*0.9;
				}
				C=C+(gsl_rng_uniform(rgen)*changec);
			}
			//cout << C << " " << logL << "\n";
		}
		
		//Gradient optimisation here - one dimension...
	
		
		//Calculate likelihood
		logL=0;
		for (int i=0;i<mltrajs.size();i++) {  //Each partial haplotype set
			
			if (include_ml[i]==1) {
				
				//Inferred i.e. mean values
				vector<double> inf;
				for (int j=0;j<mltrajs[i].size();j++) { //Each haplotype
					inf.push_back(mltrajs[i][j].m);
				}
				
				//Observed values for each time point
				for (int k=0;k<mltrajs[i][0].n.size();k++) {  //Each time point
					vector<int> obs;
					int N=0;
					for (int j=0;j<mltrajs[i].size();j++) { //Each haplotype
						obs.push_back(mltrajs[i][j].n[k]);
						N=N+mltrajs[i][j].n[k];
					}
					logL=logL+DirichletMultiCalc(N,C,obs,inf,fact_store);
				}
				
				//	cout << i << " " << logL << "\n";
				
			}
			
		}
		if (first==1) {
			logL_target=logL-5;
			first=0;
		}
		
	}
}
	
void GetInfMeanFreq(int traj, vector<str> sltrajs, vector<double>&init_freqs) {
	init_freqs.clear();
	init_freqs.push_back(sltrajs[traj].mA);
	init_freqs.push_back(sltrajs[traj].mC);
	init_freqs.push_back(sltrajs[traj].mG);
	init_freqs.push_back(sltrajs[traj].mT);
}


void GetInfStartFreq(int traj, vector<str> sltrajs, vector<double>&init_freqs) {
	init_freqs.clear();
	init_freqs.push_back(sltrajs[traj].qA[0]);
	init_freqs.push_back(sltrajs[traj].qC[0]);
	init_freqs.push_back(sltrajs[traj].qG[0]);
	init_freqs.push_back(sltrajs[traj].qT[0]);
}


void GetDetInfTDVector(int traj, vector<str> sltrajs, vector< vector<double> >& inf) {
	for (int i=0;i<sltrajs[traj].nA.size();i++) {
		vector<double> in;
		in.push_back(sltrajs[traj].qA[i]);
		in.push_back(sltrajs[traj].qC[i]);
		in.push_back(sltrajs[traj].qG[i]);
		in.push_back(sltrajs[traj].qT[i]);
		inf.push_back(in);
	}
}


double OptimiseSLTraj (int verb, run_params p, int maxrep, int traj, char sel_model, int tds, int n_times, vector<int> times, vector<int> N, double Csl_opt, vector<str> sltrajs, vector<double> fact_store, gsl_rng *rgen) {
	
	int n_haps=4;
	vector< vector<int> > obs;
	GetObsVector(traj,sltrajs,obs);
	vector< vector<double> > inf;

	double mu=p.mu;
	vector <vector <double> > mut;
	if (p.det==0) {
		GetMutationMatrix(mu,mut);
		SquareMutationMatrix(mut);
	}
	
	//Initial frequencies
	vector<double> init_freqs;
	SetInitialFreqs (init_freqs,n_haps,rgen);

	vector<double> sigs;
	SetInitialSelection(p,tds,sel_model,n_times,sigs,rgen);
	
	vector< vector<double> > hapsigs;
	SigsHapSigs(tds,sel_model,n_times,sigs,hapsigs);

	double L=-1e-10;
	double storeL=-1e10;
	double bestL=-1e10;
	double dL=-1;
	int acceptx=0;
	int accepts=0;
	int trys=0;
	int tryx=0;
	int max_it=1000000;
	double changex=0.5;
	double changes=0.1;
	double c=Csl_opt;
	int first=1;
	int movex=0;
	double a_rate=0;
	vector<double> init_freqs_store;
	vector<double> init_freqs_best;
	vector<double> sigs_store;
	vector<double> sigs_best;
	vector<double> progress;
	
	for (int rep=0;rep<maxrep;rep++) {
		if (p.verb==1) {
			cout << "Iteration # " << rep << "\n";
		}
		changex=0.5;
		changes=0.1;
		acceptx=0;
		tryx=0;
		storeL=-1e10;
		progress.clear();
	
		for (int it=0;it<max_it;it++) {
			if (first==0) {
				//Evaluate likelihood
				dL=L-storeL;
				if (dL>0) {
					storeL=L;
					init_freqs_store=init_freqs;
					sigs_store=sigs;
					
					if (L>bestL) {
						bestL=L;
						init_freqs_best=init_freqs;
						sigs_best=sigs;
					}

					if (movex==1) {
						acceptx++;
					}
					if (movex==0) {
						accepts++;
					}
					progress.push_back(L);
					if (progress.size()>25&&progress[progress.size()-21]-progress[progress.size()-1]>-0.001) {
						if (verb==1) {
							if (p.det==1) {
								DeterministicPropagation(1,times,init_freqs,hapsigs,inf);
							} else {
								Propagation(1,times,init_freqs,mut,hapsigs,inf);
							}
							cout << "Inferred freqs\n";
							for (int i=0;i<init_freqs.size();i++) {
								cout << init_freqs[i] << " ";
							}
							cout << "\n";
							cout << "Inferred sigma\n";
							for	(int i=0;i<sigs.size();i++) {
								cout << sigs[i] << "\n";
							}
							cout << "Optimal L = " << L << "\n";
						}
						break;
					}
				} else {
					init_freqs=init_freqs_store;
					sigs=sigs_store;
				}
			}
			first=0;
			//Modify parameter change variables
			if (it>10000) {
				if (it%1000==0) {
					a_rate=(acceptx+0.)/(tryx+0.);
					changex=changex*(0.95+a_rate);
					if (trys>0) {
						a_rate=(accepts+0.)/(trys+0.);
						changes=changes*(0.95+a_rate);
					}
					acceptx=0;
					tryx=0;
					accepts=0;
					trys=0;
				}
			}
		
			//Initial frequency
			if (it%2==0) {
				int j=floor(gsl_rng_uniform(rgen)*4);
				init_freqs[j]=init_freqs[j]+(gsl_rng_uniform(rgen)*changex)-(changex/2);
				if(init_freqs[j]<0) {
					init_freqs[j]=1e-20;
				}
				NormaliseFreqs (init_freqs,n_haps);
				tryx++;
				movex=1;
			}
		
			//Selection coefficients
			if (it%2==1) {
				for (int j=0;j<sigs.size();j++) {
					sigs[j]=sigs[j]+(gsl_rng_uniform(rgen)*changes)-(changes/2);
				}
				trys++;
				movex=0;
			}
		
			//Convert sigmas to haplotype sigmas
			hapsigs.clear();
			SigsHapSigs (tds,sel_model,n_times,sigs,hapsigs);
		
			//Propagate under mutation/selection
			if (p.det==1) {
				DeterministicPropagation(0,times,init_freqs,hapsigs,inf);
			} else {
				Propagation(0,times,init_freqs,mut,hapsigs,inf);
			}
		
			//Evaluate likelihood
			L=0;
			for (int i=0;i<times.size();i++) {
				double lL=DirichletMultiCalc(N[i],c,obs[i],inf[i],fact_store);
				L=L+lL;
			}
		}
	
		if (verb==1) {
			sigs=sigs_store;
			init_freqs=init_freqs_store;
			SigsHapSigs (tds,sel_model,n_times,sigs,hapsigs);
			if (p.det==1) {
				DeterministicPropagation(1,times,init_freqs,hapsigs,inf);
			} else {
				Propagation(1,times,init_freqs,mut,hapsigs,inf);
			}
			cout << "Inferred freqs\n";
			for (int i=0;i<init_freqs.size();i++) {
				cout << init_freqs[i] << " ";
			}
			cout << "\n";
			cout << "Inferred sigma\n";
			for	(int i=0;i<sigs.size();i++) {
				cout << sigs[i] << "\n";
			}
			cout << "Optimal L = " << storeL << "\n";
		}

		init_freqs=init_freqs_best;
		sigs=sigs_best;
	}
	return bestL;
}

void GetObsVector(int traj, vector<str> sltrajs, vector< vector<int> >& obs) {
	for (int i=0;i<sltrajs[traj].nA.size();i++) {
		vector<int> o;
		o.push_back(sltrajs[traj].nA[i]);
		o.push_back(sltrajs[traj].nC[i]);
		o.push_back(sltrajs[traj].nG[i]);
		o.push_back(sltrajs[traj].nT[i]);
		obs.push_back(o);
	}
}

void GetMutationMatrix (double mu, vector <vector <double> >& mut) {
	for (int i=0;i<4;i++) {
		vector<double> row;
		for (int j=0;j<4;j++) {
			if (i==j) {
				row.push_back(1-(3*mu));
			} else {
				row.push_back(mu);
			}
		}
		mut.push_back(row);
	}
}

void SquareMutationMatrix (vector<vector<double> > &m ) {
	vector<vector<double> > m2;
	for (unsigned int i=0;i<m.size();i++) {
		vector<double> m2i;
		m2i.assign(m.size(),0);
		m2.push_back(m2i);
	}
	for (unsigned int i=0;i<m.size();i++) {
		for (unsigned int j=0;j<m[i].size();j++) {
			for (unsigned int k=0;k<m[i].size();k++) {
				m2[i][j]=m2[i][j]+(m[i][k]*m[k][j]);
			}
		}
	}
	m=m2;
}

void SetInitialFreqs (vector<double>& init_freqs, int n_haps, gsl_rng *rgen) {
	for (int i=0;i<n_haps;i++) {
		double x=gsl_rng_uniform(rgen);
		init_freqs.push_back(x);
	}
	NormaliseFreqs(init_freqs,n_haps);
//	cout << "Initial freqs: ";
//	for (int i=0;i<n_haps;i++) {
//		cout << init_freqs[i] << " ";
//	}
//	cout << "\n";
}

void NormaliseFreqs (vector<double>& init_freqs, int n_haps) {
	double tot=0;
	for (int i=0;i<n_haps;i++) {
		tot=tot+init_freqs[i];
	}
	for (int i=0;i<n_haps;i++) {
		init_freqs[i]=init_freqs[i]/tot;
	}
}

void SetInitialSelection (run_params p, int tds, char sel_model, int n_times, vector<double>& sigs, gsl_rng *rgen) {
	if (sel_model!='0') {
		if (tds==0) {
			double s=0;
			if (p.det==1) {
				s=(gsl_rng_uniform(rgen)*0.05);
			} else {
				s=gsl_rng_uniform(rgen)+0.05;
			}
			sigs.push_back(s);
		} else {
			for (int t=0;t<n_times-1;t++) {
				double s=0;
				if (p.det==1) {
					s=(gsl_rng_uniform(rgen)*0.1);
				} else {
					s=gsl_rng_uniform(rgen)-0.5;
				}
				sigs.push_back(s);
			}
		}
	}
}

void SigsHapSigs (int tds, char sel_model, int n_times, vector<double> sigs, vector< vector<double> >& hapsigs) {
	char nucs[] = {'A','C','G','T'};
	for (int h=0;h<=4;h++) { //Haplotypes A C G T
		vector<double> hs;
		if (sel_model==nucs[h]) {
			if (tds==1) {
				for (int j=0;j<n_times-1;j++) {
					hs.push_back(sigs[j]);
				}
			} else {
				for (int j=0;j<n_times-1;j++) {
					hs.push_back(sigs[0]);
				}
			}
		} else {
			for (int j=0;j<n_times-1;j++) {
				hs.push_back(0);
			}
		}
		hapsigs.push_back(hs);
	}
}


void Propagation(int verb, vector<int> times, vector<double> init_freqs, vector< vector<double> > mut, vector< vector<double> > hapsigs, vector< vector<double> >& inf) {
	inf.clear();
	int start=times[0];
	int fin=times[times.size()-1];
	int index=0;
	vector<double> si;
	si=hapsigs[index];
	GetSI(index,hapsigs,si);
	vector<double> q=init_freqs;
	for (int j=0;j<=fin;j++) {
		if (j==0&&j==times[index]) {
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

void DeterministicPropagation(int verb, vector<int> times, vector<double> init_freqs, vector< vector<double> > hapsigs, vector< vector<double> >& inf) {
	inf.clear();
	vector<double> si;
	si=hapsigs[0];
	vector<double> q=init_freqs;
	if (verb==1) {
		for (int i=0;i<q.size();i++) {
			cout << q[i] << " ";
		}
		cout << "\n";
	}
	inf.push_back(q);
	int index=0;
	if (times.size()>1) {
		for (int j=1;j<times.size();j++) {
			GetSI(index,hapsigs,si);
			int dt=times[j]-times[j-1];
			for (int i=0;i<si.size();i++) {
				si[i]=si[i]*dt;
			}
			//Note here - selection acts twice per day in the with-mutation model; where selection is applied per day the inferred value is doubled...
			GrowSig(q,si);
			if (verb==1) {
				for (int i=0;i<q.size();i++) {
					cout << q[i] << " ";
				}
				cout << "\n";
			}
			inf.push_back(q);
			index++;
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
	MatMult(m,q,temp);  //Using m^2 matrix...
	q=temp;
	GrowSig(q,sigma);
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


void OptimiseLD (run_params p, double C, ld_info ld, int t, vector<double>& fact_store, ofstream& ld_file, gsl_rng *rgen) {
	//cout << "Find optimal value of LD\n";
	int n_haps=4;
	vector<double> q;
	SetInitialFreqs (q,n_haps,rgen);
	vector<double> q_store=q;
	double logL=-1e10;
	double logL_store=-1e10;
	double changex=0.05;
	int tryc=0;
	int acceptc=0;
	int max_it=100000;
	int first=1;
	vector<double> progress;
	double Dij=0;
	double D=0;
	int N=0;
	N=N+ld.n_11[t];
	N=N+ld.n_10[t];
	N=N+ld.n_01[t];
	N=N+ld.n_00[t];
	if (N>0) {
		for (int it=0;it<max_it;it++) {
			if (it==max_it-1) {
				//cout << "Warning: No convergence in LD\n";
			}
			if (first==0) {
				if (logL-logL_store>0) {
					q_store=q;
					if (p.verb==1) {
					//	cout << "q " << q[0] << " " << q[1] << " " << q[2] << " " << q[3] << " LD " << q[0]-((q[0]+q[1])*(q[2]+q[3])) << " logL " << logL << "\n";
					}
					logL_store=logL;
					progress.push_back(logL);
					acceptc++;
				} else {
					q=q_store;
				}
			}
			first=0;
			if (progress.size()>21&&progress[progress.size()-1]-progress[progress.size()-20]<0.001) {
				break;
			}
		
			if (it%100==0&&it>0) {
				double a_rate=(acceptc+0.)/(tryc+0.);
				changex=changex*(0.95+a_rate);
				acceptc=0;
				tryc=0;
			}
		
			int j=floor(gsl_rng_uniform(rgen)*n_haps);
			q[j]=q[j]+(gsl_rng_uniform(rgen)*changex)-(changex/2);
			if(q[j]<0) {
				q[j]=1e-20;
			}
			NormaliseFreqs (q,n_haps);
			tryc++;
			if (it==max_it-1) {
				q=q_store;
			}
		
			//Calculate likelihood
			logL=0;
			//Inferred values are the q^ab_ij
			vector<double> inf=q;
			vector<int> obs;
			int N=0;
			obs.push_back(ld.n_11[t]);
			obs.push_back(ld.n_10[t]);
			obs.push_back(ld.n_01[t]);
			obs.push_back(ld.n_00[t]);
			N=N+ld.n_11[t];
			N=N+ld.n_10[t];
			N=N+ld.n_01[t];
			N=N+ld.n_00[t];
			logL=logL+DirichletMultiCalc(N,C,obs,inf,fact_store);
			
			//Single-locus data i
			inf.clear();
			obs.clear();
			N=0;
			obs.push_back(ld.n_i1[t]);
			obs.push_back(ld.n_i0[t]);
			N=N+ld.n_i1[t];
			N=N+ld.n_i0[t];
			inf.push_back(q[0]+q[1]);//q^1_i
			inf.push_back(q[2]+q[3]);//q^0_i
			logL=logL+DirichletMultiCalc(N,C,obs,inf,fact_store);
			
			//Single-locus data j
			inf.clear();
			obs.clear();
			N=0;
			obs.push_back(ld.n_j1[t]);
			obs.push_back(ld.n_j0[t]);
			N=N+ld.n_j1[t];
			N=N+ld.n_j0[t];
			inf.push_back(q[0]+q[2]);//q^1_j
			inf.push_back(q[1]+q[3]);//q^0_j
			logL=logL+DirichletMultiCalc(N,C,obs,inf,fact_store);
		
			//cout << "Log L " << logL << "\n";
		
		}
		Dij=q[0]-((q[0]+q[1])*(q[0]+q[2]));
		double mx=0;
		if (Dij>0) {
			double a=(q[0]+q[1])*(q[1]+q[3]);
			double b=(q[2]+q[3])*(q[0]+q[2]);
			if (a<b) {
				mx=a;
			} else {
				mx=b;
			}
		} else {
			double a=(q[0]+q[1])*(q[0]+q[2]);
			double b=(q[2]+q[3])*(q[1]+q[3]);
			if (a<b) {
				mx=a;
			} else {
				mx=b;
			}
		}
		
		if (mx!=0) {
			D=Dij/mx;
		}
		if (pow(Dij,2)<1e-20) {
			Dij=0;
			D=0;
		}
	}
	ld_file << "Opt " << t << " " << ld.i << " " << ld.j << " " << Dij << " " << D << "\n";
	if (p.verb==1) {
		cout << "Opt " << t << " " << ld.i << " " << ld.j << " " << Dij << " " << D << "\n";
	}
	
}


