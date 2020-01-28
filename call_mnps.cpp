#include "shared.h"
#include "call_mnps.h"
#include <iostream>
#include <string>

void MaxReadLength (run_params p, int& max_r_size, vector< vector<joined> >& t_reads) {
	for (int i=0;i<t_reads.size();i++) {
		for (int j=0;j<t_reads[i].size();j++) {
			if (t_reads[i][j].seq.size()>max_r_size) {
				max_r_size=t_reads[i][j].seq.size();
			}
		}
	}
	if (p.verb==1) {
		cout << "Maximum read length " << max_r_size << "\n";
	}
}

void MaxRdLength (run_params p, int& max_r_size,  vector<joined>& t_read) {
	for (int j=0;j<t_read.size();j++) {
		if (t_read[j].seq.size()>max_r_size) {
			max_r_size=t_read[j].seq.size();
		}
	}
	if (p.verb==1) {
		cout << "Maximum read length " << max_r_size << "\n";
	}
}


void MaxReadSpan(run_params p, int max_r_size, int& max_l, vector<int> polys) {
	int l=1;
	int hit=1;
	while (hit==1) {
		hit=0;
		for (int s=0;s<polys.size()-l;s++) {
			if (hit==0) {
				for (int f=s+l;f<polys.size();f++) {
					int dist=polys[f]-polys[s];
					//cout << polys[f] << " " << polys[s] << " " << dist << "\n";
					if (dist<max_r_size) {
						if (p.verb==1) {
							cout << "Span of length " << l << " " << polys[f] << " " << polys[s] << "\n";
						}
						hit=1;
						break;
					}
				}
			}
		}
		l++;
	}
	max_l=l-1;
	if (p.verb==1) {
		cout << "Max_l " << l-1 << "\n";
	}
}

void MapSequenceToSNPLoci (run_params p, int n_times, vector<int> polys, vector< vector<joined> > t_reads, vector <vector <vector<char> > >& t_m_vars) {
	if (p.verb==1) {
		cout << "Mapping to SNP loci\n";
	}
	for (int i=0;i<n_times;i++) {
		vector< vector<char> > m_vars;
		for (int j=0;j<t_reads[i].size();j++) {
			vector<char> nucs (polys.size(),'-');
			int start=t_reads[i][j].alpos;
			int fin=start+t_reads[i][j].seq.size();
			int inc=0;
			for (int k=0;k<polys.size();k++) {
				if (polys[k]>=start&&polys[k]<=fin) {
					int pos=polys[k]-start;
					char c=t_reads[i][j].seq[pos];
					if (c=='A'||c=='C'||c=='G'||c=='T') {
						nucs[k]=c;
						inc=1;
					}
				}
			}
			if (inc==1) {
				m_vars.push_back(nucs);
			}
		}
		t_m_vars.push_back(m_vars);
	}

}

void MapSeqToSNPLoci (run_params p, int n_times, vector<int> polys, vector<joined> t_read, vector <vector<char> >& m_vars) {
	if (p.verb==1) {
		cout << "Mapping to SNP loci\n";
	}
	for (int j=0;j<t_read.size();j++) {
		vector<char> nucs (polys.size(),'-');
		int start=t_read[j].alpos;
		int fin=start+t_read[j].seq.size();
		int inc=0;
		for (int k=0;k<polys.size();k++) {
			if (polys[k]>=start&&polys[k]<=fin) {
				int pos=polys[k]-start;
				char c=t_read[j].seq[pos];
				if (c=='A'||c=='C'||c=='G'||c=='T') {
					nucs[k]=c;
					inc=1;
				}
			}
		}
		if (inc==1) {
			m_vars.push_back(nucs);
		}
	}
}



void CallReadTypesConsecutive (run_params p, int max_l, int max_r_size, vector<int> polys, vector< vector<int> >& l_combs) {
	int s=polys.size();
	if (p.verb==1) {
		cout << "Potential sets of loci spanned by reads: here consecutive sets are considered\n";
	}
	for (int i=max_l;i>0;i--) {
		for (int j=0;j<s-i+1;j++) {
			vector <int> l_set;
			if (polys[j+i-1]-polys[j]<max_r_size) {
				for (int k=j;k<j+i;k++) {
					l_set.push_back(k);
				}
				if (l_set.size()>0) {
					l_combs.push_back(l_set);
				}
			}
		}
	}
	if (p.verb==1) {
		cout << "l_combs\n";
		for	(int i=0;i<l_combs.size();i++) {
			for (int j=0;j<l_combs[i].size();j++) {
				cout << l_combs[i][j] << " ";
			}
			cout << "\n";
		}
	}
}

void CallReadTypesGaps (run_params p, int max_l, int max_r_size, vector<int> polys, vector< vector<int> >& l_combs) {
	int s=polys.size();
	if (p.verb==1) {
		cout << "Potential sets of loci spanned by reads: here we consider the full set plus a certain number of indels\n";
	}
	vector <int> l_set;
	for (int i=0;i<s;i++) {
		l_set.push_back(i);
	}
	l_combs.push_back(l_set);

	for (int i1=0;i1<s;i1++){
		l_set.clear();
		for (int i2=0;i2<s;i2++) {
			if (i2!=i1) {
				l_set.push_back(i2);
			}
		}
		l_combs.push_back(l_set);
	}
	if (p.maxgap>=2) {
		for (int i1=0;i1<s;i1++) {
			for (int i2=0;i2<s;i2++) {
				if (i1!=i2) {
					l_set.clear();
					for (int i=0;i<s;i++) {
						if (i!=i1&&i!=i2) {
							l_set.push_back(i);
						}
					}
					l_combs.push_back(l_set);
				}
			}
		}

	}
	
	if (p.verb==1) {
		cout << "l_combs\n";
		for	(int i=0;i<l_combs.size();i++) {
			for (int j=0;j<l_combs[i].size();j++) {
				cout << l_combs[i][j] << " ";
			}
			cout << "\n";
		}
	}
}



void CallMNPs2 (run_params p, vector< vector<int> >& l_combs, vector <vector<char> > m_vars, vector< vector<mpoly> >& m_pol) {
	cout << "Calling multi-locus polymorphisms...\n";
	vector<int> done (m_vars.size(),0);
	
	for (int i=0;i<l_combs.size();i++) { // Scan across potential sets of polymorphic loci
		if (p.verb==1) {
			cout << "Calling: set " << i << " of " << l_combs.size() << "\n";
		}
		vector<mpoly> mp; // Vector of type, number
		for (int k=0;k<m_vars.size();k++) {  //Sequence in set at given time
			int inc=1; //Flag for read fitting category
			vector<char> add = m_vars[k]; //Specific nucleotides being scanned for
			if (done[k]==1) {
				inc=0;
			} else {
				for (int l=0;l<l_combs[i].size();l++) {
					add[l_combs[i][l]]='-';
					if (m_vars[k][l_combs[i][l]]=='-') {
						inc=0;
					}
				}
			}
			if (inc==1) { //Included in this set
				//Check against existing polymorphisms in this category
				int a=0;
				for (int l=0;l<mp.size();l++) {
					if (m_vars[k]==mp[l].vars) {
						mp[l].count++;
						done[k]=1;
						a=1;
					}
				}
				if (a==0) {
					//	cout << "No category yet - check\n";
					mpoly m;
					//Check here - if there are additional nucleotides to those being looked for, can retain those
					int extra=0;
					for (int l=0;l<add.size();l++) {
						if (add[l]!='-') {
							extra=1;
						}
					}
					if (extra==0) { //No extra nucleotides
						m.vars=m_vars[k];
						m.count=1;
						mp.push_back(m);
						done[k]=1;
						
					} else {	//Extra nucleotides exist - split read.  Add what is being searched for to the list
							
						//Create new variant with only the nucleotide information being looked for
						vector<char> new_var = m_vars[k];
						for (int l=0;l<add.size();l++) {
							if (add[l]!='-') {
								new_var[l]='-';
							}
						}
							
						//Now check against existing polymorphisms
						int b=0;
						for (int l=0;l<mp.size();l++) {
							if (new_var==mp[l].vars) {
								mp[l].count++;
								done[k]=1;
								b=1;
							}
						}
						if (b==0) {
							m.vars=new_var;
							m.count=1;
							mp.push_back(m);
							m_vars[k]=add;
						}
					}
				}
			}
		}
		m_pol.push_back(mp);
	}
	
	//Output of dataset
	if (p.verb==1) {
		for (int j=0;j<m_pol.size();j++) {
			cout << j << " " << m_pol[j].size() << "\n";
			for (int k=0;k<m_pol[j].size();k++) {
				for (int l=0;l<m_pol[j][k].vars.size();l++) {
					cout << m_pol[j][k].vars[l];
				}
				cout << " " << m_pol[j][k].count << "\n";
			}
		}
	}
}


void FindMaxLoc(int& max_loc, vector<vector< vector<int> > > all_l_combs) {
	for (int i=0;i<all_l_combs.size();i++) {
		for (int j=0;j<all_l_combs[i].size();j++) {
			for (int k=0;k<all_l_combs[i][j].size();k++) {
				if (all_l_combs[i][j][k]>max_loc) {
					max_loc=all_l_combs[i][j][k];
					cout << "Max_loc " << max_loc << "\n";
				}
			}
		}
	}
}

void ConstructLocVec (run_params p, int max_loc, int max_l_store, vector<vector< vector<int> > > all_l_combs, vector< vector<int> >& l_vec) {
	/*cout << "Check all_l_combs\n";
	for (int i=0;i<all_l_combs.size();i++) {
		for (int j=0;j<all_l_combs[i].size();j++) {
			for (int k=0;k<all_l_combs[i][j].size();k++) {
				cout << all_l_combs[i][j][k] << " ";
			}
			cout << "\n";
		}
	}*/
	
	for (int ml=max_l_store;ml>0;ml--) {
		//cout << "ML " << ml << " " << "\n";
		for (int index=0;index<max_loc+1;index++) {
			//cout << "Index " << index << "\n";
			int found=0;
			for	(int i=0;i<all_l_combs.size();i++) {
				if (found==0) {
					for (int j=0;j<all_l_combs[i].size();j++){
						if (all_l_combs[i][j].size()==ml&&all_l_combs[i][j][0]==index) {
							//cout << "Match " << i << " " << j << " to size " << ml << " first " << index << "\n";
							if (l_vec.size()>0) {
								//cout << "Size OK\n";
								l_vec.push_back(all_l_combs[i][j]);
								found=1;
								//cout << "Found 1\n";
								break;
							} else {
								l_vec.push_back(all_l_combs[i][j]);
								found=1;
								//cout << "Found 2\n";
							break;
							}
						}
					}
				}
			}
		}
	}
	if (p.verb==1) {
		cout << "Loc_vec\n";
		for (int i=0;i<l_vec.size();i++) {
			for (int j=0;j<l_vec[i].size();j++) {
				cout << l_vec[i][j] << " ";
			}
			cout << "\n";
		}
	}
}


void ConstructLocVecGap (run_params p, vector<vector< vector<int> > > all_l_combs, vector< vector<int> >& l_vec) {
	for	(int i=0;i<all_l_combs.size();i++) {
		for (int j=0;j<all_l_combs[i].size();j++){
			int inc=0;
			for (int k=0;k<l_vec.size();k++) {
				if (all_l_combs[i][j]==l_vec[k]) {
					inc=1;
				}
			}
			if (inc==0) {
				l_vec.push_back(all_l_combs[i][j]);
			}
		}
	}
	if (p.verb==1) {
		cout << "Loc_vec\n";
		for (int i=0;i<l_vec.size();i++) {
			for (int j=0;j<l_vec[i].size();j++) {
				cout << l_vec[i][j] << " ";
			}
			cout << "\n";
		}
	}
}

		

void ConstructMPoly (vector< vector<int> >& l_vec, vector<vector< vector<int> > > all_l_combs, vector< vector< vector<mpoly> > >& m_polys) {
	vector< vector< vector<mpoly> > > new_m_polys;
	for (int i=0;i<m_polys.size();i++) {
		vector< vector<mpoly> > mps;
		for (int j=0;j<l_vec.size();j++) {
			vector<mpoly> mp;
			//Find set that corresponds to these loci
			for (int k=0;k<all_l_combs[i].size();k++) {
				if (all_l_combs[i][k]==l_vec[j]) {
					mp=m_polys[i][k];
				}
			}
			mps.push_back(mp);
		}
		new_m_polys.push_back(mps);
	}
	
	m_polys=new_m_polys;
}


void MakeAssignedList (vector< vector <vector<char> > > t_m_vars, vector< vector<int> >& done) {
	for (int i=0;i<t_m_vars.size();i++) {
		vector<int> d (t_m_vars[i].size(),0);
		done.push_back(d);
	}
}

//Need to check around here - seem to be losing some information...

void FilterMNPs (run_params p, vector< vector< vector<mpoly> > > m_polys, vector< vector< vector<mpoly> > >& m_polys_f) {
	if (p.verb==1) {
		cout << "Filtering MNPs\n";
	}
	for (int i=0;i<m_polys.size();i++) {
		vector< vector<mpoly> > mp2;
		for (int j=0;j<m_polys[i].size();j++) { //Class of read types
			vector<mpoly> mp;
			int tot=0;
			for (int k=0;k<m_polys[i][j].size();k++) { //Loop within read type
				cout << i << " " << j << " " << k << " " << m_polys[i][j][k].count << "\n";
				tot=tot+m_polys[i][j][k].count;
			}
			for (int k=0;k<m_polys[i][j].size();k++) { //Loop within read type
				if (m_polys[i][j][k].count>=p.hap_n_min&&(m_polys[i][j][k].count+0.)/(tot+0.)>p.hap_q_cut) {
					mpoly m;
					m.vars=m_polys[i][j][k].vars;
					m.count=m_polys[i][j][k].count;
					mp.push_back(m);
				}
			}
			mp2.push_back(mp);
		}
		m_polys_f.push_back(mp2);
	}
	if (p.verb==1) {
		for (int i=0;i<m_polys_f.size();i++) {
			cout << "Time point " << i << "\n";
			for (int j=0;j<m_polys_f[i].size();j++) {
				//cout << j << " " << m_polys_f[i][j].size() << "\n";
				for (int k=0;k<m_polys_f[i][j].size();k++) {
					for (int l=0;l<m_polys_f[i][j][k].vars.size();l++) {
						cout << m_polys_f[i][j][k].vars[l];
					}
					cout << " " << m_polys_f[i][j][k].count << "\n";
				}
			}
		}
	}
}

void CombineMNPsTime (run_params p, vector< vector< vector<mpoly> > > m_polys, vector< vector<mpoly> >& c_m_polys) {
	c_m_polys = m_polys[0];
	for (int i=1;i<m_polys.size();i++) {
		for (int j=0;j<m_polys[i].size();j++) { //Class of read types
			for (int k=0;k<m_polys[i][j].size();k++) { //Loop within read type
				//Check against existing set in c_m_polys
				int done=0;
				for (int l=0;l<c_m_polys[j].size();l++) {
					if (c_m_polys[j][l].vars==m_polys[i][j][k].vars) {
						c_m_polys[j][l].count=c_m_polys[j][l].count+m_polys[i][j][k].count;
						done=1;
					}
				}
				if (done==0) {
					c_m_polys[j].push_back(m_polys[i][j][k]);
				}
			}
		}
	}
	
	//Output of dataset
	if (p.verb==1) {
		cout << "Combined dataset\n";
		for (int i=0;i<c_m_polys.size();i++) {
			for (int j=0;j<c_m_polys[i].size();j++) {
				for (int k=0;k<c_m_polys[i][j].vars.size();k++) {
					cout << c_m_polys[i][j].vars[k];
				}
				cout << " " << c_m_polys[i][j].count << "\n";
			}
		}
	}
}

void MLPMeanFreqs (run_params p, vector< vector<mtr> >& mltrajs) {
	for (int i=0;i<mltrajs.size();i++) {
		vector<int> totals;
		int t_all=0;
		for (int j=0;j<mltrajs[i].size();j++) {
			int tot=0;
			for (int k=0;k<mltrajs[i][j].n.size();k++) {
				tot=tot+mltrajs[i][j].n[k];
			}
			totals.push_back(tot);
			t_all=t_all+tot;
		}
		for (int j=0;j<mltrajs[i].size();j++) {
			mltrajs[i][j].m=(totals[j]+0.)/(t_all+0.);
		}
	}
	
	if (p.verb==1) {
		cout << "Multi-locus trajectory reads and mean:\n";
		for (int i=0;i<mltrajs.size();i++) {
			cout << "Set " << i << "\n";
			for (int j=0;j<mltrajs[i].size();j++) {
				for (int k=0;k<mltrajs[i][j].seq.size();k++) {
					cout << mltrajs[i][j].seq[k];
				}
				cout << " ";
				for (int k=0;k<mltrajs[i][j].n.size();k++) {
					cout << mltrajs[i][j].n[k] << " ";
				}
				cout << " ";
				cout << mltrajs[i][j].m << "\n";
			}
		}
	}

}

void GetIncludeML (run_params p, vector<int>& include_ml, vector< vector<mtr> > mltrajs) {
	if (p.verb==1) {
		cout << "Included ML trajectories:\n";
	}
	//Filters by mean frequency - consider including points above freuqency threshold
	for (int i=0;i<mltrajs.size();i++) {
		int inc=0;
		for (int j=0;j<mltrajs[i].size();j++) {
			if (mltrajs[i][j].m>p.hap_q_cut&&mltrajs[i][j].m<(1-p.hap_q_cut)) {
				inc=1;
			}
		}
		include_ml.push_back(inc);
		if (p.verb==1) {
		cout << i << " " << include_ml[i] << "\n";
		}
	}
}

void FilterMLFreq (run_params p, vector< vector<int> >& times, vector< vector<mtr> >& mltrajs) {
	if (p.verb==1) {
		cout << "Filter ML trajectories 1\n";
	}
	
	//Filter out times at which there are not polymorphic haplotypes.
	vector< vector<mtr> > new_mltrajs;
	vector< vector<int> > new_times;
	//Get sequence depths
	for (int i=0;i<mltrajs.size();i++) {
		if (mltrajs[i].size()>0) {
			vector<int> depths;
			GetDepths (i,depths,mltrajs);
			//Filter out times at which there are no polymorphic haplotypes
			vector<int> included;
			for (int k=0;k<mltrajs[i][0].n.size();k++) {
				int inc=0;
				for (int j=0;j<mltrajs[i].size();j++) { //Cycle through partial haplotypes
					if ((mltrajs[i][j].n[k]+0.)/(depths[k]+0.)>p.hap_q_cut&&(mltrajs[i][j].n[k]+0.)/(depths[k]+0.)<(1-p.hap_q_cut)) {
						inc=1;
					}
				}
				included.push_back(inc);
			}

			vector<mtr> mt;
			vector<int> t;
			for (int k=0;k<times[i].size();k++) {
				if (included[k]==1) {
					t.push_back(times[i][k]);
				}
			}
		
			for (int j=0;j<mltrajs[i].size();j++) {
				mtr m;
				for (int k=0;k<mltrajs[i][j].n.size();k++) {
					if (included[k]==1) {
						m.n.push_back(mltrajs[i][j].n[k]);
					}
				}
				mt.push_back(m);
			}
			if (mt[0].n.size()>1) {
				new_mltrajs.push_back(mt);
				new_times.push_back(t);
			}
		}
	}
	mltrajs=new_mltrajs;
	times=new_times;
}


void FilterMLFreq2 (run_params p, vector< vector<int> >& times, vector< vector<mtr> >& mltrajs) {
	if (p.verb==1) {
		cout << "Filter ML trajectories:\n";
	}
	
	vector< vector<mtr> > new_mltrajs;
	vector<vector <int> > new_times;
	//Get sequence depths
	for (int i=0;i<mltrajs.size();i++) {
	//	cout << "Check " << i << "\n";
		vector<int> depths;
		GetDepths (i,depths,mltrajs);
		//Filter out haplotypes at which there are no polymorphisms
		vector<int> included;
		for (int j=0;j<mltrajs[i].size();j++) { //Cycle through partial haplotypes
			int inc=1;
			for (int k=0;k<mltrajs[i][0].n.size();k++) {
				if ((mltrajs[i][j].n[k]+0.)/(depths[k]+0.)<p.hap_q_cut||(mltrajs[i][j].n[k]+0.)/(depths[k]+0.)>(1-p.hap_q_cut)) {
					inc=0;
				}
			}
			included.push_back(inc);
		}
		
		vector<mtr> mt;
		for (int j=0;j<mltrajs[i].size();j++) {
			if (included[j]==1) {
				mtr m;
				for (int k=0;k<mltrajs[i][j].n.size();k++) {
					m.n.push_back(mltrajs[i][j].n[k]);
				}
				mt.push_back(m);
			}
		}
		int done_t=0;
		//cout << "Size " << i << " " << mt.size() << "\n";
		if (mt.size()>1) {
			new_mltrajs.push_back(mt);
			if (done_t==0) {
				new_times.push_back(times[i]);
				done_t=1;
			}
		}
	}
	mltrajs=new_mltrajs;
	times=new_times;
	
	
}


void FilterMLFreq3 (run_params p, vector<int>& include_ml, vector< vector<int> >& times, vector< vector<mtr> >& mltrajs) {
	for (int i=0;i<mltrajs.size();i++) {
		if (include_ml[i]==1) { //Calculate sample depths
			vector<int> depths;
			GetDepths (i,depths,mltrajs);
		
			double max_dt=0;
			for (int j=0;j<mltrajs[i].size();j++) { //Cycle through partial haplotypes
				//Assess change in frequency with time
				double dq=0;
				double dqq=0;
				double dt=0;
				for (int k=0;k<mltrajs[i][j].n.size()-1;k++) {
					dqq=((mltrajs[i][j].n[k+1]+0.)/(depths[k+1]+0.))-((mltrajs[i][j].n[k]+0.)/(depths[k]+0.));
					if (dqq<0) {dqq=-dqq;}
					dq=dq+dqq;
					dt=dt+(times[i][k+1]-times[i][k]);
				}
				dqq=dq/dt;
				if (dqq>max_dt) {max_dt=dqq;}
			}
			if (max_dt>p.dq_cut) {
				include_ml[i]=0;
			}
		}
	}
}

void GetDepths (int i, vector<int>& depths, vector< vector<mtr> >& mltrajs) {
	for (int k=0;k<mltrajs[i][0].n.size();k++) {
		int d=0;
		for (int j=0;j<mltrajs[i].size();j++) { //Cycle through partial haplotypes
			d=d+mltrajs[i][j].n[k];
		}
		depths.push_back(d);
	}
}

void CompleteMLCalls (vector< vector<char> > haps, vector< vector<mpoly> >& c_m_polys) {//Includes data from imported called haplotypes
	cout << "Edit c_m_polys\n";
	for (int i=0;i<c_m_polys.size();i++) {
		//Positions
		vector<int> pos;
		if (c_m_polys[i].size()>0) {
			for (int k=0;k<c_m_polys[i][0].vars.size();k++) {
				if (c_m_polys[i][0].vars[k]!='-') {
					pos.push_back(k);
				}
			}
			
			for (int j=0;j<haps.size();j++) { //Loop over haplotypes.  If there is no equivalent variant, create it
				vector<char> var (haps[j].size(),'-');
				for (int k=0;k<pos.size();k++) {
					var[pos[k]]=haps[j][pos[k]];
				}
				
				//Compare against all partial haplotypes
				int seen=0;
				for (int k=0;k<c_m_polys[i].size();k++) {
					if (var==c_m_polys[i][k].vars) {
						seen=1;
					}
				}
				if (seen==0) { //Add haplotype
					mpoly m;
					m.vars=var;
					m.count=0;
					c_m_polys[i].push_back(m);
				}
			}
		}
	}
}
