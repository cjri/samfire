#include "shared.h"
#include "matchpairs.h"
#include <iostream>
#include <string>
#include <sstream>

void MatchPairs (run_params p, int ireq, vector<rd>& data) { //ireq flags whether to condition on prior inclusion
	cout << "Identifying paired end reads...\n";
	int n_pairs=0;
	for (unsigned int i=0;i<data.size();i++) {
//		cout << data[i].paircode << " " << data[i].inc << "\n";
//		cout << "i " << i << "\n";
		if (ireq==0||data[i].inc!=0) {
			int j_lim=data.size();
			if (p.sorted==1) {
				j_lim=i+10;
				if (j_lim>data.size()) {
					j_lim=data.size();
				}
			}
			for (unsigned int j=i+1;j<j_lim;j++) {
				if (ireq==0||data[j].inc!=0) {
					if (data[i].paircode==data[j].paircode) {
						//cout << i << " " << j << " " << data[i].paircode << " " << data[j].paircode << "\n";
						n_pairs++;
						data[i].pairno=j;
						data[j].pairno=i;
						break;
					}
				}
			}
		}
	}
	cout << "Number of pairs identified " << n_pairs << "\n";
}

void CheckOldPairs (vector<rd>& data) {
	cout << "Filter identified pairs...\n";
	for (unsigned int i=0;i<data.size();i++) {
		if (data[i].inc==0||data[i].del==1) {
			data[i].pairno=-1;
		}
	}
}

void JoinPairs (int i, run_params p, vector<joined>& jseqs, vector<char>& qual, vector<rd> data) {
	cout << "Joining paired-end reads...\n";
	ofstream jp_file;
    ofstream qp_file;
	ostringstream convert;
	convert << i;
	string temp=convert.str();
	string name = "Joined"+temp+".out";
	string qname = "Quality"+temp+".out";
	cout << name << "\n";
	jp_file.open(name.c_str());
	if (p.print_qual>0) {
		qp_file.open(qname.c_str());
	}
	if (p.pairs==1) {
		JoinWithPairs(p,jseqs,qual,data,jp_file,qp_file);
	} else {
		JoinNoPairs(p,jseqs,qual,data,jp_file,qp_file);
	}
}

void JoinWithPairs (run_params p, vector<joined>& jseqs, vector<char>& qual, vector<rd> data, ofstream& jp_file, ofstream& qp_file) {
	int foundpair=0;
	int n_seqs=0;
	for (unsigned int i=0;i<data.size();i++) {
		if (data[i].inc!=0&&data[i].del==0) {
			n_seqs++;
			joined jj;
			if (data[i].pairno!=-1) {
				//cout << "Existing pair at " << i << "\n";
				foundpair++;
				//Have pair  !!! Note at this point that j might not be included...
				int j=data[i].pairno;
				int i2=i;
				int j2=j;
				
				if (j2>i2) {
					//Paired with higher sequence - avoid duplicates
					
					if (data[i].alpos>data[j].alpos) {  //Order sequences so that i is the first one
						int temp=i2;
						i2=j2;
						j2=temp;
					}
					int pos1=0;
					int pos2=0;
					int loc=data[i2].alpos;
					string iseq;
					string jseq;
                    string iqual;
                    string jqual;
					//Assign sequences
					
					if (data[i2].inc==1) {  //At least one direction is included for i
						iseq=data[i2].seq;
                        iqual=data[i2].qual;
					} else {
						iseq=data[i2].revseq;
                        iqual=data[i2].revqual;
					}
					
					//Account for possibility that j is not included
					if (data[j2].inc==0||data[j2].del==1) {
						jseq="";
                        jqual="";
					} else if (data[j2].inc==1) {
						jseq=data[j2].seq;
                        jqual=data[j2].qual;
					} else {
						jseq=data[j2].revseq;
                        jqual=data[j2].revqual;
					}
					
			//		cout << "i " << i2 << " " << iseq << "\n";
			//		cout << "j " << j2 << " " << jseq << "\n";
			
					
					string newstr;
                    string newqual;
					if (jseq.length()>0) {
						while (loc<data[j2].alpos+jseq.length()) {
							if (loc<data[i2].alpos+iseq.length()) {
								//	cout << "Here1 " << loc << " " << pos1 << " " << seqs[i].seq[pos1] << "\n";
								newstr.push_back(iseq[pos1]);
                                newqual.push_back(iqual[pos1]);
								pos1++;
								loc++;
								if (loc>data[j2].alpos) {
									pos2++;
								}
							} else if (loc<data[j2].alpos){
								//	cout << "Here2\n";
								newstr.push_back('-');
                                newqual.push_back('-');
								loc++;
							} else {
								//cout << "Here3 " << loc << " " << pos2 << " " << seqs[j2].seq[pos2] << " " << seqs[j2].seq[pos2-1] << "\n";
								newstr.push_back(jseq[pos2]);
                                newqual.push_back(jqual[pos2]);
								pos2++;
								loc++;
							}
						}
					} else {
						newstr=iseq;
                        newqual=iqual;
					}
					jj.alpos=data[i2].alpos;
					jj.seq=newstr;
			//		cout << "ij " << data[i].alpos << " " << newstr << "\n";
					jp_file << data[i2].alpos << " " << newstr << "\n";
                    if (p.print_qual==1) {
                        qp_file << newqual << "\n";
                    }
                    if (p.print_qual>1) {
                        TranslateQualityVector(newqual,qual,qp_file);
                    }

					jseqs.push_back(jj);
				}
			} else {
				//Don't have pair - just report single value
				joined jj;
				string iseq;
				jj.alpos=data[i].alpos;
				if (data[i].inc==1&&data[i].del==0) {
					iseq=data[i].seq;
                    if (iseq.length()>0&&p.print_qual==1) {
                        qp_file << data[i].qual << "\n";
                    }
                    if (iseq.length()>0&&p.print_qual>1) {
                        TranslateQualityVector(data[i].qual,qual,qp_file);
                    }
				} else {
					iseq=data[i].revseq;
                    if (iseq.length()>0&&p.print_qual==1) {
                        qp_file << data[i].revqual << "\n";
                    }
                    if (iseq.length()>0&&p.print_qual>1) {
                        TranslateQualityVector(data[i].revqual,qual,qp_file);
                    }

				}
				jj.seq=iseq;
				if (iseq.length()>0) {
					jp_file << data[i].alpos << " " << iseq << "\n";
					jseqs.push_back(jj);
				}
			}
		}
	}
//	cout << "Data size " << data.size() << "\n";
//	cout << "Number of pairs found " << n_seqs << "\n";
//	cout << "Number of pairs found " << foundpair << "\n";
}

void JoinNoPairs (run_params p, vector<joined>& jseqs, vector<char>& qual, vector<rd> data, ofstream& jp_file, ofstream& qp_file) {
    for (unsigned int i=0;i<data.size();i++) {
		if (data[i].inc==1&&data[i].del==0) {
			joined jj;
			jj.alpos=data[i].alpos;
			jj.seq=data[i].seq;
			jp_file << data[i].alpos << " " << data[i].seq << "\n";
			//cout << data[i].alpos << " " << data[i].seq << "\n";
            if (p.print_qual==1) {
                qp_file << data[i].qual << "\n";
            }
            if (p.print_qual>1) {
                TranslateQualityVector(data[i].qual,qual,qp_file);
            }
			jseqs.push_back(jj);
		} else if (data[i].inc==-1&&data[i].del==0) {
			joined jj;
			jj.alpos=data[i].alpos;
			jj.seq=data[i].revseq;
			jp_file << data[i].alpos << " " << data[i].revseq << "\n";
            if (p.print_qual==1) {
                qp_file << data[i].revqual << "\n";
            }
            if (p.print_qual>1) {
                TranslateQualityVector(data[i].revqual,qual,qp_file);
            }
			jseqs.push_back(jj);
		}
	}
}

void TranslateQualityVector(string s, vector<char>& qual, ofstream& qp_file) {
    for (int k=0;k<s.size();k++) {
        for (int j=0;j<=qual.size();j++) {
            if (s[k]==qual[j]) {
                qp_file << j << " ";
            }
        }
    }
    qp_file << "\n";
}
