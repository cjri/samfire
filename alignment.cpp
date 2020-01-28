#include "shared.h"
#include "alignment.h"
#include "utilities_sam.h"
#include <iostream>
#include <sstream>
#include <string>

void AlignSequencesSam (run_params& p, int s_length, vector<char> qual, rseq refseq, vector<rd>& data) {
	int q=0;
	cout << "Length " << s_length << "\n";
	string q0;
	int xalq=0;
	MinBaseQual(qual,q0); //Find code for minimum base quality
	for (int i=0;i<s_length;i++) {
		if (p.verb==1) {
			cout << "Sequences " << i << "\n";
			cout << data[i].seq << "\n";
			cout << data[i].qual << "\n";
		}
		data[i].inc=0;
		if (data[i].del==0) { //Not previously deleted
			//cout << data[i].ref << " " << refseq.name << "\n";
			if (data[i].ref==refseq.name) {
				ReadCigar(i,data); //Read CIGAR string
				q=findqual(p,i,p.min_qual,p.max_qual,qual,data);  //Assess read by median nucleotide quality
				//cout << "Median qual " << q << "\n";
				if (q>=p.min_qual) {  //Minimum sequence quality test
					//cout << "Align qual " << data[i].alq << " " << p.ali_qual << "\n";
					if (data[i].alq>=p.ali_qual) { //Minumum alignment quality
						if (p.ali_inc==1||(p.ali_inc==0&&data[i].alq!=255)) { //Minumum alignment quality
							//cout << data[i].alpos << " " << refseq.seq.size()-p.min_rlen << "\n";
							if (data[i].alpos<refseq.seq.size()-p.min_rlen) { //Sequence must be aligned to have an overlap of at least the minimum number of alleles reported by a read; this weeds out misaligned sequences at the end of the reference sequence.
								//cout << "Included\n";
								data[i].inc=1;
								int isc=0;
								int fsc=0;
								RemoveInitialSoftClipping(i,isc,data);  //Remove initial soft clipping
								FixDeletions(i,q0,data);  //Include blank data for deletions
								FixInsertions(i,data);  //Remove insertions
								if (p.trim>0) {
									CountEndSoftClipping(i,fsc,data); //Count soft clipping at end of read - required for trim;
								}
								RemoveSoftClipping(i,data);  //Remove remaining soft clipping
								ProcessReadQual (i,isc,fsc,p,qual,data); //Trim data + Process by individual nucleotide quality.
							}
						} else if (data[i].alq==255) {
							xalq++;
						}
					}
				}
			}
		}
	}
	if (xalq>200) {
		cout << "Warning: " << xalq << " sequences excluded due to unknown alignment quality\n";
	}
}


void MinBaseQual (vector<char> qual, string& q0) {
	stringstream ss;
	char c = qual[0];
	ss << c;
	ss >> q0;
}

//Read and process CIGAR string
void ReadCigar (int i, vector<rd>& data) {
	int pos=0;
	for (int j=0;j<data[i].cigar.size();j++) {
		if (data[i].cigar.compare(j,1,"M")==0) {
			int i_num=atoi(data[i].cigar.substr(pos,j).c_str());
			//			cout << "M " << i_num << "\n";
			pos=j+1;
			string app (i_num,'M');
			data[i].cig_string.append(app);
		}
		if (data[i].cigar.compare(j,1,"S")==0) {
			int i_num=atoi(data[i].cigar.substr(pos,j).c_str());
			//			cout << "S " << i_num << "\n";
			pos=j+1;
			string app (i_num,'S');
			data[i].cig_string.append(app);
		}
		
		if (data[i].cigar.compare(j,1,"H")==0) {
			pos=j+1;
		}

		if (data[i].cigar.compare(j,1,"I")==0) {
			int i_num=atoi(data[i].cigar.substr(pos,j).c_str());
			//			cout << "I " << i_num << "\n";
			pos=j+1;
			string app (i_num,'I');
			data[i].cig_string.append(app);
		}
		if (data[i].cigar.compare(j,1,"D")==0) {
			int i_num=atoi(data[i].cigar.substr(pos,j).c_str());
			//			cout << "D " << i_num << "\n";
			pos=j+1;
			string app (i_num,'D');
			data[i].cig_string.append(app);
		}
		if (data[i].cigar.compare(j,1,"P")==0) {  //Padding
			int i_num=atoi(data[i].cigar.substr(pos,j).c_str());
			//			cout << "P " << i_num << "\n";
			pos=j+1;
			string app (i_num,'P');
			data[i].cig_string.append(app);
		}
	}
	FilterCigar(i,data);
}

//Remove padding from cigar string
void FilterCigar (int i, vector<rd>& data) {
	int j=0;
	while (j<data[i].cig_string.length()) {
		if (data[i].cig_string.compare(j,1,"P")==0) {
			data[i].cig_string.erase(j,1);
		} else {
			j++;
		}
	}
}

//Check median quality of sequence.  Reduce sequence to achieve median quality if required.  Reduction is carried out from both ends; get the longest qualifying sequence
int findqual (run_params p, int i, int min_qual, int max_qual, vector<char> qual, vector<rd>& data) {
	string q=data[i].qual;
	string s=data[i].seq;
	string c=data[i].cig_string;
	data[i].rmqual.push_back(0);
	data[i].rmqual.push_back(0);
	//	cout << s << "\n";
	//	cout << q << "\n";
	vector<int> qvec;
	for (int i=0;i<s.size();i++) {
		int done=0;
		for (int j=0;j<=max_qual;j++) {
			if (q[i]==qual[j]) {
				//cout << i << " " << q.size() << " " << q[i] << " " << qual[j] << "\n";
				if (p.qlib==3) {
					if (j==0) {j=1;}
					int qq=floor(10.*log10(j));
					//cout << q[i] << " "  << i << " " << q.size() << " " << qual[j] << " " << j << " " << qq << "\n";
					qvec.push_back(qq);
					done=1;
				} else {
					qvec.push_back(j);
					done=1;
				}
				break;
			}
		}
		if (done==0) {
			cout << "Error: No match to base quality score: " << q[i] << "\n";
			qvec.push_back(0);
		}
	}
	
	int median = GetMedian(0,0,qvec);
	//cout << "Median " << median << " ";
	int qo=median;
	
	if (qvec.size()>=p.min_rlen) {
		if (median<min_qual) {  //Edit sequence to get high quality part
			int s1=999;
			int s2=999;
			//Reduce sequence by removing nucleotides from the end of the read
			for (int a=1;a<qvec.size()-p.min_rlen;a++) {
				median=GetMedian(a,0,qvec);
				if (median>=min_qual) {
					s1=a;
					break;
				}
			}
			for (int a=1;a<qvec.size()-p.min_rlen;a++) {
				median=GetMedian(0,a,qvec);
				if (median>=min_qual) {
					s2=a;
					break;
				}
			}
			int qs=q.size();
			
			if (s1<999&&s1<=s2) {
				q=q.substr(s1,qs-s1+1);
				s=s.substr(s1,qs-s1+1);
				if (c.length()>0) {
					c=c.substr(s1,qs-s1+1);
				}
				data[i].qual=q;
				data[i].seq=s;
				data[i].cig_string=c;
				qo=min_qual;
				data[i].rmqual[0]=s1;
			} else if (s2<999&&s2<s1) {
				q=q.substr(0,qs-s2+1);
				s=s.substr(0,qs-s2+1);
				if (c.length()>0) {
					c=c.substr(0,qs-s2+1);
				}
				data[i].qual=q;
				data[i].seq=s;
				data[i].cig_string=c;
				qo=min_qual;
				data[i].rmqual[1]=s2;
			} else {
				qo=-1;
			}
		}
	} else {
		qo=-1;
	}
	//	cout << "Quality " << qo << "\n";
	
	if (qo>-1) { //Check number of minimum quality reads
		int nm=0;
		for (int i=0;i<s.size();i++) {
			for (int j=min_qual;j<=max_qual;j++) {
				if (q[i]==qual[j]) {
					nm++;
					break;
				}
			}
		}
		if (nm<p.min_rlen) {
			qo=-1;
		}
	}
	return qo;
}

int GetMedian (int a, int b, vector<int> qvec) {
	vector<int> qvec2=qvec;
	sort(qvec2.begin()+a,qvec2.end()-b);
	double median;
	double size=qvec2.size();
	if (qvec2.size()==0) {
		median=(qvec2[size/2-1]+qvec2[size/2])/2;
	} else {
		median = qvec2[size/2];
	}
	return median;
}


void RemoveInitialSoftClipping (int i, int& isc, vector<rd>& data) {
	if (data[i].cig_string.compare(0,1,"S")==0) {
		string s;
		string q;
		string c;
		int j=0;
		while (data[i].cig_string.compare(j,1,"S")==0) {
			j++;
		}
		isc=j;
		data[i].seq=data[i].seq.substr(j,data[i].seq.length()-j);
		data[i].qual=data[i].qual.substr(j,data[i].qual.length()-j);
		data[i].cig_string=data[i].cig_string.substr(j,data[i].cig_string.length()-j);
	}
}

void FixDeletions (int i, string q0, vector<rd>& data) {
	for (int j=0;j<data[i].cig_string.length();j++) {
		if (data[i].cig_string.compare(j,1,"D")==0) {
			data[i].seq.insert(j,"-");
			data[i].qual.insert(j,q0);
		}
	}
}

void FixInsertions (int i, vector<rd>& data) {
	for (int j=0;j<data[i].cig_string.length();j++) {
		if (data[i].cig_string.compare(j,1,"I")==0) {
			while (data[i].cig_string.compare(j,1,"I")==0&&j<data[i].cig_string.length()) {
				data[i].seq.erase(j,1);
				data[i].qual.erase(j,1);
				data[i].cig_string.erase(j,1);
			}
		}
	}
}

void CountEndSoftClipping (int i, int& fsc, vector<rd>& data) {
	if (data[i].cig_string.size()>0) {
		int j=1;
		while (data[i].cig_string.compare(data[i].cig_string.size()-j,1,"S")==0) {
			j++;
		};
		fsc=j-1;
	}
}

void RemoveSoftClipping (int i, vector<rd>& data) {
	for (int j=0;j<data[i].cig_string.length();j++) {
		if (data[i].cig_string.compare(j,1,"S")==0) {
			while (j<data[i].cig_string.length()&&data[i].cig_string.compare(j,1,"S")==0) {
				data[i].seq.erase(j,1);
				data[i].qual.erase(j,1);
				data[i].cig_string.erase(j,1);
			}
		}
	}
}

//Process short read, removing alleles that do not have sufficient nucleotide quality
void ProcessReadQual (int i, int isc, int fsc, run_params p, vector<char> qual, vector<rd>& data) {

	//Remove ends by trimming
	if (p.trim>0) {
		if (data[i].seq.size()<(2*p.trim)) {
			data[i].inc=0;
		} else {
			int trim=p.trim-data[i].rmqual[0];
			trim=trim-isc;
			if (trim>0) {
				data[i].seq.erase(0,trim);
				data[i].alpos=data[i].alpos+trim;
			}
			//Find final soft clipping
			trim=p.trim-data[i].rmqual[1];
			trim=trim-fsc;
			if (trim>0) {
				data[i].seq.erase(data[i].seq.end()-p.trim,data[i].seq.end());
			}
		}
	}
	
	//Count nucleotides of sufficient quality
	int incl=0;
	for (int k=0;k<data[i].seq.size();k++) {
		int keep=0;
		for (int j=p.min_qual;j<=p.max_qual;j++) {
			if (data[i].qual[k]==qual[j]) {
				keep=1;
				break;
			}
		}
		if (keep==0) {
			data[i].seq.replace(k,1,"-");
		} else {
			incl++;
		}
	}
	if (incl<p.min_rlen) {
		data[i].inc=0;
	}
}





void AlignSequencesSlide (int s_length, run_params& p, vector<char> qual, rseq refseq, vector<rd>& data) {
	int q=0;
	cout << "Length " << s_length << "\n";
	for (int i=0;i</*s_length*/s_length;i++) {
		if (p.verb==1) {
			cout << "Sequence " << i << "\n";
		}
		data[i].inc=0;
		if (data[i].refno!=-1) {
			int isok=0;
			int firstpos=-1;
			q=findqual(p,i,p.min_qual,p.max_qual,qual,data);  //Check what is going on here - does this need further revision?
			if (q>=p.min_qual) {  //Minimum sequence quality
				DoForwardSlideAlignment (p,i,firstpos,isok,qual,refseq,data);
				if (isok==0) {
					isok=1;
					int firstpos=-1;
					DoReverseSlideAlignment (p,i,firstpos,isok,qual,refseq,data);
				}
			}
		}
	}
}

void DoForwardSlideAlignment (run_params p, int i, int& firstpos, int&isok, vector<char> qual, rseq refseq, vector<rd>& data) {
	string ref=refseq.seq;
	SlideAlign(data[i].seq,ref,0.95,isok,i,data,1);
	//Edit sequence to remove low quality nucleotides
	int incl=0;
	for (int k=0;k<data[i].seq.size();k++) {
		int keep=0;
		for (int j=p.min_qual;j<=p.max_qual;j++) {
			if (data[i].qual[k]==qual[j]) {
				keep=1;
				break;
			}
		}
		if (keep==0) {
			data[i].seq.replace(k,1,"-");
		} else {
			incl++;
		}
	}
	if (incl<p.min_rlen) {
		data[i].inc=0;
	}
}


//Aligns the reverse complement of the sequence
void DoReverseSlideAlignment (run_params p, int i, int& firstpos, int&isok, vector<char> qual, rseq refseq, vector<rd>& data) {
	string ref=refseq.seq;
	data[i].revseq=RevTr(data[i].seq);
//	cout << "Rev " << data[i].revseq << "\n";
//	cout << "Q " << data[i].qual << "\n";
	data[i].revqual=RevString(data[i].qual);
	SlideAlign(data[i].revseq,ref,0.95,isok,i,data,-1);
	//Edit sequence to remove low quality nucleotides
	int incl=0;
	for (int k=0;k<data[i].revseq.size();k++) {
		int keep=0;
		for (int j=p.min_qual;j<=p.max_qual;j++) {
			if (data[i].revqual[k]==qual[j]) {
				keep=1;
				break;
			}
		}
		if (keep==0) {
			data[i].revseq.replace(k,1,"-");
		} else {
			incl++;
		}
	}
	if (incl<p.min_rlen) {
		data[i].inc=0;
	}
}

//Simple alignment algorithm - ignores indels
void SlideAlign (string read, string ref, double id, int& isok, int i, vector<rd>& data, int dir) {
	int rd_len=read.length();
	int rf_len=ref.length();

	int min_id=floor(id*rd_len);
	int max_err=rd_len-min_id;
	int max_fit=0;
	int best_a=-1e5;
	for (int step=-rd_len+1;step<rf_len;step++) {
		int match=0;
		int err=0;
		int over=0;
		for (int j=0;j<rd_len;j++) {
			if (j+step<0||j+step>rf_len) {
				over++;
			} else {
				if (read.substr(j,1)==ref.substr(j+step,1)) {
					match++;
				} else {
					err++;
					if (err>max_err) {
						break;
					}
				}
			}
		}
		if (match>max_fit) {
			best_a=step;
			max_fit=match;
		}
	}
	string out;
	if (best_a!=-1e5) {
		if (best_a<0) {
			out=read.substr(-best_a,rd_len+best_a); //Trim from start
			if (max_fit>id*out.length()) {
				data[i].alpos=0;
				if (dir==1) {
					data[i].inc=1;
					data[i].seq=out;
					data[i].qual=data[i].qual.substr(-best_a,rd_len+best_a);
				} else {
					data[i].inc=-1;
					data[i].revseq=out;
					data[i].revqual=data[i].revqual.substr(-best_a,rd_len+best_a);

				}
				isok=1;
			}
	
		} else {
			string temp (best_a, ' ');
			if (dir==1) {
				if (best_a+rd_len>rf_len) {
					read=read.substr(0,rf_len-best_a);
				}
				if (max_fit>id*read.length()) {
					if (best_a+rd_len>rf_len) {
						data[i].qual=data[i].qual.substr(0,rf_len-best_a);
					}
					out=temp+read;
					data[i].alpos=best_a;
					data[i].inc=1;
					data[i].seq=read;
					isok=1;
				}
			} else {
				if (best_a+rd_len>rf_len) {
					read=read.substr(0,rf_len-best_a);
				}
				if (max_fit>id*read.length()) {
					if (best_a+rd_len>rf_len) {
						data[i].revqual=data[i].revqual.substr(0,rf_len-best_a);
					}
					out=temp+read;
					data[i].alpos=best_a;
					data[i].inc=-1;
					data[i].revseq=read;
					isok=1;
				}
			}
		}
	}
}
