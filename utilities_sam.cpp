#include "shared.h"
#include <iostream>
#include <string>

void GetRefSeqs (ifstream& ref_file, vector<refseq>& refs) {
	cout << "Getting reference sequence...\n";
	string ref;
	string s;
	int rsize;
	for (int i=1;i<=1;i++) {
		refseq r;
		ref_file >> s;
		if (s[0]=='>') {
			string rname=s.substr(1,27);
			r.name=rname;
			//cout << rname << "\n";
			getline(ref_file,s);
		}
		ref_file >> ref;
		r.seq=ref;
		//cout << ref << "\n";
		rsize=ref.size();
		//cout << rsize << "\n";
		r.size=rsize;
		getline(ref_file,s);
		refs.push_back(r);
	}
}

void ReadSamFile (ifstream& sam_file, int& s_length, vector<rd>& data, vector<refseq> refs) {
	cout << "Reading .sam file... \n";
	string line;
	const string temp="@SQ";
	string lstore;
	while (!sam_file.eof()) {
		//	cout << "Reading " << j << "\n";
		rd r;
		for (int i=1;i<=11;i++) {
			sam_file >> line;
			if (i==1) {
				lstore=line;
			}
			if (i==3) {
				r.ref=line;
				r.refno=-1;
				for (unsigned int k=0;k<refs.size();k++) {
					if (line.compare(refs[k].name)==0) {
						r.refno=k;
						//cout << k << " ";
					}
				}
				//cout << r.ref << "\n";
			}
			if (i==4) {
				r.alpos=atoi(line.c_str());
				//				cout << r.ref << "\n";
			}
			if (i==10) {
				r.seq=line;
				//cout << r.seq << "\n";
			}
			if (i==11) {
				r.qual=line;
				//cout << r.qual << "\n";
			}
			
			//			cout << line << "\n";
		}
		r.rev=0;
		r.inc=0;
		if (lstore.compare(temp)!=0) {
			s_length++;
			data.push_back(r);
		}
		//		cout << data.size() << "\n";
		getline(sam_file,line); //Remainder of line
	}
	cout << "Number of reads " << s_length << "\n";
}

void ReadSamFile2 (ifstream& sam_file, int& s_length, vector<pr>& data) {
	//cout << "Reading .sam file... \n";
	string line;
	const string temp="@SQ";
	string lstore;
	while (!sam_file.eof()) {
		//	cout << "Reading " << j << "\n";
		pr r;
		for (int i=1;i<=11;i++) {
			sam_file >> line;
		//	cout << i << " " << line << "\n";
			if (i==1) {
				lstore=line;
				string del=":";
				string tok=lstore.substr(0,lstore.find(del));
				lstore.erase(0,lstore.find(del)+1);
				tok=lstore.substr(0,lstore.find(del));
				lstore.erase(0,lstore.find(del)+1);
				tok=lstore.substr(0,lstore.find(del));
				lstore.erase(0,lstore.find(del)+1);
				tok=lstore.substr(0,lstore.find(del));
				lstore.erase(0,lstore.find(del)+1);
				tok=lstore.substr(0,lstore.find(del));
				lstore.erase(0,lstore.find(del)+1);
				tok=lstore.substr(0,lstore.find(del));
				r.s1=tok;
				lstore.erase(0,lstore.find(del)+1);
				tok=lstore.substr(0,lstore.find(del));
				r.s2=tok;
				//cout << r.s1 << " " << r.s2 << "\n";
			}
		}
		if (lstore.compare(temp)!=0) {
			s_length++;
			data.push_back(r);
		}
		//		cout << data.size() << "\n";
		getline(sam_file,line); //Remainder of line
	}
	//cout << "Number of reads " << s_length << "\n";
}



void DoForwardAlignment (int i, int min_qual, int max_qual, int& firstpos, int&isok, vector<char> qual, vector<refseq> refs, vector<rd>& data) {
	ofstream temp_fileo;
	ifstream temp_filei;
//	cout << "Here\n";
	temp_fileo.open("dat.in");
	string ref=refs[data[i].refno].seq;
	temp_fileo << ">Ref " << "\n";
	temp_fileo << ref << "\n";
	temp_fileo << ">Seq1\n";
	temp_fileo << data[i].seq << "\n";
	temp_fileo.close();
	system("../Code/muscle3.8.31_i86darwin64 -in dat.in -out dat.out -clw -quiet"); 	temp_filei.open("dat.out");
	//Read the Muscle output file
	string line;
	string aliseq ("");
	string refseq ("");
	string star;
	while (!temp_filei.eof()) {
		if (line.compare("Ref")!=0) {
			temp_filei >> line;
		}
		if (line.compare("Ref")==0) {
			temp_filei >> line;
			refseq += line;
			temp_filei >> line;
			temp_filei >> line;
			aliseq += line;
			temp_filei >> line;
			while (line.compare(0,1,"*")==0) {
				if (temp_filei.eof()) {
					break;
				}
				star += line;
				temp_filei >> line;
			}
		}
	}
//	cout << aliseq << "\n";
//	cout << star << "\n";
//	cout << "Number of stars is " << star.length() << "\n";
//	cout << "Number of nucleotides is " << data[i].seq.size() << "\n";


	int ins=0;
	int ins_tot=0;
	MeasureInsertions (ins,aliseq);
	if (ins<=1) {
		ins_tot=ins_tot+ins;
	} else {
		ins_tot=10;
	}
	MeasureInsertions (ins,refseq);
	if (ins<=1) {
		ins_tot=ins_tot+ins;
	} else {
		ins_tot=10;
	}
	//cout << "Total insertions " << ins_tot << "\n";
	int refin_s=0;
	OverlapStart (refin_s,refseq);
	data[i].seq.erase(0,refin_s);

	int refin_e=0;
	OverlapEnd (refin_e,refseq);
	data[i].seq.erase(data[i].seq.size()-refin_e,refin_e);

	
	//cout << "Modified number of nucleotides " << data[i].seq.size() << "\n";
	
	CheckOK(ins_tot,isok,star,data[i].seq);
	
	if (isok==1) {
		//Find insertions in the sequence being aligned and remove them
		vector<int> remove;
		int ii=-1;
		for (unsigned int j=0;j<refseq.size();j++) {  //Identify positions in the reference sequence where there are insertions
			if (ii==-1&&aliseq.compare(j,1,"-")!=0) {
				ii=j;
			}
			if (refseq.compare(j,1,"-")==0) {
				int k=j-refin_s;
				if (k>=0) {  //Don't consider positions before the end of the alignment
					if (k-ii<data[i].seq.size()) {
						remove.push_back(k);
	//					cout << "Adding " << k << " " << aliseq[j] << " " << ii << "\n";
					}
				}
			}
		}
		for (unsigned int j=0;j<remove.size();j++) {
	//		cout << remove[j] << "\n";
			if (remove[j]-ii<data[i].seq.size()) {  //Don't consider positions after the end of the alignment
	//			cout << "Removing " << remove[j]-ii << " " << data[i].seq.size() << "\n";
				data[i].seq.erase(remove[j]-ii,1);  //Remove insertion from the aligned data
			}
		}

		
		
		
		//Find deletions in the sequence being aligned and replace them with "-"
		vector<int> replace;
		ii=0;
		for (unsigned int j=0;j<=aliseq.size();j++) {
			if (aliseq.compare(j,1,"-")!=0) {
				ii=j; //First position in the aligned data that contains a nucleotide
				break;
			}
		}
//		cout << "ii " << ii << "\n";
		for (unsigned int j=ii;j<ii+data[i].seq.size();j++) {
			if (aliseq.compare(j,1,"-")==0) {
				int k=j-refin_s;
				if (k-ii>=0) {  //Don't consider positions before the end of the alignment
					replace.push_back(k-ii);
	//				cout << "Deletion " << k-ii << "\n";
				}
			}
		}
		for (unsigned int j=0;j<replace.size();j++) {
			data[i].seq.insert(replace[j],"-");
		}
		
	//	cout << "Include 1 " << data[i].seq << "\n";
		
		//Edit sequence to remove low quality nucleotides
		for (int k=0;k<data[i].seq.size();k++) {
			int keep=0;
			for (int j=min_qual;j<=max_qual;j++) {
				if (data[i].qual[k]==qual[j]) {
					keep=1;
					break;
				}
			}
			if (keep==0) {
//				cout << "Replace " <<
				data[i].seq.replace(k,1,"-");
			}
		}
		
	//	cout << "Include 2 " << data[i].seq << "\n";

		data[i].inc=1; //Forward
		//Find position of first nucleotide
		FindFirstPos (firstpos,refin_s,refseq,aliseq);
	//	cout << "Position " << firstpos << "\n";
		data[i].alpos=firstpos;
	//	cout << data[i].refno << " " << data[i].alpos << " " << data[i].seq << "\n";
	}
	temp_filei.close();
}

void DoReverseAlignment (int i, int min_qual, int max_qual, int& firstpos, int&isok, vector<char> qual, vector<refseq> refs, vector<rd>& data) {
	ofstream temp_filero;
	ifstream temp_fileri;
//	cout << "Reversing...\n";
	data[i].revseq=RevTr(data[i].seq);
	temp_filero.open("dat_r.in");
	string ref=refs[data[i].refno].seq;
	temp_filero << ">Ref " << "\n";
	temp_filero << ref << "\n";
	temp_filero << ">Seq1\n";
	temp_filero << data[i].revseq << "\n";
	temp_filero.close();
	system("../Code/muscle3.8.31_i86darwin64 -in dat_r.in -out dat_r.out -clw -quiet");
	temp_fileri.open("dat_r.out");
	//Read the Muscle output file
	string line;
	string aliseq ("");
	string refseq ("");
	string star;

	while (!temp_fileri.eof()) {
		if (line.compare("Ref")!=0) {
			temp_fileri >> line;
		//	cout << line << "\n";
		}
		if (line.compare("Ref")==0) {
			temp_fileri >> line;
		//	cout << line << "\n";
			refseq += line;
			temp_fileri >> line;
		//	cout << line << "\n";
			temp_fileri >> line;
		//	cout << line << "\n";
			aliseq += line;
			temp_fileri >> line;
			while (line.compare(0,1,"*")==0) {
				if (temp_fileri.eof()) {
					break;
				}
				star += line;
				temp_fileri >> line;
				//cout << line << "\n";
			}
		}
	}

	//	cout << aliseq << "\n";
	//	cout << star << "\n";
	//	cout << "Number of stars is " << star.length() << "\n";
	//	cout << "Number of nucleotides is " << data[i].seq.size() << "\n";
	
	int ins=0;
	int ins_tot=0;
	MeasureInsertions (ins,aliseq);
	if (ins<=1) {
		ins_tot=ins_tot+ins;
	} else {
		ins_tot=10;
	}
	MeasureInsertions (ins,refseq);
	if (ins<=1) {
		ins_tot=ins_tot+ins;
	} else {
		ins_tot=10;
	}
	//cout << "Total insertions " << ins_tot << "\n";
	int refin_s=0;
	OverlapStart (refin_s,refseq);
	data[i].revseq.erase(0,refin_s);
	
	int refin_e=0;
	OverlapEnd (refin_e,refseq);
	data[i].revseq.erase(data[i].revseq.size()-refin_e,refin_e);
	
	//cout << "Modified number of nucleotides " << data[i].revseq.size() << "\n";

	CheckOK(ins_tot,isok,star,data[i].revseq);
	
	if (isok==1) {
		//Find insertions in the sequence being aligned
		vector<int> remove;
		int ii=-1;
		for (unsigned int j=0;j<refseq.size();j++) {  //Identify positions in the reference sequence where there are insertions
			if (ii==-1&&aliseq.compare(j,1,"-")!=0) {
				ii=j;
			}
			if (refseq.compare(j,1,"-")==0) {
				int k=j-refin_s;
				if (k>=0) {  //Don't consider positions before the end of the alignment
					if (k-ii<data[i].revseq.size()) {//Don't consider positions after the end of the alignment
						remove.push_back(k);
			//			cout << "Insertion " << k << " " << aliseq[i] << " " << ii << "\n";
					}
				}
			}
		}
		
		for (unsigned int k=0;k<remove.size();k++) {
		//		cout << "Removing " << remove[k]-ii << " " << data[i].revseq.size() << "\n";
				data[i].revseq.erase(remove[k]-ii,1);  //Remove insertion from the aligned data
			//}
		}
		
		

		
		vector<int> replace;
		ii=0;
		for (unsigned int j=0;j<=aliseq.size();j++) {
			if (aliseq.compare(j,1,"-")!=0) {
				ii=j;
				break;
			}
		}
		for (unsigned int j=ii;j<ii+data[i].revseq.size();j++) {
			if (aliseq.compare(j,1,"-")==0) {
				int k=j-refin_s;
				if (k-ii>=0) {  //Don't consider positions before the end of the alignment
					replace.push_back(k-ii);
		//			cout << "Deletion " << k-ii << "\n";
				}
			}
		}
		for (unsigned int j=0;j<replace.size();j++) {
			data[i].revseq.insert(replace[j],"-");
		}

		
		
		
		//Edit sequence to remove low quality nucleotides
		for (int k=0;k<data[i].revseq.size();k++) {
			int keep=0;
			for (int j=min_qual;j<=max_qual;j++) {
				if (data[i].qual[k]==qual[j]) {
					keep=1;
					break;
				}
			}
			if (keep==0) {
				//				cout << "Replace " <<
				data[i].revseq.replace(k,1,"-");
			}
		}
		data[i].inc=-1; //Forward
		//Find position of first nucleotide
		FindFirstPos (firstpos,refin_s,refseq,aliseq);
	//	cout << "Position " << firstpos << "\n";
		data[i].alpos=firstpos;
	}
	temp_fileri.close();
}


void MeasureInsertions (int& ins, string aliseq) {
	int swit=0;
	ins=-1;
	for (unsigned int i=0;i<aliseq.length();i++) {
		if (swit==0&&aliseq.compare(i,1,"-")!=0) { //Have to hit some sequence before counting
			swit=1;
			ins++;
		}
		if (swit==1&&aliseq.compare(i,1,"-")==0) {
			swit=0;
		}
		if (ins>1) {
			break;
		}
	}
//	cout << "Number of breaks in align = " << ins << "\n";
}



void OverlapStart (int& refin_s, string refseq) {
	//cout << "Refseq " << refseq << "\n";
	for (unsigned int i=0;i<refseq.length();i++) {
	//	cout << i << refseq[i] << "\n";
		if (refseq.compare(i,1,"-")!=0) {
			refin_s=i;
	//		cout << i << "\n";
			break;
		}
	}
//	cout << "Insertions at start of reference sequence " << refin_s << "\n";
}

void OverlapEnd (int& refin_e, string refseq) {
	for (int i=refseq.length()-1;i>=0;i--) {
		if (refseq.compare(i,1,"-")!=0) {
			refin_e=refseq.length()-i-1;
			//cout << refseq.length()-i-1 << "\n";
			break;
		}
	}
//	cout << "Insertions at end of reference sequence " << refin_e << "\n";
}

void CheckOK (int ins, int& isok, string star, string seq) {
	if (ins>=1) {
		isok=0;
	} else if (seq.size()<50) {
		isok=0;
	} else {
		double diff=(star.size()+0.)/(seq.size()+0.);
		if (diff<0.95) {
			isok=0;
		}
	}
}

void FindFirstPos (int& firstpos, int refin_s, string refseq, string aliseq) {
	for (unsigned int i=0;i<refseq.size();i++) {
		if (refseq.compare(i,1,"-")!=0&&aliseq.compare(i,1,"-")!=0) {
			firstpos=i;
			firstpos=firstpos-refin_s;
			break;
		}
	}
}




void makequal (vector<char>& qual) {
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
}

int findqual2 (int i, int min_qual, int max_qual, vector<char> qual, vector<rd>& data) {
	string q=data[i].qual;
	string s=data[i].seq;
//	cout << s << "\n";
//	cout << q << "\n";
	vector<int> qvec;
	for (int i=0;i<s.size();i++) {
		for (int j=0;j<=max_qual;j++) {
			if (q[i]==qual[j]) {
				qvec.push_back(j);
				//cout << j << " ";
				break;
			}
		}
	}
	
	int median = GetMedian(0,0,qvec);
//	cout << "Median " << median << " ";
	int qo=median;
	
	if (median<min_qual) {  //Edit sequence to get high quality part
		int s1=999;
		int s2=999;
		//Reduce from end
		for (int a=1;a<qvec.size()-30;a++) {
			median=GetMedian(a,0,qvec);
			if (median==min_qual) {
				s1=a;
				break;
			}
		}
		for (int a=1;a<qvec.size()-30;a++) {
			median=GetMedian(0,a,qvec);
			if (median==min_qual) {
				s2=a;
				break;
			}
		}
//		cout << " S vals " << s1 << " " << s2 << " ";
		int qs=q.size();
		
		if (s1<999&&s1<=s2) {
			q=q.substr(s1,qs-s1+1);
			s=s.substr(s1,qs-s1+1);
			data[i].qual=q;
			data[i].seq=s;
//			cout << "\n" << s << "\n";
//			cout << q << "\n";
			qo=min_qual;
		} else if (s2<999&&s2<s1) {
			q=q.substr(0,qs-s2+1);
			s=s.substr(0,qs-s2+1);
			data[i].qual=q;
			data[i].seq=s;
//			cout << "\n" << s << "\n";
//			cout << q << "\n";
			qo=min_qual;
		} else {
			qo=-1;
		}
	}
//	cout << "Quality " << qo << "\n";
	//Pre-scan number of min_quality reads
	if (qo>-1) {
		int nm=0;
		for (int i=0;i<s.size();i++) {
			for (int j=min_qual;j<=max_qual;j++) {
				if (q[i]==qual[j]) {
					nm++;
					break;
				}
			}
		}
//		cout << "Quality reads " << nm << "\n";
		if (nm<30) {
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
	//cout << median << "\n";
	return median;
}

int ScoreSim (string a, string b, int p, int rsize) {
	int s=a.size();
	int fulls=rsize;
	int sim=0;
	int l=0;

	for (int i=0;i<s;i++) {
		if ((p+i-1>=0)&&(p+i-1<fulls)) {
			l++;
	//		cout << a[i] << " " << b[p+i-1] << "\n";
			if (a[i]!=b[p+i-1]) {
				sim++;
			}
		}
	}
	
	if (l<50) {  //Minimum aligned read
		sim=50;
	}
//	cout << "Align2 " << l << " " << sim << "\n";
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
