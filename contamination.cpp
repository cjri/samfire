#include "shared.h"
#include "contamination.h"
#include "io.h"
#include <iostream>
#include <string>
#include <sstream>

//Find length of reference sequence
void FindRefLen (int& len, string ref1, string ref2) {
	len=ref1.length();
	if (ref2.length()<ref1.length()) {
		len=ref2.length();
	}
	int max_l=0;
	for (int i=0;i<len;i++) {
		if (ref1.compare(i,1,"-")!=0&&ref2.compare(i,1,"-")!=0) {
			max_l=i;
		}
	}
	len=max_l;
	
}


//Get differences between reference sequences or an alternative list of sequence positions
void FindPolymorphisms (run_params p, int len, string ref1, string ref2, int& check, int& start_c, vector<int>& polys) {
	if (p.get_in==0) {
		cout << "Find differences between reference sequences...\n";
		for (int i=0;i<len;i++) {
			if (start_c<0&&ref1.compare(i,1,"-")!=0&&ref2.compare(i,1,"-")!=0) {
				start_c=i;
			}
			if (start_c>=0&&ref1.compare(i,1,ref2,i,1)!=0) {
				if (p.verb==1) {
					cout << i << " " << ref1[i] << " " << ref2[i] << "\n";
				}
//				cout << i << " " << ref1[i] << " " << ref2[i] << "\n";
				polys.push_back(i+1);
			}
		}
//		cout << start_c << "\n";
	} else {
		for (int i=0;i<len;i++) {
			if (start_c<0&&ref1.compare(i,1,"-")!=0&&ref2.compare(i,1,"-")!=0) {
				start_c=i;
			}
		}
		int check=0;
		ImportSLTLoci(p.in_file.c_str(),check,polys);
	}
}

void FindDifferences (run_params p, int len, int start_c, string ref1, string ref2, string c1, vector<int>& polys, vector<joined>& t_read, vector<int>& filter) {
	ofstream diff_file;
	diff_file.open(c1.c_str());
	diff_file << "#Read no.      #Nucleotides at specified sites      #Differences at specified sites (putative reference)      #Differences at specified sites (alternative reference)      #Nucleotides total      #Differences total (putative reference)      #Differences total (alternative reference)\n";
	
	for (int i=0;i<t_read.size();i++) {
		//Check if read potentially covers multiple sites
		int count=CountPolys(i,polys,t_read);
		if (count>=2) {
			int diff1=0;
			int diffp1=0;
			int nuc=0;
			int nucp=0;
			for (int j=0;j<t_read[i].seq.size();j++) { //Check - references are same length?
				int p_check=0;
				if (t_read[i].seq[j]!='-'&&j+t_read[i].alpos-1<len&&j+t_read[i].alpos-1>=start_c) {
					nuc++;
					for (int k=0;k<polys.size();k++) {
						if (j+t_read[i].alpos==polys[k]) { //Check -1 here?
							nucp++;
							p_check=1;
							break;
						}
					}
					if (t_read[i].seq[j]!=ref1[j+t_read[i].alpos-1]&&ref1.compare(j+t_read[i].alpos-1,0,"-")!=0) {
						diff1++;
						if (p_check==1) {
							diffp1++;
						}
					}
				}
			}
			int diff2=0;
			int diffp2=0;
			for (int j=0;j<t_read[i].seq.size();j++) {
				if (t_read[i].seq[j]!='-'&&j+t_read[i].alpos-1<len) {
					//cout << j+t_read[i].alpos-1 << " " << ref2.length() << "\n";

					if (j+t_read[i].alpos>start_c&&t_read[i].seq[j]!=ref2[j+t_read[i].alpos-1]&&ref2.compare(j+t_read[i].alpos-1,0,"-")!=0) {
						diff2++;
						for (int k=0;k<polys.size();k++) {
							if (j+t_read[i].alpos==polys[k]) {
								diffp2++;
							}
						}
					}
				}
			}
			if (p.verb==1) {
				diff_file << i << " " << nucp << " " << diffp1 << " " << diffp2 << " " << nuc << " " << diff1 << " " << diff2 << " " << "\n";
			} else {
				if (nucp>=2&&diffp1>diffp2) {
					diff_file << i << " " << nucp << " " << diffp1 << " " << diffp2 << " " << nuc << " " << diff1 << " " << diff2 << " " << "\n";
					filter.push_back(i);
				}
			}
		}
	}
	diff_file.close();
}

int CountPolys (int i, vector<int>& polys, vector<joined>& t_read) {
	int count=0;
	int end=t_read[i].alpos+t_read[i].seq.size();
	for (int j=0;j<polys.size();j++) {
		if (polys[j]>=t_read[i].alpos&&polys[j]<=end) {
			count++;
		}
	}
	return count;
}

