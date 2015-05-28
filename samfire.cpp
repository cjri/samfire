#include <iostream>
#include <vector>
#include <string>
#include <sstream>


using namespace std;

#include "shared.h"

int main(int argc, char *argv[]) {

	//Load reference sequences
	ifstream ref_file;
	ref_file.open(argv[1]);
	vector<refseq> refs;
	GetRefSeqs (ref_file,refs);

	string ref=refs[0].seq;
	
	//Load sam file data
	ifstream sam_file;
	sam_file.open(argv[2]);
	vector<rd> data;
	int s_length=0;
	ReadSamFile (sam_file,s_length,data,refs);
	
	vector<char> qual;
	makequal(qual);

	
	int min_qual=30; //Minimum sequence quality
	int max_qual=41; //Maximum sequence quality
	//Trim from both ends - find the longest...
	int q=0;

	for (int i=1;i<s_length;i++) {
		cout << "Sequence " << i << "\n";
		data[i].inc=0;
		if (data[i].refno!=-1) {
			int isok=1;
			int firstpos=-1;
			q=findqual2(i,min_qual,max_qual,qual,data);
		//	cout << "Quality " << q << "\n";
			if (q>=min_qual) {  //Minimum sequence quality
				DoForwardAlignment (i,min_qual,max_qual,firstpos,isok,qual,refs,data);
				if (isok==0) {
					isok=1;
					int firstpos=-1;
					DoReverseAlignment (i,min_qual,max_qual,firstpos,isok,qual,refs,data);
				}
			}
		}
	}
	
	//Construct individual allele frequencies
	vector<nuc> ref_counts;
	//for (int r=0;r<=7;r++) { //Reference sequences
		int r=0;
		
		nuc counts;
		for (int i=0;i<refs[r].size;i++) {
			counts.nN.push_back(0);
			counts.nA.push_back(0);
			counts.nC.push_back(0);
			counts.nG.push_back(0);
			counts.nT.push_back(0);
		}
		ref_counts.push_back(counts);
	//}

	ofstream seq_file1;
	ofstream inc_file;
	seq_file1.open("Sequences.out");
	inc_file.open("Included.out");
	for (int i=0;i<=s_length;i++) {
		if (data[i].inc==1) {
			int ref=data[i].refno;
			if (ref==0) {
				seq_file1 << data[i].alpos << " " << data[i].seq << "\n";
				inc_file << i << "\n";
			}
		}
		if (data[i].inc==-1) {
			int ref=data[i].refno;
			if (ref==0) {
				seq_file1 << data[i].alpos << " " << data[i].revseq << "\n";
				inc_file << i << "\n";
			}
		}
		if (data[i].inc==1) {
			int ref=data[i].refno;
			for (int j=0;j<data[i].seq.size();j++) {
				
				if (data[i].seq[j]=='A') {
					ref_counts[ref].nA[data[i].alpos+j]++;
					ref_counts[ref].nN[data[i].alpos+j]++;
				}
				if (data[i].seq[j]=='C') {
					ref_counts[ref].nC[data[i].alpos+j]++;
					ref_counts[ref].nN[data[i].alpos+j]++;
				}
				if (data[i].seq[j]=='G') {
					ref_counts[ref].nG[data[i].alpos+j]++;
					ref_counts[ref].nN[data[i].alpos+j]++;
				}
				if (data[i].seq[j]=='T') {
					ref_counts[ref].nT[data[i].alpos+j]++;
					ref_counts[ref].nN[data[i].alpos+j]++;
				}
			}
		}
		if (data[i].inc==-1) {
			int ref=data[i].refno;
			for (int j=0;j<data[i].revseq.size();j++) {
				
				if (data[i].revseq[j]=='A') {
					ref_counts[ref].nA[data[i].alpos+j]++;
					ref_counts[ref].nN[data[i].alpos+j]++;
				}
				if (data[i].revseq[j]=='C') {
					ref_counts[ref].nC[data[i].alpos+j]++;
					ref_counts[ref].nN[data[i].alpos+j]++;
				}
				if (data[i].revseq[j]=='G') {
					ref_counts[ref].nG[data[i].alpos+j]++;
					ref_counts[ref].nN[data[i].alpos+j]++;
				}
				if (data[i].revseq[j]=='T') {
					ref_counts[ref].nT[data[i].alpos+j]++;
					ref_counts[ref].nN[data[i].alpos+j]++;
				}
			}
		}
	}
	ofstream all_file1;
	all_file1.open("Allele_frequencies.out");

	for (int i=0;i<refs[r].size;i++) {
		all_file1 << i << " " <<  ref_counts[r].nN[i] << " " << ref_counts[r].nA[i] << " " << ref_counts[r].nC[i] << " " << ref_counts[r].nG[i] << " " << ref_counts[r].nT[i] << "\n";
	}
	
	return EXIT_SUCCESS;
}
							
 
