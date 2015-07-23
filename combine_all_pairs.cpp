//Program to process .sam file data

#include <iostream>
#include <vector>
#include <string>
#include <sstream>


using namespace std;

#include "shared.h"

int main(int argc, char *argv[]) {
	
	vector<det> seqs;
	
	//Load sequence file
	ifstream seq_file;
	seq_file.open(argv[1]);
	int i;
	string s;
	while (!seq_file.eof()) {
		det d;
		seq_file >> i;
		seq_file >> s;
		d.seq=s;
		d.start=i;
		seqs.push_back(d);
	}
	
	
	//Load inclusion listing
	ifstream inc_file;
	inc_file.open(argv[2]);
	int c=0;
	while (!inc_file.eof()) {
		inc_file >> i;
		seqs[c].sam_ref=i;
		//cout << seqs[c].start << " " << seqs[c].sam_ref << " " << seqs[c].seq << "\n";
		c++;
	}
	

	//Load pairs file
	vector<par> pairs;
	ifstream pair_file;
	pair_file.open(argv[3]);
	while (!pair_file.eof()) {
		par p;
		pair_file >> i;
		p.i1=i;
		pair_file >> i;
		p.i2=i;
		pairs.push_back(p);
	}
	
	for (unsigned int i=0;i<seqs.size();i++) {
		//cout << i << "\n";
		int j1=-1;
		int j2=-1;
		for (unsigned int j=0;j<seqs.size();j++) {
			if (seqs[j].sam_ref==pairs[i].i1) {
				j1=j;
			}
			if (seqs[j].sam_ref==pairs[i].i2) {
				j2=j;
			}
		}
		//Combine reads
		if ((j1>-1)&&(j2>-1)) {
		//	cout << "First " << seqs[j1].start << " " << seqs[j1].seq << " " << seqs[j1].seq.length() << "\n";
		//	cout << "Second " << seqs[j2].start << " " << seqs[j2].seq << " " << seqs[j2].seq.length() << "\n";
			string newstr;
			if (seqs[j1].start>seqs[j2].start) {
				int temp=j1;
				j1=j2;
				j2=temp;
			}
			int pos1=0;
			int pos2=0;
			if (seqs[j2].start < seqs[j1].start+seqs[j1].seq.length()) {
			//	pos2++;
			}
			
			int loc=seqs[j1].start;
			while (loc<seqs[j2].start+seqs[j2].seq.length()) {
				//cout << loc << " " <<seqs[j2].start+seqs[j2].seq.length() << "\n";
				if (loc<seqs[j1].start+seqs[j1].seq.length()) {
				//	cout << "Here1 " << loc << " " << pos1 << " " << seqs[j1].seq[pos1] << "\n";
					newstr.push_back(seqs[j1].seq[pos1]);
					pos1++;
					loc++;
					if (loc>seqs[j2].start) {
						pos2++;
					}
				} else if (loc<seqs[j2].start){
				//	cout << "Here2\n";
					newstr.push_back('-');
					loc++;
				} else {
					//cout << "Here3 " << loc << " " << pos2 << " " << seqs[j2].seq[pos2] << " " << seqs[j2].seq[pos2-1] << "\n";
					newstr.push_back(seqs[j2].seq[pos2]);
					pos2++;
					loc++;
				}
			}
			cout << seqs[j1].start << " " << newstr << "\n";
		}
		
	}
	
	return 0;

	
	return EXIT_SUCCESS;
}
							
 
