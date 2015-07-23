//Program to process .sam file data

#include <iostream>
#include <vector>
#include <string>
#include <sstream>


using namespace std;

#include "shared.h"

int main(int argc, char *argv[]) {
	
	//Load sam file data
	ifstream sam_file;
	sam_file.open(argv[1]);
	vector<pr> data;
	int s_length=0;
	//	sam_file.open("../../alignments/4317/4317_d_4_AEngland1952009.sam");
	ReadSamFile2 (sam_file,s_length,data);
	data.pop_back();
	data.pop_back();
	
	for (unsigned int i=0;i<data.size();i++) {
		for (unsigned int j=i+1;j<data.size();j++) {
			if ((i!=j)&&(data[i].s1==data[j].s1)) {
				string t1=data[i].s2;
				string t2=data[j].s2;
				string del="_";
				t1=t1.substr(0,t1.find(del));
				t2=t2.substr(0,t2.find(del));
				if (t1==t2) {
					cout << i << " " << j << "\n";
				}
			}
		}
	}
	
	return 0;
	
	
	
	return EXIT_SUCCESS;
}
							
 
