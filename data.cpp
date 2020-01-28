#include "shared.h"
#include "data.h"
#include <iostream>
#include <string>

void SetupAATrans (vector< vector<char> >& cod, vector<char>& aa) {
	vector<char> a;
	for (int i=0;i<3;i++) {
		a.push_back('A');
	}
	for (int i=0;i<64;i++){
		cod.push_back(a);
	}
	for (int i=0;i<64;i++) {
		if (i%4==1) {
			cod[i][2]='C';
		}
		if (i%4==2) {
			cod[i][2]='G';
		}
		if (i%4==3) {
			cod[i][2]='T';
		}
		if ((i/4)%4==1) {
			cod[i][1]='C';
		}
		if ((i/4)%4==2) {
			cod[i][1]='G';
		}
		if ((i/4)%4==3) {
			cod[i][1]='T';
		}
		if ((i/16)%4==1) {
			cod[i][0]='C';
		}
		if ((i/16)%4==2) {
			cod[i][0]='G';
		}
		if ((i/16)%4==3) {
			cod[i][0]='T';
		}
	}
	for (int i=0;i<64;i++) {
		aa.push_back('A');
	}
	aa[0]='K';
	aa[1]='N';
	aa[2]='K';
	aa[3]='N';
	aa[4]='T';
	aa[5]='T';
	aa[6]='T';
	aa[7]='T';
	aa[8]='R';
	aa[9]='S';
	aa[10]='R';
	aa[11]='S';
	aa[12]='I';
	aa[13]='I';
	aa[14]='M';
	aa[15]='I';
	aa[16]='Q';
	aa[17]='H';
	aa[18]='Q';
	aa[19]='H';
	aa[20]='P';
	aa[21]='P';
	aa[22]='P';
	aa[23]='P';
	aa[24]='R';
	aa[25]='R';
	aa[26]='R';
	aa[27]='R';
	aa[28]='L';
	aa[29]='L';
	aa[30]='L';
	aa[31]='L';
	aa[32]='E';
	aa[33]='D';
	aa[34]='E';
	aa[35]='D';
	//Defined already
	aa[40]='G';
	aa[41]='G';
	aa[42]='G';
	aa[43]='G';
	aa[44]='V';
	aa[45]='V';
	aa[46]='V';
	aa[47]='V';
	aa[48]='X';
	aa[49]='Y';
	aa[50]='X';
	aa[51]='Y';
	aa[52]='S';
	aa[53]='S';
	aa[54]='S';
	aa[55]='S';
	aa[56]='X';
	aa[57]='C';
	aa[58]='W';
	aa[59]='C';
	aa[60]='L';
	aa[61]='F';
	aa[62]='L';
	aa[63]='F';
}








