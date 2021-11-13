#include "shared.h"
#include <iostream>
#include <string>

double DirichletMultiCalc (int N, double c, vector<int> obs, vector<double> inf, vector<double>& fact_store) {
	double bin=0;
//	cout << "C " << c << "\n";
//	cout << "N " << N << "\n";
//	cout << "Obs\n";
//	 for (unsigned int i=0;i<obs.size();i++) {
//	 cout << obs[i] << " ";
//	 }
//	 cout << "\n";
//	 cout << "Inf\n";
	 for (unsigned int i=0;i<inf.size();i++) {
		 if (inf[i]<1e-20) {
			 inf[i]=1e-20;
		 }
//		cout << inf[i] << " ";
	 }
//	 cout << "\n";

	if (N>0) {
		bin=fact_store[N];
		for (unsigned int i=0;i<obs.size();i++) {
			bin=bin-fact_store[obs[i]];
		}
		vector<double> alpha;
		for (unsigned int i=0;i<inf.size();i++) {
			alpha.push_back(c*inf[i]);
		}
		double a=0;
		for (unsigned int i=0;i<alpha.size();i++) {
			a=a+alpha[i];
			bin=bin-gsl_sf_lngamma(alpha[i]);
		}
		bin=bin+gsl_sf_lngamma(a);
		a=0;
		for (unsigned int i=0;i<alpha.size();i++) {
			double b=alpha[i]+obs[i];
			a=a+b;
			bin=bin+gsl_sf_lngamma(b);
		}
		bin=bin-gsl_sf_lngamma(a);
	} else {
		bin=0;
	}
	
//	cout << "L " << bin << "\n";
	return(bin);
}









