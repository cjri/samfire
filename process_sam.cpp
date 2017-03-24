//Program to process .sam file data

#include <iostream>
#include <vector>
#include <string>
#include <sstream>


using namespace std;

#include "shared.h"
#include "alignment.h"
#include "call_snps.h"
#include "call_mnps.h"
#include "ddups.h"
#include "findld.h"
#include "io.h"
#include "make_fullhaps.h"
#include "matchpairs.h"
#include "optimisation.h"
#include "utilities_sam.h"


int main(int argc, const char **argv) {

	run_params p;
	GetOptions (p,argc,argv);
	
	int seed=p.seed;
	gsl_rng_env_setup();
	gsl_rng *rgen = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set (rgen, seed);

	//Identify method to apply
	if (argc==1) {
		PrintInstructions (0);
		return 0;
	}
	
	string method=argv[1];
	
	if (method.compare("filter")==0) {
		//Import reference sequence
		ifstream ref_file;
		ref_file.open(p.ref.c_str());
		if (p.ref=="") {
			cout << "Error: No reference sequence specified.  Use e.g. --ref reference.fa.  Note that the sequence must be contained on a single line of the FASTA file\n";
			return 0;
		}
		rseq refseq;
		GetRefSeq (ref_file,refseq);
		string ref=refseq.seq;
		
		//Import list of .sam file names
		vector<string> sam_files;
		ImportSamFileNames(p,sam_files);
		if (sam_files.size()==0) {
			cout << "Error: No input .sam files specified\n";
			return 0;
		}
		
		//Construct quality reference vector
		vector<char> qual;
		makequal(p,qual);

		//Process files one by one...
		
		for (int i=0;i<sam_files.size();i++) {
			
			//Read .sam file
			alldat a;
			ImportSamFile (p,i,refseq,sam_files,a);
			
			//Detect and remove duplicate reads
			int found_pairs=0;
			if (p.ddup>=0) {
				DelDupSequences (p,found_pairs,qual,a);
			}

			if (p.almethod==1) {
				//Alignment of reads by data in .sam file
				if (p.verb==1) {
					cout << "Use .sam file alignment\n";
				}
				AlignSequencesSam (p,a.s_length,qual,refseq,a.data);
			} else if (p.almethod==2) {
				//Alignment of short reads by simple slide method
				AlignSequencesSlide (a.s_length,p,qual,refseq,a.data);
			}
			//Output unjoined reads
			OutputAlign(i,a.data);

			//Identify pairs
			cout << "Match paired-end reads " << p.pairs << "\n";
			if (p.pairs==1) {
				if (found_pairs==0) {
					MatchPairs(p,1,a.data);
				} else {
					CheckOldPairs(a.data);
				}
			}
			
			//Join pairs if they exist
			vector<joined> jseqs;
			JoinPairs (i,p,jseqs,a.data);

		}
		
		
		return 0;
		
	} else if (method.compare("sl_traj")==0) {
		
		//Get reference sequence
		ifstream ref_file;
		ref_file.open(p.ref.c_str());
		rseq refseq;
		GetRefSeq (ref_file,refseq);
		
		//Load sam file data
		vector<string> sam_files;
		if (p.no_sam>0) {
			for (int i=0;i<p.no_sam;i++) {
				string s;
				sam_files.push_back(s);
			}
		} else {
			ImportSamFileNames(p,sam_files);
		}
		
		//Get joined and aligned sequence reads
		vector<nuc> ref_counts;
		for (int i=0;i<sam_files.size();i++) {
			vector<joined> t_read;
			InputJnData (i,t_read);
		
			//Construct individual allele frequencies
			nuc r_count;
			CountNucs (p,refseq,t_read,r_count);

			//Output variant file
			OutputVarFile (i,refseq,r_count);
			ref_counts.push_back(r_count);
		}
		
		//Identify frequencies above a given cutoff
		vector<poly> polys;
		if (p.vs_ref==1) {
			CallPolymorphismsVsRef (p,refseq,ref_counts,polys);
		} else {
			CallPolymorphisms (p,refseq,ref_counts,polys);
		}
		
		//Construct single-locus tranjectories.  Contains locus, times, four nucleotide counts
		vector<str> sltrajs;
		ConstructSLTrajs(p,polys,ref_counts,sltrajs);
		
		//Output sltrajs data
		if (p.get_out==0) {
			p.out_file="Single_locus_trajectories.out";
		}
		
		OutputSLTData(p,p.out_file.c_str(),sltrajs);
		
		return 0;
		
	} else if (method.compare("consensus")==0) {
		
		vector<string> sam_files;
		if (p.no_sam>0) {
			for (int i=0;i<p.no_sam;i++) {
				string s;
				sam_files.push_back(s);
			}
		} else {
			ImportSamFileNames(p,sam_files);
		}
		
		//Import variant data
		vector< vector< vector<int> > > all_nucs;
		InputVariantData(sam_files,all_nucs);
		
		cout << all_nucs[0].size() << "\n";

		
		//Output consensus sequence data
		OutputGlobalConsensus(all_nucs);
		
	} else if (method.compare("ef_depth")==0) {
		
		vector<string> sam_files;
		if (p.no_sam>0) {
			for (int i=0;i<p.no_sam;i++) {
				string s;
				sam_files.push_back(s);
			}
		} else {
			ImportSamFileNames(p,sam_files);
		}

		//Read in Variant data
		vector< vector<int> > all_tots;
		GetVariantTotals(sam_files,all_tots);

		//Convert depth to effective depth
		ifstream csl_file;
		double c;
		csl_file.open("Csl.out");
		csl_file >> c;
		vector< vector<int> > all_etots;
		for (int i=0;i<all_tots.size();i++) {
			vector<int> et;
			for (int j=0;j<all_tots[i].size();j++) {
				double cd=(all_tots[i][j]*(1+c))/(all_tots[i][j]+c);
				int cdi=floor(cd);
				et.push_back(cdi);
			}
			cout << i << " " << et.size() << "\n";
			all_etots.push_back(et);
		}
		
		
		for (int i=0;i<sam_files.size();i++) {
			ofstream dep_file;
			ostringstream convert;
			convert << i;
			string temp=convert.str();
			string name = "Depths"+temp+".out";
			dep_file.open(name.c_str());
			for (int j=0;j<all_tots[i].size();j++) {
				dep_file << j << " " << all_tots[i][j] << " " <<  all_etots[i][j] << "\n";
			}
		}
	
		
	} else if (method.compare("sl_noise")==0) {

		vector<str> sltrajs;
		if (p.get_in==0) {
			p.in_file="Single_locus_trajectories.out";
		}
		
		ImportSLTData(p,sltrajs);  //N.B. Times are encoded in the trajectories

		SLTFreqs(p,sltrajs); //Convert observations to frequencies

		//Remove trajectories that move by more than a cutoff amount per day
		cout << "Number of recorded time points " << sltrajs[0].times.size() << "\n";
		
		if (sltrajs[0].times.size()>1) {
			FilterSLTrajs(p,sltrajs);
		}
		
		//Remove non-polymorphic time-points from trajectories
		FilterSLTrajs2(p,sltrajs);
		
		//Calculate mean frequencies of trajectories - approximate fit under assumption of neutrality
		SLTMeanFreqs(p,sltrajs);

		//Calculate vector of log factorials
		vector<double> fact_store;
		GetFactVectorSL(sltrajs,fact_store);
		
		double Csl_opt=0;
		OptimiseCsl(p,Csl_opt,sltrajs,fact_store,rgen);
		
		ofstream csl_file;
		csl_file.open("Csl.out");
		csl_file << Csl_opt << "\n";
		csl_file.close();
		
	} else if (method.compare("sl_neutrality")==0) {
		
		vector<str> sltrajs;
		if (p.get_in==0) {
			p.in_file="Single_locus_trajectories.out";
		}
		ImportSLTData(p,sltrajs);
		
		SLTFreqs(p,sltrajs);
		
		//Calculate mean frequencies - best fit under the assumption of neutrality with no drift
		SLTMeanFreqs(p,sltrajs);

		//Calculate fact_store vector
		vector<double> fact_store;
		GetFactVectorSL(sltrajs,fact_store);
		
		double Csl_opt=0;
		if (p.specify_csl!=0) {
			Csl_opt=p.specify_csl;
		} else {
			ifstream csl_file;
			csl_file.open("Csl.out");
			csl_file >> Csl_opt;
		}
		
		//Do single-locus calculations of movement; identify selected positions
		PotentialNonNeutrality (p,Csl_opt,sltrajs,fact_store,rgen);

		//Output potentially non-neutral trajectories
		if (p.get_out==0) {
			p.out_file="Single_locus_trajectories_sl.out";
		}
		OutputSLTData(p,p.out_file.c_str(),sltrajs);

		return 0;
		
	} else if (method.compare("call_ml")==0) {  //Call multi-locus polymorphisms from a single dataset
		
		cout << p.get_in << "\n";
		for (int i=0;i<argc;i++) {
			cout << argv[i] << "\n";
		}
		
		//Read polymorphism data
		if (p.get_in==0) {
			p.in_file="Single_locus_trajectories_sl.out";
		}
		
		//Read in polymorphic loci
		vector<int> polys;
		int check=0;
		ImportSLTLoci(p.in_file.c_str(),check,polys);
		if (check==1) {
			return 0;
		}
		
		//Read in times
		vector<int> times;
		ImportTimeData(times);
		int n_times=times.size();
		
		vector< vector< vector<mpoly> > > m_polys;
		vector<vector< vector<int> > > all_l_combs;
		int max_l_store=0;
		for (int i=0;i<n_times;i++) {  //Problem exists here - require the same subsets of partial haplotypes at each time point
			vector<joined> t_read;
			//Read in files
			InputJnData (i,t_read);
			//Find maximum read length
			int max_r_size=0;
			MaxRdLength (p,max_r_size,t_read);
			//Find maximum number of loci potentially spanned by a read
			int max_l=0;
			MaxReadSpan(p,max_r_size,max_l,polys);
			//Identify nucleotides for each polymorphism for each paired end read
			vector <vector<char> > m_vars;
			MapSeqToSNPLoci(p,n_times,polys,t_read,m_vars);

			//Make collection of all types; indices of loci that might be involved in a multi-locus variant call
			vector< vector<int> > l_combs;
			CallReadTypesConsecutive(p,max_l,max_r_size,polys,l_combs);  //Note for longer genomes - could consider limit on #loci.  Are full haplotypes required?

			for (int j=0;j<l_combs.size();j++) {
				if (l_combs[j].size()>max_l_store) {
					max_l_store=l_combs[j].size();
				}
			}
			all_l_combs.push_back(l_combs);
			
			vector< vector<mpoly> > m_pol;
			CallMNPs2 (p,l_combs,m_vars,m_pol);
			m_polys.push_back(m_pol);
		}
		
		//Order polymorphisms by loci spanned
		
		cout << "Find Max Loc\n";
		int max_loc=0;
		FindMaxLoc(max_loc,all_l_combs);
		vector< vector<int> > l_vec;
		ConstructLocVec (p,max_loc,max_l_store,all_l_combs,l_vec);
		ConstructMPoly (l_vec,all_l_combs,m_polys);

		//Construct full haplotypes
		//Filter data at time-resolved level: minimum frequency and #occurrences within each dataset
		vector< vector< vector<mpoly> > > m_polys_f;
		FilterMNPs(p,m_polys,m_polys_f);
		
		//Combine data across time points
		vector< vector<mpoly> > c_m_polys;
		CombineMNPsTime(p,m_polys_f,c_m_polys);
		
		vector< vector<char> > haps;
		MakeInitialHaps (p,c_m_polys,haps);
		
		if (p.full_haps==0) {
			GenerateOutputFilesNoHaps(p,polys,m_polys_f,c_m_polys);
		} else {
			
			//Identify haplotypes either from raw data or from external file
			vector< vector<char> > haps;
			if (p.readhap==0) {
				//Construct full haplotypes from cross-time data
				MakeInitialHaps (p,c_m_polys,haps);
				ConstructFullHaps2(p,haps);
			} else {
				//Import haplotype data
				InputHaplotypeData(p,polys,haps);
				//Edit polymorphism data to span all haplotypes even if not explicitly observed
				CompleteMLCalls(haps,c_m_polys);
			}
			vector< vector<mtr> > mltrajs;
			GenerateOutputFiles(p,polys,haps,m_polys,c_m_polys,mltrajs);
		}
		return 0;
		
	} else if (method.compare("ml_noise")==0) {
		
		//Read multi-locus haplotype data.
		vector< vector<mtr> > mltrajs;
		vector< vector<int> > times;
		if (p.full_haps==0) {
			ReadMultiLocTrajs (p,times,mltrajs);
		} else {
			ReadMultiMLTrajs(p,times,mltrajs);
		}
		//Should include an alternative arrangement in which data is read from a single set of haplotypes i.e. within a single dataset...
		
		//cout << "Size " << times.size() << "\n";

		
		FilterMLFreq (p,times,mltrajs); //Edit out non-polymorphic time-points
		
		cout << "Size " << times.size() << "\n";
		
		FilterMLFreq2 (p,times,mltrajs); //Edit out haplotypes below the noise threshold
		//Calculate mean multi-locus observations over time
		MLPMeanFreqs(p,mltrajs);
		
		vector<int> include_ml;
		for (int i=0;i<mltrajs.size();i++) {
			include_ml.push_back(1);
		}
		
		//Filter observations by change in haplotype frequency - removes partial haplotype sets that change too much
		FilterMLFreq3(p,include_ml,times,mltrajs);
		
		vector<double> fact_store;
		GetFactVectorML (mltrajs,fact_store);
		//Calculate estimate for partial haplotype noise level
		
		double Cml_opt;
		OptimiseCml(p,Cml_opt,include_ml,mltrajs,fact_store,rgen);
		ofstream cml_file;
		cml_file.open("Cml.out");
		cml_file << Cml_opt << "\n";
		cml_file.close();

		if (p.conservative==1) {
			ConservativeCml (p,Cml_opt,include_ml,mltrajs,fact_store,rgen);
			ofstream cml_file;
			cml_file.open("Cml_cons.out");
			cml_file << Cml_opt << "\n";
			cml_file.close();
		}
		
		return 0;

	} else if (method.compare("ld_calc")==0) {
		if (p.get_in==0) {
			p.in_file="Single_locus_trajectories_sl.out";
		}
		
		double Cml=0;
		ifstream cml_file;
		cml_file.open("Cml.out");
		cml_file >> Cml;
		
		vector<str> sltrajs;
		ImportSLTData(p,sltrajs);
		if (p.verb==1) {
			for (int i=0;i<sltrajs.size();i++) {
				cout << sltrajs[i].locus << " " << sltrajs[i].cons << " " << sltrajs[i].nuc << "\n";
			}
		}
		
		//Read in multi-locus trajectories
		vector< vector<int> > times;
		vector< vector<mtr> > mltrajs;
		
		//Vector of information
		vector<ld_info> ld_data;
		CollectLDInfo(sltrajs,ld_data);
		
		ReadMultiLocTrajs2(p,sltrajs,ld_data);
		vector<double> fact_store;
		ProcessLDInfo (fact_store,ld_data);
		
		//Optimise Dij
		cout << "Calculating LD statistics...\n";
		
		ofstream ld_file;
		ld_file.open("LD_info.out");
		for (int s=0;s<ld_data.size();s++) {
			for (int t=0;t<8;t++) {
				OptimiseLD(p,Cml,ld_data[s],t,fact_store,ld_file,rgen);
			}
		}
	}  else {
		PrintInstructions (0);
		return 0;
	}
}
							
 




