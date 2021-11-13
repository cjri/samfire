//Program to process .sam file data

#include <iostream>
#include <vector>
#include <string>
#include <sstream>


using namespace std;

#include "shared.h"
#include "alignment.h"
#include "bioinf.h"
#include "call_snps.h"
#include "call_mnps.h"
#include "contamination.h"
#include "ddups.h"
#include "data.h"
#include "findld.h"
#include "ld_simple.h"
#include "io.h"
#include "stats.h"
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
            if (p.print_ali==1) {
                OutputAlign(i,a.data);
            }
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
			JoinPairs (i,p,jseqs,qual,a.data);

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
		
		OutputSLTData(p,1,p.out_file.c_str(),sltrajs);
		
		return 0;
		
	} else if (method.compare("stats")==0) {
		
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
		//Import variant data
		vector< vector< vector<int> > > all_nucs;
		vector< vector<char> > all_cons;
		InputVariantData(0,sam_files,all_nucs,all_cons);

		//Calculate distribution of read lengths
		if (p.calc_len==1) {
			cout << "Calculating read lengths\n";
			vector<nuc> ref_counts;
			for (int i=0;i<sam_files.size();i++) {
				vector<joined> t_read;
				InputJnData (i,t_read);
				vector<int> l_dist;
				Calc_Seq_Lengths(t_read,l_dist);
				PrintSeqLens(i,l_dist);
			}
		}
		
		//Calculate distribution of PHRED scores
		if (p.calc_qual==1) {
			cout << "Calculating PHRED score distribution\n";
			vector<ql> q_read;
			vector<char> qual;
			makequal(p,qual);
			vector<long> q_dist;
			for (int i=0;i<sam_files.size();i++) {
				InputQualData(p,i,q_read,qual);
				Calc_Base_Quality(q_read,q_dist);
				PrintSeqQual(i,q_dist);
			}
		}
		
		//Generate bootstrap samples
		vector< vector< vector< vector<int> > > > bootstrap;
		if (p.bootstrap==1) {
			MakeBootstraps (all_nucs,bootstrap,rgen);
		}
		
		//Calculate sequence diversity \pi
		if (p.calc_pi==1) {
			Calc_Pi_Diversity(p,all_nucs,bootstrap);
		}
		
		//Calculate variant composition
		if (p.calc_var_comp==1) {
			cout << "Calculate composition of low-frequency variants\n";
			Calc_Variant_Composition(p,all_nucs,all_cons,bootstrap);
		}
		
		//Calculate Hamming distances between consensus sequences
		if (p.calc_ham_cons==1) {
			Calculate_Consensus_Hamming_Distances(sam_files);
		}
		
		//Calculate Hamming distances between Joined files
		if (p.calc_ham_var==1) {
			Calculate_Variant_Hamming_Distances(all_nucs);
		}

		//Calculate sum of polymorphisms
		if (p.calc_sum_var==1) {
			Calc_Tot_Freqs (p,all_nucs,bootstrap);
		}
		
		return 0;

		
	} else if (method.compare("reading_frames")==0) {
		if (p.get_in==0) {
			p.in_file="Consensus0.fa";
		}

		Find_Reading_Frames(p);
		return 0;
		
		//Next step here - use proteins to find synonymous and non-synonymous mutations across the reference
		
	} else if (method.compare("variant_types")==0) {
		Find_Variant_Types(p);
		return 0;
		
		
	} else if (method.compare("codon_traj")==0) {

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
		
		//Setup Amino acid translation
		vector< vector<char> > cod;
		vector<char> aa;
		SetupAATrans (cod,aa);
		
		cout << "Look for codons.  Start point for translation is " << p.trans_start << "\n";

		cout << refseq.size << "\n";
		int n_cods=(refseq.size+1)/3;
		int mod=p.trans_start % 3;
		cout << "Mod " << mod << "\n";
		
		//Setup structure to contain codons by position and time point.  Vector of 64 codons
		vector< vector< vector<int> > > cc;
		vector<int> cod_count (64,0);
		for (int t=0;t<sam_files.size();t++) {
			//Each time point
			vector< vector<int> > t_cc;
			for (int i=0;i<=n_cods;i++) {
				//Position
				t_cc.push_back(cod_count);
			}
			cc.push_back(t_cc);
		}

		int min_pos=10000;
		int max_pos=0;
		
//		for (int i=0;i<sam_files.size();i++) {
		for (int t=0;t<sam_files.size();t++) {
			vector<joined> t_read;
			InputJnData (t,t_read);
			
			
			for (int i=0;i<t_read.size();i++) {
		//	for (int i=4;i<5;i++) {
			//	cout << i << " Initial " << t_read[i].alpos << " " << t_read[i].seq.size() << "\n";
				//for (int j=0;j<t_read[i].seq.size();j++) {
				//	cout << t_read[i].seq[j];
				//}
				//cout << "\n";
				int rstart=p.trans_start-t_read[i].alpos;
			//	cout << "Rstart " << rstart << "\n";

				if (rstart<0) {
			//		cout << "Fix\n";
					rstart=(3000+p.trans_start-t_read[i].alpos)%3;
				}
			//	cout << "Read from " <<  p.trans_start << "\n";
			//	cout << "Rstart " << rstart << "\n";
				for (int j=rstart;j<t_read[i].seq.size()-2;j=j+3) {
					int p1=j;
					int p2=j+1;
					int p3=j+2;
					int pos=(t_read[i].alpos+j)/3;
			//		cout << t_read[i].alpos+j << " " << t_read[i].seq[p1] << t_read[i].seq[p2] << t_read[i].seq[p3] << "\n";
					int n=FastTranslateCodonN(t_read[i].seq[p1],t_read[i].seq[p2],t_read[i].seq[p3]);
					//if (pos==461) {
						//cout << i << " " << t_read[i].alpos << " " << j << " " << t_read[i].alpos+j << " " << rstart << " " << t_read[i].seq[p1] << t_read[i].seq[p2] << t_read[i].seq[p3] << " " << n << "\n";
					//}
				//	cout << pos << " " << n << "\n";
					if (n>0) {
						if (pos<min_pos) {
							min_pos=pos;
						}
						if (pos>max_pos) {
							max_pos=pos;
						}
			//			cout << "Push " << t << " " << pos << " " << n << " " << cc.size() << " " << cc[t].size() << " " << cc[t][pos].size() << "\n";
						cc[t][pos][n]++;
					}
				}
			}
			/*cout << "Size " << cc[0].size() << "\n";
			for (int i=0;i<64;i++) {
				char c=TranslateNAA(i,aa);
				cout << c << " ";
			}
			cout << "\n";
			for (int i=0;i<cc[0].size();i++) {
				cout << i << " ";
				for (int j=0;j<cc[0][i].size();j++) {
					cout << cc[0][i][j] << " ";
				}
				cout << "\n";
			}*/
			
		}
		
		cout << min_pos << " " << max_pos << "\n";
		
		//Re-process cc data
		vector< vector< vector<int> > > cc_new;
		for (int i=0;i<cc.size();i++) {
			vector< vector<int> > cr;
			for (int j=min_pos;j<=max_pos;j++) {
				cr.push_back(cc[i][j]);
			}
			cc_new.push_back(cr);
		}
		for (int i=0;i<cc_new[0].size();i++) {
			cout << i << " ";
			for (int j=0;j<cc_new[0][i].size();j++) {
				cout << cc_new[0][i][j] << " ";
			}
			cout << "\n";
		}
		cc=cc_new;
		
		//Next step - look for polymorphisms in the data
		
		//Setup
		vector<npoly> np;
		
		for (int i=0;i<cc.size();i++) {
			for (int j=0;j<cc[i].size();j++) {
				int max=-1;
				int mpos=0;
				//Find max;
				int tot=0;
				for (int k=0;k<cc[i][j].size();k++) {
					tot=tot+cc[i][j][k];
					if (cc[i][j][k]>max) {
						max=cc[i][j][k];
						mpos=k;
					}
				}
				//Find polymorphism
				for (int k=0;k<cc[i][j].size();k++) {
					if (k!=mpos) {
						double q=(cc[i][j][k]+0.)/(max+0.);
						if (q>p.q_cut&&cc[i][j][k]>p.n_min) {
							npoly pf;
							pf.locus=j;
							pf.cod=k;
							np.push_back(pf);
							cout << "Found poly time " << i << " locus " << j << " cod " << k << "\n";
						}
					}
				}
			}
		}
		vector<int> all_loci;
		vector< vector<int> > all_cod;
		vector<int> v;
		for (int i=0;i<cc[0].size();i++) {
			all_loci.push_back(0);
			all_cod.push_back(v);
		}
		for (int i=0;i<np.size();i++) {
			all_loci[np[i].locus]=1;
			int match=0;
			for (int j=0;j<all_cod[np[i].locus].size();j++) {
				if (all_cod[np[i].locus][j]==np[i].cod) {
					match=1;
				}
			}
			if (match==0) {
				all_cod[np[i].locus].push_back(np[i].cod);
			}
		}
		
		cout << "Print polymorphic sites\n";
		for (int i=0;i<all_loci.size();i++) {
			if (all_loci[i]>0) {
				cout << i << " ";
				sort(all_cod[i].begin(),all_cod[i].end());
				for (int j=0;j<all_cod[i].size();j++) {
					cout << all_cod[i][j] << " ";
				}
				cout << "\n";
			}
		}
		
		//Calculate polymorphism data
		
		cout << "np size " << np.size() << "\n";
		
		vector<npoly> new_np;

		
		//Find initial consensus and rationalise polymorphism data
		for (int i=0;i<all_loci.size();i++) {
			if (all_loci[i]>0) {
				int cons=-1;
				int c_max=0;
				for (int j=0;j<cc[0][i].size();j++) {
					if (cc[0][i][j]>c_max) {
						c_max=cc[0][i][j];
						cons=j;
					}
				}
				npoly newp;
				newp.cons=cons;
				newp.locus=i;
				for (int j=0;j<all_cod[i].size();j++) {
					newp.cod=all_cod[i][j];
					if (newp.cod!=newp.cons) {
						new_np.push_back(newp);
					}
				}
			}
		}
		
		cout << "new_np size " << new_np.size() << "\n";
		np=new_np;
		
		//Collect polymorphism data
		for (int i=0;i<np.size();i++) {
			for (int t=0;t<cc.size();t++) {
				np[i].count.push_back(cc[t][np[i].locus][np[i].cod]);
				np[i].ccount.push_back(cc[t][np[i].locus][np[i].cons]);
				int tot=0;
				for (int j=0;j<64;j++) {
					tot=tot+cc[t][np[i].locus][j];
				}
				np[i].tot.push_back(tot);
			}
		}
		
		//Print polymorphism data
		for (int i=0;i<np.size();i++) {
			cout << np[i].locus << " " << np[i].cons << " " << np[i].cod << "\n";
			for (int j=0;j<np[i].ccount.size();j++) {
				cout << np[i].ccount[j] << " ";
			}
			cout << "\n";

			for (int j=0;j<np[i].count.size();j++) {
				cout << np[i].count[j] << " ";
			}
			cout << "\n";
		}
		
		
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
		vector< vector<char> > all_cons;
		InputVariantData(1,sam_files,all_nucs,all_cons);
		
		cout << all_nucs[0].size() << "\n";
		
		//Output consensus sequence data
		vector<char> global_cons;
		OutputGlobalConsensus(all_nucs,global_cons);
		
		//Branch here - do translation?
		
		if (p.translate==1) {
			cout << "Translate consensus sequences.  Start point is " << p.trans_start << "\n";
			cout << "Use --trans_start to change start point for translation\n";
			//Translate consensus sequences...
		
			//Output coding sequence
			OutputCodingConsensii(p,sam_files,all_cons);
			
			//Generate lookup table of nucleotides
			vector< vector<char> > cod;
			vector<char> aa;
			SetupAATrans(cod,aa);
		
			//Do translation of nucleotide sequence from specified start codon.
			for (int i=0;i<all_cons.size();i++) {
				//Repair consensus sequence if p.repair_consensus>0.
				//Aim here is to remove 'N' nucleotides so as to enable translation
				RepairConsensusSequence(p,i,all_cons,global_cons);
				vector<char> aaseq;
				ostringstream convert;
				convert << i;
				string temp=convert.str();
				string name2 = "Consensus"+temp+"_AA.fa";
				string name3 = ">Consensus_"+temp;
				int j=p.trans_start-1;
				int detect_error=0;
				vector< vector<int> > s_types;
				while (j<all_cons[i].size()) {
					vector<char> c;
					for (int k=j;k<j+3;k++) {
						c.push_back(all_cons[i][k]);
						if (detect_error==0) {
							if (all_cons[i][k]=='N') {
								cout << "Warning: Ambiguous nucleotide in sequence " << i << "\n";
								cout << "Translation of nucleotide sequence may contain errors\n";
								cout << "You could look into the --repair_consensus options.\n";
								detect_error=1;
							}
						}
					}
					char s=FastTranslateCodon(c,aa);
				
					//Branch here - Get variant types?
					if (p.get_variants==1) {
						GetVarTypes(c,s,cod,aa,s_types);
					}
				
					j=j+3;
					if (s=='X') { //End output upon getting to a stop codon
						break;
					}
					aaseq.push_back(s);
				}
				
				OutputParticularAAConsensus(name2,name3,aaseq);
				//Branch here - Get variant types?
				if (p.get_variants==1) {
					string name4 = "Variant_mask"+temp+".out";
					OutputVarMask(name4,p.trans_start,all_cons[i].size(),s_types);
				}
			}
			//Thoughts here - the translation can go wrong in a case where a start codon is missing at a particular time point; a lack of translation can then occur.  A solution would be to fix N values in the nucleotide consensus with the Consensus_all nucleotides.  Might want a --fix flag to do this, with the defalt set to 1?
			
		}
		
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
		
	} else if (method.compare("ml_load")==0) {
		
		//Get reference sequence
		ifstream ref_file;
		ref_file.open(p.ref.c_str());
		rseq refseq;
		GetRefSeq (ref_file,refseq);
		
		vector<string> sam_files;
		if (p.no_sam>0) {
			for (int i=0;i<p.no_sam;i++) {
				string s;
				sam_files.push_back(s);
			}
		} else {
			ImportSamFileNames(p,sam_files);
		}
		
		
		PerSequenceSNPs (p,refseq,sam_files);
		
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
		OutputSLTData(p,1,p.out_file.c_str(),sltrajs);

		if (p.printnl==1) {
			p.out_file="Single_locus_trajectories_nl.out";
			OutputSLTData(p,0,p.out_file.c_str(),sltrajs);
		}
		
		return 0;
		
	} else if (method.compare("contamination")==0) {
		
		//Read in consensus sequences
		ifstream ref_file1;
		ifstream ref_file2;
		ref_file1.open(p.ref1.c_str());
		ref_file2.open(p.ref2.c_str());
		if (p.ref1==""||p.ref2=="") {
			cout << "Error: Must specify two reference files.  Use --ref1 and --ref2 \n";
			return 0;
		}
		rseq refseq1;
		rseq refseq2;
		GetRefSeq (ref_file1,refseq1);
		GetRefSeq (ref_file2,refseq2);
		string ref1=refseq1.seq;
		string ref2=refseq2.seq;

		//Find reference length
		int len=0;
		FindRefLen(len,ref1,ref2);

		
		//Read in polymorphic loci
		vector<int> polys;
		int start_c=-1;
		int check;
		FindPolymorphisms(p,len,ref1,ref2,check,start_c,polys);
		if (check==1) {
			return 0;
		}
		
		//Read in Joined 1
		vector<joined> t_read1;
		vector<joined> t_read2;
		if (p.jn1==""||p.jn2=="") {
			cout << "Error: Must specify two Joined.out files.  Use --jn1 and --jn2 \n";
			return 0;
		}
		InputJnDataGeneral(p.jn1,t_read1);
		InputJnDataGeneral(p.jn2,t_read2);
		
		//Go through reads individually and find number of mismatches to each reference.
		string c1="Flagged_reads1.out";
		string c2="Flagged_reads2.out";
		vector<int> filter1;
		vector<int> filter2;
		FindDifferences(p,len,start_c,ref1,ref2,c1,polys,t_read1,filter1);
		FindDifferences(p,len,start_c,ref2,ref1,c2,polys,t_read2,filter2);
		
		cout << "Data from potentially contaminating reads in " << p.jn1 << " shown in Flagged_reads1.out\n";
		cout << "Data from potentially contaminating reads in " << p.jn2 << " shown in Flagged_reads2.out\n";
		
		if (p.decon==1) {
			OutputDeconJoined(p,t_read1,t_read2,filter1,filter2);
		}
		
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
		for (int i=0;i<n_times;i++) {
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
			if (p.multi_gap==0) {
				CallReadTypesConsecutive(p,max_l,max_r_size,polys,l_combs);  //Note for longer genomes - could consider limit on #loci.  Are full haplotypes required?
			} else {
				CallReadTypesGaps(p,max_l,max_r_size,polys,l_combs);  //Note for longer genomes - could consider limit on #loci.  Are full haplotypes required?
			}
			
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
		
		int max_loc=0;
		vector< vector<int> > l_vec;
		if (p.multi_gap==0) {
			cout << "Find Max Loc\n";
			FindMaxLoc(max_loc,all_l_combs);
			ConstructLocVec (p,max_loc,max_l_store,all_l_combs,l_vec);
		} else {
			ConstructLocVecGap (p,all_l_combs,l_vec);
		}


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
		
    } else if (method.compare("ld_simple")==0) {
        if (p.get_in==0) {
            p.in_file="Single_locus_trajectories_sl.out";
        }
       
        //Read in list of sites
        vector<str> sltrajs;
        ImportSLTData(p,sltrajs);
        if (p.verb==1) {
            for (int i=0;i<sltrajs.size();i++) {
                cout << sltrajs[i].locus << " " << sltrajs[i].cons << " " << sltrajs[i].nuc << "\n";
            }
        }
        //Get polymorphic sites
        vector<int> polys;
        for (int i=0;i<sltrajs.size();i++) {
            polys.push_back(sltrajs[i].locus);
        }
        //Read in times
        vector<int> times;
        ImportTimeData(times);
        int n_times=times.size();

        //Load sam file data
        vector<string> sam_files;
        ImportSamFileNamesFix(sam_files);
        cout << sam_files.size() << "\n";
        
        for (int t=0;t<sam_files.size();t++) { //Read in Joined file
            vector<joined> t_read;
            InputJnData (t,t_read);
            cout << "Sorting data\n";
            sort(t_read.begin(),t_read.end(),compareJoined);
            
            vector <vector<char> > m_vars;
            MapSeqToSNPLoci(p,n_times,polys,t_read,m_vars);
            
            //Count how many polymorphic sites are in each read.  Vector in n_sites
            vector< vector<int> > n_sites;
            ConstructNSites(m_vars,n_sites);
            
            int maxs=CountSNPsPerRead(n_sites);
            cout << "Max coverage of SNPs in a read = " << maxs << "\n";
            if (maxs==1) {
                cout << "Insufficient to calculate LD\n";
                return 0;
            }
            
            //Make grid of pairs
            vector< vector<int> > pairgrid;
            MakePairgrid(polys,n_sites,pairgrid);
            
            //Vector to store LD information - go through sites and find pairs
            if (p.verb==1) {
                cout << "Construct LD data\n";
            }
            vector<ld_info> ld_data;
            ConstructLDPairsData (pairgrid,ld_data);
            
            //Go through again.  Import data into ld_data.
            if (p.verb==1) {
                cout << "Collate LD data " << ld_data.size() << "\n";
            }

            CollateLDData (t,sltrajs,m_vars,n_sites,ld_data);
            ConvertLocusNumbers (sltrajs,ld_data);
            
            if (p.verb==1) {
                OutputLDRawInformation(t,ld_data);
            }
            
            //Optimise LD : Calculations
            vector<double> fact_store;
            OptimiseLDSetup (fact_store,ld_data);
            
            //Note : Check for duplicate sites
            
            if (p.noopt==0) { //Flag here - might not want to do the optimisation?
                RunLDOptimisation (p,t,fact_store,ld_data,rgen);
            }
        }
        
        
        
        
        
        
        
        
        
        
        
	} else if (method.compare("calc_distances")==0) {
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
		vector< vector<char> > all_cons;
		InputVariantData(0,sam_files,all_nucs,all_cons); //All_nucs contains A, C, G, T by position
		
		//Get variant frequencies
		vector< vector< vector<double> > > all_freqs;
		CalculateVariantFrequencies(all_nucs,all_freqs);
		
		//Calculate distances between frequency vectors.  Metric is half of the absolute distance between vectors
		cout << "Calculating distances between Variant files\n";
		vector< vector<double> > distances;
		vector< vector<double> > all_counts;
		CalculateVariantDistances(all_freqs, distances, all_counts);
		OutputDistancesCounts(distances,all_counts);
		
		if (p.sns_distances==1) {
			cout << "Calculating S and NS distances requires Variant_mask files - generated from consensus routine\n";
			//Import Variant_mask files;
			cout << "Reading in variant_mask files...\n";
			vector< vector< vector<int> > > all_types;
			ImportVariantMasks(sam_files,all_types);

			//Calculate distances for S and NS
			vector< vector<double> > distances_s;
			vector< vector<double> > distances_n;
			vector< vector<double> > all_counts_s;
			vector< vector<double> > all_counts_n;
			cout << "Calculating S and NS distances between Variant files\n";
			CalculateVariantDistancesSNS (all_freqs,all_types,distances_s,distances_n,all_counts_s,all_counts_n);
			OutputDistancesCountsSNS (distances_n,distances_s,all_counts_n,all_counts_s);
		}
		
		
		
	}  else {
		PrintInstructions (0);
		return 0;
	}
}
							
 




