//============================================================================
// Name        : bed2wig.cpp
// Author      : Wouter Van Gool (wouter@igbmc.fr)
// Version	   : 0.1
// Copyright   : GPL
// Description : Bed2Wig-writer C++, Ansi-style
//============================================================================
#include "rcpp_wigWriter.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include <map>

#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/bam_io.h>

using namespace std;
using namespace boost;
using namespace Rcpp;

long int  getChromosomeLength(string genomeBuild, string chromosomeName) {

	long int chromosomeLength = -1;

	if (genomeBuild == "mm8") {

		std::map<string, int> mm8;
        mm8["chr1"]  = 197069962;
        mm8["chr2"]  = 181976762;
        mm8["chr3"]  = 159872112;
        mm8["chr4"]  = 155029701;
        mm8["chr5"]  = 152003063;
        mm8["chr6"]  = 149525685;
        mm8["chr7"]  = 145134094;
        mm8["chr8"]  = 132085098;
        mm8["chr9"]  = 124000669;
        mm8["chr10"] = 129959148;
        mm8["chr11"] = 121798632;
        mm8["chr12"] = 120463159;
        mm8["chr13"] = 120614378;
        mm8["chr14"] = 123978870;
        mm8["chr15"] = 103492577;
        mm8["chr16"] = 98252459;
        mm8["chr17"] = 95177420;
        mm8["chr18"] = 90736837;
        mm8["chr19"] = 61321190;
        mm8["chrX"]  = 165556469;
        mm8["chrY"] = 165556469;

		if(mm8.count(chromosomeName)>0){

			chromosomeLength = mm8[chromosomeName];
		}
	}

	if (genomeBuild == "mm9") {

		std::map<string, int> mm9;
		mm9["chr1"] = 197195432;
		mm9["chr2"] = 181748087;
		mm9["chrX"] = 166650296;
		mm9["chr3"] = 159599783;
		mm9["chr4"] = 155630120;
		mm9["chr5"] = 152537259;
		mm9["chr7"] = 152524553;
		mm9["chr6"] = 149517037;
		mm9["chr8"] = 131738871;
		mm9["chr10"] = 129993255;
		mm9["chr14"] = 125194864;
		mm9["chr9"] = 124076172;
		mm9["chr11"] = 121843856;
		mm9["chr12"] = 121257530;
		mm9["chr13"] = 120284312;
		mm9["chr15"] = 103494974;
		mm9["chr16"] = 98319150;
		mm9["chr17"] = 95272651;
		mm9["chr18"] = 90772031;
		mm9["chr19"] = 61342430;
		mm9["chrY"] = 15902555;
		mm9["chrM"] = 16299;
		if(mm9.count(chromosomeName)>0){

			chromosomeLength = mm9[chromosomeName];
		}
	}

	if(genomeBuild == "hg18"){

		std::map<string, int> hg18;
        hg18["chr1"] = 247249719;
        hg18["chr2"] = 242951149;
        hg18["chr3"] = 199501827;
        hg18["chr4"] = 191273063;
        hg18["chr5"] = 180857866;
        hg18["chr6"] = 170899992;
        hg18["chr7"] = 158821424;
        hg18["chr8"] = 146274826;
        hg18["chr9"] = 140273252;
        hg18["chr10"] = 135374737;
        hg18["chr11"] = 134452384;
        hg18["chr12"] = 132349534;
        hg18["chr13"] = 114142980;
        hg18["chr14"] = 106368585;
        hg18["chr15"] = 100338915;
        hg18["chr16"] = 88827254;
        hg18["chr17"] = 78774742;
        hg18["chr18"] = 76117153;
        hg18["chr19"] = 63811651;
        hg18["chr20"] = 62435964;
        hg18["chr21"] = 46944323;
        hg18["chr22"] = 49691432;
        hg18["chrX"]  = 154913754;
        hg18["chrY"]  = 57772954;

		if(hg18.count(chromosomeName)>0){

			chromosomeLength = hg18[chromosomeName];
		}
	}

	if(genomeBuild == "hg19"){

		std::map<string, int> hg19;
		hg19["chr1"]  = 249250621;
		hg19["chr2"]  = 243199373;
		hg19["chr3"]  = 198022430;
		hg19["chr4"]  = 191154276;
		hg19["chr5"]  = 180915260;
		hg19["chr6"]  = 171115067;
		hg19["chr7"]  = 159138663;
		hg19["chr8"]  = 146364022;
		hg19["chr9"]  = 141213431;
		hg19["chr10"] = 135534747;
		hg19["chr11"] = 135006516;
		hg19["chr12"] = 133851895;
		hg19["chr13"] = 115169878;
		hg19["chr14"] = 107349540;
		hg19["chr15"] = 102531392;
		hg19["chr16"] = 90354753;
		hg19["chr17"] = 81195210;
		hg19["chr18"] = 78077248;
		hg19["chr19"] = 59128983;
		hg19["chr20"] = 63025520;
		hg19["chr21"] = 48129895;
		hg19["chr22"] = 51304566;
		hg19["chrX"]  = 155270560;
		hg19["chrY"]  = 59373566;

		if(hg19.count(chromosomeName)>0){

			chromosomeLength = hg19[chromosomeName];
		}
	}

	if(genomeBuild == "dm2"){

		std::map<string, int> dm2;
        dm2["chr3R"] = 27905053;
        dm2["chr3L"] = 23771897;
        dm2["chr2L"] = 22407834;
        dm2["chrX"] = 22224390;
        dm2["chr2R"] = 20766785;
        dm2["chrU"] = 8724946;
        dm2["chr3h"] = 2955737;
        dm2["chr2h"] = 1694122;
        dm2["chr4"] = 1281640;
        dm2["chrYh"] = 396896;
        dm2["chrXh"] = 359526;
        dm2["chr4h"] = 88110;
        dm2["chrM"] = 19517;

		if(dm2.count(chromosomeName)>0){

			chromosomeLength = dm2[chromosomeName];
		}
	}


	if(genomeBuild == "dm3"){

		std::map<string, int> dm3;
        dm3["chrUextra"] = 29004656;
        dm3["chr3R"] = 27905053;
        dm3["chr3L"] = 24543557;
        dm3["chr2L"] = 23011544;
        dm3["chrX"] = 22422827;
        dm3["chr2R"] = 21146708;
        dm3["chrU"] = 10049037;
        dm3["chr2RHet"] = 3288761;
        dm3["chr3LHet"] = 2555491;
        dm3["chr3RHet"] = 2517507;
        dm3["chr4"] = 1351857;
        dm3["chr2LHet"] = 368872;
        dm3["chrYHet"] = 347038;
        dm3["chrXHet"] = 204112;
        dm3["chrM"] = 19517;

		if(dm3.count(chromosomeName)>0){

			chromosomeLength = dm3[chromosomeName];
		}
	}

	if(genomeBuild == "ce4"){

		std::map<string, int> ce4;
        ce4["chrV"] = 20919398;
        ce4["chrX"] = 17718852;
        ce4["chrIV"] = 17493784;
        ce4["chrII"] = 15279316;
        ce4["chrI"] = 15072419;
        ce4["chrIII"] = 13783681;
        ce4["chrM"] = 13794;

		if(ce4.count(chromosomeName)>0){

			chromosomeLength = ce4[chromosomeName];
		}
	}

	if(genomeBuild == "ce6"){

		std::map<string, int> ce6;
        ce6["chrV"] = 20919568;
        ce6["chrX"] = 17718854;
        ce6["chrIV"] = 17493785;
        ce6["chrII"] = 15279323;
        ce6["chrI"] = 15072421;
        ce6["chrIII"] = 13783681;
        ce6["chrM"] = 13794;

		if(ce6.count(chromosomeName)>0){

			chromosomeLength = ce6[chromosomeName];
		}
	}

	if(genomeBuild == "danRer4"){

		std::map<string, int> danRer4;
        danRer4["chr14"] = 91717235;
        danRer4["chr7"] = 87691871;
        danRer4["chr5"] = 84656180;
        danRer4["chr3"] = 77179095;
        danRer4["chr1"] = 70589895;
        danRer4["chr6"] = 69554819;
        danRer4["chr8"] = 66798501;
        danRer4["chr16"] = 65489547;
        danRer4["chr13"] = 64258675;
        danRer4["chr20"] = 63653707;
        danRer4["chr17"] = 63411520;
        danRer4["chr2"] = 61889685;
        danRer4["chr18"] = 59765243;
        danRer4["chr12"] = 58719258;
        danRer4["chr15"] = 57214918;
        danRer4["chr21"] = 56255777;
        danRer4["chr9"] = 55712184;
        danRer4["chr10"] = 54070595;
        danRer4["chr23"] = 53215897;
        danRer4["chr11"] = 52342180;
        danRer4["chr19"] = 51715404;
        danRer4["chr22"] = 47751166;
        danRer4["chr4"] = 47249802;
        danRer4["chr24"] = 46081529;
        danRer4["chr25"] = 40315040;
        danRer4["chrM"] =  16596;

		if(danRer4.count(chromosomeName)>0){

			chromosomeLength = danRer4[chromosomeName];
		}
	}


	if(genomeBuild == "danRer6"){

		std::map<string, int> danRer6;
        danRer6["chr7"] = 76918211;
        danRer6["chr5"] = 74451498;
        danRer6["chr4"] = 71658100;
        danRer6["chr6"] = 61647013;
        danRer6["chr3"] = 60907308;
        danRer6["chr1"] = 59305620;
        danRer6["chr2"] = 58009534;
        danRer6["chr8"] = 55568185;
        danRer6["chr9"] = 54736511;
        danRer6["chr14"] = 52930158;
        danRer6["chr16"] = 51890894;
        danRer6["chr20"] = 51884995;
        danRer6["chr13"] = 50748729;
        danRer6["chr17"] = 49469313;
        danRer6["chr18"] = 49271716;
        danRer6["chr19"] = 48708673;
        danRer6["chr21"] = 47572505;
        danRer6["chr15"] = 47237297;
        danRer6["chr12"] = 46853116;
        danRer6["chr23"] = 44714728;
        danRer6["chr11"] = 44116856;
        danRer6["chr10"] = 43467561;
        danRer6["chr22"] = 41415389;
        danRer6["chr24"] = 40403431;
        danRer6["chr25"] = 38768535; 

		if(danRer6.count(chromosomeName)>0){

			chromosomeLength = danRer6[chromosomeName];
		}
	}


	if(genomeBuild == "rn4"){

		std::map<string, int> rn4;
        rn4["chr1"] = 267910886;
        rn4["chr2"] = 258207540;
        rn4["chr4"] = 187126005;
        rn4["chr5"] = 173096209;
        rn4["chr3"] = 171063335;
        rn4["chrX"] = 160699376;
        rn4["chr6"] = 147636619;
        rn4["chr7"] = 143002779;
        rn4["chr8"] = 129041809;
        rn4["chr9"] = 113440463;
        rn4["chr14"] = 112194335;
        rn4["chr13"] = 111154910;
        rn4["chr10"] = 110718848;
        rn4["chr15"] = 109758846;
        rn4["chr17"] = 97296363;
        rn4["chr16"] = 90238779;
        rn4["chr11"] = 87759784;
        rn4["chr18"] = 87265094;
        rn4["chrUn"] = 75822765;
        rn4["chr19"] = 59218465;
        rn4["chr20"] = 55268282;
        rn4["chr12"] = 46782294;

		if(rn4.count(chromosomeName)>0){

			chromosomeLength = rn4[chromosomeName];
		}
	}


	if(genomeBuild == "rn3"){

		std::map<string, int> rn3;
		rn3["chr1"] = 268121971;
		rn3["chr2"] = 258222147;
		rn3["chr4"] = 187371129;
		rn3["chr5"] = 173106704;
		rn3["chr3"] = 170969371;
		rn3["chrX"] = 160775580;
		rn3["chr6"] = 147642806;
		rn3["chr7"] = 143082968;
		rn3["chr8"] = 129061546;
		rn3["chr9"] = 113649943;
		rn3["chr14"] = 112220682;
		rn3["chr13"] = 111348958;
		rn3["chr10"] = 110733352;
		rn3["chr15"] = 109774626;
		rn3["chr17"] = 97307196;
		rn3["chr16"] = 90224819;
		rn3["chr11"] = 87800381;
		rn3["chr18"] = 87338544;
		rn3["chrUn"] = 75822765;
		rn3["chr19"] = 59223525;
		rn3["chr20"] = 55296979;
		rn3["chr12"] = 46649226;

		if(rn3.count(chromosomeName)>0){

			chromosomeLength = rn3[chromosomeName];
		}
	}
	return chromosomeLength;
}


SEXP countReads(SEXP bedIn){

	string bed = as <string>(bedIn);
	ifstream in(bed.c_str());
	string line;
	int counter = 0;
	while (getline(in, line)) {

		if( line.substr(0, 1)!="#"){

			counter++;
		}
	}

	Rcpp::NumericVector r(1);
	r[0] = counter;
	return r;
}


SEXP removeClonalReads(SEXP fileIn, SEXP readsToKeep, SEXP outputDir, SEXP format){

	string fin = as<string>(fileIn);
	//cout << "Removing clonal reads from file " << fin << "..." << endl;
	
	string dirOut = as<string>(outputDir);
	int nrReadsToKeep = as<int>(readsToKeep);
	string nrDuplicates = boost::lexical_cast<string>(nrReadsToKeep);
	//string nrDuplicates = "test";
	string ft = as<string>(format);

	std::map<string,int> readMap;
	double totalReads = 0;

	if(ft == "bed"){

		ifstream in(fin.c_str());
		string line;

		while (getline(in, line)) {

			if( line.substr(0, 1)!="#"){

				vector<string> fields;
				typedef boost::tokenizer<boost::char_separator<char> >
						tokenizer;
				boost::char_separator<char> sep("\t");
				tokenizer tokens(line, sep);
				fields.assign(tokens.begin(), tokens.end());

				if(fields.size()==3){

					if(readMap.count(line) == 0) {
						readMap[line] = 1;
					} else {
						readMap[line]++;
					}	
					totalReads++;
				}

				if(fields.size()==6){

					string chr = fields[0];
					string start = fields[1];
					string end = fields[2];
					string strand = fields[5];

					string hashKey = chr + "\t" + start + "\t" + end + "\tfilteredClonalReads\t0\t" + strand;

					if(readMap.count(hashKey) == 0) {
						readMap[hashKey] = 1;
					} else {
						readMap[hashKey]++;
					}	
					totalReads++;
				}
			}
		}
		in.close();

	}else if((ft == "bam") || (ft == "sam")){

		seqan::BamStream bamStreamIn(fin.c_str());
		seqan::BamAlignmentRecord record;

		string line;
		vector<string> fields;
		while (!atEnd(bamStreamIn)){

	  		if (!hasFlagUnmapped(record)){

				readRecord(record, bamStreamIn);
				string chr = seqan::toCString(bamStreamIn.header.sequenceInfos[record.rId].i1);
				int start = record.pos;
				string seqString = seqan::toCString(record.seq);
				int end = record.pos + seqString.length();
				string strand = "+";
				if(hasFlagRC(record)){
					strand = "-";
				};

				string hashKey = chr + "\t" + boost::lexical_cast<string>(start) + "\t" + boost::lexical_cast<string>(end) + "\tfilteredClonalReads\t0\t" + strand;
				if(readMap.count(hashKey) == 0) {
					readMap[hashKey] = 1;
				} else {
					readMap[hashKey]++;
				}	
				totalReads++;
			}
		}
		seqan::close(bamStreamIn);

	}else if(ft == "bowtie"){

		ifstream in(fin.c_str());
		string line;

		while (getline(in, line)) {

			vector<string> fields;
			typedef boost::tokenizer<boost::char_separator<char> >
					tokenizer;
			boost::char_separator<char> sep("\t");
			tokenizer tokens(line, sep);
			fields.assign(tokens.begin(), tokens.end());

			string chr = fields[2];
			long int start = atol(fields[3].c_str());
			int readLength = fields[4].length();
			long int end = start + readLength;
			string strand = fields[1];
			
			string hashKey = chr + "\t" + boost::lexical_cast<string>(start) + "\t" + boost::lexical_cast<string>(end) + "\tfilteredClonalReads\t0\t" + strand;
			if(readMap.count(hashKey) == 0) {
				readMap[hashKey] = 1;
			} else {
				readMap[hashKey]++;
			}	
			totalReads++;
		}
		in.close();

	}else if(ft == "soap"){

		ifstream in(fin.c_str());
		string line;

		while (getline(in, line)) {

			vector<string> fields;
			typedef boost::tokenizer<boost::char_separator<char> >
					tokenizer;
			boost::char_separator<char> sep("\t");
			tokenizer tokens(line, sep);
			fields.assign(tokens.begin(), tokens.end());

			string chr = fields[2];
			long int start = atol(fields[3].c_str());
			int readLength = fields[4].length();
			long int end = start + readLength;
			string strand = fields[1];
			
			string hashKey = chr + "\t" + boost::lexical_cast<string>(start) + "\t" + boost::lexical_cast<string>(end) + "\tfilteredClonalReads\t0\t" + strand;
			if(readMap.count(hashKey) == 0) {
				readMap[hashKey] = 1;
			} else {
				readMap[hashKey]++;
			}	
			totalReads++;
		}
		in.close();

	}

	string fileOut;
	size_t pos = fin.find_last_of("/");
	if(pos != string::npos)
		fileOut.assign(fin.begin() + pos + 1, fin.end());
	else
		fileOut = fin;

	int lastIndex = fileOut.find_last_of("."); 
	std::string rawName = fileOut.substr(0, lastIndex); 
	string filePath = dirOut + "/" + rawName + "_" + nrDuplicates + "-duplicate_cleaned.bed";
	ofstream fout;
	fout.open(filePath.c_str());

	int readsWritten = 0;
	for (std::map<string, int>::iterator readIter = readMap.begin(); readIter != readMap.end(); ++readIter) {
		string line = readIter->first;
        int counts = readIter->second;
		if(counts > nrReadsToKeep){

			counts = nrReadsToKeep;
		}

		for(int i=0;i<counts;i++){

			//cout << line << endl;
			fout << line << endl;
			readsWritten++;
		}
	}

	int readsRemoved = totalReads - readsWritten;
	string sRemoved = boost::lexical_cast<string>(readsRemoved);

	double percRemoved = (readsRemoved/totalReads)*100.00;
	string sPercRemoved = boost::lexical_cast<string>(percRemoved);

	string sTotal = boost::lexical_cast<string>(totalReads);
	
	//cout << sRemoved << " clonal reads removed from a total of " << sTotal << " reads (" << sPercRemoved <<"%)." << endl;

	Rcpp::CharacterVector fileName(1);
	fileName[0] = filePath;
	return fileName;
}


std::vector<string> bam2wig(string fileIn, string sampleType, string outputDir, string chromName,
		int windowSize, unsigned int fragLength, bool split, string genomeBuild,
		bool zeros) {


	std::vector<string> results(2);

	string ws = boost::lexical_cast<string>(windowSize);
	string d = boost::lexical_cast<string>(fragLength);

	long unsigned int chromLength = getChromosomeLength(genomeBuild, chromName);

	static int nrOfBins = (chromLength / windowSize) + 1;

	std::string filename;
	size_t pos = fileIn.find_last_of("/");
	if(pos != std::string::npos)
		filename.assign(fileIn.begin() + pos + 1, fileIn.end());
	else
		filename = fileIn;

	int lastIndex = filename.find_last_of("."); 
	std::string rawName = filename.substr(0, lastIndex); 

	seqan::BamStream bamStreamIn(fileIn.c_str());
    seqan::BamAlignmentRecord record;
//    if (readRecord(record, bamStreamIn) != 0)
//    {
//        std::cerr << "ERROR: Could not read record!\n";
//        //return 1;
//    }

	if (split) {
		
		std::map<long int, int> binsFw;
		std::map<long int, int> binsRev;
		//cout << "Writing " << chromName << " to strand-specific wigs..." << endl;

		string line;
		vector<string> fields;
		while (!atEnd(bamStreamIn)){

				readRecord(record, bamStreamIn);

				string chr = seqan::toCString(bamStreamIn.header.sequenceInfos[record.rId].i1);

		  		if (!hasFlagUnmapped(record)){

					if (chr == chromName) {

						string seqString = seqan::toCString(record.seq);

						int start = record.pos;
						int end = record.pos + seqString.length();

						string strand = "+";
						if(hasFlagRC(record)){
							strand = "-";
						}

						int readLength = abs(start-end);

						int distance = fragLength - readLength;
				
						if(distance < 0){

							distance = 0;
						}

						if (strand == "+") {

							end = end + distance;

							unsigned long int lower = start / windowSize;
							unsigned long int upper = end / windowSize;

							for (; lower < upper; lower++) {

								if (binsFw.count(lower) > 0){ 

									binsFw[lower]++;
								} else {
									binsFw[lower] = 1;
								}
							}
						}

						if (strand == "-") {

							start = start - distance;
							unsigned long int lower = start / windowSize;
							unsigned long int upper = end / windowSize;
							for (; lower < upper; lower++) {

								if (binsRev.count(lower) > 0) {
									binsRev[lower]++;
								} else {
									binsRev[lower] = 1;
								}
							}
						}
					}

			}
		}

		ofstream foutFw;
		ofstream foutRev;

		string fileOutFw = outputDir + "/" + chromName + "_" + sampleType + "_" + rawName + "_res-" + ws + "_dist-" + d + "_fw.wig";
		string fileOutRev = outputDir + "/" + chromName + "_" + sampleType + "_" + rawName + "_res-" + ws + "_dist-" + d + "_rev.wig";

		results[0] = fileOutFw;
		results[1] = fileOutRev;

		//cout << "writing to " << fileOutFw << endl;
		//cout << "writing to " << fileOutRev << endl;

		foutFw.open(fileOutFw.c_str());
		foutRev.open(fileOutRev.c_str());

		foutFw
				<< "track type=wiggle_0 name=\"peak intensities\" description=\"Peak intensities "
					"forward strand\"\nvariableStep chrom=" << chromName
				<< " span=" << windowSize << endl;

		foutRev
				<< "track type=wiggle_0 name=\"peak intensities\" description=\"Peak intensities "
					"reverse strand\"\nvariableStep chrom=" << chromName
				<< " span=" << windowSize << endl;

		if (zeros) {

			for (unsigned int i = 0; i < chromLength; i=i+windowSize) {

				long int positionFw = i;//(i * windowSize) ;
				if (binsFw.count(positionFw) > 0) {
					string pos = boost::lexical_cast<string>(positionFw);
					string ity = boost::lexical_cast<string>(binsFw[i]);
					foutFw << positionFw << "\t" << ity << endl;
				}else{

					foutFw << positionFw << "\t" << 0 << endl;
				}

				long int positionRev = i;  // (i * windowSize) ;
				if (binsRev.count(positionRev) > 0) {
					string pos = boost::lexical_cast<string>(positionRev);
					string ity = boost::lexical_cast<string>(binsRev[i]);
					foutRev << positionRev << "\t" << ity << endl;
				}else{

					foutRev << positionRev << "\t" << 0 << endl;
				}

			}

		} else {

			for (std::map<long int, int>::iterator iterFw = binsFw.begin(); iterFw
					!= binsFw.end(); ++iterFw) {
				long int binPosition = iterFw->first;
                                long int position = binPosition * windowSize;
				int intensity = iterFw->second;
				string pos = boost::lexical_cast<string>(position);
				string ity = boost::lexical_cast<string>(intensity);
				foutFw << pos << "\t" << ity << endl;
			}

			for (std::map<long int, int>::iterator iterRev = binsRev.begin(); iterRev
					!= binsRev.end(); ++iterRev) {
				long int binPosition = iterRev->first;
                                long int position = binPosition * windowSize;
				int intensity = iterRev->second;
				string pos = boost::lexical_cast<string>(position);
				string ity = boost::lexical_cast<string>(intensity);
				foutRev << pos << "\t" << ity << endl;
			}
		}
		foutFw.close();
		foutRev.close();

	} 

	else {

		std::map<long int, int> bins;
		//cout << "Writing " << chromName << " wig..." << endl;

		string line;
		vector<string> fields;
		while (!atEnd(bamStreamIn)){

	  		if (!hasFlagUnmapped(record)){

				readRecord(record, bamStreamIn);

				string chr = seqan::toCString(bamStreamIn.header.sequenceInfos[record.rId].i1);

				if (chr == chromName) {

					string seqString = seqan::toCString(record.seq);

					int start = record.pos;
					int end = record.pos + seqString.length();

					string strand = "+";
					if(hasFlagRC(record)){
						strand = "-";
					}

					int readLength = abs(start-end);

					int distance = fragLength - readLength;
				
					if(distance < 0){

						distance = 0;
					}

					unsigned long int lower;
					unsigned long int upper;
					if (strand == "+") {

						end = end + distance;
						lower = start / windowSize;
						upper = end / windowSize;

					}

					if (strand == "-") {

						start = start - distance;
						lower = start / windowSize;
						upper = end / windowSize;
					}

					for (; lower < upper; lower++) {

						if (bins.count(lower) > 0) {
							bins[lower]++;
						} else {
							bins[lower] = 1;
						}
					}
				}
			}
		}

		ofstream fout;
		string fileOut = outputDir + "/" + chromName + "_" + sampleType + "_" + rawName + "_res-" + ws + "_dist-" + d + "_both.wig";
		results[0] = fileOut;
		results[1] = "N/A";
		fout.open(fileOut.c_str());

		fout
				<< "track type=wiggle_0 name=\"peak intensities\" description=\"Peak intensities "
					"both strands\"\nvariableStep chrom=" << chromName
				<< " span=" << windowSize << endl;

		if (zeros) {

			for (int i = 0; i < nrOfBins; i++) {

				long int position = (i * windowSize) + 1;
				if (bins.count(position) > 0) {
					string pos = boost::lexical_cast<string>(position);
					string ity = boost::lexical_cast<string>(bins[i]);
					fout << position << "\t" << ity << endl;
				}
			}

		} else {

			for (std::map<long int, int>::iterator iter = bins.begin(); iter
					!= bins.end(); ++iter) {

				long int binPosition = iter->first;
                                long int position = binPosition * windowSize;
				int intensity = iter->second;
				string pos = boost::lexical_cast<string>(position);
				string ity = boost::lexical_cast<string>(intensity);
				fout << pos << "\t" << ity << endl;
				
			}
		}
		fout.close();
	}
	return results;
}


bool checkFormat(string fin, string format){

	bool wellformed = false;

	//"bed3"or "bed6"
	if(format=="bed"){
		int counter = 1;
		ifstream in(fin.c_str());
		string line;

		while ((getline(in, line)) && counter <=3) {

			counter++;
			if( line.substr(0, 1)!="#"){

				vector<string> fields;
				typedef boost::tokenizer<boost::char_separator<char> >
						tokenizer;
				boost::char_separator<char> sep("\t");
				tokenizer tokens(line, sep);
				fields.assign(tokens.begin(), tokens.end());

				if(fields.size()>=3){

					string chr = fields[0];
					if(chr.substr (0,3)!="chr"){

						//cout << "WARNING: Chromosome names in BED files don't seem to be properly formatted" << endl;
						wellformed == false;
						//return wellformed;
					}else{
						wellformed == true;
						//return wellformed;
					}
					if(fields.size()==4 or fields.size()==5){

						wellformed == false;
					} 

					if(fields.size()>=6){

						string strand = fields[5];
						if (strand!="+" and strand!="-"){

							//cout << "WARNING: Strand sign should either be \"+\" or \"-\"" << endl;
							wellformed = false;
							//return wellformed;
						}else{

							wellformed = true;
						}
					}

				}else{

					//cout << "WARNING: Bed file should at least contain 3 columns" << endl;
					wellformed = false;
					//return wellformed;
				}
			}
		}
	}
	return wellformed;
}


std::vector<string> map2wig(string fileIn, string sampleType, string outputDir, string chromName,
		int windowSize, unsigned int fragLength, bool split, string genomeBuild,
		bool zeros, string format) {

	std::vector<string> results(2);

	bool wellformed = checkFormat(fileIn,format);
	if(!wellformed)
	{
		//cout << "wig " << fileIn << " is not well-formed!" << endl;
		results[0] = "ILLFORMED";
		results[1] = "ILLFORMED";
	}
	else
	{
		//cout << "wig " << fileIn << " seems to be well-formed!" << endl;
		string ws = boost::lexical_cast<string>(windowSize);
		string d = boost::lexical_cast<string>(fragLength);

		long unsigned int chromLength = getChromosomeLength(genomeBuild, chromName);

		static long unsigned int nrOfBins = (chromLength / windowSize) + 1;
		//static int nrOfBins = (chromLength / windowSize) + 1;

		std::string filename;
		size_t pos = fileIn.find_last_of("/");
		if(pos != std::string::npos)
		{
			filename.assign(fileIn.begin() + pos + 1, fileIn.end());
		}
		else{
			filename = fileIn;
		}
		int lastIndex = filename.find_last_of("."); 
		std::string rawName = filename.substr(0, lastIndex); 

		ifstream in(fileIn.c_str());

		if (split) 
		{

			std::map<long int, int> binsFw;
			std::map<long int, int> binsRev;
			//cout << "Writing " << chromName << " to strand-specific wigs..." << endl;

			string line;
			vector<string> fields;
			while (getline(in, line)) 
			{

				typedef boost::tokenizer<boost::char_separator<char> >
						tokenizer;
				boost::char_separator<char> sep("\t");
				tokenizer tokens(line, sep);
				fields.assign(tokens.begin(), tokens.end());

				//TODO: convert to lower case
				if(format == "bed")
				{

					if (fields[0] == chromName) {
						long int start = atol(fields[1].c_str());
						long int end = atol(fields[2].c_str());

						int readLength = abs(start-end);
						int distance = fragLength - readLength;
						if(distance < 0){

							distance = 0;
						}

						if (fields[5] == "+") {
							//start = atol(fields[1].c_str());
							end = end + distance;

							unsigned long int lower = start / windowSize;
							unsigned long int upper = end / windowSize;

							for (; lower < upper; lower++) {

								if (binsFw.count(lower) > 0){ 
									binsFw[lower]++;
								} else {
									binsFw[lower] = 1;
								}
							}
						}

						if (fields[5] == "-") {
							start = start - distance;
							//end = atol(fields[2].c_str());
							unsigned long int lower = start / windowSize;
							unsigned long int upper = end / windowSize;
							for (; lower < upper; lower++) {

								if (binsRev.count(lower) > 0) {
									binsRev[lower]++;
								} else {
									binsRev[lower] = 1;
								}
							}
						}
					}
				}
				else
				{

					if(format == "bowtie")
					{

						if (fields[2] == chromName) {
							long int start = atol(fields[3].c_str());
							int readLength = fields[4].length();
							long int end = start + readLength;
							int distance = fragLength - readLength;
							if(distance < 0){

								distance = 0;
							}

							if (fields[1] == "+") {
								//start = atol(fields[1].c_str());
								end = end + distance;

								unsigned long int lower = start / windowSize;
								unsigned long int upper = end / windowSize;

								for (; lower < upper; lower++) {

									if (binsFw.count(lower) > 0){ 
										binsFw[lower]++;
									} else {
										binsFw[lower] = 1;
									}
								}
							}

							if (fields[1] == "-") {
								start = start - distance;
								//end = atol(fields[2].c_str());
								unsigned long int lower = start / windowSize;
								unsigned long int upper = end / windowSize;
								for (; lower < upper; lower++) {

									if (binsRev.count(lower) > 0) {
										binsRev[lower]++;
									} else {
										binsRev[lower] = 1;
									}
								}
							}
						}

					}
					else
					{

						if(format == "soap")
						{

							if (fields[7] == chromName) {
								long int start = atol(fields[8].c_str());
								int readLength = fields[1].length();
								long int end = start + readLength;
								int distance = fragLength - readLength;
								if(distance < 0){

									distance = 0;
								}

								if (fields[6] == "+") {
									//start = atol(fields[1].c_str());
									end = end + distance;

									unsigned long int lower = start / windowSize;
									unsigned long int upper = end / windowSize;

									for (; lower < upper; lower++) {

										if (binsFw.count(lower) > 0){ 
											binsFw[lower]++;
										} else {
											binsFw[lower] = 1;
										}
									}
								}

								if (fields[6] == "-") {
									start = start - distance;
									//end = atol(fields[2].c_str());
									unsigned long int lower = start / windowSize;
									unsigned long int upper = end / windowSize;
									for (; lower < upper; lower++) {

										if (binsRev.count(lower) > 0) {
											binsRev[lower]++;
										} else {
											binsRev[lower] = 1;
										}
									}
								}
							}
						}
					}
				}
			}

			ofstream foutFw;
			ofstream foutRev;

			string fileOutFw = outputDir + "/" + chromName + "_" + sampleType + "_" + rawName + "_res-" + ws + "_dist-" + d + "_fw.wig";
			string fileOutRev = outputDir + "/" + chromName + "_" + sampleType + "_" + rawName + "_res-" + ws + "_dist-" + d + "_rev.wig";

			results[0] = fileOutFw;
			results[1] = fileOutRev;

			//cout << "writing to " << fileOutFw << endl;
			//cout << "writing to " << fileOutRev << endl;

			foutFw.open(fileOutFw.c_str());
			foutRev.open(fileOutRev.c_str());

			foutFw
					<< "track type=wiggle_0 name=\"peak intensities\" description=\"Peak intensities "
						"forward strand\"\nvariableStep chrom=" << chromName
					<< " span=" << windowSize << endl;

			foutRev
					<< "track type=wiggle_0 name=\"peak intensities\" description=\"Peak intensities "
						"reverse strand\"\nvariableStep chrom=" << chromName
					<< " span=" << windowSize << endl;

			if (zeros) 
			{

				for (unsigned int i = 0; i < chromLength; i=i+windowSize) {

					long int positionFw = i;//(i * windowSize) ;
					if (binsFw.count(positionFw) > 0) {
						string pos = boost::lexical_cast<string>(positionFw);
						string ity = boost::lexical_cast<string>(binsFw[i]);
						foutFw << positionFw << "\t" << ity << endl;
					}else{

						foutFw << positionFw << "\t" << 0 << endl;
					}

					long int positionRev = i;  // (i * windowSize) ;
					if (binsRev.count(positionRev) > 0) {
						string pos = boost::lexical_cast<string>(positionRev);
						string ity = boost::lexical_cast<string>(binsRev[i]);
						foutRev << positionRev << "\t" << ity << endl;
					}else{

						foutRev << positionRev << "\t" << 0 << endl;
					}

				}

			} 
			else 
			{

				for (std::map<long int, int>::iterator iterFw = binsFw.begin(); iterFw
						!= binsFw.end(); ++iterFw) {
					long int binPosition = iterFw->first;
		                            long int position = binPosition * windowSize;
					int intensity = iterFw->second;
					string pos = boost::lexical_cast<string>(position);
					string ity = boost::lexical_cast<string>(intensity);
					foutFw << pos << "\t" << ity << endl;
				}

				for (std::map<long int, int>::iterator iterRev = binsRev.begin(); iterRev
						!= binsRev.end(); ++iterRev) {
					long int binPosition = iterRev->first;
		                            long int position = binPosition * windowSize;
					int intensity = iterRev->second;
					string pos = boost::lexical_cast<string>(position);
					string ity = boost::lexical_cast<string>(intensity);
					foutRev << pos << "\t" << ity << endl;
				}
			}
			foutFw.close();
			foutRev.close();

		} 
		else 
		{ 
			std::map<long int, int> bins;
			//cout << "Writing " << chromName << " wig..." << endl;

			string line;
			vector<string> fields;
			while (getline(in, line)) 
			{
				typedef boost::tokenizer<boost::char_separator<char> >
						tokenizer;
				boost::char_separator<char> sep("\t");
				tokenizer tokens(line, sep);
				fields.assign(tokens.begin(), tokens.end());

				//TODO: convert to lower case
				if(format == "bed")
				{

					if (fields[0] == chromName) 
					{
						long int start = atol(fields[1].c_str());
						long int end = atol(fields[2].c_str());

						int readLength = abs(start-end);
						int distance = fragLength - readLength;
						if(distance < 0)
						{

							distance = 0;
						}

						if (fields[5] == "+") 
						{
							//start = atol(fields[1].c_str());
							end = end + distance;

							unsigned long int lower = start / windowSize;
							unsigned long int upper = end / windowSize;

							for (; lower < upper; lower++) 
							{

								if (bins.count(lower) > 0)
								{ 
									bins[lower]++;
								} else {
									bins[lower] = 1;
								}
							}
						}

						if (fields[5] == "-") 
						{
							start = start - distance;
							//end = atol(fields[2].c_str());
							unsigned long int lower = start / windowSize;
							unsigned long int upper = end / windowSize;
							for (; lower < upper; lower++) {

								if (bins.count(lower) > 0) {
									bins[lower]++;
								} else {
									bins[lower] = 1;
								}
							}
						}
					}
				}
				else
				{
					if(format == "bowtie")
					{
						if (fields[2] == chromName) 
						{
							long int start = atol(fields[3].c_str());
							int readLength = fields[4].length();
							long int end = start + readLength;
							int distance = fragLength - readLength;
							if(distance < 0)
							{

								distance = 0;
							}

							if (fields[1] == "+")
							{
								//start = atol(fields[1].c_str());
								end = end + distance;

								unsigned long int lower = start / windowSize;
								unsigned long int upper = end / windowSize;

								for (; lower < upper; lower++) 
								{

									if (bins.count(lower) > 0)
									{ 
										bins[lower]++;
									} else {
										bins[lower] = 1;
									}
								}
							}
							if (fields[1] == "-")
							{
								start = start - distance;
								//end = atol(fields[2].c_str());
								unsigned long int lower = start / windowSize;
								unsigned long int upper = end / windowSize;
								for (; lower < upper; lower++) 
								{

									if (bins.count(lower) > 0) {
										bins[lower]++;
									} else {
										bins[lower] = 1;
									}
								}
							}
						}
					}
					else
					{
						if(format == "soap")
						{
							if (fields[7] == chromName) 
							{
								long int start = atol(fields[8].c_str());
								int readLength = fields[1].length();
								long int end = start + readLength;
								int distance = fragLength - readLength;
								if(distance < 0)
								{
									distance = 0;
								}

								if (fields[6] == "+") 
								{
									//start = atol(fields[1].c_str());
									end = end + distance;

									unsigned long int lower = start / windowSize;
									unsigned long int upper = end / windowSize;

									for (; lower < upper; lower++) {

										if (bins.count(lower) > 0){ 
											bins[lower]++;
										} else {
											bins[lower] = 1;
										}
									}
								}

								if (fields[6] == "-") 
								{
									start = start - distance;
									//end = atol(fields[2].c_str());
									unsigned long int lower = start / windowSize;
									unsigned long int upper = end / windowSize;
									for (; lower < upper; lower++) {

										if (bins.count(lower) > 0) {
											bins[lower]++;
										} else {
											bins[lower] = 1;
										}
									}
								}
							}
						}
					}
				}
			}
			ofstream fout;
			string fileOut = outputDir + "/" + chromName + "_" + sampleType + "_" + rawName + "_res-" + ws + "_dist-" + d + "_both.wig";
			results[0] = fileOut;
			results[1] = "N/A";
			fout.open(fileOut.c_str());
			fout
					<< "track type=wiggle_0 name=\"peak intensities\" description=\"Peak intensities "
						"both strands\"\nvariableStep chrom=" << chromName
					<< " span=" << windowSize << endl;

			if (zeros) {

				for (int i = 0; i < nrOfBins; i++) {

					long int position = (i * windowSize) + 1;
//					fout << boost::lexical_cast<string>(i) << "\t" << boost::lexical_cast<string>(position);
					if (bins.count(position) > 0) {

						string pos = boost::lexical_cast<string>(position);
						string ity = boost::lexical_cast<string>(bins[i]);
//						fout << "\t" << boost::lexical_cast<string>(bins.count(position)) << "\t" << pos  << "\t" << ity << endl;
						fout << position << "\t" << ity << endl;
					}
//					else{
//						fout << endl;
//					}
					else{

						fout << position << "\t" << 0 << endl;
					}
				}
			} 
			else 
			{

				for (std::map<long int, int>::iterator iter = bins.begin(); iter
						!= bins.end(); ++iter) {
					long int binPosition = iter->first;
		                            long int position = binPosition * windowSize;
					int intensity = iter->second;
					string pos = boost::lexical_cast<string>(position);
					string ity = boost::lexical_cast<string>(intensity);
					fout << pos << "\t" << ity << endl;
				}
			}
			fout.close();
		}
	}
	return results;
}


//std::vector<string> chromosome2wig(string fileIn, string sampleType, string outputDir, string chromName,
//		int windowSize, unsigned int fragLength, bool split, string genomeBuild,
//		bool zeros) {

//	std::vector<string> results(2);

//	string ws = boost::lexical_cast<string>(windowSize);
//	string d = boost::lexical_cast<string>(fragLength);

//	long unsigned int chromLength = getChromosomeLength(genomeBuild, chromName);

//	static int nrOfBins = (chromLength / windowSize) + 1;

//	std::string filename;
//	size_t pos = fileIn.find_last_of("/");
//	if(pos != std::string::npos)
//		filename.assign(fileIn.begin() + pos + 1, fileIn.end());
//	else
//		filename = fileIn;

//	int lastIndex = filename.find_last_of("."); 
//	std::string rawName = filename.substr(0, lastIndex); 

//	ifstream in(fileIn.c_str());

//	if (split) {
//		
//		std::map<long int, int> binsFw;
//		std::map<long int, int> binsRev;
//		cout << "Writing " << chromName << " to strand-specific wigs..."
//				<< endl;
//		//if (distance > 0) {

//		string line;
//		vector<string> fields;
//		while (getline(in, line)) {

//			typedef boost::tokenizer<boost::char_separator<char> >
//					tokenizer;
//			boost::char_separator<char> sep("\t");
//			tokenizer tokens(line, sep);
//			fields.assign(tokens.begin(), tokens.end());

//			if (fields[0] == chromName) {
//				long int start = atol(fields[1].c_str());
//				long int end = atol(fields[2].c_str());

//				int readLength = abs(start-end);
//				int distance = fragLength - readLength;
//				if(distance < 0){

//					distance = 0;
//				}

//				if (fields[5] == "+") {
//					start = atol(fields[1].c_str());
//					end = end + distance;

//					unsigned long int lower = start / windowSize;
//					unsigned long int upper = end / windowSize;

//					for (; lower < upper; lower++) {

//						//if (binsFw.count(lower > 0)) {
//						if (binsFw.count(lower) > 0){ 
//							binsFw[lower]++;

//							//cout << "check increment forward" << endl;
//						} else {
//							binsFw[lower] = 1;
//							//cout << "check new bin forward" << endl;
//						}
//					}
//				}

//				if (fields[5] == "-") {
//					start = start - distance;
//					end = atol(fields[2].c_str());
//					unsigned long int lower = start / windowSize;
//					unsigned long int upper = end / windowSize;
//					for (; lower < upper; lower++) {

//						if (binsRev.count(lower) > 0) {
//							binsRev[lower]++;
//							//cout << "check increment reverse" << endl;
//						} else {
//							binsRev[lower] = 1;
//							//cout << "check new bin forward" << endl;
//						}
//					}
//				}
//			}
//		}

//		ofstream foutFw;
//		ofstream foutRev;

//		string fileOutFw = outputDir + "/" + chromName + "_" + sampleType + "_" + rawName + "_res-" + ws + "_dist-" + d + "_fw.wig";
//		string fileOutRev = outputDir + "/" + chromName + "_" + sampleType + "_" + rawName + "_res-" + ws + "_dist-" + d + "_rev.wig";

//		results[0] = fileOutFw;
//		results[1] = fileOutRev;

//		cout << "writing to " << fileOutFw << endl;
//		cout << "writing to " << fileOutRev << endl;

//		foutFw.open(fileOutFw.c_str());
//		foutRev.open(fileOutRev.c_str());

//		foutFw
//				<< "track type=wiggle_0 name=\"peak intensities\" description=\"Peak intensities "
//					"forward strand\"\nvariableStep chrom=" << chromName
//				<< " span=" << windowSize << endl;

//		foutRev
//				<< "track type=wiggle_0 name=\"peak intensities\" description=\"Peak intensities "
//					"reverse strand\"\nvariableStep chrom=" << chromName
//				<< " span=" << windowSize << endl;

//		if (zeros) {

//			for (unsigned int i = 0; i < chromLength; i=i+windowSize) {

//				long int positionFw = i;//(i * windowSize) ;
//				if (binsFw.count(positionFw) > 0) {
//					string pos = boost::lexical_cast<string>(positionFw);
//					string ity = boost::lexical_cast<string>(binsFw[i]);
//					foutFw << positionFw << "\t" << ity << endl;
//				}else{

//					foutFw << positionFw << "\t" << 0 << endl;
//				}

//				long int positionRev = i;  // (i * windowSize) ;
//				if (binsRev.count(positionRev) > 0) {
//					string pos = boost::lexical_cast<string>(positionRev);
//					string ity = boost::lexical_cast<string>(binsRev[i]);
//					foutRev << positionRev << "\t" << ity << endl;
//				}else{

//					foutRev << positionRev << "\t" << 0 << endl;
//				}

//			}

//		} else {

//			for (std::map<long int, int>::iterator iterFw = binsFw.begin(); iterFw
//					!= binsFw.end(); ++iterFw) {
//				long int binPosition = iterFw->first;
//                                long int position = binPosition * windowSize;
//				int intensity = iterFw->second;
//				string pos = boost::lexical_cast<string>(position);
//				string ity = boost::lexical_cast<string>(intensity);
//				foutFw << pos << "\t" << ity << endl;
//			}

//			for (std::map<long int, int>::iterator iterRev = binsRev.begin(); iterRev
//					!= binsRev.end(); ++iterRev) {
//				long int binPosition = iterRev->first;
//                                long int position = binPosition * windowSize;
//				int intensity = iterRev->second;
//				string pos = boost::lexical_cast<string>(position);
//				string ity = boost::lexical_cast<string>(intensity);
//				foutRev << pos << "\t" << ity << endl;
//			}
//		}
//		foutFw.close();
//		foutRev.close();

//	} else {

//		std::map<long int, int> bins;
//		cout << "Writing " << chromName << " wig..." << endl;
//		//if (distance > 0) {

//		string line;
//		vector<string> fields;
//		while (getline(in, line)) {

//			typedef boost::tokenizer<boost::char_separator<char> >
//					tokenizer;
//			boost::char_separator<char> sep("\t");
//			tokenizer tokens(line, sep);
//			fields.assign(tokens.begin(), tokens.end());
//			if (fields[0] == chromName) {

//				long int start = atol(fields[1].c_str());
//				long int end = atol(fields[2].c_str());

//				int readLength = abs(start-end);
//				int distance = fragLength - readLength;
//				
//				if(distance < 0){

//					distance = 0;
//				}

//				unsigned long int lower;
//				unsigned long int upper;
//				if (fields[5] == "+") {

//					end = end + distance;
//					lower = start / windowSize;
//					upper = end / windowSize;

//				}

//				if (fields[5] == "-") {

//					start = start - distance;
//					lower = start / windowSize;
//					upper = end / windowSize;
//				}

//				for (; lower < upper; lower++) {

//					if (bins.count(lower) > 0) {
//						bins[lower]++;
//					} else {
//						bins[lower] = 1;
//					}
//				}
//			}
//		}

//		ofstream fout;
//		string fileOut = outputDir + "/" + chromName + "_" + sampleType + "_" + rawName + "_res-" + ws + "_dist-" + d + "_both.wig";
//		results[0] = fileOut;
//		results[1] = "N/A";
//		fout.open(fileOut.c_str());

//		fout
//				<< "track type=wiggle_0 name=\"peak intensities\" description=\"Peak intensities "
//					"both strands\"\nvariableStep chrom=" << chromName
//				<< " span=" << windowSize << endl;

//		if (zeros) {

//			for (int i = 0; i < nrOfBins; i++) {

//				long int position = (i * windowSize) + 1;
//				if (bins.count(position) > 0) {
//					string pos = boost::lexical_cast<string>(position);
//					string ity = boost::lexical_cast<string>(bins[i]);
//					fout << position << "\t" << ity << endl;
//				}			//added
//				else{

//					fout << position << "\t0" << endl;
//				}
//			}


//		} else {

//			for (std::map<long int, int>::iterator iter = bins.begin(); iter
//					!= bins.end(); ++iter) {
//				long int binPosition = iter->first;
//                                long int position = binPosition * windowSize;
//				int intensity = iter->second;
//				string pos = boost::lexical_cast<string>(position);
//				string ity = boost::lexical_cast<string>(intensity);
//				fout << pos << "\t" << ity << endl;
//			}
//		}
//		fout.close();
//	}
//	return results;
//}


SEXP wigWriter(SEXP fileName, SEXP sampleType, SEXP outputDir,SEXP chromName, SEXP resolution, SEXP elongation, SEXP split, SEXP genomeBuild, SEXP zeros, SEXP format) {

	using namespace Rcpp;
	string fn = as<string>(fileName);
	string od = as<string>(outputDir);
	string cn = as<string>(chromName);
	string st = as<string>(sampleType);
	string ft = as<string>(format);
	int re = as<int>(resolution);
	int el = as<int>(elongation);
	bool sp = as<bool>(split);
	string gb = as<string>(genomeBuild);
	bool zs = as<bool>(zeros);

	std::vector<string> results;
	if((ft == "sam") || (ft == "bam")){

		results = bam2wig(fn, st, od, cn, re, el, sp, gb, zs);	
	}

	if((ft == "bowtie") || (ft == "bed") || (ft == "eland") || (ft == "soap")){

		results = map2wig(fn, st, od, cn, re, el, sp, gb, zs, ft );
	}
	Rcpp::CharacterVector filePathVec(2);
	filePathVec[0] = results[0];
	filePathVec[1] = results[1];
	return filePathVec;
}
