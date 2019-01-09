#include <iostream>
#include <iomanip>
#include <istream>
#include <sstream>
#include <fstream>
#include "Alignment.h"



Alignment::Alignment(std::string fileName) {

    std::cout << "   * Reading data file \"" << fileName << "\"" << std::endl;

	/* open the file */
	std::ifstream seqStream(fileName.c_str());
	if (!seqStream) 
		{
		std::cerr << "Cannot open file \"" + fileName + "\"" << std::endl;
		exit(1);
		}

	std::string linestring = "";
	int line = 0;
	std::string theSequence = "";
	int taxonNum = 0;
	matrix = NULL;
	numTaxa = numChar = 0;
	bool excludeLine = false, charSetLine = false;
	bool* tempVec = NULL;
	int pid = 0;
	while( getline(seqStream, linestring).good() )
		{
		std::istringstream linestream(linestring);
		int ch;
		std::string word = "";
		int wordNum = 0;
		int siteNum = 0;
		excludeLine = false;
		charSetLine = false;
		std::string cmdString = "";
        bool foundEqualSign = false;
		do
			{
			word = "";
			linestream >> word;
			wordNum++;
			if (line == 0)
				{
				/* read the number of taxa/chars from the first line */
				int x;
				std::istringstream buf(word);
				buf >> x;
				if (wordNum == 1)
					numTaxa = x;
				else
					numChar = numSitePatterns = x;
				if (numTaxa > 0 && numChar > 0 && matrix == NULL)
					{	
					matrix = new int*[numTaxa];
					matrix[0] = new int[numTaxa * numChar];
					for (size_t i=1; i<numTaxa; i++)
						matrix[i] = matrix[i-1] + numChar;
					for (size_t i=0; i<numTaxa; i++)
						for (size_t j=0; j<numChar; j++)
							matrix[i][j] = 0;
					isExcluded = new bool[numChar];
					partitionId = new size_t[numChar];
					for (size_t i=0; i<numChar; i++)
						{
						isExcluded[i] = false;
						partitionId[i] = -1;
						}
					tempVec = new bool[numChar];
                    std::cout << "   * Analysis has " << numTaxa << " taxa and " << numChar << " characters" << std::endl;
					}
				}
			else
				{
				if (wordNum == 1)
					{
					if ( word == "exclude" )
                        {
						excludeLine = true;
                        foundEqualSign = false;
                        }
					else if ( word == "charset" )
                        {
						charSetLine = true;
                        foundEqualSign = false;
                        }
					else
						{
						taxonNames.push_back(word);
						taxonNum++;
						}
					}
				else
					{
					if (excludeLine == true || charSetLine == true)
						{
						for (int i=0; i<word.size(); i++)
							{
                            if (word.at(i) == '=')
                                foundEqualSign = true;
                            if (foundEqualSign == true)
                                {
                                if ( isdigit(word.at(i)) )
                                    cmdString += word.at(i);
                                else if ( word.at(i) == '-' )
                                    cmdString += " - ";
                                else if (word.at(i) == '\\')
                                    cmdString += " \\ ";
                                }
							}
						cmdString += " ";
						for (int i=0; i<word.size(); i++)
							{
							if ( word.at(i) == ';' )
								{
								interpretString(cmdString, tempVec, numChar);
								if (excludeLine == true)
									{
									for (int i=0; i<numChar; i++)
										if (tempVec[i] == true)
											isExcluded[i] = true;
									}
								else if (charSetLine == true)
									{
									pid++;
									for (int i=0; i<numChar; i++)
										if (tempVec[i] == true)
											partitionId[i] = pid;
									}
								}
							}
						}
					else
						{
						for (int i=0; i<word.length(); i++)
							{
							char site = word.at(i);
							matrix[taxonNum-1][siteNum++] = nucID(site);
							}
						}
					}
				}
			} while ( (ch=linestream.get()) != EOF );
			
		// NOTE: We probably do not need this bit of code.
		if (line == 0)
			{
			/* the first line should contain the number of taxa and the sequence length */
			std::istringstream buf(word);
			//buf >> genomeSize;
			}
		else
			{
			for (int i=0; i<word.length(); i++)
				{
				char site = word.at(i);
				if (tolower(site) == 'a' || tolower(site) == 'c' || tolower(site) == 'g' || tolower(site) == 't')
					theSequence += tolower(site);
				}
			}
		line++;
		}	

	delete [] tempVec;
	
	/* close the file */
	seqStream.close();
}

Alignment::~Alignment(void) {

	delete [] matrix[0];
	delete [] matrix;
	delete [] isExcluded;
	delete [] partitionId;
}

int Alignment::getNucleotide(size_t i, size_t j) {

    return matrix[i][j];
}

std::vector<std::string> Alignment::getTaxonNames(void) {

    return taxonNames;
}

/*-------------------------------------------------------------------
|
|   GetPossibleNucs: 
|
|   This function initializes a vector, nuc[MAX_NUM_STATES]. The four elements
|   of nuc correspond to the four nucleotides in alphabetical order.
|   We are assuming that the nucCode is a binary representation of
|   the nucleotides that are consistent with the observation. For
|   example, if we observe an A, then the nucCode is 1 and the 
|   function initalizes nuc[0] = 1 and the other elements of nuc
|   to be 0.
|
|   Observation    nucCode        nuc
|        A            1           1000
|        C            2           0100
|        G            4           0010
|        T            8           0001
|        R            5           1010
|        Y           10           0101
|        M            3           1100
|        K           12           0011
|        S            6           0110
|        W            9           1001
|        H           11           1101
|        B           14           0111
|        V            7           1110
|        D           13           1011
|        N - ?       15           1111
|
-------------------------------------------------------------------*/
void Alignment::getPossibleNucs (int nucCode, int* nuc) {

	if (nucCode == 1)
		{
		nuc[0] = 1;
		nuc[1] = 0;
		nuc[2] = 0;
		nuc[3] = 0;
		}
	else if (nucCode == 2)
		{
		nuc[0] = 0;
		nuc[1] = 1;
		nuc[2] = 0;
		nuc[3] = 0;
		}
	else if (nucCode == 3)
		{
		nuc[0] = 1;
		nuc[1] = 1;
		nuc[2] = 0;
		nuc[3] = 0;
		}
	else if (nucCode == 4)
		{
		nuc[0] = 0;
		nuc[1] = 0;
		nuc[2] = 1;
		nuc[3] = 0;
		}
	else if (nucCode == 5)
		{
		nuc[0] = 1;
		nuc[1] = 0;
		nuc[2] = 1;
		nuc[3] = 0;
		}
	else if (nucCode == 6)
		{
		nuc[0] = 0;
		nuc[1] = 1;
		nuc[2] = 1;
		nuc[3] = 0;
		}
	else if (nucCode == 7)
		{
		nuc[0] = 1;
		nuc[1] = 1;
		nuc[2] = 1;
		nuc[3] = 0;
		}
	else if (nucCode == 8)
		{
		nuc[0] = 0;
		nuc[1] = 0;
		nuc[2] = 0;
		nuc[3] = 1;
		}
	else if (nucCode == 9)
		{
		nuc[0] = 1;
		nuc[1] = 0;
		nuc[2] = 0;
		nuc[3] = 1;
		}
	else if (nucCode == 10)
		{
		nuc[0] = 0;
		nuc[1] = 1;
		nuc[2] = 0;
		nuc[3] = 1;
		}
	else if (nucCode == 11)
		{
		nuc[0] = 1;
		nuc[1] = 1;
		nuc[2] = 0;
		nuc[3] = 1;
		}
	else if (nucCode == 12)
		{
		nuc[0] = 0;
		nuc[1] = 0;
		nuc[2] = 1;
		nuc[3] = 1;
		}
	else if (nucCode == 13)
		{
		nuc[0] = 1;
		nuc[1] = 0;
		nuc[2] = 1;
		nuc[3] = 1;
		}
	else if (nucCode == 14)
		{
		nuc[0] = 0;
		nuc[1] = 1;
		nuc[2] = 1;
		nuc[3] = 1;
		}
	else if (nucCode == 15)
		{
		nuc[0] = 1;
		nuc[1] = 1;
		nuc[2] = 1;
		nuc[3] = 1;
		}
	else if (nucCode == 16)
		{
		nuc[0] = 1;
		nuc[1] = 1;
		nuc[2] = 1;
		nuc[3] = 1;
		}
}

void Alignment::listTaxa(void) {

	int i = 1;
	for (std::vector<std::string>::iterator p=taxonNames.begin(); p != taxonNames.end(); p++)
		std::cout << std::setw(4) << i++ << " -- " << (*p) << '\n';
}

size_t Alignment::getNumSubsets(void) {

	size_t largestId = 0;
	for (size_t i=0; i<numChar; i++)
		{
		if (partitionId[i] > largestId)
			{
			largestId = partitionId[i];
			}
		}
		
	bool* isPartHere = new bool[largestId];
	for (size_t i=0; i<largestId; i++)
		isPartHere[i] = false;
	for (size_t i=0; i<numChar; i++)
		{
		if (isExcluded[i] == false)
			{
			isPartHere[ partitionId[i]-1 ] = true;
			}
		}
	size_t numParts = 0;
	for (size_t i=0; i<largestId; i++)
		if (isPartHere[i] == true)
			numParts++;
	delete [] isPartHere;
		
	return numParts;
}

std::string Alignment::getTaxonName(int i) {

	return taxonNames[i];
}

int Alignment::getTaxonIndex(std::string ns) {

	int taxonIndex = -1;
	int i = 0;
	for (std::vector<std::string>::iterator p=taxonNames.begin(); p != taxonNames.end(); p++)
		{
		if ( (*p) == ns )
			{
			taxonIndex = i;
			break;
			}
		i++;
		}
	return taxonIndex;
}

void Alignment::interpretString(std::string s, bool* v, int n) {

	for (size_t i=0; i<n; i++)
		v[i] = false;
	int nums[3];
	int numToRemember = 0;
	//(*log) << "s = \"" << s << "\"" << std::endl;

	/* push the individual words (numbers, hyphens, or back slashes into a vector */
	std::vector<std::string> cmds;
	std::istringstream linestream(s);
	int ch;
	std::string word = "";
	do
		{
		word = "";
		linestream >> word;
		if (word != "")
			{
			cmds.push_back( word );
			}
		} while( (ch=linestream.get()) != EOF );
		
	/* loop over the vector of individual words and correctly interpret things */
	int i = 0;
	for (std::vector<std::string>::iterator p=cmds.begin(); p != cmds.end(); p++)
		{
		//(*log) << "\"" << (*p) << "\"" << std::endl;
		if ( isdigit((*p)[0]) )
			{
			
			/* we can potentially complete a sentence */
			std::string prevWord = "";
			if (i > 0)
				prevWord = cmds[i-1];
			std::string nextWord = "";
			if (i != cmds.size() - 1)
				nextWord = cmds[i+1];
			
			int x;
			std::istringstream buf(cmds[i]);
			buf >> x;

			if (prevWord == "" || isNumber(prevWord) == true)
				{
				nums[0] = x;
				numToRemember = 1;
				}
			else if (prevWord == "-")
				{
				nums[1] = x;
				numToRemember = 2;
				}
			else if (prevWord == "\\")
				{
				nums[2] = x;
				numToRemember = 3;
				}
			else
				{
				std::cerr << "ERROR: Problem interpreting string" << std::endl;
				exit(1);
				}
			
			if ( (prevWord == "" || isNumber(prevWord) == true) && (nextWord == "" || isNumber(nextWord) == true) )
				{
				v[nums[0]-1] = true;
				numToRemember = 0;
				}
			else if ( prevWord == "-" && (nextWord == "" || isNumber(nextWord) == true) )
				{
				for (int i=nums[0]-1; i<nums[1]; i++)
					v[i] = true;
				numToRemember = 0;
				}
			else if ( prevWord == "\\" )
				{
				for (int i=nums[0]-1, k=nums[2]; i<nums[1]; i++, k++)
					if ( k % nums[2] == 0 )
						v[i] = true;
				numToRemember = 0;
				}
			}
		i++;
		}
}

bool Alignment::isNumber(std::string s) {

	if (s == "")
		return false;

	bool isnum = true;
	for (size_t i=0; i<s.size(); i++)
		if ( !isdigit(s[i]) )
			isnum = false;
	return isnum;
}

/*-------------------------------------------------------------------
|
|   NucID: 
|
|   Take a character, nuc, and return an integer:
|
|       nuc        returns
|        A            1 
|        C            2     
|        G            4      
|        T U          8     
|        R            5      
|        Y           10       
|        M            3      
|        K           12   
|        S            6     
|        W            9      
|        H           11      
|        B           14     
|        V            7      
|        D           13  
|        N - ?       15       
|
-------------------------------------------------------------------*/
int Alignment::nucID(char nuc) {

	char		n;
	
	if (nuc == 'U' || nuc == 'u')
		n = 'T';
	else
		n = nuc;

	if (n == 'A' || n == 'a')
		{
		return 1;
		}
	else if (n == 'C' || n == 'c')
		{
		return 2;
		}
	else if (n == 'G' || n == 'g')
		{
		return 4;
		}
	else if (n == 'T' || n == 't')
		{
		return 8;
		}
	else if (n == 'R' || n == 'r')
		{
		return 5;
		}
	else if (n == 'Y' || n == 'y')
		{
		return 10;
		}
	else if (n == 'M' || n == 'm')
		{
		return 3;
		}
	else if (n == 'K' || n == 'k')
		{
		return 12;
		}
	else if (n == 'S' || n == 's')
		{
		return 6;
		}
	else if (n == 'W' || n == 'w')
		{
		return 9;
		}
	else if (n == 'H' || n == 'h')
		{
		return 11;
		}
	else if (n == 'B' || n == 'b')
		{
		return 14;
		}
	else if (n == 'V' || n == 'v')
		{
		return 7;
		}
	else if (n == 'D' || n == 'd')
		{
		return 13;
		}
	else if (n == 'N' || n == 'n')
		{
		return 15;
		}
	else if (n == '-')
		{
		return 15;
		}
	else if (n == '?')
		{
		return 15;
		}
	else
		return -1;
}

void Alignment::print(void) {

	int** x = matrix;
		
	std::cout << "                        ";
	for (size_t i=0; i<numTaxa; i++)
		std::cout << std::setw(3) << i;
	std::cout << '\n';
	std::cout << "------------------------";
	for (size_t i=0; i<numTaxa; i++)
		std::cout << "---";
	std::cout << '\n';
	for (size_t j=0; j<numSitePatterns; j++)
		{
		std::cout << std::setw(4) << j+1 << " -- ";
		for (size_t i=0; i<numTaxa; i++)
			{
			std::cout << std::setw(3) << x[i][j];
			}
		std::cout << '\n';
		}
		
	int sum = 0;
}


