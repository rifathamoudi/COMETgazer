
/***************************************************************************
*
* COMETgazer Software for methylation analysis
* 
* By: Emanuele Libertini
*
*
* blocks.cpp: C++ program assigning CpGs to COMETs
*
* Input : CpG-specific OM distributions across chromosomes
*
*         Chr No | Begin | End | methylation value | OM | rounded OM value
*
* Output :  Chr No | Begin | End | methylation value | OM | rounded OM value | block
*
***************************************************************************/



#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>

using namespace std;

int main(int argc, char *argv[])
{
	int NUM_lines = 0;

    	std:string line;
   	std::ifstream myfile("chr.df.txt");
    
	while (std::getline(myfile, line))
	{
	    ++NUM_lines;
	}

  	// Main variables
  	int index;
  	int chrID;
  	int start;
  	int end;
  	float score;
  	float Delt;
  	float g ;
                                                                                                                                                                                      ////////
  	// define and open output file 

  	ofstream outputFile("output.chr.df.txt");

  	if(outputFile.is_open())
  	{
		outputFile << "chr" << "\t" << "start" << "\t" << "stop" << "\t "<< "score" << "\t" << "Delt" << "\t"<< "g" << "\t" << "block" << endl;
	  
		ifstream inputfile1("chr.df.txt");
	
		if (!inputfile1)
		{
	      		cout<< "Unable to open file\n";

		}
		else 
		{   
  			for(index = 0;index < NUM_lines;index++) 
   			{     
     				int block=1;

       				//Read the file and process blocks

	  			while (inputfile1 >> chrID >> start >> end >> score >> Delt >> g) 
				{
	      				if (g == 0)
					{
	    					outputFile << chrID << "\t" << start << "\t" << end << "\t" << score << "\t" << Delt << "\t" << g << "\t" << block << "\n";
 				        	cout << chrID << "\t" << start << "\t" << end << "\t" << score << "\t" << Delt << "\t" << g << "\t" << block << "\n";
					}
	      				else 
					{ 
						block++;
                				cout << "g is " << g << " block is " << block << endl;
	    					outputFile << chrID << "\t" << start << "\t" << end << "\t" << score << "\t" << Delt << "\t" << g << "\t" << block << "\n";
            					cout << chrID << "\t" << start << "\t" << end << "\t" << score << "\t" << Delt << "\t" << g << "\t" << block << "\n";
	      				}
	  			}
   			}

       			//Close input file
        		inputfile1.close();
		}                                                                    	   

     		//Close output file
        	outputFile.close();
	}
     	else cout << "Unable to open file";
    	return(0);
}
