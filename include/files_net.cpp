
#include <fstream.h> // file manipulation
#include <string.h> // string manipulation
#include <stdlib.h> // C++ functions...
#include <math.h> // math operations
#include <iomanip.h> // classes, struct
#include <ctype.h> // toupper (GenBank sequences)
#include <stdio.h> // C library: fprintf..
#include <iostream.h> // cout..

#define MAX_DIR 61
#define LINE_LEN 81
#define MAX 10000000 // max length for nucleic acid sequences is 1 Mb
#define MAX_PROTEIN 1000000 // max length for proteins is 50.000 amino acids

//===================================================================
//
//           Allocate names to Output files
//
//===================================================================


void Allocate_names (char *OutputFile, char *RandomFile,
							char *OutputMotifs, char *OutputScores,char *GraphicalDisplay,
                     char *SimpScores, char *dir)  // Aug 99
{

//int i;
// note: GraphicalDisplay is not currently in use

// Form output file names
strcpy(OutputFile,dir);
strcat(OutputFile,"S1");
strcpy(OutputMotifs,dir);
strcat(OutputMotifs,"S2");
strcpy(OutputScores,dir);
strcat(OutputScores,"S3");
strcpy(RandomFile,dir);
strcat(RandomFile,"S4");
strcpy(GraphicalDisplay,dir);
strcat(GraphicalDisplay,".PS");
strcpy(SimpScores,dir);
strcat(SimpScores,"SIMPDATA");


// write name of the directory as part of the name of the file, Aug 99
//for (i=0; i<MAX_DIR; i++)
//	{
//   OutputFile[i]=RandomFile[i]=OutputMotifs[i]=
//   OutputScores[i]=GraphicalDisplay[i]=SimpScores[i]=dir[i];
//   if (dir[i]=='\0' || dir[i]=='\n')
//   	break;
//   }

// file S1
//OutputFile[i+1]='S'; OutputFile[i+2]='1'; OutputFile[i+3]='\0';

// file S2
//OutputMotifs[i+1]='S';OutputMotifs[i+2]='2';OutputMotifs[i+3]='\0';

// file S3
//OutputScores[i+1]='S';OutputScores[i+2]='3';OutputScores[i+3]='\0';

// file S4
//RandomFile[i+1]='S'; RandomFile[i+2]='4'; RandomFile[i+3]='\0';

// file PS (not currently in use)

//GraphicalDisplay[i+1]='.'; GraphicalDisplay[i+2]='P'; GraphicalDisplay[i+3]='S'; GraphicalDisplay[i+4]='\0';

// SimpScore (SIMPDATA) Aug 99

//SimpScores[i+1]='S';SimpScores[i+2]='I';SimpScores[i+3]='M';SimpScores[i+4]='P';
//SimpScores[i+5]='D';SimpScores[i+6]='A';SimpScores[i+7]='T';SimpScores[i+8]='A';SimpScores[i+9]='\0';


}

//===================================================================
//
//           Get sequence from file
//
//===================================================================



// ===== Function Get_sequence: get id, seq, base comp and length from file,
// transform into DNA or Purine/pyrimidine if requested

void Get_sequence (struct sequence *seq1, struct param_process *userp,
						 char InputFile[LINE_LEN])
{
char input_line[LINE_LEN], c;
int i,j;
ifstream in_file;
in_file.open(InputFile);
if (!in_file)
	{
	cerr << "sorry, can´t open input file\n";
	exit(1);
	}

// GenBank file

if ((userp->FileType==1)&&
	 ((strncmp(userp->seq_nature,"y",1)==0)||(strncmp(userp->seq_nature,"n",1)==0)))
	{
	while (in_file.getline (input_line, LINE_LEN))
		{
		if (strncmp (input_line,"LOCUS"+0,5)==0)
			{
         seq1->seq_id1 = new char[sizeof(char)*(LINE_LEN)];
			if (!seq1->seq_id1)
				{
				cerr << "Insufficient memory for seq_id1";
				exit (1);
            }
         strcpy(seq1->seq_id1, input_line);
         }
		if (strncmp (input_line,"DEFINITION"+0,10)==0)
			{
         seq1->seq_id2 = new char[sizeof(char)*(LINE_LEN)];
			if (!seq1->seq_id2)
				{
				cerr << "Insufficient memory for seq_id2";
				exit(1);
            }
			strcpy(seq1->seq_id2, input_line);
			}
		if (strncmp (input_line, "ACCESSION"+0,9)==0)
        	{
         seq1->seq_id3 = new char[sizeof(char)*(LINE_LEN)];
			if (!seq1->seq_id3)
				{
				cerr << "Insufficient memory for seq_id3";
				exit(1);
            }
			strcpy(seq1->seq_id3, input_line);
			}
		if (strncmp (input_line, "NID"+0,3)==0)
			{
         seq1->seq_id4=new char[sizeof(char)*(LINE_LEN)];
			if (!seq1->seq_id4)
				{
				cerr << "Insufficient memory for seq_id4";
				exit(1);
            }
			strcpy(seq1->seq_id4, input_line);
			}
		if (strncmp (input_line,"ORIGIN"+0,6) == 0)
			{
			seq1->seq = new char[sizeof(char)*(MAX)];
			if (!seq1->seq)
         	{
				cerr << "Not enough memory for sequence";
				exit(1);
            }
			in_file.get(c);
			seq1->length=0;
			while (c!='/')
        		{
            c=toupper(c);
            if ((c=='A')||(c=='G'))
            	{
               seq1->length++;
               if (strncmp(userp->seq_nature,"y",1)==0)
            		seq1->seq[seq1->length-1]='R'; // conversion to purines
               else seq1->seq[seq1->length-1]=c;
               }
         	else if ((c=='C')||(c=='T')||(c=='U'))
            	{
               seq1->length++;
               if (c=='U')                       // conversion U to T
               	seq1->seq[seq1->length-1]='T';
               else
               	seq1->seq[seq1->length-1]=c;
               if (strncmp(userp->seq_nature,"y",1)==0)              // conversion to pyrimidines
               	seq1->seq[seq1->length-1]='Y';
               }
				in_file.get(c);
				}
			}
      }
	}

// EMBL file

else if ((userp->FileType==2)&&
         ((strncmp(userp->seq_nature,"y",1)==0)||(strncmp(userp->seq_nature,"n",1)==0)))
	{
	while (in_file.getline (input_line, LINE_LEN))
   	{
		if (strncmp (input_line,"ID"+0,2) == 0)
			{
			seq1->seq_id1= new char[sizeof(char)*(LINE_LEN)];
			if (!seq1->seq_id1)
         	{
				cerr << "Unsufficient memory for seq_id1";
            exit(1);
            }
	    	strncpy(seq1->seq_id1, input_line, LINE_LEN);
			}
		if (strncmp (input_line,"AC"+0,2) == 0)
			{
			seq1->seq_id3 = new char[sizeof(char)*(LINE_LEN)];
			if (!seq1->seq_id3)
         	{
				cerr << "Unsufficient memory for seq_id3";
            exit(1);
            }
			strncpy(seq1->seq_id3, input_line, LINE_LEN);
			}
		if (strncmp (input_line,"NI"+0,2) == 0)
			{
			seq1->seq_id4 = new char[sizeof(char)*(LINE_LEN)];
			if (!seq1->seq_id4)
         	{
				cerr << "Unsufficient memory for seq_id4";
            exit(1);
            }
			strncpy(seq1->seq_id4, input_line, LINE_LEN);
			}
		if (strncmp (input_line,"DE"+0,2) == 0)
			{
			seq1->seq_id2 = new char[sizeof(char)*(LINE_LEN)];
			if (!seq1->seq_id2)
         	{
				cerr << "Unsufficient memory for seq_id2";
         	exit(1);
            }
			strncpy(seq1->seq_id2, input_line, LINE_LEN);
			}
		if (strncmp (input_line,"SQ"+0,2) == 0)
			{
			seq1->seq = new char[sizeof(char)*(MAX)];
			if (!seq1->seq)
         	{
				cerr << "Not enough memory for sequence";
            exit(1);
            }
			in_file.get(c);
			seq1->length=0;
			while (c!='/')
        		{
            c=toupper(c);
            if ((c=='A')||(c=='G'))
            	{
               seq1->length++;
               if (strncmp(userp->seq_nature,"y",1)==0)
            		seq1->seq[seq1->length-1]='R'; // conversion to purines
               else seq1->seq[seq1->length-1]=c;
               }
         	else if ((c=='C')||(c=='T')||(c=='U'))
            	{
               seq1->length++;
               if (c=='U')                       // conversion U to T
               	seq1->seq[seq1->length-1]='T';
               else seq1->seq[seq1->length-1]=c;
               if (strncmp(userp->seq_nature,"y",1)==0)              // conversion to pyrimidines
               	seq1->seq[seq1->length-1]='Y';
               }
				in_file.get(c);
				}
			}
		}
	}

// "Only sequence" file

else if (userp->FileType==3)
	{
	seq1->length=0;
   	seq1->seq_id1=seq1->seq_id2=seq1->seq_id3=seq1->seq_id4="..";
   	if (strncmp(userp->seq_nature,"p",1)==0)
			seq1->seq= new char[sizeof(char)*(MAX_PROTEIN)];
	else
      	seq1->seq= new char[sizeof(char)*(MAX)];
   	if (!seq1->seq)
   		{
		cerr << "Not enough memory for sequence";
		exit(1);
      		}
	/*in_file.getline(input_line,LINE_LEN);
	if (input_line[0]=='>')
		{
		in_file.getline(input_line,LINE_LEN); // for sequences with a header line
		}
	// first line
	for (i=0; i<LINE_LEN; i++) 
		{
		if (input_line[i]=='\0')
			{
			break;
			}
		c=input_line[i];
		c=toupper(c);
      		if (((strncmp(userp->seq_nature,"n",1)==0)||(strncmp(userp->seq_nature,"y",1)==0)) &&
      	 	(c=='U'))
      		c='T';                                     // Conversion of U to T
      		if ((strncmp(userp->seq_nature,"y",1)==0)&&((c=='A')||(c=='G')))  // Conversion to Pur/Pyr
		c='R';
      		if ((strncmp(userp->seq_nature,"y",1)==0)&&((c=='C')||(c=='T')))
      		c='Y';
      		for (i=0; i<userp->NbrOfElements; i++)
      			{
         		if (c==seq1->elements[i])
				{
         			seq1->seq[seq1->length]=c;
					seq1->length++;
            			}
         		}							
		}	
	while(in_file.getline(input_line,LINE_LEN))
		{
		for (i=0; i<LINE_LEN; i++) 
			{
		      	if (input_line[i]=='\0')
				{
				break;
				}
			c=input_line[i];
			c=toupper(c);
      			if (((strncmp(userp->seq_nature,"n",1)==0)||(strncmp(userp->seq_nature,"y",1)==0)) &&
      	 		(c=='U'))
      			c='T';                                     // Conversion of U to T
      			if ((strncmp(userp->seq_nature,"y",1)==0)&&((c=='A')||(c=='G')))  // Conversion to Pur/Pyr
			c='R';
      			if ((strncmp(userp->seq_nature,"y",1)==0)&&((c=='C')||(c=='T')))
      			c='Y';
      			for (i=0; i<userp->NbrOfElements; i++)
      				{
         			if (c==seq1->elements[i])
					{
         				seq1->seq[seq1->length]=c;
					seq1->length++;
            				}
         			}							
			}
		}*/
	in_file.get(c);
	if (c=='>') 
		{
		while((c!='\0')&&(c!='\n'))
			{
			in_file.get(c);
			}
		}
	else
		{
      		c=toupper(c);
      		if (((strncmp(userp->seq_nature,"n",1)==0)||(strncmp(userp->seq_nature,"y",1)==0)) &&
      	 	(c=='U'))
      		c='T';                                     // Conversion of U to T
      		if ((strncmp(userp->seq_nature,"y",1)==0)&&((c=='A')||(c=='G')))  // Conversion to Pur/Pyr
			c='R';
      		if ((strncmp(userp->seq_nature,"y",1)==0)&&((c=='C')||(c=='T')))
      		c='Y';
      		for (i=0; i<userp->NbrOfElements; i++)
      			{
         		if (c==seq1->elements[i])
				{
         			seq1->seq[seq1->length]=c;
				seq1->length++;
            			}
         		}		
						
		}	
	while (in_file.get(c))
   		{
      		c=toupper(c);
      		if (((strncmp(userp->seq_nature,"n",1)==0)||(strncmp(userp->seq_nature,"y",1)==0)) &&
      	 	(c=='U'))
      		c='T';                                     // Conversion of U to T
      		if ((strncmp(userp->seq_nature,"y",1)==0)&&((c=='A')||(c=='G')))  // Conversion to Pur/Pyr
			c='R';
      		if ((strncmp(userp->seq_nature,"y",1)==0)&&((c=='C')||(c=='T')))
      		c='Y';
      		for (i=0; i<userp->NbrOfElements; i++)
      			{
         		if (c==seq1->elements[i])
				{
         			seq1->seq[seq1->length]=c;
				seq1->length++;
            			}
         		}
		}
      }
   

// SwissProt file

else if ((userp->FileType==2)&&(strncmp(userp->seq_nature,"p",1)==0))
	{
   seq1->seq_id4="..";
   while (in_file.getline (input_line, LINE_LEN))
   	{
		if (strncmp (input_line,"ID"+0,2) == 0)
			{
			seq1->seq_id1= new char[sizeof(char)*(LINE_LEN)];
			if (!seq1->seq_id1)
         	{
				cerr << "Unsufficient memory for seq_id1";
            exit(1);
            }
	    	strncpy(seq1->seq_id1, input_line, LINE_LEN);
			}
		if (strncmp (input_line,"AC"+0,2) == 0)
			{
			seq1->seq_id3 = new char[sizeof(char)*(LINE_LEN)];
			if (!seq1->seq_id3)
         	{
				cerr << "Unsufficient memory for seq_id3";
            exit(1);
            }
			strncpy(seq1->seq_id3, input_line, LINE_LEN);
			}
		if (strncmp (input_line,"DE"+0,2) == 0)
			{
			seq1->seq_id2 = new char[sizeof(char)*(LINE_LEN)];
			if (!seq1->seq_id2)
         	{
				cerr << "Unsufficient memory for seq_id2";
         	exit(1);
            }
			strncpy(seq1->seq_id2, input_line, LINE_LEN);
			}
		if (strncmp (input_line,"SQ"+0,2) == 0)
			{
			seq1->seq = new char[sizeof(char)*(MAX_PROTEIN)];
			if (!seq1->seq)
         	{
				cerr << "Not enough memory for sequence";
            exit(1);
            }
			in_file.get(c);
			seq1->length=0;
			while (c!='/')
        		{
            c=toupper(c);
				for (i=0; i<userp->NbrOfElements; i++)
      			{
         		if (c==seq1->elements[i])
						{
         			seq1->seq[seq1->length]=c;
						seq1->length++;
            		}
         		}
            in_file.get(c);
            }
         }
      }
   }

// GenBank for protein (only one translation in the file, assumes toupper characters)
else if ((userp->FileType==1)&&(strncmp(userp->seq_nature,"p",1)==0))
   {
	while (in_file.getline (input_line, LINE_LEN))
		{
		if (strncmp (input_line,"LOCUS",5)==0)
			{
         seq1->seq_id1 = new char[sizeof(char)*(LINE_LEN)];
			if (!seq1->seq_id1)
				{
				cerr << "Insufficient memory for seq_id1";
				exit (1);
            }
         strcpy(seq1->seq_id1, input_line);
         }
		if (strncmp (input_line,"DEFINITION",10)==0)
			{
         seq1->seq_id2 = new char[sizeof(char)*(LINE_LEN)];
			if (!seq1->seq_id2)
				{
				cerr << "Insufficient memory for seq_id2";
				exit(1);
            }
			strcpy(seq1->seq_id2, input_line);
			}
		if (strncmp (input_line, "ACCESSION",9)==0)
        	{
         seq1->seq_id3 = new char[sizeof(char)*(LINE_LEN)];
			if (!seq1->seq_id3)
				{
				cerr << "Insufficient memory for seq_id3";
				exit(1);
            }
			strcpy(seq1->seq_id3, input_line);
			}
		if (strncmp (input_line, "NID",3)==0)
			{
         seq1->seq_id4=new char[sizeof(char)*(LINE_LEN)];
			if (!seq1->seq_id4)
				{
				cerr << "Insufficient memory for seq_id4";
				exit(1);
            }
			strcpy(seq1->seq_id4, input_line);
			}
		if (strncmp (input_line,"                     /translation",33)==0)
			{
         seq1->seq = new char[sizeof(char)*(MAX_PROTEIN)];
			if (!seq1->seq)
         	{
				cerr << "Not enough memory for sequence";
            exit(1);
            }
         seq1->length=0;
         for (i=35; i<LINE_LEN; i++)
         	{
            c=input_line[i];
            c=toupper(c);
				for (j=0; j<userp->NbrOfElements; j++)
      			{
         		if (c==seq1->elements[j])
						{
         			seq1->seq[seq1->length]=c;
						seq1->length++;
            		}
         		}
            }
         for (i=0; i<(MAX_PROTEIN)*2; i++)
         	{
            in_file.get(c);
            c=toupper(c);
				for (j=0; j<userp->NbrOfElements; j++)
      			{
         		if (c==seq1->elements[j])
						{
         			seq1->seq[seq1->length]=c;
						seq1->length++;
            		}
         		}
            if ((c=='"')||(seq1->length==MAX_PROTEIN))
            	break;
            }
			}
      }
	}
in_file.close();
}

//===================================================================
//
//           Output routines
//
//===================================================================

void Output_test_sequence (struct sequence *seq1, struct param_process *userp,
				 char OutputFile[], char InputFile[])
{
int i;
FILE *File;
File = fopen (OutputFile, "w");
if (File == NULL)
  {
    printf("*** ERROR. Unable to open file [%s]\n",OutputFile);
    exit(1);
  }
if (userp->FileType!=3)
	{
	fprintf(File, "\n %s", seq1->seq_id1);
	fprintf(File, "\n %s", seq1->seq_id2);
	fprintf(File, "\n %s", seq1->seq_id3);
	fprintf(File, "\n %s\n", seq1->seq_id4);
	}
fprintf(File,"\n *** RESULTS ***\n");
fprintf(File, "\n *** Sequence characteristics ***\n\n");
for (i=0; i<(userp->NbrOfElements); i++)
	{
	fprintf(File, " %c:%3d ", seq1->elements[i],seq1->count[i]);
	if (i==9)
   	fprintf(File, "\n");
   }
fprintf(File, "\n\n Length: %3d \n", seq1->length);
fprintf(File,"\n *** Calculation of simplicity ***\n");
fprintf(File, "\n Score mono-elements: %3d  Score di-elements: %3d",
        userp->score_mono, userp->score_di);
fprintf(File, "\n Score tri-elements: %3d  Score tetra-elements: %3d Size of window: %3d",
        userp->score_tri, userp->score_tetra, userp->window);
fprintf(File, "\n Simplicity was calculated from element %3d to element %3d\n",
        userp->window+1, seq1->length-userp->window-3);
fprintf(File,"\n *** Randomization of sequences ***\n");
fprintf(File,"\n Number of generated random sequences: %3d",userp->num_random);
if (userp->random_type==1)
	fprintf(File,"\n Randomization according to element frequency.\n");
else if (userp->random_type==3)
	fprintf(File,"\n Randomization according to di-element frequency.\n");
else if (userp->random_type==4)
	fprintf(File,"\n Randomization independently in three different frames.\n");
else if (userp->random_type==2)
	fprintf(File,"\n Randomization by random positioning of elements.\n");
fclose(File);
}


