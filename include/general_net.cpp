
#include <string.h> // string manipulation
#include <stdlib.h> // C++ functions...
#include <math.h> // math operations
#include <iomanip.h> // classes, struct
#include <stdio.h> // C library: fprintf..
#include <iostream.h> // cout..

#define MAX_DIR 61
#define LINE_LEN 81

#ifdef MAX
#undef MAX
#endif

#define MAX 1000000
#define TRUE 1
#define FALSE 0

//===================================================================
//
//            Get command arguments
//
//===================================================================

void Get_command_arguments(char *strptrs[], int num, char InputFile[],struct param_process *userp, char directory[])
{
	/* Initialise the parameters */
	InputFile[0] = '\0';

	/* Get the sequence file name */
	strcpy(InputFile,strptrs[1]);

	/* Get the parameters in userp*/
	strcpy(userp->seq_nature, strptrs[2]);
	userp->FileType=atoi(strptrs[3]);
	userp->score_mono=atoi(strptrs[4]);
	userp->score_di=atoi(strptrs[5]);
	userp->score_tri=atoi(strptrs[6]);
	userp->score_tetra=atoi(strptrs[7]);
	userp->window=atoi(strptrs[8]);
	userp->num_random=atoi(strptrs[9]);
	userp->random_type=atoi(strptrs[10]);
	userp->stringency=atof(strptrs[11]);
	strcpy(userp->graphical_display, strptrs[12]);
	strcpy(directory,strptrs[13]); // Aug. 99

	/* Show name of sequence file */
	if (InputFile[0] == '\0')
	{
		printf("*** ERROR. No sequence file entered\n");
		 exit(-1);
	}
	else
	{
		printf("\nParameters are: %s %s %d %d %d %d %d %d %d %d %f %s %s",InputFile,
			userp->seq_nature, userp->FileType, userp->score_mono, userp->score_di,
			userp->score_tri, userp->score_tetra, userp->window, userp->num_random,
			userp->random_type, userp->stringency, userp->graphical_display, directory);
		printf("\nNumber of parameters entered: %d", num);
	}
}

//===================================================================
//
//            Sequence composition
//
//===================================================================

// ===== Function Assign_elements: define elements of sequence =====
// This will be used when getting the elements from the sequence in the file

void Assign_elements (struct sequence *seq1, struct param_process *userp)
{
	int i;
	char *protein_elements={"GVIYWNSCPRAKFDEQTMHL"};
	char *nucleic_elements={"GATC"};
	char *purine_pyrimidine_elements={"RY"};

	if (strncmp(userp->seq_nature,"n",1)==0)    // assign elements
	{
		userp->NbrOfElements=4;
		seq1->elements=new char[sizeof(char)*(userp->NbrOfElements)];

		if (!seq1->elements)
		{
			cerr << '\n' << "Not enough memory";
			exit(1);
		}

		strncpy(seq1->elements, nucleic_elements, (userp->NbrOfElements));
	}
	else if (strncmp(userp->seq_nature,"y",1)==0)
	{
		userp->NbrOfElements=2;
		seq1->elements=new char[sizeof(char)*(userp->NbrOfElements)];

		if (!seq1->elements)
   		{
      			cerr << '\n' << "Not enough memory";
      			exit(1);
      		}

		strncpy(seq1->elements, purine_pyrimidine_elements, (userp->NbrOfElements));
	}
	else if (strncmp(userp->seq_nature,"p",1)==0)
	{
		userp->NbrOfElements=20;
		seq1->elements=new char[sizeof(char)*(userp->NbrOfElements)];

		if (!seq1->elements)
   		{	
      			cerr << '\n' << "Not enough memory";
      			exit(1);
      		}

   		strncpy(seq1->elements, protein_elements, (userp->NbrOfElements));
	}

	seq1->count=new int[sizeof(int)*(userp->NbrOfElements)]; // assign count of elements


	for (i=0; i<(userp->NbrOfElements); i++) seq1->count[i]=0;

}

// === Function Three_frames_frequency: calculates number of elements in each
// different frame ====

void Three_frames_frequency (struct sequence *seq1, int n)
{
long i;
int k,z;
long int all_frames;
float l;
l=int(float(seq1->length)/3); // divide the length in three parts
all_frames=(long int)(3.0*l); // part of the sequence that contains the three frames


seq1->count_three_frames=new int*[3];  // allocate count in three different frames
for (k=0; k<3; k++)
	{
   seq1->count_three_frames[k]=new int[n];
   }

for (k=0; k<3; k++)    // assign 0 to the count in three different frames
	{
   for (z=0; z<n; z++)
   	seq1->count_three_frames[k][z]=0;
   }

for (i=0; i<(l-1); i++)  // count elements in the part the contains
	{                              // the three frames
   for (k=0; k<3; k++)
   	{
      for (z=0; z<n; z++)
      	{
         if (seq1->seq[(i*3)+k]==seq1->elements[z])
         	seq1->count_three_frames[k][z]++;
         }
      }
   }

for (i=0; i<(seq1->length-all_frames); i++)  // end up counting elements
	{
   for (z=0; z<n; z++)
   	{
      if (seq1->seq[all_frames+i]==seq1->elements[z])
      	seq1->count_three_frames[i][z]++;
      }
   }

}

// === Function DiFrequency: calculates number of di-elements ===

void DiFrequency (struct sequence *seq1, int n)
{
long i;
int k,z;
seq1->di_count=new int*[n]; // allocate di-element matrix
for (z=0; z<n; z++)
	seq1->di_count[z]=new int[n];

for (z=0; z<n; z++) // assign 0 to di-element count
	{
 	for (k=0; k<n; k++)
   	seq1->di_count[z][k]=0;
	}

for (i=0; i<((seq1->length)-1); i++) // di-element count
	{
   for (k=0; k<n; k++)
		{
		for (z=0; z<n; z++)
			{
			if ((seq1->seq[i]==seq1->elements[k])&&
         	 (seq1->seq[i+1]==seq1->elements[z]))
         	{
				seq1->di_count[k][z]++;
            }
			}
		}
	}
}

// ==== Determine sequence characteristics ====
// This counts the elements in the sequence and the di-elements if requested
// It will be used to generate random sequences of the same composition and
// keep record of the number of different elements
// It also defines max, which is the last nucleotide to which a simplicity score
// will be assigned

void Determine_sequence_characteristics (struct sequence *seq1, int n,
													  struct param_process *userp)
{
// define max
userp->max=(seq1->length)-(userp->window);
if (userp->score_tetra!=0)
	userp->max=userp->max-3;
else if ((userp->score_tri!=0) && (userp->score_tetra==0))
	userp->max=userp->max-2;
else if ((userp->score_tri==0) && (userp->score_tetra==0) && (userp->score_di!=0))
	userp->max=userp->max-1;

// allocate variables
long i;
int j;

for (i=0; i<seq1->length; i++)   // count different elements
	{
   for (j=0; j<n; j++)
   	{
		if (seq1->seq[i]==seq1->elements[j])
      	seq1->count[j]++;
      }
   }

if (userp->random_type==3)  // calculate di-element frequency if required
   DiFrequency(seq1, n);
else
	seq1->di_count=NULL;
if (userp->random_type==4)
	Three_frames_frequency(seq1, n);
else
	seq1->count_three_frames=NULL;
}


//===================================================================
//
//            Calculations on simplicity
//
//===================================================================


// ===== Score_frequencies =====
float* FindScoreFrequency(long *ss,int maxs,long max,int win)
{
long i;
int j;
float *freq;
freq = new float[sizeof(float)*(maxs+1)];
if (!freq)
	{
   cerr << '\n' << "Not enough memory";
   exit(1);
   }
for (j=0; j<=maxs; j++)
	freq[j]=0;

for (i=win; i<max; i++)    // find frequency of scores
	{
   for (j=0; j<=maxs; j++)
   	{
      if (ss[i]==j)
         {
      	freq[j]++;
         break;
			}
      }
   }
return(freq);
}

// ==== Display simplicity along the sequence ====
// Writes in a file the position on the sequence(X) and the simplicity score(Y)

void DisplaySimplicity(long *ss, long max,int win, char *SimpScores)
{
long i;
FILE *File;
File=fopen(SimpScores,"w");
fprintf(File, "position\tss\n");  // mod Aug 99

for (i=win; i<max; i++)
	{
   fprintf(File, "%6ld %3ld\n",i+1,ss[i]);
   }
}


// ====== Maximum score ======
int MaxSimplicity(long *ss,long max,int win)
{
int maxs=0;
long i;
for (i=win; i<max; i++)
	{
   if (ss[i]>maxs)
   	maxs=ss[i];
   }
return (maxs);
}

// ===== Simplicity Factor ======
float SimpleFactor(long *ss,long max,int win)
{
float sf=0;
long i;
for (i=win; i<max; i++)
	{
   sf=sf+ss[i];
   }
sf=sf/(max-win);
return(sf);
}

// ====== Simplicity Scores ======
// In this version of SimplicityScores the repeat starting at the reference nucleotide
// to the left is taken as the reference repeat. There will be as many reference
// repeats as scores for different short repeats being different than 0
// The window on the right of the reference nucleotide will have a length of
// userp->window plus the number of nucleotides of the highest reference repeat
// minus 1 (there is no alteration for mononucleotide repeats)

long* SimplicityScores(char *simple_seq, struct param_process *userp)
{
int j,k;
long i;
long *ss;
ss= new long[sizeof(long)*(userp->max)];  // initialize ss
if (!ss)
	{
   cerr << '\n' << "Not enough memory";
   exit(1);
   }

for (i=userp->window; i<userp->max; i++)
	{
	ss[i]=0;
   if (userp->score_mono!=0)            // mono-elements
		{
      j=i+1;
      k=i-1;
		while ((j<i+1+userp->window)&&(k>=i-(userp->window)))
			{
      	if (simple_seq[i]==simple_seq[j])
      		{
         	ss[i]=ss[i]+userp->score_mono;
         	}
         if (simple_seq[i]==simple_seq[k])
         	{
            ss[i]=ss[i]+userp->score_mono;
         	}
         j=j+1;
         k=k-1;
         }
      }
   if (userp->score_di!=0)               // di-elements
   	{
      j=i+2;
      while (j<i+1+userp->window)
      	{
         if ((simple_seq[i]==simple_seq[j])&&(simple_seq[i+1]==simple_seq[j+1]))
				{
            ss[i]=ss[i]+userp->score_di;
            j=j+2;
            }
         else j=j+1;
			}
      k=i-2;
      while (k>=i-(userp->window))
      	{
         if ((simple_seq[i]==simple_seq[k])&&(simple_seq[i+1]==simple_seq[k+1]))
				{
            ss[i]=ss[i]+userp->score_di;
            k=k-2;
            }
         else
         	k=k-1;
         }
      }
   if (userp->score_tri!=0)                 // tri-elements
   	{
      j=i+3;
      while (j<i+1+userp->window)
      	{
         if ((simple_seq[i]==simple_seq[j])&&(simple_seq[i+1]==simple_seq[j+1])&&
             (simple_seq[i+2]==simple_seq[j+2]))
				{
            ss[i]=ss[i]+userp->score_tri;
            j=j+3;
            }
         else j=j+1;
			}
      k=i-3;
      while (k>=i-(userp->window))
      	{
         if ((simple_seq[i]==simple_seq[k])&&(simple_seq[i+1]==simple_seq[k+1])&&
             (simple_seq[i+2]==simple_seq[k+2]))
				{
            ss[i]=ss[i]+userp->score_tri;
            k=k-3;
            }
         else
         	k=k-1;
         }
      }
   if (userp->score_tetra!=0)                // tetra-elements
   	{
      j=i+4;
      while (j<i+1+userp->window)
      	{
         if ((simple_seq[i]==simple_seq[j])&&(simple_seq[i+1]==simple_seq[j+1])&&
             (simple_seq[i+2]==simple_seq[j+2])&&(simple_seq[i+3]==simple_seq[j+3]))
				{
            ss[i]=ss[i]+userp->score_tetra;
            j=j+4;
            }
         else
         	j=j+1;
			}
      k=i-4;
      while (k>=i-(userp->window))
      	{
         if ((simple_seq[i]==simple_seq[k])&&(simple_seq[i+1]==simple_seq[k+1])&&
             (simple_seq[i+2]==simple_seq[k+2])&&(simple_seq[i+3]==simple_seq[k+3]))
				{
            ss[i]=ss[i]+userp->score_tetra;
            k=k-4;
            }
         else
         	k=k-1;
         }
      }
   }

return(ss);
}


// ==== Calculations on simplicity of test sequence =====

void Simplicity_test_sequence (char *simple_seq, struct param_process *userp,
					                struct simple_scores *test, char *SimpScores)
{
// calculate simplicity scores
test->ss_test=SimplicityScores(simple_seq,userp);

// display simplicity scores along the sequence if required
if (strncmp(userp->graphical_display,"y",1)==0)
	DisplaySimplicity(test->ss_test,userp->max,userp->window, SimpScores);

// calculate overall simplicity
test->SimplicityFactor=SimpleFactor(test->ss_test,userp->max,userp->window);

// calculate maximum score
test->MaxScore=MaxSimplicity(test->ss_test,userp->max,userp->window);

// calculate frequency of scores
test->score_frequency=FindScoreFrequency(test->ss_test,test->MaxScore,userp->max,
													  userp->window);
}

// ==== Calculations on simplicity of random sequence =====

void Simplicity_random_sequence(char *simple_seq, struct param_process *userp,
                                float & sr, struct simple_scores *random,int maxs,double *possd)
{
	// calculate simplicity scores
	int i;
	long *ss_random;
	float *freq_random;

	ss_random = SimplicityScores(simple_seq,userp);

	// calculate overall simplicity
	sr=SimpleFactor(ss_random,userp->max,userp->window);

	// calculate frequency of scores
	freq_random=FindScoreFrequency(ss_random, maxs,userp->max,userp->window);

	// add to values in struct simple_scores random (score_frequency) and simplicity
	// factor of random sequence (sr)
	for (i=0; i<(maxs+1); i++)
	{
   		random->score_frequency[i]+=freq_random[i];
   	}

	random->SimplicityFactor+=sr;

	// ADDED SJG 2005/08/10

	double z_ss=0,z_ss2=0;

	// printf("Random sequence (length=%i):\n",userp->max);

	for(i=0;i<userp->max;i++)
    	{
		z_ss+=ss_random[i];
        	z_ss2+=ss_random[i]*ss_random[i];

		// printf("\t%i\t%ld\n",i,ss_random[i]);
    	}

    	*possd=sqrt((userp->max*z_ss2-(z_ss*z_ss))/(userp->max*userp->max));

	// printf("ss=%lf ss2=%lf SD=%lf\n",z_ss,z_ss2,*possd);

    	// END OF ADDITION

	// delete variables
	delete ss_random, freq_random;
}

//===================================================================
//
//            Output Random sequences
//
//===================================================================


void WriteOut_Random(float sr, char *ele, int num_random, int *num,
							char RandomFile[], char InputFile[], int n)
{
	int i;
	long total=0;
	for (i=0; i<n; i++) total+=num[i];

	FILE *File;

	if (num_random!=0)
		File = fopen (RandomFile, "a");
	else
   	{
		File = fopen (RandomFile, "w");
   	fprintf(File,"\n Input file: %s\n", InputFile);
	}

	fprintf(File,"\n Random sequence %3d :  simplicity : %f \n ", num_random+1, sr);

	for (i=0; i<n; i++) fprintf(File,"%3c %3d", ele[i], num[i]);

	fprintf(File,"  Total: %3d ",total);
	fclose(File);
}

// SJG Added 'ave_sd'

void WriteOut_Variance_Confidence(float var,float ave_sd,float random_mean, int n, float test_mean,
										    float significance, char OutputFile[])

{
	// calculate confidence limits
	float z1=1.96; // for 95% confidence
	float z2=2.58; // for 99% confidence
	float confidence_limit, confidence_limit_high;
	float sample_size=float(n);
	float standard_error=sqrt(var)/sqrt(sample_size);
	confidence_limit=z1*(standard_error);
	confidence_limit_high=z2*(standard_error);

	// output
	FILE *File;
	File = fopen (OutputFile, "a");
	fprintf(File,"\n *** Simplicity ***\n");
	fprintf(File,"\n Simplicity test sequence %.4f", test_mean);
	fprintf(File,"\n Average simplicity random sequences: %.4f", random_mean);
	fprintf(File,"\n Average of each sequences positional simplicity scores: %.4f", ave_sd);
	fprintf(File,"\n Ratio (test/random): %.4f", significance);
	fprintf(File, "\n\n *** Confidence limits ***\n");
	fprintf(File,"\n Variance simplicity random sequences: %.4f", var);
	fprintf(File,"\n Standard Error: %.4f", standard_error);
	fprintf(File,"\n Confidence limit 0.95: (%.4f - %.4f)", random_mean-confidence_limit,
		  random_mean+confidence_limit);
	fprintf(File,"\n Confidence limit 0.99: (%.4f - %.4f)", random_mean-confidence_limit_high,
		  random_mean+confidence_limit_high);
		  
	fprintf(File,"\n\n *** Significance of simplicity ***\n");

	if(((test_mean)/(random_mean+confidence_limit_high))>1)
	{
		fprintf(File,"\nSimplicity in sequence is significant(confidence 0.99)");
	}
	else if (((test_mean)/(random_mean+confidence_limit))>1)
	{
		fprintf(File,"\nSimplicity in sequence is significant(confidence 0.95)");
	}
	else fprintf(File,"\nSimplicity in sequence is not significant");

	fclose(File);
}

//===================================================================
//
//            Random sequences
//
//===================================================================

// *** Function GetRand3: makes randomized sequence using dinucleotide probability ***

void GetRand3(char *simple_rseq, struct sequence *seq1, int *num, int n)
{
long i;
int j,k;
int r;
float **di_ratio;
float *cumulative_di_ratio;
di_ratio = new float*[sizeof(float)*n];
for (i=0; i<n; i++)
	{
	di_ratio[i] = new float[sizeof(float)*n];
	cumulative_di_ratio = new float[sizeof(float)*n];
	}

// define range of random numbers to use
int range=1000;

// define dinucleotide ratio
for (i=0; i<n; i++)
	{
   if (seq1->elements[i]==seq1->seq[seq1->length-1]) // discount last nucleotide
   	seq1->count[i]--;
	for (j=0; j<n; j++)
   	{
	   di_ratio[i][j]=(float(seq1->di_count[i][j])/float(seq1->count[i]))*float(range);
	   }
   }

for (i=1; i<seq1->length; i++)
	{
	r=rand() %range;
   for (j=0; j<n; j++)
   	{
   	if (simple_rseq[i-1]==seq1->elements[j])
			{
         cumulative_di_ratio[j]=di_ratio[j][0];
         for (k=0; k<n; k++)
         	{
            if (r<cumulative_di_ratio[j])
					{
					simple_rseq[i]=seq1->elements[k];
					num[k]++;
               break;
					}
				else
				cumulative_di_ratio[j]+=di_ratio[j][k+1];
				}
         }
    	}
   }
}

// *** Function GetRand1: makes randomized sequence of the same
// length than the query sequence and with the same probability of bases ***

char* GetRand1(struct sequence *seq1, struct param_process *userp, int *num)
{
long i,l;
int j;
int r;
float *ratio;
float cumulative_ratio;
ratio = new float[sizeof(float)*(userp->NbrOfElements)];
char *simple_rseq;
simple_rseq = new char[sizeof(char)*(seq1->length)];
if (!simple_rseq)
	{
	cerr << '\n' << "Not enough memory";
	exit(1);
	}

// define range of random numbers to use
int range=1000;

// define length
if (userp->random_type==3)
	l=1;
else
	l=seq1->length;

// define sequence element ratio

for (i=0; i<(userp->NbrOfElements); i++)
	{
	ratio[i]=(float(seq1->count[i])/float(seq1->length))*float(range);
	}

// assign elements to random sequence
for (i=0; i<l; i++)
	{
	r=rand() %range;
   cumulative_ratio=ratio[0];
   for (j=0; j<(userp->NbrOfElements); j++)
   	{
      if (r<cumulative_ratio)
      	{
         simple_rseq[i]=seq1->elements[j];
         num[j]++;
         break;
         }
      else cumulative_ratio+=ratio[j+1];
      }
   }
if (userp->random_type==3)
	GetRand3(simple_rseq,seq1,num,userp->NbrOfElements);

return (simple_rseq);
}

// *** Function GetRand4: randomizes in three different frames according to
// nucleotide frequency ***

char* GetRand4(struct sequence *seq1, int *num, int n)
{
int i,j;
int r;
long z;
long partial_total[3];

long int all_frames;
float l;
l=int(float(seq1->length)/3); // divide the length in three parts
all_frames=(long int)(3.0*l); // part of the sequence that contains the three frames


float **ratio;
float *cumulative_ratio;
ratio = new float*[sizeof(float)*3];
for (i=0; i<3; i++)
	{
   ratio[i] = new float[sizeof(float)*n];
   }
cumulative_ratio = new float[sizeof(float)*3];

char *simple_rseq;
simple_rseq = new char[sizeof(char)*(seq1->length)];
if (!simple_rseq)
	{
	cerr << '\n' << "Not enough memory";
	exit(1);
	}

// define range of random numbers to use
int range=1000;

for (i=0; i<3; i++) // calculate partial totals (per frame)
	{
   partial_total[i]=0;
	for (j=0; j<n; j++)
		{
		partial_total[i]+=seq1->count_three_frames[i][j];
		}
   }

for (i=0; i<3; i++) // calculate ratios per frame
	{
	for (j=0; j<n; j++)
		{
		ratio[i][j]=(float(seq1->count_three_frames[i][j])/float(partial_total[i]))
      				*float(range);
		}
	}

// assign elements to random sequence
for (z=0; z<l; z++)
	{
   for (i=0; i<3; i++)
   	{
		r=rand() %range;
   	cumulative_ratio[i]=ratio[i][0];
   	for (j=0; j<n; j++)
   		{
      	if (r<cumulative_ratio[i])
      		{
         	simple_rseq[(3*z)+i]=seq1->elements[j];
         	num[j]++;
         	break;
         	}
      	else cumulative_ratio[i]+=ratio[i][j+1];
      	}
   	}
	}
for (i=0; i<(seq1->length-all_frames); i++)  // end up counting elements
	{
   r=rand() %range;
   cumulative_ratio[i]=ratio[i][0];
   for (j=0; j<n; j++)
   	{
      if (r<cumulative_ratio[i])
      	{
         simple_rseq[all_frames+i]=seq1->elements[z];
         num[j]++;
         break;
         }
      else cumulative_ratio[i]+=ratio[i][j+1];
      }
   }
return(simple_rseq);
}
// *** Function GetRand2: randomly gets nucleotides from sequence and
// assigns them to random sequence ***

char* GetRand2(struct sequence *seq1, int *num, int n)
{
	long i,l;
	int j;
	int r;
	char *simple_rseq;
	char *seq2;
	l=seq1->length;
	simple_rseq = new char[sizeof(char)*(seq1->length)];

	if (!simple_rseq)
    	{
		cerr << '\n' << "Not enough memory";
		exit(1);
    	}

	seq2 = new char[sizeof(char)*(seq1->length)];
	if (!seq2)
	{
		cerr << "Not enough memory";
		exit(1);
    	}


	strncpy(simple_rseq,seq1->seq,seq1->length);

	for (i=0; i<seq1->length-1; i++)
	{
		r=rand()%l;
	
		simple_rseq[i]^=simple_rseq[r];
		simple_rseq[r]^=simple_rseq[i];
		simple_rseq[i]^=simple_rseq[r];
        }
/*
	strncpy(seq2,seq1->seq,seq1->length);

	for (i=0; i<(seq1->length-1); i++)
	{
		r=rand()%l;
		simple_rseq[i]=seq2[r];

		strncpy(&seq2[r],&seq2[r+1],l-r);

    		l--;
    	}

	simple_rseq[seq1->length-1]=seq2[0];
*/
	delete seq2;
	return (simple_rseq);
}


// === Random_calculations: random sequences and parameters ===

void Random_sequence_calculations (struct sequence *seq1, struct param_process *userp,
 				               struct simple_scores *random, char RandomFile[],
                           char OutputFile[], char InputFile[], int maxs,
                           float sf_test)
{
	int i,j;
	float significance;
	// for calculating variance
	float sr, sr2, var, quadrat_sum;
	quadrat_sum=0;

	// initialize elements in struct simple_scores random
	random->score_frequency= new float[sizeof(float)*(maxs+1)];

	if (!random->score_frequency)
	{
   		cerr << '\n' << "Not enough memory";
   		exit(1);
   	}

	for (i=0; i<(maxs+1); i++)
	{
   		random->score_frequency[i]=0;
   	}

	random->SimplicityFactor=0;

	// initialize count elements in random sequences
	int *num;
	num = new int[sizeof(int)*(userp->NbrOfElements)];
	double sum_sd=0;

	// calculations on num_rand random sequences

	for (i=0; i<(userp->num_random); i++)
	{
   		// assign 0 to element count
		for (j=0; j<(userp->NbrOfElements); j++)
   		{
			num[j]=0;
		}

   		// get sequences and count number of elements
   		char* simple_rseq;

		if ((userp->random_type==1)||(userp->random_type==3))
    		{
			simple_rseq=GetRand1 (seq1, userp, num);
		}

		if (userp->random_type==2)
		{
			simple_rseq=GetRand2 (seq1, num, userp->NbrOfElements);
		}

   		if (userp->random_type==4)
   		{
      			simple_rseq=GetRand4 (seq1, num, userp->NbrOfElements);
      		}

   		// calculate simplicity parameters
   		sr=0;

   		// MODIFIED SJG 2005/08/10

    		double possd;

    		Simplicity_random_sequence(simple_rseq,userp,sr,random,maxs,&possd);

    		sum_sd+=possd;
    
   		// END OF MODIFICATION SJG 2005/08/10

   		// for calculating variance
		sr2=sr*sr;
		quadrat_sum=quadrat_sum + sr2;
		// Call to Output Function   ELIMINATED IN THIS VERSION
   		/* WriteOut_Random(sr,seq1->elements,i,num,RandomFile,InputFile,userp->NbrOfElements); */
   		// delete simple_rseq
		delete simple_rseq;
	}

	// Final values of struct simple_scores random (total/number random sequences)
	random->SimplicityFactor=random->SimplicityFactor/userp->num_random;

	for (i=0; i<(maxs+1); i++)
	{
   		random->score_frequency[i]=random->score_frequency[i]/userp->num_random;
   	}

	// ADDED SJG 2005/08/10

	// Average of standard devatiions

	double ave_sd=sum_sd/userp->num_random;

	// END OF ADDITION SJG 2005/08/10

	// Variance
	var=(quadrat_sum/userp->num_random)-(random->SimplicityFactor*random->SimplicityFactor);

	// Calculate significance of simplicity of test sequence
	significance=sf_test/random->SimplicityFactor;

	// Call to Output Function to plot variance and confidence intervals
	WriteOut_Variance_Confidence(var,ave_sd,random->SimplicityFactor,userp->num_random,
                             sf_test,significance,OutputFile); // SJG Added 'ave_sd' to parameter list 2005/08/1
}

//===================================================================
//
//            Output on parameters related to simplicity
//
//===================================================================

/*Output_motifs (motifs, position, count, number_elements, OutputMotifs,InputFile); */

// ==== Output Segments involved in high simplicity  =====

void Output_motifs_segments(char **motifs, struct sequence *seq1, struct param_process *userp,
							long *seq_high_score, int count, int number_elements,
                     char OutputMotifs[],char InputFile[])
{
int i,k;
long j;
int number_of_concatenations=0; // this will hold the concatenated positions with
long beginning;                 // high scores
long end;

FILE *File;
File=fopen(OutputMotifs,"w");
fprintf(File, "\n Simplicity-rich motifs and segments:");

if (count==0)
	fprintf(File,"\n No motifs with significant simplicity found!");
else
   {
   fprintf(File, "\n\n The motifs that accumulated high simplicity scores are shown.");
   fprintf(File, "\n The number in front of the motif indicates its position in the");
   fprintf(File, " sequence.");
   fprintf(File, "\n Below the motifs is the correspondent segment of the sequence where the repeats");
   fprintf(File, "\n can be found.");
	fprintf(File, "\n The motifs that accumulated high scores are marked with asterisks.");
   fprintf(File, "\n The abundance of the different motifs is also displayed at the end of the file.");
   fprintf(File, "\n\n");
	}

i=0;
while(i<count)
	{
   for (k=i+1; k<count; k++)
   	{
      if (seq_high_score[i]==(seq_high_score[k]-(k-i)))
      	{
         number_of_concatenations++;
         }
      else
      	break;
      }
   for (k=0; k<=number_of_concatenations; k++)
   	{
      fprintf(File," %5d  %s\n",seq_high_score[i+k]+1,motifs[i+k]);
      }
   // define beginning
   beginning=seq_high_score[i]-(userp->window);
	// print the part of the segment before the nucleotide/s with high score/s
   fprintf(File,"\n %6d- ",beginning+1);
   for (j=beginning; j<seq_high_score[i]; j++)
   	{
      fprintf(File,"%c",seq1->seq[j]);
      }
   // write the nucleotides with a high score
   for (k=0; k<=number_of_concatenations; k++)
   	{
      fprintf(File,"%c",seq1->seq[seq_high_score[i]+k]);
      }
   // define the end of the segment
	end=seq_high_score[i]+number_of_concatenations+(seq1->length-userp->max);
   // write the rest of the segment
   for (j=(seq_high_score[i]+number_of_concatenations+1);j<=end; j++)
   	{
      fprintf(File,"%c",seq1->seq[j]);
      }
   // write the position that corresponds to the end of this segment
   fprintf(File," -%5d\n",end+1);
   // underlie the nucleotides (repeats) with high scores
   for (j=0; j<41; j++)
   	fprintf(File,"%c",' ');
	for (j=0; j<(number_of_concatenations+1+(seq1->length)-(userp->max)-32); j++)
   	fprintf(File,"%c",'*');
   fprintf(File,"\n\n");
	// advance positions in i if there has been concatenations
   i=i+number_of_concatenations+1;
   number_of_concatenations=0;
   }

// count different motifs and keep position
int **freq_motifs;
freq_motifs=new int*[count];
for (i=0; i<count; i++)
	freq_motifs[i]=new int[2]; // the 2nd dimension holds position

int n=0;
int flag=TRUE;

for (i=0; i<count; i++)
	{
   for (j=i;j>0; j--)
   	{
      if (strncmp(motifs[i],motifs[j-1],number_elements)==0)
      	flag=FALSE;
      }
   if (flag==TRUE)
		{
      freq_motifs[n][0]=1;
   	for (k=i+1; k<count; k++)
     		{
      	if (strncmp(motifs[i],motifs[k],number_elements)==0)
      		freq_motifs[n][0]++;
      	}
   	freq_motifs[n][1]=i;
   	n++;
   	}
	flag=TRUE;
   }

// put abundance of motifs in decreasing order
long swap[2];
int flipped;

flipped=1;
while (flipped==1)
	{
   flipped=0;
   for (i=0; i<n-1; i++)
   	{
      if (freq_motifs[i][0]<freq_motifs[i+1][0])
      	{
         swap[0]=freq_motifs[i][0];
         swap[1]=freq_motifs[i][1];
         freq_motifs[i][0]=freq_motifs[i+1][0];
         freq_motifs[i][1]=freq_motifs[i+1][1];
         freq_motifs[i+1][0]=swap[0];
         freq_motifs[i+1][1]=swap[1];
         flipped=1;
         }
      }
   }

// Print out the frequency of motifs
if (count!=0)
	{
	fprintf(File, "\n\n Total  Motif\n");
	for (i=0; i<n; i++)
		{
   	fprintf(File, "\n %4d    %s",freq_motifs[i][0],motifs[freq_motifs[i][1]]);
   	}
   }

}


// =====  Frequency of scores in test and random sequences =====
void Output_score_frequencies (float *s, struct simple_scores *test,
										 struct simple_scores *random,char OutputScores[],
                               char InputFile[],float str)
{
int i;
FILE *File;
File=fopen(OutputScores,"w");
fprintf(File, "\n Input File: %s\n", InputFile);
fprintf(File, "\n Significance treshold (S=1-(freq score random/freq score test)): %5.3f",str);
fprintf(File, "\n Each position in the sequence has been assigned a simplicity score.");
fprintf(File, "\n Here scores in the test and random sequences are compared.");
fprintf(File, "\n Significant scores will correspond to positions in the test sequence");
fprintf(File, "\n associated with high simplicity.");
fprintf(File, "\n Motifs and segments of high simplicity are retrieved according to this.");
fprintf(File, "\n\n Score         Test           Random       Significant?\n");
for (i=0; i<=(test->MaxScore); i++)
	{
   fprintf(File, "\n %4d         %5.0f          %7.1f",i,test->score_frequency[i],
           random->score_frequency[i]);
	if (s[i]>=str)
   	fprintf(File, "           yes");
   else
   	fprintf(File, "           no");
   }
fclose(File);
}

//===================================================================
//
//            Calculations on significance of simplicity
//
//===================================================================



// ==== Find motifs involved in high simplicity =====
void Simple_motifs (struct sequence *seq1,struct param_process *userp,
                    long *seq_high_score,int count,char OutputMotifs[],
                    char InputFile[])
{
long i,j;
int k;
int number_elements=0;

// determine the number of elements of the repeat that have been considered
if (userp->score_tetra!=0)
	number_elements=4;
else if (userp->score_tri!=0)
	number_elements=3;
else if (userp->score_di!=0)
	number_elements=2;
else if (userp->score_mono!=0)
	number_elements=1;

// define variable that will hold the motifs
char **motifs;
motifs = new char*[sizeof(char)*(count)];
for (i=0; i<count; i++)
	{
   motifs[i] = new char[sizeof(char)*(number_elements+1)];
   }
if (!motifs)
	{
   cerr << '\n' << "Not enough memory";
   exit(1);
   }

// look for the motifs
for (i=0; i<count; i++)
	{
   for (j=userp->window; j<userp->max; j++)
   	{
      if (seq_high_score[i]==j)
      	{
         for (k=0; k<number_elements; k++)
         	{
            motifs[i][k]=seq1->seq[j+k];
            }
         motifs[i][number_elements]='\0';
         }
      }
   }
Output_motifs_segments(motifs, seq1, userp, seq_high_score, count, number_elements,
					        OutputMotifs,InputFile);
delete motifs;
}


// ==== Find significant scores, motifs and segments in test sequence =====
void Parameters_simplicity (struct sequence *seq1,struct param_process *userp,
										struct simple_scores *test,
                              struct simple_scores *random,char OutputScores[],
                              char OutputMotifs[], char InputFile[])

{
int i;
long j;
long count=0;

// In the first place we will search for those scores in the test sequence
// which are significantly high, that is which ocurr at least 10 times more
// frequently in the test sequence than in the random sequences (averaged)
// or which are only present in the test sequence
// This will be calculated through parameter "significance" (s), if the score
// is significant then the positions in the sequence with that score will be
// stored in seq_high_score

// initialize significance
float *s;
s=new float[sizeof(float)*(test->MaxScore+1)];
if (!s)
	{
   cerr << '\n' << "Not enough memory";
   exit(1);
   }

// initialize variable that will hold the positions with significant scores
long *seq_high_score;
seq_high_score=new long[sizeof(long)*(MAX)];
if (!seq_high_score)
	{
   cerr << '\n' << "Not enough memory";
   exit(1);
   }

// assign the positions with high scores to seq_high_score
for (i=0; i<=(test->MaxScore); i++)
	{
   if (test->score_frequency[i]==0)
		s[i]=0;
	else
   	s[i]=1-(random->score_frequency[i]/test->score_frequency[i]);
   if (s[i]>=(userp->stringency))  // in case it is significant look for the positions
   	{
   	for (j=userp->window; j<userp->max; j++)
      	{
         if (test->ss_test[j]==i)
            {
         	seq_high_score[count]=j;
            count++;
            }
         }
      }
   }

// put seq_high_score in order of appearance in the sequence
long swap;
int flipped;

flipped=1;
while (flipped==1)
	{
   flipped=0;
   for (i=0; i<count-1; i++)
   	{
      if (seq_high_score[i]>seq_high_score[i+1])
      	{
         swap=seq_high_score[i];
         seq_high_score[i]=seq_high_score[i+1];
         seq_high_score[i+1]=swap;
         flipped=1;
         }
      }
   }


// free test->ss
delete test->ss_test;

// Output score_frequencies  ELIMINATED IN THIS VERSION
Output_score_frequencies (s,test,random,OutputScores,InputFile,userp->stringency);

// free s and score frequencies
delete s;
delete test->score_frequency;
delete random->score_frequency;

// Find motifs which are involved in high simplicity
Simple_motifs (seq1,userp,seq_high_score,count,OutputMotifs,InputFile);

// free seq_high_score
delete seq_high_score;
}





