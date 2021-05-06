/*#############################################################################

  # to get input data for CN from Battenberg format file for BayClone2

#############################################################################*/

#include "stdio.h"
#include "assert.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"


#define MAX_LINE_LEN	50000
#define MAX_FLD_LEN		500
#define SAMPLE_OFFSET	10
#define FORMAT_INDX		8
#define VCF_CHR_NAME_INDX	0
#define VCF_POS_INDX	1
#define VCF_PASS_INDX	6


#define BB_CHR_NAME_INDX	0
#define BB_CHR_POS_STRT_INDX 1
#define BB_CHR_POS_END_INDX 2
#define BB_MAJ1A_INDX	7
#define BB_MIN1A_INDX	8
#define	BB_FRAC1A_INDX	9
#define BB_MAJ2A_INDX	10
#define BB_MIN2A_INDX	11
#define	BB_FRAC2A_INDX	12


typedef struct _SNV_CN
{
	// location info
	char chrName[10];
	int startPos;
	int endPos;

	// from mutect vcf file
	int N_tot_T; // tumor total cnt
	int n_alt_T; // tumor alt cnt
	int N_tot_N; // normal total cnt
	int n_alt_N; // tumor alt cnt

	// from BB file
	int nMaj1_A;
	int nMin1_A;
	double frac1_A;
	int nMaj2_A;
	int nMin2_A;
	double frac2_A;
	double sampleCN; // our input M for BayClone2


}SNV_CN;


int readNextLine (char *l, FILE *f)
{
	if (fgets (l, MAX_LINE_LEN, f) == NULL)
	{
		/* NULL.  EOF was read.  reset l, return */
		l[0] = '\0';
	   	return 1;
	}
	else
	{
		if (strlen (l) == MAX_LINE_LEN -1)
		{
			/* Line too long - do some malloc stuff.
	 		 * For now, just terminate and return. */
 	 		//error("Encountered line longer than %d:\n%s\n\nEnding ...\n", MAX_LINE_LEN, l);
 	 		fprintf(stdout,"Encountered line longer than %d:\n%s\n\nEnding ...\n", MAX_LINE_LEN, l);
	     	l[0] = '\0';
	     	return 2;
		}
		else
		{
			/* Proper line read. Return 0 */
	    	return 0;
	   	}
	}
}

int file_line_count(char *file_name)
{
	char c;
	int lines = 0;
	FILE *fp;
	//fprintf(stdout,"fileName: %s\n",file_name);
	fp = fopen(file_name,"r");
	if(NULL==fp)
	{
		//fprintf(stdout,"BB file open error !!\n");
		return 0;
	}
	//fprintf(stdout,"reading files...\n");
	while ((c = fgetc(fp)) != EOF)
	{
		if (c == '\n')
			++lines;
	}

	fclose(fp);
	return lines;
}

int export_N_n(char* vcf_name, SNV_CN **SC1, int* VCF_n_rcd)
{

	//Rprintf("%s:%s:%s\n",*R_vcf_name,*R_N_tot_file,*R_n_alt_file);

	//char vcf_name[100];
	//char N_tot_file[100];
	//char n_alt_file[100];

	FILE* fp;
	FILE* ft1,* ft2;
	char vcf_line[MAX_LINE_LEN], header_line[MAX_LINE_LEN];
	char* pch,*pch1;
	//char *pch2;
	int cnt1,s_cnt;
	int GT_fld;
	int AD_fld;
	char fmt_fld[MAX_FLD_LEN];
	char pass_str[MAX_FLD_LEN];
	int i,j;
	int e_flag;
	//int N_tot;
	int n_alt,N_minus_n;
	char N_n_cnt_str[MAX_FLD_LEN];
	int N_tot_samples,n_alt_samples;
	int total_line_cnt;
	int header_line_cnt = 0;
	int snv_line_cnt;

	char chrName[10];
	int posSNV;
	char* pass_pch;

	//strcpy(vcf_name,R_vcf_name);
	//strcpy(N_tot_file,R_N_tot_file);
	//strcpy(n_alt_file,R_n_alt_file);

	total_line_cnt = file_line_count(vcf_name);
	//fprintf(stdout,"%d lines\n",total_line_cnt);

	fp = fopen(vcf_name,"r");

	//ft1=fopen(N_tot_file,"w");
	//ft2=fopen(n_alt_file,"w");

	int vcf_line_cnt = 0;
	if(NULL==fp)
		//||(NULL==ft1)||(NULL==ft2))
	{
		//error("file open error !\n");
		fprintf(stdout,"VCF file open error !!\n");
		return -1;
	}
	else
	{
		while (readNextLine(vcf_line, fp) == 0)
		{
			/* If it is a comment, keep track as header line */
	     	if(vcf_line[0] == '#')
	     	{
	     		header_line_cnt++;
	     	}

	     	if (vcf_line[0] != '#')
			{
				//fprintf(stderr, "%s\n", vcf_line);
				e_flag = 0;
			// 	header_line_cnt++;
   			//  strcpy (header_line, vcf_line);
  		// 	}
			// /* If not a comment, break out of this loop */
  		// 	else
  		// 	{
    			vcf_line_cnt++;
				if(vcf_line_cnt == 1)
				{
				 	snv_line_cnt = total_line_cnt - header_line_cnt + 1;
				 	*SC1 = (SNV_CN*)(calloc(snv_line_cnt,sizeof(SNV_CN)));
				}
				// 	//fprintf(stdout,"%d SNV lines are there\n",snv_line_cnt);
				// 	// we can now find number of samples
				// 	cnt1=0;
				// 	pch = strtok (header_line,"\t");
				// 	while (pch != NULL)
				// 	{
				// 		//Rprintf("%s\n",pch);
				// 		pch = strtok (NULL, "\t");
				// 		cnt1++;
				// 	}
				// 	s_cnt = cnt1-SAMPLE_OFFSET;

				// 	N_n_cnt_str=(char**)(calloc(s_cnt,sizeof(char*)));
				// 	for(i=0;i<s_cnt;i++)
				// 		N_n_cnt_str[i] = (char*)(calloc(1,MAX_FLD_LEN*sizeof(char)));

				// 	N_tot_samples = (int*)(calloc(s_cnt,sizeof(int)));
				// 	n_alt_samples = (int*)(calloc(s_cnt,sizeof(int)));
				// }
				cnt1=0;
				j=0;
				pch = strtok(vcf_line,"\t");
				while (pch != NULL)
				{
					if(cnt1==VCF_CHR_NAME_INDX)
						strcpy(chrName,pch);

					if(cnt1==VCF_POS_INDX)
						sscanf(pch,"%d",&posSNV);

					if(cnt1==VCF_PASS_INDX)
					{
						strcpy(pass_str,pch);
						pass_pch = strstr(pass_str,"PASS");
						if(pass_pch == NULL)
							e_flag = 1;
					}
					if(cnt1==FORMAT_INDX)
						strcpy(fmt_fld,pch);

					if(cnt1==SAMPLE_OFFSET)
						strcpy(N_n_cnt_str,pch);

					pch = strtok (NULL, "\t");
					cnt1++;
				}

				cnt1=0;
				pch = strtok (fmt_fld,":");
				while (pch != NULL)
				{
					if(strcmp("GT",pch)==0)
						GT_fld=cnt1;
					if(strcmp("AD",pch)==0)
						AD_fld=cnt1;

					pch = strtok (NULL, ":");
					cnt1++;
				}


				j = 0;
				pch1 = strtok(N_n_cnt_str,":");
				while (pch1 != NULL)
				{
					if(j==AD_fld)
						sscanf(pch1,"%d,%d",&N_minus_n,&n_alt);

					pch1 = strtok (NULL, ":");
					j++;

				}

				N_tot_samples = n_alt+N_minus_n;
				n_alt_samples = n_alt;
				if(N_tot_samples == 0)
					e_flag = 1;
				//printf(stdout,"%d\t%d\n",N_tot_samples,e_flag);
				if(e_flag==0)
				{
					int k = (*VCF_n_rcd);
					//fprintf(stderr,"%s\n",chrName);
					strcpy((*SC1)[k].chrName,chrName);
					(*SC1)[k].startPos = posSNV;
					(*SC1)[k].endPos = posSNV;
					(*SC1)[k].N_tot_T = N_tot_samples;
					(*SC1)[k].n_alt_T = n_alt_samples;
					(*SC1)[k].N_tot_N = 0;
					(*SC1)[k].n_alt_N = 0;

					/*
					fprintf(stdout,"%s:%d\t",chrName,posSNV);
					for(i=0;i<s_cnt;i++)
					{
						fprintf(stdout,"%d:%d\t",N_tot_samples[i],n_alt_samples[i]);
					}
					fprintf(stdout,"\n");
					*/
					(*VCF_n_rcd)++;
				}
			}
		}

	}// else file open success

	fclose(fp);
	//fclose(ft1);
	//fclose(ft2);
	return 0;
	//Rprintf("files are generated!\n");
	//fprintf(stdout,"files are generated!\n");

}

int get_sample_copy_val_for_SNVs(char* BB_file_name,SNV_CN** SC2,int* BB_n_rcd)
{
	FILE* fp;
	char BB_line[MAX_LINE_LEN];
	//fprintf(stdout,"BB file Name: %s\n",BB_file_name);

	fp = fopen(BB_file_name,"r");
	if(NULL == fp)
	{
		fprintf(stdout,"BB input file open error!!\n");
		return -1;
	}
	int total_line_cnt = file_line_count(BB_file_name);
	//fprintf(stdout,"BB line cnt: %d\n",total_line_cnt);

	*SC2 = (SNV_CN*)(calloc(total_line_cnt,sizeof(SNV_CN)));
	int cnt1,i=0;
	char* pch;
	int nLine = 0;
	while (readNextLine(BB_line, fp) == 0)
	{
			nLine++;
			//fprintf(stdout,"%s\n",BB_line);
			/* If it is a comment, keep track as header line */
	     	if (BB_line[0] == '#' || nLine == 1)
				continue;
			else
			{
				cnt1=0;
				pch = strtok(BB_line,"\t");
				while (pch != NULL)
				{
					//Rprintf("%s\n",pch);
					switch(cnt1)
					{
						case BB_CHR_NAME_INDX:
							strcpy((*SC2)[i].chrName,pch);break;
						case BB_CHR_POS_STRT_INDX:
							sscanf(pch,"%d",&((*SC2)[i].startPos));break;
						case BB_CHR_POS_END_INDX:
							sscanf(pch,"%d",&((*SC2)[i].endPos));break;
						case BB_MAJ1A_INDX:
							sscanf(pch,"%d",&((*SC2)[i].nMaj1_A));break;
						case BB_MIN1A_INDX:
							sscanf(pch,"%d",&((*SC2)[i].nMin1_A));break;
						case BB_FRAC1A_INDX:
							sscanf(pch,"%lf",&((*SC2)[i].frac1_A));break;
						//default: fprintf(stdout,"BB file parsing error !!\n");break;
					}
					if((*SC2)[i].frac1_A < 1.0)
					{
						switch(cnt1)
						{
							case BB_MAJ2A_INDX:
								sscanf(pch,"%d",&((*SC2)[i].nMaj2_A));break;
							case BB_MIN2A_INDX:
								sscanf(pch,"%d",&((*SC2)[i].nMin2_A));break;
							case BB_FRAC2A_INDX:
								sscanf(pch,"%lf",&((*SC2)[i].frac2_A));break;
						}
					}
					else
					{
						(*SC2)[i].nMaj2_A = 0;
						(*SC2)[i].nMin2_A = 0;
						(*SC2)[i].frac2_A = 0.0;

					}
					pch = strtok (NULL, "\t");
					cnt1++;
				}
				i++;
			}

	}

	*BB_n_rcd = i;

	/*
	for (int j=0;j<(*BB_n_rcd);j++)
	{
			fprintf(stdout,"%s:%d-%d\t", (*SC2)[j].chrName,(*SC2)[j].startPos,(*SC2)[j].endPos);

			fprintf(stdout,"(%d+%d)*%lf+(%d+%d)*%lf = ",SC2[j].nMaj1_A,SC2[j].nMin1_A,SC2[j].frac1_A,
													  SC2[j].nMaj2_A,SC2[j].nMin2_A,SC2[j].frac2_A);
			SC2[j].sampleCN =  (SC2[j].nMaj1_A+SC2[j].nMin1_A)*SC2[j].frac1_A +
							(SC2[j].nMaj2_A+SC2[j].nMin2_A)*SC2[j].frac2_A;
			fprintf(stdout,"%lf\n",SC2[j].sampleCN);

	}
	*/
	fclose(fp);
	return 0;
}

int main(int argc,char* argv[])
{
	SNV_CN *SC1, *SC2;
	int BB_n_rcd = 0;
	int VCF_n_rcd = 0;
	FILE* ft;

	if(argc < 3)
	{
		fprintf(stdout,"Please enter: <prog_name> <input vcf file> <input BB file>\n");
		return -1;
	}
	else
	{

		ft = fopen("BayClone3_input.txt","w");
		if(NULL == ft)
		{
			fprintf(stdout,"output file open error!!\n");
			return -1;
		}

		int iRet = get_sample_copy_val_for_SNVs(argv[2],&SC2,&BB_n_rcd);
		if(iRet == -1)
			return -1;

		fprintf(stdout,"%d BB records\n",BB_n_rcd);

		for (int j=0;j<BB_n_rcd;j++)
		{

			//fprintf(stdout,"%s:%d-%d\t", SC2[j].chrName,SC2[j].startPos,SC2[j].endPos);

			//fprintf(stdout,"(%d+%d)*%lf+(%d+%d)*%lf = ",SC2[j].nMaj1_A,SC2[j].nMin1_A,SC2[j].frac1_A,
			//										  SC2[j].nMaj2_A,SC2[j].nMin2_A,SC2[j].frac2_A);
			SC2[j].sampleCN =  (SC2[j].nMaj1_A+SC2[j].nMin1_A)*SC2[j].frac1_A +
							   (SC2[j].nMaj2_A+SC2[j].nMin2_A)*SC2[j].frac2_A;
			///actually this involves CCF of that CNA region; so this is sampleCN in the tumor fraction
			//fprintf(stdout,"%lf\n",SC2[j].sampleCN);

		}


		iRet = export_N_n(argv[1],&SC1,&VCF_n_rcd);
		if(iRet == -1)
			return -1;

		fprintf(stdout,"%d VCF records\n",VCF_n_rcd);

		/*
		for(int j=0;j<VCF_n_rcd;j++)
		{
			fprintf(stdout,"%s:%d-%d\t",SC1[j].chrName,SC1[j].startPos,SC1[j].endPos);
			fprintf(stdout,"%d:%d\t%d:%d\n",SC1[j].N_tot_T,SC1[j].n_alt_T,SC1[j].N_tot_N,SC1[j].n_alt_N);
		}
		*/

		//fprintf(ft,"#chr\tpos\ttumor_N\ttumor_n\tsampleCN_tumor\n");

		for(int j=0;j<VCF_n_rcd;j++)
		{
			SC1[j].sampleCN = 2.0;
			//SC1[j].sampleCN = -1.0; // when only clonal CN info is avalable
			int loc1 = SC1[j].startPos;
			for(int k=0;k<BB_n_rcd;k++)
			{
				if( (strcmp(SC1[j].chrName,SC2[k].chrName)==0) && (loc1 >= SC2[k].startPos) && (loc1 <= SC2[k].endPos) )
				{
					SC1[j].sampleCN = SC2[k].sampleCN;

				}
			}
			//fprintf(stdout,"%s:%d\t",SC1[j].chrName,SC1[j].startPos);
			//fprintf(stdout,"%d\t%d\t%d\t%d\t",SC1[j].N_tot_T,SC1[j].n_alt_T,SC1[j].N_tot_N,SC1[j].n_alt_N);
			//fprintf(stdout,"%lf\n",SC1[j].sampleCN);

			fprintf(ft,"%s\t%d\t",SC1[j].chrName,SC1[j].startPos);
			//fprintf(ft,"%d\t%d\t%d\t%d\t",SC1[j].N_tot_T,SC1[j].n_alt_T,SC1[j].N_tot_N,SC1[j].n_alt_N);
			fprintf(ft,"%d\t%d\t",SC1[j].N_tot_T,SC1[j].n_alt_T);

			/// just to test new -1 input
			fprintf(ft,"%lf\n",SC1[j].sampleCN);

		}

		fprintf(stdout,"outout file is written !!\n");
	}


	fclose(ft);

	free(SC1);
	free(SC2);


	return 0;
}
