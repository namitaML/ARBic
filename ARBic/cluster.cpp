/***********************************************************************/

#include "cluster.h"
#include"struct.h"
#include "dag_longest_path.h"
#include <iostream>
using namespace std;

/************************************************************************/

static void update_colcand(bool *colcand, const char *s1)
{	
	int i=0;

	for (i=0; i< cols; i++)
		if (colcand[i] &&(s1[i]==0))
			colcand[i] = FALSE;
}
/**********************/

/*********************/

/**********************/

/**********************************************************/
static int intersect_row_LCS(const bool *colcand, char *colcand2)
/*caculate the weight of the edge*/
{
	int i;
	int cn = 0;
	for (i=0; i< cols; i++)
	{	if (colcand[i] && colcand2[i]!=0) 
			cn++;
	}
	return cn;
}

/***********************************************************/
static int seed_current(struct dyStack *genes, bool *colcand, int components,int *colsStat)
/* calculate the coverage of any row to the current consensus
 * cnt = # of valid consensus columns
 */
{
	std::cout<<"components= "<<components<<std::endl;
	int i,j;
	int cnt =0;
	int speedup_par=1;
	float col_par=0.8;
	if (components>300)
	{
		speedup_par=components/300;
	}
	if(rows>1000)col_par=1;
	//int threshold = floor(components * po->TOLERANCE)-1;//g//lxy
	//int threshold = floor(components *0.8/speedup_par)-1;//lxy
	int threshold = floor(components *po->TOLERANCE/speedup_par)*col_par-1;//lxy
	if(threshold <1)
		threshold=1;
	/*get the statistical results of each column produced by seed*/
	char *temptag = (char *)xmalloc(cols*sizeof(char));
	char *temptag2 = (char *)xmalloc(cols*sizeof(char));
	for(i=0;i<cols;i++)
	{	
		colsStat[i] = 0;
		temptag[i] = 0;
	}
	int gene_loc,gene_loc2;
	for(i=0;i<components/speedup_par-1;i++)
	{	
		gene_loc=i;
		if (components>300)
		{
			gene_loc=rand()%components;
			gene_loc2=rand()%components;
		}
		vector<int>a(0);
		get_Genes_LCS(arr_c[dsItem(genes,gene_loc)], arr_c[dsItem(genes,gene_loc+1)],temptag,a);
		vector<int>(0).swap(a);
		for(j=0;j<cols;j++)
		{
			if(temptag[j]!=0)
				colsStat[j]++;
			temptag[j]=0;
		}	
	}
	for(i=0;i<cols;i++)
	{	
		if (colsStat[i] >= threshold)
		{
			colcand[i] = TRUE; 
			cnt++;
		}
	}
	free(temptag);
	free(temptag2);
	return cnt;
}

static bool check_seed(Edge *e, Block **bb, const int block_id)
/*check whether current edge can be treat as a seed*/
{
	int profiles[rows];
	int i,b1,b2,b3;
	bool fg = FALSE;
	b1 = b2 = -1;
	for (i = 0; i < block_id; i++)
		if ( isInStack(bb[i]->genes,e->gene_one) && isInStack(bb[i]->genes, e->gene_two) ) 
			return FALSE; 
	return TRUE;
	/*
	for ( i = 0; i < rows; i++) profiles[i] = 0;
	fg = FALSE;	
	for ( i = 0; i < block_id; i++)
		if ( isInStack(bb[i]->genes, e->gene_one) ) 
		{ 
			fg = TRUE;
		       	break; 
		}
	if (fg) 
		b1 = i;
	fg = FALSE;	
	for ( i = 0; i < block_id; i++)
		if ( isInStack(bb[i]->genes, e->gene_two) ) 
		{ 
			fg = TRUE; 
			break; 
		}
	if (fg) 
		b2 = i;
	if ( (b1 == -1)||(b2 == -1) ) 
		return TRUE;
	else
	{
		for ( i = 0; i < bb[b1]->block_rows; i++)
			profiles[dsItem(bb[b1]->genes,i)]++;
		for ( i = 0; i < bb[b2]->block_rows; i++)
			profiles[dsItem(bb[b2]->genes,i)]++;
		for ( i = 0; i < rows; i++)
 			if (profiles[i] > 1) 
				return FALSE;
		b3 = MAX(bb[b1]->block_cols, bb[b2]->block_cols);
		if ( e->score <b3)// (bb[b1]->block_cols + bb[b2]->block_cols) / 2  
			return FALSE;
		else 
			return TRUE;
	}

	err("never see this message\n");
	return FALSE;
	*/
}
//###################33

int most_repeat(vector<int>a)

{
	vector<int>b(cols,0);
	int max_count=0;
	for(int i=0;i<a.size();i++)
		{
			b[a[i]]++;
		}
	for(int i=0;i<b.size();i++)
		{
			if(b[i]>max_count)
			{
				max_count=b[i];
			}
		}
	vector<int>(0).swap(b);
return max_count;
}









   
/************************************************/
static void block_init(Edge *e, Block *b, 
                     struct dyStack *genes, struct dyStack *scores,
                     bool *candidates, const int cand_threshold,
                     int *components, struct dyStack *allincluster, long double *pvalues)
{
	int i,score,j;
	int cnt = 0, cnt_all=0, pid=0;
	continuous cnt_ave=0, row_all = rows;
	long double pvalue;
	int max_cnt, max_i;
	int t0,t1;
	int *arr_rows, *arr_rows_b;
        arr_rows = (int*) malloc(sizeof(int) * rows);
        arr_rows_b = (int*) malloc(sizeof(int) * rows);
	bool *colcand;
        colcand = (bool *)malloc(sizeof(bool) * cols);
        for (i=0; i< cols; i++) 
		colcand[i] = FALSE;
	discrete *g1, *g2;
	t0=dsItem(genes,0);
	t1=dsItem(genes,1);
	g1 = arr_c[t0];
	g2 = arr_c[t1];
	
	for(i=0;i<rows;i++)	
	{
		lcs_length[i]=0;
		for(j=0;j<cols;j++)
			lcs_tags[i][j]=0;
	}
	
	/*************************calculate the lcs*********************************/
	vector<int>a;
	lcs_length[t1]=get_Genes_LCS(g1,g2,lcs_tags[t1],a); 
	//cout<<"a.size= "<<a.size()<<endl;
	for(i=0;i<cols;i++)
	{
		if(lcs_tags[t1][i]!=0)
			{
			colcand[i]=TRUE;
			//std::cout<<col_candi_matrix[0][i]<<std::endl;
			}
	}
	discrete *gm=new discrete[cols];
	for(i=0;i<a.size();i++)
	{
		gm[a[i]]=i;
	}
	vector<int>(0).swap(a);
	for(j=0;j<rows;j++)
	{	
		if (j==t1 || j==t0) continue;
		lcs_length[j]= get_Genes_LCS_with_lcs(gm,arr_c[j],lcs_tags[j],lcs_tags[t1]); 
	}
	delete [] gm;
/*	printf("\n==chose genes from MAXT:");*/
		int last_max_cnt=lcs_length[t1];
	while (*components < rows)
	{	
		max_cnt = -1;
		max_i = -1;
		(*components)++;
		int a_temp=*components;
		//std::cout<<a_temp<<"...."<<std::endl;
		cnt_all =0;
		cnt_ave = 0;
		/******************************************************/
		/*add a function of controling the bicluster by pvalue*/
		/******************************************************/
		for (i=0; i< rows; i++)
		{
			if (!candidates[i]) continue;
			if (po->IS_list && !sublist[i]) continue;
			cnt = intersect_row_LCS(colcand,lcs_tags[i]);
			
		//	cnt=get_Genes_LCS_with_lcs(gm,arr_c[i],lcs_tag[i],colcand_tag)
			cnt_all += cnt;
		//	if (cnt < 43)
			if (cnt < cand_threshold) 
			//if (cnt < last_max_cnt/10) 
				candidates[i] = FALSE;
			if (cnt > max_cnt)
			{
				max_cnt = cnt;
				max_i = i;
			}
		}
	/**	printf("\t%d,%d",lcs_rowNo[max_i],max_cnt);*******************/

		cnt_ave = cnt_all/row_all;
		//std::cout<<"cnt_ave= "<<cnt_ave<<std::endl;
		//std::cout<<"max_cnt= "<<max_cnt<<std::endl;
		pvalue = get_pvalue (cnt_ave, max_cnt);
		if (po->IS_cond)
		{
			if (max_cnt < po->COL_WIDTH || max_i < 0|| max_cnt < b->cond_low_bound) break;
		}
		else
		{
			if (max_cnt < po->COL_WIDTH || max_i < 0) break;
		}



		if (po->IS_area)
			score = *components*max_cnt;
		else
			score = MIN(*components, max_cnt);
		if (score > b->score)
			b->score = score;
		if (pvalue < b->pvalue)
			b->pvalue = pvalue;
		last_max_cnt=max_cnt;
		dsPush(genes, max_i);
		dsPush(scores,score);
		pvalues[pid++] = pvalue;
		update_colcand(colcand, lcs_tags[max_i]);
		///
		///
		candidates[max_i] = FALSE;
	}
	/*be sure to free a pointer when you finish using it*/
	free(colcand);
	free(arr_rows);
	free(arr_rows_b);
	//std::cout<<"pvalue= "<<pvalue<<std::endl;
}
/************************************************************************/
/* Core algorithm */
int cluster (FILE *fw, Edge **el, int n/*num of seeds */)
{
	int block_id = 0;
	Block **bb;
	int cnt = 0;
	int allocated = po->SCH_BLOCK;
//lxy
//        bb = (Block **)malloc(sizeof(*bb) * allocated);
	AllocArray(bb,Block **, allocated);//lxy
//
	Edge *e;
	Block *b;
	struct dyStack *genes, *scores, *b_genes, *allincluster;
	
	int i, j, k, components;
	
//	int *colsStat = (int*) malloc(sizeof(int) * cols);//lxy 
	int *colsStat;
	AllocArray(colsStat,int *,cols);	

	genes = dsNew(rows);
	scores = dsNew(rows);
	allincluster = dsNew(rows);

//        long double *pvalues = (long double *)malloc(sizeof(*pvalues) * rows);//lxy
	long double *pvalues;
	AllocArray(pvalues, long double *,rows);

        bool *candidates = (bool *)malloc(sizeof(bool) * rows);
	e = *el;

        lcs_length = (discrete *)malloc(sizeof(discrete) * rows);
        lcs_tags = (char **)malloc(sizeof(char *) * rows);

        for(i=0;i<rows;i++)
	{
          lcs_tags[i] = (char *)malloc(cols*sizeof(char));
        }	
	//

	//
	i = 0;
	while (i++ < n)
	{	
		e = *el++;
		/* check if both genes already enumerated in previous blocks */
		bool flag = TRUE;
		/* speed up the program if the rows bigger than 200 */
	        if (rows > 250)
		{ 
			if ( isInStack(allincluster,e->gene_one) && isInStack(allincluster,e->gene_two) )
				{flag = FALSE;}
			else if ((po->IS_TFname)&&(e->gene_one!= TFindex)&&(e->gene_two!=TFindex))
				flag = FALSE;
			else if ((po->IS_list)&&(!sublist[e->gene_one] || !sublist[e->gene_two]))
				flag =FALSE;
		}
		else   
		{
			flag = check_seed(e, bb, block_id);
			if ((po->IS_TFname)&&(e->gene_one!= TFindex)&&(e->gene_two!=TFindex))
				flag = FALSE;
			if ((po->IS_list)&&(!sublist[e->gene_one] || !sublist[e->gene_two]))
				flag = FALSE;
		}
		/*
		if(e->gene_one > rows || e->gene_two > rows)
			flag=FALSE;
			*/
		if (!flag) continue;
		
		/*you must allocate a struct if you want to use the pointers related to it*/
		b = (Block *)malloc(sizeof(*b));
		/*initial the b->score*/
                b->score = MIN(2, e->score);
		/*initial the b->pvalue*/
		b->pvalue = 1;
		/* initialize the stacks genes and scores */		
		int ii;		
		dsClear(genes);
		dsClear(scores);		
		for(ii = 0; ii < rows; ii ++)
		{
			dsPush(genes,-1);
			dsPush(scores,-1);
		}		
		dsClear(genes);
		dsClear(scores);
		
		dsPush(genes, e->gene_one);
		dsPush(genes, e->gene_two);

		dsPush(scores, 1);
		dsPush(scores, b->score);

		/* branch-and-cut condition for seed expansion */
		//int cand_threshold = floor(po->COL_WIDTH * po->TOLERANCE);
		int cand_threshold = floor(po->COL_WIDTH * po->TOLERANCE);/////
                if (cand_threshold < 2) 
			cand_threshold = 2;

		/* maintain a candidate list to avoid looping through all rows */		
		for (j = 0; j < rows; j++) 
			candidates[j] = TRUE;
		candidates[e->gene_one] = candidates[e->gene_two] = FALSE;
		components = 2;
		/* expansion step, generate a bicluster without noise */
		block_init(e, b, genes, scores, candidates, cand_threshold, &components, allincluster, pvalues);
		/* track back to find the genes by which we get the best score*/
		for(k = 0; k < components; k++)
		{
			if (po->IS_pvalue)
				if ((pvalues[k] == b->pvalue) &&(k >= 2) &&(dsItem(scores,k)!=dsItem(scores,k+1))) break;
			{if ((dsItem(scores,k) == b->score)&&(dsItem(scores,k+1)!= b->score)) break;}
		}
		components = k + 1;
		////
		if(components<5)
		{continue;}
	////
		std::cout<<"components= "<<components<<std::endl;
		int ki;
		for (ki=0; ki < rows; ki++)
		{	
			candidates[ki] = TRUE;
		}

		for (ki=0; ki < components - 1 ; ki++)
		{
			candidates[dsItem(genes,ki)] = FALSE;
		}

		candidates[dsItem(genes,k)] = FALSE;
		genes->top = k;
		
		bool *colcand = (bool *) malloc(sizeof(bool) * cols);
		for(ki = 0; ki < cols; ki++) 	
		{
			colcand[ki] = FALSE;
		} 
			//colcand[ki] = FALSE;
		/*       
		for(ki = 0; ki < cols; ki++) 	
		{
			colcand[ki] = FALSE;
			if(colcand[ki])
		{std::cout<<"ki "<<ki<<std::endl;}
		}
		*/
		/* get init block */ 
		cnt = seed_current(genes, colcand, components,colsStat);
		/* add some new possible genes */

		int m_ct=0;
		bool colChose = TRUE;
		for(ki=0;ki < rows;ki++)
		{
			colChose=TRUE;
			if ((po->IS_list && !sublist[ki]) || !candidates[ki])
				continue;
			//m_ct = intersect_row_LCS(colcand,lcs_tags[ki]);
//
			vector<int>a(0);
			get_Genes_LCS(arr_c[dsItem(genes,0)],arr_c[dsItem(genes,1)],lcs_tags[dsItem(genes,1)],a);
			
			discrete *gm=new discrete[cols];
       			 for(int gmi=0;gmi<a.size();gmi++)
       			 {
           			 gm[a[gmi]]=cols-gmi;
    
   			     }
        vector<int>(0).swap(a);
			char* gm_tem_tag=new char[cols];
			for(int gmi=0;gmi<cols;gmi++)
			{gm_tem_tag[gmi]=0;}
			get_Genes_LCS_with_lcs(gm,arr_c[ki],gm_tem_tag,lcs_tags[dsItem(genes,1)]);						    		
			m_ct = intersect_row_LCS(colcand,gm_tem_tag);
			delete [] gm;
			delete []gm_tem_tag;
//			m_ct=get_Genes_LCS_with_lcs(arr_c[dsItem(genes,0)],arr_c[ki],colcand);

			if (candidates[ki]&& (m_ct >= floor(cnt * po->TOLERANCE)))
//			if (candidates[ki]&& (m_ct >= floor(cnt * 0.85)-1))
			{std::cout<<"new++"<<ki<<std::endl;
				int temp;
				for(temp=0;temp<cols;temp++)
				{
					if(colcand[temp])
					{	int tmpcount = colsStat[temp];
						//int tmpcount = colsStat[temp];
						if(lcs_tags[ki][temp]!=0)
							tmpcount++;
						//if(tmpcount < floor(components *(1- po->TOLERANCE)))//g//lxy
						if(tmpcount < floor(components * 0.1)+1)//lxy
						{	
							colChose = FALSE;
							break;
						}
					}
				}
				if(colChose==TRUE)			
				{	
					dsPush(genes,ki);
			
					//std::cout<<"new+ ="<<ki<<std::endl;
					components++;
					candidates[ki] = FALSE;
					for(temp=0;temp<cols;temp++)
					{
						if(lcs_tags[ki][temp]!=0 && colcand[ki])
						{
							colsStat[temp]++;
						}
					}
				}
			}
		}
						
                b->block_rows_pre = components;
		// add genes that negative regulated to the consensus 
		char * reve_tag;
		for ( ki = 0; ki < rows; ki++)
		{
			colChose=TRUE;
			if (po->IS_list && !sublist[ki] && !candidates[ki]) continue;
				
                reve_tag = (char *)malloc(cols*sizeof(char));
                int kk = 0;
                for(kk=0;kk<cols;kk++)
				reve_tag[kk]=0;
			//if(str_intersect_r(arr_c[dsItem(genes,0)],arr_c[ki],'N') < floor(cnt * po->TOLERANCE))
			if(str_intersect_r(arr_c[dsItem(genes,0)],arr_c[ki],'N') < floor(cnt * po->TOLERANCE))
			{
				candidates[ki] = FALSE;
				continue;
			}


			vector<int>a(0);
			get_Genes_LCS(arr_c[dsItem(genes,0)],arr_c[dsItem(genes,1)],lcs_tags[dsItem(genes,1)],a);
			
			discrete *gm=new discrete[cols];
       			 for(int gmi=0;gmi<a.size();gmi++)
       			 {
           			 gm[a[gmi]]=cols-gmi;
    
   			     }
        vector<int>(0).swap(a);
			get_Genes_LCS_R_with_lcs(gm,arr_c[ki],reve_tag,lcs_tags[dsItem(genes,1)]);						    		
			m_ct = intersect_row_LCS(colcand,reve_tag);
			delete [] gm;

//			get_Genes_LCS_R_with_lcs(arr_c[dsItem(genes,0)],arr_c[ki],reve_tag,lcs_tags[dsItem(genes,1)]);
//			m_ct = intersect_row_LCS(colcand,reve_tag);

			if (candidates[ki] && (m_ct >= floor(cnt * po->TOLERANCE)))//g//lxy
			//if (candidates[ki] && (m_ct >= floor(cnt * 0.85)-1))//g//lxy
			{
				int temp;
				for(temp=0;temp<cols;temp++)
				{
					if(colcand[temp])
					{	int tmpcount = colsStat[temp];
						if(reve_tag[temp]!=0)
							tmpcount++;
						if(tmpcount < floor(components * 0.1)-1)//g//lxy
						//if(tmpcount < floor(components * 0.1)-1)//lxy
						{	
							colChose = FALSE;
							break;
						}
					}
				}
				if(colChose == TRUE)			
				{	
					dsPush(genes,ki);
					
					//std::cout<<"new- ="<<ki<<std::endl;
					components++;
					candidates[ki] = FALSE;
					for(temp=0;temp<cols;temp++)
					{
						if(reve_tag[temp]!=0 && colcand[ki])
						{
							colsStat[temp]++;
						}
					}
				}
			}
			free(reve_tag);
			reve_tag=NULL;
		}
		/* save the current cluster*/
		//std::cout<<"######"<<std::endl;
		b_genes = dsNew(b->block_rows_pre);
		for (ki = 0; ki < b->block_rows_pre; ki++)
			dsPush(b_genes, dsItem(genes,ki));
		/* store gene arrays inside block */
		b->genes = dsNew(components);
		b->conds = dsNew(cols);
	
		scan_block(colcand, b);
		free(colcand);
		colcand=NULL;
		if (b->block_cols <po->COL_WIDTH || components < 5) continue;
		b->block_rows = components;
                if (po->IS_pvalue)
			b->score = -(100*log(b->pvalue));
		else
			b->score=-(b->block_rows*log(2*rows*cols*b->block_rows)-b->block_rows*b->block_cols*log(b->block_cols/2.7)-0.5*b->block_rows*log(6.28*b->block_cols));
;
		dsClear(b->genes);
		for ( ki=0; ki < components; ki++)
			dsPush(b->genes,dsItem(genes,ki));
		for(ki = 0; ki < components; ki++)
		{
			if(!isInStack(allincluster, dsItem(genes,ki))) 
				dsPush(allincluster,dsItem(genes,ki));
		}	
		/*save the current block b to the block list bb so that we can sort the blocks by their score*/
		bb[block_id++] = b;

		/* reaching the results number limit */
		if (block_id == po->SCH_BLOCK) break;
		verboseDot();	
	}
	/* writes character to the current position in the standard output (stdout) and advances the internal file position indicator to the next position.
	 * It is equivalent to putc(character,stdout).*/
	putchar('\n');
	/* free-up the candidate list */
	free(candidates);
	free(allincluster);
	free (pvalues);
	free(colsStat);
	return report_blocks(fw, bb, block_id);
}

/************************************************************************/
static void print_params(FILE *fw)
{
	char filedesc[LABEL_LEN];
	strcpy(filedesc, "continuous");
	if (po->IS_DISCRETE) 
		strcpy(filedesc, "discrete");
	fprintf(fw, "# LCCS version %.1f output\n", VER);
	fprintf(fw, "# Datafile %s: %s type\n", po->FN, filedesc);
	fprintf(fw, "# Parameters: -k %d -f %.2f -c %.2f -o %d",
			po->COL_WIDTH, po->FILTER, po->TOLERANCE, po->RPT_BLOCK);
	if (!po->IS_DISCRETE) 
		fprintf(fw, " -q %.2f -r %d", po->QUANTILE, po->DIVIDED);
	fprintf(fw, "\n\n");
}

/************************************************************************/
static int report_blocks(FILE* fw, Block** bb, int num)
{
	//print_params(fw);
	sort_block_list(bb, num);
	
	int i, j,k;
	/*MIN MAX et al functions can be accessed in struct.h*/
        int n = MIN(num, po->RPT_BLOCK);
	bool flag;

	Block **output = (Block **)malloc(sizeof(*output) * n);

	Block **bb_ptr = output;
	Block *b_ptr;
	double cur_rows, cur_cols;
	double inter_rows, inter_cols;
        /*double proportion;*/
	
	/* the major post-processing here, filter overlapping blocks*/
	i = 0; j = 0;
	while (i < num && j < n)
	{
		b_ptr = bb[i];
		cur_rows = b_ptr->block_rows;
		cur_cols = b_ptr->block_cols;
		
		

		flag = TRUE;
		k = 0;
		while (k < j)
		{
			inter_rows = dsIntersect(output[k]->genes, b_ptr->genes);
			inter_cols = dsIntersect(output[k]->conds, b_ptr->conds);
			
			if (inter_rows*inter_cols > po->FILTER*cur_rows*cur_cols)
			{
				flag = FALSE; 
				break;
			}
                        k++;
		}
	        i++;
		if (flag)
		{
			//print_bc(fw, b_ptr, j++);//原来�?
			test_print(fw,b_ptr,j++);
			*bb_ptr++ = b_ptr;
		}

		
	}
	return j;
}
/************************************************************************/

static int block_cmpr(const void *a, const void *b)
/* compare function for qsort, descending by score */
{
	return ((*(Block **)b)->score - (*(Block **)a)->score);
}

static void sort_block_list(Block **el, int n)
{
	qsort(el, n, sizeof *el, block_cmpr);
}
/************************************************************************/
long double get_pvalue (continuous a, int b)
{
	int i =0;
	long double one = 1.0, pvalue=0.0;
	long double poisson=one/exp(a);
	for (i=0;i<b;i++)
	{
			poisson=poisson*a/(i+1);
	}
	pvalue=poisson;
//		std::cout<<"pvalue= "<<scientific<<pvalue<<std::endl;
	return pvalue;
}
/*
long double get_pvalue(continuous a, int b)
{
	continuous	pvalue=(a-pow(b,0.5)*2)*1.0/pow(b,1.0/6);
	continuous p=pvalue-0.5;
	if (p<0){p=2;}
	return p;

}
*/
/**************************************************************************/
int get_Genes_LCS(const discrete *s1, const discrete *s2, char *lcs_tg,vector<int>&a) {
  vector<int> columns;
/*
  for (int j = 0; j < cols; j++) {
    if (s1[j] != 0 && s2[j] != 0) {
      columns.push_back(j);
    }
  }
  */
  for (int j = 0; j < cols; j++) {
  
    if (s1[j]!=0&&s2[j]!=0){
 
      columns.push_back(j);
      }
  }  

  if (columns.empty()) {
    return 0;
  }


  DagGraph dag(s1, s2, columns,po->IS_const);
  vector<vector<int>> lpList;//longest path list??
  dag.getAllLongestPaths(lpList);
  for (const auto &path : lpList) {
    for (auto cc : path) {
      //lcs_tg[cc] = 1;//感觉可以�?
	lcs_tg[columns[cc]] = 1;
	a.push_back(columns[cc]);
	//cout<<"cc "<<a.back()<<endl;
    }
    break;
  }
  vector<int>(0).swap(columns);  
	int lpsize=lpList[0].size();
	vector<vector<int> >(0).swap(lpList);
  	return lpsize;
}

/**************************************************************************/
int get_Genes_LCS_length(const discrete *s1, const discrete *s2)
{
  vector<int> columns;

  for (int j = 0; j < cols; j++) {
    if(s1[j]!=0&&s2[j]!=0)
      {
      columns.push_back(j);  
      }  
  }

	
  if (columns.empty()) {
    return 0;
  }

  DagGraph dag(s1, s2, columns,po->IS_const);
  vector<vector<int>> lpList;
  dag.getAllLongestPaths(lpList);
  vector<int>(0).swap(columns);
	int lpsize=lpList[0].size();
	vector<vector<int> >(0).swap(lpList);
  	return lpsize;
}

/**************************************************************************/
int get_Genes_LCS_with_lcs(const discrete *s1, const discrete *s2,char *lcs_tg,char *lcs_seed)
{
  vector<int> columns;

  for (int j = 0; j < cols; j++) {
    //if (lcs_seed[j]!=0) {     
    if (lcs_seed[j]==0) {  //2020-7-24   
      continue;
}  
    if(s1[j]!=0&&s2[j]!=0)
    {
columns.push_back(j); 
}  
  }

  if (columns.empty()) {
    return 0;
  }

  DagGraph dag(s1, s2, columns,po->IS_const);
  vector<vector<int>> lpList;
  dag.getAllLongestPaths(lpList);
  for (const auto &path : lpList) {
    for (auto cc : path) {
      //lcs_tg[cc] = 1;//
	lcs_tg[columns[cc]] = 1;
    }
    break;
  }
  vector<int>(0).swap(columns);
	int lpsize=lpList[0].size();
	vector<vector<int> >(0).swap(lpList);
  	return lpsize;
}

/**************************************************************************/
/*原来�?
int get_Genes_LCS_R_with_lcs(const discrete *s1,const discrete *s2,char *lcs_tagR,char *lcs_seed)
{
  vector<int> columns;

  for (int j = 0; j < cols; j++) {
    if (!lcs_seed[j]) {
      continue;
    }
    if (s1[j] != 0 && s2[j] != 0) {
      columns.push_back(j);
    }
  }

  if (columns.empty()) {
    return 0;
  }
  DagGraph dag(s1, s2, columns, true);
  vector<vector<int>> lpList;
  dag.getAllLongestPaths(lpList);
  for (const auto &path : lpList) {
    for (auto cc : path) {
      lcs_tagR[columns[cc]] = 1;
    }
    break;
  }
  return lpList[0].size();

}
*/
int get_Genes_LCS_R_with_lcs(const discrete *s1,const discrete *s2,char *lcs_tagR,char *lcs_seed)

{
	int rank=po ->DIVIDED;
	short int *s3=new  short int[cols];

  vector<int> columns;
  for (int j = 0; j < cols; j++) {

    if (!lcs_seed[j]) {
      continue;
    }
    if(s1[j]!=0&&s2[j]!=0)
    {
      columns.push_back(j);
      }
  }
  if (columns.empty()) 
  	{
    	return 0;
  	}
  if(po->dataMode==0 && po->QUANTILE < 0.5)
	{
	DagGraph dag(s1, s2, columns,po->IS_const);
  	vector<vector<int>> lpList;
  	dag.getAllLongestPaths(lpList);
  	for (const auto &path : lpList) 
  	{
    for (auto cc : path) {
      lcs_tagR[columns[cc]] = 1;
    }
    break;
    } 
  	delete[] s3;
  	vector<int>(0).swap(columns);
   	return lpList[0].size();
	}
  else
	{	for(int j=0;j<cols;j++)
	{
		s3[j]=rank-s2[j]+1;
	}
  	DagGraph dag(s1, s3, columns,po->IS_const);
  	vector<vector<int>> lpList;
  	dag.getAllLongestPaths(lpList);
  	for (const auto &path : lpList) {
  	  for (auto cc : path) {
   		   lcs_tagR[columns[cc]] = 1;
   	 }
   	 break;
  	}
 	 delete[] s3;
 	 vector<int>(0).swap(columns);
	int lpsize=lpList[0].size();
	vector<vector<int> >(0).swap(lpList);
  	return lpsize;
  }
 
}
