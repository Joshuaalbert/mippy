#include <stdlib.h>
#include <stdio.h>
#include <glpk.h>
#include <time.h>

#define NO_FILE 1
#define NO_ARGS 2
#define FAIL_MOD 3
#define SIMP 4
#define MIP 5
#define INTERIOR 6

int make_mps (int argc, char **argv)
{
	FILE *file;
	char str[256];
	if (argc == 1)
	{
		printf("Usage: %s <model from mod_make.py> <model_name> <solution_type>\n",argv[0]);
		printf("solution_type:\n %d - Simplex\n %d - Mix Integer\n %d - Interior Branch Cut (dysfunct)\n", SIMP, MIP, INTERIOR);
		return NO_ARGS;
	}
	file = fopen(argv[1],"r");
	if (file == NULL)
	{
		printf("%s does not exist\n",argv[1]);
		return NO_FILE;
	}
	int solve_type;
	sscanf(argv[3],"%d",&solve_type);
	int nrows, ncols, row, col;
	double val,upper,lower,fixed,*ar;
	unsigned long count;
	int *ia, *ja, *coltype;
	glp_prob *lp;
	lp = glp_create_prob();
	glp_set_prob_name(lp,argv[2]);
	count = 1;
	while (fgets(str,256,file))
	{
		switch (str[0])
		{
			//Case of maximization
			case 'M':// number of rows, cols
				sscanf(str,"M %d %d",&nrows,&ncols);
				printf ("A is %d x %d\n",nrows,ncols);
				glp_set_obj_dir(lp,GLP_MAX); 
				glp_add_rows(lp,nrows);
				glp_add_cols(lp,ncols);
				ar = (double *)calloc(nrows*ncols + 1,sizeof(double));
				ia = (int *)calloc(nrows*ncols + 1,sizeof(int));
				ja = (int *)calloc(nrows*ncols + 1,sizeof(int));
				coltype = (int *)calloc(ncols + 1,sizeof(int));
				break;
			//Case of minimization
			case 'N':// number of rows, cols
				sscanf(str,"N %d %d",&nrows,&ncols);
				printf ("A is %d x %d\n",nrows,ncols);
				glp_set_obj_dir(lp,GLP_MIN);
				glp_add_rows(lp,nrows);
				glp_add_cols(lp,ncols);
				ar = (double *)calloc(nrows*ncols + 1,sizeof(double));
				ia = (int*)calloc(nrows*ncols + 1,sizeof(int));
				ja = (int *)calloc(nrows*ncols + 1,sizeof(int));
				coltype = (int *) calloc(ncols + 1,sizeof(int));
				break;
			//Set constraints
			case 'A':// set constraint: row, col, val
				sscanf(str,"A %d %d %lf",&row,&col,&val);
				ia[count]=row;
				ja[count]=col;
				ar[count]=val;
				++count;
				break;
			//Row constraints
			case 'L':// set row lower bound: row, lower
				sscanf(str,"L %d %lf",&row,&lower);
				glp_set_row_bnds(lp,row,GLP_LO,lower,lower);
				break;
			case 'U':// set row upper: row, upper
				sscanf(str,"U %d %lf",&row,&upper);
				glp_set_row_bnds(lp,row,GLP_UP,upper,upper);
				break;
			case 'F':// set fixed: row, fixed
				sscanf(str,"F %d %lf",&row,&fixed);
				glp_set_row_bnds(lp,row,GLP_FX,fixed,fixed);
				break;
			case 'D':// set double bounds: row, lower, upper
				sscanf(str,"D %d %lf %lf",&row,&lower,&upper);
				glp_set_row_bnds(lp,row,GLP_DB,lower,upper);
				break;
			case 'R':// set free: row
				sscanf(str,"R %d",&row);
				glp_set_row_bnds(lp,row,GLP_FR,0.0,0.0);
				break;
			//Variable bounds
			case 'l':// set lower: col, lower
				sscanf(str,"l %d %lf",&col,&lower);
				glp_set_col_bnds(lp,col,GLP_LO,lower,lower);
				break;
			case 'u':// set upper: col, upper
				sscanf(str,"u %d %lf",&col,&upper);
				glp_set_col_bnds(lp,col,GLP_UP,upper,upper);
				break;
			case 'f':// set fixed: col, val
				sscanf(str,"f %d %lf",&col,&fixed);
				glp_set_col_bnds(lp,col,GLP_FX,fixed,fixed);
				break;
			case 'd':// set double: col, lower, upper
				sscanf(str,"d %d %lf %lf",&col,&lower,&upper);
				glp_set_col_bnds(lp,col,GLP_DB,lower,upper);
				break;
			case 'r':// set free: col
				sscanf(str,"r %d",&col);
				glp_set_col_bnds(lp,col,GLP_FR,0.0,0.0);
				break;
			//Set objective value
			case 'C':// set obj: col, val
				sscanf(str,"C %d %lf",&col,&val);	
				glp_set_obj_coef(lp,col,val);
				break;
			//Variable types
			case 'b':// set var binary: col
				sscanf(str,"b %d",&col);
			//	printf("%d is BV\n",col);
				glp_set_col_kind(lp,col,GLP_IV);
				glp_set_col_bnds(lp,col,GLP_DB,0,1);
				coltype[col] = GLP_BV;
				break;
			case 'i':// set var integer: col
				sscanf(str,"i %d",&col);
			//	printf("%d in IV\n",col);
				glp_set_col_kind(lp,col,GLP_IV);
				coltype[col] = GLP_IV;
				break;
			case 'c':// set var continuous: col
				sscanf(str,"c %d",&col);
				glp_set_col_kind(lp,col,GLP_CV);
				coltype[col] = GLP_CV;
				break;

			default:
				printf("Failed to parse: %s\n",str);
		}
	}
	glp_iocp iocp;
	glp_smcp smcp;
	//glp_iptcp iptcp;
	//glp_init_iptcp(&iptcp);
	//iptcp.msg_lev = GLP_MSG_ERR;
	glp_init_smcp(&smcp);
	smcp.presolve=GLP_ON;
	smcp.msg_lev=GLP_MSG_ERR;
	smcp.meth=GLP_DUALP;
	glp_init_iocp(&iocp);
	iocp.presolve=GLP_ON;
	iocp.msg_lev=GLP_MSG_ALL;
	glp_load_matrix(lp,count-1,ia,ja,ar);
	char mod_file[50];
	sprintf(mod_file,"%s.mps",argv[2]);
	printf("MPS model in: %s\n", mod_file);
	glp_write_mps(lp,GLP_MPS_FILE,NULL,mod_file);
	int glpk_solve;
	double x;
	//int solve_type = SIMP;
	char mod_results[50];
	glpk_solve = 1;
	switch (solve_type)
	{
		case SIMP:
			glp_simplex(lp,&smcp);
		//	exact arithmetic is slow
		//	glp_exact(lp,&smcp);
		//	printf("Simplex Obj: %f\n",glp_get_obj_val(lp));
			
			sprintf(mod_results,"%s.simplex.lpsol",argv[2]);
			file=fopen(mod_results,"w+");
			for (int i=1;i<=ncols;++i)
			{
				fprintf(file,"%d %.17e\n",i,glp_get_col_prim(lp,i));
			}
			fclose(file);
			printf("Solution in: %s\n", mod_results);
			break;
		case MIP:
			glp_intopt(lp,&iocp);
		//	printf("Obj: %f\n",glp_mip_obj_val(lp));
			
			sprintf(mod_results,"%s.mip.lpsol",argv[2]);
			file=fopen(mod_results,"w+");
			for (int i=1;i<=ncols;++i)
			{
				switch (glp_get_col_kind(lp,i))
				{
					case GLP_CV:
						fprintf(file,"%d %.17e\n",i,(double) glp_mip_col_val(lp,i));	
						break;
					case GLP_IV:
						fprintf(file,"%d %d\n",i,(int) glp_mip_col_val(lp,i));
						break;
					case GLP_BV:
						fprintf(file,"%d %d\n",i,(int) glp_mip_col_val(lp,i));
						break;
					default:
						printf("No type for col %d\n",i);
				}
			}
			fclose(file);
			printf("Solution in: %s\n", mod_results);
			break;
		/*
		case INTERIOR:
			glp_interior(lp,&iptcp);
			
			sprintf(mod_results,"%s.interior.lpsol",argv[2]);
			file=fopen(mod_results,"w+");
			for (int i=1;i<=ncols;++i)
			{
				fprintf(file,"%d %.17e\n",i,glp_ipt_col_prim(lp,i));
			}
			fclose(file);
			printf("Solution in: %s", mod_results);
			break;
			*/
		default:
			printf("Not Solving with solve_type: %d\n",solve_type);
	}
	//	glp_interior(lp,NULL);
	
	free(ia),free(ja),free(ar);free(coltype);
	glp_delete_prob(lp);
	return 0;
}

int  main(int argc, char **argv)
{
	return make_mps (argc,argv);
}
