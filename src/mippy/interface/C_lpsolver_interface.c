#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_9_API_VERSION
#include <numpy/arrayobject.h>
#include "C_lpsolver_interface.h"
#include <stdlib.h>
#include <stdio.h>
#include <glpk.h>
#include <time.h>

#define _PY3_



/* ==== Set up the methods table ====================== */
static PyMethodDef C_lpsolver_interfaceMethods[] = {
	{"lpsolve", lpsolve, METH_VARARGS,"Solves the LP given"},
	{NULL, NULL,0,NULL}     /* Sentinel - marks the end of this structure */
};

#ifdef _PY3_
	/* python 3 way */
	static struct PyModuleDef C_lpsolver_interfaceModule =
	{
		PyModuleDef_HEAD_INIT,
		"C_lpsolver_interface", /* name of module */
		"Handles the interface with GLPK",          /* module documentation, may be NULL */
		-1,          /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
		C_lpsolver_interfaceMethods
	};

	PyMODINIT_FUNC PyInit_C_lpsolver_interface(void)
	{
		//printf("Loading C module\n");
		import_array();
		return PyModule_Create(&C_lpsolver_interfaceModule);
	}

//#else
	/* python 2 way */
/*
	// Module name must be _lpsolver_interface in compile and linked 
	void init_lpsolver_interface()  {
		(void) Py_InitModule("_lpsolver_interface", _lpsolver_interfaceMethods);
		import_array();  // Must be present for NumPy.  Called first after above line.
	}
*/
#endif

static PyObject *lpsolve (PyObject *self, PyObject *args)
{
	char *problem_name;
	int num_rows, num_cols,maximize, num_lines, solve_type;
	PyObject *id_obj, *lines_obj, *result_obj;
	PyObject *id_array, *lines_array, *result_array;
	double *id;
	double *lines, *result;
	/* Parse tuple (char *problem_name, int solve_type, int num_rows, int num_cols, int maximize, int *id, double *lines, double *result)*/
	if (!PyArg_ParseTuple(args, "siiiiOOO", &problem_name, &solve_type, &num_rows, &num_cols, &maximize,&id_obj, &lines_obj, &result_obj))  return NULL;
	id_array = PyArray_FROM_OTF(id_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
	lines_array = PyArray_FROM_OTF(lines_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
	result_array = PyArray_FROM_OTF(result_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);

	if (lines_array == NULL || id_array == NULL || result_array == NULL){
		printf("%ld %ld %ld\n",lines_array, id_array, result_array);
		Py_XDECREF(id_array);
		Py_XDECREF(lines_array);
		Py_XDECREF(result_array);
		printf("Failed to get arrays\n");
		return NULL;
	}

	num_lines = (int) PyArray_DIM(id_array,0);
	//printf("num lines %d\n",num_lines);
	//cid = (int *) malloc((size_t) num_lines*sizeof(int));
	id = (double *) PyArray_DATA(id_array);
	lines = (double *) PyArray_DATA(lines_array);
	result = (double *) PyArray_DATA(result_array);

	
	int row, col;
	double val,upper,lower,fixed,*ar;
	unsigned long count;
	int *ia, *ja, *coltype;
	glp_prob *lp;
	lp = glp_create_prob();
	glp_set_prob_name(lp,problem_name);
	
	printf ("A is %d x %d\n",num_rows,num_cols);

	glp_add_rows(lp,num_rows);
	glp_add_cols(lp,num_cols);
	ar = (double *)calloc(num_lines + 1,sizeof(double));
	ia = (int *)calloc(num_lines + 1,sizeof(int));
	ja = (int *)calloc(num_lines + 1,sizeof(int));
	coltype = (int *)calloc(num_cols + 1,sizeof(int));

	if (maximize){
		glp_set_obj_dir(lp,GLP_MAX); 
	} else {
		glp_set_obj_dir(lp,GLP_MIN);
	}
	//printf("reading lines\n");


	//for (int i = 0; i < num_lines; i++) printf("id = %d\n", (int)id[i]);
	count = 1;
	for (int i = 0; i < num_lines; i++){
		int line_idx = i*3;
		double slot1 = lines[line_idx];
		double slot2 = lines[line_idx + 1];
		double slot3 = lines[line_idx + 2];
		//printf("%lf %lf %lf\n",slot1, slot2, slot3);
		//printf("%d %d %d\n", (int) slot1, (int) slot2, (int) slot3);
		switch ((int)id[i])
		{	
			//Set constraints
			case A_LINE:// set constraint: row, col, val
				//printf("%d, %d %d %lf\n", (int)id[i], (int) slot1, (int) slot2, slot3);
				ia[count] = (int) slot1;
				ja[count] = (int) slot2;
				ar[count] = slot3;
				++count;
				break;
			//Row constraints
			case L_LINE:// set row lower bound: row, lower
				glp_set_row_bnds(lp,(int) slot1, GLP_LO, slot2, slot2);
				break;
			case U_LINE:// set row upper: row, upper
				glp_set_row_bnds(lp,(int) slot1 ,GLP_UP,slot2,slot2);
				break;
			case F_LINE:// set fixed: row, fixed
				glp_set_row_bnds(lp,(int) slot1,GLP_FX,slot2,slot2);
				break;
			case D_LINE:// set double bounds: row, lower, upper
				glp_set_row_bnds(lp,(int) slot1,GLP_DB,slot2,slot3);
				break;
			case R_LINE:// set free: row
				glp_set_row_bnds(lp,(int) slot1,GLP_FR,0.0,0.0);
				break;
			//Variable bounds
			case l_LINE:// set lower: col, lower
				glp_set_col_bnds(lp,(int) slot1,GLP_LO,slot2,slot2);
				break;
			case u_LINE:// set upper: col, upper
				glp_set_col_bnds(lp,(int) slot1,GLP_UP,slot2,slot2);
				break;
			case f_LINE:// set fixed: col, val
				glp_set_col_bnds(lp,(int) slot1,GLP_FX,slot2,slot2);
				break;
			case d_LINE:// set double: col, lower, upper
				glp_set_col_bnds(lp, (int)slot1, GLP_DB, slot2, slot3);
				break;
			case r_LINE:// set free: col
				glp_set_col_bnds(lp,(int) slot1,GLP_FR,0.0,0.0);
				break;
			//Set objective value
			case C_LINE:// set obj: col, val
				glp_set_obj_coef(lp,(int) slot1, slot2);
				break;
			//Variable types
			case b_LINE:// set var binary: col
				glp_set_col_kind(lp,(int) slot1,GLP_BV);
				//glp_set_col_bnds(lp,(int) slot1,GLP_DB,0,1);
				coltype[col] = GLP_BV;
				break;
			case i_LINE:// set var integer: col
				glp_set_col_kind(lp,(int) slot1,GLP_IV);
				coltype[col] = GLP_IV;
				break;
			case c_LINE:// set var continuous: col
				glp_set_col_kind(lp,(int) slot1,GLP_CV);
				coltype[col] = GLP_CV;
				break;

			default:
				printf("Failed to parse: %d\n",i);
		}
	}
	glp_iocp iocp;
	glp_smcp smcp;
	glp_init_smcp(&smcp);
	smcp.presolve=GLP_ON;
	smcp.msg_lev=GLP_MSG_ERR;
	smcp.meth=GLP_DUALP;
	glp_init_iocp(&iocp);
	iocp.presolve=GLP_ON;
	iocp.msg_lev=GLP_MSG_ALL;
	glp_load_matrix(lp,count-1,ia,ja,ar);
	char mod_file[50];
	sprintf(mod_file,"%s.mps",problem_name);
	printf("MPS model in: %s\n", mod_file);
	//glp_write_mps(lp,GLP_MPS_FILE,NULL,mod_file);
		switch (solve_type)
	{
		case SIMP:
			glp_simplex(lp,&smcp);
			for (int i=1;i<=num_cols;++i)
			{
				result[i - 1] = glp_get_col_prim(lp,i);
			}
			break;
		case MIP:
			glp_intopt(lp,&iocp);
			
			for (int i=1;i<=num_cols;++i)
			{
				result[i-1] = glp_mip_col_val(lp,i);
			}
			break;
				default:
			printf("Not Solving with solve_type: %d\n",solve_type);
	}
	
	free(ia),free(ja),free(ar);free(coltype);
	Py_DECREF(id_array);
	Py_DECREF(lines_array);
	Py_DECREF(result_array);
	glp_delete_prob(lp);

	return Py_BuildValue("i", 1);
}
