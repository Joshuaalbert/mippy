'''Does the interfacing to the c module'''

from mippy.interface.C_lpsolver_interface import lpsolve
import numpy as np

def submit_problem(problem_name, solver_type, num_rows, num_cols, maximize, id, lines):
    '''Does some type checking and submits LP job to glpk interface'''
    assert int(solver_type) in [4,5], "solver_type {} must be one of SIMP(4), MIP(5)".format(int(solver_type))
    assert isinstance(id,np.ndarray), "id must be an array"
    assert isinstance(lines,np.ndarray), "lines must be an array"
    assert id.dtype == int, "id must be integer array"
    assert lines.dtype == float, "lines must be float array"
    result = np.zeros(int(num_cols),dtype=np.double) 
    #print(id,lines)
    lpsolve(str(problem_name), int(solver_type), int(num_rows), int(num_cols), int(maximize), id.astype(np.double), lines.astype(np.double), result)
    return result
