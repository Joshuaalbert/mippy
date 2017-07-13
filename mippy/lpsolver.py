'''
Interface to linear programming solver such that the user
can stick to using familiar numpy arrays to define constraints.
Deviates from scipy which doesn't allow MIP (mixed integer programming)

author: Joshua G. Albert, albert@strw.leidenuniv.nl
'''

import mippy.interface

import numpy as np
import os

class LPSolver(object):
    """
    Class to interface standard numpy arrays to glpk solver and get the results.

    Available solvers are:
    'SIMPLEX' - Simplex algorithm
    'MIP' - Mixed Integer Programming using Branch Cut

    Given a set of matrix equalities and inequalities of the forms:

    A_eq.v = b_eq, A_lt.v < b_lt, A_gt.v > b_gt
    
    and and objective function to maximize ([True]/False), 
    
    f = c_obj.v
    
    assemble and Simplex problem with,

    l = LPSolver(A_eq, b_eq, A_lt, b_lt, A_gt, b_gt, c_obj, maximize, problem_name, 'SIMP')

    If you have mor constraints to add then read the docs on functions `add_constraints`.
    By default your variables are unconstrained. You can change this using `set_variable_type`.
    
    Then compile the problem which will create a mippy problem file using the problem_name.

    problem_file = l.compile()
    print("Name of the mippy problem file is {}".format(problem_file))

    We can then solve the problem by calling,

    """
    def __init__(self, A_eq, b_eq, A_lt, b_lt, A_gt, b_gt, c_obj, maximize=True, problem_name="Mippy_problem", solver_type='SIMP'):
        self.solver_type = solver_type
        self.problem_name = problem_name
        self.maximize = maximize
        self.num_variables = int(np.size(c_obj))
        #by default all continuous variables
        self.variable_types = ['c']*self.num_variables
        self.variable_bounds = [('*',)]*self.num_variables
        self.objective = c_obj
        self.A = []
        self.contraints = []
        self.num_constraints = 0
        if A_eq is not None:
            i = 0
            while i < A_eq.shape[0]:
                self.add_constraint(A_eq[i,:],("=",b_eq))
                i += 1
        if A_lt is not None:
            i = 0
            while i < A_lt.shape[0]:
                self.add_constraint(A_lt[i,:],("<",b_lt))
                i += 1
        if A_gt is not None:
            i = 0
            while i < A_gt.shape[0]:
                self.add_constraint(A_gt[i,:],(">",b_gt))
                i += 1
        

    @property
    def maximize(self):
        return self._maximize
    @maximize.setter
    def maximize(self,val):
        assert isinstance(val,bool), "maximize should be True or False"
        self._maximize = val
    @property
    def num_variables(self):
        return self._num_variables

    @num_variables.setter
    def num_variables(self,val):
        assert val > 0, "Setting invalid number of vaiables ({})".format(val)
        self._num_variables = int(val)
        self.objective = np.zeros(self.num_variables, dtype=np.double)

    @property
    def objective(self):
        assert self._objective is not None, "objective not allocated"
        return self._objective
    @objective.setter
    def objective(self,objective):
        assert np.size(objective) == self.num_variables
        self._objective = objective

    def add_constraint(self,row,constraint):
        """Add a constraint with row and constraint with corresponding lines in A.x ? b
        constraint is a tuple with one of:
        ('>', lower_bound),  v > lower_bound
        ('<', upper_bound),  v < upper_bound
        ('=', fixed_equality), v = fixed_equality
        ('<>', lower_bound, upper_bound), lower_bound < v < upper_bound
        ('*',), free variable no bounds
        """
        assert self.num_variables == len(row), "row size ({}) not same a variable count ({})".format(len(row),self.num_variables)
        assert not np.any(np.isnan(row)) and  not np.isnan(rhs), "NaNs in constraint"
        self.A.append(row)
        assert constraint[0] in ["<","=",">","<>","*"], "equality type ({}) not one of '<', '=', '>', '<>', or '*'".format(type)
        self.constraints.append(constraint)
        self.num_constraints += 1

    def set_objective_col(self,col, val):
        """Set the objective col to given value"""
        self.objective[col] = val

    def set_variable_type(self,col,type,bounds=('*',)):
        """Set variable at col to one of the following with given bounds,
        'c' - continuous
        'b' - binary 
        'i' - integer
        `bounds` if not None should be defined as a tuple, one of:
        ('>', lower_bound),  v > lower_bound
        ('<', upper_bound),  v < upper_bound
        ('=', fixed_equality), v = fixed_equality
        ('<>', lower_bound, upper_bound), lower_bound < v < upper_bound
        ('*',), free variable no bounds [default]
        """
        assert type in ['c','b','i'], "variable type ({}) not one of 'c', 'b', or 'i'"
        assert col < self.num_variables, "index out of range"
        self.variable_types[col] = type
        if bounds is not None:
            assert bounds[0] in ['>','<','=','<>','*'], "bound type ({}) not understood".format(bounds[0])
            self.variable_bounds[col] = bounds

    def compile(self):
        mippy_file = "{}.mippy".format(self.problem_name)
        f = open(mippy_file,"w+")
        # max or mix problem
        if self.maximize:
            f.write("M {:d} {:d}\n".format(self.num_constraints,self.num_variables)) #set minimization and nrows, ncols
        else:
            f.write("N {:d} {:d}\n".format(self.num_constraints,self.num_variables)) #set minimization and nrows, ncols
        #declare variable types and bounds
        col = 0
        while col < self.num_variables:
            f.write("{:s} {:d}\n".format(self.variable_types[col], col))
            if self.variable_bounds[col][0] == '<':
                f.write("u {:d} {:.15f}\n".format(col, self.variable_bounds[col][1]))
            if self.variable_bounds[col][0] == '>':
                f.write("l {:d} {:.15f}\n".format(col, self.variable_bounds[col][1]))
            if self.variable_bounds[col][0] == '=':
                f.write("f {:d} {:.15f}\n".format(col, self.variable_bounds[col][1]))
            if self.variable_bounds[col][0] == '<>':
                f.write("d {:d} {:.15f} {:.15f}\n".format(col, self.variable_bounds[col][1],self.variable_bounds[col][2]))
            if self.variable_bounds[col][0] == '*':
                f.write("r {:d}\n".format(col))
            col += 1
        #declare objective function
        col = 0
        while col < self.num_variables:
            f.write("C {:d} {:.15f}\n".format(col, self.objective[col]))
            col += 1
        #declare constraint matrix and bounds
        row = 0
        while row < self.num_constraints:
            #set row constraint
            if self.constraints[row][0] == '<':
                f.write("U {:d} {:.15f}\n".format(col, self.constraints[row][1]))
            if self.constraints[row][0] == '>':
                f.write("L {:d} {:.15f}\n".format(col, self.constraints[row][1]))
            if self.constraints[row][0] == '=':
                f.write("F {:d} {:.15f}\n".format(col, self.constraints[row][1]))
            if self.constraints[row][0] == '<>':
                f.write("D {:d} {:.15f} {:.15f}\n".format(col, self.constraints[row][1],self.constraints[row][2]))
            if self.constraints[row][0] == '*':
                f.write("R {:d}\n".format(col))
            # set the row
            col = 0
            while col < self.num_variables:
                f.write("A {:d} {:d} {:.15f}\n".format(row,col,self.A[row][col]))
                col += 1
        f.close()
        return mippy_file

    def submit_problem(self,mippy_file):
        '''Submit the problem that was compiled into mippy_file.
        Return the solution of each variable.'''
        if self.solver_type == 'SIMP':
            os.system("{:s}/lpsolver_interface {:s}  {:s} {:d}".format(mippy.interface.__file__, mippy_file, self.problem_name, 4))
        if self.solver_type == 'MIP':
            os.system("{:s}/lpsolver_interface {:s}  {:s} {:d}".format(mippy.interface.__file__, mippy_file, self.problem_name, 5))

if __init__ == '__main__':
    import numpy as np
    rows = 4
    cols = 5
    A_eq = np.random.normal(size=[rows,cols])
    b_eq = np.random.normal(size=cols)
    A_lt = None
    b_lt = None
    A_gt = None
    b_gt = None
    maximize = True
    problem_name = "TEST"
    c_obj = np.zeros(cols)
    l = LPSolver(A_eq, b_eq, A_lt, b_lt, A_gt, b_gt, c_obj, maximize, problem_name, 'SIMP')

