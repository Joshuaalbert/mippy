'''
Interface to GNU linear programming kit (GLPK) solver such that the user
can stick to using familiar numpy arrays to define constraints.
Deviates from scipy which doesn't allow MIP (mixed integer programming)

author: Joshua G. Albert, albert@strw.leidenuniv.nl
'''

from mippy.interface.lpsolver_interface import submit_problem

import numpy as np

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
    def __init__(self, c_obj, A_eq = None, b_eq = None, A_lt = None, b_lt = None, A_gt = None, b_gt = None,  A_db = None, b_db = None, maximize=True, problem_name="Mippy_problem", solver_type='SIMP'):
        self.solver_type = solver_type
        self.problem_name = problem_name
        self.maximize = maximize
        self.num_variables = int(np.size(c_obj))
        self.objective = c_obj
        self.variable_map = {'c':14,'b':12,'i':13,'<':7,'>':6,'=':8,'<>':9,'*':10}
        self.constraint_map = {'>':1,'<':2,'=':3,'<>':4,'*':5}

        
        #by default all continuous variables
        self.variable_types = [14]*self.num_variables#cont
        self.variable_bounds = [(10,0.,0.)]*self.num_variables#unbounded
        
        self.A = []
        self.constraints = []
        self.num_constraints = 0
        if A_eq is not None:
            assert b_eq is not None, "b_eq must be given"
            A_eq = np.array(A_eq)
            b_eq = np.array(b_eq)
            i = 0
            while i < A_eq.shape[0]:
                self.add_constraint(A_eq[i,:],('=',b_eq[i],0.))
                i += 1
        if A_lt is not None:
            assert b_lt is not None, "b_lt must be given"
            A_lt = np.array(A_lt)
            b_lt = np.array(b_lt)
            i = 0
            while i < A_lt.shape[0]:
                self.add_constraint(A_lt[i,:],('<',b_lt[i],0.))
                i += 1
        if A_gt is not None:
            assert b_gt is not None, "b_gt must be given"
            A_gt = np.array(A_gt)
            b_gt = np.array(b_gt)
            i = 0
            while i < A_gt.shape[0]:
                self.add_constraint(A_gt[i,:],('>',b_gt[i],0.))
                i += 1
        if A_db is not None:
            assert b_db is not None, "b_db must be given"
            A_db = np.array(A_db)
            b_db = np.array(b_db)
            assert b_db.shape[1] == 2, "double bounds require second dimension be 2 {}".format(b_db.shape)
            i = 0
            while i < A_db.shape[0]:
                self.add_constraint(A_db[i,:],('<>',b_db[i,0],b_db[i,1]))
                i += 1
        self.id = None
        self.lines = None
        

    @property
    def id(self):
        assert self._id is not None, "You must compile the problem first"
        return self._id
    @id.setter
    def id(self,val):
        self._id = val
    @property
    def lines(self):
        assert self._lines is not None, "You must compile the problem first"
        return self._lines
    @lines.setter
    def lines(self,val):
        self._lines = val    
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
        if val is None:
            self._num_variables = None
            return
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
        (1, lower_bound, 0.),  v > lower_bound
        (2, upper_bound, 0.),  v < upper_bound
        (3, fixed_equality, 0.), v = fixed_equality
        (4, lower_bound, upper_bound), lower_bound < v < upper_bound
        (5,0,0), free variable no bounds
        """
        assert self.num_variables == len(row), "row size ({}) not same a variable count ({})".format(len(row),self.num_variables)
        assert not np.any(np.isnan(row)), "NaNs in constraint"
        self.A.append(row)
        assert constraint[0] in ['<','>','=','<>','*'], "equality type ({}) not one of '<','>','=','<>','*'".format(constraint[0])
        self.constraints.append((self.constraint_map[constraint[0]],*constraint[1:]))
        self.num_constraints += 1

    def set_objective_col(self,col, val):
        """Set the objective col to given value"""
        self.objective[col] = val
    def set_variable_type(self,col,type,bounds=None):
        """Set variable at col to one of the following with given bounds,
        'c' 14 - continuous
        'b' 12 - binary 
        'i' 13 - integer
        `bounds` if not None should be defined as a tuple, one of:
        ('>' 6, lower_bound, 0.),  v > lower_bound
        ('<' 7, upper_bound, 0.),  v < upper_bound
        ('=' 8, fixed_equality, 0.), v = fixed_equality
        ('<>' 9, lower_bound, upper_bound), lower_bound < v < upper_bound
        ('*' 10, 0., 0,.), free variable no bounds [default]
        """
        assert type in ['c','b','i'], "variable type ({}) not one of 'b','c','i'".format(type)
        assert col < self.num_variables, "index out of range"
        self.variable_types[col] = self.variable_map[type]
        if bounds is not None:
            assert bounds[0] in ['>','<','=','<>','*'], "bound type ({}) not understood".format(bounds[0])
            if bounds[0] in ['>', '<', '=']:
                self.variable_bounds[col] = (self.variable_map[bounds[0]], bounds[1],0.)
            if bounds[0] == '<>':
                assert len(bounds) == 3, "bounds invalid {}".format(bounds)
                self.variable_bounds[col] = (self.variable_map[bounds[0]], *bounds[1:])
            if bounds[0] == '*':
                self.variable_bounds[col] = (self.variable_map[bounds[0]], 0.,0.)

    def compile(self):
        """Compile the program into a file that can be submitted.
        Returns the filename which can be used by `submit_problem`"""
        id = []
        lines = []
        line_idx = 0
        #declare variable types and bounds
        col = 0
        while col < self.num_variables:
            #variable types
            id.append(self.variable_types[col])
            lines.append((col + 1,0.,0.))
            line_idx += 1
            #variable bounds
            id.append(self.variable_bounds[col][0])
            lines.append((col+1,*self.variable_bounds[col][1:]))
            line_idx += 1
            #constraints
            if self.objective[col] != 0.:
                id.append(11)
                lines.append((col+1, self.objective[col],0.))
                line_idx += 1
                #declare constraint matrix and bounds
            col += 1
        row = 0
        while row < self.num_constraints:
            id.append(self.constraints[row][0])
            lines.append((row+1,*self.constraints[row][1:]))
            line_idx += 1
            constraint_row = self.A[row]
            col = 0
            while col < self.num_variables:
                if constraint_row[col] != 0.:
                    id.append(0)
                    lines.append((row+1, col+1, constraint_row[col]))
                    line_idx += 1
                col += 1
            row += 1
        self.id = np.array(id,dtype=int)
        #print(lines)
        self.lines = np.array(lines,dtype=np.double)
        assert not np.any(np.isnan(self.lines)), "NaNs in constraint"

    def submit_problem(self):
        '''Submit the problem that was compiled into mippy_file.
        Return the solution of each variable.'''
        if self.solver_type == 'SIMP':
            results = submit_problem(self.problem_name, 4, self.num_constraints, self.num_variables, self.maximize, self.id, self.lines)
        if self.solver_type == 'MIP':
            results = submit_problem(self.problem_name, 5, self.num_constraints, self.num_variables, self.maximize, self.id, self.lines)
        return results
