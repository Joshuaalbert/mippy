def test_lpsolver():
    import numpy as np
    from ..lpsolver import LPSolver
    rows = 5
    cols = 5
    A_eq = np.random.normal(size=[rows,cols])
    b_eq = np.random.normal(size=rows)
    sol = np.linalg.inv(A_eq).dot(b_eq)
    maximize = True
    problem_name = "TEST"
    c_obj = np.zeros(cols)
    l = LPSolver(c_obj,A_eq = A_eq, b_eq = b_eq, maximize=maximize, problem_name=problem_name, solver_type='SIMP')
    l.compile()
    result = l.submit_problem()
    print("Comparing solutions {} ? {}".format(result,sol))
    assert np.allclose(result,sol)

    A_lt = np.array([[1,2],[3,1]])
    b_lt = np.array([1,2])
    c_obj = np.array([0.6,0.5])
    l = LPSolver(c_obj,A_lt=A_lt, b_lt = b_lt, maximize=True, problem_name="TEST2", solver_type='SIMP')
    #l.set_variable_type(0,'i',('<>',0,1))
    #l.set_variable_type(1,'i',('<>',0,1))
    l.compile()
    result = l.submit_problem()
    sol = np.array([0.6, 0.2])
    print("Comparing solutions {} ? {}".format(result,sol))
    assert np.allclose(result,sol)


    rows = 2
    cols = 2
    A_gt = np.array([[1,2,3]])
    b_gt = np.array([1])
    A_lt = np.array([[2,3,1]])
    b_lt = np.array([0])
    maximize = True
    problem_name = "TEST3"
    c_obj = np.array([0,0,1])
    l = LPSolver(c_obj,A_lt=A_lt, b_lt=b_lt,A_gt= A_gt, b_gt=b_gt,  maximize=maximize, problem_name=problem_name, solver_type='MIP')
    l.set_variable_type(0,'i',('<>',-1,1))
    l.set_variable_type(1,'i',('<>',-2,2))
    l.set_variable_type(2,'i',('<>',-3,3))
    l.compile()
    result = l.submit_problem()
    #print(result)
