def test_lpsolver():
    import numpy as np
    from ..lpsolver import LPSolver
    rows = 5
    cols = 5
    A_eq = np.random.normal(size=[rows,cols])
    b_eq = np.random.normal(size=rows)
    sol = np.linalg.inv(A_eq).dot(b_eq)
    A_lt = None
    b_lt = None
    A_gt = None
    b_gt = None
    maximize = True
    problem_name = "TEST"
    c_obj = np.zeros(cols)
    l = LPSolver(A_eq, b_eq, A_lt, b_lt, A_gt, b_gt, c_obj, maximize, problem_name, 'SIMP')
    prob = l.compile()
    result = l.submit_problem(prob)
    print("Comparing solutions {} ? {}".format(result,sol))
    assert np.allclose(result,sol)

    rows = 2
    cols = 2
    A_eq = None
    b_eq = None
    A_gt = np.array([[1,2,3]])
    b_gt = np.array([1])
    A_lt = np.array([[2,3,1]])
    b_lt = np.array([0])
    maximize = True
    problem_name = "TEST"
    c_obj = np.array([0,0,1])
    l = LPSolver(A_eq, b_eq, A_lt, b_lt, A_gt, b_gt, c_obj, maximize, problem_name, 'MIP')
    l.set_variable_type(0,'i',('<>',-1,1))
    l.set_variable_type(1,'i',('<>',-2,2))
    l.set_variable_type(2,'i',('<>',-3,3))
    prob = l.compile()
    result = l.submit_problem(prob)
    print(result)
