def test_lpsolver():
    import numpy as np
    from ..lpsolver import LPSolver
    rows = 5
    cols = 5
    A_eq = np.random.normal(size=[rows,cols])
    b_eq = np.random.normal(size=cols)
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
    #l.submit_problem(prob)
