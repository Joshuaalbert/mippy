from ..lpsolver import LPSolver

def test_scheduler():
    import numpy as np
    import astropy.units as au
    import astropy.coordinates as ac
    import astropy.time as at
    #Ni pointings, Nj observation sessions, Nk 1hour time slots per observation
    
    #num of hours each pointing should be obserbed
    Nh = 4
    obs_loc = ac.EarthLocation(6.86889*au.deg,52.90889*au.deg, 0*au.m)
    Nk = 3
    offsets = np.arange(Nk)#hours of obs

    Nj = 4
    #time of first obs0 timeslot0
    obs0 = at.Time("2017-07-25T00:00:00.000",format='isot').gps
    # lst
    lst = ac.AltAz(alt=90*au.deg,az=0*au.deg,location=obs_loc,obstime=at.Time(obs0,format='gps')).transform_to(ac.ICRS)
    Ni = 3
    coords = ac.SkyCoord((lst.ra.deg+np.random.uniform(low=-10,high=10,size=Ni))*au.deg,
                            (lst.dec.deg + np.random.uniform(low=-10,high=10,size=Ni))*au.deg,frame='icrs')
    #separated by 1 day
    observations = at.Time([obs0 + 24*3600*j for j in range(Nj)], format='gps')
    
    e_ijk = np.zeros([Ni,Nj,Nk])
    hourangle_ijk = np.zeros([Ni,Nj,Nk])
   

    for k in range(Nk):
        for j in range(Nj):
            time = at.Time(observations[j].gps + offsets[k]*3600.,format='gps')
            altaz = ac.AltAz(location=obs_loc,obstime=time)
            pointings = coords.transform_to(altaz)
            lst = ac.AltAz(alt=90*au.deg,az=0*au.deg,location=obs_loc,obstime=time).transform_to(ac.ICRS).ra
            for i in range(Ni):
                e_ijk[i,j,k] = pointings[i].alt.deg
                hourangle_ijk[i,j,k] =  (lst - coords[i].ra).to(au.deg).value
    #import pylab as plt
    #plt.hist(e_ijk.flatten())
    #plt.show()

    def h(i,j,k):
        """The index for ith pointing, kth day, jth timeslot"""
        return i + Ni*(j + Nj * k)

    #v_ijk = 1 if choosing pointing i in jk window else 0
    #Constraints
    A_eq = []
    b_eq = []
    A_lt = []
    b_lt = []
    A_gt = []
    b_gt = []
    # sum_i v_ijk <= 1 for all j,k (one pointing per slot)
    for j_ in range(Nj):
        for k_ in range(Nk): 
            row = np.zeros(Ni*Nj*Nk)
            for i in range(Ni):
                idx = h(i,j_,k_)
                row[idx] = 1.
            A_lt.append(row)
            b_lt.append(1.)
    
    # sum_j,k v_ijk = Nh for all i (require Nh of each pointing)
    for i_ in range(Ni):
        row = np.zeros(Ni*Nj*Nk)
        for k in range(Nk): 
            for j in range(Nj):
                idx = h(i_,j,k)
                row[idx] = 1.
        A_eq.append(row)
        b_eq.append(Nh)

    # sum_k v_ijk <= 1 for all i,j (one pointing per day or none)
    for i_ in range(Ni):
        for j_ in range(Nj):
            row = np.zeros(Ni*Nj*Nk)
            for k in range(Nk):
                idx = h(i_,j_,k)
                row[idx] = 1.
            A_lt.append(row)
            b_lt.append(1.)
    # eijk vijk > 30 #will fail because should only test iff vijk==1 (so add a conditional)
    for i_ in range(Ni):
        for j_ in range(Nj):
            for k_ in range(Nk):
                row = np.zeros(Ni*Nj*Nk)
                idx = h(i_,j_,k_)
                row[idx] = e_ijk[i_,j_,k_]
                A_gt.append(row)
                b_gt.append(30)
    # maximize the elevation of the points c_obj = e_ijk so max e_ijk v_ijk
    c_obj = np.zeros(Ni*Nj*Nk)
    for i in range(Ni):
        for j in range(Nj):
            for k in range(Nk):
                idx = h(i,j,k)
                c_obj[idx] = e_ijk[i,j,k]
    maximize=True
    problem_name="obs_schdlr"
    A_eq = np.array(A_eq)
    b_eq = np.array(b_eq)
    A_lt = np.array(A_lt)
    b_lt = np.array(b_lt)
    A_gt = None#np.array(A_gt)
    b_gt = None#np.array(b_gt)

    l = LPSolver(A_eq, b_eq, A_lt, b_lt, A_gt, b_gt, c_obj, maximize, problem_name, 'MIP')
    for i in range(Ni):
        for j in range(Nj):
            for k in range(Nk):
                idx = h(i,j,k)
                l.set_variable_type(idx,'i',('<>',0.,1.))
    prob = l.compile()
    result = l.submit_problem(prob)

    # Printing of results and testing validity
    print("Result is:".format(result))
    num_obs = np.zeros(Ni) 
    # each obs is a block
    for j in range(Nj):
        print("---- Obs {} ----".format(j))
        s = "{:<{width}}".format("",width=4)
        #x axis is pointing
        for i in range(Ni):
            s = "{} | {:<{width}}".format(s,i,width=4)
        print(s)
        #y axis is time slot
        const3 = np.zeros(Ni)
        for k in range(Nk):
            s = "{:<{width}}".format(k,width=4)
            const1 = 0
            for i in range(Ni):
                idx = h(i,j,k)
                s = "{} | {:<{width}}".format(s,result[idx], width=4)
                if result[idx] == 1:
                    num_obs[i] += 1
                    const1 += 1
                    const3[i] += 1
                #    s = "{} | {}".format(s,"T")
                #else:
                #    s = "{} | {}".format(s,"F")
            assert const1 <= 1, "Too many pointings per slot"
            s = "{} | Sum: {:<{width}}".format(s,const1, width=4)
            print(s)    
        assert np.all(const3 <= 1), "More than one pointing per observation"
        s = "Sum ".format(k,width=4)
        for i in range(Ni):
            s = "{} | {:<{width}}".format(s,const3[i], width=4)
        print(s)    

    print("----")
    s = "Totl".format(k,width=4)
    for i in range(Ni):
        s = "{} | {:<{width}}".format(s,num_obs[i], width=4)
        assert num_obs[i] == Nh, "Not enough obs per pointing"
    print(s)    

    print("")
    print("     | pointing ...")
    print("time |")
    print("...  |")
