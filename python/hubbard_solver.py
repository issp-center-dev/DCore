import pytriqs.applications.impurity_solvers.hubbard_I

class Solver(pytriqs.applications.impurity_solvers.hubbard_I.Solver, object):
    """
       Hubbard I Solver
    """

    # initialisation:
    def __init__(self, beta, l, n_msb=1025, use_spin_orbit=False, Nmoments=5):
        super(Solver, self).__init__(beta, l, n_msb, use_spin_orbit, Nmoments)

    # Make read-only getter
    @property
    def Sigma_iw(self):
        return self.Sigma

    # Make read-only getter
    @property
    def G_iw(self):
        return self.G

    # Make read-only getter
    @property
    def G0_iw(self):
        return self.G0
