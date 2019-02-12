from .triqs_cthyb import TRIQSCTHYBSolver
from .triqs_hubbard_I import TRIQSHubbardISolver
from .alps_cthyb import ALPSCTHYBSolver
from .null_solver import NullSolver
from .pomerol import PomerolSolver

solver_classes = {
    'TRIQS/cthyb' : TRIQSCTHYBSolver,
    'TRIQS/hubbard-I' : TRIQSHubbardISolver,
    'ALPS/cthyb': ALPSCTHYBSolver,
    'null': NullSolver,
    'pomerol': PomerolSolver,
}
