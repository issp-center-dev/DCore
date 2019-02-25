from .triqs_cthyb import TRIQSCTHYBSolver
from .triqs_hubbard_I import TRIQSHubbardISolver
from .alps_cthyb import ALPSCTHYBSolver
from .null_solver import NullSolver
from .alps_cthyb_seg import ALPSCTHYBSEGSolver
from .pomerol import PomerolSolver

solver_classes = {
    'TRIQS/cthyb': TRIQSCTHYBSolver,
    'TRIQS/hubbard-I': TRIQSHubbardISolver,
    'ALPS/cthyb': ALPSCTHYBSolver,
    'null': NullSolver,
    'ALPS/cthyb-seg': ALPSCTHYBSEGSolver,
    'pomerol': PomerolSolver,
}
