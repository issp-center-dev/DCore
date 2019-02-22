from .triqs_cthyb import TRIQSCTHYBSolver
from .triqs_hubbard_I import TRIQSHubbardISolver
from .alps_cthyb import ALPSCTHYBSolver
from .null_solver import NullSolver
from .alps_ctseg import ALPSCTSEGSolver

solver_classes = {
    'TRIQS/cthyb': TRIQSCTHYBSolver,
    'TRIQS/hubbard-I': TRIQSHubbardISolver,
    'ALPS/cthyb': ALPSCTHYBSolver,
    'null': NullSolver,
    'ALPS/ctseg': ALPSCTSEGSolver,
}
