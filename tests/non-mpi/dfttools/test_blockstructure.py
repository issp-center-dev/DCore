from dcore.backend.triqs_compat.dft_tools.sumk_dft import *
from dcore.backend.triqs_compat.utility.h5diff import h5diff, compare, failures
from dcore.backend.triqs_compat.gf import *
from dcore.backend.triqs_compat.utility.comparison_tests import assert_block_gfs_are_close
from scipy.linalg import expm
from dcore.backend.triqs_compat.dft_tools.block_structure import BlockStructure
import numpy as np


def cmp(a, b, precision=1.e-15):
    compare('', a, b, 0, precision)
    if failures:
        raise AssertionError('\n'.join(failures))

SK = SumkDFT('blockstructure.in.h5', use_dft_blocks=True)

original_bs = SK.block_structure
cmp(original_bs.effective_transformation_sumk,
    [{'up': np.array([[1., 0., 0.],
                      [0., 1., 0.],
                      [0., 0., 1.]]),
      'down': np.array([[1., 0., 0.],
                        [0., 1., 0.],
                        [0., 0., 1.]])}])
cmp(original_bs.effective_transformation_solver,
    [{'up_0': np.array([[1., 0., 0.],
                        [0., 1., 0.]]),
      'up_1': np.array([[0., 0., 1.]]),
      'down_0': np.array([[1., 0., 0.],
                          [0., 1., 0.]]),
      'down_1': np.array([[0., 0., 1.]])}])

created_matrix = original_bs.create_matrix()
cmp(created_matrix,
    {'up_0': np.array([[0. + 0.j, 0. + 0.j],
                       [0. + 0.j, 0. + 0.j]]),
     'up_1': np.array([[0. + 0.j]]),
     'down_0': np.array([[0. + 0.j, 0. + 0.j],
                         [0. + 0.j, 0. + 0.j]]),
     'down_1': np.array([[0. + 0.j]])})


# check pick_gf_struct_solver
pick1 = original_bs.copy()
pick1.pick_gf_struct_solver([{'up_0': [1], 'up_1': [0], 'down_1': [0]}])

cmp(pick1.effective_transformation_sumk,
    [{'up': np.array([[0., 0., 0.],
                      [0., 1., 0.],
                      [0., 0., 1.]]),
     'down': np.array([[0., 0., 0.],
                        [0., 0., 0.],
                        [0., 0., 1.]])}])

cmp(pick1.effective_transformation_solver,
    [{'up_0': np.array([[0., 1., 0.]]),
      'up_1': np.array([[0., 0., 1.]]),
      'down_1': np.array([[0., 0., 1.]])}])

# check loading a block_structure from file
SK.block_structure = SK.load(['block_structure'], 'mod')[0]
assert SK.block_structure == pick1, 'loading SK block structure from file failed'

# check SumkDFT backward compatibility
sk_pick1 = BlockStructure(gf_struct_sumk=SK.gf_struct_sumk,
                          gf_struct_solver=SK.gf_struct_solver,
                          solver_to_sumk=SK.solver_to_sumk,
                          sumk_to_solver=SK.sumk_to_solver,
                          solver_to_sumk_block=SK.solver_to_sumk_block,
                          deg_shells=SK.deg_shells,
                          corr_to_inequiv=SK.corr_to_inequiv)
assert sk_pick1 == pick1, 'constructing block structure from SumkDFT properties failed'

cmp(pick1.effective_transformation_sumk,
    [{'up': np.array([[0., 0., 0.],
                      [0., 1., 0.],
                      [0., 0., 1.]]),
      'down': np.array([[0., 0., 0.],
                        [0., 0., 0.],
                        [0., 0., 1.]])}])

cmp(pick1.effective_transformation_solver,
    [{'up_0': np.array([[0., 1., 0.]]),
      'up_1': np.array([[0., 0., 1.]]),
      'down_1': np.array([[0., 0., 1.]])}])

# check pick_gf_struct_sumk
pick2 = original_bs.copy()
pick2.pick_gf_struct_sumk([{'up': [1, 2], 'down': [0, 1]}])

cmp(pick2.effective_transformation_sumk,
    [{'up': np.array([[0., 0., 0.],
                      [0., 1., 0.],
                      [0., 0., 1.]]),
    'down': np.array([[1., 0., 0.],
                        [0., 1., 0.],
                        [0., 0., 0.]])}])

cmp(pick2.effective_transformation_solver,
    [{'up_0': np.array([[0., 1., 0.]]),
      'up_1': np.array([[0., 0., 1.]]),
      'down_0': np.array([[1., 0., 0.],
                          [0., 1., 0.]])}])

pick3 = pick2.copy()
pick3.transformation = [np.reshape(range(9), (3, 3))]
cmp(pick3.effective_transformation_sumk,
    [{'up': np.array([[0, 0, 0],
                      [3, 4, 5],
                      [6, 7, 8]]),
      'down': np.array([[0, 1, 2],
                        [3, 4, 5],
                        [0, 0, 0]])}])

cmp(pick3.effective_transformation_solver,
    [{'up_0': np.array([[3, 4, 5]]),
      'up_1': np.array([[6, 7, 8]]),
      'down_0': np.array([[0, 1, 2],
                          [3, 4, 5]])}])

pick4 = original_bs.copy()
pick4.transformation = [np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]])]
pick4.pick_gf_struct_sumk([{'up': [1, 2], 'down': [0, 1]}])
cmp(pick2.gf_struct_sumk, pick4.gf_struct_sumk)
cmp(pick2.gf_struct_solver, pick4.gf_struct_solver)
assert pick4.sumk_to_solver == [{('up', 0): ('up_0', 0),
                                 ('up', 1): (None, None),
                                 ('up', 2): ('up_1', 0),
                                 ('down', 2): (None, None),
                                 ('down', 1): ('down_0', 1),
                                 ('down', 0): ('down_0', 0)}]
assert pick4.solver_to_sumk == [{('up_1', 0): ('up', 2),
                                 ('up_0', 0): ('up', 0),
                                 ('down_0', 0): ('down', 0),
                                 ('down_0', 1): ('down', 1)}]

# check map_gf_struct_solver
mapping = [{('down_0', 0): ('down', 0),
            ('down_0', 1): ('down', 2),
            ('down_1', 0): ('down', 1),
            ('up_0', 0): ('down_1', 0),
            ('up_0', 1): ('up_0', 0)}]
map1 = original_bs.copy()
map1.map_gf_struct_solver(mapping)

# check create_gf
G1 = original_bs.create_gf(beta=40, n_points=3)
widths = dict(up_0=1, up_1=2, down_0=4, down_1=3)
for block, gf in G1:
    gf << SemiCircular(widths[block])
original_bs.check_gf(G1)
original_bs.check_gf([G1])

# check approximate_as_diagonal
offd = original_bs.copy()
offd.approximate_as_diagonal()

# check map_gf_struct_solver
import warnings
with warnings.catch_warnings(record=True) as w:
    G2 = map1.convert_gf(G1, original_bs, beta=40, n_points=3,
                         show_warnings=True)
    assert len(w) == 1
    assert issubclass(w[-1].category, UserWarning)
    assert "Block up_1 maximum difference" in str(w[-1].message)

m2 = map1.convert_matrix(created_matrix, original_bs, show_warnings=True)
cmp(m2,
    {'down': np.array([[0. + 0.j, 0. + 0.j, 0. + 0.j],
                       [0. + 0.j, 0. + 0.j, 0. + 0.j],
                       [0. + 0.j, 0. + 0.j, 0. + 0.j]]),
     'down_1': np.array([[0. + 0.j]]),
     'up_0': np.array([[0. + 0.j]])})

# check full_structure
full = BlockStructure.full_structure(
    [{'up_0': [0, 1], 'up_1': [0], 'down_1': [0], 'down_0': [0, 1]}], None)

G_sumk = BlockGf(mesh=G1.mesh, gf_struct=original_bs.gf_struct_sumk[0])
for i in range(3):
    G_sumk['up'][i, i] << SemiCircular(1 if i < 2 else 2)
    G_sumk['down'][i, i] << SemiCircular(4 if i < 2 else 3)
G3 = original_bs.convert_gf(G_sumk,
                            None,
                            space_from='sumk',
                            beta=40,
                            n_points=3)
assert_block_gfs_are_close(G1, G3)

# check convert_gf with transformation
# np.random.seed(894892309)
H = np.random.rand(3, 3) + 1.0j * np.random.rand(3, 3)
H = H + H.conjugate().transpose()
T = expm(1.0j * H)
G_T = G_sumk.copy()
for block, gf in G_T:
    gf.from_L_G_R(T.conjugate().transpose(), gf, T)
transformed_bs = original_bs.copy()
transformed_bs.transformation = [T]
G_bT = transformed_bs.convert_gf(G_T, None, space_from='sumk',
                                 beta=40, n_points=3)
assert_block_gfs_are_close(G1, G_bT)

assert original_bs.gf_struct_sumk_list ==\
    [[('up', [0, 1, 2]), ('down', [0, 1, 2])]]
assert original_bs.gf_struct_solver_dict ==\
    [{'up_0': [0, 1], 'up_1': [0], 'down_1': [0], 'down_0': [0, 1]}]
assert original_bs.gf_struct_sumk_dict ==\
    [{'down': [0, 1, 2], 'up': [0, 1, 2]}]
assert original_bs.gf_struct_solver_list ==\
    [[('down_0', [0, 1]), ('down_1', [0]), ('up_0', [0, 1]), ('up_1', [0])]]

# check __eq__
assert full == full, 'equality not correct (equal structures not equal)'
assert pick1 == pick1, 'equality not correct (equal structures not equal)'
assert pick1 != pick2, 'equality not correct (different structures not different)'
assert original_bs != offd, 'equality not correct (different structures not different)'

if mpi.is_master_node():
    with HDFArchive('blockstructure.out.h5', 'w') as ar:
        ar['original_bs'] = original_bs
        ar['pick1'] = pick1
        ar['pick2'] = pick2
        ar['map1'] = map1
        ar['offd'] = offd
        ar['G1'] = G1
        ar['G2'] = G2
        ar['full'] = full

    # cannot use h5diff because BlockStructure testing is not implemented
    # there (and seems difficult to implement because it would mix triqs
    # and dft_tools)
    with HDFArchive('blockstructure.out.h5', 'r') as ar,\
            HDFArchive('blockstructure.ref.h5', 'r') as ar2:
        for k in ar2:
            print(k)
            if isinstance(ar[k], BlockGf):
                assert_block_gfs_are_close(ar[k], ar2[k], 1.e-6)
            else:
                assert ar[k] == ar2[k], '{} not equal'.format(k)
