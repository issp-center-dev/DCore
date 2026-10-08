"""The GPU opt-in must preserve false values from solver configuration."""
import pytest

create_solver_params = pytest.importorskip('dcore.dmft_core').create_solver_params


@pytest.mark.parametrize('literal, expected', [('False', False), ('false', False),
                                               ('True', True), ('true', True)])
def test_gpu_boolean_literal(literal, expected):
    params = create_solver_params({'gpu{bool}': literal, 'n_bath{int}': '2',
                                   'weight_threshold{float}': '1e-6'})
    assert params['gpu'] is expected
    assert params['n_bath'] == 2
    assert params['weight_threshold'] == 1e-6
