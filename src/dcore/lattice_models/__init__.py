from .wannier90_model import Wannier90Model
from .predefined_models import ChainModel, SquareModel, CubicModel, BetheModel
from .external_model import ExternalModel

def create_lattice_model(params):
    model_name =  params["model"]["lattice"]

    for c in all_lattice_models:
        if model_name == c.name():
            return c(params)

    raise RuntimeError('Unknown lattice model name: ' + model_name)


all_lattice_models = [ChainModel, SquareModel, CubicModel, BetheModel, Wannier90Model, ExternalModel]
