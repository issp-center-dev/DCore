from .wannier90_model import Wannier90Model
from .base import LatticeModel
from .predefined_models import ChainModel, SquareModel, CubicModel, BetheModel
from .external_model import ExternalModel
import sys

def create_lattice_model(params) -> LatticeModel:
    model_name =  params["model"]["lattice"]

    for c in all_lattice_models:
        if model_name == c.name():
            return c(params)

    # raise RuntimeError('Unknown lattice model name: ' + model_name)
    sys.exit(f"ERROR: Unknown lattice model name: [model] lattice={model_name!r}")


all_lattice_models = [ChainModel, SquareModel, CubicModel, BetheModel, Wannier90Model, ExternalModel]
