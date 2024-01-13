from pymatgen.core.structure import Structure

class Calculator:

  def __init__ ( self ):
    pass

  def predict_formation_energy ( self, lattice, species, positions ):
    pass



class MEGNet_Calculator (Calculator):

  model = None

  def __init__ ( self, model='Eform_MP_2019', model_fname=None ):

    if model_fname is None:
      from megnet.utils.models import load_model
      self.model = load_model(model)

    else:
      from megnet.models import MEGNetModel
      self.model = MEGNetModel.from_file(model_fname)


  def predict_formation_energy ( self, lattice, species, positions ):
    pymatgen_struct = Structure(lattice, species, positions)
    return self.model.predict_structure(pymatgen_struct).ravel()[0]
