from chem_functions import validate_molecule, molarweight, react_solver

class Compound:
    """
    Class to form a compound

    Args: 
        formula: string input of the chemical formula of the compound
        density: Density string input with unit, will be converted to kg/m3
        meltingpoint: string input with unit, will be converted to K
        boilingpoint: string input with unit, will be converted to K
    """
    def __init__(self, formula:str, density:str=None, meltingpoint:str=None, boilingpoint:str=None) -> None:
        self.formula = formula
        self.weight = molarweight(formula)
        self.density = density
        self.meltingpoint = meltingpoint
        self.boilingpoint = boilingpoint
    
    @property
    def density(self):
        return self._density
    
    @density.setter
    def density(self, data):
        if data == None:
            self._density = None
            return
        try:
            value, unit = data.split()
            value = float(value)
            conversion = {
                "kg/m3": 1,
                "g/cm3": 1000,
                "kg/L": 1000
            }
            if not unit in conversion:
                raise ValueError(f"Unit: {unit} not found in database")
            
            self._density = value * conversion[unit]
        except (ValueError, AttributeError):
            raise ValueError(f"Input should be a value(float) follow by a unit(str)")

    @property
    def meltingpoint(self):
        return self._meltingpoint
    
    @meltingpoint.setter
    def meltingpoint(self, data):
        if data == None:
            self._meltingpoint = None
            return
        try:
            value, unit = data.split()
            value = float(value)
            conversion = {
                "K": 0,
                "C": 273.15
            }
            if not unit in conversion:
                raise ValueError(f"Unit: {unit} not found in database")
            
            self._meltingpoint = value + conversion[unit]
        except (ValueError, AttributeError):
            raise ValueError(f"Input should be a value(float) follow by a unit(str)")
        
    @property
    def boilingpoint(self):
        return self._boilingpoint
    
    @boilingpoint.setter
    def boilingpoint(self, data):
        if data == None:
            self._boilingpoint = None
            return
        try:
            value, unit = data.split()
            value = float(value)
            conversion = {
                "K": 0,
                "C": 273.15
            }
            if not unit in conversion:
                raise ValueError(f"Unit: {unit} not found in database")
            
            self._boilingpoint = value + conversion[unit]
        except (ValueError, AttributeError):
            raise ValueError(f"Input should be a value(float) follow by a unit(str)")        

    def __str__(self) -> str:
        result = f'Compound: {self.formula} ({self.weight} g/mol)'
        if self.density:
            result += f"\n\tdensity: {self.density} kg/m^3"
        if self._meltingpoint:
            result += f"\n\tmelting point: {self.meltingpoint} K"
        if self._boilingpoint: 
            result += f"\n\tboiling point: {self.boilingpoint} K"
        return result

class Reaction:
    def __init__(self, reactants:list, products:list):
        if not isinstance(reactants, (tuple, list)):
            raise ValueError("reactans should be either a list or a tuple")
        map(validate_molecule, (e.formula for e in reactants))
        self.reactants = list(reactants)
        if not isinstance(reactants, (tuple, list)):
            raise ValueError("products should be either a list or a tuple")
        map(validate_molecule, (e.formula for e in products))
        self.products = list(products)
        try:
            self.solution = react_solver(list(e.formula for e in reactants), list(g.formula for g in products))
        except ValueError:
            raise ValueError("Reaction is not solvable")

    
class Mixture:
    def __init__(self, volume:float, solute:Compound=Compound("H2O", "977 kg/m3", "273.15 K", "373.15 K")) -> None:
        if not isinstance(volume, (int, float)):
            raise TypeError("")
        self.volume = volume
        self.solute = solute
        if not solute.density:
            raise ValueError("Solute needs density attribute")
        self._solute_density = solute.density
        self._contents = dict()

    @property
    def contents(self) -> list:
        return [(k.formula, v) for k, v in self._contents.items()]

    def add(self, compound:Compound, concentration:str) -> None:
        if not isinstance(compound, Compound):
            raise ValueError("Can only add compounds to a mixture")
        if compound == self.solute:
            raise ValueError("Cannot add a compound which is the same as the solute")
        try:
            value, unit = concentration.split()
            value = float(value)
            conversion= {
                "M": 1,
                "mM": 0.001,
                "kM": 1000
            }
            if unit not in conversion:
                raise ValueError("Unit not in conversion database")
            self._contents[compound] = value * conversion[unit]
        except AttributeError:
            raise ValueError("Input should be str of a value followed by a unit")