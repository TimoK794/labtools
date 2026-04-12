from chem_functions import molarweight

class Compound:
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
    
class Mixture:
    def __init__(self, volume:float, solute:Compound=Compound("H2O", "977 kg/m3", "273.15 K", "373.15 K")) -> None:
        if not isinstance(volume, (int, float)):
            raise TypeError("")
        self.volume = volume
        self.solute = solute
        if not solute.density:
            raise ValueError("Solute needs density attribute")
        self._solute_density = solute.density
        self.contents = dict()
    
    def add(self, compound:Compound, concentration:str) -> None:
        if not isinstance(compound, Compound):
            raise ValueError("Can only add compounds to a mixture")
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
            self.contents[compound] = value * conversion[unit]
        except AttributeError:
            raise ValueError("Input should be str of a value followed by a unit")

sulfuric_acid = Compound("H2SO4")
mix1 = Mixture(1)
mix1.add(sulfuric_acid, "1 M")
print(mix1.contents)