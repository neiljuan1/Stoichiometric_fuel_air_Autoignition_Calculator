class Species:
    def __init__(self, name:str, mol:float, enth_f: float, cp: float):
        self.name = name
        self.mol = mol
        self.mol_conc = 0 
        self.mol_frac = 0
        self.w = 0
        self.enth_f = enth_f
        self.cp = cp
        self.conc_gradient = 0
        
    def __str__(self)->str:
        return f"{self.name} - {self.mol} moles."
        
        