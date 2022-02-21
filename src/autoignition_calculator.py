import matplotlib.pyplot as plt
from conditions import Conditions
from species import Species
import math
# Define constants
A = 4.6e+11
ACT_ENERGY   = 15098
M    = 0.25
N    = 1.5

REFERENCE_TEMP = 298
DT = 1e-7
GAS_CONST = 8.314

class AutoignitionCalculator:
    def __init__(self):
        self.conditions = Conditions()
        self.reactant_species = {}
        self.product_species = {}
        self.t = 0
        self.temp = self.conditions.temp
        self.w_sum = 0
        self.mol_conc_sum = 0
        self.temp_gradient = 0

        self.initialize_species()
        self.calculate_mol_fraction()
        self.calculate_mol_conc()
        self.calculate_w()
        
    def initialize_species(self):
        self.reactant_species["fuel"] = Species("octane", 1, -208700.0, 431.37) # format - name, mol, enthalpy of formation, specific heat cp
        self.reactant_species["o2"] = Species("oxygen", 12.5, 0, 34.936)
        self.reactant_species["n2r"] = Species("nitrogen", 3.76*12.5, 0, 32.762)
        
        # products
        self.product_species["co2"] = Species("carbon dioxide", 8, -393546, 54.36)
        self.product_species["h2o"] = Species("water", 9, -241845, 41.315)
        self.product_species["n2p"] = Species("nitrogen", 3.76*12.5, 0, 32.762)
        
        
    def calculate_mol_fraction(self):
        reactant_sum = sum([spec.mol for spec in self.reactant_species.values()])
        product_sum = sum([spec.mol for spec in self.product_species.values()])
        
        for name, spec in self.reactant_species.items():
            spec.mol_frac = spec.mol / reactant_sum
            #print(name, spec, spec.mol_frac)
            
        for name, spec in self.product_species.items():
            spec.mol_frac = spec.mol / product_sum
            #print(name, spec, spec.mol_frac)
            

    def calculate_mol_conc(self):
        
        for name, spec in self.reactant_species.items():
            spec.mol_conc = spec.mol_frac * (self.conditions.pressure_bar * self.conditions.pressure) / (GAS_CONST * self.temp) / 1e6
            #print(name, spec, spec.mol_conc)
            self.mol_conc_sum += spec.mol_conc
            
        for name, spec in self.product_species.items():
            spec.mol_conc = spec.mol_frac * (self.conditions.pressure_bar * self.conditions.pressure) / (GAS_CONST * self.temp) / 1e6
            #print(name, spec, spec.mol_conc)
            
        self.calculate_mol_conc_sum()
    
    def calculate_mol_conc_sum(self):
        self.mol_conc_sum = 0
        for name, spec in self.get_species().items():
            self.mol_conc_sum += spec.mol_conc
            
          
            
    def calculate_w(self):
        self.reactant_species["fuel"].w = -A*math.exp(-ACT_ENERGY/self.temp)*(self.reactant_species["fuel"].mol_conc**M)*(self.reactant_species["o2"].mol_conc**N)
        self.reactant_species["o2"].w = self.reactant_species["fuel"].w * self.reactant_species["o2"].mol
        self.product_species["co2"].w = -self.reactant_species["fuel"].w * self.product_species["co2"].mol
        self.product_species["h2o"].w = -self.reactant_species["fuel"].w * self.product_species["h2o"].mol
        
        self.w_sum = -3.5 * self.reactant_species["fuel"].w
        
        #[print(a, a.w) for a in self.reactant_species.values()]
        #[print(a, a.w) for a in self.product_species.values()]
        #print(self.w_sum)
        
        
        
    def calculate_enthalpy(self, sp: Species)->float:
        """
        returns enthalpy*rxn rate at a time t
        """
        return (sp.enth_f + sp.cp * (self.temp - REFERENCE_TEMP)) * sp.w
    
    def bot(self, sp: Species)->float:
        """returns conc * cp at time i"""
        return sp.mol_conc * sp.cp
    
    def calculate_temp_gradient(self):
        
        gradient_numerator = sum([self.calculate_enthalpy(spec) for spec in self.get_species().values()])
        gradient_denom = sum([self.bot(spec) for spec in self.get_species().values()])
        
        self.temp_gradient = -(gradient_numerator / gradient_denom)
    
    def conc_gradient(self, sp: Species):
        return sp.w - sp.mol_conc * ((self.w_sum / self.mol_conc_sum) + (self.temp_gradient / self.temp))
        
    def get_species(self):
        return {**self.reactant_species, **self.product_species}
        

    

if __name__ == "__main__":
    ai = AutoignitionCalculator()
    conditions = Conditions()
    #solution_vector = {
    #    "t":    [0],
    #    "T":    [conditions.temp],
    #    "fuel": [],
    #    "o2":   [],
    #    "n2":   [],
    #    "co2":  [],
    #    "h2o":  [],
    #    "n2p":  []        
    #}
    
    sv = {
        "time":    [0],
        "T":    [conditions.temp],
        "species":{
            "fuel": [ai.get_species()["fuel"].mol_conc],
            "o2":   [ai.get_species()["o2"].mol_conc],
            "n2r":  [ai.get_species()["n2r"].mol_conc],
            "co2":  [ai.get_species()["co2"].mol_conc],
            "h2o":  [ai.get_species()["h2o"].mol_conc],
            "n2p":  [ai.get_species()["n2p"].mol_conc]}
    }
    
    i = 1
    tau = 5.8e-4
    t = 0
   
    ai.calculate_temp_gradient()
    print(t)
    
    while (t < tau):
        t += DT
        #print(t)
        sv["time"].append(t)
        ai.calculate_temp_gradient()

        
        #print("temp grad", ai.get_species()["fuel"].mol_conc, "wsum", ai.w_sum, "xsum", ai.mol_conc_sum)
        
        
        for spec in sv["species"]:
            #print("conc", ai.conc_gradient(ai.get_species()[spec]))
            ai.get_species()[spec].mol_conc = ai.get_species()[spec].mol_conc + ai.conc_gradient(ai.get_species()[spec])*DT
            #print(ai.get_species()[spec].mol_conc)
            sv["species"][spec].append(ai.get_species()[spec].mol_conc)
            #sv["species"][spec].append(sv["species"][spec][i-1] + ai.conc_gradient(ai.get_species()[spec]))
            
        # update Temperature
        ai.temp = ai.temp + ai.temp_gradient * DT
        sv["T"].append(ai.temp)
        
        
        # Update w
        ai.calculate_w()
        
        # update concentrations
        ai.calculate_mol_conc_sum()
    
        i += 1
    print(sv["species"]["fuel"])
    plt.plot(sv["time"], sv["species"]["fuel"])
    plt.show()
    plt.plot(sv["time"], sv["T"])
    plt.show()

    
    print("end")
        
    

    
