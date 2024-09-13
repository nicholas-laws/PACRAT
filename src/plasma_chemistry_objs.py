# Created by: Nicholas Laws
# Date: 2024

class Reaction:
    def __init__(self, equation, rate, species):
        self.equation = equation
        self.rate = rate
        self.species = species

class Equation:
    def __init__(self, template, reactions):
        self.template = template
        self.reactions = reactions

class Rate:
    def __init__(self, template, rates):
        self.template = template
        self.rates = rates

class Species:
    def __init__(self, reactants, products):
        self.reactants = reactants
        self.products = products