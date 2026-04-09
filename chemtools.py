import csv
import re
from collections import Counter
import numpy as np
import sympy

weight_dict = dict()
"""
Converts Ptable.csv to a readable dictionary containing chemical symbols as keys, 
associated atomic mass as the value. 
"""

try:
    with open("/home/timokneppers/Projects/Labtools/Ptable.csv") as pt:
        reader = csv.DictReader(pt)
        for row in reader:
            symbol = row["Symbol"]
            mass = float(row["AtomicMass"])
            weight_dict[symbol] = mass
except FileNotFoundError:
    print("File Ptable.csv not found in correct dir")
    exit()

def molarweight(formula):
    """
    Calculates the molar mass of a given compound

    Args:
        formula: Compound as a string input
    Return:
        Molar mass in g/mol
    Raises:
        TypeError: If the input is not a string type
        ValueError: If the compound is not a real chemical compound
    """

    if not isinstance(formula, str):
        raise TypeError("Input type should be a string")
    
    formula = formula.strip()
    
    tokens = re.findall(r"[A-Z][a-z]?|[\d]+|[()\[\]{}]", formula)

    if not "".join(tokens) == formula:
        raise ValueError("Input contains incorrect chemical data")

    stack = [0]
    buffer = 0

    for token in tokens:
        match token:
            case e if e in weight_dict.keys():
                stack[-1] += buffer
                buffer = weight_dict[e]
            case t if t.isnumeric():
                buffer *= int(t)
            case "(" | "[" | "{":
                stack[-1] += buffer
                buffer = 0
                stack.append(0)
            case ")" | "]" | "}":
                stack[-1] += buffer
                buffer = stack.pop()
            case _:
                raise ValueError(f"Onbekend symbool: {token}")
    stack[-1] += buffer
    mass = stack[0]

    return mass

def chemdict(formula):
    """
    Returns a dictionary of all elements and their amount for a given chemical compound

    Args:
        formula: Chemical formula as a string input
    Returns:
        Dictionary with the element as keys, the amount as values
    Raises:
        TypeError: For incorrect input type
        ValueError: For incorrect chemical data,
    """

    if not isinstance(formula, str):
        raise TypeError("Compound must be a string input")
    formula = formula.strip()

    tokens = re.findall(r"[A-Z][a-z]?|[\d]+|[()\[\]{}]", formula)
    
    if not "".join(tokens) == formula:
        raise ValueError("Input contains incorrect chemical data")
    
    stack = list()
    stack.append(Counter())
    buffer = Counter()
    element_set = set(re.findall(r"[A-Z][a-z]?", formula))

    if not element_set.issubset(weight_dict.keys()):
        raise ValueError("The input contains a false element")    

    for n in tokens:
        match n:
            case e if e in element_set:
                if buffer:
                    stack[-1] = stack[-1] + buffer
                buffer = Counter({e: 1})
            case "(" | "[" | "{":
                if buffer:
                    stack[-1] = stack[-1] + buffer
                stack.append(Counter())
                buffer = Counter()
            case ")" | "]" | "}":
                if buffer:
                    stack[-1] = stack[-1] + buffer
                    buffer = stack.pop()
            case t if t.isnumeric():
                buffer = Counter({m: k * int(n) for m, k in buffer.items()})
            case _:
                raise ValueError("Unknow symbol being parsed")
    if len(stack) != 1:
        raise ValueError("Brackets not properly closed")
    product = stack[-1] + buffer
    product = dict(product)
    

    return product

def react_solver(reactants: list, products: list):
    if not isinstance(reactants, (list, tuple)):
        raise TypeError("Reactants input should be either a list or tuple")
    if not isinstance(products, (list, tuple)):
        raise TypeError("Products input should be either a list or a tuple")

    react_dict_list = [chemdict(r) for r in reactants]
    prod_dict_list = [chemdict(p) for p in products]

    elements = {key for d in (react_dict_list + prod_dict_list) for key in d}
    elements_dict = {el:i for i, el in enumerate(elements)}

    react_vectors = []
    prod_vectors = []

    for react in react_dict_list:
        vector = np.zeros(len(elements_dict))
        for k, v in react.items():
            vector[elements_dict[k]] = v
        react_vectors.append(vector)
    for prod in prod_dict_list:
        vector = np.zeros(len(elements_dict))
        for k, v in prod.items():
            v *= -1
            vector[elements_dict[k]] = v
        prod_vectors.append(vector)
    
    matrix = np.array(react_vectors + prod_vectors, dtype="i")
    matrix = matrix.T
    ns_list = sympy.Matrix(matrix).nullspace()[0]
    ns_list = ns_list * (1 / ns_list[0])
    return ns_list

print(chemdict("Cr[Fe(OH)2(NO3)2Cl2]3"))
print(react_solver(["Fe(OH)2"], ["Fe3O4", "H2", "H2O"]))