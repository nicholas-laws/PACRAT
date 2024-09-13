# Created by: Nicholas Laws
# Date: 2024

## Imports ##
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx
import cantera as ct
from scipy.constants import Avogadro
import matplotlib.colors as mcolors
from networkx.drawing.nx_agraph import graphviz_layout
from netgraph import Graph
import re
from scipy import integrate
from typing import *
import time

"""
A function to collect the name of each combustion species tracked in the numerical model.
"""
def get_species(mechanism: str) -> List[str]:
    """
    Parameters:
    -----------
    mechanism : str
        A string representing the filepath to the combustion mechanism of interest.

    Returns:
    --------
    species : List[str]
        A list of strings representing the chemical names of each species tracked in the combustion mechanism.
    """
    # Leverage Cantera to list all species objects from the combustion mechanism of interest.
    S = ct.Species.list_from_file(mechanism, section='species')
    # Initialize an empty list to collect the species strings.
    species = []
    # Iterate through each species objects from the combustion mechanism of interest.
    for i in S:
        # For each iteration, collect the current species' chemical name and append to the species list.
        temp_s = i.name
        species.append(temp_s)
    # Return the completed species string list.
    return species

"""
A function to collect the reaction objects tracked in the numerical model and combustion mechanism of interest.
"""
def get_reactions(mechanism: str) -> List[ct.Reaction]:
    """
    Parameters:
    -----------
    mechanism : str
        A string representing the filepath to the combustion mechanism of interest.
    Returns:
    --------
    rxns : List[ct.Reaction]
        A list of Cantera reaction objects representing the reactions tracked in the combustion mechanism.
    """
    # Generate a placeholder gas object using Cantera and the combustion mechanism of interest.
    gas = ct.Solution(mechanism)
    # List the reaction objects from the combustion mechansim of interest using Cantera.
    rxns = ct.Reaction.list_from_file(mechanism, gas)
    # Return the reaction objects.
    return rxns

"""
A function to collect the reaction equation strings tracked in the numerical model and combustion mechanism of interest.
"""
def get_reaction_strings(mechanism: str) -> List[str]:
    """
    Parameters:
    -----------
    mechanism : str
        A string representing the filepath to the combustion mechanism of interest.

    Returns:
    --------
    reactions : List[str]
        A list of strings representing the reaction equations tracked in the combustion mechanism.
    """
    # Generate a placeholder gas object using Cantera and the combustion mechanism of interest.
    gas = ct.Solution(mechanism)
    # List the reaction objects from the combustion mechansim of interest using Cantera.
    rxns = ct.Reaction.list_from_file(mechanism, gas)
    # Initialize an empty list to collect the reaction equations list. 
    reactions = []
    # Iterate through each reaction objects from the combustion mechanism of interest.
    for j in rxns:
        # For each iteration, collect the current reaction's equation and append to the reaction equations list.
        temp_r = j.equation
        reactions.append(temp_r)
    # Return the completed reaction equations string list.
    return reactions

"""
A function to identify all reaction pathways present in a combustion reaction of interest.
"""
def pathways_identifier(reaction: ct.Reaction, reactants: List[str], products: List[str], element_of_interest: str) -> List[str]:
    """
    Parameters:
    -----------
    reaction : ct.Reaction
        An object representing the combustion reaction of interest.
    reactants : List[str]
        A list of strings representing the reactants involved in the combustion reaction of interest.
    products : List[str]
        A list of strings representing the products involved in the combustion reaction of interest.
    element_of_interest : str
        A string representing the element to be tracked in the reaction pathway analysis.

    Returns:
    --------
    pathways : List[str]
        A list of strings representing all identified reaction pathways for the element of interest.
    """
    # Initialize an empty reaction pathways list.
    pathways = []
    # Initialize an empty variable representing each possible reaction pathway, to be iteratively changed.
    temp_pathway = []
    # Iterate through each reactant, then each character in the reactant as well each product.
    for reactant in reactants:
            for char in reactant:
                # Check if the reactant character is a digit, if not skip this character.
                if char.isdigit() == False:
                    for product in products:
                        if char in product:
                            # Check if the reactant is not the product and if the element of interest is in the reactant and the product.
                            if reactant != product and element_of_interest in reactant and element_of_interest in product:
                                # Pathway identification rule.
                                if len(product) > 2 and sum(1 for element in product if element_of_interest in element) == 1:
                                    temp_pathway = f"{reactant} → {product}"
                                else:
                                    pass
                            # Verify that the pathway is not pathways list and that the pathway has been initialized to reactant -> product.
                            if temp_pathway not in pathways and temp_pathway == f"{reactant} → {product}":
                                # Append the temporary pathway to the complete pathways list.
                                pathways.append(temp_pathway)
                                # Check if the reaction is reversible, if so, reinitialize the temporary pathway to the reverse pathway.
                                if reaction.reversible == True:
                                    temp_pathway = f"{product} → {reactant}"
                                # Ensure that the temporary pathway is not in the complete pathways list and append the temporary pathway to the complete pathway list.
                                if temp_pathway not in pathways:
                                    pathways.append(temp_pathway)
    # Return the complete pathways list.
    return pathways

"""
A function to identify all reaction pathways present in a combustion reaction of interest.
"""
def find_pathways(reaction: ct.Reaction, element_of_interest: str) -> Dict[str, List[Union[ct.Reaction, List[str]]]]:
    """
    Parameters:
    -----------
    reaction : ct.Reaction
        An object representing the combustion reaction of interest.
    element_of_interest : str
        A string representing the element to be tracked in the reaction pathway analysis.

    Returns:
    --------
    pathways : Dict[str, List[Union[ct.Reaction, List[str]]]]
        A dictionary where the key is the reaction equation (as a string), and the value is a list containing:
        - The reaction object.
        - A list of all relevant pathways involving the element of interest.
    """
    # Access the data for the combustion reaction of interest to initialize variables for the reaction type and equation.
    rxn_type = reaction.reaction_type
    rxn_eq = reaction.equation
    # Check if the combustion reaction of interest is Arrhenius.
    if rxn_type == "Arrhenius":
        # Gather lists of all the reactants and products.
        reactants = list(reaction.reactants.keys())
        products = list(reaction.products.keys())
        # Construct a list of all pathways in the combustion reaction of interest.
        temp_pathways = pathways_identifier(reaction, reactants, products, element_of_interest)
        # Initialize an empty dictionary where the key is the chemical reaction string and the value is a list where 
        # the first element is the reaction object and the second element is a list of all pathways relevant to the 
        # combustion reaction of interest.
        pathways = {}
        # Verify that the reaction is not in the pathways dictionary, if not, add the key-value pair.
        if rxn_eq not in pathways:
                pathways[rxn_eq] = [reaction, temp_pathways]
        # Return the pathways dictionary.
        return pathways

"""
A function to gather the rate constants of a combustion reaction of interest.
"""
def get_rate_constants(reaction: ct.Reaction) -> List[float]:
    """
    Parameters:
    -----------
    reaction : ct.Reaction
        An object representing the combustion reaction of interest.

    Returns:
    --------
    rate_constants : List[float]
        A list containing the rate constants for the reaction:
        - A : Pre-exponential factor.
        - b : Temperature exponent.
        - Ea : Activation energy (in appropriate units).
    """
    # Initialize a variable representing the reaction type by accessing the data for the combustion reaction of interest.
    rxn_type = reaction.reaction_type
    # Check if the reaction type is Arrhenius.
    if rxn_type == "Arrhenius":
        # Initialize three variables that represent the pre-exponential factor, the temperature exponent, and the activation energy.
        A = float(reaction.rate.input_data["rate-constant"]["A"])
        b = float(reaction.rate.input_data["rate-constant"]["b"])
        Ea = float(reaction.rate.input_data["rate-constant"]["Ea"])
        # Construct list for the rate constants.
        rate_constants = [A, b, Ea]
    # Return the rate constants.
    return rate_constants

"""
A function to access the time-dependent reverse reaction rate data for a combustion reaction of interest.
"""
def reverse_rate(data: pd.DataFrame, reaction_eqn: str) -> pd.Series:
    """
    Parameters:
    -----------
    data : pd.DataFrame
        A Pandas DataFrame storing all time-dependent reverse reaction rates for all combustion reactions in the numerical model.
    reaction_eqn : str
        A string representing the combustion reaction equation of interest.

    Returns:
    --------
    pd.Series
        A Pandas Series containing the time-dependent reverse reaction rate data for the specified combustion reaction.
    """
    # Return the time-dependent reverse reaction rate data for the combustion reaction of interest.
    return data[reaction_eqn]

"""
A function to solve for the reaction rate for a pathway-reaction pair.
"""
def solve_reaction_rate(reverse_rate_coeff_data: pd.DataFrame, pathway: str, reaction_obj: ct.Reaction, 
                        rate_constants: List[float], T: np.array, R: float, P0: float, 
                        chi_data: pd.DataFrame, efficiency: float = 1) -> List[Union[np.array, int, int]]:
    """
    Parameters:
    -----------
    reverse_rate_coeff_data : pd.DataFrame
        A Pandas DataFrame storing all time-dependent reverse reaction rate data for all combustion reactions in the numerical model.
    pathway : str
        A string representing the combustion reaction pathway of interest.
    reaction_obj : ct.Reaction
        An object representing the combustion reaction of interest.
    rate_constants : List[float]
        A list representing the rate constants of the combustion reaction of interest.
    T : np.array
        A numpy array representing the time-dependent temperature of the gas mixture model.
    R : float
        A float representing the universal gas constant.
    P0 : float
        A float representing the initial pressure of the gas mixture model.
    chi_data : pd.DataFrame
        A Pandas DataFrame storing all time-dependent mass fraction data for all combustion reactions in the numerical model.
    efficiency : float, optional
        A float representing the efficiency of three-body Arrhenius reactions in the numerical model. Default is 1.

    Returns:
    --------
    [q_t, tstart, tend] : List[Union[np.array, int, int]]
        A list containing:
        - q_t : np.array
            The time-dependent reaction rate for the pathway-reaction pair.
        - tstart : int
            The start index of the time array for the interval considered.
        - tend : int
            The end index of the time array for the interval considered.
    """
    # Initialize a variable that represents the combustion reaction of interest's type.
    reaction_type = reaction_obj.reaction_type
    # Split the reaction pathway of interest into an incident and product species.
    pathways = pathway.split(" → ")
    j = pathways[0]
    k = pathways[1]
    # Initialize the time array for the mass fractions.
    time = chi_data["t"].to_numpy()
    # Determine the time boundaries for the pulse and interpulse periods.
    [tt_pulse_bounds, tt_interpulse_bounds] = find_pulse_interpulse(15, 20e3, time)
    # Find the indices in the time array for the pulse and interpulse periods.
    tstart_p, tend_p = find_time_idx(tt_pulse_bounds, time)
    tstart_ip, tend_ip = find_time_idx(tt_interpulse_bounds, time)
    # Set which period to investigate.
    tstart = tstart_p
    tend = tend_ip
    # Initialize an empty array to store reactants.
    reactants = []
    # Iterate through all reactants in the combustion reaction of interest, append each reactant to the reactants list.
    for keyr in list(reaction_obj.reactants.keys()):
        for i in range(int(reaction_obj.reactants[keyr])):
            reactants.append(keyr)
    # Initialize an empty array to store products.
    products = []
    # Iterate through all products in the combustion reaction of interest, append each product to the products list.
    for keyp in list(reaction_obj.products.keys()):
        for w in range(int(reaction_obj.products[keyp])):
            products.append(keyp)
    # Intialize a placeholder value for the reaction rate coefficient, k_t.
    k_t = 0
    # Check if the partial pathway j is in the reactants, the partial pathway k is in the products, and the reaction type is Arrhenius.
    if j in reactants and k in products and reaction_type == "Arrhenius":
        # If so, k_t is equivalent to the Arrhenius formulation for the time-dependent reaction rate coefficient.
        k_t = Arrhenius(T[tstart:tend], rate_constants, R)
    # Check if the partial pathway j is in the products, the partial pathway k is in the reactants, the combustion reaction
    # of interest is in the reverse reaction rate data file, and the reaction type is Arrhenius.
    elif j in products and k in reactants and reaction_obj.equation in reverse_rate_coeff_data.columns and reaction_type == "Arrhenius":
        # If so, k_t is equal to the reverse reaction rate data.
        k_t = reverse_rate(reverse_rate_coeff_data[tstart:tend], reaction_obj.equation).to_numpy()
        # Initialize an empty list to store the reactants.
        reactants = []
        # Iterate through all products in the combustion reaction of interest, append each product to the reactants list.
        for keyr in list(reaction_obj.products.keys()):
            for i in range(int(reaction_obj.products[keyr])):
                reactants.append(keyr)
        # Initialize an empty list to store the products.
        products = []
        # Iterate through all reactants in the combustion reaction of interest, append each reactant to the products list.
        for keyp in list(reaction_obj.reactants.keys()):
            for j in range(int(reaction_obj.reactants[keyp])):
                products.append(keyp)
    # Initialize a placeholder variable for the time-dependent reaction rate.
    q_t = 1
    # Iterate through each reactant.
    for reactant in reactants:
        # Initialize a string variable that represents the key to access the relevant time-dependent mass fraction data.
        temp_mole_frac_str = "X_" + reactant.lower()
        # Access the time-dependent mass fraction data, removing NaN and infinity values.
        temp_instantaneous_mole_frac = chi_data[temp_mole_frac_str].to_numpy()[tstart:tend]
        temp_instantaneous_mole_frac = np.nan_to_num(temp_instantaneous_mole_frac, nan=0, posinf=0, neginf=0)
        # Compute the time-dependent reaction rate data by converting the mass-fraction data to time-dependent number density.
        q_t *= y2n(T[tstart:tend], P0, temp_instantaneous_mole_frac)
    # Convert the units of the time-dependent reaction rate to 1 / (gmol*cm^3*s) for each reaction type.
    if reaction_type == "three-body-Arrhenius":
        q_t *= (efficiency*1e6*k_t)
    elif reaction_type == "Arrhenius":
        q_t *= (1e3*k_t)
    # Return the time-dependent reaction rate and the start and end bounds for the time array.
    return [q_t, tstart, tend]

"""
A function to solve for the atomic transfer multipliers.
"""
def atom_transfer_multiplier(element: str, reaction: ct.Reaction, pathway: str, third_body: str = None) -> List[float]:
    """
    Parameters:
    -----------
    element : str
        A string representing the element to be tracked in the reaction pathway analysis.
    reaction : ct.Reaction
        An object representing the combustion reaction of interest.
    pathway : str
        A string representing the combustion reaction pathway of interest.
    third_body : str, optional
        A string representing the third-body species in a three-body Arrhenius reaction. Default is None.

    Returns:
    --------
    [vAj, vAk, vAi] : List[float]
        A list of atomic transfer multipliers:
        - vAj : float
            The atomic transfer multiplier for the incident pathway node.
        - vAk : float
            The atomic transfer multiplier for the product pathway node.
        - vAi : float
            The total number of atoms of the element of interest present in the reaction.
    """
    # Intialize placeholder values for atomic transfer multipliers for the incident and product pathway nodes as well as
    # the total number of atoms of the element of interest present in the reaction.
    vAj = 0
    vAk = 0
    vAi = 0
    # Split the pathway into incident and product pathway nodes.
    j, k = pathway.split(" → ")
    # Collect the reactants and products into dictionaries where the key is the reactant/product chemical formula and 
    # the value is the number of reactants/products.
    reactants = reaction.reactants
    products = reaction.products
    # If the reaction type is three body Arrhenius, set the third body efficiency to 1.
    if reaction.reaction_type == "three-body-Arrhenius":
        reactants[third_body] = 1.0
        products[third_body] = 1.0
    # If the incident pathway node is in the reactants and the product pathway node is in the products, add the float-casted
    # number of j reactants to vAj and the float-casted number of k products to vAk.
    if j in reactants and k in products:
        vAj += float(reactants[j])
        vAk += float(products[k])
    # If the incident pathway node is the products and the product pathway node is the reactants, add the float-casted
    # number of k reactants to vAj and the float-casted number of j products to vAk.
    elif k in reactants and j in products:
        vAj += float(reactants[k])
        vAk += float(products[j])
    # Collect a list of reactants in the combustion reaction of interest.
    filtered_reactant_keys = [key for key in reactants.keys() if key is not None]
    for reactant in filtered_reactant_keys:
        # Check if the element of interest is in the reactant.
        if element in reactant:
            # Iterate through each reactant string.
            for idx in range(len(reactant)):
                # Initialize a variable that represents the character at index idx of the reactant string.
                char = reactant[idx]
                # Check if char is the element of interest and if the next idx is not the end of the current reactant.
                if char == element and idx + 1 < len(reactant):
                    # Initialize a variable that represents the next character in the reactant string.
                    next_char = reactant[idx+1]
                    # Check if the next character is a digit.
                    if next_char.isdigit():
                        # If so, cast the character to a float and multiply it by the number of the reactant of interest
                        # in the reaction.
                        vAi += float(next_char)*float(reactants[reactant])
                    else:
                        # If not, multiply the number of reactant of interest in the reaction by 1.0.
                        vAi += 1.0*float(reactants[reactant])
                # Check if char is the element of interest and if the next idx is the end of the current reactant.
                elif char == element and idx + 1 == len(reactant):
                    # If so, multiply the number of reactant of interest in the reaction by 1.0.
                    vAi += 1.0*float(reactants[reactant])
    # Collect a list of products in the combustion reaction of interest.
    filtered_product_keys = [key for key in products.keys() if key is not None]
    for product in filtered_product_keys:
        # Check if the element of interest is in the product.
        if element in product:
            # Iterate through each product string.
            for idx in range(len(product)):
                # Initialize a variable that represents the character at index idx of the product string.
                char = product[idx]
                # Check if char is the element of interest and if the next idx is not the end of the current product.
                if char == element and idx + 1 < len(product):
                    # Initialize a variable that represents the next character in the product string.
                    next_char = product[idx+1]
                    # Check if the next character is a digit.
                    if next_char.isdigit():
                        # If so, cast the character to a float and multiply it by the number of the product of interest
                        # in the reaction.
                        vAi += float(next_char)*float(products[product])
                    else:
                        # If not, multiply the number of product of interest in the reaction by 1.0.
                        vAi += 1.0*float(products[product])
                # Check if char is the element of interest and if the next idx is the end of the current product.
                elif char == element and idx + 1 == len(product):
                    # If so, multiply the number of reactant of interest in the reaction by 1.0.
                    vAi += 1.0*float(products[product])    
    # Return the atomic transfer multipliers.
    return [vAj, vAk, vAi]

"""
A function to solve for the integrated element flux.
"""
def compute_pathway_flux(multiplier_constants: List[float], q_t: np.array, t: np.array, tstart: float, tend: float) -> float:
    """
    Parameters:
    -----------
    multiplier_constants : List[float]
        A list representing the atomic transfer multipliers for the pathway j → k [vAj, vAk, vAi].
    q_t : np.array
        A numpy array representing the time-dependent reaction rate for the combustion reaction of interest.
    t : np.array
        A numpy array representing the time points corresponding to the reaction rate data.
    tstart : float
        A float representing the starting time boundary of calculation.
    tend : float
        A float representing the ending time boundary of calculation.

    Returns:
    --------
    A : float
        The integrated element flux in mol/cm^3.
    """
    # Unpack the atomic transfer multipliers.
    vAj, vAk, vAi = multiplier_constants
    # Compute the aggregate multiplier. 
    multiplier = (vAj*vAk)/vAi
    # Solve for the element flux rate in mol/(cm^3 s).
    Adot = q_t*multiplier
    # Solve for the element flux in mol/(cm^3) by integrating with respect to time.
    A = integrate.trapezoid(Adot, x=t)/(tend - tstart)
    # Return the element flux.
    return A

"""
A function to solve for the time-dependent Arrhenius reaction rate.
"""
def Arrhenius(T: np.array, rate_constants: List[float], R: float) -> np.array:
    """
    Parameters:
    -----------
    T : np.array
        A numpy array representing the time-dependent temperature of the gas mixture.
    rate_constants : List[float]
        A list representing the rate constants for the combustion reaction of interest [A, b, E].
    R : float
        A float representing the universal gas constant.

    Returns:
    --------
    k : np.array
        A numpy array representing the time-dependent reaction rate, k, computed using the Arrhenius equation.
    """
    # Unpack the rate constants for the time-dependent reaction rate.
    A, b, E = rate_constants
    # Solve for the time-dependent reaction rate using the Arrhenius formulation.
    k = A*T**b*np.exp(-E/(R*T))
    # Return the time-dependent reaction rate, k.
    return k

"""
A function to convert mass-fraction to number density.
"""
def y2n(T: np.array, P0: float, y: np.array) -> np.array:
    """
    Parameters:
    -----------
    T : np.array
        A numpy array representing the time-dependent temperature of the gas mixture.
    P0 : float
        A float representing the initial pressure of the gas mixture (in Pascals).
    y : np.array
        A numpy array representing the time-dependent mass-fraction of the combustion species of interest.

    Returns:
    --------
    n : np.array
        A numpy array representing the time-dependent number density of the species (in mol/cm^3).
    """
    # Leveraging the Boltzmann formulation of the ideal gas law, convert mass-fraction to number density.
    n = (y*P0/(1.38e-23*T)/1e6)
    n = n*(1/Avogadro)
    # Return the time-dependent number density.
    return n

"""
A function to convert strings to LaTeX format.
"""
def latexify(species_str: str) -> str:
    """
    Parameters:
    -----------
    species_str : str
        A string representing a species to be plotted in the network diagram.

    Returns:
    --------
    latex : str
        A string representing the LaTeX formatted version of the species.
    """
    # Convert the input string to LaTeX format.
    latex = r"$" + re.sub(r'(\d)', r'_\1', species_str) + "$"
    # Return the LaTeX formatted input string.
    return latex

"""
A function to collect all element flux values into a data structure that associates element flux to a reaction pathway.
"""
def get_flux_dict(reactions: List[ct.Reaction], element_of_interest: str, reverse_rate_coeffs: pd.DataFrame, 
                  T: np.array, R: float, P0: float, chi: pd.DataFrame, cutoff: float) -> Dict[str, np.array]:
    """
    Parameters:
    -----------
    reactions : List[ct.Reaction]
        A list of objects representing the combustion reactions considered in the numerical model.
    element_of_interest : str
        A string representing the element of interest to be tracked.
    reverse_rate_coeffs : pd.DataFrame
        A Pandas DataFrame storing all time-dependent reverse reaction rate coefficients for all combustion reactions considered in the numerical model.
    T : np.array
        A numpy array representing the time-dependent temperature of the gas mixture.
    R : float
        A float representing the universal gas constant.
    P0 : float
        A float representing the initial pressure of the gas mixture.
    chi : pd.DataFrame
        A Pandas DataFrame storing all time-dependent mass fraction data for all combustion reactions in the numerical model.
    cutoff : float
        A user-selected cutoff value for element flux.

    Returns:
    --------
    flux_dict : Dict[str, np.array]
        A dictionary where the key is the reaction pathway (as a string), and the value is the time-dependent element flux (as a numpy array).
    """
    # Print out starting sequence.
    print("===== Starting reaction pathways analysis of the input combustion mechanism =====")
    # Track the start time of the combustion reaction pathways analysis.
    start_time = time.time()
    # Initialize an empty dictionary to store the pathways where the key is the reaction and the value is a list of all possible pathways.
    pathways = {}
    # Iterate through each reaction object.
    for rxn in reactions:
        # Find all pathways possible for the current reaction, given an element of interest.
        temp_pathway = find_pathways(rxn, element_of_interest)
        # Check if the pathway is in the dictionary.
        if temp_pathway is not None:
            # If not, update the keys in the dictionary.
            pathways.update(temp_pathway)
        # Else, pass on the current reaction.
        else:
            pass
    # Initialize an empty dictionary to store the element flux (value) in relation to the pathway (key).
    flux_dict = {}
    # Initialize a variable representing the reaction number being analyzed.
    rxn_counter = 0
    # Initialize a variable representing the total number of reactions tracked in the combustion mechanism.
    total_reactions = len(list(pathways.keys()))
    # Iterate through all reactions in the numerical model.
    for rxn_str in list(pathways.keys()):
        # Iterate the reaction number counter.
        rxn_counter += 1
        # Compute the elapsed time.
        elapsed_time = time.time() - start_time 
        # Estimate remaining time (assuming constant processing time per reaction).
        if rxn_counter > 1:
            avg_time_per_reaction = elapsed_time / rxn_counter
            estimated_remaining_time = avg_time_per_reaction * (total_reactions - rxn_counter)
        else:
            estimated_remaining_time = 0
        # Format and display progress message.
        print(f"Processing Reaction {rxn_counter}/{total_reactions} - Estimated Remaining Time: {estimated_remaining_time / 60:.2f} minutes")
        print(f"Combustion Reaction Pathways Analysis ===== Elapsed Time: {elapsed_time:.2f} seconds =====")
        
        # Collect the reaction string.
        rxn = pathways[rxn_str][0]
        # Generate the pathways of interest.
        temp_pathways = pathways[rxn_str][1]
        # Iterate through all pathways.
        for path in temp_pathways:
            # Collect the rate constants.
            rate_constants = get_rate_constants(rxn)
            # Solve for the reaction rate.
            [q_t, tstart, tend] = solve_reaction_rate(reverse_rate_coeffs, path, rxn, rate_constants, T, R, P0, chi, efficiency=1)
            # Solve for the atomic transfer multipliers.
            multiplier_constants = atom_transfer_multiplier(element_of_interest, rxn, path)
            # Compute element flux.
            flux = compute_pathway_flux(multiplier_constants, q_t, chi["t"].to_numpy()[tstart:tend], tstart, tend)
            # Invoke the element flux cutoff. If the maximum flux is greater than the cutoff, continue.
            if np.max(flux) > cutoff:
                # Check if the pathway of interest is not in the flux dictionary.
                if path not in flux_dict:
                    # If so, add the new key-value pair.
                    flux_dict[path] = flux
                else:
                    # If not, add the flux contribution to the pathway of interest.
                    flux_dict[path] = flux_dict[path] + flux
    # Print out ending sequence.
    print("===== Ending reaction pathways analysis of the input combustion mechanism =====")
    # Return the flux dictionary.
    return flux_dict

def coupled_network_plotter(flux_dict: Dict[str, np.array]):
    """
    A function to plot the coupled plasma-combustion reaction pathways in a network diagram.

    Parameters:
    -----------
    flux_dict : Dict[str, np.array]
        A dictionary representing the element flux for the total PAC chemistry where the key is the pathway of interest and the value is the associated element flux.
    """
    # Generate a list of all pathways tracked.
    paths = list(flux_dict.keys())
    # Construct an empty list representing the intensity of the colormap where the intensity is correspondent to the element flux.
    intensity = []
    # Iterate through all paths tracked.
    for path in paths:
        # Gather and append the element flux to the intensity array.
        temp_flux = flux_dict[path]
        intensity.append(temp_flux)
    # Cast the intensity list to a numpy.array.
    intensity = np.array(intensity)
    # Generate a matplotlib figure.
    fig, ax = plt.subplots(figsize=(7.5, 5))  # Create figure and axis
    G = nx.MultiDiGraph()
    # Iterate through each pathway tracked to add an edge to the NetworkX diagram.
    for i, pathway in enumerate(paths):
        source, target = pathway.split(' → ')
        G.add_edge(source, target, intensity=intensity[i])
    # Draw the layout using Graphviz.
    pos = nx.drawing.nx_pydot.graphviz_layout(G, prog='neato')
    # Draw the nodes on the NetworkX diagram.
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color='white', node_size=100)
    # Color the edges of the NetworkX diagram.
    cmap = plt.cm.turbo
    min_intensity = np.min(intensity)
    max_intensity = np.max(intensity)
    norm = mcolors.LogNorm(vmin=min_intensity, vmax=max_intensity)
    # Iterate through each edge.
    for u, v, d in G.edges(data=True):
        edge_intensity = d['intensity']
        edge_color = cmap(norm(edge_intensity))
        # Draw the edges of the NetworkX diagram.
        nx.draw_networkx_edges(G, pos, ax=ax, edgelist=[(u, v)], arrowsize=12, edge_color=[edge_color], connectionstyle="arc3,rad=0.125")
    # Draw node labels on the NetworkX diagram.
    labels = {node: latexify(node) for node in G.nodes()}
    nx.draw_networkx_labels(G, pos, ax=ax, labels=labels, font_color='black', font_size=8, font_weight='bold')
    # Add the colorbar to represent the element flux in log-scale.
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=mcolors.LogNorm(vmin=min_intensity, vmax=max_intensity))
    sm.set_array([])
    # Create a new axis for the colorbar
    cax = fig.add_axes([0.88, 0.15, 0.03, 0.7]) # Adjust position as necessary
    cbar = plt.colorbar(sm, cax=cax, shrink=0.5)
    cbar.set_label(r'Element Flux ($\frac{mol}{cm^3}$)', labelpad=15)
    # Turn the axis off.
    ax.axis('off')
    # Present the diagram.
    plt.show()

"""
A function to collect the pulse and interpulse bounds for a given pulse number.
"""
def find_pulse_interpulse(pulse: int, PRF: float, time: np.array) -> Tuple[np.array, np.array]:
    """
    Parameters:
    -----------
    pulse : int
        An integer representing the given pulse number.
    PRF : float
        A float representing the pulse-repetition frequency (PRF).
    time : np.array
        A numpy array representing the relevant plasma time.

    Returns:
    --------
    [tt_pulse_bounds, tt_interpulse_bounds] : List[np.array, np.array]
        A list containing two numpy arrays:
        - tt_pulse_bounds : np.array
            A numpy array containing the start and end time of the pulse interval [start, end].
        - tt_interpulse_bounds : np.array
            A numpy array containing the start and end time of the interpulse interval [start, end].
    """
    # Compute the period given the PRF.
    period = 1/PRF
    # Find the initial pulse time boundary.
    begin_of_pulse = (pulse-1)*period
    # Find the ending pulse time boundary.
    end_of_pulse = begin_of_pulse + 1e-7
    # Initialize empty arrays representing the pulse and interpulse arrays.
    tt_pulse = []
    tt_interpulse = []
    # Iterate through the plasma time array to append relevant time values to each array.
    for i in time:
        if i >= begin_of_pulse and i <= end_of_pulse:
            tt_pulse.append(i)
        if i >= end_of_pulse and i <= pulse*period:
            tt_interpulse.append(i)
    # Gather the pulse and interpulse boundary values.
    tt_pulse_bounds = np.array([tt_pulse[0], tt_pulse[-1]])
    tt_interpulse_bounds = np.array([tt_interpulse[0], tt_interpulse[-1]])
    # Return the bounds.
    return [tt_pulse_bounds, tt_interpulse_bounds]


"""
A function to collect the filtered time array boundary indexes given a pulse or interpulse boundary period.
"""
def find_time_idx(tbounds: np.array, filt_time: np.array) -> Tuple[int, int]:
    """
    Parameters:
    -----------
    tbounds : np.array
        A numpy array representing the pulse or interpulse time boundary values [start, end].
    filt_time : np.array
        A numpy array representing the filtered time array for a pulse or interpulse period.

    Returns:
    --------
    [tstart, tend] : List[int]
        A list of two integers:
        - tstart : int
            The index in filt_time where the time exceeds the start boundary in tbounds.
        - tend : int
            The index in filt_time where the time exceeds the end boundary in tbounds.
    """
    # Initialize boundary index values.
    tstart = 0
    tend = 0
    # Iterate through the filtered time to find the index for the starting and ending index.
    for i in range(len(filt_time)):
        if filt_time[i] > tbounds[0]:
            tstart = i
            break
    for j in range(len(filt_time)):
        if filt_time[j] > tbounds[1]:
            tend = j
            break
    # Return the indexes.
    return [tstart, tend]
