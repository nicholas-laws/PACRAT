# Created by: Nicholas Laws
# Date: 2024

## Imports ##
from plasma_chemistry_objs import *
import numpy as np
from scipy.constants import Avogadro
from fuzzywuzzy import process
import matplotlib.pyplot as plt
import networkx as nx
import matplotlib.colors as mcolors
import re
from scipy import integrate
from typing import *
import pandas as pd
import time

"""
A function to read the plasma mechanism in .INP file format.
"""
def read_inp(file_path: str) -> str:
    """
    Parameters:
    -----------
    file_path : str
        The file path to the .INP plasma mechanism file.

    Returns:
    --------
    inp_contents : str
        The contents of the .INP plasma mechanism file as a string.
    """
    # Open and read the .INP plasma mechanism file.
    with open(file_path, 'r') as file:
        inp_contents = file.read()
    # Return the file's contents.
    return inp_contents

"""
A function to find the species in the input .INP plasma mechanism of interest file.
"""
def find_species(inp_contents: str) -> List[str]:
    """
    Parameters:
    -----------
    inp_contents : str
        A string representing the text in the input .INP plasma mechanism file.

    Returns:
    --------
    species : List[str]
        A list of species tracked by the input .INP plasma mechanism file.
    """
    # Split each line of the .INP file.
    lines_inp = inp_contents.splitlines()
    # Intialize start and end placeholder variables.
    start = 0
    end = 0
    # Iterate through each line of the .INP file.
    for i in range(len(lines_inp)):
        temp_line = lines_inp[i]
        # Locate where the species tracked by the plasma mechanism are.
        if "SPECIES" in temp_line:
            # Let the start variable equal to the index where the species tracked by the plasma mechanism begin.
            start = i
            # Iterate through each line after the start of the species tracked by the plasma mechanism begins until the
            # end sequence begins.
            counter = start
            while "END" not in lines_inp[counter]:
                counter += 1
            # Let the end variable equal to the index where the species content end sequence is located.
            end = counter
    # Initialize a variable that represents the sections of the .INP file where the species content is located.
    species_content = lines_inp[start+1:end]
    # Filter the species content for miscellaneous characters.
    filtered_species_content = [element for element in species_content if "#" not in element and element != '']
    # Initialize an empty list that will collect all species tracked by the plasma mechanism of interest.
    species = []
    # Iterate through each element in the filtered species content list.
    for element in filtered_species_content:
        # If there is a space between elements, split the line into separate elements, and append these values to the species list.
        if ' ' in element:
            species.extend(element.split())
        # Otherwise, append the singular element to the species list.
        else:
            species.append(element)
    # For uniformity purposes, search for the electron species and make this element in the species list lowercase.
    if "E" in species:
        index_electron = species.index("E")
        species[index_electron] = "e"
    # Append the ANY_NEUTRAL species to the species list.
    species.append('ANY_NEUTRAL')

    # Return the species list.
    return species

"""
A function to find the reactions in the input .INP plasma mechanism of interest file.
"""
def find_reactions(inp_contents: str) -> Dict[str, List[Union[str, Dict[str, List[str]]]]]:
    """
    Parameters:
    -----------
    inp_contents : str
        A string representing the text in the input .INP plasma mechanism file.

    Returns:
    --------
    reactions : Dict[str, List[Union[str, Dict[str, List[str]]]]]
        A dictionary where the key is the reaction string, and the value is a list containing the reaction rate 
        and a dictionary of species-specific reaction data.
    """
    # Split each line of the .INP file.
    lines_inp = inp_contents.splitlines()
    # Intialize start and end placeholder variables.
    start = 0
    end = 0
    # Iterate through each line of the .INP file.
    for i in range(len(lines_inp)):
        temp_line = lines_inp[i]
        # Locate where the reactions tracked by the plasma mechanism are.
        if "REACTIONS" in temp_line:
            # Let the start variable equal to the index where the reactions tracked by the plasma mechanism begin.
            start = i
            # Iterate through each line after the start of the reactions tracked by the plasma mechanism begins until the
            # end sequence begins.
            counter = start
            while "END" not in lines_inp[counter]:
                counter += 1
            # Let the end variable equal to the index where the reactions content end sequence is located.
            end = counter
    # Initialize a variable that represents the sections of the .INP file where the reactions content is located.
    reactions_content = lines_inp[start+1:end]
    # Filter the species content for miscellaneous characters.
    filtered_reactions_content = [element for element in reactions_content if "=>" in element or "@" in element]
    # Initialize placeholder variables for start and end counters to distinguish where a reaction is located.
    start = 0
    end = 0
    # Initialize a dictionary to collect reactions and rates.
    reactions = {}
    # Iterate through the filtered reaction content.
    for j in range(len(filtered_reactions_content)):
        temp_element = filtered_reactions_content[j]
        # Check if there is a "!" in the current line.
        if "!" in temp_element:
            # If so, split the line at the ! to separate the reaction rate from the reaction.
            split_reaction_str = temp_element.split('!')
            # Initialize a variable that represents the reaction string.
            first_str = split_reaction_str[0].strip()
            # Initialize a variable that represents the rate string.
            second_str = split_reaction_str[1].strip()
            # Let the start variable equal the current index.
            start = j
            # Locate where the end sequence of the reactions content is.
            counter = start + 1
            # Iterate through the filtered reactions content until the end sequence has been initiated.
            while counter < len(filtered_reactions_content) and ("!" not in filtered_reactions_content[counter] or "@" not in filtered_reactions_content[counter]):
                if "END" in filtered_reactions_content[counter]:
                    break
                counter += 1
            end = counter
            # Initialize a list that represents all possible reaction rates for a given reaction. 
            data = [second_str]
            # Initialize a temporary dictionary that will associate all formatted reactions to a formatted reaction rate.
            temp_data = {}
            # Iterate through the start and end content of the reactions tracked by the input plasma mechanism.
            for k in range(start+1,end):
                # Initialize a variable that represents the current line in the reaction content section, stripped of extra spaces.
                temp_data_str = filtered_reactions_content[k].strip()
                # Check if there is an @ character in the current line.
                if "@" in temp_data_str:
                    # Substitute the @ with the species of interest and split the reaction.
                    temp_key = temp_data_str[temp_data_str.index("@")] + temp_data_str[temp_data_str.index("@")+1] 
                    split_str_data = temp_data_str.split('=')
                    split_str_specific = split_str_data[1].split()
                    # Check if the temporary key is not in the data dictionary.
                    if temp_key not in temp_data:
                        # If the key isn't in the temporary data array, set the reaction as the key and the associated reaction rate as the value.
                        temp_data[temp_key] = split_str_specific
            # Append the temporary data to the complete data array.
            data.append(temp_data)
            # Apply the placeholder reaction string as the key and the value as the reaction rates dictionary.
            reactions[first_str] = data
    # Return the reactions dictionary.
    return reactions

"""
A function to create reaction objects to meld seamlessly with the Cantera reaction objects.
"""
def reaction_objects(reactions: Dict[str, List[Union[str, Dict[str, List[str]]]]], species: List[str]) -> List['Reaction']:
    """
    Parameters:
    -----------
    reactions : Dict[str, List[Union[str, Dict[str, List[str]]]]]
        A dictionary representing the reactions tracked by the input .INP plasma mechanism file.
    species : List[str]
        A list representing the species tracked by the input .INP plasma mechanism file.

    Returns:
    --------
    reaction_objs : List['Reaction']
        A list of reaction objects tracked by the input .INP plasma mechanism file.
    """
    # Initialize an empty to hold reaction objects.
    reaction_objs = []
    # Iterate through each reaction stored in the reaction dictionary.
    for key in list(reactions.keys()):
        # Define the reaction equation via the reaction dictionary key.
        eqn = key
        # Define the rate coefficient equations based on the value stored in the reaction dictionary key.
        val = reactions[key]
        rate = val[0]
        conditions = val[1]
        # Initialize two empty lists used to store the reaction equations and the reaction rate coefficients.
        eqns = []
        rates = []
        # Check if conditions are present.
        if conditions:
            # If so, iterate through each conditional reaction rate coefficient for a compounded reaction.
            conditions_keys = list(conditions.keys())
            for key in conditions_keys:
                # Check if the current iteration reaction rate condition equation is within the reaction equation.
                if key in eqn:
                    # Check if the eqns list is empty, if so, fill the list with the respective reaction equation for the conditional reaction rate coefficients.
                    if len(eqns) == 0:
                        eqns = [eqn.replace(key, value) for value in conditions[key]]
                    # If not, iterate through each reaction equation and replace the template with the respective reaction equation variables for the conditional reaction rate coefficient.
                    else:
                        for idx, rxn in enumerate(eqns):
                            temp_eqn = rxn.replace(key, conditions[key][idx])
                            eqns[idx] = temp_eqn
                # Check if the current iteration reaction rate condition is within the rate list.
                elif key in rate:
                    # Check if the rates list is empty, if so, fill the list with the conditional reaction rate coefficients.
                    if len(rates) == 0:
                        rates = [rate.replace(key, value) for value in conditions[key]]
                    # If not, iterate through each conditional reaction rate coefficient and replace the template with the respective reaction equation variables for the conditional reaction rate coefficient.
                    else:
                        for idx, rate in enumerate(rates):
                            temp_rate = rate.replace(key, conditions[key][idx])
                            rates[idx] = temp_rate
            # Define the reactant and product data for the .INP file of interest.
            reactant_data, product_data = format_species(eqns, species)
        # Check if the reaction of interest is not a compound reaction.
        else:
            # If so, set the conditions and compound reaction equations and rate coefficients to None.
            conditions = None
            eqns = None
            rates = None
            reactant_data, product_data = format_species([eqn], species)
        # Format the plasma chemistry objects.
        temp_eqn_obj = Equation(eqn, eqns)
        temp_rate_obj = Rate(rate, rates)
        temp_species_obj = Species(reactant_data, product_data)
        temp_reaction_obj = Reaction(temp_eqn_obj, temp_rate_obj, temp_species_obj)
        # Append the parent plasma chemistry object to the overall list.
        reaction_objs.append(temp_reaction_obj)
    # Return the reaction objects list.
    return reaction_objs


"""
A function to collect the number of different reactants and products in each reaction, formatting them in a dictionary 
where the key is the reactant/product species and the value is the number of occurrences for each species.
"""
def format_species(rxns: List[str], species: List[str]) -> List[List[Union[Dict[str, int], List[Dict[str, int]]]]]:
    """
    Parameters:
    -----------
    rxns : List[str]
        A list representing the reactions tracked by the input .INP plasma mechanism file.
    species : List[str]
        A list representing the species tracked by the input .INP plasma mechanism file.

    Returns:
    --------
    [reactant_data, product_data] : List[List[Union[Dict[str, int], List[Dict[str, int]]]]]
        A list containing two lists: one for reactant data and one for product data. Each list contains dictionaries 
        where keys are species names and values are the count of each species in the reaction.
    """
    # Initialize empty array for reactant and product data for each species.
    reactant_data = []
    product_data = []
    # Iterate through the reactions of interest.
    for rxn in rxns:
        # For each reaction, initialize empty dictionaries for reactant and product data.
        temp_reactant = {}
        temp_product = {}
        # Iterate through each species tracked by the input .INP plasma mechanism of interest.
        for s in species:
            # Check if the iterated species is in the reaction.
            if s in rxn:
                # If so, split the iterated reaction into a string array.
                components = rxn.split()
                # Find the index of the forward-directioned arrow.
                index = components.index('=>')
                # Count how many of the iterated species are in the reactant and product sides of the iterated reaction, respectively.
                reactant_count = sum(1 for component in components[:index] if component.strip() == s)
                product_count = sum(1 for component in components[index+1:] if component.strip() == s)
                # Check if the number of reactants and products are greater than 0, if so, add the reactant/product and associated count to the respective dictionary.
                if reactant_count != 0:
                    temp_reactant[s] = reactant_count
                if product_count != 0:
                    temp_product[s] = product_count
        # Append the temporary reactant and product counts to the overall lists.
        reactant_data.append(temp_reactant)
        product_data.append(temp_product)
    # Return the reactant and product data.
    return [reactant_data, product_data]

"""
A function to identify all reaction pathways present in a plasma reaction of interest.
"""
def pathways_identifier(reaction: 'Reaction', element_of_interest: str) -> Dict[str, List[List[Union[str, Dict[str, int]]]]]:
    """
    Parameters:
    -----------
    reaction : Reaction
        An object representing a reaction tracked in the input .INP plasma mechanism.
    element_of_interest : str
        A string representing the element to be tracked in the reaction pathway analysis.

    Returns:
    --------
    pathways : Dict[str, List[List[Union[str, Dict[str, int]]]]]
        A dictionary where the key is the reaction equation and the value is the list of relevant reaction pathways data identified in the reaction. 
        Each list consists of the pathway string, a dictionary of reactants, and a dictionary of products.
    """
    # Initialize an empty array representing the storage for reaction pathways present in the input .INP plasma mechanism of interest.
    pathways = {}
    # Initialize a mutable temporary pathway value. 
    temp_pathway = 0
    # Access the reactants, products, and reactions objects.
    reactants = reaction.species.reactants
    products = reaction.species.products
    reactions = reaction.equation.reactions
    # Check if the reactions object is not empty.
    if reactions != None:
        # Iterate through all associated reaction strings for the reaction object of interest.
        for idx in range(len(reactions)):
            # Iterate through each reactant in the reaction of interest.
            for reactant in list(reactants[idx].keys()):
                    # Iterate through each character of the iterated reactant.
                    for char in reactant:
                        # Check if the iterated character of the reactant is not a number.
                        if char.isdigit() == False:
                            # Iterate through each product in the reaction of interest.
                            for product in list(products[idx].keys()):
                                # Check if the iterated character is not a number and is in the product.
                                if char in product:
                                    # Locate the indices of the character present in both the iterated reactant and product.
                                    index_char_reactant = reactant.index(char)
                                    index_char_product = product.index(char)
                                    # Check if the reactant/product is not a neutral, the reactant is not the product, the user-selected element of interest to be tracked is in both the reactant and product, and that the index before the corresponding character in both the reactant/product is not a "(".
                                    if reactant != "ANY_NEUTRAL" and product != "ANY_NEUTRAL" and reactant != product and element_of_interest in reactant and element_of_interest in product and reactant[index_char_reactant-1] != "(" and product[index_char_product-1] != "(":
                                        # Check if the length of the product is greater than 2 and that the number of the tracked element of interest in the selected product is equal to 1.
                                        if len(product) > 2 and sum(1 for element in product if element_of_interest in element) == 1:
                                            # Define the pathway.
                                            temp_pathway = f"{reactant} → {product}"
                                        # If not, pass.
                                        else:
                                            pass
                                    # Check if the iterated reaction is not already in the pathways dictionary and if the temporary pathway has formatted properly.
                                    if reactions[idx] not in pathways and temp_pathway == f"{reactant} → {product}":
                                        # Insert the first iterated reaction into the pathways dictionary where the key is the iterated reaction and the value is a list where the first element is the pathway, the second element is the reactants associated with the reaction, and the third element is the products associated with the reaction.
                                        pathways[reactions[idx]] = [[temp_pathway, reactants[idx], products[idx]]]
                                    # Check if the iterated reaction is already in the pathways dictionary, the temporary pathway has formatted properly, and the data list associated with the pathway is not already in the dictionary.
                                    elif reactions[idx] in pathways and temp_pathway == f"{reactant} → {product}" and [temp_pathway, reactants[idx], products[idx]] not in pathways[reactions[idx]]:
                                        # Extend the pathways dictionary for the iterated reaction with the Nth pathway and the associated reactants/products data. 
                                        pathways[reactions[idx]].extend([[temp_pathway, reactants[idx], products[idx]]])
    # Check if the reaction object is not a coupled reaction.
    elif reactions == None:
        # If so, initialize the reaction equation as the template.
        eqn = reaction.equation.template
        # Iterate through each reactant in the reaction of interest.
        for reactant in list(reactants[0].keys()):
                # Iterate through each character of the iterated reactant.
                for char in reactant:
                    # Check if the iterated character of the reactant is not a number.
                    if char.isdigit() == False:
                        # Iterate through each product in the reaction of interest.
                        for product in list(products[0].keys()):
                            # Check if the iterated character is not a number and is in the product.
                            if char in product:
                                # Locate the indices of the character present in both the iterated reactant and product.
                                index_char_reactant = reactant.index(char)
                                index_char_product = product.index(char)
                                # Check if the reactant/product is not a neutral, the reactant is not the product, the user-selected element of interest to be tracked is in both the reactant and product, and that the index before the corresponding character in both the reactant/product is not a "(".
                                if reactant != "ANY_NEUTRAL" and product != "ANY_NEUTRAL" and reactant != product and element_of_interest in reactant and element_of_interest in product and reactant[index_char_reactant-1] != "(" and product[index_char_product-1] != "(":
                                    # Define the pathway.
                                    temp_pathway = f"{reactant} → {product}"
                                # Check if the iterated reaction is not already in the pathways dictionary and if the temporary pathway has formatted properly.
                                if eqn not in pathways and temp_pathway == f"{reactant} → {product}":
                                    # Insert the first iterated reaction into the pathways dictionary where the key is the iterated reaction and the value is a list where the first element is the pathway, the second element is the reactants associated with the reaction, and the third element is the products associated with the reaction.
                                    pathways[eqn] = [[temp_pathway, reactants[0], products[0]]]
                                # Check if the iterated reaction is already in the pathways dictionary, the temporary pathway has formatted properly, and the data list associated with the pathway is not already in the dictionary.
                                elif eqn in pathways and temp_pathway == f"{reactant} → {product}" and [temp_pathway, reactants[0], products[0]] not in pathways[eqn]:
                                    # Extend the pathways dictionary for the iterated reaction with the Nth pathway and the associated reactants/products data. 
                                    pathways[eqn].extend([[temp_pathway, reactants[0], products[0]]])
    # Return the pathways for the reaction of interest.
    return pathways

"""
A function to identify all reaction pathways present in the input .INP plasma mechanism of interest.
"""
def collect_pathways(reactions: List['Reaction'], element_of_interest: str) -> Dict[str, List[List[Union[str, Dict[str, int]]]]]:
    """
    Parameters:
    -----------
    reactions : List[Reaction]
        A list of reaction objects representing the reactions tracked in the input .INP plasma mechanism.
        
    element_of_interest : str
        A string representing the element to be tracked in the reaction pathway analysis.

    Returns:
    --------
    pathways : Dict[str, List[List[Union[str, Dict[str, int]]]]]
        A dictionary where the key is the reaction equation (a string), and the value is a list of pathways.
        Each pathway includes the pathway string (reactant → product), a dictionary of reactants with species names as keys
        and their counts as values, and a dictionary of products with species names as keys and their counts as values.
    """
    # Initialize an empty pathways dictionary.
    pathways = {}
    # Iterate through each reaction in the list of reaction objects.
    for rxn in reactions:
        # Identify the pathways within the reaction of interest.
        temp_pathways = pathways_identifier(rxn, element_of_interest)
        # Update the overall pathways dictionary.
        pathways.update(temp_pathways)
    # Return the pathways for the input .INP plasma mechanism of interest.
    return pathways

"""
A function that iterates through a dataset, replacing NaN and infinity values with the most recent non-NaN value.
"""
def replace_nan_with_recent(array: np.array) -> np.array:
    """
    Parameters:
    -----------
    array : np.array
        A numpy array of data with arbitrary typing, which may contain NaN or infinity values.

    Returns:
    --------
    array : np.array
        The input array with NaN and infinity values replaced by the most recent non-NaN value.
    """
    # Initialize a placeholder value for the most recent non-NaN value.
    recent_non_nan = None
    # Iterate through the dataset.
    for i in range(len(array)):
        val = array[i]
        # Check if the current value is not NaN or infinity.
        if not np.isnan(val) and not np.isinf(val):
            # If so, define the placeholder variable as val.
            recent_non_nan = val
        # Check if the first value in array is NaN or infinity.
        elif recent_non_nan is None:
            # If so, define the first value in array as 0.
            array[i] = 0
        # Check if the most recent non-NaN value is not empty.
        elif recent_non_nan is not None:
            # If so, define the current index array to recent_non_nan.
            array[i] = recent_non_nan
    # Return the updated dataset.
    return array

"""
A function that iterates through the species tracked in the input .INP plasma mechanism to locate and collect all neutral species.
"""
def find_neutrals(species: List[str]) -> List[str]:
    """
    Parameters:
    -----------
    species : List[str]
        A list of strings representing the species tracked by the input .INP plasma mechanism.

    Returns:
    --------
    neutrals : List[str]
        A list of strings representing all neutral species.
    """
    # Initialize an empty array representing the storage for neutral species strings.
    neutrals = []
    # Iterate through each species.
    for s in species:
        # Check if the species is not charged, is not the ANY_NEUTRAL species, and is not an electron.
        if "^+" not in s and s != "ANY_NEUTRAL" and s != "e":
            # If so, append the species to the neutrals list.
            neutrals.append(s)
    # Return the list of neutrals.
    return neutrals


"""
A function that locates the start and end times, as well as the associated indices, for the pulse period of a given pulse number.
"""
def find_pulse_range(time: np.array, pulse_num: int, PRF: float, sigma: float, EN: np.array) -> List[Union[float, float, int, int]]:
    """
    Parameters:
    -----------
    time : np.array
        A numpy array representing the relevant plasma time.
    pulse_num : int
        An integer representing the given pulse number.
    PRF : float
        A float representing the pulse-repetition frequency (PRF).
    sigma : float
        A float representing the standard deviation of the Gaussian plasma pulse.
    EN : np.array
        A numpy array representing the reduced electric field corresponding to time.

    Returns:
    --------
    [tp_start, tp_end, tp_start_idx, tp_end_idx] : List[Union[float, float, int, int]]
        A list containing the following:
        - tp_start : float
            The start time of the pulse period, offset by 3 times the standard deviation.
        - tp_end : float
            The end time of the pulse period, offset by 3 times the standard deviation.
        - tp_start_idx : int
            The index in the time array corresponding to tp_start.
        - tp_end_idx : int
            The index in the time array corresponding to tp_end.
    """
    # Solve for the pulse and interpulse boundaries.
    [tt_pulse_bounds, tt_interpulse_bounds] = find_pulse_interpulse(pulse_num, PRF, time)
    [tstart, tend] = find_time_idx(tt_pulse_bounds, time)
    # Locate the relevant portions of the pulse period of interest in the time and EN arrays.
    tt_loc = time[tstart:tend]
    EN_loc = EN[tstart:tend]
    # Offset the start and ends of the pulse period by 3 times the standard deviation on either side.
    tp_start = tt_loc[np.argmax(EN_loc)] - 3*sigma
    tp_end = tt_loc[np.argmax(EN_loc)] + 3*sigma
    # Iterate through the time array.
    for idx_s, i in enumerate(time):
        # Check if the time array value is greater than the selected tp_start value.
        if i > tp_start:
            # If so, define tp_start_idx as the iterated index in the time array, then break.
            tp_start_idx = idx_s
            break
    # Iterate through the time array.
    for idx_e, j in enumerate(time):
        # Check if the time array value is greater than the selected tp_end value.
        if j > tp_end:
            # If so, define tp_end_idx as the iterated index in the time array, then break.
            tp_end_idx = idx_e
            break
    # Return the start and end times, as well as the associated indices, for the pulse period of a given pulse number.
    return [tp_start, tp_end, tp_start_idx, tp_end_idx]

"""
A function that locates the start and end times, as well as the associated indices, for the interpulse period of a given pulse number.
"""
def find_interpulse_range(time: np.array, pulse_num: int, PRF: float, sigma: float, EN: np.array) -> List[Union[float, float, int, int]]:
    """
    Parameters:
    -----------
    time : np.array
        A numpy array representing the relevant plasma time.
    pulse_num : int
        An integer representing the given pulse number.
    PRF : float
        A float representing the pulse-repetition frequency (PRF).
    sigma : float
        A float representing the standard deviation of the Gaussian plasma pulse.
    EN : np.array
        A numpy array representing the reduced electric field corresponding to time.

    Returns:
    --------
    [tip_start, tip_end, tip_start_idx, tip_end_idx] : List[Union[float, float, int, int]]
        A list containing the following:
        - tip_start : float
            The start time of the interpulse period, offset by 3 times the standard deviation.
        - tip_end : float
            The end time of the interpulse period, offset by 3 times the standard deviation.
        - tip_start_idx : int
            The index in the time array corresponding to tip_start.
        - tip_end_idx : int
            The index in the time array corresponding to tip_end.
    """
    # Define the pulse period start and end times, as well as associated time indices for the selected pulse number.
    [tp_start, tp_end, tp_start_idx, tp_end_idx] = find_pulse_range(time, pulse_num, PRF, sigma, EN)
    # Define the pulse period start and end times, as well as associated time indices for the selected pulse number + 1.
    [tp_start_next, tp_end_next, tp_start_idx_next, tp_end_idx_next] = find_pulse_range(time, pulse_num+1, PRF, sigma, EN)
    # Using the boundaries between the selected pulse number and the next pulse, define the interpulse period start and end times, as well as associated time indices.
    tip_start = tp_end
    tip_end = tp_start_next
    tip_start_idx = tp_end_idx
    tip_end_idx = tp_start_idx_next
    # Return the start and end times, as well as the associated indices, for the interpulse period of a given pulse number. 
    return [tip_start, tip_end, tip_start_idx, tip_end_idx]

"""
A function that solves for the reaction rate of a given plasma reaction.
"""
def solve_reaction_rate(rate_str: str, reactants: List[str], plasma_data: pd.DataFrame, neutrals: List[str], 
                        element_of_interest: str, plot_boolean: bool, sigma: float) -> List[Union[np.array, int, int]]:
    """
    Parameters:
    -----------
    rate_str : str
        A string representing the reaction rate coefficient for a given plasma reaction of interest.
    reactants : List[str]
        A list of reactant species corresponding to the plasma reaction of interest.
    plasma_data : pd.DataFrame
        A pandas DataFrame storing all time-dependent number densities and reaction rate coefficients.
    neutrals : List[str]
        A list of neutral species tracked by the input .INP plasma mechanism.
    element_of_interest : str
        A string representing the element to be tracked in the reaction pathway analysis.
    plot_boolean : bool
        A boolean indicating whether the pulse (True) or interpulse (False) period should be considered.
    sigma : float
        A float representing the standard deviation of the Gaussian plasma pulse.

    Returns:
    --------
    [q_t, tstart, tend] : List[Union[np.array, int, int]]
        A list containing:
        - q_t : np.array
            The time-dependent reaction rate for the plasma reaction in units of (mol/cm^3 * s).
        - tstart : int
            The start index in the time array for the time interval considered in the calculation.
        - tend : int
            The end index in the time array for the time interval considered in the calculation.
    """
    # By leveraging the input rate_str, find the closest matching column in plasma_data for the reaction rate coefficient.
    closest_match, similarity_score = process.extractOne(rate_str, plasma_data.columns.tolist())
    # Define the time corresponding to the ZDPlasKin solver.
    time = plasma_data["PDtt"].to_numpy()
    # Define the time-dependent reduced-electric field.
    EN = plasma_data["PDtt"].to_numpy()
    # Define the time corresponding to the Cantera solver.
    time_ct = plasma_data["PDtt_ct"].to_numpy()
    # Check whether plot_boolean is True or False.
    if plot_boolean == True:
        # If True, consider the pulse period.
        [tstart_val, tend_val, tstart, tend] = find_pulse_range(time, 15, 20e3, sigma, EN)
    elif plot_boolean == False:
        # If False, consider the interpulse period.
        [tstart_val, tend_val, tstart, tend] = find_interpulse_range(time, 15, 20e3, sigma, EN)
    # Check if the ANY_NEUTRAL wildcard species is not in the reactants.
    if "ANY_NEUTRAL" not in reactants:
        # If not, collect the time-dependent reaction rate coefficient data between the start and end times.
        rate_coeff_data = plasma_data[closest_match].to_numpy()[tstart:tend]
        # Replace all NaN and infinity values in the time-dependent reaction rate coefficient data.
        k_t = replace_nan_with_recent(rate_coeff_data)
    # Check if the ANY_NEUTRAL wildcard species is in the reactants.
    elif "ANY_NEUTRAL" in reactants:
        # If so, gather the time-dependent electron temperature, replacing all NaN and infinity values.
        Te = replace_nan_with_recent(plasma_data["Te"].to_numpy()[tstart:tend])
        # Hard-coded time-dependent reaction rate coefficient.
        k_t = (6e-27)*((300/Te)**1.5)
    # Initialize a placeholder variable for the time-dependent reaction rate.
    q_t = 1
    # Iterate through all reactants.
    for r in reactants:
        # Check if the iterated reactant is not the ANY_NEUTRAL wildcard species.
        if r != "ANY_NEUTRAL":
            # If so, find the closest matching number density for the iterated reactant in plasma_data, replacing NaN and infinity values.
            closest_match_reactant, similarity_score_reactant = process.extractOne("n_" + r, plasma_data.columns.tolist())
            temp_n_data = plasma_data[closest_match_reactant].to_numpy()[tstart:tend]
            temp_n_data = replace_nan_with_recent(temp_n_data)
            temp_n = temp_n_data
            # Multiply the overall time-dependent reaction rate by the number density for the iterated reactant.
            q_t *= temp_n
        # Check if the iterated reactant is the ANY_NEUTRAL wildcard specis.
        elif r == "ANY_NEUTRAL":
            # Initialize a placeholder variable for the time_dependent number density for all neutral species.
            temp_n = 0
            # Iterate through all neutral species.
            for s in neutrals:
                # Find the closest matching number density for the iterated neutral reactant in plasma_data, replacing NaN and infinity values.
                closest_match_reactant, similarity_score_reactant = process.extractOne("n_" + s, plasma_data.columns.tolist())
                temp_n_data = plasma_data[closest_match_reactant].to_numpy()[tstart:tend]
                temp_n_data = replace_nan_with_recent(temp_n_data)
                # Add the number density for the iterated neutral species to the total number density for all neutral specis.
                temp_n += temp_n_data
            # Multiply the overall time-dependent number density for neutral species by the overall time-dependent reaction rate.
            q_t *= temp_n
    # Convert the time-dependent reaction rate to units of (mol/cm^3 * s).
    q_t = (k_t*q_t)/Avogadro
    # Return the time-dependent reaction rate for the plasma reaction of interest and the start and end times considered in the calculation.
    return [q_t, tstart, tend]

"""
A function to solve for the atomic transfer multipliers.
"""
def atom_transfer_multiplier(reactants: Dict[str, float], products: Dict[str, float], element: str, pathway: str) -> List[float]:
    """
    Parameters:
    -----------
    reactants : Dict[str, float]
        A dictionary representing the reactants and their respective counts in the plasma reaction of interest.
    products : Dict[str, float]
        A dictionary representing the products and their respective counts in the plasma reaction of interest.
    element : str
        A string representing the element to be tracked in the reaction pathway analysis.
    pathway : str
        A string representing the combustion reaction pathway of interest.

    Returns:
    --------
    [vAj, vAk, vAi] : List[float]
        A list of atomic transfer multipliers:
        - vAj : float
            The atomic transfer multiplier for the incident node of the reaction pathway.
        - vAk : float
            The atomic transfer multiplier for the product node of the reaction pathway.
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
    # Check if the incident pathway node is in the reactants and the product pathway node is in the products.
    if j in reactants and k in products:
        # Initialize placeholder multipliers for the reactants and products.
        reactant_multiplier = 1
        product_multiplier = 1
        # Iterate through the characters of the indicent reactant node.
        for char_r_idx in range(len(j)):
            # Initialize a temporary variable representing the iterated character.
            char_r = j[char_r_idx]
            # Check if the character is a number and the character prior to the current character is the element of interest.
            if char_r.isdigit() and j[char_r_idx-1] == element:
                # If so, define the reactant multiplier as the float-casted iterated digit character.
                reactant_multiplier = float(char_r)
        # Iterate through the characters of the product pathway node.
        for char_p_idx in range(len(k)):
            # Initialize a temporary variable representing the iterated character.
            char_p = k[char_p_idx]
            # Check if the character is a number and the character prior to the current character is the element of interest.
            if char_p.isdigit() and k[char_p_idx-1] == element:
                # If so, define the product multiplier as the float-casted iterated digit character.
                product_multiplier = float(char_p)
        # Sum the product of the number of reactants/products and the reactant/product multiplier for each atomic transfer coefficient.
        vAj += float(reactants[j])*reactant_multiplier
        vAk += float(products[k])*product_multiplier
    # Check if the product pathway node is in the reactants and the incident pathway node is in the products.
    elif k in reactants and j in products:
        # Initialize placeholder multipliers for the reactants and products.
        reactant_multiplier = 1
        product_multiplier = 1
        # Iterate through the characters of the product pathway node.
        for char_p_idx in range(len(k)):
            # Initialize a temporary variable representing the iterated character.
            char_p = k[char_p_idx]
            # Check if the character is a number and the character prior to the current character is the element of interest.
            if char_p.isdigit() and k[char_p_idx-1] == element:
                # If so, define the product multiplier as the float-casted iterated digit character.
                reactant_multiplier = float(char_p)
        # Iterate through the characters of the indicent reactant node.
        for char_r_idx in range(len(j)):
            # Initialize a temporary variable representing the iterated character.
            char_r = j[char_r_idx]
            # Check if the character is a number and the character prior to the current character is the element of interest.
            if char_r.isdigit() and j[char_r_idx-1] == element:
                # If so, define the product multiplier as the float-casted iterated digit character.
                product_multiplier = float(char_r)    
        # Sum the product of the number of reactants/products and the reactant/product multiplier for each atomic transfer coefficient.
        vAj += float(reactants[k])*reactant_multiplier
        vAk += float(products[j])*product_multiplier
    # Initialize a list of reactants for the plasma reaction of interest.
    filtered_reactant_keys = [key for key in reactants.keys() if key is not None]
    # Iterate through the reactants.
    for reactant in filtered_reactant_keys:
        # Check if the element of interest is in the reactant. 
        if element in reactant:
            # If so, iterate through the characters of the reactant.
            for idx in range(len(reactant)):
                # Initialize a temporary variable representing the iterated character.
                char = reactant[idx]
                # Check if the iterated character is the element of interest and if the iterated character is not the final character.
                if char == element and idx + 1 < len(reactant):
                    # Initialize a variable representing the next character.
                    next_char = reactant[idx+1]
                    # Check if the next character is a number.
                    if next_char.isdigit():
                        # If so, add the product of the float-casted character and the float-casted number of the reactant of interest to vAi.
                        vAi += float(next_char)*float(reactants[reactant])
                    else:
                        # If not, add the float-casted number of the reactant of interest to vAi.
                        vAi += 1.0*float(reactants[reactant])
                # Check if the iterated character is the element of interest and if the iterated character is the final character.
                elif char == element and idx + 1 == len(reactant):
                    # If so, add the float-casted number of the reactant of interest to vAi.
                    vAi += 1.0*float(reactants[reactant])
    # Initialize a list of products for the plasma reaction of interest.
    filtered_product_keys = [key for key in products.keys() if key is not None]
    # Iterate through the products.
    for product in filtered_product_keys:
        # Check if the element of interest is in the product. 
        if element in product:
            # If so, iterate through the characters of the products.
            for idx in range(len(product)):
                # Initialize a temporary variable representing the iterated character.
                char = product[idx]
                # Check if the iterated character is the element of interest and if the iterated character is not the final character.
                if char == element and idx + 1 < len(product):
                    # Initialize a variable representing the next character.
                    next_char = product[idx+1]
                    # Check if the next character is a number.
                    if next_char.isdigit():
                        # If so, add the product of the float-casted character and the float-casted number of the reactant of interest to vAi.
                        vAi += float(next_char)*float(products[product])
                    else:
                        # If not, add the float-casted number of the reactant of interest to vAi.
                        vAi += 1.0*float(products[product])
                # Check if the iterated character is the element of interest and if the iterated character is the final character.
                elif char == element and idx + 1 == len(product):
                    # If so, add the float-casted number of the reactant of interest to vAi.
                    vAi += 1.0*float(products[product])
    # Return the atomic transfer coefficients.     
    return [vAj, vAk, vAi]

"""
A function to solve for the integrated, time-dependent element flux.
"""
def compute_pathway_flux(multiplier_constants: List[float], q_t: np.array, t: np.array, tstart: float, tend: float) -> float:
    """
    Parameters:
    -----------
    multiplier_constants : List[float]
        A list of floats representing the atomic transfer coefficients [vAj, vAk, vAi].
    q_t : np.array
        A numpy array representing the time-dependent reaction rate for a plasma reaction of interest.
    t : np.array
        A numpy array representing the time corresponding to the plasma reaction of interest.
    tstart : float
        A float representing the start value in the time array for the time interval considered in the calculation.
    tend : float
        A float representing the end value in the time array for the time interval considered in the calculation.

    Returns:
    --------
    A : float
        The integrated element flux for the plasma reaction of interest.
    """
    # Unpack the atomic transfer coefficients.
    vAj, vAk, vAi = multiplier_constants
    # Initialize a value representing the atomic transfer multiplier.
    multiplier = (vAj*vAk)/vAi
    # Solve for the element flux rate.
    Adot = q_t*multiplier
    # Solve for the element flux by integrating the element flux rate using trapezoidal Riemann sum.
    A = integrate.trapezoid(Adot, x=t)/(tend - tstart)
    # Return the element flux.
    return A

"""
A function to locate the reaction rate coefficient string utilized in the input plasma data.
"""
def find_rate_str(reaction_objs: List['Reaction'], reaction: 'Reaction') -> str:
    """
    Parameters:
    -----------
    reaction_objs : List[Reaction]
        A list of reaction objects representing the reactions tracked in the input .INP plasma mechanism.
    reaction : Reaction
        An object representing the plasma reaction of interest.

    Returns:
    --------
    rate_str : str
        A string representing the reaction rate coefficient for the given reaction. If the reaction is a coupled reaction,
        the corresponding rate coefficient for the matching reaction is returned.
    """
    # Initialize a placeholder variable for the reaction rate coefficient string.
    rate_str = 0
    # Iterate through all reaction objects.
    for reaction_obj in reaction_objs:
        # Initialize a variable representing the plasma reaction's equation
        temp_rxn_template = reaction_obj.equation.template
        # Check if the input reaction is the reaction template and that the reaction template is not a coupled reaction.
        if reaction == temp_rxn_template and "@" not in temp_rxn_template:
            # If so, access the reaction rate template.
            rate_str = reaction_obj.rate.template
            # Check if the string "BOLSIG" is not in the reaction rate coefficient string.
            if "BOLSIG" not in rate_str:
                # If so, define the reaction rate coefficient string as the reaction.
                rate_str = reaction
            # Return the selected reaction rate coefficient string.
            return rate_str
        # Check if the reaction rate object is a coupled reaction.
        elif reaction_obj.equation.reactions != None:
            # Define the list of reactions and associated reaction rate coefficient strings for the coupled plasma reaction.
            temp_rxn_strs = reaction_obj.equation.reactions
            temp_rate_strs = reaction_obj.rate.rates
            # Check if the length of the reaction rate coefficient strings equals 0.
            if len(temp_rate_strs) == 0:
                # If so, all reactions in the coupled plasma reaction have the same reaction rate; thus, define the reaction rate coefficient string as the template.
                rate_str = reaction_obj.rate.template
                # Check if the string "BOLSIG" is not in the reaction rate coefficient string.
                if "BOLSIG" not in rate_str:
                    # If so, define the reaction rate coefficient string as the reaction. 
                    rate_str = reaction
                # Return the selected reaction rate coefficient string.
                return rate_str
            else:
                # If not, iterate through the reactions in the coupled plasma reaction.
                for i in range(len(temp_rxn_strs)):
                    # Initialize a temporary variable representing the iterated reaction.
                    temp_rxn = temp_rxn_strs[i]
                    # Check if the temporary reaction is the input reaction.
                    if temp_rxn == reaction:
                        # If so, define the reaction rate coefficient string as the current iterated index. 
                        rate_str = temp_rate_strs[i]
                        # Check if the string "BOLSIG" is not in the reaction rate coefficient string.
                        if "BOLSIG" not in rate_str:
                            # If so, define the reaction rate coefficient string as the reaction. 
                            rate_str = reaction
                        # Return the selected reaction rate coefficient string.
                        return rate_str

"""
A function to collect all element flux values into a data structure that associates element flux to a reaction pathway for the plasma mechanism of interest.
"""
def get_plasma_flux_dict(pathways_dict: Dict[str, List[List[Union[str, Dict[str, int]]]]], reaction_objs: List['Reaction'], 
                         element_of_interest: str, plasma_data: pd.DataFrame, cutoff: float, neutrals: List[str], 
                         plot_boolean: bool, sigma: float) -> Dict[str, np.array]:
    """
    Parameters:
    -----------
    pathways_dict : Dict[str, List[List[Union[str, Dict[str, int]]]]]
        A dictionary where the key is the reaction equation (a string), and the value is a list of pathways.
        Each pathway includes the pathway string (reactant → product), a dictionary of reactants with species names as keys
        and their counts as values, and a dictionary of products with species names as keys and their counts as values.
    reaction_objs : List[Reaction]
        A list of reaction objects representing the reactions tracked in the input .INP plasma mechanism.
    element_of_interest : str
        A string representing the element to be tracked in the reaction pathway analysis.
    plasma_data : pd.DataFrame
        A pandas DataFrame storing all time-dependent number densities and reaction rate coefficients.
    cutoff : float
        A float representing the cutoff value for element flux.
    neutrals : List[str]
        A list of neutral species tracked by the input .INP plasma mechanism.
    plot_boolean : bool
        A boolean indicating whether the pulse (True) or interpulse (False) period should be considered.
    sigma : float
        A float representing the standard deviation of the Gaussian plasma pulse.

    Returns:
    --------
    flux_dict : Dict[str, np.array]
        A dictionary where the key is the pathway string and the value is the computed time-dependent element flux 
        for each reaction pathway that exceeds the cutoff value.
    """
    # Print out starting sequence.
    print("===== Starting reaction pathways analysis of the input plasma mechanism =====")
    # Track the start time of the combustion reaction pathways analysis.
    start_time = time.time()
    # Initialize an empty dictionary to store all element flux values associated with pathways.
    flux_dict = {}
    # Initialize a variable representing the reaction number being analyzed.
    rxn_counter = 0
    # Initialize a variable representing the total number of reactions tracked in the combustion mechanism.
    total_reactions = len(list(pathways_dict.keys()))
    # Initialize a variable that represents the reactions tracked by the input .INP plasma mechanism.
    reactions = list(pathways_dict.keys())
    # Iterate through all reactions.
    for rxn in reactions:
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
        print(f"Plasma Reaction Pathways Analysis ===== Elapsed Time: {elapsed_time:.2f} seconds =====")
        # Find the rate string corresponding to the iterated reaction.
        temp_rate_str = find_rate_str(reaction_objs, rxn)
        # Iterate through the pathway data associated with the iterated reaction of interest.
        for pathway_arr in pathways_dict[rxn]:
            # Unpack the pathway data: pathway, reactants dictionary, and products dictionary.
            pathway = pathway_arr[0]
            reactants = pathway_arr[1]
            products = pathway_arr[2]
            # Solve for the time-dependent reaction rate coefficient.
            [q_t, tstart, tend] = solve_reaction_rate(temp_rate_str, reactants, plasma_data, neutrals, element_of_interest, plot_boolean, sigma)
            # Solve for the atomic transfer coefficients.
            multiplier_constants = atom_transfer_multiplier(reactants, products, element_of_interest, pathway)
            # Compute the element flux for the iterated pathway.
            flux = compute_pathway_flux(multiplier_constants, q_t, plasma_data["PDtt"].to_numpy()[tstart:tend], tstart, tend)
            # Check if the maximum element flux is greater than the cutoff value.
            if np.max(flux) > cutoff:
                # If so, check if the pathway is not already in the overall element flux dictionary.
                if pathway not in flux_dict:
                    # If so, initialize a new key-value pair in flux_dict, where the key is the pathway and the value is the time-dependent element flux.
                    flux_dict[pathway] = flux
                else:
                    # If not, add the computed time-dependent element flux to the existing key-value pair.
                    flux_dict[pathway] = flux_dict[pathway] + flux
    # Print out ending sequence.
    print("===== Ending reaction pathways analysis of the input plasma mechanism =====")
    # Return the element flux dictionary.
    return flux_dict

"""
A function to combine pathways corresponding to the combustion and plasma mechanisms into a single data structure.
"""
def coupled_pathways(combustion_flux_dict: Dict[str, np.array], plasma_flux_dict: Dict[str, np.array]) -> Dict[str, np.array]:
    """
    Parameters:
    -----------
    combustion_flux_dict : Dict[str, np.array]
        A dictionary representing all pathways corresponding to the combustion mechanism, where the key is the pathway string 
        and the value is the computed time-dependent element flux for each reaction pathway.
    plasma_flux_dict : Dict[str, np.array]
        A dictionary representing all pathways corresponding to the plasma mechanism, where the key is the pathway string 
        and the value is the computed time-dependent element flux for each reaction pathway.

    Returns:
    --------
    combined_dict : Dict[str, np.array]
        A dictionary combining the element flux for both the combustion and plasma mechanisms, where the key is the pathway string
        and the value is the sum of the fluxes for common pathways or the individual flux for unique pathways.
    """
    # Initialize an empty dictionary to store the pathways for both the combustion and plasma mechanisms.
    combined_dict = {}

    # Combine the element flux for both the combustion and plasma mechanisms for common pathways.
    for key in set(combustion_flux_dict.keys()) | set(plasma_flux_dict.keys()):
        if key in combustion_flux_dict and key in plasma_flux_dict:
            combined_dict[key] = combustion_flux_dict[key] + plasma_flux_dict[key]
        elif key in combustion_flux_dict:
            combined_dict[key] = combustion_flux_dict[key]
        elif key in plasma_flux_dict:
            combined_dict[key] = plasma_flux_dict[key]
    # Return the combined pathways dictionary.
    return combined_dict

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