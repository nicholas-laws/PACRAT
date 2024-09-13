# Created by: Nicholas Laws
# Date: 2024

## Imports ##
from plasma_pathways import *
from plasma_chemistry_objs import *
import pandas as pd
import time
import os

"""
A function to adjust the combustion data to be tracked on the ZDPlasKin timescale.
"""
def process_plasma_data(plasma_data: pd.DataFrame, comb_data: pd.DataFrame, file_name: str) -> pd.DataFrame:
    """
    Parameters:
    -----------
    plasma_data : pd.DataFrame
        A pandas DataFrame storing all time-dependent number densities and reaction rate coefficients.
    comb_data : pd.DataFrame
        A Pandas DataFrame storing all time-dependent mass fraction data for all combustion reactions in the numerical model.
    file_name : str
        A user-selected string that represents the output filename for the adjusted combustion data.

    Returns:
    --------
    comb_data_adj : pd.DataFrame
        A Pandas DataFrame containing the adjusted combustion data aligned with the ZDPlasKin timescale.
    """
    # Get the directory of the currently running script.
    dir = os.path.dirname(os.path.abspath(__file__))
    # Initialize a variable that represents the boolean for Cantera tracked time-intervals.
    ct_bool = plasma_data["PDtt_ct"].to_numpy()
    # Initialize temporary placeholder variables for the adjusted slices of data and the combustion index.
    comb_data_adj_slices = []
    comb_idx = 0
    # Iterate through the length of the Cantera boolean.
    for idx in range(len(ct_bool)):
        # Check if ct_bool at index idx is equal to 0.
        if ct_bool[idx] == 0:
            # If so, take the row of data at comb_idx.
            comb_data_adj_slices.append(comb_data.iloc[comb_idx:comb_idx+1])
        else:
            # If not, iterate comb_idx, and take the data from the next row.
            comb_idx += 1
            comb_data_adj_slices.append(comb_data.iloc[comb_idx:comb_idx+1])
    # Concatenate the data from each slice into a single DataFrame.
    comb_data_adj = pd.concat(comb_data_adj_slices, ignore_index=True)
    # Add the plasma time column to the DataFrame.
    comb_data_adj['t'] = plasma_data['PDtt'].values
    # Write the adjusted combustion data to a .CSV.
    comb_data_adj.to_csv(os.path.join(dir, f'data/output/{file_name}.csv'), index=False)
    # Return the adjusted combustion data.
    return comb_data_adj

# Get the directory of the currently running script.
dir = os.path.dirname(os.path.abspath(__file__))
# Define the filepath of the plasma kinetic mechanism.
mech_path = os.path.join(dir, "data/plasma/kinet_CH4_PAC_V5_nrg.inp") 
# Read the user-selected plasma data. 
plasma_data = pd.read_csv(os.path.join(dir, "data/plasma/result_a00_b05_c00104_PD.csv"))
# Read the user-selected combustion data. 
comb_data = pd.read_csv(os.path.join(dir, "data/combustion/result_a00_b05_c00104.csv") )
# Initialize an output filename.
file_name = "result_a00_b05_c00104_adj"
# Generate the time-adjusted combustion data.
comb_data_adj = process_plasma_data(plasma_data, comb_data, file_name)
