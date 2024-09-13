# Created by: Nicholas Laws
# Date: 2023

## Imports ##
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import os
import glob
from typing import *

"""
A function to take the time-dependent derivative of a 0D temperature profile.
"""
def get_deriv(temp: np.array, time: np.array, PRF: Union[float, int]) -> List[np.array]:
    """
    Parameters:
    -----------
    temp : np.array
        A numpy array representing the temperature as a function of time for the 0D profile.
    time : np.array
        A numpy array representing the time associated with the 0D temperature profile.
    PRF : Union[float, int]
        A float or integer value representing the pulse repetition frequency of the 0D profile.

    Returns:
    --------
    [time_filtered, temp_filtered, deriv_filtered] : List[np.array]
        - time_filtered: Filtered time array corresponding to the processed temperature data.
        - temp_filtered: Filtered temperature array after processing.
        - deriv_filtered: Time-dependent derivative (dT/dt) of the filtered temperature array.
    """
    # Take the approximate derivative, delta(T)/delta(t).
    deriv = np.diff(temp)/np.diff(time)
    # Define relevant values.
    FWHM = 20e-9
    sigma = FWHM/(2*np.sqrt(2*np.log(2)))
    duration_pulse = 10*sigma
    N = 0
    period = 1/PRF
    step = 1e-9
    # Create empty lists to represent the filtered time, filtered derivative, and filtered temperature arrays.
    time_filtered = []
    deriv_filtered = []
    temp_filtered = []
    for j in range(len(deriv)):
        # Temporary values of the unfiltered time, derivative, and temperature arrays.
        temp_time = time[j]
        temp_deriv = deriv[j]
        temp_temp = temp[j]
        # Check if the current time exceeds the number of pulses plus one times the period. If so, one pulse has occurred.
        if temp_time >= (N + 1)*period:
            N += 1
        # Start looking at the unfiltered arrays after the tenth value.
        if j >= 10:
            if temp[j] - temp[j-10] > 500 and temp_time <= period:
                time_filtered.append(temp_time)
                deriv_filtered.append(temp_deriv)
                temp_filtered.append(temp_temp)
            # Check if the time falls between [N*T - pulse_duration, N*T + pulse_duration]. If so, pass.
            elif temp_time >= N*period - duration_pulse*400 and temp_time <= N*period + duration_pulse*400 or np.isnan(temp_deriv):
                pass
            # If not, append the value.
            else:
                time_filtered.append(temp_time)
                deriv_filtered.append(temp_deriv)
                temp_filtered.append(temp_temp)
    return [time_filtered, temp_filtered, deriv_filtered]

"""
A function to compute the ignition delay time (IDT).
"""
def get_idt(temp: np.array, time: np.array, PRF: Union[float, int]) -> float:
    """
    Parameters:
    -----------
    temp : np.array
        A numpy array representing the temperature as a function of time for the 0D profile.
    time : np.array
        A numpy array representing the time associated with the 0D temperature profile.
    PRF : Union[float, int]
        A float or integer value representing the pulse repetition frequency of the 0D profile.

    Returns:
    --------
    idt : float
        The ignition delay time (IDT), or NaN if ignition does not occur.
    """
    # Generate filtered time, temperature, and derivative arrays.
    filt_time, filt_temp, filt_deriv = get_deriv(temp, time, PRF)
    # Find maximum derivative value; base IDT index on this maximum. 
    idt_idx = np.argmax(filt_deriv)
    idt_idx_min = np.argmin(filt_deriv)
    idt_temp = filt_temp[idt_idx]
    max_deriv = filt_deriv[idt_idx]
    min_deriv = filt_deriv[idt_idx_min]
    # Check if the maximum derivative falls below 5e4. If so, ignition does not happen.
    if max_deriv < 5e4:
        idt = math.nan
        idt_temp = math.nan
        #print("Ignition does not occur.")
    # Fail-safe check.
    elif min_deriv < -5e4:
        idt_temp = filt_temp[idt_idx_min]
        idt = filt_time[idt_idx_min]
    # If not, ignition does happen.
    else:
        idt = filt_time[idt_idx]
    return idt

"""
A function to compute the maximum temperature achieved just after ignition.
"""
def get_max_ignition_temp(filtered_time: np.array, filtered_temp: np.array, 
                          filtered_deriv: np.array, idt: float) -> List[float]:
    """
    Parameters:
    -----------
    filtered_temp : np.array
        A numpy array representing the filtered temperature as a function of time for the 0D profile.
    filtered_time : np.array
        A numpy array representing the filtered time associated with the 0D temperature profile.
    filtered_deriv : np.array
        A numpy array representing the filtered time-dependent temperature derivative associated with the 0D profile.
    idt : float
        A float value representing the ignition delay time (IDT).

    Returns:
    --------
    [ignition_temp, ignition_temp_time] : List[float]
        - ignition_temp: The maximum temperature achieved after ignition.
        - ignition_temp_time: The time at which this maximum temperature is achieved.
    """
    # Check ignition happens. If so, find maximum temperature based on IDT.
    if np.isnan(idt) == False:
        # Generate time, temperature, and derivative arrays based on IDT.
        index_after_idt = filtered_time.index(idt)
        time_after_idt = filtered_time[index_after_idt:]
        temp_after_idt = filtered_temp[index_after_idt:]
        deriv_after_idt = filtered_deriv[index_after_idt:]
        # Find the magnitude of the derivative's maximum.
        mag = int(math.log10(abs(np.max(filtered_deriv))))
        # Find when the maximum derivative is less than the maximum derivative's order of magnitude.
        ignition_temp_idx = 0
        for k in range(len(deriv_after_idt)):
            temp_deriv = deriv_after_idt[k]
            if temp_deriv < mag:
                ignition_temp_idx = k
                break
        ignition_temp = temp_after_idt[ignition_temp_idx]
        ignition_temp_time = time_after_idt[ignition_temp_idx]
    else:
        time_after_idt = filtered_time
        temp_after_idt = filtered_temp
        deriv_after_idt = filtered_deriv
        # Find the magnitude of the derivative's maximum.
        mag = int(math.log10(abs(np.max(filtered_deriv))))
        # Find when the maximum derivative is less than the maximum derivative's order of magnitude.
        ignition_temp_idx = 0
        for k in range(len(deriv_after_idt)):
            temp_deriv = deriv_after_idt[k]
            if temp_deriv < mag:
                ignition_temp_idx = k
                break
        ignition_temp = temp_after_idt[ignition_temp_idx]
        ignition_temp_time = time_after_idt[ignition_temp_idx]
    return [ignition_temp, ignition_temp_time]