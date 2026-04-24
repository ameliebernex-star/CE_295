import pandas as pd

#-----------------------------------------------
"""REAL-TIME PRICE"""

def energy_cost(p_k, rt_price, alpha = 1):
  # p_k: power injection (MW)
  # rt_price: real time prices ($/MWh)
  # alpha: time step (hour)
  total = 0
  for k in range(len(rt_price)):
    total += rt_price[k] * p_k * alpha
  return total
#-----------------------------------------------
"""OPERATING COSTS"""

def operating_cost(p_k, opex_params, ):
    raise NotImplementedError
#-----------------------------------------------
"""DEGRADATION COSTS"""

def degradation_cost(p_k, battery_params):
    raise NotImplementedError
#-----------------------------------------------
"""OVERALL OBJECTIVE FXN"""

def total_cost(dispatch, prices, battery_params, opex_params):
    return energy_cost(dispatch, prices) + \
           degradation_cost(dispatch, battery_params) + \
           operating_cost(dispatch, opex_params)
