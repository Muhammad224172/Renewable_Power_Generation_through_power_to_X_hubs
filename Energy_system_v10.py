# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 06:01:28 2024

@author: muham
"""
#import libraries and functions 
from gurobipy import GRB,Model,quicksum #*
import pandas as pd
def Output(m):  
    # Print the result
    status_code = {1:'LOADED', 2:'OPTIMAL', 3:'INFEASIBLE', 4:'INF_OR_UNBD', 5:'UNBOUNDED'} #this is how a 'dictionary' 
                                                                                            #is defined in Python
    status = m.status
    
    print('The optimization status is ' + status_code[status])
    if status == 2:    
        # Retrieve variables value
        print('Optimal solution:')
        for v in m.getVars():
            print(str(v.varName) + " = " + str(v.x))    
        print('Optimal objective value: ' + str(m.objVal) + "\n")
        #
################
#import data and parameters 
Electrical_system = pd.read_excel("C:/Users/muham/Desktop/PURE/Electrical system.xlsx")
####power system 
# Electricity demand (in MWh)
elec_demand = Electrical_system["Electricity demand (MWh)"] 
# Hydrogen demand (in MWh)
#h2_demand = Hydrogen_demand['Hydrogen demand (MWh)']
# Availability factors (Wind and solar)
avail_factor_wind = Electrical_system["Wind availability (%)"]
avail_factor_solar = Electrical_system["PV availability (%)"]
#Annualized Investment costs 
ann_inv_cost_wind = 107643
ann_inv_cost_solar = 95730
ann_inv_cost_biomass = 243140
#fixed and operational cost 
fix_om_cost_wind = 30551
fix_om_cost_solar = 18695
fix_om_cost_biomass = 91459
#production cost 
ammonia_cost = 1000
methanol_cost = 1000
#unmet demand cost
unmet = 1000000#at high price to avoid unmet demand 

###Batteries 

#annaual costs 
ann_inv_cost_b1 = 42099
ann_inv_cost_b2 = 51319
#fixed and operational costs 
fix_om_cost_b1 = 22145
fix_om_cost_b2 = 11519
# charging and discharging efficienciies 
charge_b1 = 0.9
charge_b2 = 0.9
discharge_b1 = 0.9
discharge_b2 = 0.9

# Energy to power ratio
e2p_b1 = 4
e2p_b2 = 2


### Hydrogen system

#efficiencies 
electrolizer_1_eff = 0.64
electrolizer_2_eff = 0.80
electrolizer_3_eff = 0.72

#annual investment costs 
electrolizer_1_ann = 148005
electrolizer_2_ann = 579150
electrolizer_3_ann = 92950.02288

### fixed operational and maintainance costs 
electrolizer_1_op = 71875
electrolizer_2_op = 281250
electrolizer_3_op = 45138.9
# charging and discaharging of hydrogen stroage 
h2_Storage_eff_ch = 0.84
h2_Storage_eff_ds = 0.84
#Investement cost 
h2_ann_investment = 5868.72
#fixed and maintaince cost 
h2_fixed_cost = 2850
##eff of fuels 
fuel1_eff = 0.5
fuel2_eff = 0.6 
#Annualized 
fuel1_ann = 133848
fuel2_ann = 339768
#fixed and operational costs
fuel1_op = 65000
fuel2_op = 165000
# Hydrogen sales price (in €/MWh)
h2_sales_price = 75 #considering the cost for each Mwh
#electricity import price 
elec_import_price = 10000
###########################
###########################
#Modelling 
# Create the model
model = Model("Energy_System")
# Set parameters
model.setParam('OutputFlag',True)
# time period 
time_periods = len(Electrical_system["Electricity demand (MWh)"])

#######
#decision variables 

###Capacity of each Technology 
#Note ----> find demand for year 
#max half of peak  demand
# 4500 caapcity  maximum capapcity 
#Using Biomass as an additional technology ---- peak 
#Annualized --> 243140 
#fixed -----> 9149
#production -----> 100 euros per megawatt of electricity since its different from other renewable 
#Biomass acts as a thermal power unit which further reduces the heat requirement for the plan
#	Find the hourly electricity peak demand for the whole year (almost 9,000), and max bound of biomass availability is 4,500
#for biomass
#annualized capital cost (243,140)
#fixed op and maintenance (91,459)

# Capacity of each technology (wind, solar, individual batteries, electrolyzers, hydrogen storage, individual fuel cells)
cap_wind = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="cap_wind")
cap_solar = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="cap_solar")
cap_biomass = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="biomass")

# Individual battery capacities
cap_battery_1 = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="cap_battery_1")
cap_battery_2 = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="cap_battery_2")

# Individual electrolyzer capacities
cap_electrolyzer_1 = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="cap_electrolyzer_1")
cap_electrolyzer_2 = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="cap_electrolyzer_2")
cap_electrolyzer_3 = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="cap_electrolyzer_3")

# Hydrogen storage capacity
cap_h2_storage = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="cap_h2_storage")

# Individual fuel cell capacities
cap_fuelcell_1 = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="cap_fuelcell_1")
cap_fuelcell_2 = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="cap_fuelcell_2")

### Charging power to each electricity storage tecnology 

# Charging power to each electricity storage technology in each hourly time period
chg_pwr_battery_1 = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="chg_pwr_battery_1")
chg_pwr_battery_2 = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="chg_pwr_battery_2")
# discharging power from each electricity storage technology in each hourly time period
dchg_pwr_battery_1 = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="dchg_pwr_battery_1")
dchg_pwr_battery_2 = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="dchg_pwr_battery_2")
#state of energy in each hourly time period 
soe_battery_1 = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="soe_battery_1")
soe_battery_2 = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="soe_battery_2")


###Electricity Imports 
import_elec = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="import_elec")
#Unmet electricity demand 
unmet_elec = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="unmet_elec")


####Electricity production 
prod_wind = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="prod_wind")
prod_solar = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="prod_solar")
prod_biomass = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="prod_biomass")


# Electricity input to each electrolyzer in each hourly time period
elec_input_elec_1 = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="elec_input_elec_1")
elec_input_elec_2 = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="elec_input_elec_2")
elec_input_elec_3 = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="elec_input_elec_3")

# Hydrogen production from each electrolyzer in each hourly time period
h2_prod_elec_1 = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="h2_prod_elec_1")
h2_prod_elec_2 = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="h2_prod_elec_2")
h2_prod_elec_3 = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="h2_prod_elec_3")


# Direct Hydrogen input to all fuel cells in each time period (without being stored)
dir_h2_fuelcell_1 = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="dir_h2_fuelcell_1")
dir_h2_fuelcell_2 = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="dir_h2_fuelcell_2")

# Direct Hydrogen directed to cover hydrogen demand (without being stored)
dir_h2_demand = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="dir_h2_demand")

# discharging from hydrogen storage to fuel cells and demand coverage 
dchg_h2_storage_fuelcell_1 = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="chg_h2_storage_fuelcell_1")
dchg_h2_storage_fuelcell_2 = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="chg_h2_storage_fuelcell_2")
dchg_h2_storage_demand = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="chg_h2_storage_demand")

# Hydrogen input to each fuel cell in each hourly time period
h2_input_fuelcell_1 = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="h2_input_fuelcell_1")
h2_input_fuelcell_2 = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="h2_input_fuelcell_2")


# Electricity production from each fuel cell in each hourly time period
elec_prod_fuelcell_1 = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="elec_prod_fuelcell_1")
elec_prod_fuelcell_2 = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="elec_prod_fuelcell_2")


# Charging hydrogen to each hydrogen storage technology in each hourly time period
chg_h2_storage_hourly = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="chg_h2_storage_hourly")
# Discharging hydrogen to each hydrogen storage technology in each hourly time period
dchg_h2_storage_hourly = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="dchg_h2_storage_hourly")

# State-of-energy of each hydrogen storage technology in each hourly time period
soe_h2_storage = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="soe_h2_storage")
#Total demand 
h2_ammonia = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="h2_ammonia")
h2_methanol = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="h2_methanol")
h2_output = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="h2_output")


ammonia = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="ammonia")
elec_ammonia = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="elec_ammonia")
n2 = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="n2")

methanol = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="methanol")
elec_methanol = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="elec_methanol")
co2 = model.addVars(time_periods, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="co2")

##########################
##########################

# Objective function components

# 1. Cost of Capacity Investment
capacity_cost = (
    cap_wind * ann_inv_cost_wind +
    cap_solar * ann_inv_cost_solar +
    cap_biomass * ann_inv_cost_biomass +
    cap_battery_1 * ann_inv_cost_b1 +
    cap_battery_2 * ann_inv_cost_b2 +
    cap_electrolyzer_1 * electrolizer_1_ann +
    cap_electrolyzer_2 * electrolizer_2_ann +
    cap_electrolyzer_3 * electrolizer_3_ann +
    cap_h2_storage * h2_ann_investment +
    cap_fuelcell_1 * fuel1_ann +
    cap_fuelcell_2 * fuel2_ann
)

# 2. Fixed Operational and Maintenance Costs
fixed_om_cost = (
    cap_wind * fix_om_cost_wind +
    cap_solar * fix_om_cost_solar +
    cap_biomass * fix_om_cost_biomass +
    cap_battery_1 * fix_om_cost_b1 +
    cap_battery_2 * fix_om_cost_b2 +
    cap_electrolyzer_1 * electrolizer_1_op +
    cap_electrolyzer_2 * electrolizer_2_op +
    cap_electrolyzer_3 * electrolizer_3_op +
    cap_h2_storage * h2_fixed_cost +
    cap_fuelcell_1 * fuel1_op +
    cap_fuelcell_2 * fuel2_op
)

# 3. Electricity Imports Cost
import_cost = quicksum(import_elec[t] * elec_import_price for t in range(time_periods))

# 4. Hydrogen Sales Revenue 
hydrogen_revenue = quicksum((h2_output[t]) * h2_sales_price for t in range(time_periods))

# 5. Unmet demand 
unmet_demand = quicksum(unmet_elec[t] * unmet for t in range(time_periods))

##6. Biomass production 

biomass_cost = quicksum(prod_biomass[t] * 100 for t in range(time_periods))

#synthesis 
ammonia_Revenue = quicksum(ammonia[t] * 800 for t in range(time_periods))
ammonia_production_cost = quicksum( n2[t] * 1 + ammonia[t] * 50 for t in range(time_periods))

methanol_Revenue = quicksum(methanol[t] * 400 for t in range(time_periods))
methanol_production_cost = quicksum( co2[t] * 50 + methanol[t] * 50 for t in range(time_periods))

##total cost 
total_cost =   (capacity_cost + fixed_om_cost + import_cost + unmet_demand + biomass_cost + ammonia_production_cost + methanol_production_cost)
- hydrogen_revenue - ammonia_Revenue - methanol_Revenue

# Objective: Minimize Total Cost
model.setObjective(total_cost, GRB.MINIMIZE)
###############################################
###############################################
# 1. Electricity demand balance
#############
for t in range(time_periods):
    model.addConstr(
             import_elec[t] + prod_solar[t] + prod_wind[t] + prod_biomass[t] + dchg_pwr_battery_1[t] + dchg_pwr_battery_2[t] +
             unmet_elec[t] + elec_prod_fuelcell_1[t] + elec_prod_fuelcell_2[t] == 
             elec_demand[t] + chg_pwr_battery_1[t] + chg_pwr_battery_2[t] + elec_input_elec_1[t] + elec_input_elec_2[t]
             + elec_input_elec_3[t] + (0.34*elec_ammonia[t]) + (0.169 * elec_methanol[t]) 
    )
    model.addConstr(
        prod_solar[t] + prod_wind[t] + prod_biomass[t] + import_elec[t] + dchg_pwr_battery_1[t] + dchg_pwr_battery_2[t] == 
        elec_input_elec_1[t] + elec_input_elec_2[t]
        + elec_input_elec_3[t] 
        + (0.34*elec_ammonia[t]) + (0.169 * elec_methanol[t]) 
        )


######################
############
# 2. Eleticity Production from each RES tecnology 
for t in range(time_periods):
    model.addConstr( 
    prod_solar[t] == cap_solar * avail_factor_solar[t])
    model.addConstr(
     prod_wind[t] == cap_wind * avail_factor_wind[t]
        )
    model.addConstr(
        prod_biomass[t] == cap_biomass * 4500  
        )
############
############
#3. Electricity Storage modelling 

#Initial State 
soe_battery_1[0] = 0 
soe_battery_1[0] = 0
for t in range(1 , time_periods):
    #State of energy of batteries 
    model.addConstr(
        soe_battery_1[t]  == soe_battery_1[t-1] + (chg_pwr_battery_1[t]*charge_b1) - (dchg_pwr_battery_1[t]/discharge_b1)
        )
    model.addConstr(
      soe_battery_2[t] ==  soe_battery_2[t-1] + (chg_pwr_battery_2[t]*charge_b2) - (dchg_pwr_battery_2[t]/discharge_b2)
        )
    #Charging power capapcity 
    model.addConstr(
        chg_pwr_battery_1[t] <= cap_battery_1/e2p_b1 
        )
    model.addConstr(
       chg_pwr_battery_2[t] <= cap_battery_2/e2p_b2
        )
    #discharging capapcity 
    model.addConstr(
        dchg_pwr_battery_1[t] <= cap_battery_1/e2p_b1 
        )
    model.addConstr(
        dchg_pwr_battery_2[t] <= cap_battery_2/e2p_b2
        )
    
    #Constrain on capacity of storage systems 
    model.addConstr(
        soe_battery_1[t]  <= cap_battery_1 
        )
    model.addConstr(
        soe_battery_2[t] <= cap_battery_2
        )
#########################################################
#########################################################
# 4. Electrolyzer modelling 
   
for t in range(time_periods):
    #Hydrogen production from each electrolizer 
    model.addConstr(
           h2_prod_elec_1[t]  == (elec_input_elec_1[t] * electrolizer_1_eff)) 
    model.addConstr(
           h2_prod_elec_2[t]  == (elec_input_elec_2[t] * electrolizer_2_eff))  
    model.addConstr(
           h2_prod_elec_3[t] == (elec_input_elec_3[t] * electrolizer_3_eff))
    
    #Electrolizer input capacity contrain
    model.addConstr(
        elec_input_elec_1[t]  <= cap_electrolyzer_1 
        )
    model.addConstr(
        elec_input_elec_2[t] <= cap_electrolyzer_2
        )
    model.addConstr(
        elec_input_elec_3[t] <= cap_electrolyzer_3
        )

#########################################################
#########################################################
#5.  hydrogen supply balanace 
for t in range(time_periods):
    model.addConstr(
        h2_prod_elec_1[t] + h2_prod_elec_2[t] + h2_prod_elec_3[t] ==
        dir_h2_fuelcell_1[t] + dir_h2_fuelcell_2[t] + dir_h2_demand[t] 
        + chg_h2_storage_hourly[t]        
        )
    #Additional
    
    model.addConstr(
        chg_h2_storage_hourly[t] <= cap_h2_storage - soe_h2_storage[t]
        )
    model.addConstr(
        dchg_h2_storage_demand[t] + dchg_h2_storage_fuelcell_1[t] + dchg_h2_storage_fuelcell_2[t] <= soe_h2_storage[t]
        )

##########################
##########################
#6.  Hydrogen storage modelling 
for t in range(1,time_periods):
    soe_h2_storage[0] = 0 
    model.addConstr(
        soe_h2_storage[t] == soe_h2_storage[t-1] + 
        (chg_h2_storage_hourly[t]*h2_Storage_eff_ch) -
        (dchg_h2_storage_demand[t]/h2_Storage_eff_ds) -
        ((dchg_h2_storage_fuelcell_1[t] + dchg_h2_storage_fuelcell_2[t])/h2_Storage_eff_ds)
        )
    #cpapcity storage
    model.addConstr(
        soe_h2_storage[t] <= cap_h2_storage
        )  
##########################
##########################
# 7. Fuel cell modelling 

for t in range(time_periods):

    #electriicty production constrains 
    model.addConstr(
        elec_prod_fuelcell_1[t]  == ((dir_h2_fuelcell_1[t] +  dchg_h2_storage_fuelcell_1[t]) * fuel1_eff) 
        )
    model.addConstr(
        elec_prod_fuelcell_2[t] ==  ((dir_h2_fuelcell_2[t] +  dchg_h2_storage_fuelcell_2[t])*fuel2_eff)
        )
    #Electricity production Capacity Constrains 
    model.addConstr(
        elec_prod_fuelcell_1[t] <= cap_fuelcell_1
        ) 
    model.addConstr(
        elec_prod_fuelcell_2[t] <= cap_fuelcell_2
        ) 
######################


####################################
####################################
# synthesis 
for t in range(time_periods):
    model.addConstr(
        h2_output[t] + (0.192 * h2_methanol[t]) + (0.18*h2_ammonia[t]) == dir_h2_demand[t]  + dchg_h2_storage_demand[t]
        )
    #production of methanol and ammonia 
    model.addConstr(
        ammonia[t] == (0.34*elec_ammonia[t]) + (0.18*h2_ammonia[t]) + (0.84*n2[t])
        )
    model.addConstr(
        methanol[t] == (0.169 * elec_methanol[t]) + (0.192 * h2_methanol[t]) + (1.37*co2[t])
        )


    

        
######################
######################
#solve and output solutions 
model.optimize()
Output(model)    
# print the LP file
model.write('Energy_System.lp')
# print the sol file
model.write('Energy_System.sol')

#################################
#################################
#################################

data = {
    'Time Period': list(range(time_periods))
}

# Add time-varying variables to the dictionary

data['import_elec'] = [import_elec[t].X for t in range(time_periods)]

data['h2_prod_elec'] = [
    h2_prod_elec_1[t].X + h2_prod_elec_2[t].X + h2_prod_elec_3[t].X for t in range(time_periods)
]
data['total_energy_prod'] = [
    prod_wind[t].X + prod_solar[t].X  for t in range(time_periods)
]
data['prod_solar'] = [
     prod_solar[t].X  for t in range(time_periods)
]
data['prod_wind'] = [
    prod_wind[t].X  for t in range(time_periods)]
    
data['prod_biomass'] = [
    prod_biomass[t].X  for t in range(time_periods)
]
data['Ammonia'] = [
    ammonia[t].X  for t in range(time_periods)
]
data['methanol'] = [
    methanol[t].X  for t in range(time_periods)
]
data['h2_output'] = [
    h2_output[t].X for t in range(time_periods)
]


# Add constant variables (capacities) to the dictionary
data['cap_wind'] = [cap_wind.X] * time_periods
data['cap_solar'] = [cap_solar.X] * time_periods
data['cap_biomass'] = [cap_biomass.X] * time_periods
data['cap_battery_1'] = [cap_battery_1.X] * time_periods
data['cap_battery_2'] = [cap_battery_2.X] * time_periods
data['cap_electrolyzer_1'] = [cap_electrolyzer_1.X] * time_periods
data['cap_electrolyzer_2'] = [cap_electrolyzer_2.X] * time_periods
data['cap_electrolyzer_3'] = [cap_electrolyzer_3.X] * time_periods
data['cap_h2_storage'] = [cap_h2_storage.X] * time_periods
data['cap_fuelcell_1'] = [cap_fuelcell_1.X] * time_periods
data['cap_fuelcell_2'] = [cap_fuelcell_2.X] * time_periods

# Convert dictionary to DataFrame
df = pd.DataFrame(data)

# Save DataFrame to CSV
df.to_csv('C:/Users/muham/Desktop/PURE/Energy_system_v12.csv', index=False)


################################################
################################################ 
#Analytical values 
print("The capacity of the wind power is " , cap_wind)
print("The capacity of the solar power is " , cap_solar)
print("The capacity of the biomass is " , cap_biomass)
print("The capacity of the battery 1 is " , cap_battery_1)
print("The capacity of the battery 2 is " , cap_battery_2)
print("The capacity of the electrolizer 1 is " , cap_electrolyzer_1)
print("The capacity of the electrolizer 2 is " , cap_electrolyzer_2)
print("The capacity of the electrolizer 3 is " , cap_electrolyzer_3)
print("The capacity of the fuel cell 1 is " , cap_fuelcell_1)
print("The capacity of the fuel cell 2 is " , cap_fuelcell_2)








































