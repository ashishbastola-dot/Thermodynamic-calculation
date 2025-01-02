# =============================================================================
# Import required packages
# =============================================================================

from reaktoro import *
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Database selection
# Ensure that the database file is in the working directory
# =============================================================================

db = PhreeqcDatabase().load("CEMDATA18-31-03-2022-phaseVol.dat")

# =============================================================================
# System definition and activity model selection
# =============================================================================

aqueous = AqueousPhase(speciate("Al Si Ca K Mg Na S O H"))

system = ChemicalSystem(db, aqueous)

aqueous.set(ActivityModelPitzer())

# =============================================================================
# Define specification to setup pH constraint
# =============================================================================
specs = EquilibriumSpecs(system)
specs.temperature()
specs.pressure()
specs.pH()

# =============================================================================
# Import results saved in a file
# =============================================================================

stab_a = pd.read_excel("result.xlsx",sheet_name="stabA")    # ICP result
stab_b = pd.read_excel("result.xlsx",sheet_name="stabB")    # ICP result
type_3 = pd.read_excel("result.xlsx",sheet_name="type3")    # ICP result
csa = pd.read_excel("result.xlsx",sheet_name="csa")         # ICP result
strength = pd.read_excel("result.xlsx",sheet_name="strength")   # UCS result

# =============================================================================
# Create  chemical state
# =============================================================================

df = stab_b                         # select the ICP result to be analyzed
strength=strength["stab_b"]         # select corresponding strength data   

# ----------------------------------------------------------------------------
# Cheimcal state for 0 min

t = 0

state = ChemicalState(system)
state.set("H2O",  1, "kg")
state.set("Al+3",  df['Al'][t], "mol")
state.set("Ca+2",  df['Ca'][t], "mol")
state.set("K+",  df['K'][t], "mol")
state.set("Na+",  df['Na'][t], "mol")
state.set("SO4-2",  df['S'][t], "mol")
state.set("SiO2",  df['Si'][t], "mol")

state_0=state

conditions_0 = EquilibriumConditions(specs)
conditions_0.temperature(25.0, "celsius")
conditions_0.pressure(1.0, "atm")
conditions_0.pH(df['pH'][t])

# ----------------------------------------------------------------------------

# Chemical state for 15 min

t = 1
state = ChemicalState(system)
state.set("H2O",  1, "kg")
state.set("Al+3",  df['Al'][t], "mol")
state.set("Ca+2",  df['Ca'][t], "mol")
state.set("K+",  df['K'][t], "mol")
state.set("Na+",  df['Na'][t], "mol")
state.set("SO4-2",  df['S'][t], "mol")
state.set("SiO2",  df['Si'][t], "mol")

state_1=state

conditions_1 = EquilibriumConditions(specs)
conditions_1.temperature(25.0, "celsius")
conditions_1.pressure(1.0, "atm")
conditions_1.pH(df['pH'][t])

# ----------------------------------------------------------------------------

# Chemical state for 30 min

t = 2
state = ChemicalState(system)
state.set("H2O",  1, "kg")
state.set("Al+3",  df['Al'][t], "mol")
state.set("Ca+2",  df['Ca'][t], "mol")
state.set("K+",  df['K'][t], "mol")
state.set("Na+",  df['Na'][t], "mol")
state.set("SO4-2",  df['S'][t], "mol")
state.set("SiO2",  df['Si'][t], "mol")

state_2=state

conditions_2 = EquilibriumConditions(specs)
conditions_2.temperature(25.0, "celsius")
conditions_2.pressure(1.0, "atm")
conditions_2.pH(df['pH'][t])

# ----------------------------------------------------------------------------

# Chemical state for 1 hr

t = 3
state = ChemicalState(system)
state.set("H2O",  1, "kg")
state.set("Al+3",  df['Al'][t], "mol")
state.set("Ca+2",  df['Ca'][t], "mol")
state.set("K+",  df['K'][t], "mol")
state.set("Na+",  df['Na'][t], "mol")
state.set("SO4-2",  df['S'][t], "mol")
state.set("SiO2",  df['Si'][t], "mol")

state_3=state

conditions_3 = EquilibriumConditions(specs)
conditions_3.temperature(25.0, "celsius")
conditions_3.pressure(1.0, "atm")
conditions_3.pH(df['pH'][t])

# ----------------------------------------------------------------------------

# Chemical state for 1 d

t = 4
state = ChemicalState(system)
state.set("H2O",  1, "kg")
state.set("Al+3",  df['Al'][t], "mol")
state.set("Ca+2",  df['Ca'][t], "mol")
state.set("K+",  df['K'][t], "mol")
state.set("Na+",  df['Na'][t], "mol")
state.set("SO4-2",  df['S'][t], "mol")
state.set("SiO2",  df['Si'][t], "mol")

state_4=state

conditions_4 = EquilibriumConditions(specs)
conditions_4.temperature(25.0, "celsius")
conditions_4.pressure(1.0, "atm")
conditions_4.pH(df['pH'][t])

# ----------------------------------------------------------------------------

# Chemical state for 7 d

t = 5
state = ChemicalState(system)
state.set("H2O",  1, "kg")
state.set("Al+3",  df['Al'][t], "mol")
state.set("Ca+2",  df['Ca'][t], "mol")
state.set("K+",  df['K'][t], "mol")
state.set("Na+",  df['Na'][t], "mol")
state.set("SO4-2",  df['S'][t], "mol")
state.set("SiO2",  df['Si'][t], "mol")

state_5=state

conditions_5 = EquilibriumConditions(specs)
conditions_5.temperature(25.0, "celsius")
conditions_5.pressure(1.0, "atm")
conditions_5.pH(df['pH'][t])

# =============================================================================
# Solve for the equilibirum system
# =============================================================================

# Define solver
solver = EquilibriumSolver(specs)

# Solve chemical state for 0 min
state = state_0; solver.solve(state, conditions_0)
aprops_0 = AqueousProps(state)

# Solve chemical state for 15 min
state = state_1; solver.solve(state, conditions_1)
aprops_1 = AqueousProps(state)

# Solve chemical state for 30 min
state = state_2; solver.solve(state, conditions_2)
aprops_2 = AqueousProps(state)

# Solve chemical state for 1hr
state = state_3; solver.solve(state, conditions_3)
aprops_3 = AqueousProps(state)

# Solve chemical state for 1d
state = state_4; solver.solve(state, conditions_4)
aprops_4 = AqueousProps(state)

# Solve chemical state for 7d
state = state_5; solver.solve(state, conditions_5)
aprops_5 = AqueousProps(state)

# =============================================================================
# Satuation index of the phases of interest
# =============================================================================

csh=np.array([])
ettrin=np.array([])
gyp=np.array([])
port=np.array([])
gib=np.array([])
strat=np.array([])
mono=np.array([])

csh_eff=np.array([])
ettrin_eff=np.array([])
gyp_eff=np.array([])
port_eff=np.array([])
gib_eff=np.array([])
strat_eff=np.array([])
mono_eff=np.array([])


# At 0 min
aprops = aprops_0

csh=np.append(csh,max(aprops.saturationIndex("CSH3T-T2C"),aprops.saturationIndex("CSH3T-T5C"),
        aprops.saturationIndex("CSH3T-TobH"),aprops.saturationIndex("CSHQ-JenD"),
        aprops.saturationIndex("CSHQ-JenH"),aprops.saturationIndex("CSHQ-TobD"),
        aprops.saturationIndex("CSHQ-TobH"))[0])
ettrin = np.append(ettrin,max(aprops.saturationIndex("ettringite"),aprops.saturationIndex("ettringite13"),
              aprops.saturationIndex("Ettringite13_des"),aprops.saturationIndex("ettringite30"),
              aprops.saturationIndex("ettringite9"),aprops.saturationIndex("Ettringite9_des"))[0])
gyp=np.append(gyp,aprops.saturationIndex("Gp")[0])
gib=np.append(gib,aprops.saturationIndex("Gbs")[0])
port=np.append(port,aprops.saturationIndex("Portlandite")[0])
strat=np.append(strat,max(aprops.saturationIndex("straetlingite"),aprops.saturationIndex("straetlingite5_5"),
          aprops.saturationIndex("straetlingite7"))[0])
mono=np.append(mono,max(aprops.saturationIndex("monosulphate10_5"),aprops.saturationIndex("monosulphate12"),
          aprops.saturationIndex("monosulphate1205"),aprops.saturationIndex("monosulphate14"),
          aprops.saturationIndex("monosulphate16"),aprops.saturationIndex("monosulphate9"))[0])

csh_eff=np.append(csh_eff,max(aprops.saturationIndex("CSH3T-T2C")/5.5,aprops.saturationIndex("CSH3T-T5C")/5,
        aprops.saturationIndex("CSH3T-TobH")/4.5,aprops.saturationIndex("CSHQ-JenD")/5.167,
        aprops.saturationIndex("CSHQ-JenH")/4.999,aprops.saturationIndex("CSHQ-TobD")/3.166825,
        aprops.saturationIndex("CSHQ-TobH")/3.0001)[0])
gyp_eff=np.append(gyp_eff,aprops.saturationIndex("Gp")[0]/2)
gib_eff=np.append(gib_eff,aprops.saturationIndex("Gbs")[0]/2)
port_eff=np.append(port_eff,aprops.saturationIndex("Portlandite")[0]/3)
strat_eff=np.append(strat_eff,max(aprops.saturationIndex("straetlingite"),aprops.saturationIndex("straetlingite5_5"),
          aprops.saturationIndex("straetlingite7"))[0]/7)
ettrin_eff = np.append(ettrin_eff,max(aprops.saturationIndex("ettringite"),aprops.saturationIndex("ettringite13"),
              aprops.saturationIndex("Ettringite13_des"),aprops.saturationIndex("ettringite30"),
              aprops.saturationIndex("ettringite9"),aprops.saturationIndex("Ettringite9_des"))[0]/15)
mono_eff=np.append(mono_eff,max(aprops.saturationIndex("monosulphate10_5"),aprops.saturationIndex("monosulphate12"),
          aprops.saturationIndex("monosulphate1205"),aprops.saturationIndex("monosulphate14"),
          aprops.saturationIndex("monosulphate16"),aprops.saturationIndex("monosulphate9"))[0]/11)

# At 15 min

aprops = aprops_1
csh=np.append(csh,max(aprops.saturationIndex("CSH3T-T2C"),aprops.saturationIndex("CSH3T-T5C"),
        aprops.saturationIndex("CSH3T-TobH"),aprops.saturationIndex("CSHQ-JenD"),
        aprops.saturationIndex("CSHQ-JenH"),aprops.saturationIndex("CSHQ-TobD"),
        aprops.saturationIndex("CSHQ-TobH"))[0])
ettrin = np.append(ettrin,max(aprops.saturationIndex("ettringite"),aprops.saturationIndex("ettringite13"),
              aprops.saturationIndex("Ettringite13_des"),aprops.saturationIndex("ettringite30"),
              aprops.saturationIndex("ettringite9"),aprops.saturationIndex("Ettringite9_des"))[0])
gyp=np.append(gyp,aprops.saturationIndex("Gp")[0])
gib=np.append(gib,aprops.saturationIndex("Gbs")[0])
port=np.append(port,aprops.saturationIndex("Portlandite")[0])
strat=np.append(strat,max(aprops.saturationIndex("straetlingite"),aprops.saturationIndex("straetlingite5_5"),
          aprops.saturationIndex("straetlingite7"))[0])
mono=np.append(mono,max(aprops.saturationIndex("monosulphate10_5"),aprops.saturationIndex("monosulphate12"),
          aprops.saturationIndex("monosulphate1205"),aprops.saturationIndex("monosulphate14"),
          aprops.saturationIndex("monosulphate16"),aprops.saturationIndex("monosulphate9"))[0])

csh_eff=np.append(csh_eff,max(aprops.saturationIndex("CSH3T-T2C")/5.5,aprops.saturationIndex("CSH3T-T5C")/5,
        aprops.saturationIndex("CSH3T-TobH")/4.5,aprops.saturationIndex("CSHQ-JenD")/5.167,
        aprops.saturationIndex("CSHQ-JenH")/4.999,aprops.saturationIndex("CSHQ-TobD")/3.166825,
        aprops.saturationIndex("CSHQ-TobH")/3.0001)[0])
gyp_eff=np.append(gyp_eff,aprops.saturationIndex("Gp")[0]/2)
gib_eff=np.append(gib_eff,aprops.saturationIndex("Gbs")[0]/2)
port_eff=np.append(port_eff,aprops.saturationIndex("Portlandite")[0]/3)
strat_eff=np.append(strat_eff,max(aprops.saturationIndex("straetlingite"),aprops.saturationIndex("straetlingite5_5"),
          aprops.saturationIndex("straetlingite7"))[0]/7)
ettrin_eff = np.append(ettrin_eff,max(aprops.saturationIndex("ettringite"),aprops.saturationIndex("ettringite13"),
              aprops.saturationIndex("Ettringite13_des"),aprops.saturationIndex("ettringite30"),
              aprops.saturationIndex("ettringite9"),aprops.saturationIndex("Ettringite9_des"))[0]/15)
mono_eff=np.append(mono_eff,max(aprops.saturationIndex("monosulphate10_5"),aprops.saturationIndex("monosulphate12"),
          aprops.saturationIndex("monosulphate1205"),aprops.saturationIndex("monosulphate14"),
          aprops.saturationIndex("monosulphate16"),aprops.saturationIndex("monosulphate9"))[0]/11)


# At 30 min

aprops = aprops_2
csh=np.append(csh,max(aprops.saturationIndex("CSH3T-T2C"),aprops.saturationIndex("CSH3T-T5C"),
        aprops.saturationIndex("CSH3T-TobH"),aprops.saturationIndex("CSHQ-JenD"),
        aprops.saturationIndex("CSHQ-JenH"),aprops.saturationIndex("CSHQ-TobD"),
        aprops.saturationIndex("CSHQ-TobH"))[0])
ettrin = np.append(ettrin,max(aprops.saturationIndex("ettringite"),aprops.saturationIndex("ettringite13"),
              aprops.saturationIndex("Ettringite13_des"),aprops.saturationIndex("ettringite30"),
              aprops.saturationIndex("ettringite9"),aprops.saturationIndex("Ettringite9_des"))[0])
gyp=np.append(gyp,aprops.saturationIndex("Gp")[0])
gib=np.append(gib,aprops.saturationIndex("Gbs")[0])
port=np.append(port,aprops.saturationIndex("Portlandite")[0])
strat=np.append(strat,max(aprops.saturationIndex("straetlingite"),aprops.saturationIndex("straetlingite5_5"),
          aprops.saturationIndex("straetlingite7"))[0])
mono=np.append(mono,max(aprops.saturationIndex("monosulphate10_5"),aprops.saturationIndex("monosulphate12"),
          aprops.saturationIndex("monosulphate1205"),aprops.saturationIndex("monosulphate14"),
          aprops.saturationIndex("monosulphate16"),aprops.saturationIndex("monosulphate9"))[0])

csh_eff=np.append(csh_eff,max(aprops.saturationIndex("CSH3T-T2C")/5.5,aprops.saturationIndex("CSH3T-T5C")/5,
        aprops.saturationIndex("CSH3T-TobH")/4.5,aprops.saturationIndex("CSHQ-JenD")/5.167,
        aprops.saturationIndex("CSHQ-JenH")/4.999,aprops.saturationIndex("CSHQ-TobD")/3.166825,
        aprops.saturationIndex("CSHQ-TobH")/3.0001)[0])
gyp_eff=np.append(gyp_eff,aprops.saturationIndex("Gp")[0]/2)
gib_eff=np.append(gib_eff,aprops.saturationIndex("Gbs")[0]/2)
port_eff=np.append(port_eff,aprops.saturationIndex("Portlandite")[0]/3)
strat_eff=np.append(strat_eff,max(aprops.saturationIndex("straetlingite"),aprops.saturationIndex("straetlingite5_5"),
          aprops.saturationIndex("straetlingite7"))[0]/7)
ettrin_eff = np.append(ettrin_eff,max(aprops.saturationIndex("ettringite"),aprops.saturationIndex("ettringite13"),
              aprops.saturationIndex("Ettringite13_des"),aprops.saturationIndex("ettringite30"),
              aprops.saturationIndex("ettringite9"),aprops.saturationIndex("Ettringite9_des"))[0]/15)
mono_eff=np.append(mono_eff,max(aprops.saturationIndex("monosulphate10_5"),aprops.saturationIndex("monosulphate12"),
          aprops.saturationIndex("monosulphate1205"),aprops.saturationIndex("monosulphate14"),
          aprops.saturationIndex("monosulphate16"),aprops.saturationIndex("monosulphate9"))[0]/11)

# At 1 hr

aprops = aprops_3
csh=np.append(csh,max(aprops.saturationIndex("CSH3T-T2C"),aprops.saturationIndex("CSH3T-T5C"),
        aprops.saturationIndex("CSH3T-TobH"),aprops.saturationIndex("CSHQ-JenD"),
        aprops.saturationIndex("CSHQ-JenH"),aprops.saturationIndex("CSHQ-TobD"),
        aprops.saturationIndex("CSHQ-TobH"))[0])
ettrin = np.append(ettrin,max(aprops.saturationIndex("ettringite"),aprops.saturationIndex("ettringite13"),
              aprops.saturationIndex("Ettringite13_des"),aprops.saturationIndex("ettringite30"),
              aprops.saturationIndex("ettringite9"),aprops.saturationIndex("Ettringite9_des"))[0])
gyp=np.append(gyp,aprops.saturationIndex("Gp")[0])
gib=np.append(gib,aprops.saturationIndex("Gbs")[0])
port=np.append(port,aprops.saturationIndex("Portlandite")[0])
strat=np.append(strat,max(aprops.saturationIndex("straetlingite"),aprops.saturationIndex("straetlingite5_5"),
          aprops.saturationIndex("straetlingite7"))[0])
mono=np.append(mono,max(aprops.saturationIndex("monosulphate10_5"),aprops.saturationIndex("monosulphate12"),
          aprops.saturationIndex("monosulphate1205"),aprops.saturationIndex("monosulphate14"),
          aprops.saturationIndex("monosulphate16"),aprops.saturationIndex("monosulphate9"))[0])

csh_eff=np.append(csh_eff,max(aprops.saturationIndex("CSH3T-T2C")/5.5,aprops.saturationIndex("CSH3T-T5C")/5,
        aprops.saturationIndex("CSH3T-TobH")/4.5,aprops.saturationIndex("CSHQ-JenD")/5.167,
        aprops.saturationIndex("CSHQ-JenH")/4.999,aprops.saturationIndex("CSHQ-TobD")/3.166825,
        aprops.saturationIndex("CSHQ-TobH")/3.0001)[0])
gyp_eff=np.append(gyp_eff,aprops.saturationIndex("Gp")[0]/2)
gib_eff=np.append(gib_eff,aprops.saturationIndex("Gbs")[0]/2)
port_eff=np.append(port_eff,aprops.saturationIndex("Portlandite")[0]/3)
strat_eff=np.append(strat_eff,max(aprops.saturationIndex("straetlingite"),aprops.saturationIndex("straetlingite5_5"),
          aprops.saturationIndex("straetlingite7"))[0]/7)
ettrin_eff = np.append(ettrin_eff,max(aprops.saturationIndex("ettringite"),aprops.saturationIndex("ettringite13"),
              aprops.saturationIndex("Ettringite13_des"),aprops.saturationIndex("ettringite30"),
              aprops.saturationIndex("ettringite9"),aprops.saturationIndex("Ettringite9_des"))[0]/15)
mono_eff=np.append(mono_eff,max(aprops.saturationIndex("monosulphate10_5"),aprops.saturationIndex("monosulphate12"),
          aprops.saturationIndex("monosulphate1205"),aprops.saturationIndex("monosulphate14"),
          aprops.saturationIndex("monosulphate16"),aprops.saturationIndex("monosulphate9"))[0]/11)

# At 1d
aprops = aprops_4
csh=np.append(csh,max(aprops.saturationIndex("CSH3T-T2C"),aprops.saturationIndex("CSH3T-T5C"),
        aprops.saturationIndex("CSH3T-TobH"),aprops.saturationIndex("CSHQ-JenD"),
        aprops.saturationIndex("CSHQ-JenH"),aprops.saturationIndex("CSHQ-TobD"),
        aprops.saturationIndex("CSHQ-TobH"))[0])
ettrin = np.append(ettrin,max(aprops.saturationIndex("ettringite"),aprops.saturationIndex("ettringite13"),
              aprops.saturationIndex("Ettringite13_des"),aprops.saturationIndex("ettringite30"),
              aprops.saturationIndex("ettringite9"),aprops.saturationIndex("Ettringite9_des"))[0])
gyp=np.append(gyp,aprops.saturationIndex("Gp")[0])
gib=np.append(gib,aprops.saturationIndex("Gbs")[0])
port=np.append(port,aprops.saturationIndex("Portlandite")[0])
strat=np.append(strat,max(aprops.saturationIndex("straetlingite"),aprops.saturationIndex("straetlingite5_5"),
          aprops.saturationIndex("straetlingite7"))[0])
mono=np.append(mono,max(aprops.saturationIndex("monosulphate10_5"),aprops.saturationIndex("monosulphate12"),
          aprops.saturationIndex("monosulphate1205"),aprops.saturationIndex("monosulphate14"),
          aprops.saturationIndex("monosulphate16"),aprops.saturationIndex("monosulphate9"))[0])

csh_eff=np.append(csh_eff,max(aprops.saturationIndex("CSH3T-T2C")/5.5,aprops.saturationIndex("CSH3T-T5C")/5,
        aprops.saturationIndex("CSH3T-TobH")/4.5,aprops.saturationIndex("CSHQ-JenD")/5.167,
        aprops.saturationIndex("CSHQ-JenH")/4.999,aprops.saturationIndex("CSHQ-TobD")/3.166825,
        aprops.saturationIndex("CSHQ-TobH")/3.0001)[0])
gyp_eff=np.append(gyp_eff,aprops.saturationIndex("Gp")[0]/2)
gib_eff=np.append(gib_eff,aprops.saturationIndex("Gbs")[0]/2)
port_eff=np.append(port_eff,aprops.saturationIndex("Portlandite")[0]/3)
strat_eff=np.append(strat_eff,max(aprops.saturationIndex("straetlingite"),aprops.saturationIndex("straetlingite5_5"),
          aprops.saturationIndex("straetlingite7"))[0]/7)
ettrin_eff = np.append(ettrin_eff,max(aprops.saturationIndex("ettringite"),aprops.saturationIndex("ettringite13"),
              aprops.saturationIndex("Ettringite13_des"),aprops.saturationIndex("ettringite30"),
              aprops.saturationIndex("ettringite9"),aprops.saturationIndex("Ettringite9_des"))[0]/15)
mono_eff=np.append(mono_eff,max(aprops.saturationIndex("monosulphate10_5"),aprops.saturationIndex("monosulphate12"),
          aprops.saturationIndex("monosulphate1205"),aprops.saturationIndex("monosulphate14"),
          aprops.saturationIndex("monosulphate16"),aprops.saturationIndex("monosulphate9"))[0]/11)


# At 7d
aprops = aprops_5
csh=np.append(csh,max(aprops.saturationIndex("CSH3T-T2C"),aprops.saturationIndex("CSH3T-T5C"),
        aprops.saturationIndex("CSH3T-TobH"),aprops.saturationIndex("CSHQ-JenD"),
        aprops.saturationIndex("CSHQ-JenH"),aprops.saturationIndex("CSHQ-TobD"),
        aprops.saturationIndex("CSHQ-TobH"))[0])
ettrin = np.append(ettrin,max(aprops.saturationIndex("ettringite"),aprops.saturationIndex("ettringite13"),
              aprops.saturationIndex("Ettringite13_des"),aprops.saturationIndex("ettringite30"),
              aprops.saturationIndex("ettringite9"),aprops.saturationIndex("Ettringite9_des"))[0])
gyp=np.append(gyp,aprops.saturationIndex("Gp")[0])
gib=np.append(gib,aprops.saturationIndex("Gbs")[0])
port=np.append(port,aprops.saturationIndex("Portlandite")[0])
strat=np.append(strat,max(aprops.saturationIndex("straetlingite"),aprops.saturationIndex("straetlingite5_5"),
          aprops.saturationIndex("straetlingite7"))[0])
mono=np.append(mono,max(aprops.saturationIndex("monosulphate10_5"),aprops.saturationIndex("monosulphate12"),
          aprops.saturationIndex("monosulphate1205"),aprops.saturationIndex("monosulphate14"),
          aprops.saturationIndex("monosulphate16"),aprops.saturationIndex("monosulphate9"))[0])

csh_eff=np.append(csh_eff,max(aprops.saturationIndex("CSH3T-T2C")/5.5,aprops.saturationIndex("CSH3T-T5C")/5,
        aprops.saturationIndex("CSH3T-TobH")/4.5,aprops.saturationIndex("CSHQ-JenD")/5.167,
        aprops.saturationIndex("CSHQ-JenH")/4.999,aprops.saturationIndex("CSHQ-TobD")/3.166825,
        aprops.saturationIndex("CSHQ-TobH")/3.0001)[0])
gyp_eff=np.append(gyp_eff,aprops.saturationIndex("Gp")[0]/2)
gib_eff=np.append(gib_eff,aprops.saturationIndex("Gbs")[0]/2)
port_eff=np.append(port_eff,aprops.saturationIndex("Portlandite")[0]/3)
strat_eff=np.append(strat_eff,max(aprops.saturationIndex("straetlingite"),aprops.saturationIndex("straetlingite5_5"),
          aprops.saturationIndex("straetlingite7"))[0]/7)
ettrin_eff = np.append(ettrin_eff,max(aprops.saturationIndex("ettringite"),aprops.saturationIndex("ettringite13"),
              aprops.saturationIndex("Ettringite13_des"),aprops.saturationIndex("ettringite30"),
              aprops.saturationIndex("ettringite9"),aprops.saturationIndex("Ettringite9_des"))[0]/15)
mono_eff=np.append(mono_eff,max(aprops.saturationIndex("monosulphate10_5"),aprops.saturationIndex("monosulphate12"),
          aprops.saturationIndex("monosulphate1205"),aprops.saturationIndex("monosulphate14"),
          aprops.saturationIndex("monosulphate16"),aprops.saturationIndex("monosulphate9"))[0]/11)

# %%
# =============================================================================
# Plot the result
# =============================================================================

x=np.arange(1,7,1)

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial"]
plt.rcParams["font.size"]=10
plt.rcParams["font.style"]="italic"
plt.rcParams["font.weight"]="bold"
plt.rcParams["lines.linewidth"]=0.75

# define axis behavior

fig , ax1 = plt.subplots(figsize=[4.7,3])
ax2 = ax1.twinx()
ax2.spines["right"].set_color("red")
ax2.tick_params(axis='y', colors='red') 
ax3 = ax1.twinx()
ax3.spines["right"].set_position(("axes", 1.22))
ax3.spines["right"].set_color("blue")
ax3.tick_params(axis='y', colors='blue') 

# plot result

line1=ax1.plot(x,csh_eff,'ko--',markersize=3,fillstyle='none',markeredgewidth=0.5,label="C-S-H",mfc='k')
line2=ax1.plot(x,ettrin_eff,'kd--',markersize=3,fillstyle='none',markeredgewidth=0.5,label="Ettringite",mfc='k')
line3=ax1.plot(x,strat_eff,'kx--',markersize=3,label="Stratlingite")

line4=ax1.plot(x,gyp_eff,'o--',color='0.5',markersize=3,label="Gypsum")
line5=ax1.plot(x,gib_eff,'v--',color='0.5',markersize=3,label="Gibbsite")
line6=ax1.plot(x,port_eff,'^--',color='0.5',markersize=3,label="Portlandite")
line7=ax1.plot(x,mono_eff,'x-',color='0.5',markersize=3,label="Monosulfate")

ax1.set_xlabel("Curing duration",weight="bold")
ax1.set_xlim([0.5,6.5])
ax1.set_xticks(np.arange(1,7,1),labels=np.array(['0 min','15 min','30 min','1 h','1 d','7 d']))
ax1.set_ylim([-6,2])
ax1.set_ylabel("Effective saturation index",weight="bold")
# ax1.legend(ncol=2,fontsize=8,loc="lower right",frameon=False,borderpad=0.1)

ax1.grid(lw=0.1,alpha=.7)

x2 = np.array([1,4,5,6])
line8=ax2.plot(x2,strength*100,'r--.',label="Strength",markeredgewidth=0.5,fillstyle='none')
ax2.set_ylim([0,300])
ax2.set_ylabel("Relative UCS (%)",weight="bold",color='r')
# ax2.legend(ncol=2,fontsize=8,loc="lower right",frameon=False,borderpad=0.1)

line9=ax3.plot(x,df['pH'],'b--.',fillstyle='none',markeredgewidth=0.5,label="pH")
ax3.set_ylabel("pH",weight="bold",color='b')
ax3.set_ylim([4,12])
# ax3.legend(ncol=2,fontsize=8,loc="right",frameon=False,borderpad=0.1)


line10=ax1.plot([4],[ettrin_eff[3]],'rd',markersize=3,label="Observed")

line = line1 + line2 + line3 + line4 + line5+ line6 + line7 + line8 + line9 + line10
labels = [l.get_label() for l in line]
ax1.legend(line, labels, ncol=2,fontsize=8,loc="lower right",frameon=False,borderpad=0.1)
# ax1.grid(lw=0.1,alpha=0.8)
plt.tight_layout()

# plt.savefig("stab_b.svg")