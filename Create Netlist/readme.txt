THIS PROGRAM DOES NOT GENERATE DATA ABOUT THE JUNCTIONS OR FLUX... YET 
THE MUTUAL INDUCTANCE INPUT FROM USER IS 1/Ma SO IN THE NETLIST WRITE 1/(USER_INPUT) TO ACCOUNT FOR THIS

-------------------------------------------------------------------------------------------------------------------------------

When entering mutually inducted loops enter in the form '1,2;5' where loop 1 and 2 are coupled with mutaul inductance 5μA/Φ𝜊

When defining components in a loop seperate by line and use the following naming conventions:
Resistor names should start with R		e.g. "R1", "Ra", "Resistor1"
Capacitor names should start with C		e.g. "C1", "Ca", "Capacitor1"
Inductor names should start with L		e.g. "L1", "La", "L_Inductor1"
Josephson Junction names should start with J	e.g. "J1", "Ja", "Josephson1"

-------------------------------------------------------------------------------------------------------------------------------

Do not include units or prefixes when defining component parameters.     <--- Use of prefixes may be included in future version

-------------------------------------------------------------------------------------------------------------------------------

The file loops+components.csv contains the structure of the circuit in terms of loops where each loop is separated by line and 
each component of a loop is separated by commas.

For example:
Loop 0,Ib,J1		<--- Loop 0 (Loop Ib) contains components Ib and J1
Loop 1,J1,J2		<--- Loop 1 contains components J1 and J2
Loop 2,J2,Ra,La		<--- Loop 3 contains components J2, Ra and La

-------------------------------------------------------------------------------------------------------------------------------

The file component_parameters.txt contains data about each relevant component and it's parameters, examples for components are 
shown below.

Data for RESITORS is stored in the form:
	RESISTOR_NAME: Any["RESISTANCE"]

Data for CAPACITORS is stored in the form:
	CAPACITOR_NAME: Any["CAPACITANCE"]

Data for INDUCTORS is stored in the form:
	INDUCTOR_NAME: Any["INDUCTANCE"]

Data for JOSEPHSON JUNCTIONS is stored in the form:
	JOSEPHSON_JUNCTION_NAME: Any["CRITICAL_CURRENT, SHUNT_RESISTANCE, STEWART-MCCUMBER_PARAMETER, HALF_LOOP_INDUCTANCE"]

-------------------------------------------------------------------------------------------------------------------------------