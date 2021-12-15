# SQUID simulation netlist compiler, reads netlist file in and produces system of differential equations.
using DelimitedFiles
function netlist_read(name, path)
    f_nl = open(path * '/' * name * "_nl.csv", "r")               #open file containing TOPOLOGY PARAMETERS OF SQUID APF CIRCUIT
    f_cp = open(path * '/' * name * "_cp.txt", "r")  #open file containing COMPONENT PARAMETERS OF SQUID APF CIRCUIT PARAMETER
    f_lc = open(path * '/' * name * "_lc.csv", "r")      #open file containing LOOPS AND COMPONENTS CONTAINED IN EACH LOOP

    nl = readdlm(f_nl, ',', String, '\n'); close(f_nl)
    cp = readdlm(f_cp, ':', String, '\n'); close(f_cp)
    lc = readdlm(f_lc, ',', String, '\n'); close(f_lc)
    return cp


    #L̂ = map(x->(v = tryparse(Float64,x); isnothing(v) ? 0.0 : get(v)),a[2:end,2:end])
end
