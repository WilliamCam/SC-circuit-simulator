mutable struct Resistor
    R::Float64
end

struct Capacitor
    C::Float64
end

struct Inductor
    L::Float64
end

struct RSCJ
    R::Float64
    C::Float64
    I₀::Float64
end

Components = Union{Resistor, Capacitor, Inductor, RSCJ}

struct Loop
    components::Vector{Components}
