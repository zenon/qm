# load with include("qm.jl")

# infrastructure
using Plots, UnicodePlots, StatsPlots, Distributions, Primes



using AlgebraicNumbers

# As AlgebraicNumbers loads Nemo.jl, and this prints a greeting
# we clear the screen.
# (No, I don't love this either.)
# This is a POSIX command, btw.
print("\e[1;1H\e[2J")
# print( "\033[2J" ) # clears without moving the cursor.


n(v) = real(v'*v)

# Pauli matrices

zero = AlgebraicNumber(0)
# can't call it one, as that is a base function
uno  = AlgebraicNumber(1)
i    = sqrt(AlgebraicNumber(-1))

# Pauli Matrices
σ1 = [[zero,  uno],
      [uno , zero]]

σ2 = [[zero,  -i],
      [i   , zero]]

σ3 = [[uno,   zero],
      [zero , -uno]]


function printSingleQB(q)
    α  = q1' * q 
    β  = q2' * q
    lu = '╔'
    u  = '═'
    mu = '╤'
    ru = '╗'
    l  = '║'
    m  = '│'
    r  = l
    ld = '╚'
    d  = u
    md = '╧'
    rd = '╝'
    upper = lu*u*u*u*mu*u*u*u*mu*u*u*u*ru
    middle = l*' '*'0'*' '*m*' '*'1'*' '*m*' '*'?'*' '*r
    lower = ld*d*d*d*md*d*d*d*md*d*d*d*rd
    println(upper)
    println(middle)
    println(lower)
end

#=
'─': Unicode U+2500
'│': Unicode U+2502

'┌': Unicode U+250c
'┐': Unicode U+2510
'└': Unicode U+2514
'┘': Unicode U+2518
'├': Unicode U+251c
'┤': Unicode U+2524
'┬': Unicode U+252c
'┴': Unicode U+2534
'┼': Unicode U+253c
'┴': Unicode U+2534
'┼': Unicode U+253c
'═': Unicode U+2550
'║': Unicode U+2551
═║╒╓╔╕╖╗╘╙╚╛╜╝╞╟
'╟': Unicode U+255f
'╠': Unicode U+2560
╠╡╢╣╤╥╦╧╨╩╪╫╬
'╬': Unicode U+256c

'○': Unicode U+25cb
'◌': Unicode U+25cc
'●': Unicode U+25cf

=#


# override show(io,mime::MIME...,x) with the most useful interactive output with a given multimedia capability (MIME…). E.g. for the REPL version use mime::MIME"text/plain". This will be used for display.

# hcat([1 2; 3 4], [5 6 7; 8 9 10], [11, 12])
# vcat([1 2; 3 4], [5 6; 7 8; 9 10], [11 12])



import Base.show


# needs Primes
# factorizing might not be the most efficient way to do this
# returns (rootOfSquarePart, squareFreePart)
function squarePart(n :: BigInt)
    squarePart = big(1)
    squareFree = big(1)
    for (p,e) in factor(n)
        ef = mod(e,2)
        es = div(e - ef, 2)
#        println(p," ",e," ",es," ",ef)
        if ef != 0
            squareFree = squareFree*p
        end
        if es != 0
            squarePart = squarePart*(p^es)
        end        
    end
    (squarePart, squareFree)
end


# better
function show(io::IO, x::Rational{BigInt})
    show(io, numerator(x))
    if denominator(x) != 1
        print(io, "/")
        show(io, denominator(x))
    end
end

function show(io::IO, an::AlgebraicNumber)
    if length(an.coeff) > 3
print(io,"≈")
#ndigits = max(10, round(Int,ceil(convert(Float64,log(an.prec)/log(10)))))
show(io,convert(Complex{Float64},an.apprx))
#print(io,"...")
    else
        if length(an.coeff) < 2        ## should not happen
            show(io,"unexpected coeffs: $(an.coeff)")
        elseif length(an.coeff) == 2   ## rational number
            if an.coeff[2] == 1
                show(io, -an.coeff[1])
            else
                show(io, -(an.coeff[1])//an.coeff[2])
            end
        else                           ## surd
            #
            #
            #
            c, b, a = an.coeff
            d = b*b - (4*a*c)
            ds, df = squarePart(d)
            dfSign = sign(df)
            dfAbs = abs(df)
            r = -b // 2a
            if r!=0
                show(io,r)
            end
            f = ds//(2a)
            # which of the solutions is intended?
            x  = sqrt(dfAbs)*f*(if dfSign == -1; im; else; 1; end)
            x1 = x+r
            x2 = -x+r
            #            println("a=$a, b=$b, c=$c, d=$d, ds=$ds, df=$df, dfSign=$dfSign, dfAbs=$dfAbs, r=$r, f=$f, x=$x, x1=$x1, x2=$x2")
            # open: f can be negative, and this switches the signs.
            if abs(x1-an.apprx) < an.prec
                if r != 0
                    print(io, "+")
                end
            elseif abs(x2-an.apprx) <= an.prec # danger. prec can be zero! e.g. i+i
                print(io, "-")
            else
                println("\n could not identify root.")
            end
            if f!= 1
                show(io, f)
                print(io,"*")
            end
            printedSurd = false
            if dfAbs != 1
                print(io, "√")
                show(io,dfAbs)
                printedSurd = true
            end
            if dfSign == -1
                if printedSurd
                    print(io, "*i")
                else
                    print(io, "i")
                end
            end
        end
    end
end


#=

Bugs:

1+i -> 2

=#




