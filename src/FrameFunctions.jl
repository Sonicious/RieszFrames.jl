"""
Held(x,[,m])

generates the Held-Wavelet with the smoothness m

The Wavelet is constructed according to the paper
S. Held, M. Storath, P. Massopust and B. Forster, "Steerable Wavelet Frames Based on the Riesz Transform," in IEEE Transactions on Image Processing, vol. 19, no. 3, pp. 653-667, March 2010, doi: 10.1109/TIP.2009.2036713.
https://ieeexplore.ieee.org/document/5339191
"""
function Held(x, m=4)
    return cos(2 * pi * q(x, m)) * ((1 / 8 < x <= 1 / 4)) + sin(2 * pi * q(1 / 2 * x, m)) * ((1 / 4 < x <= 1 / 2))
end
# Support function for held wavelet
function q(x, m)
    # polynomial for the wavelet
    # gives values of the supporting function q in dependency of m
    if m == 0
        return 1 / 8 * 4 - 1 / 8 * (4 * 4) * x
    elseif m == 1
        return -1 / 4 * 4 + 3 / 2 * (4 * 4) * x - 9 / 4 * (4 * 4 * 4) * (x * x) + 1 * (4 * 4 * 4 * 4) * (x * x * x)
    elseif m == 2
        return 2 * 4 - 15 * 4^2 * x + 45 * 4^3 * x .^ 2 - 65 * 4^4 * x .^ 3 + 45 * 4^5 * x .^ 4 - 12 * 4^6 * x .^ 5
    elseif m == 3
        return return -13 * 4 + 140 * 4^2 * x - 630 * 4^3 * x .^ 2 + 1540 * 4^4 * x .^ 3 - 2205 * 4^5 * x .^ 4 + 1848 * 4^6 * x .^ 5 - 840 * 4^7 * x .^ 6 + 160 * 4^8 * x .^ 7
    elseif m == 4
        return 92 * 4 - 1260 * 4^2 * x + 7560 * 4^3 * x .^ 2 - 26040 * 4^4 * x .^ 3 + 56700 * 4^5 * x .^ 4 - 80892 * 4^6 * x .^ 5 + 75600 * 4^7 * x .^ 6 - 44640 * 4^8 * x .^ 7 + 15120 * 4^9 * x .^ 8 - 2240 * 4^10 * x .^ 9
    else
        throw(DomainError(m, "the smoothness factor m must be 0 <= m <= 4"))
    end
end

"""
Papadakis(x)

generates the Papadakis Wavelet Frame
"""
function Papadakis(x)
    return sqrt((1 + sin(2 * pi * 5 * x)) / 2) * (3 / 20 < x <= 1 / 4) + (1 / 4 < x <= 3 / 10) + sqrt((1 - sin(pi * 5 * x)) / 2) * (3 / 10 < x <= 1 / 2)
end

"""
Shannon(x)

generates theShannon Wavelet Frame
"""
function Shannon(x)
    return (1 / 4 < x <= 1 / 2)
end

"""
Simoncelli(x)

generates the Simoncelli Wavelet Frame
"""
function Simoncelli(x)
    return (1 / 8 < x <= 1 / 2) ? cos(pi / 2 * log2(4 * x)) : zero(x)
end