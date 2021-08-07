# critical_ra.jl
#
# Analytic expression for the critical Rayleigh (2D linear convection) as a function of horizontal wavenumber and aspect ratio
# 
# Note: smallest critical Ra is achieved with vertical wavenumber m=1

using Plots
using LaTeXStrings

m=1;
nx=100;

# Function of aspect ratio, a
a = LinRange(0.1, 10.0, nx)
n = [1.0, 3.0, 10.0]
nn = size(n)[1]
Ra_c_1 = fill(0.0, nn, nx)

for i=1:1:nn
    Ra_c_1[i,:] = (pi./a).^4 .* (n[i]^2 .+ (a.*m).^2).^3 ./ n[i]^2
end

l = @layout [a;b]

p1 = plot(a, 
    Ra_c_1',
    yaxis=:log,
    ylims=(1,10^10), 
    label=[L"n=1" L"n=3" L"n=10"], lw=2)
xlabel!(L"a")
ylabel!(L"Ra_c")

# Function of horizontal wavenumber, n
n = range(1,10,length=10)
a = [0.1, sqrt(2), 5, 10]
na = size(a)[1]
Ra_c_2 = fill(0.0, na, 10)

for j=1:1:10
    for i=1:1:na
        Ra_c_2[i,j] = (pi/a[i])^4 * (n[j]^2 + (a[i]*m)^2)^3 / n[j]^2
    end
end

# Plot
p2 = plot(n,
    Ra_c_2',
    xticks=n,
    ylim=(1,10^10),
    yaxis=:log,
    label=[L"a=0.1" L"a=\sqrt{2}" L"a=5" L"a=10"],
    lw=2,
    marker = (:square, 8, 0.6, :green, stroke(2, 0.2, :black, :dot)))
xlabel!(L"n")
ylabel!(L"Ra_c")

p3 = plot(p1,p2,layout=l,size=(720,512))


# savefig(p3,string(@__DIR__, "/docs/figures/critical_ra.png"))