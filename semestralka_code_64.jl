# 1:59
using CairoMakie, LinearAlgebra, Trapz

heaviside(x::AbstractFloat) = ifelse(x < 0, zero(x), ifelse(x > 0, one(x), oftype(x,0.5))); #heavisadeovy funkce
jama(x, sirka, vyska) = -oftype(x,vyska).*(heaviside.(x+oftype(x,0.5*sirka))-heaviside.(x-oftype(x,0.5*sirka))); #pravouhla jama dane sirky a vysky pro pole souradnic x
N = 2000; # počet dílků
rozsah = 3;
x = LinRange(-rozsah, rozsah, N); # vytvoří pole v daném rozsahu s N prvky
# Pole je převedeno na typ Float32 pro úsporu výpočetního výkonu
dx = x[2] - x[1];
ħ = 1.0;
m = 1.0;

sirkaJamy = 1.7e0
hloubkaJamy = 80e0

U = jama.(x, sirkaJamy, hloubkaJamy); # vytvoří pole potenciálu
Δ = Tridiagonal(ones(N-1), -2*ones(N), ones(N-1) )/(dx.^2); # 1DLaplaceúv operátor
# Tridiagonal umístí příslušné vektory na hlavní a první vedlejší diagonály
H = (-0.5*ħ^2/m).*Δ + diagm(U) # Hamiltonián - matice je vypsána

energie, funkce = eigen(H) # najde vlastní honoty a vlastní vektory

fig = Figure(backgroundcolor = RGBf(1, 1, 1),
	resolution = (500, 500))
ax1 = Axis(fig[1, 1], ylabel = L"Vlnové funkce $\Psi_{i}$")
ax2 = Axis(fig[2, 1], ylabel = L"Hustoty pravděpodobnosti $\left|\Psi_i\right|^2$")
for i=1:3
	text = "E_$i = " * string(round(energie[i], digits=3))
	lines!(ax1, x, (funkce[:, i]), label = text)
	lines!(ax2, x, (funkce[:, i]).^2, label = text)
end
axislegend.([ax1, ax2])
save("figure.pdf", fig)
save("figure.png", fig, px_per_unit = 3)

grJamy = Figure(backgroundcolor = RGBf(1, 1, 1), resolution = (16*30, 9*30))
ax = Axis(grJamy[1,1], ylabel = L"V(x)", xlabel = L"x")
lines!(ax, x, U);
save("jama.pdf", grJamy)
save("jama.png", grJamy, px_per_unit = 3)

# stavKVykresleni = 40
# fig1 = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
# 	resolution = (1000, 1000))
# ax11 = Axis(fig1[1, 1], ylabel = L"Vlnové funkce $\Psi_{i}$")
# ax12 = Axis(fig1[2, 1], ylabel = L"Hustoty pravděpodobnosti $\left|\Psi_i\right|^2$")
# lines!(ax11, x, (funkce[:, stavKVykresleni]), label = "E = " * string(energie[stavKVykresleni]))
# lines!(ax12, x, (funkce[:, stavKVykresleni]).^2, label = "E = " * string(energie[stavKVykresleni]))
# axislegend.([ax11, ax12])
# fig1

