### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 2762bbd1-b4c9-45e2-883b-9839258e00ba
begin
	using Pkg
	Pkg.activate(".")
end

# ╔═╡ 188e11bf-20e0-44bb-b9fc-2694f4b024da
using CairoMakie, Trapz

# ╔═╡ 81c6679a-ecff-42a1-a858-01d0c8a3cede
using LinearAlgebra

# ╔═╡ b5004efe-cb98-4488-9f4a-15f680355f19
using PlutoUI

# ╔═╡ 637fb152-463f-11ec-3f28-2bdda77d61e9
md"""
# Sample codes from the material

See: Garcia, R., Zozulya, A., &#38; Stickney, J. (2007). _MATLAB codes for teaching quantum physics: Part 1_.

## Code 1: Numerical integration and derivation, plotting

"""

# ╔═╡ 343cb553-2aa9-4ae5-95bb-6f2fdb80e3bb
begin
	N = 2000
	L = 20
	x = Float32.(LinRange(-L, L, N))
	dx = x[2] - x[1]
	#y = exp(-x.^2/16)
	A = 0.5
	y = (sin.(x*A).^2)./(pi*(A*x.^2)) # representing delta function
	y = sin.(x*A).^2 
	y1 = exp.(-x.^(2))
	lines(x, y) 
	lines!(x, real(y))
	lines!(x, imag(y))
	lines!(x, y1)
	current_figure()
end

# ╔═╡ 69d156b9-976a-431c-a53f-8420ff8e3fd3
integration = sum(y)*dx

# ╔═╡ 879a869a-4d91-4efa-93ff-3641d9b5cdae
anotherApproach = trapz(x, y)

# ╔═╡ 3c55c19c-f338-44b9-af4b-78101d01080a
md"""
## Representing differential operators with matrieces

Trying represent operator $\hat D$ acting on ket (__column__) vector $\left|f\right>$ with an matrix.

The MATLAB commands as `ones` or `diag` should be used.

"""

# ╔═╡ dc8e7b1a-44dd-476f-a652-146f6fecea76
ones(Float32, 2,4)

# ╔═╡ ddaf1397-a5bf-41a4-ad04-3b9a3c216dea
begin
	vec = [1, 2, 3, 4]
	Diagonal(vec)
end

# ╔═╡ 949af4b7-f6e9-4f65-aa7c-fc1abb1a93e1
md"""
> An exercise we suggest is for students to verify that the derivative matrix is not Hermitian while the derivative matrix times the imaginary number $i$ is. This can be very valuable for promoting student understanding if one in conjunction with the proof usually given for the differential operator

For the Julia implemetation the `Bidiagonal` command could be helpful.
"""

# ╔═╡ 4fd25b31-304d-4a7c-84cd-cca026fc1065
difMatrix = Bidiagonal(ones(Int16, N), -ones(Int16, N-1), :U)/dx

# ╔═╡ cff32a65-7856-486d-bbc1-0373986a8d2d
centreDif = Tridiagonal(ones(Int16, N-1), zeros(Int16, N), -ones(Int16, N-1), )/(2*dx)

# ╔═╡ 87846b6e-860c-4ecf-af17-79416865bb18
Δ = Tridiagonal(ones(N-1), -2*ones(N), ones(N-1) )/(dx.^2)

# ╔═╡ 0db82ac5-5888-4b48-89ec-68399f49256b
begin
	lines(x, y)
	lines!(x[1:N-1], (difMatrix*y[:])[1:N-1])
	lines!(x[2:N-1], (centreDif*y[:])[2:N-1])
	lines!(x[2:N-1], (Δ*y[:])[2:N-1])
	current_figure()
end

# ╔═╡ 76fe08e0-b157-41b9-b34b-5fab71f91606
md"""

## Infinite square well

Means to solve the Schrödinger equation 

$$\hat{H} \left|\Psi\right> = E\left|\Psi\right>$$

This means search for eigenvalues $E$ and eigenvectors $\left|\Psi\right>$ of matrix $\hat{H}$.

1. Let's write the matrix representation of $\hat{H}$
$$H = -\frac{\hbar^2}{2m}\Delta$$
"""

# ╔═╡ 2e44bb7c-1955-40ec-bfbb-a4db31420867
begin
	ħ = 1;
	m = 1;
	H = Float32.(Matrix(-0.5*ħ^2/m .*Δ))
end

# ╔═╡ 2b6046d5-3d6e-4131-9bcf-53277c10adf8
vals, vecs = eigen(H);

# ╔═╡ 2c2259c7-b9d1-4893-bf0c-918d2d8243ef
begin
	for i=1:6
		if i==1
			lines(x, vecs[:, i])
		else
			lines!(x, vecs[:, i])
		end
	end
	current_figure()
end

# ╔═╡ b6a4a876-78e7-4b0e-8bbb-f7407aae8eb5
sum(vecs[:,4].*vecs[:, 2])

# ╔═╡ 722b3778-0b1d-44db-a12b-a50817530b5a
md"""
arbitary potentials to come 🙂
"""

# ╔═╡ d7746175-97ec-43e1-b3b8-cc1d9fa48af4
md"""
## Arbitary potentials

The program finds eigenmodes and eigenenergies for Schrödinger equation 

$$- \frac{\hbar^2}{2m} \frac{\mathrm{d^2}}{\mathrm{d}x^2} \Psi_m(x) + U(x) \Psi_m(x) = E\Psi_m(x)$$

Where $\Psi_m(x)$ are the eigenmodes and $E$ are the eigenenergies.

"""

# ╔═╡ 5bb93743-9b96-47e7-84e0-ab47942f7363
begin 
	pom = rand(100)
	indi = sortperm(pom)
	pomSorted = sort(pom)
	indi
end

# ╔═╡ 686b0e60-c68f-48f7-821e-7d4206a8460e
md"""
Slider pro zvolení stavu k vykreslení na následujícím grafu.

$(@bind stav PlutoUI.Slider(1:1:N, show_value = true))
"""

# ╔═╡ 5e49a80b-754c-45d6-91fe-8b9915ae27ba
md"""
### Same figure for probability density
"""

# ╔═╡ 34764fa5-adab-4c8c-91dd-9f8c9e1e5b43
md"""
### Vytvoření heaviside funkce  jámy s pravoúhlým potenciálem
"""

# ╔═╡ ebab8409-587f-4056-863b-b49adb569ef5
heaviside(x::AbstractFloat) = ifelse(x < 0, zero(x), ifelse(x > 0, one(x), oftype(x,0.5)));

# ╔═╡ a7244b82-b9ed-423f-b199-1ceff58d50bf
function jama(x, sirka, vyska)
	return -vyska.*(heaviside.(x+sirka)-heaviside.(x-sirka))
end

# ╔═╡ a4908d66-6d18-4968-b659-70033014fea0
begin
	#U = 1/2*x.^4 # defining the potential
	U = jama.(x, 15, 1000)
	e = ones(N,1)
	Hpot = -1/2*ħ^2/m*Δ + diagm(U)
	eng, modes = eigen(Hpot)
	ind = sortperm(eng) # vytvoří pole indexů uspořádání jedn
	engSorted = eng[ind]
	modesSorted = modes[ind, :]
end

# ╔═╡ 3332c780-e8f8-4590-807b-deb67a210182
lines(x, modes[:, stav])

# ╔═╡ ec13678e-0630-419f-b158-c7508f2faf66
begin
	for i=1:3
		if i==1
			lines(x, modes[:, i], label = "1. stav")
		else
			text = string(i) * ". stav" 
			lines!(x, modes[:, i], label = text)
		end
	end
	axislegend()
	current_figure()
end

# ╔═╡ a6d8c178-edc4-4570-8363-84630f4d66d3
begin
	for i=1:3
		if i==1
			lines(x, (modes[:, i]).^2, label = "1. stav")
		else
			text = string(i) * ". stav" 
			lines!(x, (modes[:, i]).^2, label = text)
		end
	end
	axislegend()
	current_figure()
end

# ╔═╡ deddb979-7a5d-42ec-b2de-e7e7c767b258
modes[:,1]

# ╔═╡ 38a00213-3cbb-4cc4-84fc-2bac8c55ed6f
@bind pocet PlutoUI.Slider(2:1:2000, show_value = true)

# ╔═╡ b7efbd75-536b-4401-b0c2-545db82e0e9c
@bind sir PlutoUI.Slider(0:0.001:1, show_value = true)

# ╔═╡ c0366d0d-8fb8-4c30-8adf-570f2f2bccd8
begin
	promenna = Float32.(LinRange( -1, 1, pocet))
	lines(promenna, heaviside.(promenna))
	lines(promenna, jama.(promenna, sir, 3))
end

# ╔═╡ 8907d5f0-bf10-435c-bce9-a2dc60e4970b
md"""
### Chování částice s v jámě s pravoúhlým potenciálem
"""

# ╔═╡ Cell order:
# ╟─637fb152-463f-11ec-3f28-2bdda77d61e9
# ╠═2762bbd1-b4c9-45e2-883b-9839258e00ba
# ╠═188e11bf-20e0-44bb-b9fc-2694f4b024da
# ╠═343cb553-2aa9-4ae5-95bb-6f2fdb80e3bb
# ╠═69d156b9-976a-431c-a53f-8420ff8e3fd3
# ╠═879a869a-4d91-4efa-93ff-3641d9b5cdae
# ╟─3c55c19c-f338-44b9-af4b-78101d01080a
# ╠═81c6679a-ecff-42a1-a858-01d0c8a3cede
# ╠═dc8e7b1a-44dd-476f-a652-146f6fecea76
# ╠═ddaf1397-a5bf-41a4-ad04-3b9a3c216dea
# ╟─949af4b7-f6e9-4f65-aa7c-fc1abb1a93e1
# ╠═4fd25b31-304d-4a7c-84cd-cca026fc1065
# ╠═cff32a65-7856-486d-bbc1-0373986a8d2d
# ╠═87846b6e-860c-4ecf-af17-79416865bb18
# ╠═0db82ac5-5888-4b48-89ec-68399f49256b
# ╟─76fe08e0-b157-41b9-b34b-5fab71f91606
# ╠═2e44bb7c-1955-40ec-bfbb-a4db31420867
# ╠═2b6046d5-3d6e-4131-9bcf-53277c10adf8
# ╠═2c2259c7-b9d1-4893-bf0c-918d2d8243ef
# ╠═b6a4a876-78e7-4b0e-8bbb-f7407aae8eb5
# ╟─722b3778-0b1d-44db-a12b-a50817530b5a
# ╟─d7746175-97ec-43e1-b3b8-cc1d9fa48af4
# ╠═5bb93743-9b96-47e7-84e0-ab47942f7363
# ╠═a4908d66-6d18-4968-b659-70033014fea0
# ╠═b5004efe-cb98-4488-9f4a-15f680355f19
# ╟─686b0e60-c68f-48f7-821e-7d4206a8460e
# ╠═3332c780-e8f8-4590-807b-deb67a210182
# ╠═ec13678e-0630-419f-b158-c7508f2faf66
# ╟─5e49a80b-754c-45d6-91fe-8b9915ae27ba
# ╟─a6d8c178-edc4-4570-8363-84630f4d66d3
# ╠═deddb979-7a5d-42ec-b2de-e7e7c767b258
# ╟─34764fa5-adab-4c8c-91dd-9f8c9e1e5b43
# ╠═ebab8409-587f-4056-863b-b49adb569ef5
# ╠═a7244b82-b9ed-423f-b199-1ceff58d50bf
# ╠═38a00213-3cbb-4cc4-84fc-2bac8c55ed6f
# ╠═b7efbd75-536b-4401-b0c2-545db82e0e9c
# ╠═c0366d0d-8fb8-4c30-8adf-570f2f2bccd8
# ╟─8907d5f0-bf10-435c-bce9-a2dc60e4970b