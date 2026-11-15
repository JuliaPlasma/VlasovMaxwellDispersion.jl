@testitem "Maxwell-Juttner rejects damped ω at k∥ = 0" begin
    species = NormalizedPopulation(-1.0, 0.2, MaxwellJuttner(; theta = 1 / 1.0e4))
    @test_throws ArgumentError contribution(species, 0.8 - 0.01im, Wavevector(0.1, 0.0))
end

@testitem "Maxwell-Juttner tends to Maxwellian for large mu" begin
    let mu = 1.0e4, k = Wavevector(0.1, 0.2), ω = 0.8 + 0.01im
        χrel = contribution(NormalizedPopulation(-1.0, 0.2, MaxwellJuttner(; theta = 1 / mu)), ω, k)
        χcoldish = contribution(NormalizedPopulation(-1.0, 0.2, BiMaxwellian(; vth = sqrt(2 / mu))), ω, k)
        @test χrel ≈ χcoldish rtol = 1.0e-3
    end
end

@testitem "Maxwell-Juttner Bessel K stays exact at large |z|" begin
    using SpecialFunctions: besselkx
    kx = VlasovMaxwellDispersion._besselkx
    for ν in (2, 3), z in (150.0 + 40im, 1.0e7 + 1.0e6im, 1.0e9)
        @test kx(ν, z) ≈ besselkx(ν, z) rtol = 1.0e-14
    end
    @test kx(2, 1.0e12) ≈ sqrt(π / 2.0e12) * (1 + 15 / 8.0e12)  # AMOS throws here
end
