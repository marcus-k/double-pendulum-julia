using DifferentialEquations
using ProgressBars
using Printf
using Plots

g = 9.81    # acceleration due to gravity, in m/s^2
l₁ = 1.0    # length of pendulum 1 in m
l₂ = 2.0    # length of pendulum 2 in m
m₁ = 1.0    # mass of pendulum 1 in kg
m₂ = 1.0    # mass of pendulum 2 in kg

function deriv!(du, u, p, t)
    θ₁, ω₁, θ₂, ω₂ = u

    # From solving the lagrangian:
    Δθ = θ₂ - θ₁
    f₁ = (l₂ / l₁) * (m₂ / (m₁ + m₂)) * ω₂ * ω₂ * sin(Δθ) - (g / l₁) * sin(θ₁)
    f₂ = -(l₁ / l₂) * ω₁ * ω₁ * sin(Δθ) - (g / l₂) * sin(θ₂)
    α₁ = (l₂ / l₁) * (m₂ / (m₁ + m₂)) * cos(Δθ)
    α₂ = (l₁ / l₂) * cos(Δθ)
    
    du[1] = ω₁
    du[2] = (f₁ - α₁ * f₂) / (1 - α₁ * α₂)
    du[3] = ω₂
    du[4] = (-α₂ * f₁ + f₂) / (1 - α₁ * α₂)
end

# Initial conditions
θ₁ = 120.0
ω₁ = 0.0
θ₂ = -10.0
ω₂ = 0.0

t_range = (0.0, 100.0)
u = deg2rad.([θ₁, ω₁, θ₂, ω₂])

# Setup and solve ODE system
prob = ODEProblem(deriv!, u, t_range)
sol = solve(prob, reltol=1e-6)

# Create a gif
println("Starting creation.")
anim = @animate for t in tqdm(0:0.025:50)
    frame = sol(t)
    x₁ = l₁ * sin(frame[1])
    y₁ = -l₁ * cos(frame[1])
    x₂ = l₂ * sin(frame[3]) + x₁
    y₂ = -l₂ * cos(frame[3]) + y₁

    plot(
        [0, x₁, x₂], [0, y₁, y₂], 
        xlim = (-3, 3), 
        ylim = (-3, 3),
        label = @sprintf("t = %.2f", t),
        aspect_ratio = :equal,
        market = :circle,
        size = (500, 500)
    )

    scatter!([0, x₁, x₂], [0, y₁, y₂], label="")
end

gif(anim, "output.gif", fps=30)
println("Finished creation.")
