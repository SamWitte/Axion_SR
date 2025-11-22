using SciMLSensitivity
using DifferentialEquations, Flux, Random, Plots
using IterTools: ncycle

rng = Random.default_rng()
layersN = 3
nodesN = 20

ann = Chain(Dense(layersN, nodesN, tanh), Dense(nodesN, layersN, tanh))
θ, re = Flux.destructure(ann)

function dudt_(u, p, t)
    re(p)(u)[1] .* u
end

function predict_adjoint(time_batch)
    _prob = remake(prob, u0 = u0, p = θ)
    Array(solve(_prob, Tsit5(), saveat = time_batch))
end

function loss_adjoint(batch, time_batch)
    pred = predict_adjoint(time_batch)
    sum(abs2, batch - pred)#, pred
end

u0 = Float32[200.0]
datasize = 30
tspan = (0.0f0, 3.0f0)

t = range(tspan[1], tspan[2], length = datasize)
true_prob = ODEProblem(true_sol, u0, tspan)
ode_data = Array(solve(true_prob, Tsit5(), saveat = t))

prob = ODEProblem{false}(dudt_, u0, tspan, θ)

k = 10
train_loader = Flux.Data.DataLoader((ode_data, t), batchsize = k)

for (x, y) in train_loader
    @show x
    @show y
end

numEpochs = 300
losses = []
function cb()
    begin
        l = loss_adjoint(ode_data, t)
        push!(losses, l)
        @show l
        pred = predict_adjoint(t)
        pl = scatter(t, ode_data[1, :], label = "data", color = :black, ylim = (150, 200))
        scatter!(pl, t, pred[1, :], label = "prediction", color = :darkgreen)
        display(plot(pl))
        false
    end
end

opt = Adam(0.05)
Flux.train!(loss_adjoint, Flux.params(θ), ncycle(train_loader, numEpochs), opt,
    cb = Flux.throttle(cb, 10))

#Now lets see how well it generalizes to new initial conditions

starting_temp = collect(10:30:250)
true_prob_func(u0) = ODEProblem(true_sol, [u0], tspan)
color_cycle = palette(:tab10)
pl = plot()
for (j, temp) in enumerate(starting_temp)
    ode_test_sol = solve(ODEProblem(true_sol, [temp], (0.0f0, 10.0f0)), Tsit5(),
        saveat = 0.0:0.5:10.0)
    ode_nn_sol = solve(ODEProblem{false}(dudt_, [temp], (0.0f0, 10.0f0), θ))
    scatter!(pl, ode_test_sol, var = (0, 1), label = "", color = color_cycle[j])
    plot!(pl, ode_nn_sol, var = (0, 1), label = "", color = color_cycle[j], lw = 2.0)
end
display(pl)
title!("Neural ODE for Newton's Law of Cooling: Test Data")
xlabel!("Time")
ylabel!("Temp")
