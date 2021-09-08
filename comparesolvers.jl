using DifferentialEquations,OrdinaryDiffEq, DiffEqDevTools, Sundials, ParameterizedFunctions,Plots, LinearAlgebra

rober = @ode_def begin
  dy₁ = -k₁*y₁+k₃*y₂*y₃
  dy₂ =  k₁*y₁-k₂*y₂^2-k₃*y₂*y₃
  dy₃ =  k₂*y₂^2
end k₁ k₂ k₃

prob = ODEProblem(rober,[1.0,0.0,0.0],(0.0,1e5),(1e-2,2e6,3e8)) # function, y0, tspan, k0
sol = solve(prob,CVODE_BDF(),abstol=1/10^14,reltol=1/10^14)
sol_1 = solve(prob,AutoVern7(Rodas5(),nonstifftol = 11/10),abstol=1/10^10,reltol=1/10^6)
sol_2 = solve(prob,AutoVern7(Rodas5(),nonstifftol = 11/10))
sol_3 = solve(prob,AutoTsit5(Rosenbrock23(autodiff=false),nonstifftol=9/10),abstol=1/10^6,reltol=1/10^3)
test_sol = TestSolution(sol)

abstols = 1.0 / 10.0^10
reltols = 1.0 / 10.0^6
wp_1 = WorkPrecision(prob,AutoVern7(Rodas5(),nonstifftol = 11/10),abstols,reltols;name="AutoVern7(Rodas5(),nonstifftol = 11/10)",
                      appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
abstols = 1.0 / 10.0^6
reltols = 1.0 / 10.0^3
wp_2 = WorkPrecision(prob,AutoTsit5(Rosenbrock23(autodiff=false),nonstifftol=9/10),abstols,reltols;name="AutoTsit5(Rosenbrock23(autodiff=false),nonstifftol=9/10)",
                      appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
#wp=WorkPrecision(prob,alg,abstols,reltols,dts=nothing;
#                       name=nothing,appxsol=nothing,error_estimate=:final,numruns=20,seconds=2,kwargs...)
