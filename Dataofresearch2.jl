using OrdinaryDiffEq, DiffEqDevTools, Sundials, ParameterizedFunctions,Plots, LinearAlgebra

rober = @ode_def begin
  dy₁ = -k₁*y₁+k₃*y₂*y₃
  dy₂ =  k₁*y₁-k₂*y₂^2-k₃*y₂*y₃
  dy₃ =  k₂*y₂^2
end k₁ k₂ k₃

prob = ODEProblem(rober,[1.0,0.0,0.0],(0.0,1e5),(1e6,0.003,4e3)) # function, y0, tspan, k0
sol = solve(prob,CVODE_BDF(),abstol=1/10^14,reltol=1/10^14)
test_sol = TestSolution(sol)
abstols = 1.0 ./ 10.0 .^ (5:10)
reltols = 1.0 ./ 10.0 .^ (1:6)
#plot(sol,labels=["y1","y2","y3"])
#wp=WorkPrecision(prob,alg,abstols,reltols,dts=nothing;
#                       name=nothing,appxsol=nothing,error_estimate=:final,numruns=20,seconds=2,kwargs...)
#AutoTsit5
wp_T5KC3 = WorkPrecision(prob,AutoTsit5(KenCarp3(),nonstifftol = 14/10),abstols,reltols;name="AutoTsit5(KenCarp3(),nonstifftol = 14/10))",
                      appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_T5KC3_1 = WorkPrecision(prob,AutoTsit5(KenCarp3(),nonstifftol = 11/10),abstols,reltols;name="AutoTsit5(KenCarp3(),nonstifftol = 11/10))",
                      appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_T5KC3_11 = WorkPrecision(prob,AutoTsit5(KenCarp3(),nonstifftol = 11/10),abstols,reltols;name="AutoTsit5(KenCarp3(),nonstifftol = 11/10))",
                      appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_T5KC3_2= WorkPrecision(prob,AutoTsit5(KenCarp3(),nonstifftol = 9/10),abstols,reltols;name="AutoTsit5(KenCarp3(),nonstifftol = 9/10))",
                      appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_T5KC3_3 = WorkPrecision(prob,AutoTsit5(KenCarp3(),nonstifftol = 6/10),abstols,reltols;name="AutoTsit5(KenCarp3(),nonstifftol = 6/10))",
                      appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_T5KC3_4= WorkPrecision(prob,AutoTsit5(KenCarp3(),nonstifftol = 3/10),abstols,reltols;name="AutoTsit5(KenCarp3(),nonstifftol = 3/10))",
                      appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)

wp_T5KC4 = WorkPrecision(prob,AutoTsit5(KenCarp4(),nonstifftol = 14/10),abstols,reltols;name="AutoTsit5(KenCarp4(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_T5KC4_1 = WorkPrecision(prob,AutoTsit5(KenCarp4(),nonstifftol = 11/10),abstols,reltols;name="AutoTsit5(KenCarp4(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_T5KC4_11 = WorkPrecision(prob,AutoTsit5(KenCarp4(),nonstifftol = 11/10),abstols,reltols;name="AutoTsit5(KenCarp4(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_T5KC4_2= WorkPrecision(prob,AutoTsit5(KenCarp4(),nonstifftol = 9/10),abstols,reltols;name="AutoTsit5(KenCarp4(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_T5KC4_3 = WorkPrecision(prob,AutoTsit5(KenCarp4(),nonstifftol = 6/10),abstols,reltols;name="AutoTsit5(KenCarp4(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_T5KC4_4= WorkPrecision(prob,AutoTsit5(KenCarp4(),nonstifftol = 3/10),abstols,reltols;name="AutoTsit5(KenCarp4(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)

wp_T5KC5 = WorkPrecision(prob,AutoTsit5(KenCarp5(),nonstifftol = 14/10),abstols,reltols;name="AutoTsit5(KenCarp5(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_T5KC5_1 = WorkPrecision(prob,AutoTsit5(KenCarp5(),nonstifftol = 11/10),abstols,reltols;name="AutoTsit5(KenCarp5(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_T5KC5_11 = WorkPrecision(prob,AutoTsit5(KenCarp5(),nonstifftol = 11/10),abstols,reltols;name="AutoTsit5(KenCarp5(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_T5KC5_2= WorkPrecision(prob,AutoTsit5(KenCarp5(),nonstifftol = 9/10),abstols,reltols;name="AutoTsit5(KenCarp5(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_T5KC5_3 = WorkPrecision(prob,AutoTsit5(KenCarp5(),nonstifftol = 6/10),abstols,reltols;name="AutoTsit5(KenCarp5(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_T5KC5_4= WorkPrecision(prob,AutoTsit5(KenCarp5(),nonstifftol = 3/10),abstols,reltols;name="AutoTsit5(KenCarp5(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)


wp_T5Tz = WorkPrecision(prob,AutoTsit5(Trapezoid(),nonstifftol = 14/10),abstols,reltols;name="AutoTsit5(Trapezoid(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_T5Tz_1 = WorkPrecision(prob,AutoTsit5(Trapezoid(),nonstifftol = 11/10),abstols,reltols;name="AutoTsit5(Trapezoid(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_T5Tz_11 = WorkPrecision(prob,AutoTsit5(Trapezoid(),nonstifftol = 11/10),abstols,reltols;name="AutoTsit5(Trapezoid(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_T5Tz_2= WorkPrecision(prob,AutoTsit5(Trapezoid(),nonstifftol = 9/10),abstols,reltols;name="AutoTsit5(Trapezoid(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_T5Tz_3 = WorkPrecision(prob,AutoTsit5(Trapezoid(),nonstifftol = 6/10),abstols,reltols;name="AutoTsit5(Trapezoid(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_T5Tz_4= WorkPrecision(prob,AutoTsit5(Trapezoid(),nonstifftol = 3/10),abstols,reltols;name="AutoTsit5(Trapezoid(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)


wp_T5Rb23 = WorkPrecision(prob,AutoTsit5(Rosenbrock23(),nonstifftol = 14/10),abstols,reltols;name="AutoTsit5(Rosenbrock23(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_T5Rb23_1 = WorkPrecision(prob,AutoTsit5(Rosenbrock23(),nonstifftol = 11/10),abstols,reltols;name="AutoTsit5(Rosenbrock23(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_T5Rb23_11 = WorkPrecision(prob,AutoTsit5(Rosenbrock23(),nonstifftol = 11/10),abstols,reltols;name="AutoTsit5(Rosenbrock23(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_T5Rb23_2= WorkPrecision(prob,AutoTsit5(Rosenbrock23(),nonstifftol = 9/10),abstols,reltols;name="AutoTsit5(Rosenbrock23(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_T5Rb23_3 = WorkPrecision(prob,AutoTsit5(Rosenbrock23(),nonstifftol = 5/10),abstols,reltols;name="AutoTsit5(Rosenbrock23(),nonstifftol = 5/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
# nonstifftol=6/10不出结果
wp_T5Rb23_4= WorkPrecision(prob,AutoTsit5(Rosenbrock23(),nonstifftol = 3/10),abstols,reltols;name="AutoTsit5(Rosenbrock23(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)


wp_T5Rd3 = WorkPrecision(prob,AutoTsit5(Rodas3(),nonstifftol = 14/10),abstols,reltols;name="AutoTsit5(Rodas3(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_T5Rd3_1 = WorkPrecision(prob,AutoTsit5(Rodas3(),nonstifftol = 11/10),abstols,reltols;name="AutoTsit5(Rodas3(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_T5Rd3_11 = WorkPrecision(prob,AutoTsit5(Rodas3(),nonstifftol = 11/10),abstols,reltols;name="AutoTsit5(Rodas3(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_T5Rd3_2 = WorkPrecision(prob,AutoTsit5(Rodas3(),nonstifftol = 9/10),abstols,reltols;name="AutoTsit5(Rodas3(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_T5Rd3_3 = WorkPrecision(prob,AutoTsit5(Rodas3(),nonstifftol = 6/10),abstols,reltols;name="AutoTsit5(Rodas3(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_T5Rd3_4 = WorkPrecision(prob,AutoTsit5(Rodas3(),nonstifftol = 3/10),abstols,reltols;name="AutoTsit5(Rodas3(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)


wp_T5Rd4 = WorkPrecision(prob,AutoTsit5(Rodas4(),nonstifftol = 14/10),abstols,reltols;name="AutoTsit5(Rodas4(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_T5Rd4_1 = WorkPrecision(prob,AutoTsit5(Rodas4(),nonstifftol = 11/10),abstols,reltols;name="AutoTsit5(Rodas4(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_T5Rd4_11 = WorkPrecision(prob,AutoTsit5(Rodas4(),nonstifftol = 11/10),abstols,reltols;name="AutoTsit5(Rodas4(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_T5Rd4_2= WorkPrecision(prob,AutoTsit5(Rodas4(),nonstifftol = 9/10),abstols,reltols;name="AutoTsit5(Rodas4(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_T5Rd4_3 = WorkPrecision(prob,AutoTsit5(Rodas4(),nonstifftol = 6/10),abstols,reltols;name="AutoTsit5(Rodas4(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_T5Rd4_4 = WorkPrecision(prob,AutoTsit5(Rodas4(),nonstifftol = 4/10),abstols,reltols;name="AutoTsit5(Rodas4(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)


wp_T5Rd5 = WorkPrecision(prob,AutoTsit5(Rodas5(),nonstifftol = 14/10),abstols,reltols;name="AutoTsit5(Rodas5(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_T5Rd5_1 = WorkPrecision(prob,AutoTsit5(Rodas5(),nonstifftol = 11/10),abstols,reltols;name="AutoTsit5(Rodas5(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_T5Rd5_11 = WorkPrecision(prob,AutoTsit5(Rodas5(),nonstifftol = 11/10),abstols,reltols;name="AutoTsit5(Rodas5(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_T5Rd5_2= WorkPrecision(prob,AutoTsit5(Rodas5(),nonstifftol = 9/10),abstols,reltols;name="AutoTsit5(Rodas5(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_T5Rd5_3 = WorkPrecision(prob,AutoTsit5(Rodas5(),nonstifftol = 6/10),abstols,reltols;name="AutoTsit5(Rodas5(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_T5Rd5_4= WorkPrecision(prob,AutoTsit5(Rodas5(),nonstifftol = 3/10),abstols,reltols;name="AutoTsit5(Rodas5(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)

# AutoDP5
wp_DP5KC3 = WorkPrecision(prob,AutoDP5(KenCarp3(),nonstifftol = 14/10),abstols,reltols;name="AutoDP5(KenCarp3(),nonstifftol = 14/10))",
                      appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_DP5KC3_1 = WorkPrecision(prob,AutoDP5(KenCarp3(),nonstifftol = 11/10),abstols,reltols;name="AutoDP5(KenCarp3(),nonstifftol = 11/10))",
                      appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_DP5KC3_11 = WorkPrecision(prob,AutoDP5(KenCarp3(),nonstifftol = 11/10),abstols,reltols;name="AutoDP5(KenCarp3(),nonstifftol = 11/10))",
                      appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_DP5KC3_2= WorkPrecision(prob,AutoDP5(KenCarp3(),nonstifftol = 9/10),abstols,reltols;name="AutoDP5(KenCarp3(),nonstifftol = 9/10))",
                      appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_DP5KC3_3 = WorkPrecision(prob,AutoDP5(KenCarp3(),nonstifftol = 6/10),abstols,reltols;name="AutoDP5(KenCarp3(),nonstifftol = 6/10))",
                      appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_DP5KC3_4= WorkPrecision(prob,AutoDP5(KenCarp3(),nonstifftol = 3/10),abstols,reltols;name="AutoDP5(KenCarp3(),nonstifftol = 3/10))",
                      appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)

wp_DP5KC4 = WorkPrecision(prob,AutoDP5(KenCarp4(),nonstifftol = 14/10),abstols,reltols;name="AutoDP5(KenCarp4(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_DP5KC4_1 = WorkPrecision(prob,AutoDP5(KenCarp4(),nonstifftol = 11/10),abstols,reltols;name="AutoDP5(KenCarp4(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_DP5KC4_11 = WorkPrecision(prob,AutoDP5(KenCarp4(),nonstifftol = 11/10),abstols,reltols;name="AutoDP5(KenCarp4(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_DP5KC4_2= WorkPrecision(prob,AutoDP5(KenCarp4(),nonstifftol = 9/10),abstols,reltols;name="AutoDP5(KenCarp4(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_DP5KC4_3 = WorkPrecision(prob,AutoDP5(KenCarp4(),nonstifftol = 6/10),abstols,reltols;name="AutoDP5(KenCarp4(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_DP5KC4_4= WorkPrecision(prob,AutoDP5(KenCarp4(),nonstifftol = 3/10),abstols,reltols;name="AutoDP5(KenCarp4(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)

wp_DP5KC5 = WorkPrecision(prob,AutoDP5(KenCarp5(),nonstifftol = 14/10),abstols,reltols;name="AutoDP5(KenCarp5(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_DP5KC5_1 = WorkPrecision(prob,AutoDP5(KenCarp5(),nonstifftol = 11/10),abstols,reltols;name="AutoDP5(KenCarp5(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_DP5KC5_11 = WorkPrecision(prob,AutoDP5(KenCarp5(),nonstifftol = 11/10),abstols,reltols;name="AutoDP5(KenCarp5(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_DP5KC5_2= WorkPrecision(prob,AutoDP5(KenCarp5(),nonstifftol = 9/10),abstols,reltols;name="AutoDP5(KenCarp5(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_DP5KC5_3 = WorkPrecision(prob,AutoDP5(KenCarp5(),nonstifftol = 6/10),abstols,reltols;name="AutoDP5(KenCarp5(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_DP5KC5_4= WorkPrecision(prob,AutoDP5(KenCarp5(),nonstifftol = 3/10),abstols,reltols;name="AutoDP5(KenCarp5(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)


wp_DP5Tz = WorkPrecision(prob,AutoDP5(Trapezoid(),nonstifftol = 14/10),abstols,reltols;name="AutoDP5(Trapezoid(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_DP5Tz_1 = WorkPrecision(prob,AutoDP5(Trapezoid(),nonstifftol = 11/10),abstols,reltols;name="AutoDP5(Trapezoid(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_DP5Tz_11 = WorkPrecision(prob,AutoDP5(Trapezoid(),nonstifftol = 11/10),abstols,reltols;name="AutoDP5(Trapezoid(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_DP5Tz_2= WorkPrecision(prob,AutoDP5(Trapezoid(),nonstifftol = 9/10),abstols,reltols;name="AutoDP5(Trapezoid(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_DP5Tz_3 = WorkPrecision(prob,AutoDP5(Trapezoid(),nonstifftol = 6/10),abstols,reltols;name="AutoDP5(Trapezoid(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_DP5Tz_4 = WorkPrecision(prob,AutoDP5(Trapezoid(),nonstifftol = 3/10),abstols,reltols;name="AutoDP5(Trapezoid(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)


wp_DP5Rb23 = WorkPrecision(prob,AutoDP5(Rosenbrock23(),nonstifftol = 14/10),abstols,reltols;name="AutoDP5(Rosenbrock23(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_DP5Rb23_1 = WorkPrecision(prob,AutoDP5(Rosenbrock23(),nonstifftol = 11/10),abstols,reltols;name="AutoDP5(Rosenbrock23(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_DP5Rb23_11 = WorkPrecision(prob,AutoDP5(Rosenbrock23(),nonstifftol = 11/10),abstols,reltols;name="AutoDP5(Rosenbrock23(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_DP5Rb23_2= WorkPrecision(prob,AutoDP5(Rosenbrock23(),nonstifftol = 9/10),abstols,reltols;name="AutoDP5(Rosenbrock23(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_DP5Rb23_3 = WorkPrecision(prob,AutoDP5(Rosenbrock23(),nonstifftol = 6/10),abstols,reltols;name="AutoDP5(Rosenbrock23(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_DP5Rb23_4= WorkPrecision(prob,AutoDP5(Rosenbrock23(),nonstifftol = 3/10),abstols,reltols;name="AutoDP5(Rosenbrock23(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)


wp_DP5Rd3 = WorkPrecision(prob,AutoDP5(Rodas3(),nonstifftol = 14/10),abstols,reltols;name="AutoDP5(Rodas3(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_DP5Rd3_1 = WorkPrecision(prob,AutoDP5(Rodas3(),nonstifftol = 11/10),abstols,reltols;name="AutoDP5(Rodas3(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_DP5Rd3_11 = WorkPrecision(prob,AutoDP5(Rodas3(),nonstifftol = 11/10),abstols,reltols;name="AutoDP5(Rodas3(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_DP5Rd3_2= WorkPrecision(prob,AutoDP5(Rodas3(),nonstifftol = 9/10),abstols,reltols;name="AutoDP5(Rodas3(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_DP5Rd3_3 = WorkPrecision(prob,AutoDP5(Rodas3(),nonstifftol = 6/10),abstols,reltols;name="AutoDP5(Rodas3(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_DP5Rd3_4= WorkPrecision(prob,AutoDP5(Rodas3(),nonstifftol = 3/10),abstols,reltols;name="AutoDP5(Rodas3(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)


wp_DP5Rd4 = WorkPrecision(prob,AutoDP5(Rodas4(),nonstifftol = 14/10),abstols,reltols;name="AutoTsit5(Rodas4(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_DP5Rd4_1 = WorkPrecision(prob,AutoDP5(Rodas4(),nonstifftol = 11/10),abstols,reltols;name="AutoTsit5(Rodas4(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_DP5Rd4_11 = WorkPrecision(prob,AutoDP5(Rodas4(),nonstifftol = 11/10),abstols,reltols;name="AutoTsit5(Rodas4(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_DP5Rd4_2= WorkPrecision(prob,AutoDP5(Rodas4(),nonstifftol = 9/10),abstols,reltols;name="AutoTsit5(Rodas4(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_DP5Rd4_3 = WorkPrecision(prob,AutoDP5(Rodas4(),nonstifftol = 6/10),abstols,reltols;name="AutoTsit5(Rodas4(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_DP5Rd4_4= WorkPrecision(prob,AutoDP5(Rodas4(),nonstifftol = 3/10),abstols,reltols;name="AutoTsit5(Rodas4(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)


wp_DP5Rd5 = WorkPrecision(prob,AutoDP5(Rodas5(),nonstifftol = 14/10),abstols,reltols;name="AutoDP5(Rodas5(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_DP5Rd5_1 = WorkPrecision(prob,AutoDP5(Rodas5(),nonstifftol = 11/10),abstols,reltols;name="AutoDP5(Rodas5(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_DP5Rd5_11 = WorkPrecision(prob,AutoDP5(Rodas5(),nonstifftol = 11/10),abstols,reltols;name="AutoDP5(Rodas5(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_DP5Rd5_2= WorkPrecision(prob,AutoDP5(Rodas5(),nonstifftol = 9/10),abstols,reltols;name="AutoDP5(Rodas5(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_DP5Rd5_3 = WorkPrecision(prob,AutoDP5(Rodas5(),nonstifftol = 6/10),abstols,reltols;name="AutoDP5(Rodas5(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_DP5Rd5_4= WorkPrecision(prob,AutoDP5(Rodas5(),nonstifftol = 3/10),abstols,reltols;name="AutoDP5(Rodas5(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)

# AutoVern6
wp_V6KC3 = WorkPrecision(prob,AutoVern6(KenCarp3(),nonstifftol = 14/10),abstols,reltols;name="AutoVern6(KenCarp3(),nonstifftol = 14/10))",
                      appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V6KC3_1 = WorkPrecision(prob,AutoVern6(KenCarp3(),nonstifftol = 11/10),abstols,reltols;name="AutoVern6(KenCarp3(),nonstifftol = 11/10))",
                      appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V6KC3_11 = WorkPrecision(prob,AutoVern6(KenCarp3(),nonstifftol = 11/10),abstols,reltols;name="AutoVern6(KenCarp3(),nonstifftol = 11/10))",
                      appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_V6KC3_2= WorkPrecision(prob,AutoVern6(KenCarp3(),nonstifftol = 9/10),abstols,reltols;name="AutoVern6(KenCarp3(),nonstifftol = 9/10))",
                      appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V6KC3_3 = WorkPrecision(prob,AutoVern6(KenCarp3(),nonstifftol = 6/10),abstols,reltols;name="AutoVern6(KenCarp3(),nonstifftol = 6/10))",
                      appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V6KC3_4= WorkPrecision(prob,AutoVern6(KenCarp3(),nonstifftol = 3/10),abstols,reltols;name="AutoVern6(KenCarp3(),nonstifftol = 3/10))",
                      appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)

wp_V6KC4 = WorkPrecision(prob,AutoVern6(KenCarp4(),nonstifftol = 14/10),abstols,reltols;name="AutoVern6(KenCarp4(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V6KC4_1 = WorkPrecision(prob,AutoVern6(KenCarp4(),nonstifftol = 11/10),abstols,reltols;name="AutoVern6(KenCarp4(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V6KC4_11 = WorkPrecision(prob,AutoVern6(KenCarp4(),nonstifftol = 11/10),abstols,reltols;name="AutoVern6(KenCarp4(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_V6KC4_2= WorkPrecision(prob,AutoVern6(KenCarp4(),nonstifftol = 9/10),abstols,reltols;name="AutoVern6(KenCarp4(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V6KC4_3 = WorkPrecision(prob,AutoVern6(KenCarp4(),nonstifftol = 6/10),abstols,reltols;name="AutoVern6(KenCarp4(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V6KC4_4= WorkPrecision(prob,AutoVern6(KenCarp4(),nonstifftol = 3/10),abstols,reltols;name="AutoVern6(KenCarp4(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)

wp_V6KC5 = WorkPrecision(prob,AutoVern6(KenCarp5(),nonstifftol = 14/10),abstols,reltols;name="AutoVern6(KenCarp5(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V6KC5_1 = WorkPrecision(prob,AutoVern6(KenCarp5(),nonstifftol = 11/10),abstols,reltols;name="AutoVern6(KenCarp5(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V6KC5_11 = WorkPrecision(prob,AutoVern6(KenCarp5(),nonstifftol = 11/10),abstols,reltols;name="AutoVern6(KenCarp5(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_V6KC5_2= WorkPrecision(prob,AutoVern6(KenCarp5(),nonstifftol = 9/10),abstols,reltols;name="AutoVern6(KenCarp5(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V6KC5_3 = WorkPrecision(prob,AutoVern6(KenCarp5(),nonstifftol = 6/10),abstols,reltols;name="AutoVern6(KenCarp5(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V6KC5_4= WorkPrecision(prob,AutoVern6(KenCarp5(),nonstifftol = 3/10),abstols,reltols;name="AutoVern6(KenCarp5(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)


wp_V6Tz = WorkPrecision(prob,AutoVern6(Trapezoid(),nonstifftol = 14/10),abstols,reltols;name="AutoVern6(Trapezoid(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V6Tz_1 = WorkPrecision(prob,AutoVern6(Trapezoid(),nonstifftol = 11/10),abstols,reltols;name="AutoVern6(Trapezoid(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V6Tz_11 = WorkPrecision(prob,AutoVern6(Trapezoid(),nonstifftol = 11/10),abstols,reltols;name="AutoVern6(Trapezoid(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_V6Tz_2= WorkPrecision(prob,AutoVern6(Trapezoid(),nonstifftol = 9/10),abstols,reltols;name="AutoVern6(Trapezoid(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V6Tz_3 = WorkPrecision(prob,AutoVern6(Trapezoid(),nonstifftol = 6/10),abstols,reltols;name="AutoVern6(Trapezoid(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V6Tz_4= WorkPrecision(prob,AutoVern6(Trapezoid(),nonstifftol = 3/10),abstols,reltols;name="AutoVern6(Trapezoid(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)


wp_V6Rb23 = WorkPrecision(prob,AutoVern6(Rosenbrock23(),nonstifftol = 14/10),abstols,reltols;name="AutoVern6(Rosenbrock23(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V6Rb23_1 = WorkPrecision(prob,AutoVern6(Rosenbrock23(),nonstifftol = 11/10),abstols,reltols;name="AutoVern6(Rosenbrock23(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V6Rb23_11 = WorkPrecision(prob,AutoVern6(Rosenbrock23(),nonstifftol = 11/10),abstols,reltols;name="AutoVern6(Rosenbrock23(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_V6Rb23_2= WorkPrecision(prob,AutoVern6(Rosenbrock23(),nonstifftol = 9/10),abstols,reltols;name="AutoVern6(Rosenbrock23(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V6Rb23_3 = WorkPrecision(prob,AutoVern6(Rosenbrock23(),nonstifftol = 6/10),abstols,reltols;name="AutoVern6(Rosenbrock23(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V6Rb23_4= WorkPrecision(prob,AutoVern6(Rosenbrock23(),nonstifftol = 3/10),abstols,reltols;name="AutoVern6(Rosenbrock23(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)


wp_V6Rd3 = WorkPrecision(prob,AutoVern6(Rodas3(),nonstifftol = 14/10),abstols,reltols;name="AutoVern6(Rodas3(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V6Rd3_1 = WorkPrecision(prob,AutoVern6(Rodas3(),nonstifftol = 11/10),abstols,reltols;name="AutoVern6(Rodas3(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V6Rd3_11 = WorkPrecision(prob,AutoVern6(Rodas3(),nonstifftol = 11/10),abstols,reltols;name="AutoVern6(Rodas3(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_V6Rd3_2= WorkPrecision(prob,AutoVern6(Rodas3(),nonstifftol = 9/10),abstols,reltols;name="AutoVern6(Rodas3(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V6Rd3_3 = WorkPrecision(prob,AutoVern6(Rodas3(),nonstifftol = 6/10),abstols,reltols;name="AutoVern6(Rodas3(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V6Rd3_4= WorkPrecision(prob,AutoVern6(Rodas3(),nonstifftol = 3/10),abstols,reltols;name="AutoVern6(Rodas3(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)


wp_V6Rd4 = WorkPrecision(prob,AutoVern6(Rodas4(),nonstifftol = 14/10),abstols,reltols;name="AutoVern6(Rodas4(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V6Rd4_1 = WorkPrecision(prob,AutoVern6(Rodas4(),nonstifftol = 11/10),abstols,reltols;name="AutoVern6(Rodas4(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V6Rd4_11 = WorkPrecision(prob,AutoVern6(Rodas4(),nonstifftol = 11/10),abstols,reltols;name="AutoVern6(Rodas4(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_V6Rd4_2= WorkPrecision(prob,AutoVern6(Rodas4(),nonstifftol = 9/10),abstols,reltols;name="AutoVern6(Rodas4(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V6Rd4_3 = WorkPrecision(prob,AutoVern6(Rodas4(),nonstifftol = 6/10),abstols,reltols;name="AutoVern6(Rodas4(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V6Rd4_4= WorkPrecision(prob,AutoVern6(Rodas4(),nonstifftol = 3/10),abstols,reltols;name="AutoVern6(Rodas4(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)


wp_V6Rd5 = WorkPrecision(prob,AutoVern6(Rodas5(),nonstifftol = 14/10),abstols,reltols;name="AutoVern6(Rodas5(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V6Rd5_1 = WorkPrecision(prob,AutoVern6(Rodas5(),nonstifftol = 11/10),abstols,reltols;name="AutoVern6(Rodas5(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V6Rd5_11 = WorkPrecision(prob,AutoVern6(Rodas5(),nonstifftol = 11/10),abstols,reltols;name="AutoVern6(Rodas5(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_V6Rd5_2= WorkPrecision(prob,AutoVern6(Rodas5(),nonstifftol = 9/10),abstols,reltols;name="AutoVern6(Rodas5(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V6Rd5_3 = WorkPrecision(prob,AutoVern6(Rodas5(),nonstifftol = 6/10),abstols,reltols;name="AutoVern6(Rodas5(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V6Rd5_4= WorkPrecision(prob,AutoVern6(Rodas5(),nonstifftol = 3/10),abstols,reltols;name="AutoVern6(Rodas5(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)

# AutoVern7
wp_V7KC3 = WorkPrecision(prob,AutoVern7(KenCarp3(),nonstifftol = 14/10),abstols,reltols;name="AutoVern7(KenCarp3(),nonstifftol = 14/10))",
                      appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V7KC3_1 = WorkPrecision(prob,AutoVern7(KenCarp3(),nonstifftol = 11/10),abstols,reltols;name="AutoVern7(KenCarp3(),nonstifftol = 11/10))",
                      appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V7KC3_11 = WorkPrecision(prob,AutoVern7(KenCarp3(),nonstifftol = 11/10),abstols,reltols;name="AutoVern7(KenCarp3(),nonstifftol = 11/10))",
                      appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_V7KC3_2= WorkPrecision(prob,AutoVern7(KenCarp3(),nonstifftol = 9/10),abstols,reltols;name="AutoVern7(KenCarp3(),nonstifftol = 9/10))",
                      appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V7KC3_3 = WorkPrecision(prob,AutoVern7(KenCarp3(),nonstifftol = 6/10),abstols,reltols;name="AutoVern7(KenCarp3(),nonstifftol = 6/10))",
                      appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V7KC3_4= WorkPrecision(prob,AutoVern7(KenCarp3(),nonstifftol = 3/10),abstols,reltols;name="AutoVern7(KenCarp3(),nonstifftol = 3/10))",
                      appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)

wp_V7KC4 = WorkPrecision(prob,AutoVern7(KenCarp4(),nonstifftol = 14/10),abstols,reltols;name="AutoVern7(KenCarp4(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V7KC4_1 = WorkPrecision(prob,AutoVern7(KenCarp4(),nonstifftol = 11/10),abstols,reltols;name="AutoVern7(KenCarp4(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V7KC4_11 = WorkPrecision(prob,AutoVern7(KenCarp4(),nonstifftol = 11/10),abstols,reltols;name="AutoVern7(KenCarp4(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_V7KC4_2= WorkPrecision(prob,AutoVern7(KenCarp4(),nonstifftol = 9/10),abstols,reltols;name="AutoVern7(KenCarp4(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V7KC4_3 = WorkPrecision(prob,AutoVern7(KenCarp4(),nonstifftol = 6/10),abstols,reltols;name="AutoVern7(KenCarp4(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V7KC4_4= WorkPrecision(prob,AutoVern7(KenCarp4(),nonstifftol = 3/10),abstols,reltols;name="AutoVern7(KenCarp4(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)

wp_V7KC5 = WorkPrecision(prob,AutoVern7(KenCarp5(),nonstifftol = 14/10),abstols,reltols;name="AutoVern7(KenCarp5(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V7KC5_1 = WorkPrecision(prob,AutoVern7(KenCarp5(),nonstifftol = 11/10),abstols,reltols;name="AutoVern7(KenCarp5(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V7KC5_11 = WorkPrecision(prob,AutoVern7(KenCarp5(),nonstifftol = 11/10),abstols,reltols;name="AutoVern7(KenCarp5(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_V7KC5_2= WorkPrecision(prob,AutoVern7(KenCarp5(),nonstifftol = 9/10),abstols,reltols;name="AutoVern7(KenCarp5(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V7KC5_3 = WorkPrecision(prob,AutoVern7(KenCarp5(),nonstifftol = 6/10),abstols,reltols;name="AutoVern7(KenCarp5(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V7KC5_4= WorkPrecision(prob,AutoVern7(KenCarp5(),nonstifftol = 3/10),abstols,reltols;name="AutoVern7(KenCarp5(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)


wp_V7Tz = WorkPrecision(prob,AutoVern7(Trapezoid(),nonstifftol = 14/10),abstols,reltols;name="AutoVern7(Trapezoid(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V7Tz_1 = WorkPrecision(prob,AutoVern7(Trapezoid(),nonstifftol = 11/10),abstols,reltols;name="AutoVern7(Trapezoid(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V7Tz_11 = WorkPrecision(prob,AutoVern7(Trapezoid(),nonstifftol = 11/10),abstols,reltols;name="AutoVern7(Trapezoid(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_V7Tz_2= WorkPrecision(prob,AutoVern7(Trapezoid(),nonstifftol = 9/10),abstols,reltols;name="AutoVern7(Trapezoid(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V7Tz_3 = WorkPrecision(prob,AutoVern7(Trapezoid(),nonstifftol = 6/10),abstols,reltols;name="AutoVern7(Trapezoid(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V7Tz_4= WorkPrecision(prob,AutoVern7(Trapezoid(),nonstifftol = 3/10),abstols,reltols;name="AutoVern7(Trapezoid(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)


wp_V7Rb23 = WorkPrecision(prob,AutoVern7(Rosenbrock23(),nonstifftol = 14/10),abstols,reltols;name="AutoVern7(Rosenbrock23(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V7Rb23_1 = WorkPrecision(prob,AutoVern7(Rosenbrock23(),nonstifftol = 11/10),abstols,reltols;name="AutoVern7(Rosenbrock23(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V7Rb23_11 = WorkPrecision(prob,AutoVern7(Rosenbrock23(),nonstifftol = 11/10),abstols,reltols;name="AutoVern7(Rosenbrock23(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_V7Rb23_2= WorkPrecision(prob,AutoVern7(Rosenbrock23(),nonstifftol = 9/10),abstols,reltols;name="AutoVern7(Rosenbrock23(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V7Rb23_3 = WorkPrecision(prob,AutoVern7(Rosenbrock23(),nonstifftol = 6/10),abstols,reltols;name="AutoVern7(Rosenbrock23(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V7Rb23_4= WorkPrecision(prob,AutoVern7(Rosenbrock23(),nonstifftol = 3/10),abstols,reltols;name="AutoVern7(Rosenbrock23(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)


wp_V7Rd3 = WorkPrecision(prob,AutoVern7(Rodas3(),nonstifftol = 14/10),abstols,reltols;name="AutoVern7(Rodas3(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V7Rd3_1 = WorkPrecision(prob,AutoVern7(Rodas3(),nonstifftol = 11/10),abstols,reltols;name="AutoVern7(Rodas3(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V7Rd3_11 = WorkPrecision(prob,AutoVern7(Rodas3(),nonstifftol = 11/10),abstols,reltols;name="AutoVern7(Rodas3(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_V7Rd3_2= WorkPrecision(prob,AutoVern7(Rodas3(),nonstifftol = 9/10),abstols,reltols;name="AutoVern7(Rodas3(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V7Rd3_3 = WorkPrecision(prob,AutoVern7(Rodas3(),nonstifftol = 6/10),abstols,reltols;name="AutoVern7(Rodas3(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V7Rd3_4= WorkPrecision(prob,AutoVern7(Rodas3(),nonstifftol = 3/10),abstols,reltols;name="AutoVern7(Rodas3(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)


wp_V7Rd4 = WorkPrecision(prob,AutoVern7(Rodas4(),nonstifftol = 14/10),abstols,reltols;name="AutoVern7(Rodas4(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V7Rd4_1 = WorkPrecision(prob,AutoVern7(Rodas4(),nonstifftol = 11/10),abstols,reltols;name="AutoVern7(Rodas4(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V7Rd4_11 = WorkPrecision(prob,AutoVern7(Rodas4(),nonstifftol = 11/10),abstols,reltols;name="AutoVern7(Rodas4(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_V7Rd4_2= WorkPrecision(prob,AutoVern7(Rodas4(),nonstifftol = 9/10),abstols,reltols;name="AutoVern7(Rodas4(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V7Rd4_3 = WorkPrecision(prob,AutoVern7(Rodas4(),nonstifftol = 6/10),abstols,reltols;name="AutoVern7(Rodas4(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V7Rd4_4= WorkPrecision(prob,AutoVern7(Rodas4(),nonstifftol = 3/10),abstols,reltols;name="AutoVern7(Rodas4(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)


wp_V7Rd5 = WorkPrecision(prob,AutoVern7(Rodas5(),nonstifftol = 14/10),abstols,reltols;name="AutoVern7(Rodas5(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V7Rd5_1 = WorkPrecision(prob,AutoVern7(Rodas5(),nonstifftol = 11/10),abstols,reltols;name="AutoVern7(Rodas5(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V7Rd5_11 = WorkPrecision(prob,AutoVern7(Rodas5(),nonstifftol = 11/10),abstols,reltols;name="AutoVern7(Rodas5(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_V7Rd5_2= WorkPrecision(prob,AutoVern7(Rodas5(),nonstifftol = 9/10),abstols,reltols;name="AutoVern7(Rodas5(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V7Rd5_3 = WorkPrecision(prob,AutoVern7(Rodas5(),nonstifftol = 6/10),abstols,reltols;name="AutoVern7(Rodas5(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V7Rd5_4= WorkPrecision(prob,AutoVern7(Rodas5(),nonstifftol = 3/10),abstols,reltols;name="AutoVern7(Rodas5(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)


# AutoVern8
wp_V8KC3 = WorkPrecision(prob,AutoVern8(KenCarp3(),nonstifftol = 14/10),abstols,reltols;name="AutoVern8(KenCarp3(),nonstifftol = 14/10))",
                      appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V8KC3_1 = WorkPrecision(prob,AutoVern8(KenCarp3(),nonstifftol = 11/10),abstols,reltols;name="AutoVern8(KenCarp3(),nonstifftol = 11/10))",
                      appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V8KC3_11 = WorkPrecision(prob,AutoVern8(KenCarp3(),nonstifftol = 11/10),abstols,reltols;name="AutoVern8(KenCarp3(),nonstifftol = 11/10))",
                      appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_V8KC3_2= WorkPrecision(prob,AutoVern8(KenCarp3(),nonstifftol = 9/10),abstols,reltols;name="AutoVern8(KenCarp3(),nonstifftol = 9/10))",
                      appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V8KC3_3 = WorkPrecision(prob,AutoVern8(KenCarp3(),nonstifftol = 6/10),abstols,reltols;name="AutoVern8(KenCarp3(),nonstifftol = 6/10))",
                      appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V8KC3_4= WorkPrecision(prob,AutoVern8(KenCarp3(),nonstifftol = 3/10),abstols,reltols;name="AutoVern8(KenCarp3(),nonstifftol = 3/10))",
                      appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)

wp_V8KC4 = WorkPrecision(prob,AutoVern8(KenCarp4(),nonstifftol = 14/10),abstols,reltols;name="AutoVern8(KenCarp4(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V8KC4_1 = WorkPrecision(prob,AutoVern8(KenCarp4(),nonstifftol = 11/10),abstols,reltols;name="AutoVern8(KenCarp4(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V8KC4_11 = WorkPrecision(prob,AutoVern8(KenCarp4(),nonstifftol = 11/10),abstols,reltols;name="AutoVern8(KenCarp4(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_V8KC4_2= WorkPrecision(prob,AutoVern8(KenCarp4(),nonstifftol = 9/10),abstols,reltols;name="AutoVern8(KenCarp4(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V8KC4_3 = WorkPrecision(prob,AutoVern8(KenCarp4(),nonstifftol = 6/10),abstols,reltols;name="AutoVern8(KenCarp4(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V8KC4_4= WorkPrecision(prob,AutoVern8(KenCarp4(),nonstifftol = 3/10),abstols,reltols;name="AutoVern8(KenCarp4(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)

wp_V8KC5 = WorkPrecision(prob,AutoVern8(KenCarp5(),nonstifftol = 14/10),abstols,reltols;name="AutoVern8(KenCarp5(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V8KC5_1 = WorkPrecision(prob,AutoVern8(KenCarp5(),nonstifftol = 11/10),abstols,reltols;name="AutoVern8(KenCarp5(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V8KC5_11 = WorkPrecision(prob,AutoVern8(KenCarp5(),nonstifftol = 11/10),abstols,reltols;name="AutoVern8(KenCarp5(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_V8KC5_2= WorkPrecision(prob,AutoVern8(KenCarp5(),nonstifftol = 9/10),abstols,reltols;name="AutoVern8(KenCarp5(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V8KC5_3 = WorkPrecision(prob,AutoVern8(KenCarp5(),nonstifftol = 6/10),abstols,reltols;name="AutoVern8(KenCarp5(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V8KC5_4= WorkPrecision(prob,AutoVern8(KenCarp5(),nonstifftol = 3/10),abstols,reltols;name="AutoVern8(KenCarp5(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)


wp_V8Tz = WorkPrecision(prob,AutoVern8(Trapezoid(),nonstifftol = 14/10),abstols,reltols;name="AutoVern8(Trapezoid(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V8Tz_1 = WorkPrecision(prob,AutoVern8(Trapezoid(),nonstifftol = 11/10),abstols,reltols;name="AutoVern8(Trapezoid(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V8Tz_11 = WorkPrecision(prob,AutoVern8(Trapezoid(),nonstifftol = 11/10),abstols,reltols;name="AutoVern8(Trapezoid(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_V8Tz_2= WorkPrecision(prob,AutoVern8(Trapezoid(),nonstifftol = 9/10),abstols,reltols;name="AutoVern8(Trapezoid(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V8Tz_3 = WorkPrecision(prob,AutoVern8(Trapezoid(),nonstifftol = 6/10),abstols,reltols;name="AutoVern8(Trapezoid(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V8Tz_4= WorkPrecision(prob,AutoVern8(Trapezoid(),nonstifftol = 3/10),abstols,reltols;name="AutoVern8(Trapezoid(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)


wp_V8Rb23 = WorkPrecision(prob,AutoVern8(Rosenbrock23(),nonstifftol = 14/10),abstols,reltols;name="AutoVern8(Rosenbrock23(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V8Rb23_1 = WorkPrecision(prob,AutoVern8(Rosenbrock23(),nonstifftol = 11/10),abstols,reltols;name="AutoVern8(Rosenbrock23(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V8Rb23_11 = WorkPrecision(prob,AutoVern8(Rosenbrock23(),nonstifftol = 11/10),abstols,reltols;name="AutoVern8(Rosenbrock23(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_V8Rb23_2= WorkPrecision(prob,AutoVern8(Rosenbrock23(),nonstifftol = 9/10),abstols,reltols;name="AutoVern8(Rosenbrock23(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V8Rb23_3 = WorkPrecision(prob,AutoVern8(Rosenbrock23(),nonstifftol = 6/10),abstols,reltols;name="AutoVern8(Rosenbrock23(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V8Rb23_4= WorkPrecision(prob,AutoVern8(Rosenbrock23(),nonstifftol = 3/10),abstols,reltols;name="AutoVern8(Rosenbrock23(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)


wp_V8Rd3 = WorkPrecision(prob,AutoVern8(Rodas3(),nonstifftol = 14/10),abstols,reltols;name="AutoVern8(Rodas3(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V8Rd3_1 = WorkPrecision(prob,AutoVern8(Rodas3(),nonstifftol = 11/10),abstols,reltols;name="AutoVern8(Rodas3(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V8Rd3_11 = WorkPrecision(prob,AutoVern8(Rodas3(),nonstifftol = 11/10),abstols,reltols;name="AutoVern8(Rodas3(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_V8Rd3_2= WorkPrecision(prob,AutoVern8(Rodas3(),nonstifftol = 9/10),abstols,reltols;name="AutoVern8(Rodas3(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V8Rd3_3 = WorkPrecision(prob,AutoVern8(Rodas3(),nonstifftol = 6/10),abstols,reltols;name="AutoVern8(Rodas3(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V8Rd3_4= WorkPrecision(prob,AutoVern8(Rodas3(),nonstifftol = 3/10),abstols,reltols;name="AutoVern8(Rodas3(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)


wp_V8Rd4 = WorkPrecision(prob,AutoVern8(Rodas4(),nonstifftol = 14/10),abstols,reltols;name="AutoVern8(Rodas4(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V8Rd4_1 = WorkPrecision(prob,AutoVern8(Rodas4(),nonstifftol = 11/10),abstols,reltols;name="AutoVern8(Rodas4(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V8Rd4_11 = WorkPrecision(prob,AutoVern8(Rodas4(),nonstifftol = 11/10),abstols,reltols;name="AutoVern8(Rodas4(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_V8Rd4_2= WorkPrecision(prob,AutoVern8(Rodas4(),nonstifftol = 9/10),abstols,reltols;name="AutoVern8(Rodas4(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V8Rd4_3 = WorkPrecision(prob,AutoVern8(Rodas4(),nonstifftol = 6/10),abstols,reltols;name="AutoVern8(Rodas4(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V8Rd4_4= WorkPrecision(prob,AutoVern8(Rodas4(),nonstifftol = 3/10),abstols,reltols;name="AutoVern8(Rodas4(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)


wp_V8Rd5 = WorkPrecision(prob,AutoVern8(Rodas5(),nonstifftol = 14/10),abstols,reltols;name="AutoVern8(Rodas5(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V8Rd5_1 = WorkPrecision(prob,AutoVern8(Rodas5(),nonstifftol = 11/10),abstols,reltols;name="AutoVern8(Rodas5(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V8Rd5_11 = WorkPrecision(prob,AutoVern8(Rodas5(),nonstifftol = 11/10),abstols,reltols;name="AutoVern8(Rodas5(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_V8Rd5_2= WorkPrecision(prob,AutoVern8(Rodas5(),nonstifftol = 9/10),abstols,reltols;name="AutoVern8(Rodas5(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V8Rd5_3 = WorkPrecision(prob,AutoVern8(Rodas5(),nonstifftol = 6/10),abstols,reltols;name="AutoVern8(Rodas5(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V8Rd5_4= WorkPrecision(prob,AutoVern8(Rodas5(),nonstifftol = 3/10),abstols,reltols;name="AutoVern8(Rodas5(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)


# AutoVern9
wp_V9KC3 = WorkPrecision(prob,AutoVern9(KenCarp3(),nonstifftol = 14/10),abstols,reltols;name="AutoVern9(KenCarp3(),nonstifftol = 14/10))",
                      appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V9KC3_1 = WorkPrecision(prob,AutoVern9(KenCarp3(),nonstifftol = 11/10),abstols,reltols;name="AutoVern9(KenCarp3(),nonstifftol = 11/10))",
                      appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V9KC3_11 = WorkPrecision(prob,AutoVern9(KenCarp3(),nonstifftol = 11/10),abstols,reltols;name="AutoVern9(KenCarp3(),nonstifftol = 11/10))",
                      appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_V9KC3_2= WorkPrecision(prob,AutoVern9(KenCarp3(),nonstifftol = 9/10),abstols,reltols;name="AutoVern9(KenCarp3(),nonstifftol = 9/10))",
                      appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V9KC3_3 = WorkPrecision(prob,AutoVern9(KenCarp3(),nonstifftol = 6/10),abstols,reltols;name="AutoVern9(KenCarp3(),nonstifftol = 6/10))",
                      appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V9KC3_4= WorkPrecision(prob,AutoVern9(KenCarp3(),nonstifftol = 3/10),abstols,reltols;name="AutoVern9(KenCarp3(),nonstifftol = 3/10))",
                      appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)

wp_V9KC4 = WorkPrecision(prob,AutoVern9(KenCarp4(),nonstifftol = 14/10),abstols,reltols;name="AutoVern9(KenCarp4(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V9KC4_1 = WorkPrecision(prob,AutoVern9(KenCarp4(),nonstifftol = 11/10),abstols,reltols;name="AutoVern9(KenCarp4(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V9KC4_11 = WorkPrecision(prob,AutoVern9(KenCarp4(),nonstifftol = 11/10),abstols,reltols;name="AutoVern9(KenCarp4(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_V9KC4_2= WorkPrecision(prob,AutoVern9(KenCarp4(),nonstifftol = 9/10),abstols,reltols;name="AutoVern9(KenCarp4(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V9KC4_3 = WorkPrecision(prob,AutoVern9(KenCarp4(),nonstifftol = 6/10),abstols,reltols;name="AutoVern9(KenCarp4(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V9KC4_4= WorkPrecision(prob,AutoVern9(KenCarp4(),nonstifftol = 3/10),abstols,reltols;name="AutoVern9(KenCarp4(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)

wp_V9KC5 = WorkPrecision(prob,AutoVern9(KenCarp5(),nonstifftol = 14/10),abstols,reltols;name="AutoVern9(KenCarp5(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V9KC5_1 = WorkPrecision(prob,AutoVern9(KenCarp5(),nonstifftol = 11/10),abstols,reltols;name="AutoVern9(KenCarp5(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V9KC5_11 = WorkPrecision(prob,AutoVern9(KenCarp5(),nonstifftol = 11/10),abstols,reltols;name="AutoVern9(KenCarp5(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_V9KC5_2= WorkPrecision(prob,AutoVern9(KenCarp5(),nonstifftol = 9/10),abstols,reltols;name="AutoVern9(KenCarp5(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V9KC5_3 = WorkPrecision(prob,AutoVern9(KenCarp5(),nonstifftol = 6/10),abstols,reltols;name="AutoVern9(KenCarp5(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V9KC5_4= WorkPrecision(prob,AutoVern9(KenCarp5(),nonstifftol = 3/10),abstols,reltols;name="AutoVern9(KenCarp5(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)


wp_V9Tz = WorkPrecision(prob,AutoVern9(Trapezoid(),nonstifftol = 14/10),abstols,reltols;name="AutoVern9(Trapezoid(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V9Tz_1 = WorkPrecision(prob,AutoVern9(Trapezoid(),nonstifftol = 11/10),abstols,reltols;name="AutoVern9(Trapezoid(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V9Tz_11 = WorkPrecision(prob,AutoVern9(Trapezoid(),nonstifftol = 11/10),abstols,reltols;name="AutoVern9(Trapezoid(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_V9Tz_2= WorkPrecision(prob,AutoVern9(Trapezoid(),nonstifftol = 9/10),abstols,reltols;name="AutoVern9(Trapezoid(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V9Tz_3 = WorkPrecision(prob,AutoVern9(Trapezoid(),nonstifftol = 6/10),abstols,reltols;name="AutoVern9(Trapezoid(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V9Tz_4= WorkPrecision(prob,AutoVern9(Trapezoid(),nonstifftol = 3/10),abstols,reltols;name="AutoVern9(Trapezoid(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)


wp_V9Rb23 = WorkPrecision(prob,AutoVern9(Rosenbrock23(),nonstifftol = 14/10),abstols,reltols;name="AutoVern9(Rosenbrock23(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V9Rb23_1 = WorkPrecision(prob,AutoVern9(Rosenbrock23(),nonstifftol = 11/10),abstols,reltols;name="AutoVern9(Rosenbrock23(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V9Rb23_11 = WorkPrecision(prob,AutoVern9(Rosenbrock23(),nonstifftol = 11/10),abstols,reltols;name="AutoVern9(Rosenbrock23(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_V9Rb23_2= WorkPrecision(prob,AutoVern9(Rosenbrock23(),nonstifftol = 9/10),abstols,reltols;name="AutoVern9(Rosenbrock23(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V9Rb23_3 = WorkPrecision(prob,AutoVern9(Rosenbrock23(),nonstifftol = 6/10),abstols,reltols;name="AutoVern9(Rosenbrock23(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V9Rb23_4= WorkPrecision(prob,AutoVern9(Rosenbrock23(),nonstifftol = 3/10),abstols,reltols;name="AutoVern9(Rosenbrock23(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)


wp_V9Rd3 = WorkPrecision(prob,AutoVern9(Rodas3(),nonstifftol = 14/10),abstols,reltols;name="AutoVern9(Rodas3(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V9Rd3_1 = WorkPrecision(prob,AutoVern9(Rodas3(),nonstifftol = 11/10),abstols,reltols;name="AutoVern9(Rodas3(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V9Rd3_11 = WorkPrecision(prob,AutoVern9(Rodas3(),nonstifftol = 11/10),abstols,reltols;name="AutoVern9(Rodas3(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_V9Rd3_2= WorkPrecision(prob,AutoVern9(Rodas3(),nonstifftol = 9/10),abstols,reltols;name="AutoVern9(Rodas3(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V9Rd3_3 = WorkPrecision(prob,AutoVern9(Rodas3(),nonstifftol = 6/10),abstols,reltols;name="AutoVern9(Rodas3(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V9Rd3_4= WorkPrecision(prob,AutoVern9(Rodas3(),nonstifftol = 3/10),abstols,reltols;name="AutoVern9(Rodas3(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)


wp_V9Rd4 = WorkPrecision(prob,AutoVern9(Rodas4(),nonstifftol = 14/10),abstols,reltols;name="AutoVern9(Rodas4(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V9Rd4_1 = WorkPrecision(prob,AutoVern9(Rodas4(),nonstifftol = 11/10),abstols,reltols;name="AutoVern9(Rodas4(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V9Rd4_11 = WorkPrecision(prob,AutoVern9(Rodas4(),nonstifftol = 11/10),abstols,reltols;name="AutoVern9(Rodas4(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_V9Rd4_2= WorkPrecision(prob,AutoVern9(Rodas4(),nonstifftol = 9/10),abstols,reltols;name="AutoVern9(Rodas4(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V9Rd4_3 = WorkPrecision(prob,AutoVern9(Rodas4(),nonstifftol = 6/10),abstols,reltols;name="AutoVern9(Rodas4(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V9Rd4_4= WorkPrecision(prob,AutoVern9(Rodas4(),nonstifftol = 3/10),abstols,reltols;name="AutoVern9(Rodas4(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)


wp_V9Rd5 = WorkPrecision(prob,AutoVern9(Rodas5(),nonstifftol = 14/10),abstols,reltols;name="AutoVern9(Rodas5(),nonstifftol = 14/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V9Rd5_1 = WorkPrecision(prob,AutoVern9(Rodas5(),nonstifftol = 11/10),abstols,reltols;name="AutoVern9(Rodas5(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V9Rd5_11 = WorkPrecision(prob,AutoVern9(Rodas5(),nonstifftol = 11/10),abstols,reltols;name="AutoVern9(Rodas5(),nonstifftol = 11/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e6),numruns=10)
wp_V9Rd5_2= WorkPrecision(prob,AutoVern9(Rodas5(),nonstifftol = 9/10),abstols,reltols;name="AutoVern9(Rodas5(),nonstifftol = 9/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V9Rd5_3 = WorkPrecision(prob,AutoVern9(Rodas5(),nonstifftol = 6/10),abstols,reltols;name="AutoVern9(Rodas5(),nonstifftol = 6/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
wp_V9Rd5_4= WorkPrecision(prob,AutoVern9(Rodas5(),nonstifftol = 3/10),abstols,reltols;name="AutoVern9(Rodas5(),nonstifftol = 3/10))",
                        appxsol=test_sol, error_estimate=:l2,maxiters=Int(1e5),numruns=10)
