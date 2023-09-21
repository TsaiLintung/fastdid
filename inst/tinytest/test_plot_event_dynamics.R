

dt <- sim_did(100, 10)[["dt"]]
dynamic_est <- generate_est(dt, "dynamic")
expect_silent(dynamic_est |> plot_event_dynamics(), info = "no error in single outcome 2 stratify")