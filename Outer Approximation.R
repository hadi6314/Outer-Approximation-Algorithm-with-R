library(nloptr)
library(ROI)
library(ROI.plugin.cplex)
library(ompr) 
library(ompr.roi)
library(dplyr)
library(numDeriv)

zu = c(Inf)
zl = c(-Inf)
opt.x = c()
k.set = data.frame()
nlpr = function() 
{
  eval_f = function(x) {return(list("objective" = -x[1] - x[2] - x[3], "gradient" = c(-1,-1,-1,0)))}
  
  eval_g_ineq = function(x) 
  {
    return(rbind(c((x[2]-0.5)^2 + (x[3] - 0.5)^2 - 0.25),
                 c(x[1] - x[2]),
                 c(x[1] + x[3] + x[4] - 2)))
  }
  
  eval_jac_g0 = function(x) {return(rbind(c(0,2*(x[2] - 0.5), 2*(x[3] - 0.5),0), c(1,-1,0,0), c(1,0,1,1)))
  }
  x0 = c(0,0,0,0)
  lb = c(0,0,0,0)
  ub = c(1,Inf,Inf,5)
  local_opts = list("algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-7)
  opts = list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-7, "maxeval" = 10000,"local_opts" = local_opts)
  res = nloptr(x0 = x0, eval_f = eval_f, lb = lb, ub = ub, eval_g_ineq = eval_g_ineq, eval_jac_g_ineq = eval_jac_g0, opts = opts)
  return(res$solution)
}

master = function(crr, upper, lower) 
{
  f.grad = apply(crr, 1, function(y) {grad(function(x){-x[1]-x[2]},(y))})
  c.grad = apply(crr, 1, function(y) {grad(function(x){(x[1]-0.5)^2 + (x[2]-0.5)^2},(y))})
  coef1 = c(1,-1,0,0)
  coef2 = c(1,0,1,1)
  model = MIPModel()
  model = add_variable(model, x, type = "binary")
  model = add_variable(model, y[i], i = 1:2, type = "continuous")
  model = add_variable(model, z, type = "integer")
  model = add_variable(model, eta, type = "continuous")
  model = set_bounds(model, y[i], i = 1:2, lb = 0, ub = Inf) 
  model = set_bounds(model, z, lb = 0, ub = 5)
  model = set_bounds(model, eta, lb = -Inf, ub = Inf)
  model = set_objective(model, eta - x - z, "min")
  model = add_constraint(model, coef1[1]*x + coef1[2]*y[1] <= 0)
  model = add_constraint(model, x + y[2] + z <= 2)
  model = add_constraint(model, eta <= upper - .Machine$double.eps)
  model = add_constraint(model, eta >= lower)
  for (i in 1:nrow(crr)) 
  {
    model = add_constraint(model, (-crr[i,1]-crr[i,2]) + f.grad[1,i]*(y[1]-crr[i,1]) + f.grad[2,i]*(y[2]-crr[i,2]) <= eta)
    model = add_constraint(model, ((crr[i,1]-0.5)^2 + (crr[i,2]-0.5)^2) - 0.25 + c.grad[1,i]*(y[1]-crr[i,1]) + c.grad[2,i]*(y[2]-crr[i,2]) <= 0)
  }
  
  results = solve_model(model, with_ROI("cplex", verbose = TRUE))
  if (solver_status(results) == 'infeasible') {return(c('infeasible'))}
  return(c(as.numeric(results$objective_value), as.numeric(results$solution[c(2,5)]))) #solution[c(1,4)
}
#x.int = c(1,3)
nlp = function(x.int) 
{
  eval_f = function(x) {return(list("objective" = -x.int[1] - x[1] - x[2] - x.int[2], "gradient" = c(-1,-1)))}
  eval_g_ineq = function(x) 
  {
    return(rbind(c((x[1]-0.5)^2 + (x[2] - 0.5)^2 - 0.25),
                 c(x.int[1] - x[1]),
                 c(x.int[1] + x[2] + x.int[2] - 2)))
  }
  
  eval_jac_g0 = function(x) {return(rbind(c(2*(x[1] - 0.5), 2*(x[2] - 0.5)), c(-1,0), c(1,0)))}
  x0 = c(0,0)
  lb = c(0,0)
  ub = c(Inf,Inf)
  local_opts = list("algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-7)
  opts = list( "algorithm" = "NLOPT_LD_MMA", "xtol_abs" = .Machine$double.eps, "maxeval" = 1000000,"local_opts" = local_opts)
  res = nloptr(x0 = x0, eval_f = eval_f, lb = lb, ub = ub, 
               eval_g_ineq = eval_g_ineq, eval_jac_g_ineq = eval_jac_g0, opts = opts)
  #if ()
  return(res$solution) 
}
nlpf = function(x.int) 
{
  eval_f = function(x) {return(list("objective" = (x[1]-0.5)^2 + (x[2] - 0.5)^2 - 0.25,
                                    "gradient" = c(2*(x[1] - 0.5), 2*(x[2] - 0.5))))}
  
  eval_g_ineq = function(x) {return(rbind(c(x.int[1] - x[1]), c(x.int[1] + x[2] + x.int[2] - 2)))}
  
  eval_jac_g0 = function(x) {return(rbind(c(-1,0), c(1,0)))}
  x0 = c(0,0)
  lb = c(0,0)
  ub = c(Inf,Inf)
  local_opts = list("algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-7)
  opts = list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-7, "maxeval" = 10000,"local_opts" = local_opts)
  res = nloptr(x0 = x0, eval_f = eval_f, lb = lb, ub = ub, 
               eval_g_ineq = eval_g_ineq, eval_jac_g_ineq = eval_jac_g0, opts = opts)
  return(res$solution)
}

feasibility = function(ints, cons) 
{
  return (((cons[1]-0.5)^2 + (cons[2]-0.5)^2 <= 0.25) & (ints[1] - cons[1] <= 0) & 
            (ints[1] + cons[2] + ints[2] <= 2) & (ints[1] <= 1) & (ints[2] <= 5) & (ints[1] >= 0) & 
            (ints[2] >= 0) & (cons[1] >= 0) & (cons[2] >= 0))
}

init.x = c(0,0,0,0)
#k.set = rbind(k.set,init.x)
while (zu[length(zu)] - zl[length(zl)] > 1.0e-3) 
{
  nlp.solution = nlp(c(init.x[1],init.x[4]))
  zu = c(zu, -init.x[1]-nlp.solution[1]-nlp.solution[2]-init.x[4])
  k.set = rbind(k.set,nlp.solution)
  
  master.solution = master(k.set,zu[length(zu)],zl[length(zl)])
  if (master.solution[1] != 'infeasible') {zl = c(zl, master.solution[1])}
  else {break()}
  current.solution = c(init.x[1],nlp.solution[1],nlp.solution[2],init.x[4])
  init.x = c(master.solution[2],nlp.solution[1],nlp.solution[2],master.solution[3])
  
  #nlp.solution = nlp(master.solution[-1])
}
print(paste(c('The optimal solution is',current.solution), collapse=" "))
print(paste(c("The optimal objective function is",zu[length(zu)-1]), collapse = " "))
