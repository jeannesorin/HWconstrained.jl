module HWconstrained

	using JuMP, NLopt, DataFrames, Ipopt
	using LinearAlgebra

	export data, table_NLopt, table_JuMP, max_JuMP, obj, constr, max_NLopt


	function data(a=0.5)
		na = 3
		p = ones(3)
		e = [2.0, 0.0, 0.0]
		z = zeros(16,3)
		z_2 = [0.72, 0.92, 1.12, 1.32]
		z_3 = [0.86, 0.96, 1.06, 1.16]
		z[:,1] .= 1.0
		for i=1:4
		z[1+(4*(i-1)):4+(4*(i-1)), 2] .= z_2[i]
		end
		for j = 1:4
			z[j,3] = z_3[j]
			z[j+4,3] = z_3[j]
			z[j+4*2,3] = z_3[j]
			z[j+4*3,3] = z_3[j]
		end
		ns = length(z_2)
		nss = length(z_2) * length(z_3)
		pi = 1/nss
		nc = 1
		return Dict("a"=>a,"na"=>na,"nc"=>nc,"ns"=>ns,"nss"=>nss,"e"=>e,"p"=>p,"z"=>z,"pi"=>pi)
	end



# DOES NOT WORK
	function max_JuMP(a=0.5)
		m = Model(with_optimizer(Ipopt.Optimizer))
		a = 0.5
		N = 3
		@variable(m, omega[1:N] >= 0 )

		#	@variable(m, c >= 0)

		@objective(m, Max,
			-exp(-a*c) + data()["pi"] * sum(-exp.(-a .*
			sum(data()["z"][:,i] * x[i] for i = 1:3))))

		@constraint(m, c + (sum(data()["p"] .* (X .- data()["e"]))) <= 0)

		print(m)

		status = optimize!(m)


		return Dict("obj"=>objective_value(m),"c"=>value(c),"omegas"=>[value(omega[i]) for i in 1:length(omega)])
	end




	function table_JuMP()
		d = DataFrame(a=[0.5;1.0;5.0],c = zeros(3),omega1=zeros(3),omega2=zeros(3),omega3=zeros(3),fval=zeros(3))
		for i in 1:nrow(d)
			xx = max_JuMP(d[i,:a])
			d[i,:c] = xx["c"]
			d[i,:omega1] = xx["omegas"][1]
			d[i,:omega2] = xx["omegas"][2]
			d[i,:omega3] = xx["omegas"][3]
			d[i,:fval] = xx["obj"]
		end
		return d
	end






## Does not work
	function obj(x::Vector,grad::Vector,data::Dict)
		if length(grad) > 0
			# grad wrt x1
			a = 0.5
			grad[1] = (-1) * sum(data["pi"] * (-a * data["z"][:,1]) * (- exp(-a * sum(data["z"][:,1] * x[1] +
																				data["z"][:,2] * x[2] +
																				data["z"][:,3] * x[3]))))
			# grad wrt x2
			grad[2] = (-1) * sum(data["pi"] * (-a * data["z"][:,2]) * (- exp(-a * sum(data["z"][:,1] * x[1] +
																				data["z"][:,2] * x[2] +
																				data["z"][:,3] * x[3]))))
			# grad wrt x3
			grad[3] = (-1) * sum(data["pi"] * (-a * data["z"][:,3]) * (- exp(-a * sum(data["z"][:,1] * x[1] +
																				data["z"][:,2] * x[2] +
																				data["z"][:,3] * x[3]))))
			# grad wrt c
			grad[4] = (-1) *  (-a) * (- exp(-a * c))
		end
		return (-1) * (-exp(-a*c) + data["pi"] * (- exp(-a * sum(data["z"][:,1] * x[1] +
										data["z"][:,1] * x[2] +
										data["z"][:,1] * x[3]))))
	end



	function constr(x::Vector,grad::Vector,data::Dict)
		if length(grad) > 0
			grad[1] = data["p"][1]
			grad[2] = data["p"][2]
			grad[3] = data["p"][3]
			#grad[4] = 1
		end
		return (c +
			data["p"][1] * (x[1] - data["e"][1]) +
			data["p"][2] * (x[2] - data["e"][2]) +
			data["p"][3] * (x[3] - data["e"][3]))

	end



### BOUND ISSUE --> problem on this function
	function max_NLopt(a=0.5)

		opt = Opt(:LD_MMA, 3)
		lower_bounds!(opt, [0., 0., 0.])
		xtol_rel!(opt,1e-4)
		min_objective!(opt,(x,grad) -> obj(x,grad,data()))
		inequality_constraint!(opt, (x,grad) -> constr(x,grad))
		ftol_rel!(opt,1e-9)
		(minfunc,minx,ret) = NLopt.optimize(opt, [0.0, 0.0, 0.0])
end
#Bounds error



	function table_NLopt()
		d = DataFrame(a=[0.5;1.0;5.0],c = zeros(3),omega1=zeros(3),omega2=zeros(3),omega3=zeros(3),fval=zeros(3))
		for i in 1:nrow(d)
			xx = max_NLopt(d[i,:a])
			for j in 2:ncol(d)-1
				d[i,j] = xx[2][j-1]
			end
			d[i,end] = xx[1]
		end
		return d
	end
table_NLopt()


end # module
