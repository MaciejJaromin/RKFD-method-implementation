function rk_step(param, h, tn, yn, f)
	s = size(param)[1]-1
	K = []
	for i in 1:s
		t = tn + param[i,1]*h
		y = yn
		for j in 2:i
			y = y + h*param[i,j]*K[j-1]
		end
		push!(K,f(t,y))
	end
	yn1 = yn
	for i in 1:s
		yn1 = yn1 + h*param[s+1,i+1]*K[i]
	end
	return yn1
end

function rkfd_step(param, h, tn, yn, f)
	s = size(param)[1]-4
	Y = []
	for i in 1:s
		ci = param[i,1]
		sum = BigFloat(0)
		for j in 2:i
			sum = sum + param[i,j] * f(tn + param[j-1,1]*h, Y[j-1])
		end
		push!(Y,yn[1] + ci*h*yn[2] + ci^2 *h^2 * yn[3]/2 + ci^3 * h^3 * yn[4]/6 + h^4 * sum)
	end
	yn1 = BigFloat[0, 0, 0, 0];
	yn1[1] = yn[1] + h*yn[2] + h^2 * yn[3]/2 + h^3 * yn[4]/6
	yn1[2] = yn[2] + h*yn[3] + h^2 * yn[4]/2
	yn1[3] = yn[3] + h*yn[4]
	yn1[4] = yn[4]
	for i in 1:s
		aux = f(tn + param[i,1]*h, Y[i])
		yn1[1] = yn1[1] + h^4 * param[s+1,i+1] * aux
		yn1[2] = yn1[2] + h^3 * param[s+2,i+1] * aux
		yn1[3] = yn1[3] + h^2 * param[s+3,i+1] * aux
		yn1[4] = yn1[4] + h * param[s+4,i+1] * aux
	end
	return yn1
end

d(a,b) = BigFloat(a)/BigFloat(b)

param_RK4 = BigFloat[
[0 0 0 0 0];
[d(1,2) d(1,2) 0 0 0];
[d(1,2) 0 d(1,2) 0 0];
[1 0 0 1 0];
[0 d(1,6) d(1,3) d(1,3) d(1,6)]]

param_RKFD4 = BigFloat[
[0 0 0 0];
[d(4,11) d(-1,5) 0 0];
[d(17,20) d(19,125) d(19,125) 0];
[0 d(17,200) d(-7,75) d(1,20)];
[0 d(1,18) d(209,1926) d(5,1926)];
[0 d(47,408) d(847,2568) d(100,1819)];
[0 d(47,408) d(1331,2568) d(2000,5457)]]	

true_y(t) = sin(t)
y0 = [BigFloat(0), BigFloat(1), BigFloat(0), BigFloat(-1)]



function calc_one_step_error_array(param, h, f, y0, t0, y_orig)
	y1 = rk_step(param,h,t0,y0,f)
	real_value = y_orig(t0+h)
	return abs(real_value-y1[1])
end

function test_error_RK4(n)
	f(t,y) = [y[2], y[3], y[4], y[1]^2 + cos(t)^2 + sin(t) - 1]
	h_space = logspace(0, -n, n+1)
	prev_err = BigFloat(0)
	for h in h_space
		err = calc_one_step_error_array(param_RK4, BigFloat(h), f, y0, BigFloat(0), true_y)
		print("h = ", h, "\n")
		print("Error in one step: ", err, "\n")
		print("Error reduced by a factor of: ", prev_err/err, "\n")
		print("\n")
		prev_err = err
	end
end

function calc_one_step_error_RKFD(param, h, f, y0, t0, y_orig)
	y1 = rkfd_step(param,h,t0,y0,f)
	real_value = y_orig(t0+h)
	return abs(real_value-y1[1])
end

function test_error_RKFD4(n)
	f(t,y) = y^2 + cos(t)^2 + sin(t) - 1
	h_space = logspace(0, -n, n+1)
	prev_err = BigFloat(0)
	for h in h_space
		err = calc_one_step_error_RKFD(param_RKFD4, BigFloat(h), f, y0, BigFloat(0), true_y)
		print("h = ", h, "\n")
		print("Error in one step: ", err, "\n")
		print("Error reduced by a factor of: ", prev_err/err, "\n")
		print("\n")
		prev_err = err
	end
end
	
function benchmark_RK4(h)
	f(t,y) = [y[2], y[3], y[4], y[1]^2 + cos(t)^2 + sin(t) - 1]
	t = BigFloat(0)
	t_end = BigFloat(10)
	h = BigFloat(h)
	yn = [BigFloat(0), BigFloat(1), BigFloat(0), BigFloat(-1)]
	tic()
	while t < t_end
		yn = rk_step(param_RK4, h, t, yn, f)
		t = t+h
	end
	toc()
	print("Error: ", abs(sin(t)-yn[1]),"\n")
end

function benchmark_RKFD4(h)
	f(t,y) = y^2 + cos(t)^2 + sin(t) - BigFloat(1)
	t = BigFloat(0)
	t_end = BigFloat(10)
	h = BigFloat(h)
	yn = [BigFloat(0), BigFloat(1), BigFloat(0), BigFloat(-1)]
	tic()
	while t < t_end
		yn = rkfd_step(param_RKFD4, h, t, yn, f)
		t = t+h
	end
	toc()
	print("Error: ", abs(sin(t)-yn[1]),"\n")
end
function benchmark_RKFD4(h)
	f(t,y) = y^2 + cos(t)^2 + sin(t) - BigFloat(1)
	t = BigFloat(0)
	t_end = BigFloat(10)
	h = BigFloat(h)
	yn = [BigFloat(0), BigFloat(1), BigFloat(0), BigFloat(-1)]
	tic()
	while t < t_end
		yn = rkfd_step(param_RKFD4, h, t, yn, f)
		t = t+h
	end
	toc()
	print("Error: ", abs(sin(t)-yn[1]),"\n")
end