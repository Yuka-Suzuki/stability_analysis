using DelimitedFiles, Random, Distributions, StatsBase, Statistics, LinearAlgebra, LightGraphs

function Initial(n_sp,n_com,initJ)
	#=
	pinit: the initial Pik for all i and k. pinit*n_sp has to be <= 1.
	The example they used in the paper was n_sp = 20, n_com = 20
	initJ: the number of occupied sites in each community at initial state
	=#
	Q = zeros(n_sp, n_com) .+ initJ;
	return Q
end;

function ReadMat(filename,d)
	#=
	symmetric matrix (undirected) filled with 0 or 1
	d: per capita emigration rate. sum(mat,dims=2) = d
	elements (average growth rates) in the file are between 0.5 and 1.5
	=#
	mat = readdlm(filename, ' ', Float64);
	mat = mat ./ sum(mat,dims=2);
	mat = mat .* d;
	return mat
end;

function BFS(A)
        # g: simple directed graph
	g = DiGraph(A);
        n_com = nv(g);
        q = zeros(Integer,0); # empty vector for queuing nodes
        dist = zeros(n_com); # vector to contain the distances of nodes from the starting node
        # start_v = 1; # starting node
        cent = closeness_centrality(g);
        start_v = findall(cent.==maximum(cent))[1]; # choose the node with highest closeness centrality as the starting node
        dist[start_v] = 0; # the distance at starting node is zero
        discov = zeros(0); # vector to contain the nodes that have been discovered
        append!(discov,start_v);
        append!(q,start_v);
        while length(q) > 0
                v = q[1]; # choose one node from the queue
                deleteat!(q,1);
                if length(discov) == n_com
                        return dist
                else
                        for w in neighbors(g,v)
                                if !(w in discov)
                                        append!(discov,w);
                                        append!(q,w);
                                        dist[w] = dist[v] + 1;
                                end;
                        end;
                end;
        end;
end;

function sort_env(dist,env_mat)
	env = zeros(size(env_mat));
	sortednodes = sortperm(dist);
	for (i,node) in enumerate(sortednodes)
		env[node] = env_mat[i]
	end;
	return env
end;

function noises(rho_eta_fluc,L,r_std)
	noises = zeros(L);
	dnorm = Normal(0,r_std);
	prenoise = rand(Truncated(dnorm,-1,1));
	noises[1] = prenoise;
	for j in 2:L
		cnoise = rho_eta_fluc .* prenoise .+ (1 - rho_eta_fluc) .* rand(dnorm);
		noises[j] = cnoise;
		prenoise = cnoise;
	end;
	return noises
end;

function Update(Q,A,r,c,K,rho_eta_fluc,n_com,n_sp,h,dist,meanK,r_std)
	#=
	<Variables>
	Matrix Q (Qik): the number of sites occupied by species i in community k. Each site can take a single individual only. !! Different from relative species abundance
	<Parameters>
	Matrix r (rik): the intrinsic growth rate of species i in community l.
	Matrix K (Kik): the carrying capacity for species i in community k.
	Matrix  A (alk): the dispersal rate from community l to k. diag(A) = 0.
	=#
	t = 0.05; # step size in Euler's method
	for i in 1:Int(1 ./ t)
		grow = r.*Q.*(1 .- (Q .+ c*Q)./K); # growth ### check this
		disp = -sum(A,dims=2)'.*Q .+ Q*A; # dispersal (-outbound + inbound)
		Q = Q .+ t.*(grow .+ disp);
        	if any(Q.<0)
            		for indices in findall(x->x<0,Q)
                		Q[indices[1],indices[2]] = 0;
            		end;
        	end;
	end;
	return Q 
end;


function Diversity(P, n_com)
	#=
	Jost 2006,2007 Multiplicative partition of the effective number of species of 
    the Gini-Simpson index, with equal weights for all communities.
    0 < gamma/alpha < 1, 1 < beta < n_com
	=#
	rP = P ./ sum(P, dims=1);
	# the above division gives nan for all elms in row i when sp.i went extinct from all communities. Replace this row with np.zeros(n_com) in that case.
	for i in (1:size(rP)[1])
		if all(broadcast(isnan,rP)) && sum(P[i,:]) == 0
			rP[i,:] = zeros(n_com);
		end;
	end;
	w = 1/n_com;
	lam_g = (sum(rP,dims=2)' * sum(rP,dims=2))[1,1] * w^2;
	lam_a = sum(rP.*rP) ./ n_com;
	lam_b = lam_g ./ lam_a;
	gamma = 1 ./ lam_g;
	alpha = 1 ./ lam_a;
	beta  = 1 ./ lam_b;
	return alpha, beta, gamma
end;

function Variability2(Q_ls,n_sp,n_com)
#= Wang2019stability:
The method in Wang 2014 is for the special case where local communities all have the same number and abundance of species.
=#
	# Variabilities at local or regional and population or community level
	mu_ik = mean(Q_ls,dims=3); # temporal mean biomass of species i in patch k
	mu = mean(sum(Q_ls,dims=[1,2])); # temporal mean of total metacommunity biomass
	nu = var(sum(Q_ls,dims=[1,2])); # temporal variance of total metacommunity biomass
	CVcr = sqrt(nu) ./ mu; # variability at community and regional level
	
	nu_kk = var(sum(Q_ls,dims=1),dims=3); # temporal variance of each local community (ignoring species difference)
	CVcl = sum(broadcast(sqrt,nu_kk)) ./ mu; # sum(CVck .* (mu_k ./ mu)); # local-scale average community variability

	nu_ii = var(sum(Q_ls,dims=2),dims=3);
	CVsr = sum(broadcast(sqrt,nu_ii)) ./ mu; # regional-scale average species variability

	nu_iikk = var(Q_ls,dims=3); # temporal variance of biomass of sp i in com k
	CVsl = sum(broadcast(sqrt,nu_iikk)) ./ mu; # local-scale average species variability

	# synchrony measures
	phi_r = CVcr ./ CVsr;
	phi_l = CVcl ./ CVsl;
	phi_c = CVcr ./ CVcl;
	phi_s = CVsr ./ CVsl;
	return CVcr,CVcl,CVsr,CVsl,phi_r,phi_l,phi_c,phi_s
end;

function env2K(env_mat,h,dist,n_com)
	X = collect(0:1/19:1);

	nE = reshape(sort_env(dist,env_mat),1,n_com);
	nD = (nE .- X).^2
	nF = broadcast(exp,-nD ./ (2*h));
	K = nF .* 2; 
	return K
end;

function Immigration(Q,initJ)
	# constant immigration into all patches equally
	Q .+= initJ * 0.001;
	return Q
end;

# Parameter and directory setting
maxtime = parse(Int64,ARGS[7]);
r_std = parse(Float64,ARGS[8]);
statrange = 2000;
n_sp = 20;
pinit = 1 ./ n_sp;
n_com = parse(Int64,ARGS[11]);
n_diffenv = parse(Int64,ARGS[12]);
name = ARGS[1];
netdir = ARGS[2];
rho_eta_hetero = parse(Float64,ARGS[3]); # spatial autocorrelation in environmental condition affecting growth rate
rho_eta_fluc = parse(Float64,ARGS[6]); # spatial autocorrelation in stochastic environmental responses through growth rate
d = parse(Float64,ARGS[4]);
h = parse(Float64,ARGS[9]);
outfile = ARGS[13];

nRep = 100;
filename = string("~/Archive/",netdir,"/",name); # set the path to the Archive directory
Random.seed!(1234);
A = ReadMat(filename,d);
dist = BFS(A);
env_mat = readdlm(string("~/Archive/environments/N",n_com,"/env_rhoeta",rho_eta_hetero,"_e",n_diffenv,".txt"),' ');
env_sorted = reshape(sort_env(dist,env_mat),1,n_com);
X = collect(0:1/(n_sp-1):1);
D = (env_sorted .- X).^2;
F = broadcast(exp,-D ./ (2*h));
meanK = F .* 2; 
r = parse(Float64,ARGS[10]);
c = parse(Float64,ARGS[5]);
if c > 0
	c = zeros(n_sp,n_sp) .+ c;
	c = c .- (zeros(size(c)) +  Matrix{Float64}(I,size(c))).*diag(c); # interspecific interaction
end;

Random.seed!(0);
for rep in 1:nRep 
	for neu in [false]
		initJ = meanK ./ 2;
		Q_ls = zeros(n_sp,n_com,0);
		Q = Initial(n_sp,n_com,initJ);
		Preg = vec(sum(Q,dims=2) ./ Float64(sum(Q)));
		adiv = zeros(0);
		bdiv = zeros(0);
		gdiv = zeros(0);
		tot_biomass = zeros(0);
		env_noises = 0;
		current_env = env_mat;
		for g in (1:maxtime)
			if g % 20 == 1
				env_noises = noises(rho_eta_fluc,length(env_mat),r_std);
				current_env = current_env .+ env_noises;
				if !all(current_env .<= 2)
					current_env[current_env .> 2] .= 1.8;
				end;
				if !all(current_env .>= -1)
					current_env[current_env .< -1] .= -0.8;
				end;
			end;
			K = env2K(current_env,h,dist,n_com);
			Q = Update(Q,A,r,c,K,rho_eta_fluc,n_com,n_sp,h,dist,meanK,r_std);
			Q = Immigration(Q,initJ);
			if sum(Q) <= 0
				println("metapopulation extinction ",g);
				append!(tot_biomass,-1)
				break
			# calculate diversity and variability if simulation has reached (maxtime-statrange+1)th generations
			elseif g >= (maxtime-statrange+1)
				alpha,beta,gamma = Diversity(Q,n_com);
				append!(adiv,alpha);
				append!(bdiv,beta);
				append!(gdiv,gamma);
				biomass = sum(Q,dims=1);
				Q_ls = cat(Q_ls,Q,dims=3);
				append!(tot_biomass,sum(biomass));
			end;
		end;
		if sum(Q_ls) == 0 # extinction event occurred before simulation runs long enough to calculate variance
			CVcr,CVcl,CVsr,CVsl,phi_r,phi_l,phi_c,phi_s = "NaN","NaN","NaN","NaN","NaN","NaN","NaN","NaN";
		else
			CVcr,CVcl,CVsr,CVsl,phi_r,phi_l,phi_c,phi_s = Variability2(Q_ls,n_sp,n_com);
		end;
		if tot_biomass[end] == -1
			open(string(outfile),"a") do io
				write(io, string(name, " ", CVcr," ",CVcl," ",CVsr," ",CVsl," ",phi_r," ",phi_l," ",phi_c," ",phi_s," ","extinction"," ",mean(adiv)," ",mean(gdiv)," ",mean(bdiv)," ",var(adiv)," ",var(gdiv)," ",var(bdiv)," ",sum(meanK)," ",sum(c)," ",sum(env_mat)," ",rho_eta_hetero," ",rho_eta_fluc," ",d," ",r_std," ",neu,"\n"));
			end;
		else
			open(string(outfile),"a") do io
				write(io, string(name, " ", CVcr," ",CVcl," ",CVsr," ",CVsl," ",phi_r," ",phi_l," ",phi_c," ",phi_s," ",tot_biomass[end]," ",mean(adiv)," ",mean(gdiv)," ",mean(bdiv)," ",var(adiv)," ",var(gdiv)," ",var(bdiv)," ",sum(meanK)," ",sum(c)," ",sum(env_mat)," ",rho_eta_hetero," ",rho_eta_fluc," ",d," ",r_std," ",neu,"\n"));
			end;
		end;
	end;
end;
