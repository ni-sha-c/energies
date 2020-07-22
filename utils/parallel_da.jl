for i = 1:2
	filename = string("../data/rijke_test_", string(i), ".jld") 
	put!(RemoteChannel(i), filename) 
	@spawnat i assimilate_parameter_and_trajectory(filename)
end

	
