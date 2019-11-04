CXXFLAGS=-fopenmp -O2 -Iinclude

rt: dt_kick dt_nokick
	$(CXX) $(CXXFLAGS) ray_trace/ray_trace.cpp -o rt

dt_kick:
	$(CXX) $(CXXFLAGS) delta_tracking/kicking/delta_tracking.cpp -o dt_kick

dt_nokick:
	$(CXX) $(CXXFLAGS) delta_tracking/no_kicking/delta_tracking.cpp -o dt_nokick

clean:
	rm dt_kick dt_nokick rt
