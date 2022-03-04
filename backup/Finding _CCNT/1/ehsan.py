import matlab.engine
eng = matlab.engine.start_matlab()
eng.nsga2(nargout=0)
