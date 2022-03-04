import matlab.engine
eng = matlab.engine.start_matlab()
zz = eng.helix(2,3)
print(zz)
