import matlab.engine
eng = matlab.engine.start_matlab()
zz = eng.helix(float(1),float(2),float(2),float(2),float(2),float(7))
print(zz)
