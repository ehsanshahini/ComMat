
import nsga_single as ns
import joblib
from numpy import savetxt

R_input = input("Enter mean diameter of the coil: ")

R_input = float(R_input)/2

theta_input = float(input("Enter pitch angle of the coil: "))

D_input = float(input("Enter diameter of the carbon nanotube: "))


# predicting the stress from Random Forest algorithm
regressor_stress = joblib.load("./RF_EL_STRESS.joblib")

stress = regressor_stress.predict([[R_input,theta_input,D_input]])

print("The predicted stress is", stress[0], "GPa")


# predicting the strain from Random Forest algorithm
regressor_strain = joblib.load("./RF_EL_STRAIN.joblib")

strain = regressor_strain.predict([[R_input,theta_input,D_input]])

print("The predicted strain is", strain[0])

with open("Predicted_results.txt", "w") as f:
    f.write("The predicted results based on random forest algorithm is as follow:\n")
    f.write("Stress = " + str(stress[0]) + " GPa\n")
    f.write("Strain =  " + str(strain[0]))

####-----------------optimization and getting the XYZ coordinate of CNC-----

Question = input("Do you want to find XYZ coordinate of the CNC? Type 'y' if you like to (it might take a long time to produce the result).")
if Question == ("y"):
    index = ns.nsga(R_input, theta_input, D_input)
    print("The best answer has indices of ", index)
    XYZ = ns.position(index)
    savetxt('CNC_coordinates.xyz', XYZ, delimiter=' ')