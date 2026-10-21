import uproot
import numpy as np


# Ensure the file path is correct
file_path = "Landau2_LB5_trgCH4_VCB.root"

try:
    with uproot.open(file_path) as file:
        # Try different methods to list contents
        print("File keys:", file.keys())
        print("File items:", list(file.items()))
        print("File classes:", list(file.classes()))
except Exception as e:
    print("Error reading file:", e)






gain=1.136449e+06





with uproot.open("Landau2_LB5_trgCH4_VCB.root") as file:    
  print(file.keys())
  tree = file["PMT1data;29"] 
  data = tree.arrays("corrected").tolist()
  
  # Extracting values from the list of dictionaries
  corrected_values = np.array([item["corrected"] for item in data])
  
  # Filtering out values smaller than -0.001 and larger than 0.2
  filtered_values = corrected_values[(corrected_values >= -0.001) & (corrected_values <= 0.2)]
  
  filtered_values=filtered_values/gain

  sum=np.sum(np.array(filtered_values))
  print("the sum of the values is:", format(sum, '.2e'))

