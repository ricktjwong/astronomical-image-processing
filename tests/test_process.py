import numpy as np
import main

mock_data = np.array([0, 2, 0], [0, 4, 1], [0, 6, 1])
brightest = main.get_brightest(mock_data)
print(brightest)