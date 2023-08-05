import torch
import numpy as np

data = [[1,2], [3,4]]
x_data = torch.tensor(data)

print(x_data)

tensor = torch.ones(4, 4)

tensor[:,1] = 0

print(tensor)