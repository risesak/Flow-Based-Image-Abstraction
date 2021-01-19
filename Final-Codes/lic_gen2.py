import numpy as np;
import matplotlib.pyplot as plt;
import scipy.io as sio;

from lic_cython import lic_internal;

print("We have assumed that the name of variable in mat file is 'ETF_output'.")
print("Enter file name: ")
file_name = input();
# print(file_name);
data = sio.loadmat(file_name); # load mat file as dictionary
data = np.array(data["ETF_output"]); # load variable having vactor field

texture = np.random.rand(data.shape[1],data.shape[0]).astype(np.float32)

kernellen=31
kernel = np.sin(np.arange(kernellen)*np.pi/kernellen)
kernel = kernel.astype(np.float32)

lic_array = lic_internal.line_integral_convolution(data.astype(np.float32),texture,kernel);
plt.imshow(lic_array,origin='upper',cmap='gray');
plt.show();