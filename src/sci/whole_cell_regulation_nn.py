import numpy as np
import pandas as pd
# In[38]:
expression = pd.read_csv('../../input/nn/data.csv')
del expression["Unnamed: 0"] # This holds the gene names just so you know...
meta = pd.read_csv('../../input/nn/meta.csv')
del meta["Unnamed: 0"]
# In[39]:
# expression
stress_one_hot = pd.get_dummies(meta["Stress"])
stress_one_hot.shape
# In[40]:
expression = expression.transpose()
expression.shape
# In[55]:
from __future__ import print_function

import numpy as np
import tflearn

target = expression
target = np.array(expression, dtype='float32')
data = stress_one_hot
data = np.array(data, dtype='float32')

# Build neural network
net = tflearn.input_data(shape=[None, 9])
net = tflearn.fully_connected(net, 32)
net = tflearn.fully_connected(net, 6692)
net = tflearn.regression(net)

# Define model
model = tflearn.DNN(net)
# Start training (apply gradient descent algorithm)
model.fit(data, target, n_epoch=10, batch_size=16, show_metric=True)

# Let's create some data for DiCaprio and Winslet
