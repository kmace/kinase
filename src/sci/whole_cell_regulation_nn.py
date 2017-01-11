import numpy as np
import pandas as pd
import tflearn
import tensorflow as tf

meta = pd.read_csv('../../input/nn/meta.csv')
experiments = meta["Unnamed: 0"]
del meta["Unnamed: 0"]

stress_one_hot = pd.get_dummies(meta["train"])

kinase_deletion = pd.get_dummies(meta["Strain"])
kinase_one_hot = 1 - kinase_deletion

expression = pd.read_csv('../../input/nn/data.csv')
genes = expression["Unnamed: 0"]
del expression["Unnamed: 0"] # This holds the gene names just so you know...

expression = expression.transpose()

# Set up data for tensorflow
# Gene expression
target = expression
target = np.array(expression, dtype='float32')
target_mean = target.mean(axis=0, keepdims=True) # I had this wrong before, the time it worked... weird axis was 1
target_std = target.std(axis=0, keepdims=True)
target = target - target_mean
target = target / target_std

# Stress information
data1 = stress_one_hot
data1 = np.array(data1, dtype='float32')
data_mean = data1.mean(axis=0, keepdims=True)
data_std = data1.std(axis=0, keepdims=True)
data1 = data1 - data_mean
data1 = data1 / data_std

# Kinase information
data2 = kinase_one_hot
data2 = np.array(data2, dtype='float32')

# For Reference
# data1.shape
# #(301, 10)
# data2.shape
# #(301, 29)


# Build the Neural Network

num_stresses = 10
num_kinase = 29
num_transcription_factors = 200
num_genes = 6692

# Build neural network
# Input variables (10)
# Which Node to dropout (32)
stress = tflearn.input_data(shape=[None, num_stresses])
kinase_deletion = tflearn.input_data(shape=[None, num_kinase])

# This is the layer that I want to perform selective dropout on,
# I should be able to specify which of the 32 nodes should output zero
# based on a 1X32 vector of ones and zeros.
kinase = tflearn.fully_connected(stress, num_kinase, activation='relu')
kinase_dropout = tf.mul(kinase, kinase_deletion)

transcription_factor = tflearn.fully_connected(kinase_dropout, num_transcription_factors, activation='relu')

gene = tflearn.fully_connected(transcription_factor, num_genes, activation='linear')

adam = tflearn.Adam(learning_rate=0.00001, beta1=0.99)

regression = tflearn.regression(gene, optimizer=adam, loss='mean_square', metric='R2')

# Define model
model = tflearn.DNN(regression, tensorboard_verbose=1)

# Start training (apply gradient descent algorithm)
model.fit([data1, data2], target, n_epoch=20000, show_metric=True, shuffle=True)#,validation_set=0.05)


# Plot results

import matplotlib.pyplot as plt # side-stepping mpl backend
import seaborn as sns



# In[ ]:

get_ipython().magic('matplotlib inline')
df = pd.DataFrame()
df['nn_beta'] = model.get_weights(gene.b)
df['gene_mean'] = np.mean(target,0)
sns.lmplot('nn_beta','gene_mean',data=df)


# In[ ]:

get_ipython().magic('matplotlib inline')
stress_to_kinase = pd.DataFrame(model.get_weights(kinase.W))
stress_to_kinase.rows = stress_one_hot.columns

#sns.clustermap(stress_to_kinase)#, xticklabels=stress_to_kinase.rows)
#
pd.melt(stress_to_kinase)
test = stress_to_kinase.T
test.columns = stress_one_hot.columns
sns.clustermap(test)


# Let us now plot the weight matrix between kinase outputs and transcription factors, each column is a TF (100) and each row is a kinase (32)

# In[ ]:

get_ipython().magic('matplotlib inline')
kinase_to_tf = model.get_weights(transcription_factor.W)
sns.clustermap(kinase_to_tf)


# Now lets look at the TF to Gene relationship. this is more interseting, although underwhelming. It seems that there are some clusters but the dimensionality is way lower than expected. maybe there are too many tfs?
#
#

# In[ ]:

get_ipython().magic('matplotlib inline')
tf_to_gene = model.get_weights(gene.W)
sns.clustermap(tf_to_gene)


# In[ ]:

model.save('mymodel_relu_cocktail_without_val.tflearn')


# In[ ]:

model.load('mymodel.tflearn')


# In[ ]:

pred = model.predict(data)


# In[ ]:

get_ipython().magic('matplotlib inline')
import matplotlib.pyplot
import pylab

for n in range(1,301):
    f, ax = matplotlib.pyplot.subplots()
    ax.scatter(target[n],pred[n])

    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
        np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
    ]

    # now plot both limits against eachother
    ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
    ax.set_aspect('equal')
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_title(meta["train"][n])
    matplotlib.pyplot.show()


# In[ ]:

meta["Sample_Name"][5]


# In[ ]:

get_ipython().magic('matplotlib inline')
sns.clustermap(target)


# In[ ]:

get_ipython().magic('matplotlib inline')
sns.clustermap(pred)


# In[ ]:

get_ipython().magic('matplotlib inline')
matplotlib.pyplot.scatter(
    np.array(pred).min(axis=0),
    target.min(axis=0)
)


# In[ ]:

get_ipython().magic('matplotlib inline')
w = model.get_weights(gene.b)
matplotlib.pyplot.hist(w)


# In[ ]:

np.sum(np.abs(w[1]))


# In[ ]:

data.shape


# In[ ]:
