
# coding: utf-8

# what if i were to encode a filtering matrix, to only allow TFs to act on particular genes, from yeastract?

# In[2]:

import numpy as np
import pandas as pd
import tflearn
import tensorflow as tf

meta = pd.read_csv('../../input/nn/meta.csv')
experiments = meta["Unnamed: 0"]
del meta["Unnamed: 0"]

stress_one_hot = pd.get_dummies(meta["Condition"])

kinase_deletion = pd.get_dummies(meta["Strain"])
kinase_one_hot = 1 - kinase_deletion

regulation = pd.read_csv('../../input/nn/RegulationMatrix_Documented_2013927.csv', sep=';')
regulation_tf = regulation["6732"]
del regulation["6732"]
reg_genes = regulation.columns
reg_genes = pd.np.array(reg_genes)

expression = pd.read_csv('../../input/nn/data.csv')
exp_genes = expression["Unnamed: 0"]
exp_genes = pd.np.array(exp_genes)
del expression["Unnamed: 0"] # This holds the gene names just so you know...
expression = expression.transpose()
expression.columns = exp_genes

print('Reg Shape: {}'.format(regulation.shape))
print('Exp Shape: {}'.format(expression.shape))


common_genes = [val for val in exp_genes if val in reg_genes]
print('Reg Genes: {}'.format(len(reg_genes)))
print('Exp Genes: {}'.format(len(exp_genes)))
print('Common Genes: {}'.format(len(common_genes)))

expression = expression[common_genes]
regulation = regulation[common_genes]

print('Reg Shape: {}'.format(regulation.shape))
print('Exp Shape: {}'.format(expression.shape))


target = expression
target = np.array(expression, dtype='float32')
target_mean = target.mean(axis=0, keepdims=True) # I had this wrong before, the time it worked... weird axis was 1
target_std = target.std(axis=0, keepdims=True)
target = target - target_mean
target = target / target_std

data1 = stress_one_hot
data1 = np.array(data1, dtype='float32')
data_mean = data1.mean(axis=0, keepdims=True)
data_std = data1.std(axis=0, keepdims=True)
data1 = data1 - data_mean
data1 = data1 / data_std

data2 = kinase_one_hot
data2 = np.array(data2, dtype='float32')

data3 = regulation
data3 = np.array(data3, dtype='float32')


# In[3]:

data1.shape
#(301, 10)
data2.shape
#(301, 29)


# In[4]:

num_stresses = 10
num_kinase = 29
num_transcription_factors = 313
num_genes = 6223

# Build neural network
# Input variables (10)
# Which Node to dropout (32)
stress = tflearn.input_data(shape=[None, num_stresses], name='Stress')
kinase_deletion = tflearn.input_data(shape=[None, num_kinase], name='Inactive_Kinase')
tf_connectivity = tf.constant(data3, shape=[num_transcription_factors, num_genes], name='TF_to_Gene_Connections')

# This is the layer that I want to perform selective dropout on,
# I should be able to specify which of the 29 nodes should output zero
# based on a 1X29 vector of ones and zeros.
kinase = tflearn.fully_connected(stress, num_kinase, activation='relu', name='Kinase')
with tf.name_scope('Active_Kinases'):
    kinase_dropout = tf.mul(kinase, kinase_deletion)

transcription_factor = tflearn.fully_connected(kinase_dropout, num_transcription_factors, activation='relu', name='Transcription_Factor')

with tf.name_scope('TF_Targets'):
    W = tf.Variable(tf.random_uniform([num_transcription_factors, num_genes]), name='Weights')
    b = tf.Variable(tf.zeros([num_genes]), name='biases')
#gene = tf.matmul( tf.mul(tf_connectivity, W), transcription_factor) + b
    transcription_factor_targets = tf.mul(tf_connectivity, W)

with tf.name_scope('Gene'):
    gene = tf.matmul(transcription_factor, transcription_factor_targets) + b

#gene = tflearn.fully_connected(transcription_factor, num_genes, activation='linear')

adam = tflearn.Adam(learning_rate=0.00001, beta1=0.99)

regression = tflearn.regression(gene, optimizer=adam, loss='mean_square', metric='R2')

# Define model
model = tflearn.DNN(regression, tensorboard_verbose=1)


# In[ ]:

# Start training (apply gradient descent algorithm)
#model.fit([data1, data2], target, n_epoch=200, show_metric=True, shuffle=True, validation_set=0.10)
model.fit([data1, data2], target, n_epoch=200, show_metric=True, shuffle=True, batch_size=20, validation_set=0.10)
model.save('mymodel_relu_cocktail_with_val_and_dropout_and_reg_test.tflearn')


# In[ ]:

tf_filter


# In[ ]:

transcription_factor


# In[ ]:

transcription_factor * tf_filter


# In[ ]:

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
stress_to_kinase.columns = kinase_one_hot.columns
test = np.abs(stress_to_kinase.T)
test.columns = stress_one_hot.columns
cg = sns.clustermap(test, metric="correlation")
plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
plt.show()


# In[ ]:

b = model.get_weights(kinase.b)
b


# Let us now plot the weight matrix between kinase outputs and transcription factors, each column is a TF (100) and each row is a kinase (32)

# In[ ]:

get_ipython().magic('matplotlib inline')
kinase_to_tf = pd.DataFrame(model.get_weights(transcription_factor.W))
kinase_to_tf.columns = regulation_tf
kinase_to_tf = kinase_to_tf.T
kinase_to_tf.columns = kinase_one_hot.columns
kinase_to_tf = np.abs(kinase_to_tf)
cg = sns.clustermap(kinase_to_tf)
plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
kinase_to_tf.shape


# In[ ]:




# Now lets look at the TF to Gene relationship. this is more interseting, although underwhelming. It seems that there are some clusters but the dimensionality is way lower than expected. maybe there are too many tfs?
#
#

# In[ ]:

get_ipython().magic('matplotlib inline')
tf_to_gene = model.get_weights(gene.W)
sns.clustermap(tf_to_gene)


# In[ ]:

model.save('mymodel_relu_cocktail_without_val_and_dropout.tflearn')


# In[ ]:

model.load('mymodel.tflearn')


# In[ ]:

pred = model.predict([data1, data2])


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
    ax.set_title(meta["Condition"][n])
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

kinase_to_tf = pd.DataFrame(model.get_weights(transcription_factor.W))
kinase_to_tf.shape


# In[ ]:

regulation = pd.read_csv('../../input/nn/RegulationMatrix_Documented_2013927.csv', sep=';')
regulation_tf = regulation["6732"]
del regulation["6732"]
genes = regulation.columns
# this is problematic becuase the gene names are used here instead of the yorf names


# In[ ]:

w = model.get_weights(W)
data3
wf = (w*data3)


# In[ ]:

get_ipython().magic('matplotlib inline')
w = model.get_weights(W)
data3
wf = (w*data3)
wf = np.abs(wf)
sns.heatmap(wf)
#plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
wf.shape


# In[ ]:

1+1


# In[ ]:
