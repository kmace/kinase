import numpy as np
import pandas as pd
# In[38]:
expression = pd.read_csv('../../input/nn/data.csv')
genes = expression["Unnamed: 0"]
del expression["Unnamed: 0"] # This holds the gene names just so you know...
meta = pd.read_csv('../../input/nn/meta.csv')
experiments = meta["Unnamed: 0"]
del meta["Unnamed: 0"]
# In[39]:
# expression
stress_one_hot = pd.get_dummies(meta["Stress"])
stress_one_hot.shape
# In[40]:
expression = expression.transpose()
expression.shape

import tflearn

target = expression
target = np.array(expression, dtype='float32')
data = stress_one_hot
data = np.array(data, dtype='float32')


num_stresses = 9
num_kinase = 32
num_transcription_factors = 100
num_genes = 6692

# Build neural network
stress = tflearn.input_data(shape=[None, num_stresses])
kinase = tflearn.fully_connected(stress, num_kinase)
transcription_factor = tflearn.fully_connected(kinase, num_transcription_factors)
gene = tflearn.fully_connected(transcription_factor, num_genes)

adam = tflearn.Adam(learning_rate=0.001, beta1=0.99)

regression = tflearn.regression(gene, optimizer=adam, loss='mean_square')


# Define model
model = tflearn.DNN(regression, tensorboard_verbose=3)
# Start training (apply gradient descent algorithm)
model.fit(data, target, n_epoch=1000, batch_size=441, show_metric=True)

df = pd.DataFrame()
df['nn_beta'] = model.get_weights(gene.b)
df['gene_mean'] = np.mean(expression,0)
fig1 = plt.figure()
sns.lmplot('nn_beta','gene_mean',data=df)
fig1.show()

stress_to_kinase = model.get_weights(kinase.W)
fig2 = plt.figure()
sns.clustermap(stress_to_kinase)
fig2.show()

kinase_to_tf = model.get_weights(transcription_factor.W)
fig3 = plt.figure()
sns.clustermap(kinase_to_tf)
fig3.show()

tf_to_gene = model.get_weights(gene.W)
fig4 = plt.figure()
sns.clustermap(tf_to_gene)
fig4.show()
