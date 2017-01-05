import tensorflow as tf
import numpy as np
import numpy as np
import pandas as pd
data = pd.read_csv('../../input/nn/data.csv')
del data["Unnamed: 0"] # This holds the gene names just so you know...
meta = pd.read_csv('../../input/nn/meta.csv')
del meta["Unnamed: 0"]
stress_one_hot = pd.get_dummies(meta["Stress"])
stress_one_hot.shape
data = data.transpose()
data.shape
import tensorflow.contrib.learn.python.learn as learn
from sklearn import datasets, metrics, preprocessing
x = preprocessing.StandardScaler().fit_transform(stress_one_hot)
feature_columns = learn.infer_real_valued_columns_from_input(x)
regressor = learn.LinearRegressor(feature_columns=feature_columns)
regressor.fit(x, data[1], steps=200, batch_size=32)
predictions = list(regressor.predict(x, as_iterable=True))
score = metrics.mean_squared_error(predictions, data[1])
print ("MSE: %f" % score)
from sklearn import linear_model
clf = linear_model.MultiTaskLasso(alpha=0.1)
clf.fit(x,data)
print(clf.coef_)
c = (clf.coef_)
c.shape
import matplotlib
matplotlib.use('GTK')
import seaborn as sns

ax = sns.heatmap(c)
