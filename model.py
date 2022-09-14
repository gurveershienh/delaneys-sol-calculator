import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn import linear_model

df = pd.read_csv('delaney-descriptors.csv')

X = df.drop(['logS'], axis=1)
Y = df['logS']

model = linear_model.LinearRegression()
model.fit(X, Y)

Y_pred = model.predict(X)
Y_series = pd.Series(Y_pred, name='logS_pred')

df = pd.concat( [X, Y, Y_series], axis=1 )

plt.figure(figsize=(5,5))
plt.scatter(x=Y, y=Y_pred, c='blue', alpha=0.4)

z = np.polyfit(Y, Y_pred, 1)
p = np.poly1d(z)

plt.plot(Y, p(Y), 'orange')
plt.ylabel('Predicted LogS')
plt.xlabel('Experimental LogS')
plt.show()

##pickle.dump(model, open('solubility_model.pkl', 'wb'))