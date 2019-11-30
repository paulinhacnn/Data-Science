import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
import keras
from keras.datasets import mnist
from keras.models import Sequential
from keras.layers import Dense, Dropout


test_size=0.2
batch_size = 1024
epochs = 45


data = pd.read_csv('genes.csv')
data_array=data.iloc[:,1:].values
X = data_array[:,:1503]
Y = data_array[:,-1]

x_train, x_test, y_train, y_test = train_test_split(X, Y, test_size=test_size)
# modelo
model = Sequential()
model.add(Dense(512, activation='relu', input_shape=(X.shape[1],)))
model.add(Dropout(0.2))
model.add(Dense(512, activation='relu'))
model.add(Dropout(0.2))
model.add(Dense(1, activation='tanh'))
model.summary()

model.summary()

model.compile(loss='mean_squared_error',
              optimizer=keras.optimizers.Adam(),
              metrics=['mae'])

history = model.fit(x_train, y_train,
                    batch_size=batch_size,
                    epochs=epochs,
                    verbose=1,
                    validation_data=(x_test, y_test))
score = model.evaluate(x_test, y_test, verbose=0)
print('Test loss:', score[0])
print('Test accuracy:', score[1])

y_pred=model.predict(x_test)

print(y_test)

np.corrcoef(y_test.flatten(), y_pred.flatten())

# Resultado do modelo de regressão linear

from sklearn.linear_model import LinearRegression

reg = LinearRegression()
reg.fit(x_train, y_train)

np.corrcoef(reg.predict(x_test).flatten(), y_test)

