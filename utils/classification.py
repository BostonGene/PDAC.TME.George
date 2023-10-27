import copy

import pandas as pd
from sklearn.neighbors import KNeighborsClassifier

from utils.utils import median_scale


class KNeighborsClusterClassifier(object):
    """
    Clustering using KNeighbors algorithm.
    """

    def __init__(self, norm=True, algorithm='auto', clip=None, scale=False, k=25):
        self.norm = norm
        self.median = 0
        self.mad = 1
        self.X = None
        self.y = None
        self.algorithm = algorithm
        self.model = None
        self.clip = clip
        self.scale = scale
        self.k = k

    # some workaround pickle.load segmentation fault
    def __getstate__(self):
        self.model = None
        return self.__dict__

    def __setstate__(self, state):
        for attr, val in state.items():
            self.__setattr__(attr, val)
        self.model = KNeighborsClassifier(algorithm=self.algorithm, n_neighbors=self.k).fit(self.X, self.y)

    def check_is_fitted(self):
        return (self.X is not None) and (self.y is not None) and (self.model is not None)

    def preprocess_data(self, X):
        """
        Preprocessing data

        :param X: pd.DataFrame, index - samples, columns - features
        """
        x = copy.deepcopy(X)
        if self.scale:
            x = median_scale(x)

        x = (x - self.median) / self.mad
        if self.clip is not None and self.clip > 0:
            x = x.clip(-1 * self.clip, self.clip)

        return x

    def check_columns(self, X):
        """
        Checking columns for X (matching with self.X.columns)
        """
        if hasattr(self.X, 'columns'):
            try:
                return X[self.X.columns]
            except KeyError:

                raise Exception('Columns do not match')
        return X

    def check_rows(self, Y):
        """
        Checking indexes for Y (matching with self.X.index and vice versa)
        """
        if hasattr(self.X, 'index'):
            try:
                return Y.loc[self.X.index]
            except KeyError:
                raise Exception('Indexes do not match')
        return Y

    def fit(self, X, y):
        """
        Fit model - calculate centroids

        :param X: pd.DataFrame, RNA data, columns - features, index - samples
        :param y: array-like, cluster labels
        """
        if X.shape[0] != len(y):
            raise Exception('Indexes do not match')

        if self.norm:
            self.median = X.median()
            self.mad = X.mad()

        self.X = self.preprocess_data(X)
        self.y = self.check_rows(Y=copy.deepcopy(y))

        self.model = KNeighborsClassifier(algorithm=self.algorithm, n_neighbors=self.k).fit(self.X, self.y)

        return self

    def predict(self, X):
        """
        Predict - return cluster labels

        :param X: pd.DataFrame, RNA data, columns - features, index - samples
        :return: pd.Series, predicted cluster labels
        """
        if X.shape[1] != self.X.shape[1]:
            raise Exception('Сolumns do not match')

        x_scaled = self.preprocess_data(self.check_columns(X))
        # Here self.model.predict is used in order to mimic its' way to select the
        #  class in case of equal probabilities
        return pd.Series(self.model.predict(x_scaled), index=x_scaled.index)

    def predict_proba(self, X):
        """
        :param X: pd.DataFrame, RNA data, columns - features, index - samples
        :return: pd.DataFrame, probabilities for each cluster. Index - samples, columns - clusters
        """

        if X.shape[1] != self.X.shape[1]:
            raise Exception('Columns do not match')

        x_scaled = self.preprocess_data(self.check_columns(X))

        return pd.DataFrame(self.model.predict_proba(x_scaled).astype(float), index=x_scaled.index,
                            columns=self.model.classes_)