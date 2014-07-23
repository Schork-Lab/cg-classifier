'''
Created on Jul 23, 2014

@author: Kunal Bhutani
'''

import cPickle
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression as LogisticClassifier

from processing.variants import VariantFile

class Classifier(object):
    '''
    Classifier class that handles wraps the scikit learn models. 
    Bare bones currently, but will be nice once we try to calculate more complex
    statistics and integrate with figures and other potential uses.
    
    Also again, abstracted into super class and subclasses to handle
    different classifiers.
    
    '''
    def __init__(self, clf=None, from_pickle=False):
        if from_pickle:
            self.clf = cPickle.load(open(clf))
        else:
            self.clf = clf

    def fit(self, VariantFile, truth):
        features = VariantFile.features
        self.clf.fit(features, truth)
        return

    def predict(self, VariantFile):
        features = VariantFile.features
        return self.clf.predict(features)

    def score(self, VariantFile, truth):
        features = VariantFile.features
        return self.clf.score(features, truth)

class RandomForest(Classifier):

    def __init__(self, default=True, *args, **kwargs):
        if default:
            clf = RandomForestClassifier(n_estimators=500,
                                         n_jobs=8,
                                         oob_score=True,
                                         criterion='entropy',
                                         min_samples_split=4,
                                         min_samples_leaf=2)
        else:
            clf = RandomForestClassifier(*args, **kwargs)

        super(RandomForest, self).__init__(clf)


class LogisticRegression(Classifier):

    def __init__(self, default=True, *args, **kwargs):
        if default:
            clf = LogisticClassifier(penalty='l2', C=1)
        else:
            clf = LogisticClassifier(*args, **kwargs)

        super(LogisticRegression, self).__init__(clf)
