#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import argparse
import numpy as np
from RhineAMP.utils import save_file, draw_plot, calculate_prediction_metrics
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold,GridSearchCV


def RF_Classifier(X, y, fold=5, n_trees=100):
    """
    Parameters:
    ----------
    :param X: 2-D ndarray
    :param y: 1-D ndarray
    :param indep: 2-D ndarray, the first column is labels and the rest are feature values
    :param fold: int, default 5
    :param n_trees: int, number of trees, default: 5
    :param out:
    :return:
        info: str, the model parameters
        cross-validation result: list with element is ndarray
        independent result: ndarray, the first column is labels and the rest are prediction scores.
    """
    classes = sorted(list(set(y)))

    prediction_result_cv = []

    param_dict = [
        {"n_estimators" : [i for i in range(50,260,10)],
         "max_depth" : [i for i in range(5,25)]}
    ]
    grid_search = GridSearchCV(estimator=RandomForestClassifier(),
                               param_grid = param_dict,
                               cv=5,
                               n_jobs=4)
    grid_search.fit(X,y)
    print("Grid_Search_done")
    print(grid_search.best_params_)
    folds = StratifiedKFold(fold).split(X, y)
    for i, (trained, valided) in enumerate(folds):
        train_y, train_X = y[trained], X[trained]
        valid_y, valid_X = y[valided], X[valided]
        model = RandomForestClassifier(n_estimators=grid_search.best_params_["n_estimators"],
                                       max_depth=grid_search.best_params_["max_depth"],
                                       bootstrap=False)
        rfc = model.fit(train_X, train_y)
        scores = rfc.predict_proba(valid_X)
        tmp_result = np.zeros((len(valid_y), len(classes) + 1))
        tmp_result[:, 0], tmp_result[:, 1:] = valid_y, scores
        prediction_result_cv.append(tmp_result)
        # independent
    header = 'n_trees: %d' % n_trees
    return header, prediction_result_cv
def RF(data,
       label,
       n_trees:int = 100,
       fold:int = 5,
       out:str = "RF_output"):
    X, y = data, label
    para_info, cv_res= RF_Classifier(X, y,fold=fold, n_trees=n_trees)

    classes = sorted(list(set(y)))
    if len(classes) == 2:
        save_file.save_CV_result_binary(cv_res, '%s_CV.txt' % out, para_info)
        mean_auc = draw_plot.plot_roc_cv(cv_res, '%s_ROC_CV.png' % out, label_column=0, score_column=2)
        mean_auprc = draw_plot.plot_prc_CV(cv_res, '%s_PRC_CV.png' % out, label_column=0, score_column=2)
        cv_metrics = calculate_prediction_metrics.calculate_metrics_cv(cv_res, label_column=0, score_column=2, )
        save_file.save_prediction_metrics_cv(cv_metrics, '%s_metrics_CV.txt' % out)

    if len(classes) > 2:
        save_file.save_CV_result(cv_res, classes, '%s_CV.txt' % out, para_info)
        cv_metrics = calculate_prediction_metrics.calculate_metrics_cv_muti(cv_res, classes, label_column=0)
        save_file.save_prediction_metrics_cv_muti(cv_metrics, classes, '%s_metrics_CV.txt' % out)
