import math

import numpy as np
import pandas
from sklearn import svm
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from RhineAMP.utils import  save_file, draw_plot, calculate_prediction_metrics

def SVM_Classifier(X,
                   y,
                   indep=None,
                   fold=5,
                   batch=None,
                   auto=False,
                   kernel='rbf',
                   degree=3,
                   gamma='auto',
                   coef0=0.0,
                   C=1.0):
    default_params = {'degree': degree,
                      'gamma': gamma,
                      'coef0': coef0,
                      'C': C}
    if auto:
        data = np.zeros((X.shape[0], X.shape[1] + 1))
        data[:, 0] = y
        data[:, 1:] = X
        np.random.shuffle(data)
        X1 = data[:, 1:]
        y1 = data[:, 0]
        parameters = {'kernel': ['linear'],
                      'C': [1, 15]} \
            if kernel == 'linear' \
            else {'kernel': [kernel], 'C': [1, 15], 'gamma': 2.0 ** np.arange(-10, 4)}
        optimizer = GridSearchCV(svm.SVC(probability=True), parameters)
        optimizer = optimizer.fit(X1[0:math.ceil(batch * X1.shape[0]), ],
                                  y1[0:math.ceil(batch * y1.shape[0]), ]) if batch else optimizer.fit(X, y)
        params = optimizer.best_params_
        default_params['C'] = params['C']
        if kernel != 'linear':
            default_params['gamma'] = params['gamma']

    classes = sorted(list(set(y)))
    svms = []
    cvs = np.zeros((X.shape[0], len(classes) + 1))
    folds = StratifiedKFold(fold).split(X, y)
    # todo ind part
    prediction_result_cv = []
    for trained, valided in folds:
        train_y, train_X = y[trained], X[trained]
        valid_y, valid_X = y[valided], X[valided]
        model = svm.SVC(C=default_params['C'],
                        kernel=kernel,
                        degree=default_params['degree'],
                        gamma=default_params['gamma'],
                        coef0=default_params['coef0'],
                        probability=True,
                        random_state=1)
        svc = model.fit(train_X, train_y)
        svms.append(svc)
        proba_ = svc.predict_proba(valid_X)
        cvs[valided, 0], cvs[valided, 1:] = valid_y, proba_
        # save the sample label and prediction result to ndarray
        tmp_result = np.zeros((len(valid_y), len(classes) + 1))
        tmp_result[:, 0], tmp_result[:, 1:] = valid_y, proba_
        prediction_result_cv.append(tmp_result)

        # todo inde
        # # independent
        # if indep.shape[0] != 0:
        #     inds[:, 1:] += svc.predict_proba(indep[:, 1:])

    header = 'C=%f\tgamma=%s' % (default_params['C'], default_params['gamma'])
    # todo

    # if indep.shape[0] != 0:
    #     inds[:, 1:] /= fold
    return header, prediction_result_cv


def SVM(X,
        y,
        kernel: str = "rbf",
        auto_opt: bool = True,
        batch: float = 1.0,
        degree: int = 3,
        gamma=None,
        coef: float = 0.0,
        cost: float = 1.0,
        fold: int = 5,
        out:str = "SVM_output"):

    if gamma is None:
        gamma = 'auto'

    independent =np.array([])
    """
    X is the training dataset containing both Positive and Negative DATA
    y is the label of data in X
    """

    # todo 什么是indep？
    # if indep:
    #     ind_X, ind_y = read_code_ml.read_code(args.indep, format='%s' % args.format)
    #     independent = np.zeros((ind_X.shape[0], ind_X.shape[1] + 1))
    #     independent[:, 0], independent[:, 1:] = ind_y, ind_X

    para_info, cv_res = SVM_Classifier(X, y, indep=independent, fold=fold, batch=batch,
                                                auto=auto_opt, kernel=kernel, degree=degree,
                                                gamma=gamma, coef0=coef, C=cost)

    classes = sorted(list(set(y)))
    if len(classes) == 2:
        save_file.save_CV_result_binary(cv_res, '%s_CV.txt' % out, para_info)
        mean_auc = draw_plot.plot_roc_cv(cv_res, '%s_ROC_CV.png' % out, label_column=0, score_column=2)
        mean_auprc = draw_plot.plot_prc_CV(cv_res, '%s_PRC_CV.png' % out, label_column=0, score_column=2)
        cv_metrics = calculate_prediction_metrics.calculate_metrics_cv(cv_res, label_column=0, score_column=2, )
        save_file.save_prediction_metrics_cv(cv_metrics, '%s_metrics_CV.txt' % out)

        # if args.indep:
        #     save_file.save_IND_result_binary(ind_res, '%s_IND.txt' % args.out, para_info)
        #     ind_auc = draw_plot.plot_roc_ind(ind_res, '%s_ROC_IND.png' % args.out, label_column=0, score_column=2)
        #     ind_auprc = draw_plot.plot_prc_ind(ind_res, '%s_PRC_IND.png' % args.out, label_column=0, score_column=2)
        #     ind_metrics = calculate_prediction_metrics.calculate_metrics(ind_res[:, 0], ind_res[:, 2])
        #     save_file.save_prediction_metrics_ind(ind_metrics, '%s_metrics_IND.txt' % args.out)
    if len(classes) > 2:
        save_file.save_CV_result(cv_res, classes, '%s_CV.txt' % out, para_info)
        cv_metrics = calculate_prediction_metrics.calculate_metrics_cv_muti(cv_res, classes, label_column=0)
        save_file.save_prediction_metrics_cv_muti(cv_metrics, classes, '%s_metrics_CV.txt' % out)

        # if args.indep:
        #     save_file.save_IND_result(ind_res, classes, '%s_IND.txt' % args.out, para_info)
        #     ind_metrics = calculate_prediction_metrics.calculate_metrics_ind_muti(ind_res, classes, label_column=0)
        #     save_file.save_prediction_metrics_ind_muti(ind_metrics, classes, '%s_metrics_IND.txt' % args.out)
