
import numpy
import numpy as np
from RhineAMP.utils import save_file, draw_plot, calculate_prediction_metrics
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import StratifiedKFold


def KNN_Classifier_binary(X,
                          y,
                          fold=5,
                          n_neighbors=3):
    classes = sorted(list(set(y)))
    prediction_result_cv = []

    folds = StratifiedKFold(fold).split(X, y)

    for i, (trained, valided) in enumerate(folds):
        train_y, train_X = y[trained], X[trained]
        valid_y, valid_X = y[valided], X[valided]
        model = KNeighborsClassifier(n_neighbors=n_neighbors).fit(train_X, train_y)
        scores = model.predict_proba(valid_X)
        tmp_result = np.zeros((len(valid_y), len(classes) + 1))
        tmp_result[:, 0], tmp_result[:, 1:] = valid_y, scores
        prediction_result_cv.append(tmp_result)

        # independent
    header = 'n_neighbors: %d\n' % n_neighbors
    return header, prediction_result_cv


def KNN(data: numpy.ndarray,
        label: list or numpy.ndarray,
        k: int = 3,
        fold: int = 5,
        out: str = "KNN_output"):
    X, y = data, label

    para_info, cv_res = KNN_Classifier_binary(X,
                                              y,
                                              fold=fold,
                                              n_neighbors=k)

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
