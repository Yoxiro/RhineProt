
import numpy as np
from RhineAMP.utils import read_code_ml, save_file, draw_plot, calculate_prediction_metrics
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import StratifiedKFold, GridSearchCV
import re


def MLP_Classifier(X,
                   y,
                   fold=5,
                   hidden_layer_size=(32, 32),
                   lost='lbfgs',
                   activation='relu',
                   lr=0.001,
                   epochs=200):
    classes = sorted(list(set(y)))
    folds = StratifiedKFold(fold).split(X, y)
    prediction_result_cv = []

    for trained, valided in folds:
        train_y, train_X = y[trained], X[trained]
        valid_y, valid_X = y[valided], X[valided]
        model = MLPClassifier(activation=activation,
                              alpha=1e-05,
                              batch_size='auto',
                              beta_1=0.9,
                              beta_2=0.999,
                              early_stopping=False,
                              epsilon=1e-08,
                              hidden_layer_sizes=hidden_layer_size,
                              learning_rate='constant',
                              learning_rate_init=lr,
                              max_iter=epochs,
                              momentum=0.9,
                              nesterovs_momentum=True,
                              power_t=0.5,
                              random_state=1,
                              shuffle=True,
                              solver=lost,
                              tol=0.0001,
                              validation_fraction=0.1,
                              verbose=False,
                              warm_start=False)
        model.fit(train_X, train_y)
        scores = model.predict_proba(valid_X)
        tmp_result = np.zeros((len(valid_y), len(classes) + 1))
        tmp_result[:, 0], tmp_result[:, 1:] = valid_y, scores
        prediction_result_cv.append(tmp_result)
    header = 'Hidden layer size: ' + str(hidden_layer_size)
    return header, prediction_result_cv

def MLP(data,
        label,
        hidden = None,
        lost:str="lbfgs",
        activation:str = "relu",
        epochs:int = 20000,
        lr:float = 0.0001,
        fold:int = 5,
        out:str= "MLP_output"):
    if hidden is None:
        hidden = "32:32"
    X,y = data,label
    hidden_layer_size = eval('(' + re.sub(':', ',', hidden) + ')')
    para_info, cv_res= MLP_Classifier(X,
                                      y,
                                      fold=fold,
                                      epochs=epochs,
                                      hidden_layer_size=hidden_layer_size,
                                      lost=lost,
                                      activation=activation,
                                      lr=lr)
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
