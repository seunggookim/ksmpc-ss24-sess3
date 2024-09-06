# Plot CV curve
# REF: https://autoreject.github.io/auto_examples/plot_global_reject.html#sphx-glr-auto-examples-plot-global-reject-py


import argparse
parser = argparse.ArgumentParser(description='Creates cross-validation curve')
reqoptions = parser.add_argument_group('Required arguments')
reqoptions.add_argument('-i', '-in', dest="fn_epochs", required=True,
                        help='.SET filename for EPOCHs from EEGLAB')
args = parser.parse_args()

## DEFINE figure filename
#import os
#path = os.path.dirname(args.fn_epochs)
#basename = os.path.basename(args.fn_epochs)
#path_png = path + "/mne.autoreject"
#if not os.path.isdir(path_png):
#    os.mkdir(path_png)
#fn_png = path_png + "/" + basename + "_CVcurv.png"
fn_png = args.fn_epochs + "_CVcurv.png"

## READ data:
import mne
epochs = mne.io.read_epochs_eeglab(args.fn_epochs)

## Plot learning curve
import numpy as np  # noqa
param_range = np.linspace(40e-6, 500e-6, 30)

from autoreject import validation_curve  # noqa
from autoreject import get_rejection_threshold  # noqa

_, test_scores, param_range = validation_curve(
    epochs, param_range=param_range, cv=5, return_param_range=True, n_jobs=1)

test_scores = -test_scores.mean(axis=1)
best_thresh = param_range[np.argmin(test_scores)]

reject2 = get_rejection_threshold(epochs, random_state=None, cv=5)

import matplotlib.pyplot as plt  # noqa
from autoreject import set_matplotlib_defaults  # noqa
set_matplotlib_defaults(plt)


#human_thresh = 80e-6  # this is a threshold determined visually by a human

unit = r'$\mu$V'
scaling = 1e6
plt.figure(figsize=(8, 5))
plt.tick_params(axis='x', which='both', bottom='off', top='off')
plt.tick_params(axis='y', which='both', left='off', right='off')
colors = ['#E24A33', '#348ABD', '#988ED5', 'k']
plt.plot(scaling * param_range, scaling * test_scores,
         'o-', markerfacecolor='w',
         color=colors[0], markeredgewidth=2, linewidth=2,
         markeredgecolor=colors[0], markersize=8, label='CV scores')
plt.ylabel('RMSE (%s)' % unit)
plt.xlabel('Threshold (%s)' % unit)
plt.xlim((scaling * param_range[0] * 0.9, scaling * param_range[-1] * 1.1))
plt.axvline(scaling * best_thresh, label='auto global', color=colors[2],
            linewidth=2, linestyle='--')
plt.axvline(scaling * reject2['eeg'], label='bayes opt', color=colors[3],
            linewidth=2, linestyle='--')
# plt.axvline(scaling * human_thresh, label='manual', color=colors[1],
#             linewidth=2, linestyle=':')
plt.legend(loc='upper right')
plt.title(args.fn_epochs, fontsize=11)
plt.tight_layout()
plt.draw()
plt.pause(0.001)
plt.savefig(fn_png)
plt.close()