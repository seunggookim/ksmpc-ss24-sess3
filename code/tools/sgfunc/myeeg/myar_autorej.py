# RUN autoreject
# REF: https://autoreject.github.io/auto_examples/plot_global_reject.html#sphx-glr-auto-examples-plot-global-reject-py


import argparse
parser = argparse.ArgumentParser(description='Creates cross-validation curve')
reqoptions = parser.add_argument_group('Required arguments')
reqoptions.add_argument('-i', '-in', dest="fn_epochs", required=True,
                        help='.SET filename for EPOCHs from EEGLAB')
args = parser.parse_args()

## READ data:
import mne
epochs = mne.io.read_epochs_eeglab(args.fn_epochs)

## RUN auto-reject:
from autoreject import (AutoReject, set_matplotlib_defaults)
import numpy as np
n_interpolates = np.array([1, 4, 32])
consensus_percs = np.linspace(0, 1.0, 11)
ar = AutoReject(n_interpolates, consensus_percs, thresh_method='random_search',
                verbose=False)
ar.fit(epochs)

## EXPORT metainfo:
loss = ar.loss_['eeg'].mean(axis=-1)
idx, jdx = np.unravel_index(loss.argmin(), loss.shape)
consensus_perc_opt = consensus_percs[idx]
n_interp_opt = n_interpolates[jdx]
log = ar.get_reject_log(epochs)
#log.labels       # 0=good 1=bad/uninterp, 2=bad/interp
#log.bad_epochs   # logical 1-D array True=BAD, False=GOOD
from scipy.io import savemat
import os
savemat(os.path.splitext(args.fn_epochs)[0] + "_ar.mat", dict(
    labels_epochs_by_channels=log.labels, bad_epochs=log.bad_epochs,
    consensus_perc_opt=consensus_perc_opt, n_interp_opt=n_interp_opt))


## SAVE cleaned epochs:
epochs_clean = ar.transform(epochs)
epochs_clean.save(os.path.splitext(args.fn_epochs)[0] + "_clean_epo.fif",
                  overwrite=True)

## PLOT
evoked = epochs.average()
evoked_clean = epochs_clean.average()
for badchan in ['LHEOG', 'RHEOG', 'LVEOG']:
    evoked.info['bads'] = [badchan]
    evoked_clean.info['bads'] = [badchan]
    
    import matplotlib.pyplot as plt  # noqa
    set_matplotlib_defaults(plt)
    
    fig, axes = plt.subplots(2, 1, figsize=(6, 6))
    for ax in axes:
        ax.tick_params(axis='x', which='both', bottom='off', top='off')
        ax.tick_params(axis='y', which='both', left='off', right='off')
    
    ylim = dict(grad=(-170, 200))
    evoked.plot(exclude=[], axes=axes[0], ylim=ylim, show=False)
    axes[0].set_title('Before autoreject: ' + badchan, fontsize=11)
    evoked_clean.plot(exclude=[], axes=axes[1], ylim=ylim, show=False)
    axes[1].set_title('After autoreject: ' + badchan, fontsize=11)
    plt.tight_layout()
    plt.draw()
    plt.pause(0.001)
    fig.savefig(
        os.path.splitext(args.fn_epochs)[0] + "_autorej_" + badchan + ".png")
    plt.close(fig)

