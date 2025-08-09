from matplotlib.ticker import MaxNLocator
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--period', type=float)
parser.add_argument('-d', '--dc', type=float)
parser.add_argument('-a', '--amplitude', type=float)
parser.add_argument('-p', '--phase', type=float)
parser.add_argument('-e', '--tolerance', type=float)
parser.add_argument('-f', '--function', type=str)
args = parser.parse_args()
t0 = args.period
dc = args.dc
amp = args.amplitude
phi = args.phase
func = args.function
tol = args.tolerance

sine = lambda t: dc + amp * np.sin(2.0 * np.pi * t / t0 + phi)

def square(t):
    rad = np.floor(2.0 * t / t0 + phi / np.pi)
    return np.where(np.mod(rad, 2), dc - amp, dc + amp)

def sawtooth(t):
    rad = t / t0 + phi / np.pi * 0.5
    return dc + amp * 2.0 * (rad - np.floor(rad + 0.5))

def triangle(t):
    aux = t - 0.25 * t0 + phi * t0 * 0.5 / np.pi
    return dc + amp * (4.0 * np.abs( (0.5 * t0 - np.mod(aux, t0)) / t0) - 1.0)

os.chdir(f"../build/{func}")

exp = pd.read_csv("data/exp.dat")
trig = pd.read_csv("data/trig.dat")
hist = pd.read_csv("data/hist.dat")

max_frames = min(exp.shape[0], 51)

n = np.arange(1, exp.shape[0] + 1)
t = np.linspace(0.0, 3.0 * t0, 300000, endpoint=False) # [s]
    
# For creating a dense representation of the approximations
tt, nn = np.meshgrid(t, n[:max_frames-1]) # [:-1] since it is the number of harmonics and not coefficients (remove DC)

# Exponential approximation as a function of the number of harmonics
exp_matrix = 2.0 * np.real((exp.Real[1:max_frames] + 1.j * exp.Imag[1:max_frames]).values.reshape(-1, 1) * np.exp(2.j * np.pi * nn * tt / t0))
exp_matrix = exp_matrix.cumsum(axis=0)
def exp_approx(k):
    if k == 0:
        return np.zeros_like(t)
    elif k== 1: 
        return np.full((t.size,), exp.Real[0])
    elif k > 1: 
        return exp.Real[0] + exp_matrix[k - 2]
    else:
        raise ValueError

# Trigonometric approximation as a function of the number of harmonics
trig_matrix = trig.Amplitude[1:max_frames].values.reshape(-1, 1) * np.cos(2. * np.pi * nn * tt / t0 + trig.Phase[1:max_frames].values.reshape(-1, 1))
trig_matrix = trig_matrix.cumsum(axis=0)
def trig_approx(k):
    if k == 0:
        return np.zeros_like(t)
    elif k == 1:
        return np.full((t.size,), trig.Amplitude[0])
    elif k > 1:
        return trig.Amplitude[0] + trig_matrix[k - 2]
    else:
        raise ValueError

# Original signal
if func == "SINE":
    signal = sine(t)
elif func == "SQUARE":
    signal = square(t)
elif func == "SAWTOOTH":
    signal = sawtooth(t)
elif func == "TRIANGLE":
    signal = triangle(t)
else:
    pass

# Plot MSE and RMSE history
fig1, ax1 = plt.subplots(2, 1, figsize=(8, 4))
xmarg = (n.max() - n.min()) * ax1[-1].margins()[0]
xlim = [n.min() - xmarg, n.max() + xmarg]
for i, l in enumerate(hist.columns): # MSE and RMSE in different columns
    ax1[i].xaxis.set_major_locator(MaxNLocator(integer=True))
    ax1[i].plot(n, hist[l], label=l)
    ax1[i].set_ylabel(l)
    ax1[i].set_xlim(xlim)
    # ax1[i].set_xticks(n[::np.where(n.size > 5, 5, 1)])
ax1[-1].hlines([tol], xlim[0], xlim[1], lw=1.0, ls="-.", color="k")
ax1[-1].set_xlabel("# of Coefficients")
fig1.suptitle("Loss x # of Coefficients")
fig1.legend(loc="upper right")
fig1.tight_layout()
fig1.savefig("figures/loss.png")
fig1.show()

# Plot amplitude of the harmonics
fig2, ax2 = plt.subplots(figsize=(9, 3))
ax2.xaxis.set_major_locator(MaxNLocator(integer=True))
ax2.stem(n[1:] - 1, np.where((trig.Phase[1:] > 0.5 * np.pi) & (trig.Phase[1:] < 1.5 * np.pi), 
                             -trig.Amplitude[1:], 
                              trig.Amplitude[1:]))
ax2.set_xlabel("n")
ax2.set_ylabel("Amplitude")
if n.size == 2: ax2.set_xticks([1])
fig2.suptitle(f"Amplitude of the n-th {t0**-1:.2g}[Hz] Harmonic")
fig2.tight_layout()
fig2.savefig("figures/amplitude.png")
fig2.show()

# Plot approximation and squared error
fig3, axes3 = plt.subplots(2, 1, figsize=(10, 5))
ax3, ax4 = axes3
title_str = "Original Signal x Fourier Series\n# of Harmonics: {:0{}}\nMSE = {:.7f}\nRMSE = {:.7f}"
title3 = fig3.suptitle(title_str)
# Approximation
ax3.plot(t, signal, label="Original signal")
fs_line = ax3.plot(t, np.zeros_like(t), label="Fourier Series Approximation")[0]
ax3.set_ylabel("Amplitude")
ax3.set_ylim([dc - amp * 1.3, dc + amp * 1.3])
ax3.yaxis.set_inverted(False)
# Squared Error
se = (signal - trig_approx(1))**2
mse_line = ax4.plot(t, se, label="Squared error", c="C2")[0]
ax4.set_ylim([-ax4.margins()[1], se.max() + ax4.margins()[1]])
ax4.set_xlabel("Time [s]")
ax4.set_ylabel("Squared Error")

zpad = min(int(np.log10(n.max())) + 1, 2) # Zero-pad string
def update_frame(frame):
    y = trig_approx(frame + 1)
    se = (signal - y)**2
    mse, rmse = hist.MSE[frame], hist.RMSE[frame]

    fs_line.set_ydata(y)
    mse_line.set_ydata(se)
    title3.set_text(title_str.format(frame, zpad, mse, rmse))
    return (fs_line, mse_line, title3)

# Only animate the 50 first harmonics.
fps = 5
fig3.legend(loc="upper right")
fig3.tight_layout()
ani3 = animation.FuncAnimation(fig=fig3, func=update_frame, frames=max_frames, interval=1000 / fps, repeat_delay=250)
ani3.save("figures/approximation.gif", writer=animation.writers["pillow"](fps=fps))
fig3.savefig("figures/func.png")
plt.show()