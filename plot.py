#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# -----------------------------
# Load data
# -----------------------------
A = np.load("Output/cpoints.npy")   # shape (T, N, 2)
T, N, _ = A.shape

# -----------------------------
# Define particle colors (length N)
# -----------------------------
# Example 1: fixed named colors
# colors = ["red", "blue", "green", "orange", ...]  # must be length N

# Example 2: random RGB colors
colors = N*["black"] 
colors = np.array(colors)

# -----------------------------
# Compute global limits ignoring NaNs
# -----------------------------
xmin = np.nanmin(A[:, :, 0])
xmax = np.nanmax(A[:, :, 0])
ymin = np.nanmin(A[:, :, 1])
ymax = np.nanmax(A[:, :, 1])

# Small manual padding
ymin -= 1
ymax += 1

# 5% automatic padding
pad_frac = 0.05
xrange = xmax - xmin
yrange = ymax - ymin

xmin -= pad_frac * xrange
xmax += pad_frac * xrange
ymin -= pad_frac * yrange
ymax += pad_frac * yrange

# -----------------------------
# Setup figure
# -----------------------------
fig, ax = plt.subplots()
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.set_aspect('equal', adjustable='box')

# -----------------------------
# Initial frame
# -----------------------------
pts0 = A[0]
valid0 = ~np.isnan(pts0).any(axis=1)
pts0 = pts0[valid0]

scat = ax.scatter(
        pts0[:, 0],
        pts0[:, 1],
        c=colors[valid0],
        s=30,
        alpha=1.
        )

# -----------------------------
# Update function
# -----------------------------
def update(frame):
    pts = A[frame]
    valid = ~np.isnan(pts).any(axis=1)
    pts = pts[valid]

    scat.set_offsets(pts)
    scat.set_color(colors[valid])

    return scat,

# -----------------------------
# Create animation
# -----------------------------
ani = FuncAnimation(
        fig,
        update,
        frames=T,
        interval=2000 // T,
        blit=True
        )

ani.save("animation.mp4", writer="ffmpeg", fps=30)


plt.show()
