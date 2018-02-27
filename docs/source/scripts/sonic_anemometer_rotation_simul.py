import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from fluxer.eddycov.flux import rotate_coordinates

# Basic settings
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams["text.usetex"] = True
mpl.rcParams["text.latex.preamble"] = [r"\usepackage{amsmath}"]


class Arrow3D(FancyArrowPatch):
    """Subclass extending FancyArrowPatch to 3D"""
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        FancyArrowPatch.draw(self, renderer)


def plot_rotation(uvw, theta, rotation_axis, fig_axis):
    """Plot rotation of vectors given angle and rotation axis

    Parameters
    ----------
    vectors : array_like
        An nx3 array of vectors with their `x`, `y`, `z` components
    theta : numeric
        The angle (degrees) by which to perform the rotation.
    rotation_axis : int
        Axis around which to perform the rotation (`x` = 0; `y` = 1; `z` =
        2)
    fig_axis : Axes
        Axes object where plotting is performed

    Returns
    -------
    Axes
    The modified fig_axis object with plot

    """
    uvw_theta = rotate_coordinates(uvw, theta=theta, axis=rotation_axis,
                                   active=False)

    # First set up unit vectors to define the starting coordinate system
    xyz = np.array([[1, 0, 0],
                    [0, 1, 0],
                    [0, 0, 1]])
    origin = np.zeros((3, 3))

    # Rotated coordinate system to place the uvw_theta coordinates; we
    # multiply the starting coordinate system by the transpose of R_z
    xyz_theta = rotate_coordinates(xyz, theta=theta, axis=rotation_axis,
                                   active=True)

    ax = fig_axis
    arrow_opts = dict(mutation_scale=20, color="k", shrinkA=0, shrinkB=0)
    txt_opts = dict(horizontalalignment="center",
                    verticalalignment="center", fontsize=12)
    axlims = [-1, 1]
    axlabels = "xyz"
    title = (r"Rotation around ${0}$-axis by " +
             r"$\theta={1:.1f}^\circ$").format(axlabels[rotation_axis],
                                               theta)
    ax.set_title(title)
    ax.set_xlim3d(axlims)
    ax.set_ylim3d(axlims)
    ax.set_zlim3d(axlims)
    # Below only as guide while developing
    ax.set_xlabel(axlabels[0])
    ax.set_ylabel(axlabels[1])
    ax.set_zlabel(axlabels[2])

    # Draw the original system, providing a matrix of tuples with the 4
    # coordinates for the ends of the axes, to iterate over rows
    xyz = np.zeros((3, 3, 2))
    np.fill_diagonal(xyz[:, :, 1], 1)
    arrow_opts.update({"arrowstyle": "->"})
    for i in np.arange(3):
        a = Arrow3D(xyz[i, 0],
                    xyz[i, 1],
                    xyz[i, 2], **arrow_opts)
        ax.add_artist(a)
        # Add axis labels
        ax.text(1.1 * xyz[i, 0, 1],
                1.1 * xyz[i, 1, 1],
                1.1 * xyz[i, 2, 1],
                r"${}_0$".format(axlabels[i]), **txt_opts)

    # Draw the new coordinate system.  Set up arrows representing the new
    # system
    arrow_opts.update({"color": "0.75", "linewidth": "1.5"})
    # For each axis of the new coordinate system
    for i in np.arange(3):
        # 2x3 matrix with 2 points on the 3D line, passing through the
        # origin, and up to the calculated end of the axis
        ax1crd = np.vstack((origin[i, :], xyz_theta[i, :]))
        # Get the coordinates of the mean(center) point for singular value
        # decomposition
        ax1crd_mean = ax1crd.mean(axis=0)
        _, _, eigv1 = np.linalg.svd(ax1crd - ax1crd_mean)
        # The first row of eigv1 (eigv[0]) is the first principal
        # component, i.e. the direction vector of the line.  We use it to
        # compute the coordinates of 2 points defining a line with known
        # length
        ax1pts = np.transpose(eigv1[0] * np.mgrid[-0.5:0.5:2j][:, np.newaxis] +
                              ax1crd_mean)
        a = Arrow3D(ax1pts[0, :],
                    ax1pts[1, :],
                    ax1pts[2, :], **arrow_opts)
        ax.add_artist(a)
        # Add axis labels
        ax.text(1.1 * ax1pts[0, 1],
                1.1 * ax1pts[1, 1],
                1.1 * ax1pts[2, 1],
                r"${}_1$".format(axlabels[i]), **txt_opts)
        # Add theta angle labels
        arc = np.linspace(0, np.radians(theta))
        arc_offset = 0.5
        # Plot a surface as reference
        xsurf, ysurf = np.meshgrid(np.linspace(-1, 1, 9),
                                   np.linspace(-1, 1, 9))
        zeros = np.zeros((9, 9))
        surf_opts = dict(rstride=1, cstride=1, color="pink", alpha=0.3)
        if i != rotation_axis:
            m = 0.45 * ((xyz[i, :, 1] + ax1pts[:, 1]) / 2.0)
            ax.text(m[0], m[1], m[2], r'$\theta$', **txt_opts)
        if rotation_axis == 0:      # yz plane
            p = np.array([arc * 0, np.cos(arc), np.sin(arc)]) * arc_offset
            ax.plot(p[0, :], p[1, :], p[2, :], ":", c="0.75")
            p = np.array([arc * 0, -np.sin(arc), np.cos(arc)]) * arc_offset
            ax.plot(p[0, :], p[1, :], p[2, :], ":", c="0.75")
            ax.plot_surface(zeros, xsurf, ysurf, **surf_opts)
            ax.view_init(elev=25, azim=-10)
        elif rotation_axis == 1:    # xz plane
            p = np.array([np.sin(arc), arc * 0, np.cos(arc)]) * arc_offset
            ax.plot(p[0, :], p[1, :], p[2, :], ":", c="0.75")
            p = np.array([np.cos(arc), arc * 0, -np.sin(arc)]) * arc_offset
            ax.plot(p[0, :], p[1, :], p[2, :], ":", c="0.75")
            ax.plot_surface(xsurf, zeros, ysurf, **surf_opts)
            ax.view_init(elev=10, azim=50)
        else:                       # xy plane
            p = np.array([np.cos(arc), np.sin(arc), arc * 0]) * arc_offset
            ax.plot(p[0, :], p[1, :], p[2, :], ":", c="0.75")
            p = np.array([-np.sin(arc), np.cos(arc), arc * 0]) * arc_offset
            ax.plot(p[0, :], p[1, :], p[2, :], ":", c="0.75")
            ax.plot_surface(xsurf, ysurf, zeros, **surf_opts)
            ax.view_init(elev=30, azim=-20)

    # Plot each vector
    arrow_opts.update({"arrowstyle": "-|>"})
    vctcols = np.array(["b", "g", "m"])
    for i in np.arange(3):
        # In original coordinate system
        arrow_opts.update({"color": vctcols[i]})
        a = Arrow3D((0, uvw[i, 0]),
                    (0, uvw[i, 1]),
                    (0, uvw[i, 2]), **arrow_opts)
        ax.add_artist(a)
        # Label vector
        txt_opts.update({"verticalalignment": "center"})  # reset
        ax.text(1.1 * uvw[i, 0],
                1.1 * uvw[i, 1],
                1.1 * uvw[i, 2],
                r"$\vec v_{{{0},0}}$".format(i + 1), **txt_opts)
        # In rotated coordinate system
        a = Arrow3D((0, uvw_theta[i, 0]),
                    (0, uvw_theta[i, 1]),
                    (0, uvw_theta[i, 2]), alpha=0.2, **arrow_opts)
        ax.add_artist(a)
        # Label vector
        txt_opts.update({"verticalalignment": "bottom"})
        ax.text(1.1 * uvw_theta[i, 0],
                1.1 * uvw_theta[i, 1],
                1.1 * uvw_theta[i, 2],
                r"$\vec v'_{{{0},1}}$".format(i + 1), **txt_opts)

    ax.set_aspect("equal")
    return ax


# Simple matrix with row vectors x, y, z, representing samples of wind
# vector coordinates u, v, w, as measured by a 3D sonic anemometer.  The
# resulting vectors are lined up with the cartesian (sonic) coordinate
# frame, except that the vector aligned with the z-axis is reversed to make
# it easy to visualize.
uvw = np.array([[0.75, 0, 0],
                [0, 0.75, 0],
                [0, 0, -0.75]])

# Plotting
theta_deg, rotation_axis = 36, 2
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection="3d")
plot_rotation(uvw, theta=theta_deg, rotation_axis=rotation_axis,
              fig_axis=ax)
# fname_fmt = ("../../build/images/savefig/"
#              "coordinate_rotation_{0:.0f}_{1}.png")
# fname = fname_fmt.format(theta_deg, rotation_axis)
# fig.savefig(fname, bbox_inches="tight", frameon=False)
