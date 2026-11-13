#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Copyright 2015-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

from PyQt5 import QtCore, QtGui, QtWidgets

from PyQt5.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QComboBox

import os
import copy
import math
import matplotlib

import matplotlib.pyplot

from matplotlib.animation import FuncAnimation
from matplotlib.patches import Rectangle

from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar,
)

from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

from collections import OrderedDict
from math import asin, atan2

from .mc_log_types import LineStyle, PlotSide, PlotType

from .mc_log_utils import InitDialogWithOkCancel
from .mc_log_utils import ShowErrorDialog

import numpy as np

import sys

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
# This might fail if mc_log_ui is imported in a headless environment
try:
    matplotlib.use("Qt5Agg")
except:  # noqa: E722
    pass


if sys.version_info[0] > 2:
    unicode = str


class GenerateRangeDialog(QtWidgets.QDialog):
    @InitDialogWithOkCancel(Layout=QtWidgets.QFormLayout, apply_=False)
    def __init__(self, parent):
        self.setWindowTitle("Generate time-range(s) based on data range")
        self.data = parent.data
        self.selectedData = QtWidgets.QComboBox(self)
        keys = list(self.data.keys())
        keys.sort()
        for k in keys:
            self.selectedData.addItem(k)
        self.layout.addRow("Data range", self.selectedData)
        self.stringSelector = QtWidgets.QComboBox(self)
        self.stringSelectorLabel = QtWidgets.QLabel("Select value")
        self.setStringSelector()
        self.layout.addRow(self.stringSelectorLabel, self.stringSelector)
        self.outputName = QtWidgets.QLineEdit()
        self.setDefaultName()
        self.layout.addRow("Name", self.outputName)
        self.selectedData.currentTextChanged.connect(self.setDefaultName)
        self.selectedData.currentTextChanged.connect(self.setStringSelector)
        self.stringSelector.currentTextChanged.connect(self.setDefaultName)

    def setStringSelector(self):
        key = self.selectedData.currentText()
        data = self.data[key]
        self.stringSelector.clear()
        if type(data[0]) is unicode:
            self.stringSelectorLabel.show()
            self.stringSelector.show()
            self.stringSelector.addItem("All")
            self.stringSelector.insertSeparator(1)
            keys = set(data)
            for k in keys:
                self.stringSelector.addItem(k)
        else:
            self.stringSelectorLabel.hide()
            self.stringSelector.hide()

    def setDefaultName(self):
        self.outputName.setText("t_" + self.selectedData.currentText())
        if self.stringSelector.currentIndex() > 0:
            self.outputName.setText(
                self.outputName.text() + "_" + self.stringSelector.currentText()
            )

    def accept(self):
        key = self.selectedData.currentText()
        data = self.data[key]

        def make_key(i, ranges):
            if len(ranges) == 1:
                return self.outputName.text()
            zeros = len(str(len(ranges)))
            text = self.outputName.text() + "_{}".format(str(i).zfill(zeros))
            i += 1
            return text

        if type(data[0]) is unicode:
            strIdx = self.stringSelector.currentIndex()
            if strIdx == 0:

                def is_valid(idx):
                    return len(data[idx]) != 0

            else:
                string = self.stringSelector.currentText()

                def is_valid(idx):
                    return data[idx] == string

        else:

            def is_valid(idx):
                return not np.isnan(data[idx])

        def get_range(i0):
            assert i0 < len(data)
            if i0 == len(data) - 1:
                return None
            while i0 < len(data) and not is_valid(i0):
                i0 += 1
            if i0 == len(data):
                return None
            iN = i0
            while iN < len(data) and is_valid(iN):
                iN += 1
            if iN == len(data):
                iN = len(data) - 1
            return (i0, iN)

        ranges = [get_range(0)]
        while ranges[-1] is not None:
            ranges.append(get_range(ranges[-1][1]))
        ranges = ranges[:-1]

        def make_range(i0, iN):
            out = copy.deepcopy(self.data["t"])
            out[:i0] = np.nan
            out[iN + 1 :] = np.nan
            return out

        i = 1
        for r in ranges:
            k = make_key(i, ranges)
            self.data[k] = make_range(r[0], r[1])
            i += 1
        self.data.notify_update()
        super(GenerateRangeDialog, self).accept()


def rpyFromMat(E):
    """Same as mc_rbdyn::rpyFromMat."""
    roll = atan2(E[1][2], E[2][2])
    pitch = -asin(E[0][2])
    yaw = atan2(E[0][1], E[0][0])
    return [roll, pitch, yaw]


def rpyFromQuat(quat):
    """Same as mc_rbdyn::rpyFromQuat."""
    import eigen

    return rpyFromMat(list(eigen.Quaterniond(*quat).toRotationMatrix()))


class PlotPolygonAxis(object):
    def __init__(self, parent, axis):
        self.figure = parent
        self._axis = axis.twinx()
        self._axis.set_zorder(-1)
        self._axis.set_yticks([])
        self._axis.set_ylim(0, 1)
        self._axis.get_yaxis().set_visible(False)
        self.data = {}
        self.plots = OrderedDict()
        self.colors = OrderedDict()

    def _legend_fontsize(self):
        return self.figure._labels_legend_fontsize

    def _legend_top_offset(self):
        return self.figure._labels_legend_top_offset

    def _plot_string(self, x, y, y_label, style):
        if y_label in self.plots:
            return False
        self.plots[y_label] = []
        self.data[y_label] = [x, y]
        self.colors[y_label] = {}
        i = 0
        i0 = 0
        label = y[i0]
        while i < len(y):
            if y[i] == label and i + 1 < len(y):
                # Continue if the y value is same as that of previous timestep
                i += 1
                continue
            # Add new color if the y value is new
            new_label = False
            if label not in self.colors[y_label]:
                new_label = True
                self.colors[y_label][label] = self.figure._next_poly_color()
            color = self.colors[y_label][label]
            # Determine the last timestep of this label
            xi = x[i]
            if np.isnan(xi):
                xi = x[i - 1]
            # Make rectangle to draw
            if len(label) > 0:
                if not new_label:
                    # Do not display labels in legends from the second time
                    label = "_" + label
                self.plots[y_label].append(
                    self._axis.add_patch(
                        Rectangle(
                            (x[i0], 0), xi - x[i0], 1, label=label, facecolor=color
                        )
                    )
                )
            # Store the last value. This is the initial value for the next label.
            i0 = i
            label = y[i0]
            i += 1
        return True

    def legend(self):
        if not len(self.plots):
            self._axis.clear()
            return
        xL = 1.0
        if len(self.figure._right()):
            xL = 1.025
        self._axis.legend(
            bbox_to_anchor=(xL, self._legend_top_offset(), 0.2, 0.1),
            mode="expand",
            borderaxespad=0.5,
            fontsize=self._legend_fontsize(),
        )

    def remove_plot(self, y):
        if y not in self.plots:
            return
        for plt in self.plots[y]:
            plt.remove()
        del self.data[y]
        del self.plots[y]
        del self.colors[y]
        self.figure.draw()

    def update_x(self, x):
        keys = self.data.keys()
        for y_label in keys:
            x, y = self.data[y_label]
            self.remove_plot(y_label)
            self._plot_string(x, y, y_label, None)


class PlotYAxis(object):
    def __init__(self, parent, x_axis=None, poly=None, _3D=False):
        self.figure = parent
        self._3D = _3D
        if x_axis is None:
            self.is_left = True
            if not _3D:
                self._axis = parent.fig.add_subplot(111)
            else:
                self._axis = parent.fig.add_subplot(projection="3d")
                if self._axis is None:
                    # Fallback for older versions of matplotlib
                    self._axis = parent.fig.gca(projection="3d")
                self._axis.set_xlabel("X")
                self._axis.set_ylabel("Y")
                self._axis.set_zlabel("Z")
            self._axis.autoscale(enable=True, axis="both", tight=False)
            self._x_axis = self._axis
        else:
            self.is_left = False
            if not _3D:
                self._axis = x_axis.twinx()
            else:
                self._axis = x_axis
            self._x_axis = x_axis
        if poly is None and not _3D:
            self._polyAxis = PlotPolygonAxis(parent, self._axis)
        else:
            self._polyAxis = poly
        try:
            self._axis.set_facecolor((1, 1, 1, 0))
        except AttributeError:
            pass
        self._axis.autoscale_view(False, True, True)
        self._axis.format_coord = parent.format_coord
        self._axis.get_yaxis().set_visible(False)
        self.grid = LineStyle(linestyle="--")
        self.plots = OrderedDict()
        self._label_fontsize = 10
        self._z_label_fontsize = 10
        self._legend_ncol = 3
        self.source = {}
        self.data = {}
        self.filtered = {}

    def __len__(self):
        return len(self.plots)

    def _data(self):
        return self.figure.data

    def _legend_fontsize(self):
        return self.figure._legend_fontsize

    def _labelpad(self):
        return self.figure._labelpad

    def _x_label_fontsize(self):
        return self.figure._x_label_fontsize

    def axis(self):
        return self._axis

    def drawGrid(self):
        if len(self.plots):
            self._axis.grid(
                color=self.grid.color,
                linestyle=self.grid.linestyle,
                linewidth=self.grid.linewidth,
                visible=self.grid.visible,
                which="both",
            )
        if not self.is_left and len(self.figure.axes[PlotSide.LEFT].plots) == 0:
            self.figure.axes[PlotSide.LEFT]._axis.grid(
                color=self.grid.color,
                linestyle=self.grid.linestyle,
                linewidth=self.grid.linewidth,
                visible=self.grid.visible,
                which="both",
                axis="x",
            )

    def legend(self):
        if not len(self.plots):
            return
        if self.is_left:
            loc = 3
            top_anchor = 1.02
        else:
            loc = 2
            top_anchor = -0.14
            if len(self.figure.x_label()):
                top_anchor = -0.175
        self._axis.legend(
            bbox_to_anchor=(0.0, top_anchor, 1.0, 0.102),
            loc=loc,
            ncol=self._legend_ncol,
            mode="expand",
            borderaxespad=0.5,
            fontsize=self._legend_fontsize(),
        )
        if self._polyAxis is not None:
            self._polyAxis.legend()

    def legendNCol(self, n=None):
        if n is None:
            return self._legend_ncol
        self._legend_ncol = n
        self.legend()

    def legendRows(self):
        return math.ceil(len(self.plots) / float(self._legend_ncol))

    def legendOffset(self, offset, sign):
        if self.legendRows() > 3:
            offset = offset + sign * 0.035 * (self.legendRows() - 3)
        return offset

    def getLimits(self, frame0, frame, idx):
        if not len(self):
            return None
        if frame0 is None:
            frame0, frame = self.getFrameRange()
        values = [
            value
            for value in self.data.values()
            if not np.all(np.isnan(value[idx][frame0:frame]))
        ]
        if len(values) == 0:
            return None
        value0 = values[0]
        min_ = np.nanmin(value0[idx][frame0:frame])
        max_ = np.nanmax(value0[idx][frame0:frame])
        for value in values:
            data = value[idx][frame0:frame]
            min_ = min(np.nanmin(data), min_)
            max_ = max(np.nanmax(data), max_)
        return min_, max_

    def getFrameRange(self):
        if not len(self):
            return None, None

        def getRange(value, i0=None, iN=None):
            first, last = 0, 0
            if self._3D:
                check = np.nonzero(np.isfinite(value[0] + value[1] + value[2]))
            else:
                check = np.nonzero(np.isfinite(value[0] + value[1]))
            if len(check[0]):
                first = check[0][0]
                last = check[0][-1]
            if i0 is None:
                return first, last
            return min(i0, first), max(iN, last)

        values = iter(self.data.values())
        i0, iN = getRange(next(values))
        for value in values:
            i0, iN = getRange(value, i0, iN)
        return i0, iN

    def setLimits(self, xlim=None, ylim=None, frame0=None, frame=None, zlim=None):
        if not len(self):
            return xlim, ylim
        if frame is None:
            frame0, frame = self.getFrameRange()
        dataLim = self._axis.dataLim.get_points()

        def setMargin(min_, max_):
            return [
                min_ - 0.01 * (max_ - min_),
                max_ + 0.01 * (max_ - min_),
            ]

        def setLimit(lim, idx, set_lim):
            if lim is not None:
                min_, max_ = lim
            elif frame is not None:
                limits = self.getLimits(frame0, frame, idx)
                # Ignore limits if all data is nan
                if limits is None:
                    return None
                min_, max_ = setMargin(limits[0], limits[1])
            else:
                min_, max_ = setMargin(dataLim[0][idx], dataLim[1][idx])
            set_lim([min_, max_])
            return min_, max_

        if self._3D:
            frame = -1
            setLimit(zlim, 2, self._axis.set_zlim)
        return setLimit(xlim, 0, self._x_axis.set_xlim), setLimit(
            ylim, 1, self._axis.set_ylim
        )

    def _label(self, get_label, set_label, label, size):
        if label is None:
            label = get_label()
        set_label(label, fontsize=size, labelpad=self._labelpad())

    def _label_property(self, get_label, set_label, label=None):
        if label is None:
            return get_label()
        set_label(label)

    def _x_label(self, label=None):
        self._label(
            self.x_label, self._x_axis.set_xlabel, label, self._x_label_fontsize()
        )

    def x_label(self, label=None):
        return self._label_property(self._x_axis.get_xlabel, self._x_label, label)

    def _y_label(self, label=None):
        self._label(self.y_label, self._axis.set_ylabel, label, self._label_fontsize)

    def y_label(self, label=None):
        return self._label_property(self._axis.get_ylabel, self._y_label, label)

    def _z_label(self, label=None):
        self._label(self.z_label, self._axis.set_zlabel, label, self._z_label_fontsize)

    def z_label(self, label=None):
        return self._label_property(self._axis.get_zlabel, self._z_label, label)

    def y_label_fontsize(self, fontsize=None):
        if fontsize is None:
            return self._label_fontsize
        self._label_fontsize = fontsize
        self._y_label()

    def z_label_fontsize(self, fontsize=None):
        if fontsize is None:
            return self._z_label_fontsize
        self._z_label_fontsize = fontsize
        self._z_label()

    def animate(self, frame0, frame):
        for y_label in self.plots.keys():
            self.plots[y_label].set_data(
                self.data[y_label][0][frame0:frame], self.data[y_label][1][frame0:frame]
            )
            if self._3D:
                self.plots[y_label].set_3d_properties(
                    self.data[y_label][2][frame0:frame]
                )
        return list(self.plots.values())

    def update_x(self, x):
        styles = {}
        for y_label in self.data.keys():
            filter_ = self.data[y_label][-1]
            if filter_ is None:
                self.data[y_label][0] = x
            else:
                filter_ = x
            styles[y_label] = self.style(y_label)
        self.clear()
        for y_label in self.data.keys():
            x = self.data[y_label][0]
            y = self.data[y_label][1]
            z = self.data[y_label][2]
            self._plot(x, y, y_label, styles[y_label], filter_=filter_, z=z)

    def _filter(self, data, filter_):
        if data is None:
            return None
        out = copy.deepcopy(data)
        for i, d in enumerate(filter_):
            if np.isnan(d):
                out[i] = np.nan
        return out

    def _plot(self, x, y, y_label, style=None, filter_=None, z=None, source=None):
        if type(y[0]) is unicode:
            if filter_ is not None:
                return False
            return self._polyAxis._plot_string(x, y, y_label, style)
        if style is None:
            return self._plot(
                x,
                y,
                y_label,
                LineStyle(color=self.figure._next_color()),
                filter_=filter_,
                z=z,
                source=source,
            )
        if y_label in self.plots:
            return False
        self._axis.get_yaxis().set_visible(True)
        self.data[y_label] = [x, y, z, filter_]
        if filter_ is not None:
            x = self._filter(x, filter_)
            y = self._filter(y, filter_)
            z = self._filter(z, filter_)
            self.filtered[y_label] = [x, y, z]
        else:
            self.filtered[y_label] = None
        self.source[y_label] = source
        if z is None:
            self.plots[y_label] = self._axis.plot(
                x,
                y,
                label=y_label,
                color=style.color,
                linestyle=style.linestyle,
                linewidth=style.linewidth,
            )[0]
        else:
            self.plots[y_label] = self._axis.plot(
                x,
                y,
                z,
                label=y_label,
                color=style.color,
                linestyle=style.linestyle,
                linewidth=style.linewidth,
            )[0]
        self.legend()
        return True

    def startAnimation(self, i0):
        for y_label in self.plots.keys():
            plotStyle = self.style(y_label)
            self.plots[y_label].remove()
            if not self._3D:
                self.plots[y_label] = self._axis.plot(
                    self.data[y_label][0][i0], self.data[y_label][1][i0], label=y_label
                )[0]
            else:
                self.plots[y_label] = self._axis.plot(
                    [self.data[y_label][0][i0]],
                    [self.data[y_label][1][i0]],
                    [self.data[y_label][2][0]],
                    label=y_label,
                )[0]
            self.style(y_label, plotStyle)

    def stopAnimation(self):
        for y_label in self.plots.keys():
            plotStyle = self.style(y_label)
            self.plots[y_label].remove()
            if self.filtered[y_label] is not None:
                if not self._3D:
                    self.plots[y_label] = self._axis.plot(
                        self.filtered[y_label][0],
                        self.filtered[y_label][1],
                        label=y_label,
                    )[0]
                else:
                    self.plots[y_label] = self._axis.plot(
                        self.filtered[y_label][0],
                        self.filtered[y_label][1],
                        self.filtered[y_label][2],
                        label=y_label,
                    )[0]
            else:
                self.plots[y_label] = self._axis.plot(
                    self.data[y_label][0], self.data[y_label][1], label=y_label
                )[0]
            self.style(y_label, plotStyle)

    def add_plot(self, x, y, y_label, style=None):
        return self._plot(self._data()[x], self._data()[y], y_label, style, source=y)

    def add_plot_xy(self, x, y, y_label, t, style=None):
        return self._plot(
            self._data()[x],
            self._data()[y],
            y_label,
            style,
            filter_=self._data()[t],
            source=[x, y, y_label],
        )

    def add_plot_xyz(self, x, y, z, y_label, t, style=None):
        return self._plot(
            self._data()[x],
            self._data()[y],
            y_label,
            style,
            filter_=self._data()[t],
            z=self._data()[z],
            source=[x, y, z, y_label],
        )

    def _add_diff_plot(self, x, data, y_label):
        dt = self._data()[x][1] - self._data()[x][0]
        return self._plot(self._data()[x][1:], np.diff(data) / dt, y_label)

    def add_diff_plot(self, x, y, y_label):
        return self._add_diff_plot(x, self._data()[y], y_label)

    def _get_rpy_data(self, y, idx):
        assert idx >= 0 and idx <= 2, "index must be 0, 1 or 2"
        rpy_label = ["roll", "pitch", "yaw"]
        y_label = "{}_{}".format(y, rpy_label[idx])
        fmt = ""
        if "{}_qw".format(y) in self._data().keys():
            fmt = "q"
        qw, qx, qy, qz = [
            self._data()[k]
            for k in ["{}_{}{}".format(y, fmt, ax) for ax in ["w", "x", "y", "z"]]
        ]
        data = [rpyFromQuat([w, x, y, z])[idx] for w, x, y, z in zip(qw, qx, qy, qz)]
        return data, y_label

    def _add_rpy_plot(self, x_label, y, idx):
        data, y_label = self._get_rpy_data(y, idx)
        return self._plot(self._data()[x_label], data, y_label)

    def _add_rpy_diff_plot(self, x_label, y, idx):
        data, y_label = self._get_rpy_data(y, idx)
        y_label += " diff"
        return self._add_diff_plot(x_label, data, y_label)

    def add_roll_plot(self, x, y):
        return self._add_rpy_plot(x, y, 0)

    def add_roll_diff_plot(self, x, y):
        return self._add_rpy_diff_plot(x, y, 0)

    def add_pitch_plot(self, x, y):
        return self._add_rpy_plot(x, y, 1)

    def add_pitch_diff_plot(self, x, y):
        return self._add_rpy_diff_plot(x, y, 1)

    def add_yaw_plot(self, x, y):
        return self._add_rpy_plot(x, y, 2)

    def add_yaw_diff_plot(self, x, y):
        return self._add_rpy_diff_plot(x, y, 2)

    def add_rpy_plot(self, x, y):
        r = self.add_roll_plot(x, y)
        p = self.add_pitch_plot(x, y)
        y = self.add_yaw_plot(x, y)
        return r or p or y

    def add_rpy_diff_plot(self, x, y):
        r = self.add_roll_diff_plot(x, y)
        p = self.add_pitch_diff_plot(x, y)
        y = self.add_yaw_diff_plot(x, y)
        return r or p or y

    def remove_plot(self, y):
        if y not in self.plots:
            self._polyAxis.remove_plot(y)
            for y_label, source in self.source.items():
                if source == y:
                    return self.remove_plot(y_label)
            return
        self.plots[y].remove()
        del self.plots[y]
        del self.data[y]
        del self.filtered[y]
        del self.source[y]
        if len(self.plots):
            self._axis.relim()
            self.legend()
        else:
            self._axis.get_yaxis().set_visible(False)
            self._axis.clear()

    def clear(self):
        self.plots = {}
        self._axis.clear()

    # Get or set the style of a given plot
    def style(self, y, style=None):
        if y not in self.plots:
            raise KeyError("No plot named {}".format(y))
        plt = self.plots[y]
        if style is None:
            return LineStyle(
                plt.get_color(),
                plt.get_linestyle(),
                plt.get_linewidth(),
                label=plt.get_label(),
            )
        plt.set_color(style.color)
        plt.set_linestyle(style.linestyle)
        plt.set_linewidth(style.linewidth)
        if len(style.label):
            plt.set_label(style.label)


class PlotFigure(object):
    def __init__(self, type_=PlotType.TIME, init_canvas=False):
        self.fig = matplotlib.pyplot.figure(figsize=(5, 4), dpi=100)
        if init_canvas:
            self.canvas = FigureCanvas(self.fig)
        self.axes = {}
        self._3D = type_ is PlotType._3D
        self._XY = type_ is PlotType.XY
        self.axes[PlotSide.LEFT] = PlotYAxis(self, _3D=self._3D)
        if not self._3D:
            self.axes[PlotSide.RIGHT] = PlotYAxis(
                self, self._left().axis(), self._left()._polyAxis, _3D=self._3D
            )
        else:
            self.axes[PlotSide.RIGHT] = None
        self.animation = None

        self._title_fontsize = 12
        self._x_label_fontsize = 10
        self._labelpad = 10
        self._tick_labelsize = 10
        self._legend_fontsize = 10
        self._labels_legend_fontsize = 10
        self._labels_legend_top_offset = 0.9
        self._top_offset = 0.9
        self._bottom_offset = 0.1
        if self._3D:
            self._bottom_offset = 0

        self.data = None
        self.computed_data = {}

        self.color = 0
        cm = matplotlib.cm.Set1
        self.Ncolor = min(cm.N, 12)
        cm2rgb = (np.array(cm(x)[0:3]) for x in np.linspace(0, 1, self.Ncolor))
        self.colors = [
            "#%02x%02x%02x" % tuple((255 * rgb).astype(int)) for rgb in cm2rgb
        ]

        self.polyColor = 0
        cm = matplotlib.cm.Pastel1
        cm2rgb = (
            np.array(cm(x)[0:3] + (0.5,)) for x in np.linspace(0, 1, min(cm.N, 32))
        )
        self.polyColors = [rgb for rgb in cm2rgb]

        self.x_data = "t"

    # Helper function to call something on all axes, call expects an axis argument
    def _axes(self, call):
        call(self._left())
        if self._right() is not None:
            call(self._right())

    # Shortcut to get the polygon axis
    def _polygons(self):
        return self._left()._polyAxis

    # Shortcut to the left axis
    def _left(self):
        return self.axes[PlotSide.LEFT]

    # Shortcut to the right axis
    def _right(self):
        return self.axes[PlotSide.RIGHT]

    def _drawGrid(self):
        self._axes(lambda axis: axis.drawGrid())

    def _legend(self):
        self._axes(lambda axis: axis.legend())

    def draw(
        self, x_limits=None, y1_limits=None, y2_limits=None, frame0=None, frame=None
    ):
        if x_limits is None:
            # No limit specified for XY/time plots, compute the best min-max x limits fitting both the left and right curves
            x1_limits = None
            x2_limits = None
            if self._left() is not None:
                x1_limits = self._left().getLimits(frame0, frame, 0)
            if self._right() is not None:
                x2_limits = self._right().getLimits(frame0, frame, 0)
            if x1_limits is None and x2_limits is not None:
                x_limits = x2_limits
            elif x2_limits is None and x1_limits is not None:
                x_limits = x1_limits
            elif x1_limits is None and x2_limits is None:
                return
            else:
                x_limits = [
                    min(x1_limits[0], x2_limits[0]),
                    max(x1_limits[1], x2_limits[1]),
                ]
            if x_limits is None:
                return
            range_ = x_limits[1] - x_limits[0]
            x_limits = [x_limits[0] - range_ * 0.01, x_limits[1] + range_ * 0.01]

        if self._3D:
            x_limits, _ = self._left().setLimits(
                x_limits, y1_limits, frame0=frame0, frame=frame, zlim=y2_limits
            )
        else:
            x_limits, y1_new_limits = self._left().setLimits(
                x_limits, y1_limits, frame0=frame0, frame=frame
            )
            x_limits, y2_new_limits = self._right().setLimits(
                x_limits, y2_limits, frame0=frame0, frame=frame
            )
            # Limits are not pre-set and there is actual data on both axis
            if (
                y1_limits is None
                and y2_limits is None
                and y1_new_limits is not None
                and y2_new_limits is not None
            ):
                y1_mid = (y1_new_limits[0] + y1_new_limits[1]) / 2
                y2_mid = (y2_new_limits[0] + y2_new_limits[1]) / 2
                y1_y2_ratio = y1_mid / y2_mid
                if y1_y2_ratio > 0.2 and y1_y2_ratio < 5:
                    y_min = min(y1_new_limits[0], y2_new_limits[0])
                    y_max = max(y1_new_limits[1], y2_new_limits[1])
                    self._left()._axis.set_ylim(y_min, y_max)
                    self._right()._axis.set_ylim(y_min, y_max)

        self._legend()
        self._drawGrid()
        top_offset = self._left().legendOffset(self._top_offset, -1)
        if self._right() is not None:
            bottom_offset = self._right().legendOffset(self._bottom_offset, 1)
        else:
            bottom_offset = self._bottom_offset
        left_offset = 0.125
        right_offset = 0.9
        if self._left()._polyAxis is not None and len(self._left()._polyAxis.plots):
            left_offset, right_offset = 0.05, right_offset - (left_offset - 0.05)
        self.fig.subplots_adjust(
            left=left_offset, right=right_offset, top=top_offset, bottom=bottom_offset
        )

    def animate(self, frame0, frame, x_limits=None, y1_limits=None, y2_limits=None):
        if self.onlyShowCheckbox.isChecked():
            # Handle only show the last N seconds
            frame0 = max(
                frame0, frame - int(self.onlyShowSpinBox.value() / (self.dt / 1000))
            )
        ret = self._left().animate(frame0, frame)
        if self._right():
            ret.extend(self._right().animate(frame0, frame))
        PlotFigure.draw(
            self,
            x_limits=x_limits,
            y1_limits=y1_limits,
            y2_limits=y2_limits,
            frame0=frame0,
            frame=frame,
        )
        return ret

    def stopAnimation(self):
        self._axes(lambda a: a.stopAnimation())

    def setData(self, data):
        self.data = data

    def setColors(self, colors):
        self.colors = colors
        self.Ncolor = len(self.colors)

    def setPolyColors(self, colors):
        self.polyColors = colors

    def show(self):
        self.fig.show()

    def top_offset(self, off=None):
        if off is None:
            return self._top_offset
        else:
            self._top_offset = off

    def bottom_offset(self, off=None):
        if off is None:
            return self._bottom_offset
        else:
            self._bottom_offset = off

    def title(self, title=None):
        if title is None:
            if self.fig._suptitle is None:
                return ""
            return self.fig._suptitle.get_text()
        self.fig.suptitle(title)

    def title_fontsize(self, fontsize=None):
        if fontsize is None:
            return self._title_fontsize
        self._title_fontsize = fontsize
        self.fig.suptitle(self.title(), fontsize=self._title_fontsize)

    def tick_fontsize(self, size=None):
        if size is None:
            return self._tick_labelsize
        self._tick_labelsize = size
        self._axes(lambda a: a.axis().tick_params(labelsize=self._tick_labelsize))

    def labelpad(self, pad=None):
        if pad is None:
            return self._labelpad
        self._labelpad = pad
        self._left()._x_label()
        self._axes(lambda axis: axis._y_label())

    def _x_label(self, label=None):
        self._left()._x_label(label)

    def _y1_label(self, label=None):
        self._left()._y_label(label)

    def _y2_label(self, label=None):
        self._right()._y_label(label)

    def x_label(self, label=None):
        return self._left().x_label(label)

    def x_label_fontsize(self, fontsize=None):
        if fontsize is None:
            return self._x_label_fontsize
        self._x_label_fontsize = fontsize
        self._x_label()

    def y1_label(self, label=None):
        return self._left().y_label(label)

    def y1_label_fontsize(self, fontsize=None):
        return self._left().y_label_fontsize(fontsize)

    def y2_label(self, label=None):
        if self._3D:
            return self._left().z_label(label)
        return self._right().y_label(label)

    def y2_label_fontsize(self, fontsize=None):
        if self._3D:
            return self._left().z_label_fontsize(fontsize)
        return self._right().y_label_fontsize(fontsize)

    def _next_poly_color(self):
        self.polyColor += 1
        return self.polyColors[(self.polyColor - 1) % len(self.polyColors)]

    def _next_color(self):
        self.color += 1
        return self.colors[(self.color - 1) % self.Ncolor]

    def legend_fontsize(self, size=None):
        if size is None:
            return self._legend_fontsize
        self._legend_fontsize = size
        self._legend()

    def labels_legend_fontsize(self, size=None):
        if size is None:
            return self._labels_legend_fontsize
        self._labels_legend_fontsize = size
        self._legend()

    def labels_legend_top_offset(self, offset=None):
        if offset is None:
            return self._labels_legend_top_offset
        self._labels_legend_top_offset = offset
        self._legend()

    def y1_legend_ncol(self, n=None):
        return self._left().legendNCol(n)

    def y2_legend_ncol(self, n=None):
        if self._right() is not None:
            return self._right().legendNCol(n)
        else:
            return 0

    def add_plot_left(self, x, y, y_label, style=None):
        return self._left().add_plot(x, y, y_label, style)

    def add_plot_left_xy(self, x, y, y_label, style=None):
        return self._left().add_plot_xy(x, y, y_label, self.x_data, style)

    def add_plot_left_xyz(self, x, y, z, y_label, style=None):
        return self._left().add_plot_xyz(x, y, z, y_label, self.x_data, style)

    def add_plot_right(self, x, y, y_label, style=None):
        return self._right().add_plot(x, y, y_label, style)

    def add_plot_right_xy(self, x, y, y_label, style=None):
        return self._right().add_plot_xy(x, y, y_label, self.x_data, style)

    def add_plot_right_xyz(self, x, y, z, y_label, style=None):
        return self._right().add_plot_xyz(x, y, z, y_label, self.x_data, style)

    def add_diff_plot_left(self, x, y, y_label):
        return self._left().add_diff_plot(x, y, y_label)

    def add_diff_plot_right(self, x, y, y_label):
        return self._right().add_diff_plot(x, y, y_label)

    def add_roll_plot_left(self, x, y):
        return self._left().add_roll_plot(x, y)

    def add_pitch_plot_left(self, x, y):
        return self._left().add_pitch_plot(x, y)

    def add_yaw_plot_left(self, x, y):
        return self._left().add_yaw_plot(x, y)

    def add_roll_diff_plot_left(self, x, y):
        return self._left().add_roll_diff_plot(x, y)

    def add_pitch_diff_plot_left(self, x, y):
        return self._left().add_pitch_diff_plot(x, y)

    def add_yaw_diff_plot_left(self, x, y):
        return self._left().add_yaw_diff_plot(x, y)

    def add_roll_plot_right(self, x, y):
        return self._right().add_roll_plot(x, y)

    def add_pitch_plot_right(self, x, y):
        return self._right().add_pitch_plot(x, y)

    def add_yaw_plot_right(self, x, y):
        return self._right().add_yaw_plot(x, y)

    def add_roll_diff_plot_right(self, x, y):
        return self._right().add_roll_diff_plot(x, y)

    def add_pitch_diff_plot_right(self, x, y):
        return self._right().add_pitch_diff_plot(x, y)

    def add_yaw_diff_plot_right(self, x, y):
        return self._right().add_yaw_diff_plot(x, y)

    def add_rpy_plot_left(self, x, y):
        return self._left().add_rpy_plot(x, y)

    def add_rpy_diff_plot_left(self, x, y):
        return self._left().add_rpy_diff_plot(x, y)

    def add_rpy_plot_right(self, x, y):
        return self._right().add_rpy_plot(x, y)

    def add_rpy_diff_plot_right(self, x, y):
        return self._right().add_rpy_diff_plot(x, y)

    def _remove_plot(self, SIDE, y_label):
        self.axes[SIDE].remove_plot(y_label)
        if len(self._left()) == 0 and self._right() and len(self._right()) == 0:
            self.color = 0

    def remove_plot_left(self, y_label):
        self._remove_plot(PlotSide.LEFT, y_label)

    def remove_plot_right(self, y_label):
        self._remove_plot(PlotSide.RIGHT, y_label)

    def format_coord(self, x, y):
        if self._right() is not None:
            display_coord = self.axes[PlotSide.RIGHT].axis().transData.transform((x, y))
            inv = self.axes[PlotSide.LEFT].axis().transData.inverted()
            ax_coord = inv.transform(display_coord)
            if len(self._left()) and len(self._right()):
                return "x: {:.3f}    y1: {:.3f}    y2: {:.3f}".format(x, ax_coord[1], y)
            elif len(self._left()):
                return "x: {:.3f}    y1: {:.3f}".format(x, ax_coord[1])
            elif len(self._right()):
                return "x: {:.3f}    y2: {:.3f}".format(x, y)
            else:
                return "x: {:.3f}".format(x)
        else:
            axis = self.axes[PlotSide.LEFT].axis()
            return type(axis).format_coord(axis, x, y)

    def clear_all(self):
        self.color = 0
        self._axes(lambda a: a.clear())

    def style_left(self, y, styleIn=None):
        return self._left().style(y, styleIn)

    def style_right(self, y, styleIn=None):
        return self._right().style(y, styleIn)


class SimpleAxesDialog(QtWidgets.QDialog):
    def __init__(self, parent):
        QtWidgets.QDialog.__init__(self, parent)
        self.setWindowTitle("Edit axes limits")
        self.setModal(True)
        self.layout = QtWidgets.QGridLayout(self)
        self.layout.addWidget(QtWidgets.QLabel("Min"), 0, 1)
        self.layout.addWidget(QtWidgets.QLabel("Max"), 0, 2)

        self.layout.addWidget(QtWidgets.QLabel("X"), 1, 0)
        x_limits = parent.x_limits
        if x_limits is None:
            x_limits = parent._left().axis().get_xlim()
        self.x_min = QtWidgets.QLineEdit(str(x_limits[0]))
        self.x_min.setValidator(QtGui.QDoubleValidator())
        self.layout.addWidget(self.x_min, 1, 1)
        self.x_max = QtWidgets.QLineEdit(str(x_limits[1]))
        self.x_max.setValidator(QtGui.QDoubleValidator())
        self.x_init = [float(self.x_min.text()), float(self.x_max.text())]
        self.layout.addWidget(self.x_max, 1, 2)

        self.layout.addWidget(QtWidgets.QLabel("Y1"), 2, 0)
        y1_limits = parent.y1_limits
        if y1_limits is None:
            y1_limits = parent._left().axis().get_ylim()
        self.y1_min = QtWidgets.QLineEdit(str(y1_limits[0]))
        self.y1_min.setValidator(QtGui.QDoubleValidator())
        self.layout.addWidget(self.y1_min, 2, 1)
        self.y1_max = QtWidgets.QLineEdit(str(y1_limits[1]))
        self.y1_max.setValidator(QtGui.QDoubleValidator())
        self.y1_init = [float(self.y1_min.text()), float(self.y1_max.text())]
        self.layout.addWidget(self.y1_max, 2, 2)

        y2_label = "Y2"
        if parent._3D:
            y2_label = "Z"
        self.layout.addWidget(QtWidgets.QLabel(y2_label), 3, 0)
        y2_limits = parent.y2_limits
        if y2_limits is None:
            if parent._3D:
                y2_limits = parent._left().axis().get_zlim()
            else:
                y2_limits = parent._right().axis().get_ylim()
        self.y2_min = QtWidgets.QLineEdit(str(y2_limits[0]))
        self.y2_min.setValidator(QtGui.QDoubleValidator())
        self.layout.addWidget(self.y2_min, 3, 1)
        self.y2_max = QtWidgets.QLineEdit(str(y2_limits[1]))
        self.y2_max.setValidator(QtGui.QDoubleValidator())
        self.y2_init = [float(self.y2_min.text()), float(self.y2_max.text())]
        self.layout.addWidget(self.y2_max, 3, 2)

        confirmLayout = QtWidgets.QHBoxLayout()
        okButton = QtWidgets.QPushButton("Ok", self)
        confirmLayout.addWidget(okButton)
        okButton.clicked.connect(self.accept)
        applyButton = QtWidgets.QPushButton("Apply", self)
        confirmLayout.addWidget(applyButton)
        applyButton.clicked.connect(self.apply)
        if not parent._3D:
            alignButton = QtWidgets.QPushButton("Align left/right", self)
            confirmLayout.addWidget(alignButton)
            alignButton.clicked.connect(self.align)
        cancelButton = QtWidgets.QPushButton("Cancel", self)
        confirmLayout.addWidget(cancelButton)
        cancelButton.clicked.connect(self.reject)
        self.layout.addLayout(confirmLayout, 4, 0, 1, 3)

    def get_limits(self):
        x_limits = [float(self.x_min.text()), float(self.x_max.text())]
        y1_limits = [float(self.y1_min.text()), float(self.y1_max.text())]
        y2_limits = [float(self.y2_min.text()), float(self.y2_max.text())]
        return x_limits, y1_limits, y2_limits

    def apply(self):
        changed = False
        x_limits, y1_limits, y2_limits = self.get_limits()
        if x_limits != self.x_init:
            changed = True
            self.parent().x_locked.setChecked(True)
            self.parent().x_limits = x_limits
        if y1_limits != self.y1_init:
            changed = True
            self.parent().y1_locked.setChecked(True)
            self.parent().y1_limits = y1_limits
        if y2_limits != self.y2_init:
            changed = True
            self.parent().y2_locked.setChecked(True)
            self.parent().y2_limits = y2_limits
        if changed:
            self.parent().draw()

    def align(self):
        _, y1_limits, y2_limits = self.get_limits()
        y1_limits[0] = min(y1_limits[0], y2_limits[0])
        y1_limits[1] = max(y1_limits[1], y2_limits[1])
        self.parent().y1_locked.setChecked(True)
        self.parent().y1_limits = y1_limits
        self.parent().y2_locked.setChecked(True)
        self.parent().y2_limits = y1_limits
        self.parent().draw()

    def accept(self):
        QtWidgets.QDialog.accept(self)
        self.apply()


class SaveAnimationDialog(QtWidgets.QDialog):
    def __init__(self, parent):
        super(QtWidgets.QDialog, self).__init__(parent)
        self.setModal(True)
        self.parent = parent
        self.layout = QtWidgets.QVBoxLayout(self)
        label = QtWidgets.QLabel("Saving animation, please wait...", self)
        self.layout.addWidget(label)
        self.progress = QtWidgets.QProgressBar(self)
        self.layout.addWidget(self.progress)

    def closeEvent(self, event):
        event.ignore()


class PlotCanvasWithToolbar(PlotFigure, QWidget):
    def __init__(self, parent=None, mode=PlotType.TIME):
        PlotFigure.__init__(self, mode, True)
        QWidget.__init__(self, parent)

        self.canvas.mpl_connect("draw_event", self.on_draw)
        self.canvas.mpl_connect("motion_notify_event", self.on_hover)
        self.canvas.mpl_connect("axes_leave_event", self.on_axes_leave)
        self.canvas.mpl_connect("scroll_event", self.on_scroll)
        self.toolbar = NavigationToolbar(self.canvas, self)

        self._hover_tooltip = QtWidgets.QLabel(self.canvas)
        self._hover_tooltip.setStyleSheet(
            "QLabel { background-color: rgba(40, 40, 40, 220); color: white;"
            " border: 1px solid #888; border-radius: 4px; padding: 4px 6px;"
            " font-size: 11px; }"
        )
        self._hover_tooltip.setAlignment(QtCore.Qt.AlignLeft)
        self._hover_tooltip.setAttribute(QtCore.Qt.WA_TransparentForMouseEvents)
        self._hover_tooltip.hide()
        self._hover_pick_radius = 8  # pixels

        # Zoom state: independent of the user-facing lock buttons.
        # None means no zoom override is active for that axis.
        self._zoom_x_limits = None
        self._zoom_y1_limits = None
        self._zoom_y2_limits = None

        self.layout = QVBoxLayout(self)
        self.layout.addWidget(self.canvas)
        self.layout.addWidget(self.toolbar)

        self.setupLockButtons()
        self.setupAnimationButtons()

        self.animation_path = QtCore.QStandardPaths.writableLocation(
            QtCore.QStandardPaths.DocumentsLocation
        )
        self.savingAnimation = False
        self.saveAnimationDialog = SaveAnimationDialog(self)

    def setupLockButtons(self):
        layout = QHBoxLayout()

        self.xSelector = QComboBox(self)
        self.xSelector.activated.connect(
            lambda: self.parent().on_xSelector_activated(
                self, self.xSelector.currentText()
            )
        )
        sizePolicy = QtWidgets.QSizePolicy(
            QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.xSelector.sizePolicy().hasHeightForWidth())
        self.xSelector.setSizePolicy(sizePolicy)
        layout.addWidget(self.xSelector)

        self.genRangeButton = QtWidgets.QPushButton("Make range")
        self.genRangeButton.released.connect(lambda: GenerateRangeDialog(self).exec_())
        layout.addWidget(self.genRangeButton)

        self.x_locked = QtWidgets.QPushButton("Lock X", self)
        self.x_locked.setCheckable(True)
        layout.addWidget(self.x_locked)
        self.x_locked.toggled.connect(self.x_locked_changed)
        self.x_limits = None

        self.y1_locked = QtWidgets.QPushButton("Lock Y1", self)
        self.y1_locked.setCheckable(True)
        layout.addWidget(self.y1_locked)
        self.y1_locked.toggled.connect(self.y1_locked_changed)
        self.y1_limits = None

        self.y2_locked = QtWidgets.QPushButton("Lock Y2", self)
        if self._3D:
            self.y2_locked.setText("Lock Z")
        self.y2_locked.setCheckable(True)
        layout.addWidget(self.y2_locked)
        self.y2_locked.toggled.connect(self.y2_locked_changed)
        self.y2_limits = None

        if not self._3D:
            self.align_y1_y2 = QtWidgets.QPushButton("Align Y1/Y2", self)
            layout.addWidget(self.align_y1_y2)
            self.align_y1_y2.released.connect(self.align_y1_y2_released)

        self.layout.addLayout(layout)

    def setupAnimationButtons(self):
        animationLayout = QtWidgets.QHBoxLayout()
        self.animation = None
        self.animationButton = QtWidgets.QPushButton("Start animation")
        self.animationButton.setCheckable(True)
        self.animationButton.toggled.connect(self.startStopAnimation)
        animationLayout.addWidget(self.animationButton)
        self.lockAxesButton = QtWidgets.QPushButton("Lock axes")
        self.lockAxesButton.released.connect(self.lockAxes)
        animationLayout.addWidget(self.lockAxesButton)
        self.saveAnimationButton = QtWidgets.QPushButton("Save animation")
        self.saveAnimationButton.released.connect(self.saveAnimation)
        animationLayout.addWidget(self.saveAnimationButton)
        self.onlyShowCheckbox = QtWidgets.QCheckBox("Only show: ")
        animationLayout.addWidget(self.onlyShowCheckbox)
        self.onlyShowSpinBox = QtWidgets.QDoubleSpinBox()
        self.onlyShowSpinBox.setMinimum(0)
        self.onlyShowSpinBox.setSuffix("s")
        self.onlyShowSpinBox.setValue(10)
        animationLayout.addWidget(self.onlyShowSpinBox)
        self.layout.addLayout(animationLayout)

    def setData(self, data):
        super(PlotCanvasWithToolbar, self).setData(data)
        self.update_x_selector()

    def update_x_selector(self):
        self.xSelector.clear()
        self.xSelector.addItems(sorted(self.data.keys()))
        idx = self.xSelector.findText(self.x_data)
        if idx != -1:
            self.xSelector.setCurrentIndex(idx)

    def startStopAnimation(self):
        if self.animationButton.isChecked():
            if self.startAnimation():
                self.animationButton.setText("Stop animation")
            else:
                self.animationButton.setChecked(False)
        else:
            if self.animation is not None:
                self.stopAnimation()
            self.animationButton.setText("Start animation")

    def restartAnimation(self):
        if self.animationButton.isChecked():
            self.stopAnimation()
            self.startAnimation()

    def update_x(self):
        self._zoom_x_limits = None
        self._zoom_y1_limits = None
        self._zoom_y2_limits = None
        self._axes(lambda a: a.update_x(self.data[self.x_data]))
        if self._polygons():
            self._polygons().update_x(self.data[self.x_data])
        self.restartAnimation()
        self.draw()

    def getAxisFrameRange(self, x_data, yValues):
        """Get range from first non-nan data to last non-nan data of a single y axis
        For XY plots and 3D plots data is considered valid if all components are non-nan
        """
        min_i0 = len(x_data) - 1
        max_iN = 0
        if len(yValues) == 0:
            return min_i0, max_iN

        def isNull(ydata, idx):
            if self._3D:
                # data is [x, y, z, t]
                return (
                    lambda ydata, idx: (ydata[0] is None or np.isnan(ydata[0][idx]))
                    or (ydata[1] is None or np.isnan(ydata[1][idx]))
                    or (ydata[2] is None or np.isnan(ydata[2][idx]))
                )
            elif self._XY:
                # data is [x, y, None, t]
                return (ydata[0] is None or np.isnan(ydata[0][idx])) or (
                    ydata[1] is None or np.isnan(ydata[1][idx])
                )
            else:
                # data is [t, y]
                return ydata[1] is None or np.isnan(ydata[1][idx])

        for y_data in yValues:
            # look for first non-nan value
            i0 = 0
            while i0 < len(x_data) and (np.isnan(x_data[i0]) or isNull(y_data, i0)):
                i0 += 1
            iN = len(x_data) - 1
            while iN > i0 and (np.isnan(x_data[iN]) or isNull(y_data, iN)):
                iN -= 1
            min_i0 = min(min_i0, i0)
            max_iN = max(max_iN, iN)
        return min_i0, max_iN

    def getFrameRange(self):
        """Returns the index corresponding to the time from the first non-nan element to the last non-nan element in the non-nan region of the selected time range"""
        if self.data is None or len(self.data) == 0:
            return 0, 0
        if self._left() is None and self.right() is None:
            return 0, 0

        x_data = self.data[self.x_data]  # time
        if self._left() is not None and self._right() is None:
            return self.getAxisFrameRange(x_data, self._left().data.values())
        if self._right() is not None and self._left() is None:
            return self.getAxisFrameRange(x_data, self._right().data.values())

        y1_data = self._left().data.values()
        y2_data = self._right().data.values()
        if len(y1_data) == 0 and len(y2_data) == 0:
            return 0, 0

        y1_range = self.getAxisFrameRange(x_data, y1_data)
        y2_range = self.getAxisFrameRange(x_data, y2_data)
        range = [min(y1_range[0], y2_range[0]), max(y1_range[1], y2_range[1])]
        assert range[1] > range[0] and range[0] < len(x_data), "Strange frame range"
        return range

    def lockAxes(self):
        i0, iN = self.getFrameRange()
        if i0 == iN:
            return
        if self.x_limits is None:
            self.x_limits = self._left().getLimits(
                i0, iN, 0
            ) or self._right().getLimits(i0, iN, 0)
        if self.x_limits is None:
            return
        self.x_locked.setChecked(True)
        if self.y1_limits is None:
            self.y1_limits = self._left().getLimits(i0, iN, 1)
        if self.y1_limits is not None:
            self.y1_locked.setChecked(True)
        if self.y2_limits is None:
            if self._3D:
                self.y2_limits = self._left().getLimits(i0, iN, 2)
            else:
                self.y2_limits = self._right().getLimits(i0, iN, 1)
        if self.y2_limits is not None:
            self.y2_locked.setChecked(True)
        return i0, iN

    def startAnimation(self):
        interval = 50  # ms
        i0, iN = self.getFrameRange()
        if i0 == iN:
            return False
        x_data = self.data[self.x_data]
        self.dt = (x_data[i0 + 1] - x_data[i0]) * 1000  # dt in ms
        step = int(math.ceil(interval / self.dt))
        self.frame0 = i0
        self.saveAnimationDialog.progress.setRange(i0, iN)
        self.animation = FuncAnimation(
            self.fig, self.animate, frames=range(i0 + 1, iN, step), interval=interval
        )
        self._axes(lambda a: a.startAnimation(i0))
        self.draw()
        return True

    def animate(self, frame):
        # Matplotlib's animation.save() call is blocking, and cannot be called from a thread.
        # However, it calls this animate function for each frame. Thus we manually call Qt's event loop to keep rendering the GUI
        # and allow the saving dialog to display, and the plot to be displayed on screen
        if self.savingAnimation:
            self.saveAnimationDialog.progress.setValue(frame)
            QtWidgets.QApplication.processEvents()
        return PlotFigure.animate(
            self, self.frame0, frame, self.x_limits, self.y1_limits, self.y2_limits
        )

    def stopAnimation(self):
        self.animation.event_source.stop()
        PlotFigure.stopAnimation(self)
        self.draw()

    def saveAnimation(self):
        if self.animationButton.isChecked():
            self.animationButton.toggle()
        if not self.animationButton.isChecked():
            self.startAnimation()
            if not self.animation:
                ShowErrorDialog("Can't save an empty animation")
                return
            self.stopAnimation()

        fpath = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Output animation",
            self.animation_path + "/output-animation.mp4",
            filter="Video (*.mp4)",
        )[0]
        if not len(fpath):
            return
        self.animation_path = os.path.split(fpath)[0]
        filename, extension = os.path.splitext(fpath)
        if not extension:
            fpath = fpath + ".mp4"

        self.savingAnimation = True
        self.saveAnimationDialog.show()
        self.animation.save(fpath)
        self.savingAnimation = False
        self.saveAnimationDialog.hide()

    def axesDialog(self):
        SimpleAxesDialog(self).exec_()

    def on_axes_leave(self, event):
        """Hide tooltip and marker when the cursor leaves the plot area."""
        self._hover_tooltip.hide()
        self._clear_hover_marker()

    def _get_data_bounds(self):
        """Return (x_bounds, y1_bounds, y2_bounds) as (min, max) tuples or None.
        These are the full extents of the actual plotted data, used to clamp zoom-out.
        A 1 % margin is added on each side to match the padding PlotFigure.draw() and
        setMargin() apply to the initial auto-fit view."""

        _MARGIN = 0.01  # must match the 0.01 used in PlotFigure.draw() / setMargin()

        def axis_bounds(plot_axis, data_idx):
            lo, hi = float("inf"), float("-inf")
            for data in plot_axis.data.values():
                arr = data[data_idx]
                if arr is None:
                    continue
                try:
                    finite = arr[np.isfinite(arr)]
                    if len(finite):
                        lo = min(lo, float(finite.min()))
                        hi = max(hi, float(finite.max()))
                except Exception:
                    pass
            if lo > hi:
                return None
            pad = _MARGIN * (hi - lo)
            return (lo - pad, hi + pad)

        x_bounds = None
        if self._left() and len(self._left()):
            x_bounds = axis_bounds(self._left(), 0)
        if x_bounds is None and self._right() and len(self._right()):
            x_bounds = axis_bounds(self._right(), 0)
        # Merge x bounds from both axes (they share x)
        if self._right() and len(self._right()):
            rb = axis_bounds(self._right(), 0)
            if rb is not None:
                x_bounds = (
                    min(x_bounds[0], rb[0]) if x_bounds else rb[0],
                    max(x_bounds[1], rb[1]) if x_bounds else rb[1],
                ) if x_bounds else rb

        y1_bounds = axis_bounds(self._left(), 1) if self._left() and len(self._left()) else None
        y2_bounds = None
        if not self._3D and self._right() and len(self._right()):
            y2_bounds = axis_bounds(self._right(), 1)

        return x_bounds, y1_bounds, y2_bounds

    def on_scroll(self, event):
        """Zoom in/out with the scroll wheel, respecting locked axes and data bounds.
        Uses a separate zoom-limit state so the user-facing Lock buttons are never touched."""
        if event.inaxes is None or self._3D:
            return

        # scroll up (step > 0) → zoom in; scroll down → zoom out
        zoom = 0.85 ** event.step  # < 1 means zoom in when step > 0

        x_bounds, y1_bounds, y2_bounds = self._get_data_bounds()

        def zoomed(cur_lo, cur_hi, focus, factor, bounds):
            """Zoom [cur_lo, cur_hi] by factor around focus, clamped to bounds."""
            new_lo = focus - (focus - cur_lo) * factor
            new_hi = focus + (cur_hi - focus) * factor
            if bounds is not None:
                b_lo, b_hi = bounds
                b_range = b_hi - b_lo
                # Clamp zoom-out: never exceed the full data range
                new_lo = max(new_lo, b_lo)
                new_hi = min(new_hi, b_hi)
                # Clamp zoom-in: prevent collapsing to a point
                min_range = b_range * 1e-5 if b_range > 0 else 1e-9
                if new_hi - new_lo < min_range:
                    mid = (new_lo + new_hi) / 2.0
                    new_lo = mid - min_range / 2.0
                    new_hi = mid + min_range / 2.0
            return [new_lo, new_hi]

        ax_left = self._left().axis()

        # ── X axis (skip if user has it locked) ──────────────────────────────
        if not self.x_locked.isChecked():
            cur = self._zoom_x_limits if self._zoom_x_limits is not None else list(ax_left.get_xlim())
            x_focus = event.xdata if event.xdata is not None else (cur[0] + cur[1]) / 2.0
            self._zoom_x_limits = zoomed(cur[0], cur[1], x_focus, zoom, x_bounds)

        # ── Y1 axis (skip if user has it locked) ─────────────────────────────
        if self._left() and len(self._left()) and not self.y1_locked.isChecked():
            if event.inaxes is ax_left:
                y1_focus = event.ydata
            else:
                try:
                    disp = self._right().axis().transData.transform((event.xdata or 0, event.ydata or 0))
                    y1_focus = ax_left.transData.inverted().transform(disp)[1]
                except Exception:
                    y1_focus = None
            cur = self._zoom_y1_limits if self._zoom_y1_limits is not None else list(ax_left.get_ylim())
            focus = y1_focus if y1_focus is not None else (cur[0] + cur[1]) / 2.0
            self._zoom_y1_limits = zoomed(cur[0], cur[1], focus, zoom, y1_bounds)

        # ── Y2 axis (skip if user has it locked) ─────────────────────────────
        if not self._3D and self._right() and len(self._right()) and not self.y2_locked.isChecked():
            ax_right = self._right().axis()
            if event.inaxes is ax_right:
                y2_focus = event.ydata
            else:
                try:
                    disp = ax_left.transData.transform((event.xdata or 0, event.ydata or 0))
                    y2_focus = ax_right.transData.inverted().transform(disp)[1]
                except Exception:
                    y2_focus = None
            cur = self._zoom_y2_limits if self._zoom_y2_limits is not None else list(ax_right.get_ylim())
            focus = y2_focus if y2_focus is not None else (cur[0] + cur[1]) / 2.0
            self._zoom_y2_limits = zoomed(cur[0], cur[1], focus, zoom, y2_bounds)

        self.draw()

    @staticmethod
    def _point_to_segment_dist(px, py, x1, y1, x2, y2):
        """Return (distance, t) where t in [0,1] is the parameter along segment (x1,y1)-(x2,y2)
        closest to point (px, py). All coordinates are in pixel space."""
        dx, dy = x2 - x1, y2 - y1
        seg_len_sq = dx * dx + dy * dy
        if seg_len_sq < 1e-12:
            t = 0.0
        else:
            t = max(0.0, min(1.0, ((px - x1) * dx + (py - y1) * dy) / seg_len_sq))
        cx = x1 + t * dx
        cy = y1 + t * dy
        dist = math.hypot(px - cx, py - cy)
        return dist, t

    def on_hover(self, event):
        """Show a tooltip with parameter name and x/y value when hovering near a plot line."""
        if event.inaxes is None or self._3D:
            self._hover_tooltip.hide()
            self._clear_hover_marker()
            return

        # Collect all axes we need to check
        axes_to_check = []
        if self._left() is not None and len(self._left()):
            axes_to_check.append((self._left(), "y1"))
        if self._right() is not None and len(self._right()):
            axes_to_check.append((self._right(), "y2"))

        best_dist = self._hover_pick_radius
        best_label = None
        best_x = None
        best_y = None
        best_ax = None
        best_color = None

        for plot_axis, side in axes_to_check:
            ax = plot_axis.axis()
            for label, line in plot_axis.plots.items():
                try:
                    xy_data = line.get_xydata()
                except Exception:
                    continue
                if len(xy_data) < 1:
                    continue
                # Filter out NaN rows
                valid = np.isfinite(xy_data[:, 0]) & np.isfinite(xy_data[:, 1])
                xy_valid = xy_data[valid]
                if len(xy_valid) == 0:
                    continue
                # Transform data coords → display (pixel) coords
                try:
                    xy_pixels = ax.transData.transform(xy_valid)
                except Exception:
                    continue

                mx, my = event.x, event.y

                # Check every segment between consecutive valid points
                for i in range(len(xy_pixels) - 1):
                    x1p, y1p = xy_pixels[i]
                    x2p, y2p = xy_pixels[i + 1]
                    dist, t = self._point_to_segment_dist(mx, my, x1p, y1p, x2p, y2p)
                    if dist < best_dist:
                        best_dist = dist
                        best_label = label
                        # Interpolate data-space coordinates along the segment
                        best_x = xy_valid[i, 0] * (1 - t) + xy_valid[i + 1, 0] * t
                        best_y = xy_valid[i, 1] * (1 - t) + xy_valid[i + 1, 1] * t
                        best_ax = ax
                        best_color = line.get_color()

                # Also check the last point in isolation (single-point lines)
                if len(xy_pixels) == 1:
                    dist = math.hypot(xy_pixels[0, 0] - mx, xy_pixels[0, 1] - my)
                    if dist < best_dist:
                        best_dist = dist
                        best_label = label
                        best_x = xy_valid[0, 0]
                        best_y = xy_valid[0, 1]
                        best_ax = ax
                        best_color = line.get_color()

        if best_label is None:
            self._hover_tooltip.hide()
            self._clear_hover_marker()
            return

        # Draw / update the marker dot
        self._update_hover_marker(best_ax, best_x, best_y, best_color)

        # Format tooltip text
        tooltip_text = "<b>{}</b><br/>x: {:.4g}<br/>y: {:.4g}".format(
            best_label, best_x, best_y
        )
        self._hover_tooltip.setText(tooltip_text)
        self._hover_tooltip.setTextFormat(QtCore.Qt.RichText)
        self._hover_tooltip.adjustSize()

        # Position tooltip near cursor, keeping it inside the canvas
        canvas_w = self.canvas.width()
        canvas_h = self.canvas.height()
        tip_w = self._hover_tooltip.width()
        tip_h = self._hover_tooltip.height()
        # event.x / event.y are in display coords (y flipped vs widget coords)
        wx = int(event.x) + 12
        wy = canvas_h - int(event.y) - tip_h - 8
        if wx + tip_w > canvas_w:
            wx = int(event.x) - tip_w - 12
        if wy < 0:
            wy = canvas_h - int(event.y) + 8
        self._hover_tooltip.move(wx, wy)
        self._hover_tooltip.show()
        self._hover_tooltip.raise_()

    def _update_hover_marker(self, ax, x, y, color):
        """Draw or move the hover marker dot on the given axes."""
        marker = getattr(self, "_hover_marker", None)
        if marker is not None and marker.axes is ax:
            marker.set_data([x], [y])
        else:
            self._clear_hover_marker()
            (marker,) = ax.plot(
                [x],
                [y],
                marker="o",
                markersize=8,
                color=color,
                markeredgecolor="white",
                markeredgewidth=1.5,
                linestyle="None",
                zorder=10,
                animated=False,
            )
            self._hover_marker = marker
        self.canvas.draw_idle()

    def _clear_hover_marker(self):
        """Remove the hover marker dot if present."""
        marker = getattr(self, "_hover_marker", None)
        if marker is not None:
            try:
                marker.remove()
            except Exception:
                pass
            self._hover_marker = None
            self.canvas.draw_idle()

    def on_draw(self, event):
        if self.x_limits is not None:
            self.x_limits = self._left().axis().get_xlim()
        elif self._zoom_x_limits is not None:
            self._zoom_x_limits = list(self._left().axis().get_xlim())
        if self.y1_limits is not None:
            self.y1_limits = self._left().axis().get_ylim()
        elif self._zoom_y1_limits is not None:
            self._zoom_y1_limits = list(self._left().axis().get_ylim())
        if self.y2_limits is not None:
            if self._3D:
                self.y2_limits = self._left().axis().get_zlim()
            else:
                self.y2_limits = self._right().axis().get_ylim()
        elif not self._3D and self._zoom_y2_limits is not None:
            self._zoom_y2_limits = list(self._right().axis().get_ylim())

    def draw(self):
        # Merge: user lock limits take priority, zoom limits act as soft override,
        # None falls through to PlotFigure auto-fit.
        x_lim = self.x_limits if self.x_limits is not None else self._zoom_x_limits
        y1_lim = self.y1_limits if self.y1_limits is not None else self._zoom_y1_limits
        y2_lim = self.y2_limits if self.y2_limits is not None else self._zoom_y2_limits
        PlotFigure.draw(self, x_lim, y1_lim, y2_lim)
        self.canvas.draw()

    def _y_lock_changed(self, name, cbox, get_lim):
        if cbox.isChecked():
            cbox.setText("Unlock {}".format(name))
            return get_lim()
        else:
            cbox.setText("Lock {}".format(name))
            return None

    def x_locked_changed(self, status):
        self.x_limits = self._y_lock_changed(
            "X", self.x_locked, self._left().axis().get_xlim
        )
        if self.x_limits is None:
            self._zoom_x_limits = None  # clear scroll-zoom so auto-fit resumes
            self.draw()

    def y1_locked_changed(self, status):
        self.y1_limits = self._y_lock_changed(
            "Y1", self.y1_locked, self._left().axis().get_ylim
        )
        if self.y1_limits is None:
            self._zoom_y1_limits = None
            self.draw()

    def y2_locked_changed(self, status):
        if self._3D:
            self.y2_limits = self._y_lock_changed(
                "Z", self.y2_locked, self._left().axis().get_zlim
            )
        else:
            self.y2_limits = self._y_lock_changed(
                "Y2", self.y2_locked, self._right().axis().get_ylim
            )
        if self.y2_limits is None:
            self._zoom_y2_limits = None
            self.draw()

    def align_y1_y2_released(self):
        i0, iN = self.getFrameRange()
        if i0 == iN:
            return
        if self.y1_limits is None:
            self.y1_limits = self._left().getLimits(i0, iN, 1)
        if self.y2_limits is None:
            self.y2_limits = self._right().getLimits(i0, iN, 1)
        if self.y1_limits is not None and self.y2_limits is not None:
            self.y1_locked.setChecked(True)
            self.y2_locked.setChecked(True)
            self.y1_limits = (
                min(self.y1_limits[0], self.y2_limits[0]),
                max(self.y1_limits[1], self.y2_limits[1]),
            )
            self.y2_limits = self.y1_limits
            self.draw()

    def show(self):
        return QWidget.show(self)
