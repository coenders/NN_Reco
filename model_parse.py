#!/usr/bin/env python
# coding: utf-8

"""This file is part of DeepIceLearning
DeepIceLearning is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import tables
import numpy as np
import importlib
import os
import sys
from collections import OrderedDict


def prepare_io_shapes(inputs, outputs, exp_file):
    inp_transformations = OrderedDict()
    inp_shapes = OrderedDict()
    out_transformations = OrderedDict()
    out_shapes = OrderedDict()
    # open example file
    inp_file = tables.open_file(exp_file)

    for br in inputs:
        inp_shapes[br] = {}
        inp_transformations[br] = {}
        for var, tr in zip(inputs[br]["variables"],
                           inputs[br]["transformations"]):
            test_arr = np.array(inp_file._get_node("/" + var)[0])
            # eval("inp_file.root.{}".format(var))[0]
            res_shape = np.shape(tr(test_arr))
            inp_shapes[br][var] = res_shape
            inp_transformations[br][var] = tr
        inp_shapes[br]["general"] = res_shape[:-1] + (len(inputs[br]["variables"]),)

    for br in outputs:
        out_shapes[br] = {}
        out_transformations[br] = {}
        for var, tr in zip(outputs[br]["variables"],
                           outputs[br]["transformations"]):
            test_arr = np.array(inp_file._get_node("/reco_vals").col(var)[0])
            # eval("inp_file.root.reco_vals.cols.{}".format(var))[0])
            tr_applied = tr(test_arr)
            res_shape = tr_applied.shape if not isinstance(tr_applied,
                                                           np.float) else 1
            out_shapes[br][var] = res_shape
            out_transformations[br][var] = tr
        out_shapes[br]["general"] = len(outputs[br]["variables"])
    return inp_shapes, inp_transformations, out_shapes, out_transformations


def parse_functional_model(cfg_file, exp_file):
    try:
        # fancy relative imports..
        sys.path.append(os.path.dirname(cfg_file))
        mname = os.path.splitext(os.path.basename(cfg_file))[0]
        func_model_def = importlib.import_module(mname)
        sys.path.pop()
    except Exception:
        raise Exception('Import of model.py failed: {}'.format(cfg_file))
    inputs = func_model_def.inputs
    outputs = func_model_def.outputs

    in_shapes, in_trans, out_shapes, out_trans = \
        prepare_io_shapes(inputs, outputs, exp_file)
    print out_shapes
    print in_shapes
    model = func_model_def.model(in_shapes, out_shapes) 
    return model, in_shapes, in_trans, out_shapes, out_trans