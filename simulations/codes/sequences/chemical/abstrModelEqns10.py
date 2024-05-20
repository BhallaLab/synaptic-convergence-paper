"""
This code has been adapted from https://github.com/ModelDBRepository/227318

It was originally a part of Upinder Singh Bhalla (2017) Synaptic input sequence discrimination on behavioral timescales mediated by reaction-diffusion chemistry in dendrites eLife 6:e25827

https://doi.org/10.7554/eLife.25827

"""

import re
import moose


def parseExpr(expr, params, hasCa):
    if hasCa:
        expr = expr.replace("Ca", "x0")
        expr = expr.replace("A", "x1")
        expr = expr.replace("B", "x2")
    else:
        expr = expr.replace("Ca", "x0")  # happens in the negFF model
        expr = expr.replace("A", "x0")  # This is the usual case.
        expr = expr.replace("B", "x1")

    parts = re.split("k", expr)
    ret = parts[0]
    for i in parts[1:]:
        ret += str(params["k" + i[:2]])
        ret += i[2:]

    if hasCa:
        return "x3*( " + ret + ")"
    else:
        return "x2*( " + ret + ")"


def makeChemProto(name, Aexpr, Bexpr, params):
    sw = params["stimWidth"]
    diffLength = params["diffusionLength"]
    dca = params["diffConstA"] * diffLength * diffLength
    dcb = params["diffConstB"] * diffLength * diffLength

    # Objects
    chem = moose.Neutral("/library/" + name)
    compt = moose.CubeMesh("/library/" + name + "/" + name)
    A = moose.Pool(compt.path + "/A")
    B = moose.Pool(compt.path + "/B")
    Z = moose.BufPool(compt.path + "/Z")
    Ca = moose.BufPool(compt.path + "/Ca")
    phase = moose.BufPool(compt.path + "/phase")
    extra = moose.BufPool(compt.path + "/extra")
    ampl = moose.BufPool(compt.path + "/ampl")
    Adot = moose.Function(A.path + "/Adot")
    Bdot = moose.Function(B.path + "/Bdot")
    CaStim = moose.Function(Ca.path + "/CaStim")
    A.diffConst = dca
    B.diffConst = dcb

    # Equations

    Adot.expr = parseExpr(Aexpr, params, True)
    Bdot.expr = parseExpr(Bexpr, params, False)
    s = "exp( -((x{0} - t)^2)/(2*" + str(sw * sw) + ") )"
    # CaStim.expr = 'x2 * exp( -((x0 - t)^2)/(2* ' + str(sw*sw) + ') )'
    CaStim.expr = "x2 * ({0} + {1})".format(s.format(0), s.format(1))

    # print Adot.expr
    # print Bdot.expr
    # print CaStim.expr

    # Connections
    Adot.x.num = 4
    moose.connect(Ca, "nOut", Adot.x[0], "input")
    moose.connect(A, "nOut", Adot.x[1], "input")
    moose.connect(B, "nOut", Adot.x[2], "input")
    moose.connect(Z, "nOut", Adot.x[3], "input")
    moose.connect(Adot, "valueOut", A, "increment")

    Bdot.x.num = 3
    if name[:5] == "negFF":
        moose.connect(Ca, "nOut", Bdot.x[0], "input")
        print("Doing special msg")
    else:
        moose.connect(A, "nOut", Bdot.x[0], "input")
    moose.connect(B, "nOut", Bdot.x[1], "input")
    moose.connect(Z, "nOut", Bdot.x[2], "input")
    moose.connect(Bdot, "valueOut", B, "increment")

    CaStim.x.num = 3
    moose.connect(phase, "nOut", CaStim.x[0], "input")
    moose.connect(extra, "nOut", CaStim.x[1], "input")
    moose.connect(ampl, "nOut", CaStim.x[2], "input")
    moose.connect(CaStim, "valueOut", Ca, "setN")

    return compt


def makeBis(args):
    params = {
        "k0a": 0.1,  # Constant
        "k1a": -5.0,  # Coeff for A
        "k2a": 5.0,  # Coeff for A^2
        "k3a": -1.0,  # Coeff for A^3
        "k4a": 10.0,  # turnon of A by A and Ca
        "k5a": -5.0,  # Turnoff of A by B
        "k1b": 0.01,  # turnon of B by A
        "k2b": -0.01,  # Decay rate of B
        "diffusionLength": 1.0e-6,  # Diffusion characteristic length, used as voxel length too.
        "dendDiameter": 10e-6,  # Diameter of section of dendrite in model
        "dendLength": 100e-6,  # Length of section of dendrite in model
        "diffConstA": 5.0,  # Diffusion constant of A
        "diffConstB": 2.0,  # Diffusion constant of B
        "stimWidth": 1.0,  # Stimulus width in seconds
        "stimAmplitude": 1.0,  # Stimulus amplitude, arb units. From FHN review
        "blankVoxelsAtEnd": 10,  # of voxels to leave blank at end of cylinder
        "preStimTime": 10.0,  # Time to run before turning on stimulus.
        "postStimTime": 40.0,  # Time to run after stimulus. ~3x decay time
        "settleTime": 20.0,  # Settling time of response, after stimulus.
        # To include in analysis of total response over
        # entire dendrite.
        "fnumber": 1,  # Number to append to fname
    }

    for i in args:
        params[i] = args[i]

    makeChemProto(
        "bis",
        "k0a + k1a*A + k2a*A*A + k3a*A*A*A + k4a*Ca*A/(1+A+10*B) + k5a*A*B",
        "k1b*A*A + k2b*B",
        params,
    )
    return params
