"""Output the refinement process internals.

This module provides output routines in two different ways:
(1) html output or (2) text output. Inside each one, three levels are
supplied: (1) for the start of the fitting, (2) for each iteration
and (3) for the final convergence achievement.
"""

import datetime

import numpy as np

import consts

# variable that controls the amount of output to be produced
# values accepted are 0-2.
default_verbosity = 1


def iteration(values, increment):
    return "\n".join(f"{val} {inc}" for val, inc in zip(values, increment))


def html_start(**kwargs):
    """Output the initial state prior to fitting start.

    Parameters:
        titration (sequence, optional): the raw data
        verbosity (int, optional): A number 0-2 representing the amount of
             verbosity for the output.
    """
    template = """<hr><h2>Fitting potentiometric data</h2>
    <h3>Fitting Parameters</h3> {}
    <h3>Data Input</h3> {}
    <p>Starting fit on potentiometric data: {}</p>
    """
    timestamp = ("Starting fit on potentiometric data: " +
                 str(datetime.datetime.now()))

    verbosity = kwargs.get('verbosity', default_verbosity)

    # TODO avoid importing numpy here
    text = "<h3>Fitting potentiometric data</h3>"
    if verbosity >= 1:
        timestamp = ("Starting fit on potentiometric data: " +
                     str(datetime.datetime.now()))
    else:
        timestamp = ""

    if verbosity == 2 and 'titration' in kwargs:
        emf = (_['emf'] for _ in kwargs['titration'])
        output_data = ('<p>Data input: %d datasets. ' % len(kwargs['titration']) +
                       "; ".join("Set #%d: %d points (%d ignored)" %
                                 (n, d.size, np.count_nonzero(d.mask))
                                 for n, d in enumerate(emf)) + '</p>')
    else:
        output_data = ""

    return text + timestamp + output_data


def html_lm_iteration(iteration, xi, dx, chisq, verbosity=default_verbosity):
    """Produce html output with information for this iteration.

    Parameters:
        iteration (int): The number of the iteration
        xi (:class:`numpy.ndarray`): A 1D array of floats with the initial
            values for the constants at the beginning of the iteration
        dx (:class:`numpy.ndarray`): A 1D array of floats with the variation
            in the values for the constants in this iteration
        chisq (float): The χ² value for this iteration
        verb (int, optional): A number 0-2 representing the amount of
             verbosity for the output.
    Return:
        str: HTML code
    """
    row_txt = "<tr>" + 3*"<td>{:10.3f}</td>" + "</tr>"
    rows = "\n".join(row_txt.format(bi, i, bo)
                     for bi, i, bo in zip(xi, dx, (xi+dx)))
    txtout = ("<h3>ITERATION #{}</h3>\n"
              "<table style=\"border: 1px solid black;\"><thead><tr><th>in log&beta;</th>"
              "<th>&Delta;log&beta;</th><th>out log&beta;</th></tr></thead>"
              "<tbody>{}</tbody>"
              "\n</table><p>&chi;<sup>2</sup> = {:.2f}</p>").format(iteration, rows, chisq)
    return txtout


def html_nm_iteration(iteration, x, chisq, verbosity=default_verbosity):
    """Produce html output with information for this iteration.

    Parameters:
        iteration (int): The number of the iteration
        xi (:class:`numpy.ndarray`): A 1D array of floats with the initial
            values for the constants at the beginning of the iteration
        chisq (float): The χ² value for this iteration
        verbosity (int, optional): A number 0--2 representing the amount of
             verbosity for the output.
    Return:
        str: HTML code
    """
    # import consts
    # txtout = ("<h3>ITERATION #{}</h3>"
    #           "<table><tr><th>parameter</th></tr>").format(iteration)
    # for b in x:
    #     txtout += "<tr><td>{:.3f}</td></tr>".format(b/consts.LOGK)
    # txtout += "</table><br>chisq = {:.2f}".format(chisq)
    row_txt = "<tr><td>{:.3f}</td></tr>"
    rows = "\n".join(row_txt.format(_x) for _x in x)
    txtout = ("<h3>ITERATION #{}</h3>\n"
              "<table>\n<tr><th>parameter</th></tr>\n{}"
              "\n</table><p>chisq = {:.2f}</p>").format(iteration, rows, chisq)
    return txtout


def html_finalconvg(params, errors, verbosity=default_verbosity, **kwargs):
    # TODO Use a different technique to construct HTML output
    # Do NOT use + to concatenate strings (slooooow)
    if 'iterations' in kwargs:
        ins = " in {} iterations".format(kwargs['iterations'])
    else:
        ins = ""
    t = "<p><i>Convergence achieved{}</i></p><h2>Final Parameters</h2>".format(ins)
    t += "<table><tr><th>#</th><th>log&beta;</th><th>error</th></tr>"
    for i, (b, e) in enumerate(zip(params, errors)):
        t += "<tr><td>{}</td>".format(i) + \
                  2*"<td>%.3f</td>" % (b, e) + "</tr>"
    t += "</table>"

    if verbosity >= 1 and 'correlation' in kwargs:
        correlation_matrix = kwargs['correlation']
        t += "<h2>Correlation matrix</h2>"
        t += '<table>\n<tr>{}</tr>\n</table>'.format(
            '</tr>\n<tr>'.join(
                ''.join('<td>{:.3f}</td>'.format(i) for i in row)
                for row in correlation_matrix))

    return t  # '<p>' + t + '</p>'


def spyq_start(beta, bflags, stoichiometry, titrations, weights,
               verbosity=default_verbosity):
    def _flag(x):
        if x == 0:
            return 'CONSTANT'
        elif x == 1:
            return 'REFINE'
        else:
            return 'CONSTRAINT %d' % (x-1)
    # raise NotImplementedError
    E = len(stoichiometry)
    S = len(stoichiometry[0])
    print("(i) model")
    print("lg beta   stoichiometry")
    for i in range(E):
        print("%6.2f" % beta[i], end="")
        print("".join(["%5d" % stoichiometry[i][j] for j in range(S)]),
              end="   ")
        print(_flag(bflags[i]))

    print("\n(ii) data")
    V = (d['V'] for d in titrations)
    emf = [d['emf'] for d in titrations]
    # import pudb
    # pudb.set_trace()
    point = 0
    for i, (v, e, w) in enumerate(zip(V, emf, weights)):
        print("curve ", i+1)
        if verbosity == 0:
            print(len(v), " points")
        else:
            print("point     titre     emf        weight")

            for n, it in enumerate(zip(v, e, w)):
                if verbosity == 2 or ((verbosity == 1) and (n % 5 == 0)):
                    print("%6d %8.3f %10.2f %8.2f" % (point, *it))
                point += 1
        print()


def spyq_lm_iteration(iteration, x_start, x_increment, chisq,
                      verbosity=default_verbosity):
    """Print output for :program:`supyquad`.

    Parameters:
        iteration (int): The number of the iteration
        x_start (:class:`numpy.ndarray`): A 1D array of floats with the initial
            values for the constants at the beginning of the iteration
        x_increment (:class:`numpy.ndarray`): A 1D array of floats with the
            variation in the values for the constants in this iteration
        chisq (float): The χ² value for this iteration
        verbosity (int, optional): A number 0-2 representing the amount of
             verbosity for the output.
    """
    def _conv_(array):
        return "".join("{:8.2f}".format(i/consts.LOGK) for i in array)

    print("ITERATION #", iteration)
    print("in param =  ", _conv_(x_start))
    print("increment = ", _conv_(x_increment))
    print("out param = ", _conv_(x_start+x_increment))
    print("chisq = %.2f" % chisq)
    print()


def spyq_nm_iteration(iteration, x, chisq, verbosity=1):
    """Print output for :program:`supyquad`.

    Parameters:
        iteration (int): The number of the iteration
        x_start (:class:`numpy.ndarray`): A 1D array of floats with the initial
            values for the constants at the beginning of the iteration
        chisq (float): The χ² value for this iteration
        verbosity (int, optional): A number 0-2 representing the amount of
             verbosity for the output.
    """

    def _conv_(array):
        return "".join("{:8.2f}".format(i/consts.LOGK) for i in array)

    print("ITERATION #", iteration)
    print("best params =  ", _conv_(x))
    print("chisq = %.2f" % chisq)
    print()


def spyq_finalconvg(params, errors, verbosity=1, **kwargs):
    """Print final stated for :program:`supyquad`.

    Parameters:
        params (:class:`numpy.ndarray`): Refined parameters' final values.
        errors (:class:`numpy.ndarray`): Errors for the refined parameters.
        verbosity (int, optional): A number 0-2 representing the amount of
             verbosity for the output.
        correlation_matrix (:class:`numpy.ndarray`, optional): Correlation
            matrix.
    """
    print("CONVERGENCE ACHIEVED\n")

    # import numpy as np

    if verbosity > 0:
        titrations = kwargs['titrations']
        residuals = kwargs['residuals']
        accumlen = 0
        for n_curve, (curve, resid) in enumerate(zip(titrations, residuals)):
            lenc = len(resid)
            print("Curve ", n_curve + 1)
            print("Sigma:", np.sqrt(np.sum(resid**2)/lenc))
            print("point   EMF     residual")
            for j, (y, r) in enumerate(zip(curve['emf'], resid)):
                print("%3d   %8.2f %8.2f" % (1+j+accumlen, y, r))
            print()
            accumlen += lenc

    print("Final parameters:")
    print("param.\tlogB\terror")
    for count, (param, error) in enumerate(zip(params, errors)):
        print("%d\t%6.3f\t%6.3f" % (count, param, error))
    print()
    if verbosity == 1 and 'correlation' in kwargs:
        correlation_matrix = kwargs['correlation']
        print("Correlation matrix")
        print('\n'.join(''.join(['{:8.3f}'.format(i) for i in b])
                        for b in correlation_matrix))


report_function = {
    consts.METHOD_LM: html_lm_iteration,
    consts.METHOD_NM: html_nm_iteration
}
