

class Bridge():
    def __init__(self, *datawidgets):
        self.functions = []
        self.functions_size = []
        self.variables = []
        self.variable_size = []
        self.variable_flags = []

        for dw in datawidgets:
            vnames, vsizes, vflags = dw.variables()
            self.variables.extend(vnames)
            self.variable_size.extend(vsizes)
            self.variable_flags.extend(vflags)

            fnames, fsize = dw.functions()
            self.functions.extend(fnames)
            self.functions_size.extend(fsize)

    def generate_jacobian(self):
        ...

    def generate_fobj(self):
        ...

    def weights(self):
        ...
