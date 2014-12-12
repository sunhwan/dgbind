from colvars import Colvars

class Controller:
    def __init__(self, conf):
        colvars = Colvars(conf)

        if conf['simulations']['method'] == 'REMD':
            from REMD.controller import Controller
            self.controller = Controller(conf)
            self.controller.colvars = colvars

    def create(self):
        self.controller.create()