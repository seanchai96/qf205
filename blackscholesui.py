import sys 
from PyQt5.QtWidgets import QMainWindow, QApplication
from PyQt5 import uic
from blackscholespde import blackscholes

qtCreatorFile = "blackscholes.ui"
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)

class Main(QMainWindow, Ui_MainWindow):   
    def __init__(self):
        super().__init__()
        self.setupUi(self)
        self.explicit_fdm.setChecked(True)
        self.Calculate.clicked.connect(self.blackscholescalc)

    def blackscholescalc(self):
        S = float(self.S.text())
        t = float(self.t.text())   
        K = float(self.K.text())
        r = float(self.r.text())
        sigma = float(self.sigma.text())
        q = float(self.q.text())
        M = int(self.M.text())
        N = int(self.N.text())
             
        method = blackscholes(S, t, K, r, sigma, q)

        if self.explicit_fdm.isChecked():
            result = method.explicitfdm(M, N)
        else:
            result = method.implicitfdm(M, N)
        
        self.call.setText(str(result[0]))
        self.put.setText(str(result[1]))

if __name__ == '__main__':
    app = QApplication(sys.argv)
    main = Main()
    main.show()
    sys.exit(app.exec_())