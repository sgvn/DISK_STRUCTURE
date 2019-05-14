import sys
from PyQt5.QtWidgets import QApplication, QWidget, QMainWindow, QLabel, QToolTip, QPushButton
from PyQt5.QtGui import QPainter, QColor, QPen, QIcon, QFont  
from PyQt5.QtCore import Qt
import random
 
class App(QWidget):
    def __init__(self):
        super().__init__()
        self.title = 'PyGrain interface'
        self.left = 400
        self.top = 250
        self.width = 640
        self.height = 480
        self.initUI()
 
    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        #button
        QToolTip.setFont(QFont('Tahoma', 10))
        self.setToolTip('Run all the models')
        btn = QPushButton('RUN', self)
        btn.setToolTip('Run all the models')
        btn.resize(btn.sizeHint())
        btn.move(400, 400)       

        # Set window background color
        self.setAutoFillBackground(True)
        p = self.palette()
        p.setColor(self.backgroundRole(), QColor("#71A7A7"))
        self.setPalette(p)
        
        # Add paint widget and paint
        self.m = PaintWidget(self)
        self.m.move(0,0)
        self.m.resize(self.width,self.height)
        self.show()


class PaintWidget(QWidget):
    def paintEvent(self, event):
        qp = QPainter(self)
        
        #qp.setPen(Qt.black)
        #size = self.size()
        
        # Colored rectangles
        #qp.setBrush(QColor("#000"))
        #qp.drawRect(0, 0, 100, 100)
        
        
 
if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())