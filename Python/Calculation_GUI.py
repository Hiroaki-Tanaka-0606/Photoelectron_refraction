# Photoelectron_refraction
# GUI

from pyqtgraph.Qt import QtGui, QtCore, QtWidgets
import pyqtgraph as pg

import numpy as np
import math
import lib
import h5py
from datetime import datetime

import Config

class MainWindow(QtGui.QMainWindow):
    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        self.setWindowTitle("Photoelectron refraction")

        font=QtGui.QFont()
        font.setFamilies(Config.fontFamilies)
        font.setPixelSize(Config.fontSize_normal)

        bFont=QtGui.QFont(font)
        bFont.setBold(True)

        vbox=QtGui.QVBoxLayout()
        vbox.setContentsMargins(*Config.ContentsMargins)
        vbox.setAlignment(QtCore.Qt.AlignTop)

        mainWidget=QtGui.QWidget()
        mainWidget.setLayout(vbox)
        self.setCentralWidget(mainWidget)

        # Row 1: potential
        row1=QtGui.QHBoxLayout()
        row1.setAlignment(QtCore.Qt.AlignLeft)
        vbox.addLayout(row1)

        label1A=QtGui.QLabel("Work function (eV)")
        label1A.setFont(bFont)
        row1.addWidget(label1A)

        self.Wtext=QtGui.QLineEdit("5")
        row1.addWidget(self.Wtext)

        label1B=QtGui.QLabel("Inner potential (eV)")
        label1B.setFont(bFont)
        row1.addWidget(label1B)

        self.V0text=QtGui.QLineEdit("12")
        row1.addWidget(self.V0text)

        # Row 2: initial state
        row2=QtGui.QHBoxLayout()
        row2.setAlignment(QtCore.Qt.AlignLeft)
        vbox.addLayout(row2)

        label2A=QtGui.QLabel("Parabola paramter")
        label2A.setFont(bFont)
        row2.addWidget(label2A)

        self.atext=QtGui.QLineEdit("-1")
        row2.addWidget(self.atext)

        label2B=QtGui.QLabel("Energy of the top from EF (eV)")
        label2B.setFont(bFont)
        row2.addWidget(label2B)
        self.V1text=QtGui.QLineEdit("0")
        row2.addWidget(self.V1text)

        label2C=QtGui.QLabel("Coordinate of the top (Ang^-1)")
        label2C.setFont(bFont)
        row2.addWidget(label2C)

        self.k0xtext=QtGui.QLineEdit("0")
        row2.addWidget(self.k0xtext)
        self.k0ytext=QtGui.QLineEdit("0")
        row2.addWidget(self.k0ytext)
        self.k0ztext=QtGui.QLineEdit("10")
        row2.addWidget(self.k0ztext)

        # Row 3: reciprocal space
        row3=QtGui.QHBoxLayout()
        row3.setAlignment(QtCore.Qt.AlignLeft)
        vbox.addLayout(row3)

        self.kPlaneChoice=QtGui.QButtonGroup()
        self.kFlat=QtGui.QRadioButton("Flat k plane")
        self.kFlat.setChecked(True)
        self.kFlat.setFont(bFont)
        row3.addWidget(self.kFlat)
        self.kPlaneChoice.addButton(self.kFlat)
        label3A=QtGui.QLabel("( kz = ")
        row3.addWidget(label3A)
        self.kFlat_kz=QtGui.QLineEdit("10")
        row3.addWidget(self.kFlat_kz)
        label3B=QtGui.QLabel(" Ang^-1)")
        row3.addWidget(label3B)

        self.kCurved=QtGui.QRadioButton("Curved k plane")
        self.kCurved.setFont(bFont)
        row3.addWidget(self.kCurved)
        self.kPlaneChoice.addButton(self.kCurved)
        label3C=QtGui.QLabel("( |k| = ")
        row3.addWidget(label3C)
        self.kCurved_k=QtGui.QLineEdit()
        row3.addWidget(self.kCurved_k)
        label3D=QtGui.QLabel(" Ang^-1)")
        row3.addWidget(label3D)

        # Row 4: surface
        row4=QtGui.QHBoxLayout()
        row4.setAlignment(QtCore.Qt.AlignLeft)
        vbox.addLayout(row4)

        self.surfaceChoice=QtGui.QButtonGroup()
        self.surfaceConst=QtGui.QRadioButton("Constant surface")
        self.surfaceConst.setFont(bFont)
        self.surfaceConst.setChecked(True)
        row4.addWidget(self.surfaceConst)
        self.surfaceChoice.addButton(self.surfaceConst)
        label4A=QtGui.QLabel("( theta, phi = ")
        row4.addWidget(label4A)
        self.surfaceConst_theta=QtGui.QLineEdit("0")
        row4.addWidget(self.surfaceConst_theta)
        self.surfaceConst_phi=QtGui.QLineEdit("0")
        row4.addWidget(self.surfaceConst_phi)
        label4B=QtGui.QLabel(" (deg))")
        row4.addWidget(label4B)

        self.surfaceRandom=QtGui.QRadioButton("Random surface")
        self.surfaceRandom.setFont(bFont)
        row4.addWidget(self.surfaceRandom)
        self.surfaceChoice.addButton(self.surfaceRandom)
        label4C=QtGui.QLabel("( samples = ")
        row4.addWidget(label4C)
        self.surfaceRandom_samples=QtGui.QLineEdit()
        row4.addWidget(self.surfaceRandom_samples)
        label4D=QtGui.QLabel(")")
        row4.addWidget(label4D)

        # Row 5: grid
        row5=QtGui.QHBoxLayout()
        row5.setAlignment(QtCore.Qt.AlignLeft)
        vbox.addLayout(row5)

        label5A=QtGui.QLabel("kx range")
        label5A.setFont(bFont)
        row5.addWidget(label5A)

        label5B=QtGui.QLabel("Min and Max (Ang^-1)")
        row5.addWidget(label5B)
        self.kxMintext=QtGui.QLineEdit("-0.5")
        row5.addWidget(self.kxMintext)
        self.kxMaxtext=QtGui.QLineEdit("0.5")
        row5.addWidget(self.kxMaxtext)
        label5C=QtGui.QLabel("Count")
        row5.addWidget(label5C)
        self.kxCounttext=QtGui.QLineEdit("51")
        row5.addWidget(self.kxCounttext)

        label5D=QtGui.QLabel("ky range")
        label5D.setFont(bFont)
        row5.addWidget(label5D)

        label5E=QtGui.QLabel("Min and Max (Ang^-1)")
        row5.addWidget(label5E)
        self.kyMintext=QtGui.QLineEdit("-0.5")
        row5.addWidget(self.kyMintext)
        self.kyMaxtext=QtGui.QLineEdit("0.5")
        row5.addWidget(self.kyMaxtext)
        label5F=QtGui.QLabel("Count")
        row5.addWidget(label5F)
        self.kyCounttext=QtGui.QLineEdit("51")
        row5.addWidget(self.kyCounttext) 

        
        label5G=QtGui.QLabel("E range")
        label5G.setFont(bFont)
        row5.addWidget(label5G)

        label5H=QtGui.QLabel("Min and Max from EF (eV)")
        row5.addWidget(label5H)
        self.eMintext=QtGui.QLineEdit("-2")
        row5.addWidget(self.eMintext)
        self.eMaxtext=QtGui.QLineEdit("0.5")
        row5.addWidget(self.eMaxtext)
        label5I=QtGui.QLabel("Count")
        row5.addWidget(label5I)
        self.eCounttext=QtGui.QLineEdit("51")
        row5.addWidget(self.eCounttext) 

        # Row 6: boradening
        row6=QtGui.QHBoxLayout()
        row6.setAlignment(QtCore.Qt.AlignLeft)
        vbox.addLayout(row6)

        label6A=QtGui.QLabel("k broadening (Ang^-1)")
        label6A.setFont(bFont)
        row6.addWidget(label6A)
        self.sigmaktext=QtGui.QLineEdit("0.01")
        row6.addWidget(self.sigmaktext)

        label6B=QtGui.QLabel("E broadening (eV)")
        label6B.setFont(bFont)
        row6.addWidget(label6B)
        self.sigmaetext=QtGui.QLineEdit("0.05")
        row6.addWidget(self.sigmaetext)

        # Row 7: start calculation and plot
        row7=QtGui.QHBoxLayout()
        row7.setAlignment(QtCore.Qt.AlignLeft)
        vbox.addLayout(row7)
        self.startCalc=QtGui.QPushButton("Start calculation")
        row7.addWidget(self.startCalc)

        self.plotChoice=QtGui.QButtonGroup()
        self.plotDisp1=QtGui.QRadioButton("Original")
        self.plotDisp1.setFont(bFont)
        self.plotChoice.addButton(self.plotDisp1)
        row7.addWidget(self.plotDisp1)
        self.plotDisp2=QtGui.QRadioButton("Refracted")
        self.plotDisp2.setFont(bFont)
        self.plotChoice.addButton(self.plotDisp2)
        row7.addWidget(self.plotDisp2)

        self.export=QtGui.QPushButton("Export")
        row7.addWidget(self.export)

        self.importH5=QtGui.QPushButton("Import")
        row7.addWidget(self.importH5)

        # Row 8: kx, ky, and e index
        row8=QtGui.QHBoxLayout()
        row8.setAlignment(QtCore.Qt.AlignLeft)
        vbox.addLayout(row8)

        label8A=QtGui.QLabel("kx index")
        row8.addWidget(label8A)
        self.kxIndex=QtGui.QSpinBox()
        self.kxIndex.setSingleStep(1)
        self.kxIndex.setMinimum(0)
        row8.addWidget(self.kxIndex)

        self.kxValue=QtGui.QLabel()
        row8.addWidget(self.kxValue)

        label8B=QtGui.QLabel("ky index")
        row8.addWidget(label8B)
        self.kyIndex=QtGui.QSpinBox()
        self.kyIndex.setSingleStep(1)
        self.kyIndex.setMinimum(0)
        row8.addWidget(self.kyIndex)

        self.kyValue=QtGui.QLabel()
        row8.addWidget(self.kyValue)

        label8C=QtGui.QLabel("Energy index")
        row8.addWidget(label8C)
        self.eIndex=QtGui.QSpinBox()
        self.eIndex.setSingleStep(1)
        self.eIndex.setMinimum(0)
        row8.addWidget(self.eIndex)
        self.eValue=QtGui.QLabel()
        row8.addWidget(self.eValue)

        # Row 9: dispersion
        row9=QtGui.QHBoxLayout()
        row9.setAlignment(QtCore.Qt.AlignLeft)
        vbox.addLayout(row9)

        self.plot3=pg.GraphicsLayoutWidget()
        self.plot3.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        row9.addWidget(self.plot3)

        cmap=pg.colormap.get("CET-L9")
        
        labelStyle={"font-size":str(Config.fontSize_normal)+"px", "color": "white"}

        self.plotEy=self.plot3.addPlot()
        self.plotxy=self.plot3.addPlot()
        self.plot3.nextRow()
        self.plot3.nextColumn()
        self.plotEx=self.plot3.addPlot()

        self.imgEx=pg.ImageItem()
        self.imgEy=pg.ImageItem()
        self.imgxy=pg.ImageItem()
        self.plotEx.addItem(self.imgEx)
        self.plotEy.addItem(self.imgEy)
        self.plotxy.addItem(self.imgxy)
        self.barEx=pg.ColorBarItem(colorMap=cmap, values=(0,1))
        self.barEy=pg.ColorBarItem(colorMap=cmap, values=(0,1))
        self.barxy=pg.ColorBarItem(colorMap=cmap, values=(0,1))
        self.barEx.setImageItem(self.imgEx)
        self.barEy.setImageItem(self.imgEy)
        self.barxy.setImageItem(self.imgxy)
        
        self.plotEx.getAxis("bottom").setStyle(tickFont=font,tickLength=Config.tickLength)
        self.plotEx.getAxis("left").setStyle(tickFont=font,tickLength=Config.tickLength)
        self.plotEx.getAxis("bottom").setPen((255,255,255))
        self.plotEx.getAxis("left").setPen((255,255,255))
        self.plotEx.getAxis("bottom").setTextPen((255,255,255))
        self.plotEx.getAxis("left").setTextPen((255,255,255))
        self.plotEx.getAxis("bottom").setLabel(**labelStyle)
        self.plotEx.getAxis("left").setLabel(**labelStyle)
        self.plotEx.showGrid(x=True, y=True, alpha=1.0)
        self.plotEx.getAxis("bottom").setZValue(1)
        self.plotEx.getAxis("left").setZValue(1)

        
        self.plotEy.getAxis("bottom").setStyle(tickFont=font,tickLength=Config.tickLength)
        self.plotEy.getAxis("left").setStyle(tickFont=font,tickLength=Config.tickLength)
        self.plotEy.getAxis("bottom").setPen((255,255,255))
        self.plotEy.getAxis("left").setPen((255,255,255))
        self.plotEy.getAxis("bottom").setTextPen((255,255,255))
        self.plotEy.getAxis("left").setTextPen((255,255,255))
        self.plotEy.getAxis("bottom").setLabel(**labelStyle)
        self.plotEy.getAxis("left").setLabel(**labelStyle)
        self.plotEy.showGrid(x=True, y=True, alpha=1.0)
        self.plotEy.getAxis("bottom").setZValue(1)
        self.plotEy.getAxis("left").setZValue(1)

        
        self.plotxy.getAxis("bottom").setStyle(tickFont=font,tickLength=Config.tickLength)
        self.plotxy.getAxis("left").setStyle(tickFont=font,tickLength=Config.tickLength)
        self.plotxy.getAxis("bottom").setPen((255,255,255))
        self.plotxy.getAxis("left").setPen((255,255,255))
        self.plotxy.getAxis("bottom").setTextPen((255,255,255))
        self.plotxy.getAxis("left").setTextPen((255,255,255))
        self.plotxy.getAxis("bottom").setLabel(**labelStyle)
        self.plotxy.getAxis("left").setLabel(**labelStyle)
        self.plotxy.showGrid(x=True, y=True, alpha=1.0)
        self.plotxy.getAxis("bottom").setZValue(1)
        self.plotxy.getAxis("left").setZValue(1)
        
        self.plotEx.setLabel(axis="left", text="E-EF (eV)")
        self.plotEx.setLabel(axis="bottom", text="kx (ang^-1)")
        
        self.plotEy.setLabel(axis="bottom", text="E-EF (eV)")
        self.plotEy.setLabel(axis="left", text="ky (ang^-1)")

        self.plotxy.setLabel(axis="left", text="ky (ang^-1)")
        self.plotxy.setLabel(axis="bottom", text="kx (ang^-1)")

        
        self.vLineEx=pg.InfiniteLine(angle=90, movable=False, pen=Config.pen1)
        self.plotEx.addItem(self.vLineEx, ignoreBounds=True)
        self.hLineEx=pg.InfiniteLine(angle=0, movable=False, pen=Config.pen1)
        self.plotEx.addItem(self.hLineEx, ignoreBounds=True)

        self.vLineEy=pg.InfiniteLine(angle=90, movable=False, pen=Config.pen1)
        self.plotEy.addItem(self.vLineEy, ignoreBounds=True)
        self.hLineEy=pg.InfiniteLine(angle=0, movable=False, pen=Config.pen1)
        self.plotEy.addItem(self.hLineEy, ignoreBounds=True)

        self.vLinexy=pg.InfiniteLine(angle=90, movable=False, pen=Config.pen1)
        self.plotxy.addItem(self.vLinexy, ignoreBounds=True)
        self.hLinexy=pg.InfiniteLine(angle=0, movable=False, pen=Config.pen1)
        self.plotxy.addItem(self.hLinexy, ignoreBounds=True)

        self.bLineEx=pg.InfiniteLine(angle=0, movable=False, pen=Config.pen2)
        self.plotEx.addItem(self.bLineEx, ignoreBounds=True)
        self.bLineEy=pg.InfiniteLine(angle=90, movable=False, pen=Config.pen2)
        self.plotEy.addItem(self.bLineEy, ignoreBounds=True)

        def changeKXYIndices(e):
            if e.key()==QtCore.Qt.Key_Down:
                self.kyIndex.setValue(self.kyIndex.value()-1)
            elif e.key()==QtCore.Qt.Key_Up:
                self.kyIndex.setValue(self.kyIndex.value()+1)
            elif e.key()==QtCore.Qt.Key_Right:
                self.kxIndex.setValue(self.kxIndex.value()+1)
            elif e.key()==QtCore.Qt.Key_Left:
                self.kxIndex.setValue(self.kxIndex.value()-1)
            elif e.key()==QtCore.Qt.Key_PageUp:
                self.eIndex.setValue(self.eIndex.value()+1)
            elif e.key()==QtCore.Qt.Key_PageDown:
                self.eIndex.setValue(self.eIndex.value()-1)
            else:
                return

        
        self.plot3.keyPressEvent=changeKXYIndices


        
app=QtGui.QApplication([])
win=MainWindow()
font=QtGui.QFont()
font.setPixelSize(Config.fontSize_normal)
font.setFamilies(Config.fontFamilies)
win.setFont(font)

dispCube1=None
dispCube2=None
def startCalc():
    print("----")
    print("Check input parameters...")
    try:
        W_eV=float(win.Wtext.text())
        W=W_eV/Config.Eh_eV
        print(("{0:32s} = {1:.2f} eV = {2:.2f} Eh").format("Work function", W_eV, W))
        V0_eV=float(win.V0text.text())
        V0=V0_eV/Config.Eh_eV
        print(("{0:32s} = {1:.2f} eV = {2:.2f} Eh").format("Inner potential", V0_eV, V0))

        a=float(win.atext.text())
        V1_eV=float(win.V1text.text())
        V1=V1_eV/Config.Eh_eV

        k0_ang=np.zeros((3))
        k0_ang[0]=float(win.k0xtext.text())
        k0_ang[1]=float(win.k0ytext.text())
        k0_ang[2]=float(win.k0ztext.text())
        k0=k0_ang*Config.Bohr_ang
        print(("{0:32s} = (k - ({1:.2f}, {2:.2f}, {3:.2f}))^2 * {4:.2f} /2 - {5:.2f} + {6:.2f} (Eh)").format("Initial state dispersion", k0[0], k0[1], k0[2], a, W, V1))

        kFlat=True
        if win.kFlat.isChecked()==True:
            pass
        elif win.kCurved.isChecked()==True:
            kFlat=False
        else:
            print("Error: kFlat or kCurved should be checked")
            return

        kFlat_kz=0
        kCurved_k=0
        if kFlat==True:
            kFlat_kz_ang=float(win.kFlat_kz.text())
            kFlat_kz=kFlat_kz_ang*Config.Bohr_ang
            print(("{0:32s} = flat, kz = {1:.2f} Ang^-1 = {2:.2f} Bohr^-1").format("k plane", kFlat_kz_ang, kFlat_kz))
        else:
            kCurved_k_ang=float(win.kCurved_k.text())
            kCurved_k=kCurved_k_ang*Config.Bohr_ang
            print(("{0:32s} = curved, |k| = {1:.2f} Ang^-1 = {2:.2f} Bohr^-1").format("k plane", kCurved_k_ang, kCurved_k))

        surfaceConst=True
        if win.surfaceConst.isChecked()==True:
            pass
        elif win.surfaceRandom.isChecked()==True:
            surfaceConst=False
        else:
            print("Error: surfaceConst or surfaceRandom should be checked")
            return

        surfaceConst_theta=0
        surfaceConst_phi=0
        surfaceRandom_samples=0

        if surfaceConst==True:
            surfaceConst_theta_deg=float(win.surfaceConst_theta.text())
            surfaceConst_theta=math.radians(surfaceConst_theta_deg)
            surfaceConst_phi_deg=float(win.surfaceConst_phi.text())
            surfaceConst_phi=math.radians(surfaceConst_phi_deg)
            print(("{0:32s} = constant, theta = {1:.1f} deg = {2:.2f} rad, phi = {3:.1f} deg = {4:.2f} rad").format("Surface orientation", surfaceConst_theta_deg, surfaceConst_theta, surfaceConst_phi_deg, surfaceConst_phi))
        else:
            surfaceRandom_samples=int(win.surfaceRandom_samples.text())
            print(("{0:32s} = random, samples = {1:d}").format("Surface orientation", surfaceRandom_samples))

        kxMin_ang=float(win.kxMintext.text())
        kxMin=kxMin_ang*Config.Bohr_ang
        kxMax_ang=float(win.kxMaxtext.text())
        kxMax=kxMax_ang*Config.Bohr_ang
        kxCount=int(win.kxCounttext.text())

        print(("{0:32s} = {1:.2f} to {2:.2f}, {3:d} points").format("kx range (Ang^-1)", kxMin_ang, kxMax_ang, kxCount))
        print(("{0:32s} = {1:.2f} to {2:.2f}, {3:d} points").format("kx range (Bohr^-1)", kxMin, kxMax, kxCount))

        kyMin_ang=float(win.kyMintext.text())
        kyMin=kyMin_ang*Config.Bohr_ang
        kyMax_ang=float(win.kyMaxtext.text())
        kyMax=kyMax_ang*Config.Bohr_ang
        kyCount=int(win.kyCounttext.text())

        print(("{0:32s} = {1:.2f} to {2:.2f}, {3:d} points").format("ky range (Ang^-1)", kyMin_ang, kyMax_ang, kyCount))
        print(("{0:32s} = {1:.2f} to {2:.2f}, {3:d} points").format("ky range (Bohr^-1)", kyMin, kyMax, kyCount))

        eMin_eV=float(win.eMintext.text())
        eMin=eMin_eV/Config.Eh_eV-W
        eMax_eV=float(win.eMaxtext.text())
        eMax=eMax_eV/Config.Eh_eV-W
        eCount=int(win.eCounttext.text())

        print(("{0:32s} = {1:.2f} to {2:.2f}, {3:d} points").format("E range from EF (eV)", eMin_eV, eMax_eV, eCount))
        print(("{0:32s} = {1:.2f} to {2:.2f}, {3:d} points").format("E range from Vacuum (eV)", eMin*Config.Eh_eV, eMax*Config.Eh_eV, eCount))
        print(("{0:32s} = {1:.2f} to {2:.2f}, {3:d} points").format("E range from Vacuum (Eh)", eMin, eMax, eCount))

        sigmak_ang=float(win.sigmaktext.text())
        sigmak=sigmak_ang*Config.Bohr_ang
        sigmae_eV=float(win.sigmaetext.text())
        sigmae=sigmae_eV/Config.Eh_eV

        print(("{0:32s} = {1:.2f} Ang^-1 = {2:.2f} Bohr^-1").format("k broadening", sigmak_ang, sigmak))
        print(("{0:32s} = {1:.2f} eV = {2:.3f} Eh").format("E broadening", sigmae_eV, sigmae))

        dkx=(kxMax-kxMin)/(kxCount-1)
        dky=(kyMax-kyMin)/(kyCount-1)
        de=(eMax-eMin)/(eCount-1)
        print(("{0:32s} = {1:.3f} Ang^-1, {2:.3f} Ang^-1, {3:.3f} eV").format("dkx, dky, de", dkx/Config.Bohr_ang, dky/Config.Bohr_ang, de*Config.Eh_eV))
        print(("{0:32s} = {1:.3f} Bohr^-1, {2:.3f} Bohr^-1, {3:.3f} Eh").format("dkx, dky, de", dkx, dky, de))       

    except Exception as e:
        print(e)
        return

    global dispCube1
    dispCube1=np.zeros((kxCount, kyCount, eCount))
    lib.calc1(W, V0, k0, a, V1, kFlat, kFlat_kz, kCurved_k, kxMin, kxMax, kxCount, dkx, kyMin, kyMax, kyCount, dky, eMin, eMax, eCount, de, sigmak, sigmae, dispCube1)

    global dispCube2
    dispCube2=np.zeros(dispCube1.shape)
    lib.calc2(W, V0, k0, a, V1, kFlat, kFlat_kz, kCurved_k, surfaceConst, surfaceConst_theta, surfaceConst_phi, surfaceRandom_samples, kxMin, kxMax, kxCount, dkx, kyMin, kyMax, kyCount, dky, eMin, eMax, eCount, de, sigmak, sigmae, dispCube2)

    win.kxIndex.setMaximum(kxCount-1)
    win.kyIndex.setMaximum(kyCount-1)
    win.eIndex.setMaximum(eCount-1)
    print("Calculation finished")

    plotDisp()

def plotDisp():
    global dispCube1
    global dispCube2

    if dispCube1 is None or dispCube2 is None:
        return

    kxIndex=win.kxIndex.value()
    kyIndex=win.kyIndex.value()
    eIndex=win.eIndex.value()

    try:    
        kxMin=float(win.kxMintext.text())
        kxMax=float(win.kxMaxtext.text())
        kxCount=int(win.kxCounttext.text())

        kyMin=float(win.kyMintext.text())
        kyMax=float(win.kyMaxtext.text())
        kyCount=int(win.kyCounttext.text())

        eMin=float(win.eMintext.text())
        eMax=float(win.eMaxtext.text())
        eCount=int(win.eCounttext.text())

        dkx=(kxMax-kxMin)/(kxCount-1)
        dky=(kyMax-kyMin)/(kyCount-1)
        de=(eMax-eMin)/(eCount-1)

    except Exception as e:
        print(e)
        return

    if dispCube1.shape[0]!=kxCount or dispCube1.shape[1]!=kyCount or dispCube1.shape[2]!=eCount or dispCube2.shape[0]!=kxCount or dispCube2.shape[1]!=kyCount or dispCube2.shape[2]!=eCount:
        print("Plot error: size mismatch")
        return 

    kxValue=kxMin+dkx*kxIndex
    kyValue=kyMin+dky*kyIndex
    eValue=eMin+de*eIndex

    win.kxValue.setText(("({0:.3f})").format(kxValue))
    win.kyValue.setText(("({0:.3f})").format(kyValue))
    win.eValue.setText(("({0:.3f})").format(eValue))

    win.vLineEx.setPos(kxValue)
    win.hLineEx.setPos(eValue)
    win.vLineEy.setPos(eValue)
    win.hLineEy.setPos(kyValue)
    win.vLinexy.setPos(kxValue)
    win.hLinexy.setPos(kyValue)

    if win.plotDisp1.isChecked()==True:
        Ex=dispCube1[:,kyIndex,:]
        Ey=dispCube1[kxIndex,:,:]
        xy=dispCube1[:,:,eIndex]
    elif win.plotDisp2.isChecked()==True:
        Ex=dispCube2[:,kyIndex,:]
        Ey=dispCube2[kxIndex,:,:]
        xy=dispCube2[:,:,eIndex]
    else:
        return


    tr_Ex=QtGui.QTransform()
    tr_Ex.translate(kxMin-dkx/2,eMin-de/2)
    tr_Ex.scale(dkx, de)
    win.imgEx.setTransform(tr_Ex)

    tr_Ey=QtGui.QTransform()
    tr_Ey.translate(eMin-de/2, kyMin-dky/2)
    tr_Ey.rotate(-90)
    tr_Ey.scale(-dky, de)
    win.imgEy.setTransform(tr_Ey)

    tr_xy=QtGui.QTransform()
    tr_xy.translate(kxMin-dkx/2, kyMin-dky/2)
    tr_xy.scale(dky, dky)
    win.imgxy.setTransform(tr_xy)

    win.imgEx.setImage(Ex)
    win.imgEy.setImage(Ey)
    win.imgxy.setImage(xy)

def importDisp():
    global dispCube1
    global dispCube2

    selectedFile, _filter=QtGui.QFileDialog.getOpenFileName(caption="Open file")
    if selectedFile!="":
        with h5py.File(selectedFile, "r") as f:
            dispCube1=np.array(f["Original"])
            dispCube2=np.array(f["Refracted"])
            offset=np.array(f.attrs["Offset"])
            delta=np.array(f.attrs["Delta"])
            size=np.array(f.attrs["Size"])

            win.kxMintext.setText(("{0:f}").format(offset[0]))
            win.kyMintext.setText(("{0:f}").format(offset[1]))
            win.eMintext.setText(("{0:f}").format(offset[2]))

            win.kxMaxtext.setText(("{0:f}").format(offset[0]+delta[0]*size[0]))
            win.kyMaxtext.setText(("{0:f}").format(offset[1]+delta[1]*size[1]))
            win.eMaxtext.setText(("{0:f}").format(offset[2]+delta[2]*size[2]))

            win.kxCounttext.setText(("{0:d}").format(size[0]))
            win.kyCounttext.setText(("{0:d}").format(size[1]))
            win.eCounttext.setText(("{0:d}").format(size[2]))

            win.Wtext.setText(("{0:f}").format(f.attrs["W"]))
            win.V0text.setText(("{0:f}").format(f.attrs["V0"]))
            win.V1text.setText(("{0:f}").format(f.attrs["V1"]))
            win.atext.setText(("{0:f}").format(f.attrs["a"]))

            k0=np.array(f.attrs["k0"])
            win.k0xtext.setText(("{0:f}").format(k0[0]))
            win.k0ytext.setText(("{0:f}").format(k0[1]))
            win.k0ztext.setText(("{0:f}").format(k0[2]))
            
            win.sigmaktext.setText(("{0:f}").format(f.attrs["sigmak"]))
            win.sigmaetext.setText(("{0:f}").format(f.attrs["sigmae"]))

            kPlane=f.attrs["kPlane"]
            if kPlane=="Flat":
                win.kFlat.setChecked(True)
                win.kCurved.setChecked(False)
                win.kFlat_kz.setText(("{0:f}").format(f.attrs["kPlane_kz"]))
            elif kPlane=="Curved":
                win.kFlat.setChecked(False)
                win.kCurved.setChecked(True)
                win.kCurved_k.setText(("{0:f}").format(f.attrs["kPlane_k"]))
            else:
                print(("Error: invalid kPlane {0:s}").format(kPlane))

            surfaceConst=f.attrs["Surface"]
            if surfaceConst=="Constant":
                win.surfaceConst.setChecked(True)
                win.surfaceRandom.setChecked(False)
                win.surfaceConst_theta.setText(("{0:f}").format(f.attrs["Surface_theta"]))
                win.surfaceConst_phi.setText(("{0:f}").format(f.attrs["Surface_phi"]))
            elif surfaceConst=="Random":
                win.surfaceConst.setChecked(False)
                win.surfaceRandom.setChecked(True)
                win.surfaceRandom_samples.setText(("{0:d}").format(f.attrs["Surface_samples"]))


    print("Import finished")
    plotDisp()
                               

    
def exportDisp():
    global dispCube1
    global dispCube2
    
    try:    
        kxMin=float(win.kxMintext.text())
        kxMax=float(win.kxMaxtext.text())
        kxCount=int(win.kxCounttext.text())

        kyMin=float(win.kyMintext.text())
        kyMax=float(win.kyMaxtext.text())
        kyCount=int(win.kyCounttext.text())

        eMin=float(win.eMintext.text())
        eMax=float(win.eMaxtext.text())
        eCount=int(win.eCounttext.text())

        dkx=(kxMax-kxMin)/(kxCount-1)
        dky=(kyMax-kyMin)/(kyCount-1)
        de=(eMax-eMin)/(eCount-1)

        
        W=float(win.Wtext.text())
        V0=float(win.V0text.text())

        a=float(win.atext.text())
        V1=float(win.V1text.text())

        k0=np.zeros((3))
        k0[0]=float(win.k0xtext.text())
        k0[1]=float(win.k0ytext.text())
        k0[2]=float(win.k0ztext.text())

        kFlat=True
        if win.kFlat.isChecked()==True:
            pass
        elif win.kCurved.isChecked()==True:
            kFlat=False
        else:
            return

        kFlat_kz=0
        kCurved_k=0
        if kFlat==True:
            kFlat_kz=float(win.kFlat_kz.text())
        else:
            kCurved_k=float(win.kCurved_k.text())

        surfaceConst=True
        if win.surfaceConst.isChecked()==True:
            pass
        elif win.surfaceRandom.isChecked()==True:
            surfaceConst=False
        else:
            return

        surfaceConst_theta=0
        surfaceConst_phi=0
        surfaceRandom_samples=0

        if surfaceConst==True:
            surfaceConst_theta=float(win.surfaceConst_theta.text())
            surfaceConst_phi=float(win.surfaceConst_phi.text())
        else:
            surfaceRandom_samples=int(win.surfaceRandom_samples.text())

            
        sigmak=float(win.sigmaktext.text())
        sigmae=float(win.sigmaetext.text())


    except Exception as e:
        print(e)
        return


    selectedFile, _filter=QtGui.QFileDialog.getSaveFileName(caption="Open file")
    if selectedFile!="":
        with h5py.File(selectedFile, "w") as f:
            f.attrs.create("Datetime", datetime.now().isoformat(" "))
            f.create_dataset("Original", data=dispCube1)
            f.create_dataset("Refracted", data=dispCube2)
            f.attrs.create("Offset", [kxMin, kyMin, eMin])
            f.attrs.create("Delta", [dkx, dky, de])
            f.attrs.create("Size", [kxCount, kyCount, eCount])
            f.attrs.create("W", W)
            f.attrs.create("V0", V0)
            f.attrs.create("V1", V1)
            f.attrs.create("a", a)
            f.attrs.create("k0", k0)
            f.attrs.create("sigmak", sigmak)
            f.attrs.create("sigmae", sigmae)

            if kFlat==True:
                f.attrs.create("kPlane", "Flat")
                f.attrs.create("kPlane_kz", kFlat_kz)
            else:
                f.attrs.create("kPlane", "Curved")
                f.attrs.create("kPlane_k", kCurved_k)

            if surfaceConst==True:
                f.attrs.create("Surface", "Constant")
                f.attrs.create("Surface_theta", surfaceConst_theta)
                f.attrs.create("Surface_phi", surfaceConst_phi)
            else:
                f.attrs.create("Surface", "Random")
                f.attrs.create("Surface_samples", surfaceRandom_samples)


    

win.startCalc.clicked.connect(startCalc)
win.plotDisp1.clicked.connect(plotDisp)
win.plotDisp2.clicked.connect(plotDisp)
win.export.clicked.connect(exportDisp)
win.importH5.clicked.connect(importDisp)
win.kxIndex.valueChanged.connect(plotDisp)
win.kyIndex.valueChanged.connect(plotDisp)
win.eIndex.valueChanged.connect(plotDisp)

pg.setConfigOptions(antialias=True)
win.show()
app.exec_()
