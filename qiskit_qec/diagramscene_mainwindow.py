from qiskit_qec.diagramscene import DiagramTextItem, DiagramScene

import math
import sys

from PySide6.QtCore import (QLineF, QPointF, QRect, QRectF, QSize, QSizeF, Qt,
                            Signal)
from PySide6.QtGui import (QAction, QColor, QFont, QIcon, QIntValidator,
                           QPainter, QPainterPath, QPen, QPixmap, QPolygonF,
                           QBrush)
from PySide6.QtWidgets import (
    QApplication, QButtonGroup, QComboBox, QFontComboBox, QGraphicsAnchorLayout,
    QGraphicsItem, QGraphicsLineItem, QGraphicsPolygonItem, QGraphicsTextItem,
    QGraphicsScene, QGraphicsView, QGridLayout, QHBoxLayout, QLabel,
    QMainWindow, QMenu, QMessageBox, QSizePolicy, QToolBox, QToolButton,
    QWidget)

import qiskit_qec.diagramscene_rc


class MainWindow(QMainWindow):
    insert_text_button = 10

    def __init__(self):
        super().__init__()

        self.create_actions()
        self.create_menus()
        self.create_tool_box()

        self.scene = DiagramScene(self._item_menu)
        self.scene.setSceneRect(QRectF(0, 0, 5000, 5000))
#        self.scene.item_inserted.connect(self.item_inserted)
        self.scene.text_inserted.connect(self.text_inserted)
        self.scene.item_selected.connect(self.item_selected)

        self.create_toolbars()

        layout = QHBoxLayout()
        layout.addWidget(self._tool_box)
        self.view = QGraphicsView(self.scene)
        layout.addWidget(self.view)

        self.widget = QWidget()
        self.widget.setLayout(layout)

        self.setCentralWidget(self.widget)
        self.setWindowTitle("Diagramscene")

    def background_button_group_clicked(self, button):
        buttons = self._background_button_group.buttons()
        for myButton in buttons:
            if myButton != button:
                button.setChecked(False)

        text = button.text()
        if text == "Blue Grid":
            self.scene.setBackgroundBrush(
                QBrush(QPixmap(':/images/background1.png')))
        elif text == "White Grid":
            self.scene.setBackgroundBrush(
                QBrush(QPixmap(':/images/background2.png')))
        elif text == "Gray Grid":
            self.scene.setBackgroundBrush(
                QBrush(QPixmap(':/images/background3.png')))
        else:
            self.scene.setBackgroundBrush(
                QBrush(QPixmap(':/images/background4.png')))

        self.scene.update()
        self.view.update()

    def button_group_clicked(self, idx):
        buttons = self._button_group.buttons()
        for button in buttons:
            if self._button_group.button(idx) != button:
                button.setChecked(False)

        if idx == self.insert_text_button:
            self.scene.set_mode(DiagramScene.InsertText)
        else:
            self.scene.set_mode(DiagramScene.InsertItem)

    def pointer_group_clicked(self, i):
        self.scene.set_mode(self._pointer_type_group.checkedId())

    def bring_to_front(self):
        if not self.scene.selectedItems():
            return

        selected_item = self.scene.selectedItems()[0]
        overlap_items = selected_item.collidingItems()

        z_value = 0
        # for item in overlap_items:
        #     if (item.zValue() >= z_value and isinstance(item, Qubit)):
        #         z_value = item.zValue() + 0.1
        selected_item.setZValue(z_value)

    def send_to_back(self):
        if not self.scene.selectedItems():
            return

        selected_item = self.scene.selectedItems()[0]
        overlap_items = selected_item.collidingItems()

        z_value = 0
        # for item in overlap_items:
        #    if (item.zValue() <= z_value and isinstance(item, Qubit)):
        #       z_value = item.zValue() - 0.1
        selected_item.setZValue(z_value)

    def item_inserted(self, item):
        self._pointer_type_group.button(DiagramScene.MoveItem).setChecked(True)
        #self.scene.set_mode(self._pointer_type_group.checkedId())
        #xself._button_group.button(item).setChecked(False)

    def text_inserted(self, item):
        self._button_group.button(self.insert_text_button).setChecked(False)
        self.scene.set_mode(self._pointer_type_group.checkedId())

    def current_font_changed(self, font):
        self.handle_font_change()

    def font_size_changed(self, font):
        self.handle_font_change()

    def scene_scale_changed(self, scale):
        new_scale = int(scale[:-1]) / 100.0
        old_matrix = self.view.transform()
        self.view.resetTransform()
        self.view.translate(old_matrix.dx(), old_matrix.dy())
        self.view.scale(new_scale, new_scale)

    def text_color_changed(self):
        self._text_action = self.sender()
        self._font_color_tool_button.setIcon(
            self.create_color_tool_button_icon(':/images/textpointer.png',
                                               QColor(
                                                   self._text_action.data())))
        self.text_button_triggered()

    def item_color_changed(self):
        self._fill_action = self.sender()
        self._fill_color_tool_button.setIcon(
            self.create_color_tool_button_icon(':/images/floodfill.png',
                                               QColor(
                                                   self._fill_action.data())))
        self.fill_button_triggered()

    def line_color_changed(self):
        self._line_action = self.sender()
        self._line_color_tool_button.setIcon(
            self.create_color_tool_button_icon(':/images/linecolor.png',
                                               QColor(
                                                   self._line_action.data())))
        self.line_button_triggered()

    def text_button_triggered(self):
        self.scene.set_text_color(QColor(self._text_action.data()))

    def fill_button_triggered(self):
        self.scene.set_item_color(QColor(self._fill_action.data()))

    def line_button_triggered(self):
        self.scene.set_line_color(QColor(self._line_action.data()))

    def handle_font_change(self):
        font = self._font_combo.currentFont()
        font.setPointSize(int(self._font_size_combo.currentText()))
        if self._bold_action.isChecked():
            font.setWeight(QFont.Bold)
        else:
            font.setWeight(QFont.Normal)
        font.setItalic(self._italic_action.isChecked())
        font.setUnderline(self._underline_action.isChecked())

        self.scene.set_font(font)

    def item_selected(self, item):
        font = item.font()
        color = item.defaultTextColor()
        self._font_combo.setCurrentFont(font)
        self._font_size_combo.setEditText(str(font.pointSize()))
        self._bold_action.setChecked(font.weight() == QFont.Bold)
        self._italic_action.setChecked(font.italic())
        self._underline_action.setChecked(font.underline())

    def about(self):
        QMessageBox.about(
            self, "About Diagram Scene",
            "The <b>Diagram Scene</b> example shows use of the graphics framework."
        )

    def create_tool_box(self):
        self._button_group = QButtonGroup()
        self._button_group.setExclusive(False)
        self._button_group.idClicked.connect(self.button_group_clicked)

        layout = QGridLayout()
        #layout.addWidget(self.create_cell_widget(Qubit.name, Qubit), 0, 0)

        text_button = QToolButton()
        text_button.setCheckable(True)
        self._button_group.addButton(text_button, self.insert_text_button)
        text_button.setIcon(
            QIcon(QPixmap(':/images/textpointer.png').scaled(30, 30)))
        text_button.setIconSize(QSize(50, 50))

        text_layout = QGridLayout()
        text_layout.addWidget(text_button, 0, 0, Qt.AlignHCenter)
        text_layout.addWidget(QLabel("Text"), 1, 0, Qt.AlignCenter)
        text_widget = QWidget()
        text_widget.setLayout(text_layout)
        layout.addWidget(text_widget, 1, 1)

        layout.setRowStretch(3, 10)
        layout.setColumnStretch(2, 10)

        item_widget = QWidget()
        item_widget.setLayout(layout)

        self._background_button_group = QButtonGroup()
        self._background_button_group.buttonClicked.connect(
            self.background_button_group_clicked)

        background_layout = QGridLayout()
        background_layout.addWidget(
            self.create_background_cell_widget("Blue Grid",
                                               ':/images/background1.png'), 0,
            0)
        background_layout.addWidget(
            self.create_background_cell_widget("White Grid",
                                               ':/images/background2.png'), 0,
            1)
        background_layout.addWidget(
            self.create_background_cell_widget("Gray Grid",
                                               ':/images/background3.png'), 1,
            0)
        background_layout.addWidget(
            self.create_background_cell_widget("No Grid",
                                               ':/images/background4.png'), 1,
            1)

        background_layout.setRowStretch(2, 10)
        background_layout.setColumnStretch(2, 10)

        background_widget = QWidget()
        background_widget.setLayout(background_layout)

        self._tool_box = QToolBox()
        self._tool_box.setSizePolicy(
            QSizePolicy(QSizePolicy.Maximum, QSizePolicy.Ignored))
        self._tool_box.setMinimumWidth(item_widget.sizeHint().width())
        self._tool_box.addItem(item_widget, "Basic Flowchart Shapes")
        self._tool_box.addItem(background_widget, "Backgrounds")

    def create_actions(self):
        self._to_front_action = QAction(QIcon(':/images/bringtofront.png'),
                                        "Bring to &Front",
                                        self,
                                        shortcut="Ctrl+F",
                                        statusTip="Bring item to front",
                                        triggered=self.bring_to_front)

        self._send_back_action = QAction(QIcon(':/images/sendtoback.png'),
                                         "Send to &Back",
                                         self,
                                         shortcut="Ctrl+B",
                                         statusTip="Send item to back",
                                         triggered=self.send_to_back)

        self._exit_action = QAction("E&xit",
                                    self,
                                    shortcut="Ctrl+X",
                                    statusTip="Quit Scenediagram example",
                                    triggered=self.close)

        self._bold_action = QAction(QIcon(':/images/bold.png'),
                                    "Bold",
                                    self,
                                    checkable=True,
                                    shortcut="Ctrl+B",
                                    triggered=self.handle_font_change)

        self._italic_action = QAction(QIcon(':/images/italic.png'),
                                      "Italic",
                                      self,
                                      checkable=True,
                                      shortcut="Ctrl+I",
                                      triggered=self.handle_font_change)

        self._underline_action = QAction(QIcon(':/images/underline.png'),
                                         "Underline",
                                         self,
                                         checkable=True,
                                         shortcut="Ctrl+U",
                                         triggered=self.handle_font_change)

        self._about_action = QAction("A&bout",
                                     self,
                                     shortcut="Ctrl+B",
                                     triggered=self.about)

    def create_menus(self):
        self._file_menu = self.menuBar().addMenu("&File")
        self._file_menu.addAction(self._exit_action)

        self._item_menu = self.menuBar().addMenu("&Item")
        self._item_menu.addSeparator()
        self._item_menu.addAction(self._to_front_action)
        self._item_menu.addAction(self._send_back_action)

        self._about_menu = self.menuBar().addMenu("&Help")
        self._about_menu.addAction(self._about_action)

    def create_toolbars(self):
        self._edit_tool_bar = self.addToolBar("Edit")
        self._edit_tool_bar.addAction(self._to_front_action)
        self._edit_tool_bar.addAction(self._send_back_action)

        self._font_combo = QFontComboBox()
        self._font_combo.currentFontChanged.connect(self.current_font_changed)

        self._font_size_combo = QComboBox()
        self._font_size_combo.setEditable(True)
        for i in range(8, 30, 2):
            self._font_size_combo.addItem(str(i))
        validator = QIntValidator(2, 64, self)
        self._font_size_combo.setValidator(validator)
        self._font_size_combo.currentIndexChanged.connect(
            self.font_size_changed)

        self._font_color_tool_button = QToolButton()
        self._font_color_tool_button.setPopupMode(QToolButton.MenuButtonPopup)
        self._font_color_tool_button.setMenu(
            self.create_color_menu(self.text_color_changed, Qt.black))
        self._text_action = self._font_color_tool_button.menu().defaultAction()
        self._font_color_tool_button.setIcon(
            self.create_color_tool_button_icon(':/images/textpointer.png',
                                               Qt.black))
        self._font_color_tool_button.setAutoFillBackground(True)
        self._font_color_tool_button.clicked.connect(self.text_button_triggered)

        self._fill_color_tool_button = QToolButton()
        self._fill_color_tool_button.setPopupMode(QToolButton.MenuButtonPopup)
        self._fill_color_tool_button.setMenu(
            self.create_color_menu(self.item_color_changed, Qt.white))
        self._fill_action = self._fill_color_tool_button.menu().defaultAction()
        self._fill_color_tool_button.setIcon(
            self.create_color_tool_button_icon(':/images/floodfill.png',
                                               Qt.white))
        self._fill_color_tool_button.clicked.connect(self.fill_button_triggered)

        self._line_color_tool_button = QToolButton()
        self._line_color_tool_button.setPopupMode(QToolButton.MenuButtonPopup)
        self._line_color_tool_button.setMenu(
            self.create_color_menu(self.line_color_changed, Qt.black))
        self._line_action = self._line_color_tool_button.menu().defaultAction()
        self._line_color_tool_button.setIcon(
            self.create_color_tool_button_icon(':/images/linecolor.png',
                                               Qt.black))
        self._line_color_tool_button.clicked.connect(self.line_button_triggered)

        self._text_tool_bar = self.addToolBar("Font")
        self._text_tool_bar.addWidget(self._font_combo)
        self._text_tool_bar.addWidget(self._font_size_combo)
        self._text_tool_bar.addAction(self._bold_action)
        self._text_tool_bar.addAction(self._italic_action)
        self._text_tool_bar.addAction(self._underline_action)

        self._color_tool_bar = self.addToolBar("Color")
        self._color_tool_bar.addWidget(self._font_color_tool_button)
        self._color_tool_bar.addWidget(self._fill_color_tool_button)
        self._color_tool_bar.addWidget(self._line_color_tool_button)

        pointer_button = QToolButton()
        pointer_button.setCheckable(True)
        pointer_button.setChecked(True)
        pointer_button.setIcon(QIcon(':/images/pointer.png'))
        line_pointer_button = QToolButton()
        line_pointer_button.setCheckable(True)
        line_pointer_button.setIcon(QIcon(':/images/linepointer.png'))

        self._pointer_type_group = QButtonGroup()
        self._pointer_type_group.addButton(pointer_button,
                                           DiagramScene.MoveItem)
        self._pointer_type_group.addButton(line_pointer_button,
                                           DiagramScene.InsertLine)
        self._pointer_type_group.idClicked.connect(self.pointer_group_clicked)

        self._scene_scale_combo = QComboBox()
        self._scene_scale_combo.addItems(["50%", "75%", "100%", "125%", "150%"])
        self._scene_scale_combo.setCurrentIndex(2)
        self._scene_scale_combo.currentTextChanged.connect(
            self.scene_scale_changed)

        self._pointer_toolbar = self.addToolBar("Pointer type")
        self._pointer_toolbar.addWidget(pointer_button)
        self._pointer_toolbar.addWidget(line_pointer_button)
        self._pointer_toolbar.addWidget(self._scene_scale_combo)

    def create_background_cell_widget(self, text, image):
        button = QToolButton()
        button.setText(text)
        button.setIcon(QIcon(image))
        button.setIconSize(QSize(50, 50))
        button.setCheckable(True)
        self._background_button_group.addButton(button)

        layout = QGridLayout()
        layout.addWidget(button, 0, 0, Qt.AlignHCenter)
        layout.addWidget(QLabel(text), 1, 0, Qt.AlignCenter)

        widget = QWidget()
        widget.setLayout(layout)

        return widget

    def create_cell_widget(self, text, class_type):
        item = class_type(self._item_menu)
        icon = QIcon(item.image())

        button = QToolButton()
        button.setIcon(icon)
        button.setIconSize(QSize(50, 50))
        button.setCheckable(True)
        self._button_group.addButton(button)

        layout = QGridLayout()
        layout.addWidget(button, 0, 0, Qt.AlignHCenter)
        layout.addWidget(QLabel(text), 1, 0, Qt.AlignCenter)

        widget = QWidget()
        widget.setLayout(layout)

        return widget

    def create_color_menu(self, slot, defaultColor):
        colors = [Qt.black, Qt.white, Qt.red, Qt.blue, Qt.yellow]
        names = ["black", "white", "red", "blue", "yellow"]

        color_menu = QMenu(self)
        for color, name in zip(colors, names):
            action = QAction(self.create_color_icon(color),
                             name,
                             self,
                             triggered=slot)
            action.setData(QColor(color))
            color_menu.addAction(action)
            if color == defaultColor:
                color_menu.setDefaultAction(action)
        return color_menu

    def create_color_tool_button_icon(self, imageFile, color):
        pixmap = QPixmap(50, 80)
        pixmap.fill(Qt.transparent)
        painter = QPainter(pixmap)
        image = QPixmap(imageFile)
        target = QRect(0, 0, 50, 60)
        source = QRect(0, 0, 42, 42)
        painter.fillRect(QRect(0, 60, 50, 80), color)
        painter.drawPixmap(target, image, source)
        painter.end()

        return QIcon(pixmap)

    def create_color_icon(self, color):
        pixmap = QPixmap(20, 20)
        painter = QPainter(pixmap)
        painter.setPen(Qt.NoPen)
        painter.fillRect(QRect(0, 0, 20, 20), color)
        painter.end()

        return QIcon(pixmap)


if __name__ == '__main__':
    app = QApplication(sys.argv)

    main_window = MainWindow()
    main_window.setGeometry(100, 100, 800, 500)
    main_window.show()

    sys.exit(app.exec())
