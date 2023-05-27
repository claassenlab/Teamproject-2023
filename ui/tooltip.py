from tkinter import *


class ToolTip(object):

    def __init__(self, widget):
        self.widget = widget
        self.tip_window = None

    def show(self, text):
        self.text = text
        if self.tip_window or not self.text:
            return
        x, y, cx, cy = self.widget.bbox("insert")
        x = x + self.widget.winfo_rootx()
        y = y + cy + self.widget.winfo_rooty() - 20
        self.tip_window = tw = Toplevel(self.widget)
        tw.wm_overrideredirect(1)
        tw.wm_geometry("+%d+%d" % (x, y))
        label = Label(tw, text=self.text, justify=LEFT,
                      background="#ffffe0", relief=SOLID, borderwidth=1,
                      font=("tahoma", "8", "normal"))
        label.pack(ipadx=1)

    def hide(self):
        tw = self.tip_window
        self.tip_window = None
        if tw:
            tw.destroy()


def CreateToolTip(widget, text):
    tool_tip = ToolTip(widget)

    def enter(event):
        tool_tip.show(text)

    def leave(event):
        tool_tip.hide()

    widget.bind('<Enter>', enter)
    widget.bind('<Leave>', leave)
