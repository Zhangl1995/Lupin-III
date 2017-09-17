from tkinter import *
import PIL
from PIL import Image
from PIL import ImageTk


i=0

def showimage_prev():
    global i
    global var
    global label
    global R6
    if(i!=0):
        i=i-1
        img=Image.open(imagelist[i].strip())
#        canvas=Canvas(top,height=100,width=100)
        image=ImageTk.PhotoImage(img,top)
        imgs=canvas.create_image(50,50,image=image)
        canvas.pack()
        
#        var = StringVar()
#        label = Message( top, textvariable=var, relief=RAISED )
        var.set(imagelist[i][11:].strip())
        label.pack()
        R6.select()
        top.mainloop()
    else:
#        var = StringVar()
#        label = Message( top, textvariable=var, relief=RAISED )
        var.set('The Begining!')
        label.pack()
        top.mainloop()

def showimage_next():
    global i
    global var
    global label
    global R6
    if(i!=(len(imagelist)-1)):
        i=i+1
        img=Image.open(imagelist[i].strip())
#        canvas=Canvas(top,height=100,width=100)
        image=ImageTk.PhotoImage(img,top)
        imgs=canvas.create_image(50,50,image=image)
        canvas.pack()

#        var = StringVar()
#        label = Message( top, textvariable=var, relief=RAISED )
        var.set(imagelist[i][11:].strip())
        label.pack()
        R6.select()
        top.mainloop()
    else:
#        var = StringVar()
#        label = Message( top, textvariable=var, relief=RAISED )
        var.set('The Ending!')
        label.pack()
        top.mainloop()

def sel1():
    class_result=open('c:/users/a/classfication.txt','a+')
    class_result.write(imagelist[i][11:].strip() + '   Class 1\n')
    class_result.close()
def sel2():
    class_result=open('c:/users/a/classfication.txt','a+')
    class_result.write(imagelist[i][11:].strip() + '   Class 2\n')
    class_result.close()
def sel3():
    class_result=open('c:/users/a/classfication.txt','a+')
    class_result.write(imagelist[i][11:].strip() + '   Class 3\n')
    class_result.close()
def sel4():
    class_result=open('c:/users/a/classfication.txt','a+')
    class_result.write(imagelist[i][11:].strip() + '   Class 4\n')
    class_result.close()
def sel5():
    class_result=open('c:/users/a/classfication.txt','a+')
    class_result.write(imagelist[i][11:].strip() + '   Class 5\n')
    class_result.close()

imagefile=open('c:/users/a/imagelist.txt','r')
imagelist=imagefile.readlines()
imagefile.close()
top=Tk()
img=Image.open(imagelist[i].strip())
canvas=Canvas(top,height=100,width=100)
image=ImageTk.PhotoImage(img,top)
imgs=canvas.create_image(50,50,image=image)

var = StringVar()
label = Message( top, textvariable=var, relief=RAISED )
var.set(imagelist[i][11:].strip())

var_r = IntVar()
R1 = Radiobutton(top, text='Class 1', variable=var_r, value=1, command=sel1)
R2 = Radiobutton(top, text='Class 2', variable=var_r, value=2, command=sel2)
R3 = Radiobutton(top, text='Class 3', variable=var_r, value=3, command=sel3)
R4 = Radiobutton(top, text='Class 4', variable=var_r, value=4, command=sel4)
R5 = Radiobutton(top, text='Class 5', variable=var_r, value=5, command=sel5)
R6 = Radiobutton(top, text='Please choose a class above!', variable=var_r, value=6)
R1.pack()
R2.pack()
R3.pack()
R4.pack()
R5.pack()
R6.pack()

button_prev=Button(top,text='Prev',command=showimage_prev)
button_next=Button(top,text='Next',command=showimage_next)

button_prev.pack()
button_next.pack()

canvas.pack()
label.pack()

top.mainloop()
