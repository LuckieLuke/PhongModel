from tkinter import *
import math
import numpy as np


def get_z(x, y, R):
    return math.sqrt(R**2 - (x-450)**2 - (y-350)**2)


def illuminate(kd, ks, n, x_l, y_l, z_l, R_b):
    for x in range(int(450-R_b), int(450+R_b)):
        for y in range(int(350-R_b), int(350+R_b)):
            if (x-450)**2 + (y-350)**2 > R_b**2:
                pass
            else:
                z = get_z(x, y, R_b)
                N = np.array([x, y, z]) / np.linalg.norm(np.array([x, y, z]))
                L = np.array([x_l - x, y_l - y, z_l - z]) / \
                    np.linalg.norm(np.array([x_l - x, y_l - y, z_l - z]))
                beta = math.acos(N[0] * L[0] + N[1] * L[1] + N[2] * L[2] / (math.sqrt(
                    N[0]**2 + N[1]**2 + N[2]**2) * math.sqrt(L[0]**2 + L[1]**2 + L[2]**2)))
                R = L - 2 * N * \
                    math.cos(beta) / np.linalg.norm(L -
                                                    2 * N * math.cos(beta))
                V = np.array([0, 0, 1])
                I = 0.25 + (1 * (kd * N.dot(L) + ks * (R.dot(V)**n)))

                color = [90, 0, 0]
                for i in range(3):
                    new_value = int(color[i] + I * 255)
                    if new_value > 255:
                        new_value = 255
                    elif new_value < 0:
                        new_value = 0
                    color[i] = new_value
                new_color = "#%02x%02x%02x" % tuple(color)

                my_canvas.create_line(x, y, x+1, y+1, fill=new_color)


root = Tk()
root.geometry('1000x800')
my_canvas = Canvas(root, width=900, height=700, bg='white')
my_canvas.pack(pady=10)

R = 100

illuminate(kd=0.8, ks=0.1, n=400, x_l=420, y_l=320, z_l=-200, R_b=R)
# illuminate(kd=0.53, ks=0.92, n=100, x_l=450, y_l=350, z_l=-200, R_b=R)
# illuminate(kd=0.85, ks=0.15, n=100, x_l=450, y_l=350, z_l=-200, R_b=R)
# illuminate(kd=0.44, ks=0.71, n=900, x_l=450, y_l=350, z_l=-200, R_b=R)

root.mainloop()
