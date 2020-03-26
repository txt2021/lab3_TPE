import math

x1_min = 10
x1_max = 60
x2_min = -25
x2_max = 10
x3_min = 10
x3_max = 15
y_min = 200 + int((x1_min + x2_min + x3_min) / 3)
y_max = 200 + int((x1_max + x2_max + x3_max) / 3)
m = 3
N = 4
p = 0.95

def genermat():
    from random import randrange
    matrix_with_y = [[randrange(y_min, y_max) for y in range(3)] for x in range(4)]
    return matrix_with_y


def aver_y(lst):
    average = []
    for k in range(len(lst)):
        average.append(sum(lst[k]) / len(lst[k]))
    return average


def aver_x(lst):
    average = [0, 0, 0]
    for k in range(4):
        average[0] += lst[k][0] / 4
        average[1] += lst[k][1] / 4
        average[2] += lst[k][2] / 4
    return average


def det(a):
    from numpy.linalg import det
    return det(a)


class Critical:

    def get_cohren_value(size_of_selections, qty_of_selections, significance):
        from _pydecimal import Decimal
        from scipy.stats import f
        size_of_selections += 1
        partResult1 = significance / (size_of_selections - 1)
        params = [partResult1, qty_of_selections, (size_of_selections - 1 - 1) * qty_of_selections]
        fisher = f.isf(*params)
        result = fisher / (fisher + (size_of_selections - 1 - 1))
        return Decimal(result).quantize(Decimal('.0001')).__float__()

    def get_student_value(f3, significance):
        from _pydecimal import Decimal
        from scipy.stats import t
        return Decimal(abs(t.ppf(significance / 2, f3))).quantize(Decimal('.0001')).__float__()


    def get_fisher_value(f3, f4, significance):
        from _pydecimal import Decimal
        from scipy.stats import f
        return Decimal(abs(f.isf(significance, f4, f3))).quantize(Decimal('.0001')).__float__()

print("\tРівняння регресії: \n\tŷ = b0 + b1*X1 + b2*X2 + b3*X3")
print("\tМатриця планування експеременту")
matrix_pfe = [[1, -1, -1, -1],
              [1, -1, 1, 1],
              [1, 1, -1, 1],
              [1, 1, 1, -1]]

for i in range(len(matrix_pfe)):
    print(" ", end=" ")
    for j in range(len(matrix_pfe[i])):
        print(matrix_pfe[i][j], end=" ")
    print(" ")


a = True
while a:
    y_matrix = genermat()
    x_matrix = [[x1_min, x2_min, x3_min], [x1_min, x2_max, x3_max], [x1_max, x2_min, x3_max], [x1_max, x2_max, x3_min]]
    matrix = []
    average_y = aver_y(y_matrix)
    average_x = aver_x(x_matrix)
    a1, a2, a3, a11, a22, a33, a12, a13, a23 = 0, 0, 0, 0, 0, 0, 0, 0, 0
    for i in range(4):
        a1 += x_matrix[i][0] * average_y[i] / 4
        a2 += x_matrix[i][1] * average_y[i] / 4
        a3 += x_matrix[i][2] * average_y[i] / 4
        a11 += x_matrix[i][0] ** 2 / 4
        a22 += x_matrix[i][1] ** 2 / 4
        a33 += x_matrix[i][2] ** 2 / 4
        a12 += x_matrix[i][0] * x_matrix[i][1] / 4
        a13 += x_matrix[i][0] * x_matrix[i][2] / 4
        a23 += x_matrix[i][1] * x_matrix[i][2] / 4

    a21 = a12
    a31 = a13
    a32 = a23
    my = sum(average_y) / len(average_y)
    b0_numerator = [[my, average_x[0], average_x[1], average_x[2]], [a1, a11, a12, a13], [a2, a21, a22, a23],
                    [a3, a31, a32, a33]]
    b1_numerator = [[1, my, average_x[1], average_x[2]], [average_x[0], a1, a12, a13], [average_x[1], a2, a22, a23],
                    [average_x[2], a3, a32, a33]]
    b2_numerator = [[1, average_x[0], my, average_x[2]], [average_x[0], a11, a1, a13], [average_x[1], a21, a2, a23],
                    [average_x[2], a31, a3, a33]]
    b3_numerator = [[1, average_x[0], average_x[1], my], [average_x[0], a11, a12, a1], [average_x[1], a21, a22, a2],
                    [average_x[2], a31, a32, a3]]
    b_denominator = [[1, average_x[0], average_x[1], average_x[2]], [average_x[0], a11, a12, a13],
                     [average_x[1], a21, a22, a23], [average_x[2], a31, a32, a33]]
    b0 = det(b0_numerator) / det(b_denominator)
    b1 = det(b1_numerator) / det(b_denominator)
    b2 = det(b2_numerator) / det(b_denominator)
    b3 = det(b3_numerator) / det(b_denominator)
    f1 = m - 1
    f2 = N
    q = 1 - p
    dispersion_y = [0, 0, 0, 0]
    for i in range(m):
        dispersion_y[0] += ((y_matrix[0][i] - average_y[0]) ** 2) / 3
        dispersion_y[1] += ((y_matrix[1][i] - average_y[1]) ** 2) / 3
        dispersion_y[2] += ((y_matrix[2][i] - average_y[2]) ** 2) / 3
        dispersion_y[3] += ((y_matrix[3][i] - average_y[3]) ** 2) / 3

    Gp = max(dispersion_y) / sum(dispersion_y)
    print("\tПеревірка за критерієм Кохрена")
    Gt = Critical.get_cohren_value(f2, f1, q)
    if Gt > Gp:
        print("\tДисперсія однорідна при рівні значимості {:.2f}\n\t".format(q))
        a = False
    else:
        print("\tДисперсія не однорідна при рівні значимості {:.2f}!".format(q))
        m += 1

for i in range(N):                                    #4 -кількість експериментів (рядків матриці планування), замість 4 можна написати N
    matrix.append(x_matrix[i] + y_matrix[i])
print("\tМатриця з натуральних значень факторів")
print("  X1 X2 X3 Y1  Y2  Y3 ")
for i in range(len(matrix)):
    print(" ", end=" ")
    for j in range(len(matrix[i])):
        print(matrix[i][j], end=" ")
    print(" ")

print("\tРівняння регресії")
print("{:.3f} + {:.3f} * X1 + {:.3f} * X2 + {:.3f} * X3 = ŷ".format(b0, b1, b2, b3))
print("\tПеревірка")
print("{:.3f} + {:.3f} * {:.3f} + {:.3f} * {:.3f} + {:.3f} * {:.3f} = ".format(b0, b1, x1_min, b2, x2_min, b3, x3_min)
      + str(b0 + b1 * x1_min + b2 * x2_min + b3 * x3_min))
print("{:.3f} + {:.3f} * {:.3f} + {:.3f} * {:.3f} + {:.3f} * {:.3f} = ".format(b0, b1, x1_min, b2, x2_max, b3, x3_max)
      + str(b0 + b1 * x1_min + b2 * x2_max + b3 * x3_max))
print("{:.3f} + {:.3f} * {:.3f} + {:.3f} * {:.3f} + {:.3f} * {:.3f} = ".format(b0, b1, x1_max, b2, x2_min, b3, x3_max)
      + str(b0 + b1 * x1_max + b2 * x2_min + b3 * x3_max))
print("{:.3f} + {:.3f} * {:.3f} + {:.3f} * {:.3f} + {:.3f} * {:.3f} = ".format(b0, b1, x1_max, b2, x2_max, b3, x3_min)
      + str(b0 + b1 * x1_max + b2 * x2_max + b3 * x3_min))

print("\tПеревірка за критерієм Стьюдента")
f3 = f1 * f2
S_2b = sum(dispersion_y) / (N * N * m)
S_b = math.sqrt(S_2b)
beta_0 = (average_y[0] + average_y[1] + average_y[2] + average_y[3]) / N
beta_1 = (-average_y[0] - average_y[1] + average_y[2] + average_y[3]) / N
beta_2 = (-average_y[0] + average_y[1] - average_y[2] + average_y[3]) / N
beta_3 = (-average_y[0] + average_y[1] + average_y[2] - average_y[3]) / N
t_0 = math.fabs(beta_0) / S_b
t_1 = math.fabs(beta_1) / S_b
t_2 = math.fabs(beta_2) / S_b
t_3 = math.fabs(beta_3) / S_b
Tt = Critical.get_student_value(f1 * f2, q)
t_lst = [t_0, t_1, t_2, t_3]
b_lst = [b0, b1, b2, b3]
for i in range(4):                
    if t_lst[i] > Tt:
        continue
    else:
        t_lst[i] = 0
for j in range(4):
    if t_lst[j] != 0:
        continue
    else:
        b_lst[j] = 0
print("\t\tЗначимі коефіцієнти: ")
yj1 = b_lst[0] + b_lst[1] * x1_min + b_lst[2] * x2_min + b_lst[3] * x3_min
yj2 = b_lst[0] + b_lst[1] * x1_min + b_lst[2] * x2_max + b_lst[3] * x3_max
yj3 = b_lst[0] + b_lst[1] * x1_max + b_lst[2] * x2_min + b_lst[3] * x3_max
yj4 = b_lst[0] + b_lst[1] * x1_max + b_lst[2] * x2_max + b_lst[3] * x3_min
print("{:.3f} + {:.3f} * {:.3f} + {:.3f} * {:.3f} + {:.3f} * {:.3f} = "
      "".format(b_lst[0], b_lst[1], x1_min, b_lst[2], x2_min, b_lst[3],
                x3_min) + str(yj1))
print("{:.3f} + {:.3f} * {:.3f} + {:.3f} * {:.3f} + {:.3f} * {:.3f} = "
      "".format(b_lst[0], b_lst[1], x1_min, b_lst[2], x2_max, b_lst[3],
                x3_max) + str(yj2))
print("{:.3f} + {:.3f} * {:.3f} + {:.3f} * {:.3f} + {:.3f} * {:.3f} = "
      "".format(b_lst[0], b_lst[1], x1_max, b_lst[2], x2_min, b_lst[3],
                x3_max) + str(yj3))
print("{:.3f} + {:.3f} * {:.3f} + {:.3f} * {:.3f} + {:.3f} * {:.3f} = "
      "".format(b_lst[0], b_lst[1], x1_max, b_lst[2], x2_max, b_lst[3],
                x3_min) + str(yj4))
print("\tПеревірка за критерієм Фішера")
for i in range(3):
    if b_lst[i] == 0:
        del b_lst[i]

while True:
 d = len(b_lst)
 f4 = N - d
 S_2ad = m * ((yj1 - average_y[0]) ** 2 + (yj2 - average_y[1]) ** 2 + (yj3 - average_y[2]) ** 2 + (
             yj4 - average_y[3]) ** 2) / f4
 Fp = S_2ad / S_2b
 Ft = Critical.get_fisher_value(f1 * f2, f4, q)
 print(Fp)
 print(Ft)
 if Fp > Ft:
     print("Рівняння регресії неадекватно оригіналу при рівні значимості {:.2f}".format(q))
     m=m+1       #але все одно рівняння ніколи не буде адекватним, адже при збільшенні m Fp теж збільшується
 else:
     print("Рівняння регресії адекватно оригіналу при рівні значимості {:.2f}".format(q))
     break         
