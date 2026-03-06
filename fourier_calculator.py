from sympy import *
import numpy as np


def calculate_fourier_coefficients(func_str: str, lower_bound: float, upper_bound: float, degree: int, speed: int):
    x = symbols('x')
    try:
        f = sympify(func_str)
    except SympifyError as e:
        raise ValueError(f"Error parsing function string: '{func_str}'") from e

    def hermite_interpolation(node_info):
        nodes, divided_diff = [], []
        for key in node_info:
            nodes += [key] * len(node_info[key])
        n = len(nodes)
        for i in range(n):
            divided_diff.append([S.Zero]*(n-i))  
            divided_diff[i][0] = node_info[nodes[i]][0]        
        
        for j in range(1, n):
            for i in range(n - j):
                same_nodes = all(Eq(nodes[i], nodes[i+k]) for k in range(j+1))
                if same_nodes:
                    divided_diff[i][j] = node_info[nodes[i]][j] / factorial(j)
                else:
                    numerator = divided_diff[i+1][j-1] - divided_diff[i][j-1]
                    denominator = nodes[i+j] - nodes[i]  
                    divided_diff[i][j] = numerator / denominator

        P, product = Integer(0), Integer(1)
        for j in range(n):
            P += divided_diff[0][j] * product
            if j < n-1:
                product *= (x - nodes[j])  
                
        return P

    def simpson38(func, a, b, n_segments):
        if n_segments % 3 != 0:
            raise ValueError("n_segments must be a multiple of 3")
        h = (b - a) / n_segments
        total = 0.0
        for i in range(0, n_segments, 3):
            x0 = a + i*h
            x1 = x0 + h
            x2 = x0 + 2*h
            x3 = x0 + 3*h
            total += func(x0) + 3*func(x1) + 3*func(x2) + func(x3)
        return 3*h/8 * total
    
    T = 1 / 16
    left = -T / 4
    delta = (upper_bound + lower_bound) / 2
    z = (x - delta) * 2 * (upper_bound - lower_bound) / T
    f = f.subs(x, z)
    f = 0.5 * f

    mid, right = left + 1 * T / 2, left + T
    N = degree

    points = {mid: [], right: []}
    ω = Integer(2) * pi / T
    for i in range(speed):
        f_ = diff(f, x, i)
        points[mid].append(f_.subs(x, mid))
        points[right].append(f_.subs(x, left))

    P = hermite_interpolation(points)

    num_f = lambdify(x, f, "numpy")     
    num_P = lambdify(x, P, "numpy")     
    N_SEG = 300
    N_SEG2 = 300
    I1 = simpson38(num_f, float(left), float(mid), N_SEG)
    I2 = simpson38(num_P, float(mid), float(right), N_SEG2)
    a0_num = (I1 + I2) * (1.0 / float(T))

    num_an = [a0_num]
    num_bn = [0]

    for ni in range(1, N+1):
        ωi = float(ω) * ni
        c1 = simpson38(lambda t: num_f(t) * np.cos(ωi * t),
                float(left), float(mid), N_SEG)
        s1 = simpson38(lambda t: num_f(t) * np.sin(ωi * t),
                    float(left), float(mid), N_SEG)
        c2 = simpson38(lambda t: num_P(t) * np.cos(ωi * t),
                    float(mid),  float(right), N_SEG2)

        s2 = simpson38(lambda t: num_P(t) * np.sin(ωi * t),
                    float(mid),  float(right), N_SEG2)

        an_i = 2.0/float(T) * (c1 + c2)
        bn_i = 2.0/float(T) * (s1 + s2)
        num_an.append(an_i)
        num_bn.append(bn_i)

    result = [complex(an, -bn) for an, bn in zip(num_an, num_bn)]
    return result
