from flask import Flask, render_template, request
from sympy import symbols, sympify
import numpy as np
import matplotlib.pyplot as plt

app = Flask(__name__)

def temp(T):
    T[0][1] = 976
    T[0][2] = 904
    T[0][3] = 784
    T[0][4] = 616
    T[0][5] = 400


@app.route('/')
def index():
    return render_template('form.html')

@app.route('/calculate', methods=['POST'])
def calculate():
    a_2 = float(request.form['a'])
    n = int(request.form['n'])
    m = int(request.form['m'])
    l = float(request.form['l'])
    t0 = float(request.form['t0'])

    h = l / n
    tau = t0 / m

    sigma = a_2 / (tau * pow(h, 2))
    a = sigma
    b = sigma
    c = 1 + 2 * sigma

    t = symbols('t')
    x = symbols('x')

    phi_expr = sympify(request.form['phi'])
    psi_expr = sympify(request.form['psi'])
    f_expr = sympify(request.form['f'])

    T = np.zeros((m + 1, n + 1))

    #крайові та граничні умови
    for i in range(m + 1):
        T[i][0] = phi_expr.subs(t, tau*i)
        T[i][n] = psi_expr.subs(t, tau*i)
    for i in range(n + 1):
        T[0][i] = f_expr.subs(x, i*h)



    # прогонка коефіціентів
    for k in range(1, m + 1):
        alpha = [0]
        beta = [T[k][0]]
        for i in range(1, n + 1):
            alpha_i = sigma / (1 + (2 - alpha[i - 1]) * sigma)
            alpha.append(alpha_i)
            beta_i = (sigma * beta[i - 1] + T[k - 1][i - 1]) / (1 + (2 - alpha[i - 1]) * sigma)
            beta.append(beta_i)
        # зворотній хід
        for i in range(n - 1, 0, -1):
            T[k][i] = alpha[i + 1] * T[k][i + 1] + beta[i + 1]
            pass

    plot(l,t0,n,m,T)

    plt.savefig('static/plot.png')
    return render_template('result.html', result=T)

# def start_cond(Matrix,phi,psi,f)

def plot(l,t0,n,m,T):
    # Створення сітки
    x = np.linspace(0, l, n + 1)
    t = np.linspace(0, t0, m + 1)
    X, T_func = np.meshgrid(x, t)

    # Побудова графіка
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, T_func, T, cmap='viridis')

    # Налаштування відображення
    ax.set_xlabel('x')
    ax.set_ylabel('t')
    ax.set_zlabel('T')
    ax.set_title('Approximation of T(x, t)')


if __name__ == '__main__':
    app.run(debug=True)


