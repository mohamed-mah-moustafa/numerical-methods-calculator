import sympy as sp
import math

tol = 0.01

# --- HELPER FUNCTIONS ---

def func():
    """Scans a function from the user and returns a machine (f) and blueprint (expr)."""
    x = sp.symbols('x')
    user_input = input("Enter your function f(x) (use ** for power, exp(x) for e): ")
    expr = sp.sympify(user_input)
    f = sp.lambdify(x, expr)
    return f, expr

def get_df(expr):
    """Calculates and returns the numerical derivative function."""
    x = sp.symbols('x')
    deriv_expr = sp.diff(expr, x)
    return sp.lambdify(x, deriv_expr)

def get_points():
    """Scans a set of (x, y) coordinates from the user."""
    n = int(input("Enter the number of points: "))
    x_pts = []
    y_pts = []
    for i in range(n):
        x_pts.append(float(input(f"Enter x{i}: ")))
        y_pts.append(float(input(f"Enter y{i}: ")))
    return x_pts, y_pts


# --- PART A: ROOT FINDING METHODS ---

def bisection(): 
    f, _ = func()
    a = float(input("enter value of a: "))
    b = float(input("enter value of b: "))
    c = (a + b) / 2
    
    if f(a) * f(b) >= 0:
        print("cant do operation")
        return
    
    print(f"\n{'a':<10} | {'b':<10} | {'c':<10} | {'f(a)':<10}| {'f(b)':<10}| {'f(c)':<10}")
    print("-" * 74)
    
    while abs(f(c)) > tol:
        print(f"{a:<10.4f} | {b:<10.4f} | {c:<10.4f} | {f(a):<10.4f}| {f(b):<10.4f}| {f(c):<10.4f}")
        if f(a) * f(c) > 0:
            a = c
        else:
            b = c
        c = (a + b) / 2
    print(f"the root is approximately = ({c:.4f}, {f(c):.4f})")

def sec(num):
    f, _ = func()
    a = float(input("enter value of a: "))
    b = float(input("enter value of b: "))
    c = b - f(b) * ((b - a) / (f(b) - f(a))) 
    print(f"\n{'a':<10} | {'b':<10} | {'c':<10} | {'f(a)':<10}| {'f(b)':<10}| {'f(c)':<10}")
    print("-" * 74)
    
    while abs(f(c)) > tol:
        print(f"{a:<10.4f} | {b:<10.4f} | {c:<10.4f} | {f(a):<10.4f}| {f(b):<10.4f}| {f(c):<10.4f}")
        if num == 1:
            a = b
            b = c
        else: # Modified Secant
            if f(a) * f(c) > 0:
                a = c
            else:
                b = c
                
        c = b - f(b) * ((b - a) / (f(b) - f(a))) 
    print(f"{a:<10.4f} | {b:<10.4f} | {c:<10.4f} | {f(a):<10.4f}| {f(b):<10.4f}| {f(c):<10.4f}")
    print(f"the root is approximately = ({c:.4f}, {f(c):.4f})")
    
def newton():
    f, expr = func()
    df = get_df(expr)
    x = float(input("enter value of x: "))
    print(f"\n{'x':<10} | {'f(x)':<10} | {'df(x)':<10}")
    print("-" * 36)

    while abs(f(x)) > tol:
        print(f"{x:<10.4f} | {f(x):<10.4f} | {df(x):<10.4f}")
        x = x - (f(x) / df(x))
                
    print(f"{x:<10.4f} | {f(x):<10.4f} | {df(x):<10.4f}")
    print(f"the root is approximately = ({x:.4f}, {f(x):.4f})")

def fixed():
    f, _ = func()
    x = float(input("enter value of x: "))
    print(f"\n{'x':<10} | {'f(x)':<10} | {'g(x)':<10}")
    print("-" * 36)
    
    while abs(f(x)) > tol:
        fx = f(x)
        gx = fx + x
        print(f"{x:<10.4f} | {f(x):<10.4f} | {gx:<10.4f}")
        x = gx

    print(f"{x:<10.4f} | {f(x):<10.4f} | {gx:<10.4f}")
    print(f"the root is approximately = ({x:.4f}, {f(x):.4f})")


# --- PART B: INTERPOLATION & REGRESSION ---

def lagrange():
    x_points, y_points = get_points()
    n = len(x_points)
    x = sp.symbols('x')
    final_poly = 0
    
    for i in range(n):
        li = 1
        for j in range(n):
            if i != j:
                li = li * (x - x_points[j]) / (x_points[i] - x_points[j])
        print(f"L{i}(x) = {sp.simplify(li)}")
        final_poly += y_points[i] * li

    print("\n--- Final Results ---")
    simplified_poly = sp.simplify(final_poly)
    print(f"Final Polynomial P(x) = {simplified_poly}")
    
    target = float(input("\nEnter an x value to interpolate: "))
    result = simplified_poly.subs(x, target)
    print(f"P({target}) = {result:.4f}")

def newton_divided_diff():
    x_pts, y_pts = get_points()
    n = len(x_pts)
    table = [[0] * n for _ in range(n)]
    
    for i in range(n):
        table[i][0] = y_pts[i]
        
    for j in range(1, n):
        for i in range(n - j):
            table[i][j] = (table[i+1][j-1] - table[i][j-1]) / (x_pts[i+j] - x_pts[i])
            
    print("\n--- Divided Difference Table ---")
    header = "xi      | f[xi]   " + " ".join([f"Order {i:<2}" for i in range(1, n)])
    print(header)
    print("-" * len(header) * 2)
    for i in range(n):
        row = f"{x_pts[i]:<7.2f} | "
        for j in range(n - i):
            row += f"{table[i][j]:<9.4f} "
        print(row)

    target_x = float(input("\nEnter x value to interpolate: "))
    res = table[0][0]
    product_term = 1
    for i in range(1, n):
        product_term *= (target_x - x_pts[i-1])
        res += table[0][i] * product_term
        
    print(f"Interpolated value at x={target_x} is {res:.6f}")

def linear_regression():
    x_pts, y_pts = get_points()
    n = len(x_pts)
    
    sum_x = sum(x_pts)
    sum_y = sum(y_pts)
    sum_xy = sum(x * y for x, y in zip(x_pts, y_pts))
    sum_x2 = sum(x**2 for x in x_pts)
    
    b = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x**2)
    a = (sum_y - b * sum_x) / n
    
    print(f"\nEquation: y = {b:.4f}x + {a:.4f}")
    print(f"Slope (b): {b:.4f}, Intercept (a): {a:.4f}")
    
    target_x = float(input("Enter x to predict y: "))
    print(f"For x = {target_x}, y = {b * target_x + a:.6f}")


# --- PART C: CALCULUS ---

def numerical_integration():
    f, _ = func()
    a = float(input("Enter lower bound a: "))
    b = float(input("Enter upper bound b: "))
    
    print("\n1. Simple Method")
    print("2. Composite Method")
    choice = input("Select integration type: ")

    if choice == '1':
        mid = (b - a) * f((a + b) / 2)
        trap = ((b - a) / 2) * (f(a) + f(b))
        simp = ((b - a) / 6) * (f(a) + f(b) + 4 * f((a + b) / 2))
        
        print(f"\n{'Method':<15} | {'Result':<10}")
        print("-" * 28)
        print(f"{'Midpoint':<15} | {mid:<10.6f}")
        print(f"{'Trapezoidal':<15} | {trap:<10.6f}")
        print(f"{'Simpson':<15} | {simp:<10.6f}")

    elif choice == '2':
        h = float(input("Enter step size h: "))
        n = int((b - a) / h)
        
        mid_sum = sum(f(a + (i + 0.5) * h) for i in range(n))
        mid_comp = h * mid_sum

        trap_sum = f(a) + f(b) + 2 * sum(f(a + i * h) for i in range(1, n))
        trap_comp = (h / 2) * trap_sum

        simp_sum = f(a) + f(b)
        for i in range(1, n):
            if i % 2 == 0:
                simp_sum += 2 * f(a + i * h)
            else:
                simp_sum += 4 * f(a + i * h)
        simp_comp = (h / 3) * simp_sum

        print(f"\nCalculated n = {n}")
        print(f"{'Method':<15} | {'Result':<10}")
        print("-" * 28)
        print(f"{'Midpoint':<15} | {mid_comp:<10.6f}")
        print(f"{'Trapezoidal':<15} | {trap_comp:<10.6f}")
        print(f"{'Simpson':<15} | {simp_comp:<10.6f}")

def numerical_diff():
    x_pts, y_pts = get_points()
    n = len(x_pts)
    h = x_pts[1] - x_pts[0]
    
    diff_table = [[0] * n for _ in range(n)]
    for i in range(n): 
        diff_table[i][0] = y_pts[i]
        
    for j in range(1, n):
        for i in range(n - j):
            diff_table[i][j] = diff_table[i+1][j-1] - diff_table[i][j-1]

    print("\n--- Unified Difference Table ---")
    header = f"{'x':<8} | {'y=f(x)':<10} | " + " | ".join([f"Diff^{i:<2}" for i in range(1, n)])
    print(header)
    print("-" * len(header))
    for i in range(n):
        row = f"{x_pts[i]:<8.2f} | {diff_table[i][0]:<10.4f} | "
        for j in range(1, n - i):
            row += f"{diff_table[i][j]:<9.4f} | "
        print(row)

    print(f"\n--- Forward Differentiation at x0 = {x_pts[0]:.2f} ---")
    dy_fwd = sum(((-1)**(k-1) / k) * diff_table[0][k] for k in range(1, n)) * (1 / h)
    print(f"First Derivative (dy/dx):   {dy_fwd:.6f}")
    
    if n >= 3:
        d2y_fwd = diff_table[0][2] 
        if n > 3: d2y_fwd -= diff_table[0][3]
        if n > 4: d2y_fwd += (11/12) * diff_table[0][4]
        if n > 5: d2y_fwd -= (5/6) * diff_table[0][5]
        d2y_fwd *= (1 / (h**2))
        print(f"Second Derivative (d2y/dx2): {d2y_fwd:.6f}")

    print(f"\n--- Backward Differentiation at xn = {x_pts[-1]:.2f} ---")
    dy_bwd = sum((1 / k) * diff_table[n - 1 - k][k] for k in range(1, n)) * (1 / h)
    print(f"First Derivative (dy/dx):   {dy_bwd:.6f}")
    
    if n >= 3:
        d2y_bwd = diff_table[n-3][2]
        if n > 3: d2y_bwd += diff_table[n-4][3]
        if n > 4: d2y_bwd += (11/12) * diff_table[n-5][4]
        if n > 5: d2y_bwd += (5/6) * diff_table[n-6][5]
        d2y_bwd *= (1 / (h**2))
        print(f"Second Derivative (d2y/dx2): {d2y_bwd:.6f}")


# --- MAIN MENU ---

if __name__ == "__main__":
    while True:
        print("\n" + "="*35)
        print("      NUMERICAL METHODS MENU")
        print("="*35)
        print("1. Bisection")
        print("2. Secant")
        print("3. Modified Secant")
        print("4. Newton-Raphson")
        print("5. Fixed Point")
        print("6. Lagrange Interpolation")
        print("7. Newton Divided Difference")
        print("8. Linear Regression")
        print("9. Numerical Integration")
        print("10. Numerical Differentiation")
        print("11. Exit")
        print("="*35)
            
        choice = input("Select an option: ")
            
        if choice == '1': bisection()
        elif choice == '2': sec(1)
        elif choice == '3': sec(2)
        elif choice == '4': newton()
        elif choice == '5': fixed()
        elif choice == '6': lagrange()    
        elif choice == '7': newton_divided_diff()
        elif choice == '8': linear_regression()
        elif choice == '9': numerical_integration()
        elif choice == '10': numerical_diff()
        elif choice == '11': 
            print("Exiting...")
            break
        else: 
            print("Invalid choice! Please select a valid number.")