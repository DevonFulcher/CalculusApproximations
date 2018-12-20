from sympy import *
from sympy.parsing.sympy_parser import parse_expr

x, y, n = symbols('x y n')

'''
returns the approximation of a differential equation (de) at a given x value
differential = de that will be approximated
estimate_at = the x value of the equation that is to be approximated
step_size = the change of the x value at each step
initial_x = initial given x value of the equation
initial_y = initial given y value of the equation
'''
def eulers_method(differential, estimate_at, step_size, initial_x, initial_y):
    #step_of_x is the x value at any given step
    step_of_x = initial_x
    #step_of_x is the y value at any given step
    step_of_y = initial_y
    while step_of_x < estimate_at:
        #the the x and y value of each step are subbed into the de
        subbed_in_differential = differential.subs(x, step_of_x)
        subbed_in_differential = subbed_in_differential.subs(y, step_of_y)
        
        #step of x and y are updated
        step_of_y = step_of_y + step_size * subbed_in_differential
        step_of_x = step_of_x + step_size
    return step_of_y

'''
returns the approximation of a convergent alternating series
f = a function that would produce a convergent alternating series
num_decimal_points = expected precision of approximation
'''
def approximate_alternating_series(f, num_decimal_points):
    #get rid of (-1)**n in the numerator
    f_without_neg_one_to_the_n = f / sympify((-1)**n)
    #sub in n+1 for all n in the function
    f_of_n_plus_one = f_without_neg_one_to_the_n.subs(n, n+1)
    precision = 1
    for _ in range(num_decimal_points):
        precision = precision / 10
    terms = 0
    check_num_terms = f_of_n_plus_one.subs(n, terms)
    while check_num_terms > precision:
        terms = terms + 1
        check_num_terms = f_of_n_plus_one.subs(n, terms)
    approximation = 0
    for term in range(terms):
        approximation = approximation + f.subs(n, term)
    return approximation.evalf()

'''
returns the approximation of an integral
f = the function whos integral is to be approximated
lower_bound = the lower bound of the integral
upper_bound = the upper bound of the integral
num_intervals = the number of intervales to calculate underneath the curve
'''
def simpsons_rule(f, lower_bound, upper_bound, num_intervals):
    approximation = f.subs(x, lower_bound)
    delta_x = (upper_bound - lower_bound) / num_intervals
    this_x_value = lower_bound + delta_x
    #the second coefficient of the formula is 4
    coefficient = 4
    for _ in range(num_intervals - 1):
        approximation = approximation + coefficient * f.subs(x, this_x_value)
        this_x_value = this_x_value + delta_x
        #the coefficients of the formula then alternate between 2 and 4 except the last one
        if coefficient == 4:
            coefficient = 2
        else:
            coefficient = 4
    return (delta_x / 3 * (approximation + f.subs(x, upper_bound))).evalf()


def main():
    while true:
        print("\nWould you like to approximate the value of a function using Euler's method,\n"
              "approximate an alternating series that starts at n=1,\n"
              "or approximate the integral of a function using Simpson's rule?")
        print("1) Euler's method")
        print("2) alternating series")
        print("3) Simpson's rule")
        response = input()
        if response == "1":
            differential = parse_expr(input("(e.g. -x-y) y' = "))
            estimate_at = float(input("estimate y at: "))
            step_size = float(input("step size: "))
            initial_x = float(input("initial x: "))
            initial_y = float(input("initial y: "))
            print("answer: ", eulers_method(differential, estimate_at, step_size, initial_x, initial_y))
        elif response == "2":
            f = parse_expr(input("enter a function that would produce a convergent alternating series "
                                 "(e.g. ((-1)**n * n)/3**n): f(n)= "))
            num_decimal_points = int(input("how many decimal points do you want the approximation to be within? "))
            print("answer: ", approximate_alternating_series(f, num_decimal_points))
        else:
            f = parse_expr(input("(e.g. sin(x)**2) f(x)= "))
            lower_bound = float(input("lower bound: "))
            upper_bound = float(input("upper bound: "))
            num_intervals = int(input("number of intervals: "))
            print("answer: ", simpsons_rule(f, lower_bound, upper_bound, num_intervals))


if __name__ == '__main__':
    main()
