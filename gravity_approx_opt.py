import matplotlib.pyplot as plt 

# this program determines the best possible approximations for gravity
# by simulating the following:
# find the force of attraction of a particle to two other particles, and add them
# then, find the force of attraction of a particle to the center of mass of those particles
# determine what is the highest distance from the center of mass given a certain "n"
# where "n" is half of the side length of the square on the "approximation grid"
# therefore, for grids of l = 100 and T = 0.90, the particle needs to be at about 300 squares away 
# this seems to have a linear relationship with a slope of 3 for T = 0.90

T = 0.90  # accuracy threshold

def f(x: int) -> float:
    return 1 / x**2

def g(x: int, n: int) -> float: 
    return (2 * f(x)) / (f(x - n) + f(x + n))

if __name__ == "__main__":
    test_n = [1, 5, 10, 20, 50, 100]
    ans = []
    for n in test_n:
        for j in range(n + 1, 1000):
            if g(j, n // 2) >= T:
                break
        ans.append(j)
    print(*zip(test_n, ans), sep="\n")
    plt.scatter(test_n, ans)
    plt.show()
