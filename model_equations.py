
def Adot(A,R,G,C,D,B):
    return 1-A+G*A**2*(1-(A/C))-R*A
def Rdot (A,R,G,C,D,B):
    return D-D*R+B*A*R
def model (y,t,G,C,D,B):
    A,R = y
    Adot = 1-A+G*A**2*(1-(A/C))-R*A
    Rdot = D-D*R+B*A*R
    return [Adot,Rdot]

def A_nullcline(x,C,G,D,B):
    return (1+G*(x**2)*(1-(x/C)))/(x)
def R_nullcline(x,C,G,D,B):
    return D/(D-B*x)

