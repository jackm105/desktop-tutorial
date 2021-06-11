import numpy as np
import scipy.linalg as la
import scipy.linalg
import math
import sys
#jack: this is just for the first question the matrix is 13x13 and is in another file
A = np.genfromtxt('bridge_matrix.csv', delimiter=',') ### Don't delete this line

################################################################################


testA = np.array([[2,1,1],[4,3,3],[8,7,9]])
testb = np.array([[1],[1],[-1]])
testx = scipy.linalg.solve(testA,testb)

### Problem 1
### A is already initialized above from a separate file (don't delete line 4).
### Initialize the data (right hand side) b.
b = np.array([[0 , 0 ,  0 , 0 , 0 , 0 , 0 , 0 , 3 , 0 , np.e**2 , 0 , np.pi]])
b = np.reshape(b , (-1 , 1))



### Solve for the force vector and save it as A1

#jack this is the main function that takes A and gives P,L,U for the PA=LU form
P, L, U = scipy.linalg.lu(A)
#jack this is because python actually does A=PLU but P is just a reordered identity matrix to the inverse is just the transpose so that transposes it so PA=LU works
P = P.T
#jack sometimes they wand Ux=y and this gives you y
y = scipy.linalg.solve_triangular(L , P @ b ,lower=True)
#jack this gives you the x which is what you susually want
x = scipy.linalg.solve_triangular(U , y)
c = U @ x
#jack this is how i imput answers
A1 = x.copy()


### Compute the PA = LU factorization of A
### You may want to use some of these variables later on so don't forget to
### use .copy() wherever appropriate
### Save L AS A2, and c as A3.
A2 = L.copy()
A3 = c.copy()


### Create a loop that breaks when one of the forces is greater than 20 tons
### Save A4 as the weight of the truck in position 8
### Save A5 as the entry of the force vector that exceeds 20 tons

#jack so basically it wants me to change the 8th entry in b and add 0.001 to it until any of the entries
#in x are greater than 20 and then it wants the value of b that broke it and what the entry of x that was greater than 20
#so it takes b[8] which is the 8th entry in b and gets rid of it and then adds 3 + i * 0.001 and then it does the whole
#function stuff again until a value in x exceeds 20 but for some reason it doesn't work
for i in range(17005) :
    b[8] = (b[8] * 0) + 3 + (i * 0.001)
    P, L, U = scipy.linalg.lu(A)
    P = P.T
    y = scipy.linalg.solve_triangular(L, P @ b, lower=True)
    x = scipy.linalg.solve_triangular(U, y)
    if np.max(abs(x)) > 20 :
        break

A4 = b[8].copy()-0.001
A5 = np.argmax(x) + 1



### Problem 2
### Initialize, alpha, omega, and A, and compute the PA = LU

#jack this is the one that im so sure i did right and isnt correct. locationx y and d are vectors full of zeros 1001 long
#then i set the first values of each according to the problem, then i set A and b which is called X0
#Then theres a function that takes a 2d vector b and gives you the 2d vector x like before
#The first for loop makes a temp 2d vector of the first value being the last entry in tjhe big ass locationx vector and
#then the second value as the last value of the location y vector and then assigns the current location x value as the
#first value of the new 2d vector that the function outputs when you input the temp matrix and the same for the locationy value
#but once again for some reason this doesnt work
total = 1001
locationx = np.zeros(total)
locationy = np.zeros(total)
locationd = np.zeros(total)
locationx[0] = 1
locationy[0] = -1
alpha = -0.002
omega = 0.06
X0 = np.array([[1] , [-1]])
A2 = np.array([[1-alpha , -omega] , [omega , 1-alpha]])


def zoro(location):
    P2, L2, U2 = scipy.linalg.lu(A2)
    P2 = P2.T
    y2 = scipy.linalg.solve_triangular(L2, P2 @ location, lower=True)
    x2 = scipy.linalg.solve_triangular(U2, y2)
    return x2


for j in range(1 , total):
    tempa = locationx[j-1]
    tempb = locationy[j-1]
    temp2 = np.array([[tempa] , [tempb]])
    locationx[j] = zoro(temp2)[0]
    locationy[j] = zoro(temp2)[1]


A6 = locationx.copy()
A7 = locationy.copy()


### The initializations can get a little tricky so definitely ask for help
### if you're stuck.
### Initialize a matrix made up of the position vector at each time
### Set the first x and y coordinates at time = 0 in your matrix
### to the values instructed in the assignment file.
### Create a loop that loops through each time given in the assignment file.
### Compute the new right hand side c using P, L, and/or U.
### You may need to recall that the inverse of P is P transpose
### Solve for the position by solving the Ux = c equation.
### Save all x coordinates as A6
### Save all y coordinates as A7
### Save the distance from the origin as A8

#jack this is part of the same problem and this loop just takes whatever value of locationd that were at and makes it
#equal to the square root of the sum of the x and y values squared or the distance essentially from the origin so
#now location d is just the distance at each value for the x and y from the origin
#and once again this doesnt work
for j in range(total):
    locationd[j] = math.sqrt((locationx[j] ** 2) + (locationy[j] ** 2))

A8 = locationd.copy()
np.set_printoptions(threshold=sys.maxsize)
print(A8)

### Initialize a position vector
### Initialize a distance variable
### Initialize a time variable
### Create a loop that breaks when the distance from the origin is
### less than 0.06.
### In the loop compute the position using P, L, and/or U and
### compute the distance from the origin.
### Iterate time at each iteration of the loop.
### Save the time the loop breaks as A9.
### Save the distance from the origin as A10.

#jack this just checks if any of hte values in locationd are less than 0.06 and gives the entry number but none of them,
#are less than 0.06 so obviously i did something wrong earlier
for j in range(total):
    if locationd[j] < 0.06 :
        break

print(j)

### Problem 3
### Create a function here for the rotation matrix that
### takes an input in radians and returns the matrix.




    
    
### Save A11 as R(pi/8)
### Rotate the vector given in the assignment file and save it as A12.






### Find the vector x that was rotated to give you vector b.
### Save the vector x as A13






### Invert the R(3*pi/4) and save it as A14.
### Find the angle theta that would give you this inverse
### without having to do matrix operations, and save the angle
### as A15.



