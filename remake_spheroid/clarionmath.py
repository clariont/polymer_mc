import math

#This function takes the derivative of a function, whose y-values are passed in as a list.
#The list must have equally-spaced x-values.  Uses the forward difference method.
#Note: The returned list will be shorter by one element.
def fwdDifference(yvalues, spacing):
	derivValues = []
	counter = 0
	while (counter < (len(yvalues)-1)):
		deriv = (yvalues[counter+1]-yvalues[counter])/float(spacing)
		derivValues.append(deriv)
		counter = counter+1
	
	return derivValues
	

def fwdDifference1(yvalues, xvalues):
	derivValues = []
	counter = 0
	while (counter < (len(yvalues)-1)):
		deriv = (yvalues[counter+1]-yvalues[counter])/float(xvalues[counter+1]-xvalues[counter])
		derivValues.append(deriv)
		counter = counter+1
	
	return derivValues
	
#xvalues contains the x values of the function to be integrated.
#yvalues contains the y values of the function to be integrated.
#lo is the lower bound of the integral.
#hi is the upper bound of the integral.
def trapIntegrate(xvalues, yvalues, lo, hi):
	area = 0
	counter = xvalues.index(lo)
	end = xvalues.index(hi)
	while (counter < end):
		area = area + (xvalues[counter+1]-xvalues[counter])*(yvalues[counter+1] + yvalues[counter])/2
		counter = counter + 1
	return area


#Returns the index of the greatest element in the list.
def greatestIndex(myList):
	biggest = myList[0]
	bigIndex = 0
	counter = 0
	for element in myList:
		if (element > biggest):
			biggest = element
			bigIndex = counter
		counter = counter + 1
	return bigIndex


# Calculates the distance between two points in three dimensions.
def distance(x1,y1,z1, x2,y2,z2):
    dsq = (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)
    return math.sqrt(dsq)
