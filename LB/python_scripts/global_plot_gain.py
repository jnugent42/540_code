import matplotlib.pyplot as plt
import numpy as np

#always give data in centre, left, right, top, bottom  
Square_mean = [129.1, 129.1, 128.9, 128.6, 128.6]
#note: between right and centre: we have cut the fibre
Square_gain = [0.4057, 0.4094, 0.4134, 0.405,0.404 ]
Square_gain_err = [0.0007, 0.0019, 0.0014, 0.0010, 0.0010]

Circle_mean = [118.9, 118.9, 118.9, 119.00, 118.9]
Circle_gain = [0.4187, 0.4217, 0.4227, 0.4196, 0.4135]
Circle_gain_err = [0.0011, 0.0011, 0.0012, 0.0012, 0.0012]

Losange_mean = [124.4, 124.4, 124.4, 124.4, 124.4]
Losange_gain = [0.4176, 0.4117, 0.4170, 0.4179, 0.4037]
Losange_gain_err = [0.0015, 0.0011, 0.0010, 0.0011, 0.0010]

Triangle_mean = [131.8, 131.7, 131.7, 131.8, 131.7]
Triangle_gain = [0.4218, 0.4179, 0.4166, 0.4053, 0.4149]
Triangle_gain_err = [0.0010, 0.0013, 0.0014, 0.0012, 0.0011]


#note: between right and centre: we have cut the fibre
Square_gain = [1.635e6, 1.848e6, 1.587e6, 1.616e6, 1.627e6 ]
Square_gain_err = [1.342e3, 1.018e5, 2.410e3, 2.077e3, 7.398e3]

Circle_gain = [9.461e5, 8.574e5, 8.730e5, 9.003e5, 9.348e5]
Circle_gain_err = [3.305e4, 1.876e4, 2.528e4, 2.454e4, 2.578e4]

Losange_gain = [9.188e5, 9.251e5, 1.092e6, 1.118e6, 9.406e5]
Losange_gain_err = [2.467e3, 2.356e3, 1.622e4, 1.809e4, 2.161e3]

Triangle_gain = [1.242e6, 1.213e6, 1.100e6, 1.354e6, 1.078e6]
Triangle_gain_err = [1.342e3, 2.119e3, 2.410e3, 2.077e3, 7.398e3]






position = ["centre", "left", "right", "top", "bottom"]


#here the plotting
plt.subplot(2,1, 1)
plt.errorbar(position[:len(Square_gain)], Square_gain, yerr = Square_gain_err, marker = 's', label = 'Square PMT', color = 'blue')
plt.plot([position[0], position[len(Square_gain)-1]],[np.array(Square_gain).mean(), np.array(Square_gain).mean()], '--', color = 'blue')

plt.errorbar(position[:len(Circle_gain)], Circle_gain, yerr = Circle_gain_err, marker = 'o', label = "Circle PMT", color = 'red')
plt.plot([position[0], position[len(Circle_gain)-1]],[np.array(Circle_gain).mean(), np.array(Circle_gain).mean()], '--', color = 'red')

plt.errorbar(position[:len(Losange_gain)], Losange_gain, yerr = Losange_gain_err, marker = 'D', label = "Losange PMT", color = 'green')
plt.plot([position[0], position[len(Losange_gain)-1]],[np.array(Losange_gain).mean(), np.array(Losange_gain).mean()], '--', color = 'green')

plt.errorbar(position[:len(Triangle_gain)], Triangle_gain, yerr = Triangle_gain_err, marker = '^', label = 'Triangle PMT', color = 'black')
plt.plot([position[0], position[len(Triangle_gain)-1]],[np.array(Triangle_gain).mean(), np.array(Triangle_gain).mean()], '--', color = 'black')



##Here the printing of the values 
Square_gain = np.array(Square_gain)
largest_pos_deviation_Square = (Square_gain - Square_gain.mean()).max()/Square_gain.mean()
largest_neg_deviation_Square = (Square_gain - Square_gain.mean()).min()/Square_gain.mean()

if abs(largest_pos_deviation_Square) > abs(largest_neg_deviation_Square):
    largest_deviation_Square = 100 * largest_pos_deviation_Square
else:
    largest_deviation_Square = 100 * largest_neg_deviation_Square

print("\nMean Gain accros the %i positions for the Square PMT: %.3e, largest deviation from the mean %.2f percent"%(len(Square_gain), Square_gain.mean(), largest_deviation_Square))

Circle_gain = np.array(Circle_gain)
largest_pos_deviation_Circle = (Circle_gain - Circle_gain.mean()).max()/Circle_gain.mean()
largest_neg_deviation_Circle = (Circle_gain - Circle_gain.mean()).min()/Circle_gain.mean()

if abs(largest_pos_deviation_Circle) > abs(largest_neg_deviation_Circle):
    largest_deviation_Circle = 100 * largest_pos_deviation_Circle
else:
    largest_deviation_Circle = 100 * largest_neg_deviation_Circle

print("Mean Gain accros the %i positions for the Circle PMT: %.3e, largest deviation from the mean %.2f percent"%(len(Circle_gain), Circle_gain.mean(), largest_deviation_Circle))

Losange_gain = np.array(Losange_gain)
largest_pos_deviation_Losange = (Losange_gain - Losange_gain.mean()).max()/Losange_gain.mean()
largest_neg_deviation_Losange = (Losange_gain - Losange_gain.mean()).min()/Losange_gain.mean()

if abs(largest_pos_deviation_Losange) > abs(largest_neg_deviation_Losange):
    largest_deviation_Losange = 100 * largest_pos_deviation_Losange
else:
    largest_deviation_Losange = 100 * largest_neg_deviation_Losange

print("Mean Gain accros the %i positions for the Losange PMT: %.3e, largest deviation from the mean %.2f percent"%(len(Losange_gain), Losange_gain.mean(), largest_deviation_Losange))

Triangle_gain = np.array(Triangle_gain)
largest_pos_deviation_Triangle = (Triangle_gain - Triangle_gain.mean()).max()/Triangle_gain.mean()
largest_neg_deviation_Triangle = (Triangle_gain - Triangle_gain.mean()).min()/Triangle_gain.mean()

if abs(largest_pos_deviation_Triangle) > abs(largest_neg_deviation_Triangle):
    largest_deviation_Triangle = 100 * largest_pos_deviation_Triangle
else:
    largest_deviation_Triangle = 100 * largest_neg_deviation_Triangle

print("Mean Gain accros the %i positions for the Triangle PMT: %.3e, largest deviation from the mean %.2f percent\n"%(len(Triangle_gain),Triangle_gain.mean(), largest_deviation_Triangle))



#plt.subplot(2,1,1)
plt.title("Gain of the four PMTs", fontsize=15, weight = 'bold')
plt.ylabel("Gain", fontsize=15)
plt.legend()
plt.grid()


#here - plot with respect to the mean
plt.subplot(2,1, 2)
#Square_gain =np.sort(Square_gain)
#Circle_gain =np.sort(Circle_gain)
#Losange_gain =np.sort(Losange_gain)
#Triangle_gain =np.sort(Triangle_gain)
#position = [0, 1,2, 3, 4]

plt.errorbar(position[:len(Square_gain)], Square_gain/np.array(Square_gain).mean(), yerr = Square_gain_err/np.array(Square_gain).mean(), marker = 's', color = 'Blue')
plt.errorbar(position[:len(Circle_gain)], Circle_gain/np.array(Circle_gain).mean(), yerr = Circle_gain_err/np.array(Circle_gain).mean(), marker = 'o', color = 'Red')
plt.errorbar(position[:len(Losange_gain)], Losange_gain/np.array(Losange_gain).mean(), yerr = Losange_gain_err/np.array(Losange_gain).mean(), marker = 'D', color = 'Green')
plt.errorbar(position[:len(Triangle_gain)], Triangle_gain/np.array(Triangle_gain).mean(), yerr = Triangle_gain_err/np.array(Triangle_gain).mean(), marker = '^', color = 'black')
plt.plot([position[0], position[len(Triangle_gain)-1]],[1, 1], '--', label = 'Mean PMT Gain', color = 'gray')
plt.ylabel("Gain/Mean PMT Gain", fontsize=15)
plt.xlabel("Position on the PMT face", fontsize=15)
plt.legend()
plt.grid()

plt.show()
