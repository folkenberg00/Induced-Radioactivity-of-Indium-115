#!/usr/bin/env python3
'''
	26.09.2022
	Folkenberg Siro, OPTEX
	Nuclear Engineering and Thermal Physics, National Research Nuclear University - IATE, Russia
	Geospatial(Optical Remote Sensing), Technical University of Kenya 2019

	==================================== INDUCED RADIOACTIVITY OF INDIUM =========================================

	N/B: time(minutes)
	     Confidence level(0,95); you can choose your preferred value and compare the test criterion to the table value at the intersection
'''

from __future__ import absolute_import;
import matplotlib.pyplot;
import numpy;
from scipy import stats;
import pandas;
#from scipy.optimize import curve_fit;
#from sklearn.metrics import r2_score;

conf = 0.95; #confidence level

def fit(time, a, b):
                return a*numpy.exp(b*time);
excel_file = pandas.read_excel("radioactivity.data.xlsx"); #reads the xlsx file format
excel_file.to_csv ("radioactivity.data.csv", index = None, header=True); #converts xlsx to csv
data = pandas.read_csv("radioactivity.data.csv"); #retrieves data from csv file
time = data['TIME'].values; #required column
pulses = data['PULSES'].values;
n = len(pulses); #number of observations of pulses
measured_bg = data['MEASURED_BACKGROUND'];
measured_bg = measured_bg[numpy.isfinite]; #eliminate  NaN rows
y = numpy.log(pulses - (numpy.sum(measured_bg))/3.0); #the natural logarithm of number of pulses after background subtraction
p = 2; #number of parameters in the adopted equation(corresponds to y = mx + c)
nu = n - p; #number of degrees of freedom

'''least square method preprocessing variables'''
t_sum = numpy.sum(time); #sum of the elements in the time array
tsq_sum = numpy.sum(time**2);
y_sum = numpy.sum(y); #sum of the elements in the y array
ysq_sum = numpy.sum(y**2);
ytime_sum = numpy.sum(time*y);

''' decay constant, theoretical values comparable to the natural logarithm of the
measured number of pulses after background subtraction, Half-life '''
lambd = -((n*ytime_sum - t_sum*y_sum)/(n*tsq_sum - t_sum**2)); #decay constant
c = (tsq_sum*y_sum - t_sum*ytime_sum)/(n*tsq_sum - t_sum**2); # y intercept
y_theoret = c-lambd*time;
half_life = numpy.log(2)/lambd; #half-life

''' average quadratic errors in coefficients lambd and c '''
sigmasq_0 = ysq_sum/(n - 2) - y_sum**2/(n*(n - 2)) - (n*ytime_sum - t_sum*y_sum)**2/(n*(n - 2)*(n*tsq_sum - t_sum**2));
sigmasq_lambd = n*sigmasq_0/(n*tsq_sum - t_sum**2); #dispersion of coefficient lambda
sigmasq_c = sigmasq_0*tsq_sum/(n*tsq_sum - t_sum**2); #dispersion of coefficient c

''' absolute uncertainties'''
abs_lambd = 2*numpy.sqrt(sigmasq_lambd); #absolute uncertainty in decay constant
abs_c = 2*numpy.sqrt(sigmasq_c); #absolute uncertainty in c
abs_halflife = ((numpy.log(2))*abs_lambd)/lambd**2; #absolute uncertainty in half-life
abs_y = ((numpy.sqrt(pulses + (numpy.sum(measured_bg))/3.0))/(pulses - numpy.sum(measured_bg)/3.0));

''' the degree of reliability of parameters c and lambd '''
M_i = (y - y_theoret)**2/abs_y**2;
M = numpy.sum(M_i); #the degree of reliability of the obtained parameters, c and lambd

''' quick report on CLI '''
def quickcli_report():
	print ('DECAY CONSTANT:', lambd,'(',abs_lambd,')/min'); #decay constant
	print ('HALF-LIFE:', half_life,'(',abs_halflife,')min'); #half life
	print ('DEGREES OF FREEDOM:', nu ); #degrees of freedom
	print ('TEST CRITERION, M:', M); #test criterion
	print ('CONFIDENCE LEVEL:', conf); #confidence level
	print ('THEORETICAL Y:', y_theoret); #retrieved from the generated equationn of a straight line
	print ('UNCERTAINTY IN Y:', 2*abs_y); #absolute uncertainties in y
grad, intercept, r_val, p_val, std_err = stats.linregress(time, y);
mn = 0; #minimum value on abscissa extent
mx = numpy.max(time);
x1 = numpy.linspace(mn, mx, 500);
y1 = numpy.y1 = grad*x1+intercept;

fig, ax = matplotlib.pyplot.subplots(facecolor=("grey"), figsize=(8, 5));
ax.grid(True);
matplotlib.pyplot.rc('axes', labelsize=8);
matplotlib.pyplot.rc('axes', titlesize=11);
ax.set_facecolor('black');
ax.grid(which='major', linewidth=0.7);
ax.grid(which='minor', color='#EEEEEE', linestyle=':', linewidth=0.5);
ax.minorticks_on();
ax.text(0.99, 0.05, 'folkenberg plots', transform=ax.transAxes, fontsize=10, color='gray', alpha=0.7, ha='right', va='bottom', rotation='0');
sep_scale = ax.twinx();
def plot():
	sigma_y = 2*abs_y; #absolute uncertainty in  pulse and background measurements
	numpy.savetxt("./output/report.out.csv", numpy.dstack((numpy.arange(1, sigma_y.size+1),time**2, y**2, sigma_y, time*y,\
	y+sigma_y, y-sigma_y, y_theoret))[0],"%d, %.2f, %.7f, %.7f, %.7f, %.7f, %.7f, %.7f",header="ID, TIME_SQUARED, Y_SQUARED,\
	UNCERTAINTIES, TIME*Y, YFLAGS_UPPER, YFLAG_LOWER, THEORETICAL_Y"); #saves array to CSV
	plot_1, = ax.plot(time, y, 'og', linewidth=1);
	plot_2, = ax.plot(x1, y1, '--r', linewidth=1.5);
	ax.errorbar(time, y, yerr = sigma_y, fmt="o", color = 'g', linewidth=1.2); #error bars
	ax.set_xlabel("Time (minutes)");
	ax.set_ylabel('ln | N$_{i}$ - N$_{bg}$ |');
	ax.set_title("Neutron Activation of Indium - 115 (Induced Radioactivity)");
	ax.text(0.90, 0.68, 'ln | N$_{i}$ - N$_{bg}$ | = %.4f - %.3f t'%(c, lambd), transform=ax.transAxes, fontsize=9, fontweight="bold",\
		color='b', ha='right', va='bottom', rotation='0', bbox=dict(facecolor='k'));
	ax.text(0.24, 0.40, 'DECAY CONSTANT: ', transform=ax.transAxes, fontsize=8, color='r', ha='right', va='bottom', rotation='0',\ 
		bbox=dict(facecolor='k'));
	ax.text(0.44, 0.40, '%.4f ( %.4f)/min'%(lambd, abs_lambd), transform=ax.transAxes, fontsize=8, color='g', ha='right', va='bottom',\
		rotation='0', bbox=dict(facecolor='k'));
	ax.text(0.169, 0.35, 'HALF-LIFE: ', transform=ax.transAxes, fontsize=8, color='r', ha='right', va='bottom', rotation='0',\
		bbox=dict(facecolor='k'));
	ax.text(0.32, 0.35, '%.2f ( %.2f)min'%(half_life, abs_halflife), transform=ax.transAxes, fontsize=8, color='g', ha='right',\
		va='bottom', rotation='0', bbox=dict(facecolor='k'));
	ax.text(0.29, 0.30, 'DEGREES OF FREEDOM: ', transform=ax.transAxes, fontsize=8, color='r', ha='right', va='bottom', rotation='0',\
		bbox=dict(facecolor='k'));
	ax.text(0.32, 0.30, '%d'%(nu), transform=ax.transAxes, fontsize=8, color='g', ha='right', va='bottom', rotation='0',\
		bbox=dict(facecolor='k'));
	ax.text(0.229, 0.25, 'TEST CRITERION: ', transform=ax.transAxes, fontsize=8, color='r', ha='right', va='bottom', rotation='0',\
		bbox=dict(facecolor='k'));
	ax.text(0.30, 0.25, '%.2f'%(M), transform=ax.transAxes, fontsize=8, color='g', ha='right', va='bottom', rotation='0',\
		bbox=dict(facecolor='k'));
	ax.text(0.26, 0.20, 'CONFIDENCE LEVEL: ', transform=ax.transAxes, fontsize=8, color='r', ha='right', va='bottom',\
		rotation='0', bbox=dict(facecolor='k'));
	ax.text(0.31, 0.20, '%.2f'%(conf), transform=ax.transAxes, fontsize=8, color='g', ha='right', va='bottom', rotation='0',\
		bbox=dict(facecolor='k'));
	ax.legend([plot_1, plot_2], ["Data Points", 'Linear Fit, R = {:.3f}'.format(r_val)], borderpad=2, prop = {'size' : 6.5},\
		  loc = 'upper right', shadow = True, facecolor="w", title = "Legend");
	matplotlib.pyplot.savefig('./output/plot.png');
	matplotlib.pyplot.show();
	exit();
if __name__ == "__main__":
	plot();
