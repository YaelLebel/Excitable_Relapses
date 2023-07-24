import numpy as np

def weeks_to_trimesters (values_weeks,weeks):
    rate_trimester_to_week = 13
    values_trimesters = np.zeros(int(len(values_weeks)/rate_trimester_to_week))
    trimesters = []
    for i in range(int(len(values_weeks)/rate_trimester_to_week)):
        for j in range(rate_trimester_to_week):
            ind = i*rate_trimester_to_week+j
            values_trimesters[i] += values_weeks[ind]
        values_trimesters[i] = values_trimesters[i]/(j+1)
        trimesters.append(i)
    trimesters = np.array(trimesters)-int(abs(weeks[0])/rate_trimester_to_week)
    return trimesters,values_trimesters

def det_values_for_fit (trimesters_lymp,values_trimesters_lymp,
                        trimesters_rr,values_trimesters_rr):
    lymp_values_for_fit = []
    freq_relapses_for_fit = []
    for v,t in zip(values_trimesters_rr,trimesters_rr):
        print(v,t)
        find = [abs(tr-t) for tr in trimesters_lymp]
        if min(find)<1:
            ind = np.argmin(find)
            lymp_values_for_fit.append(values_trimesters_lymp[ind])
            freq_relapses_for_fit.append(v)
    print(lymp_values_for_fit)
    lymp_baseline = lymp_values_for_fit[-2]
    rr_baseline = freq_relapses_for_fit[-2]

    return lymp_values_for_fit,freq_relapses_for_fit,lymp_baseline, rr_baseline

def fun(x,q):
    return np.exp(x*q)

from scipy.optimize import curve_fit
        
def create_fit (lymp_values_for_fit,lymp_baseline,freq_relapses_for_fit,rr_baseline):
    alphas_for_fit = [v/lymp_baseline-1 for v in lymp_values_for_fit]
    popt,pcov = curve_fit(fun,alphas_for_fit,np.array(freq_relapses_for_fit)/rr_baseline)
    return popt,pcov,alphas_for_fit


def compute_chi2_correlation_coeff (popt,alphas_for_fit,freq_relapses_for_fit,rr_baseline):
    from scipy import stats
    slope, intercept, r_value, p_value, std_err = stats.linregress([fun(a,*popt) for a in alphas_for_fit],np.array(freq_relapses_for_fit)/rr_baseline)
    chi_square_test_statistic1 = 0
    expected_data = [fun(a,*popt) for a in alphas_for_fit]
    observed_data = np.array(freq_relapses_for_fit)/rr_baseline
    for i in range(len(observed_data)):
        chi_square_test_statistic1 = chi_square_test_statistic1 + \
            (np.square(observed_data[i]-expected_data[i]))/expected_data[i]
            
    print('chi square value determined by formula : ' +
          str(chi_square_test_statistic1))
      
    # find Chi-Square critical value
    pval_chi2 = stats.chi2.cdf(chi_square_test_statistic1, len(expected_data)-1)
    print('p value : ',pval_chi2)
    red_chi2 = chi_square_test_statistic1/(len(expected_data)-1)
    print('reduced chi square : ',red_chi2)
    from scipy.stats import pearsonr

    # X and Y are the two variables to be correlated
    corr, pval = pearsonr(expected_data, observed_data)

    print("Correlation coefficient:", corr)
    print("p-value:", pval)
    
    return chi_square_test_statistic1, pval_chi2,red_chi2,corr,pval

import matplotlib.pyplot as plt

def create_figure_rr_from_data (weeks,rr_weeks):
    fig,ax = plt.subplots()
    ax.plot(weeks,rr_weeks,color = '#2c7bb6',label = 'data frpm paper')
    ax.fill_between(weeks,rr_weeks+0.15,rr_weeks-0.15, color = '#2c7bb6',alpha = 0.2)
    ax.set_xlim(min(weeks),max(weeks))
    ax.spines[['top','right']].set_visible(False)
    ax.set_xlabel('Weeks from delivery')
    ax.set_ylabel('Relapse Rate [1/yr]')
    ax.set_box_aspect(1)
    ax.set_xlim(min(weeks),max(weeks))
    return fig,ax

def create_figure_lymp_count (weeks,lymp_weeks):
    fig, ax = plt.subplots()
    ax.set_xlabel('Weeks from delivery')
    ax.set_ylabel('Lymphocyte count')
    ax.plot(weeks,lymp_weeks, color = '#fdae61')
    ax.spines[['top','right']].set_visible(False)
    ax.set_box_aspect(1)
    ax.set_xlim(weeks[0],weeks[-1])
    return fig,ax

def create_figure_rr_from_fitted_model (weeks,lymp_baseline,lymp_weeks,popt,rr_baseline,fig = None,ax = None):
    alphas_weeks = [l/lymp_baseline-1 for l in lymp_weeks]
    if (fig is None) and (ax is None):
        fig,ax = plt.subplots()
    ax.set_xlabel('Weeks from delivery')
    ax.set_ylabel('Relapse Rate [1/yr]')
    ax.plot(weeks,[fun(a,*popt)*rr_baseline for a in alphas_weeks],color = 'black', label = 'computed days to relapse \n from lymphocyte rate (Clalit db)')
    ax.spines[['top','right']].set_visible(False)
    ax.set_box_aspect(1)
    ax.set_xlim(min(weeks),max(weeks))
    return fig,ax

def create_figure_overlay (fig_data,ax_data,weeks,lymp_baseline,lymp_weeks,popt,rr_baseline):
    fig,ax = create_figure_rr_from_fitted_model(weeks, lymp_baseline, lymp_weeks, popt, rr_baseline,fig_data,ax_data)
    return fig, ax
import pandas as pd

if __name__=="__main__":
    filename1 = 'relapse_freq_from_delivery_trimesters.csv'
    filename2 = 'LYMP.abs.csv'
    relapse_data = pd.read_csv(filename1, header = None)
    trimesters_rr = [round(m)+1 for m in relapse_data[0]]
    values_trimesters_rr = relapse_data[1]
    lymp = pd.read_csv(filename2)
    lymp_weeks = lymp['mean']
    weeks = np.arange(-59,81)
    trimesters_lymp, values_trimesters_lymp = weeks_to_trimesters(lymp['mean'], weeks)
    
    lymp_values_for_fit,freq_relapses_for_fit,lymp_baseline, rr_baseline = det_values_for_fit(trimesters_lymp, values_trimesters_lymp, trimesters_rr, values_trimesters_rr)
    popt,pcov,alphas_for_fit = create_fit(lymp_values_for_fit,lymp_baseline,freq_relapses_for_fit,rr_baseline)
    chi_square_test_statistic1, pval_chi2,red_chi2,corr,pval_corr = compute_chi2_correlation_coeff(popt, alphas_for_fit, freq_relapses_for_fit, rr_baseline)
    
    filename3 = 'relapse_freq_from_delivery.csv'
    relapse_data_weeks = pd.read_csv(filename3,header=None)
    weeks_ = relapse_data_weeks[0]
    weeks_rr_from_data = relapse_data_weeks[1]

    create_figure_lymp_count(weeks, lymp_weeks)
    fig_data,ax_data = create_figure_rr_from_data(weeks_, weeks_rr_from_data)
    
    create_figure_rr_from_fitted_model(weeks, lymp_baseline, lymp_weeks, popt, rr_baseline)
    
    create_figure_overlay(fig_data, ax_data, weeks, lymp_baseline, lymp_weeks, popt, rr_baseline)