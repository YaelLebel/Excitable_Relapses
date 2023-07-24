import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_default_settings import load_def_settings


def rfp_hazard (arr,bins):
    #compute relapse free probability and hazard curves from data
    count, bins_count = np.histogram(arr, bins=bins)
    pdf = count / sum(count)
    cdf = np.cumsum(pdf)
    rfp = 1- cdf
    log_rfp = np.log(rfp)
    dt = bins_count[1]-bins_count[0]
    hazard = [(log_rfp[i]-log_rfp[i+1])/dt for i in range(len(log_rfp)-2)]
    return rfp,hazard,bins_count

#bootstrapping
def bootstrapping (num_bins,diffs_all,num_samples,size_samples):
    bins = np.linspace(min(diffs_all),max(diffs_all),num_bins+1)
    samples_s = np.zeros(shape = (num_samples,num_bins)) #survival
    samples_h = np.zeros(shape = (num_samples,num_bins-2)) #hazard
    for i in range(num_samples):
        y = np.random.choice(diffs_all,size_samples,replace = True)
        s,h,b = rfp_hazard(y,bins)
        samples_s[i] = s
        samples_h[i] = h
    stds_s = np.std(samples_s,axis=0)
    mean_s = np.mean(samples_s,axis = 0)
    s_top = [m+s for m,s in zip(mean_s,stds_s)]
    s_bottom = [m-s for m,s in zip(mean_s,stds_s)]
    stds_h = np.ma.masked_invalid(samples_h).std(axis=0)
    mean_h = np.ma.masked_invalid(samples_h).mean(axis=0)
    h_top = [m+s for m,s in zip(mean_h+0.5,stds_h)]
    h_bottom = [m-s for m,s in zip(mean_h+0.5,stds_h)]
    return bins,stds_s,mean_s,s_top,s_bottom,stds_h,mean_h,h_top,h_bottom


def create_rfp_hazard_curves (filename_diffs,bs_num_bins,bs_num_samples,bs_size_samples):
    #create relapse free probability and hazard curves
    df_diffs = pd.read_csv(filename_diffs)
    diffs_all = df_diffs.DAYSDIFF.values

    bins, stds_s,mean_s,s_top,s_bottom,stds_h,mean_h,h_top,h_bottom = bootstrapping(bs_num_bins, diffs_all, bs_num_samples, bs_size_samples)
    load_def_settings()
    color_s = '#fc8d59'
    color_h = '#91bfdb'
    fig,ax = plt.subplots()
    ax.fill_between(bins[1:],s_top,s_bottom,alpha = 0.4,color = color_s)
    ax.plot(bins[1:], mean_s, color_s,label="Relapse-free \n probability")

    ax.fill_between(bins[2:-1],h_top,h_bottom,alpha = 0.4,color = color_h)
    ax.plot(bins[2:-1], mean_h+0.5, color = color_h,label="Relapse hazard")
    ax.spines[['top','right']].set_visible(False)
    ax.set_xlabel('Time [Days]')
    ax.set_ylabel('Probability')
    ax.legend()
    
if __name__ == "__main__":    
    filename_diffs = 'Data_Diffs_only_exact_dates.csv'
    #bootstrapping parameters
    bs_num_bins = 50
    bs_num_samples = 1000
    bs_size_samples = 200
    create_rfp_hazard_curves(filename_diffs, bs_num_bins, bs_num_samples, bs_size_samples)