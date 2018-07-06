import soundfile as sf
import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
%matplotlib inline

import gammatone.filters as gt_filters
import scipy.signal as dsp

def Hz2Mel(freq):
    # the convert equation is not the regular one
    return 1000*np.log10(1+freq/800.0)
def Mel2Hz(mel):
    return (10**(mel/1000)-1)*800

def MelFreqBand(N,low,high):
    # get cutoff frequency of filterbanks whith center frequency equally distribute in MEL scale
    low_mel = Hz2Mel(low)
    high_mel = Hz2Mel(high)
    # the upper cutoff frequency of band i equals the lower cutoff frequency of band i+1
    bandwidth_mel = (high_mel-low_mel)/(N+1)
    
    cutoff_mel = np.arange(bandwidth_mel,high_mel+bandwidth_mel,bandwidth_mel)
    
    cutoff = np.zeros((N,2))
    cutoff[:,0] = Mel2Hz(cutoff_mel[0:N])
    cutoff[:,1] = Mel2Hz(cutoff_mel[1:N+1])
    cf = np.mean(cutoff,axis=1)
    
    return [cf,cutoff]
                           

def MakeFilterBank(fs,N):
    """buterworth filter bank"""
    # fs: sample frequency
    # N: frequency channel number
    
    cf,cutoff = MelFreqBand(N=N,low=0,high=fs/2)
    # filter coefficience
    B = []
    A = []

    is_plot = False
    if is_plot:
        plt.figure()
        
    for n in xrange(N-1):
        b,a = dsp.butter(N=3,Wn=cutoff[n]/fs*2,btype='bandpass')
        B.append(b)
        A.append(a)
        if is_plot:
            w,h = dsp.freqz(a=a,b=b)
            plt.plot(w,10*np.log10(np.abs(h)+1e-10))
    # the last filter cover the rest frequency range
    b,a = dsp.butter(N=6,Wn=cutoff[n,0]/fs*2,btype='highpass')
    B.append(b)
    A.append(a)
    
    return [B,A]

def FB_filter(x,B,A):
    if len(x.shape)>1 and x.shape[2]>1:
        raise Exception('Only input with one channel is supported in FB filter')
    x_len = x.shape[0]
    
    if len(B)!= len(A):
        raise Exception('filter coefficiences do not match in FB filter')
    channel_num = len(B)
    
    x_band = np.zeros((x_len,channel_num))
    for n in xrange(channel_num):
        x_band[:,n] = dsp.lfilter(b=B[n],a=A[n],x=x) 
        
    return x_band

def GetEvelope(x,fs):
    # simple implement
    # half-wave rectification/suquare + decimate 

    ratio = fs/4000 # equivalent to resample to 4kHz

    # half-wave rectification
    envelope_opt = 'abs'
    if envelope_opt == 'abs':
        x = np.abs(x)
    elif envelope_opt == 'square':
        x = x**2

    return dsp.decimate(x,ratio,axis=0,ftype='fir')

def GetFreqConvertMatrix(start_bin,end_bin,fft_len,channel_num):
    # Triangular window, centers are normally distributed in the given range
    # window length is twice of the distance of two adjacent centers
    step = np.float64((end_bin-start_bin))/(channel_num+1)
    cf_bin = step*np.arange(1,channel_num+1)# center frequency in bins
    matrix = np.zeros((fft_len,channel_num))
    for n in xrange(channel_num):
        bin_start = np.int16(np.ceil(cf_bin[n]-step))
        bin_end = np.int16(np.floor(cf_bin[n]))
        for freq_bin in xrange(bin_start,bin_end+1):
            matrix[freq_bin,n] = 1-(cf_bin[n]-freq_bin)/step
        
        bin_start = np.int16(np.ceil(cf_bin[n]))
        bin_end = np.int16(np.floor(cf_bin[n]+step))
        for freq_bin in xrange(bin_start,bin_end+1):
            matrix[freq_bin,n] = 1-(freq_bin-cf_bin[n])/step
    return matrix
    
def GetAMS(wav_path):
    
    # filter bank
    wav,fs = sf.read(wav_path)

    # frequency decomposition
    N = 22
    B,A = MakeFilterBank(fs=fs,N=N)# filter bank
    
    wav_band = FB_filter(x=wav,B=B,A=A)
    
    # get envelope
    env_band = GetEvelope(wav_band,fs)
    env_len = env_band.shape[0]

    # calculation of AMS feature
    feature_len = 15
    frame_len = 128
    overlap = 64
    frame_num = np.int64(np.floor((env_len-frame_len)/overlap)+1)
    
    AMS_matrix = GetFreqConvertMatrix(start_bin=0,end_bin=26,fft_len=frame_len,channel_num=feature_len)
    
    AMS = np.zeros((N,frame_num,feature_len))
    win = np.hanning(frame_len)

    for i in xrange(frame_num):
        for n in xrange(N):
            frame = env_band[i*overlap:i*overlap+frame_len,n]
            fft_amp = np.abs(np.fft.fft(np.multiply(frame,win),n=2*frame_len))[:frame_len]

            AMS[n,i,:] = np.dot(AMS_matrix.T,fft_amp)

    return AMS