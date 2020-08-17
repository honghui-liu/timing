from stingray import Lightcurve, Powerspectrum, AveragedPowerspectrum
import subprocess as sp
import numpy as np
import os


def getgti(infile):
    """
    Get the GTIs from the light curve fits file
    
    Parameter:
    ------
    infile: str
        The name of the input light curve
    """
    with fits.open(infile) as hdul:
 
        gti_row = hdul[2].data
        gti = np.empty(shape=(len(gti_row['START']),2))
        for i in range(len(gti_row['START'])):
            gti[i,0] = gti_row['START'][i]
            gti[i,1] = gti_row['STOP'][i]
 
        time = hdul[1].data['TIME']
        interval = np.array([time[0]])
        for k in range(len(time)-1):
            if (time[k+1]-time[k]>50):
                interval = np.append(interval, [time[k]], axis=0)
                interval = np.append(interval, [time[k+1]], axis=0)
 
        interval = np.append(interval, [time[len(time)-1]], axis=0)
        interval = interval.reshape((-1, 2)) 
 
        if (interval.shape[0]<gti.shape[0]):
            gti_corr = interval
        else:
            gti_corr = gti 
   
    return gti_corr


def psd2xsp(fname, outname, direct_save=False, freq=None, noise=None):
    """
    Save the averaged PSD from Stingray in a format readable by XSPEC. 
    This function is mainly a reproduce of the 'save_to_xspec' module of HENDRICS:
    https://hendrics.readthedocs.io/en/latest/_modules/hendrics/save_as_xspec.html#save_as_xspec
    
    Parameters
    ------------
    fname: str
        Input Stingray psd file name
    outname: str
        Output file name without extension (e.g. maxij1535)
    cwd: str
        Current working directory
    freq: float
        The frequency above which can be considerred as pure white noise
    noise: float
        The level of white noise
    direct_save: bool
        If True: call flx2xsp to produce the output .pha and .rsp files.
        If False: flx2xsp has to be called from the user
    ------------
    The 'flx2xsp' tool:
    https://casdc.china-vo.org/mirror/AstroSoft/HEASoft/lheasoft6.10/source/ftools/heasarc/src/flx2xsp/flx2xsp.txt
    """
    
    # Set the environment to avoid error like this:
    # 'Unable to redirect prompts to the /dev/tty (at headas_stdio.c:...)'
    # see https://heasarc.gsfc.nasa.gov/lheasoft/scripting.html
    os.environ['HEADASNOQUERY'] = ''
    os.environ['HEADASPROMPT'] = '/dev/null'

    # calculate the level of white noise by averaging the power in high frequency range.
    if freq is None:
        if noise is None:
            white_noise = 0.0
            print('psd2xsp: no white noise is considerred!')
        else:
            white_noise = noise
    else:
        if noise is None:
            white_noise = np.mean(fname.power[np.where(fname.freq>freq)])
        else:
            white_noise = np.mean(fname.power[np.where(fname.freq>freq)])
            print('Both \'freq\' and \'noise\' are defined! \n Will consider frequency as the criteria!')

    flo = fname.freq - fname.df/2
    fhi = fname.freq + fname.df/2
    power = (fname.power - white_noise) * fname.df # white noise subtracted
    power_err = fname.power_err * fname.df
    
    np.savetxt(outname+'.txt', np.transpose([flo, fhi, power, power_err]))
    
    if direct_save:
        sp.check_call('flx2xsp {0} {1}.pha {1}.rsp'.format(outname+'.txt', outname).split())
    
    return True, white_noise

#sp.check_call('flx2xsp infile=temp.txt phafil=temp1.pha rspfil=temp.rsp'.split(), cwd='/Users/honghui/projects/hxmt/timing')