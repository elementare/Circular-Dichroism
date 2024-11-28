import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
from os import listdir
from os.path import isfile, join
from matplotlib import font_manager
import os
import seaborn as sns
from LumiCD.spectrum import Spectrum
from matplotlib.colors import LinearSegmentedColormap

font_path = 'Fontes/Ruda/Ruda-VariableFont_wght.ttf'
font_manager.fontManager.addfont(font_path)
prop = font_manager.FontProperties(fname=font_path)

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = prop.get_name()

def smooth(y):
    
    window_size = 11  # Tamanho da janela (deve ser um número ímpar)
    poly_order = 3    # Ordem do polinômio
    y_smooth = savgol_filter(y, window_size, poly_order)
    
    return y_smooth

def getFiles(data):
    for caminho in data:
        path = f"Dados CD/{caminho}/"
        desiredMetada = data[caminho]
        
        for name in desiredMetada:
            
            for sample in desiredMetada[name]["samples"]:
                idSample = sample["id"]
                
                # print(data[caminho][name]["samples"])
                sample["filesTXT"] = [join(path, f) for f in listdir(path) if isfile(join(path, f)) and f.startswith(f"{idSample}_") and os.path.splitext(f)[1] == ".txt"]
                sample["filesGEN"] = [join(path, f) for f in listdir(path) if isfile(join(path, f)) and f.startswith(f"{idSample}_") and os.path.splitext(f)[1] == ".gen"]
            
        
    return data



def proCDAll(metaData, fileType = "filesTXT", buffer = False):
    metaDataProc = {}
    for sampleIdentity in metaData:
        if sampleIdentity != buffer:
            metaDataProc[sampleIdentity] = procCD(metaData, sampleIdentity, fileType=fileType, buffer=buffer)
        
    return metaDataProc
       

def getAbs(abs_actual):    
    avg = abs_actual.mean()
    std = abs_actual.std()
    normalized = (abs_actual - avg) / std
    standard_error = std / np.sqrt(len(abs_actual))

    abs_actual = normalized.copy()
    abs_actual = smooth(abs_actual)
    
    return abs_actual, standard_error

def getAbsBuffer(buffer):
    files = buffer["filesTXT"]

    for i in range(1, len(files)+1):
        
        df = pd.read_csv(files[i-1],names=['WL','CD Abs'],sep=' ')

        if i == 1: cd_abs_actual = df['CD Abs']
        else: cd_abs_actual = np.add(np.array(df['CD Abs']),cd_abs_actual)
    
    standard_error = 0
    # cd_abs_actual, standard_error = getAbs(cd_abs_actual)
    
    cd_abs_actual = smooth(cd_abs_actual)

    return np.array(df['WL']), cd_abs_actual, standard_error


def procCD(metaData, sampleName, fileType="filesTXT", buffer=False):
    
    dataFrames = {}
    
    desiredSamples = metaData[sampleName]["samples"]

    for sample in desiredSamples:
        idSample = sample["id"]
        for i in range(1, len(sample[fileType])+1):
            df = pd.read_csv(sample[fileType][i-1] ,names=['WL','CD Abs'],sep=' ')

            if i == 1: cd_abs_actual = df['CD Abs']
            else: cd_abs_actual = np.add(np.array(df['CD Abs']),cd_abs_actual) 
            
        cd_abs_actual = [x/6 for x in cd_abs_actual]
        cd_abs_actual = np.array(cd_abs_actual)

        if buffer:
            bufferData = metaData[buffer]
            bufferDesired = None
            for sampleBuffer in bufferData["samples"]:
                
                if sampleBuffer["distance"] == sample["distance"] and sampleBuffer["PMT"] == sample["PMT"] and sampleBuffer["size"] == sample["size"] and sampleBuffer["id"] != sample["id"]:
                    # print(sampleBuffer, sample)
                    print(sampleBuffer["id"], sample["id"])
                    bufferDesired = sampleBuffer
                    break
            else:
                raise FileNotFoundError(f"I wasn't able to find a blank for {sampleName} at distance {sample['distance']} and PMT {sample['PMT']}")
            
            wl_buffer, abs_buffer, se_buffer = getAbsBuffer(bufferDesired)

            if len(list((df['WL']))) > len(wl_buffer):
                cd_abs_actual = cd_abs_actual[list(df['WL']).index(260):list(df['WL']).index(180)]

            cd_abs_actual = np.subtract(cd_abs_actual,abs_buffer)

        # cd_abs_actual, standard_error = getAbs(cd_abs_actual)
        standard_error = 0
        
        cd_abs_actual = smooth(cd_abs_actual)
        size = sample["size"] if "size" in  sample.keys() else None
        distance = sample["distance"] if "distance" in  sample.keys() else None
        pmt = sample["PMT"] if "PMT" in sample.keys() else None
        if buffer:

            dataFrames[idSample] = {"WL": list(range(180,261))[::-1], "cd_abs_actual": cd_abs_actual, "standard_error": standard_error, "size": size, "distance": distance, "pmt": pmt}
        else:
            dataFrames[idSample] = {"WL": np.array(df['WL']), "cd_abs_actual": cd_abs_actual, "standard_error": standard_error, "size": size, "distance": distance, "pmt": pmt}
    return dataFrames
    

def cdproc(index,caminho,water=False):
    
    for i in range(1,7):
        
        df = pd.read_csv(f'{caminho}/{index}_converted_{i}.txt',names=['WL','CD Abs'],sep=' ')

        if i == 1: cd_abs_actual = df['CD Abs']
        else: cd_abs_actual = np.add(np.array(df['CD Abs']),cd_abs_actual)
    
    cd_abs_actual = [x/6 for x in cd_abs_actual]
    cd_abs_actual = np.array(cd_abs_actual)

    if water:

        wl_water, abs_water, se_water = cdproc(3,'CD_data/water_data')

        if len(list((df['WL']))) > len(wl_water):
            cd_abs_actual = cd_abs_actual[list(df['WL']).index(251):list(df['WL']).index(190)]

        cd_abs_actual = np.subtract(cd_abs_actual,np.abs(abs_water))

    avg = cd_abs_actual.mean()
    std = cd_abs_actual.std()
    normalized = (cd_abs_actual - avg) / std
    standard_error = std / np.sqrt(len(cd_abs_actual))

    cd_abs_actual = normalized.copy()
    cd_abs_actual = smooth(cd_abs_actual)

    # return np.array(df['WL']), cd_abs_actual, standard_error

    if water:
        return list(range(190,251))[::-1], cd_abs_actual, standard_error
    else:
        return np.array(df['WL']), cd_abs_actual, standard_error

# cd_plot(water)

def normalizeSamples(reference_data, target_samples, normalization_type="min-max"):
    """
    Normalizes target sample spectra based on the reference data.
    
    Parameters:
        reference_data (list or array): Reference data used for normalization (e.g., y-values from a reference spectrum).
        target_samples (dict): Samples to normalize. Assumes each sample has a "spectrum" key with the data to normalize.
        normalization_type (str): Type of normalization ('min-max', 'z-score').
    
    Returns:
        dict: Updated target_samples with an additional "normalized_spectrum" key for each sample.
    """
    reference_data = reference_data.spectrum["cd_abs"]
    if normalization_type == "min-max":
        y_min = min(reference_data)
        y_max = max(reference_data)
        
        def normalize(data):
            return (data - y_min) / (y_max - y_min)
    
    elif normalization_type == "z-score":
        mean = sum(reference_data) / len(reference_data)
        std_dev = (sum((x - mean) ** 2 for x in reference_data) / len(reference_data)) ** 0.5
        
        def normalize(data):
            return (data - mean) / std_dev
    
    else:
        raise ValueError(f"Unsupported normalization type: {normalization_type}")
    
    # Normalize each sample's spectrum
    for sample_id in target_samples:
        spectrum = target_samples[sample_id]["spectrum"].spectrum["cd_abs"]
        normalized_spectrum = normalize(spectrum)
        target_samples[sample_id]["spectrum"].spectrum["cd_abs_normalized"] = normalized_spectrum
    
    return target_samples

def linearAdjustSamples(reference_data, target_samples):
    """
    Adjusts the target sample spectra to match the min and max of the reference data.
    
    Parameters:
        reference_data (array-like): Reference spectrum used to define the min and max.
        target_samples (dict): Samples to adjust. Assumes each sample has a "spectrum" key.
    
    Returns:
        dict: Updated target_samples with an additional "adjusted_spectrum" key for each sample.
    """
    # Calculate min and max of the reference data
    reference_data = reference_data.spectrum["cd_abs"]
    ref_min = min(reference_data)
    ref_max = max(reference_data)
    
    def adjust(data):
        """Adjusts data to stretch to the reference range."""
        data_min = min(data)
        data_max = max(data)
        # Apply linear scaling
        return ref_min + (data - data_min) * (ref_max - ref_min) / (data_max - data_min)
    
    # Adjust each sample's spectrum
    for sample_id in target_samples:
        spectrum = target_samples[sample_id]["spectrum"].spectrum["cd_abs"]
        adjusted_spectrum = adjust(spectrum)
        target_samples[sample_id]["spectrum"].spectrum["cd_abs_normalized"] = adjusted_spectrum
    
    return target_samples



def plotSamples(data, samples = [], show=True, save = False, output_path = "output.png", config = None, style = None, custom_colors=None, create_fig = True):
    
    default_config = {
        "title": None,
        "title_fontsize": 16,
        "xlim": (180, 260),
        "ylim": None,
        "grid": True,
        "grid_linewidth": 0.5,
        "xlabel": "Wavelength (nm)",
        "ylabel": "CD Absorbance (millidegrees)",
        "tick_direction": "in",
        "tick_fontsize": 18, 
        "figsize": (10, 8),
        "line_width": 3,
        "dpi": 500
    }
    
    use_normalized = config.get("use_normalized_spectrum", False) if config else False
    # Merge user config with default config
    if config:
        default_config.update(config)
        
    if create_fig: plt.figure(figsize=default_config["figsize"])
        
    if style:
        if isinstance(style, dict):  # If it's a dictionary, update rcParams
            plt.rcParams.update(style)
        elif isinstance(style, str):  # If it's a string, use the named style
            plt.style.use(style)
    
    if custom_colors:
        if len(custom_colors) == len(samples):
            palette = custom_colors  # Use colors as-is
        else:
            # Generate gradient colors
            cmap = LinearSegmentedColormap.from_list("custom_cmap", custom_colors, N=len(samples))
            palette = [cmap(i) for i in np.linspace(0, 1, len(samples))]
    else:
        # Default to Seaborn's color palette
        palette = sns.color_palette("hls", len(samples))

    for sampleID, color  in zip(samples, palette):
        spectrum =  data[sampleID]["spectrum"]
        # Use the specified color if provided, else use default
        spectrum.plot(color=color, linewidth=default_config["line_width"], use_normalized=use_normalized)  # Assume the spectrum object has a plot method with color support

    # Apply customizations
    plt.title(default_config["title"], fontsize=default_config["title_fontsize"])
    plt.xlabel(default_config["xlabel"])
    plt.ylabel(default_config["ylabel"])
    plt.xlim(default_config["xlim"])
    if default_config["ylim"]:
        plt.ylim(default_config["ylim"])
    plt.tick_params(axis='both', direction=default_config["tick_direction"], labelsize=default_config["tick_fontsize"])
    if default_config["grid"]:
        plt.grid(linewidth=default_config["grid_linewidth"])
    plt.tight_layout()

    # Save or show plot
    if save:
        plt.savefig(output_path, dpi=default_config["dpi"])
    if show:
        plt.show()
        
        
def cdplot(wls,ramp_plot,l):
    sns.color_palette("hls", 8)
    plt.plot(wls,ramp_plot,label=l)

    plt.xlabel('Wavelength (nm)',fontsize=14),plt.ylabel('CD Absorbance (millidegrees)',fontsize=14)
    plt.legend(fontsize=12)    

def rampproc(caminho,temp_min,arqinicio,arqfim,water=False):
    
    actual_temp = temp_min
    wls,ramp_plot,temp,error = [],[],[],[]

    for i in range(arqinicio,arqfim):
        
        df = pd.read_csv(f'{caminho}/1_converted_{i}.txt',names=['WL','CD Abs'],sep=' ')

        if i % 6 == 0 and i != 1 and i != 96:
                        
            cd_abs_actual = [x/6 for x in cd_abs_actual]
            cd_abs_actual = cd_abs_actual[list(df['WL']).index(250):list(df['WL']).index(200)]
            cd_abs_actual = np.array(cd_abs_actual)

            if water:
                # Read Water Before!
                cd_abs_actual = np.subtract(cd_abs_actual,cd_abs_w[30:])

            avg = cd_abs_actual.mean()
            std = cd_abs_actual.std()
            normalized = (cd_abs_actual - avg) / std
            standard_error = std / np.sqrt(len(cd_abs_actual))

            # Descomentar para normalizar
            # cd_abs_actual = smooth(normalized)
            cd_abs_actual = smooth(cd_abs_actual)

            wls.append(list(df['WL']))
            temp.append(actual_temp)
            ramp_plot.append(cd_abs_actual)
            error.append(standard_error)

            actual_temp += 5

        if i % 6 == 0 or i == 1 or i == 96: cd_abs_actual = list(df['CD Abs'])
        else: cd_abs_actual = np.add(np.array(df['CD Abs']),cd_abs_actual)

    return wls[0], temp, np.array(ramp_plot), error

def boltzmann(x, A1, A2, x0, dx):
    return A2 + (A1 - A2) / (1 + np.exp((x - x0) / dx))

def wlspectra(wl,wls,temperaturas,cd_abs,error,TM=False):

    index = wls.index(wl)
    wl_cd = []

    for i in range(len(cd_abs)):
        wl_cd.append(cd_abs[i][index])

    plt.errorbar(temperaturas,wl_cd,yerr=error,capsize=3,fmt="r--o",color="#E0218A",ecolor = "black")
    plt.xlabel('Temperature ºC'),plt.ylabel('CD Abs')

    if TM:
        popt, pcov = curve_fit(boltzmann, temperaturas, wl_cd, p0=[0, 10, 30, 10])
        plt.plot(range(20,100), boltzmann(range(20,100), *popt), color='gray', linestyle='--', label=f'Boltzmann')
        plt.legend()
        print(popt)

def rampplot(wls,temp,ramp_plot):

    from scipy import interpolate

    wls = wls[wls.index(250):wls.index(200)]
    # wls = wls[wls.index(250):wls.index(190)]

    X_b1, Y_b1 = np.meshgrid(wls, temp)
    Z_b1 = np.array([ramp_plot[i] for i in range(len(temp))])

    Q1 = np.percentile(Z_b1, 25, axis=None)
    Q3 = np.percentile(Z_b1, 75, axis=None)
    IQR = Q3 - Q1

    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR

    outliers_mask = (Z_b1 < lower_bound) | (Z_b1 > upper_bound)

    Z_b1_cleaned = Z_b1.copy()
    Z_b1_cleaned[outliers_mask] = np.nan

    X_valid = X_b1[~outliers_mask]
    Y_valid = Y_b1[~outliers_mask]
    Z_valid = Z_b1_cleaned[~outliers_mask]

    interp_func = interpolate.griddata((X_valid, Y_valid), Z_valid, (X_b1, Y_b1), method='linear')

    Z_b1 = np.where(outliers_mask, interp_func, Z_b1)

    fig = plt.figure(figsize=(15,10))
    ax = fig.add_subplot(111, projection='3d')

    im = ax.plot_surface(X_b1, Y_b1, Z_b1, cmap='viridis')

    ax.set_xlabel('\nWavelength (nm)\n',fontsize=18)
    ax.set_ylabel('Temperature (ºC)',fontsize=18)
    # ax.set_zlabel('CD Absorbance (millidegrees)',fontsize=14)

    ax.set_xlim(200,250)

    cbar = fig.colorbar(im,shrink=0.5,pad=0.01,orientation='horizontal')
    cbar.set_label('CD Absorbance (millidegrees)', fontsize=18)
    cbar.ax.tick_params(labelsize=12)

    plt.tight_layout()

    return ax

def rampplot2D(wls,temp,ramp_plot):

    wls = wls[wls.index(250):wls.index(200)]

    X_b1, Y_b1 = np.meshgrid(wls, temp)
    
    min_value = np.min(ramp_plot)
    max_value = np.max(ramp_plot)
    ramp_plot = (ramp_plot - min_value) / (max_value - min_value)

    sc = plt.pcolormesh(X_b1,Y_b1,ramp_plot,cmap='viridis',shading='auto')
    
    plt.xlabel('Wavelength (nm)',fontsize=14)
    plt.ylabel('Temperature (ºC)',fontsize=14)

    # cb = plt.colorbar(sc)
    # cb.set_label('CD Absorbance (millidegrees)', fontsize=14)

    return sc

def pontualspectra(caminho,l,c):

    df = pd.read_csv(f'{caminho}',sep='\t',names=['Temperature','CD Absorbance'])

    temp, cdabs = df['Temperature'],df['CD Absorbance']

    # cdabs = (cdabs - min(cdabs)) / (max(cdabs) - min(cdabs))
    # cdabs = smooth(cdabs)

    plt.plot(temp, cdabs,'-o',label=l,color=c,markersize=4)
