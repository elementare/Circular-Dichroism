import pandas as pd
import numpy as np
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
import seaborn as sns
from .visualization import Visualization

class Spectrum:
    def __init__(self, sampleData, name, fileType="filesTXT", buffer = False):
        self.sampleData = sampleData
        self.name = name
        self.fileType = fileType
        self.cd_abs = None
        self.wavelengths = None
        self.standard_error = None
        self.buffer = buffer
        self.spectrum = self.process()
        
    def process(self):

        self.cd_abs = self._aggregate_cd_abs()

        if self.buffer:

            buffer_result = self._find_matching_buffer(buffer_samples=self.buffer)
            if buffer_result:
                wl_buffer, abs_buffer, se_buffer = buffer_result
                self._apply_buffer_correction(wl_buffer, abs_buffer)

        
        self.cd_abs = self._smooth(self.cd_abs)
        size = self.sampleData["size"] if "size" in  self.sampleData.keys() else None
        distance = self.sampleData["distance"] if "distance" in  self.sampleData.keys() else None
        pmt = self.sampleData["PMT"] if "PMT" in self.sampleData.keys() else None
        return {
            "WL": self.wavelengths,
            "cd_abs": self.cd_abs,
            "standard_error": self.standard_error,
            "size": size, 
            "distance": distance, 
            "pmt": pmt,
        }
        

    def _aggregate_cd_abs(self):

        cd_abs_actual = None
        # print(f"SampleData: {self.sampleData}")
        if not self.sampleData[self.fileType]:
            raise FileNotFoundError(f"I wasn't able to find the sample files.\nSample: {self.sampleData} ")
        
        for i, file_path in enumerate(self.sampleData[self.fileType], start=1):
            
            df = pd.read_csv(file_path, names=['WL', 'CD Abs'], sep=' ')
            if cd_abs_actual is None:
                cd_abs_actual = df['CD Abs']
            else:
                cd_abs_actual = np.add(cd_abs_actual, np.array(df['CD Abs']))

        self.wavelengths = list(df['WL'])
        return [x / len(self.sampleData[self.fileType]) for x in cd_abs_actual]

    def _find_matching_buffer(self, buffer_samples):
    
        for buffer_sample in buffer_samples:
            if (buffer_sample["distance"] == self.sampleData["distance"] and 
                buffer_sample["PMT"] == self.sampleData["PMT"] and 
                buffer_sample["size"] == self.sampleData["size"] and 
                buffer_sample["id"] != self.sampleData["id"]):
                return self._get_abs_buffer(buffer_sample)
        
        raise FileNotFoundError(f"I wasn't able to find a buffer for {self.sampleData['id']} "
                                f"at distance {self.sampleData['distance']} and PMT {self.sampleData['PMT']}")

    def _apply_buffer_correction(self, wl_buffer, abs_buffer):
        if len(self.wavelengths) > len(wl_buffer):
            self.cd_abs = self.cd_abs[:len(wl_buffer)]
        
        self.cd_abs = np.subtract(self.cd_abs, abs_buffer)

    def _get_abs_buffer(self, buffer_sample):

        df = pd.read_csv(buffer_sample[self.fileType][0], names=['WL', 'CD Abs'], sep=' ')
        wl_buffer = list(df['WL'])
        abs_buffer = np.array(df['CD Abs'])
        se_buffer = 0 
        
        return wl_buffer, abs_buffer, se_buffer
        
    def _smooth(self, data):
        window_size = 11  # Tamanho da janela (deve ser um número ímpar)
        poly_order = 3    # Ordem do polinômio
        data_smooth = savgol_filter(data, window_size, poly_order)
        
        return data_smooth
    
    def plot(self, mode="cdplot", color = None, linewidth = None, use_normalized=False):

        if mode == "cdplot":
            name = "$"+self.name+"$"+ f" - Size {self.spectrum['size']}" if self.spectrum["size"] is not None else self.name
            # name = name + f" - Distance {self.spectrum['distance']}" if self.spectrum['distance'] is not None else name
            # name = name + f" - PMT {self.spectrum['pmt']}" if self.spectrum['pmt'] is not None else name
            cd_abs = self.cd_abs if not use_normalized else self.spectrum["cd_abs_normalized"]
            Visualization.cdplot(self.wavelengths, cd_abs, label=name, color=color, linewidth=linewidth)
    

    @staticmethod
    def extract(metaData, fileType = "filesTXT", buffers = False):
        metaDataProc = {}
        if isinstance(buffers, str):
            buffers = {"default": buffers,
                       "buffersList": [buffers]}
        elif isinstance(buffers, dict):
            buffers["buffersList"] = []
            for key in buffers.keys():
                if  buffers[key] not in buffers["buffersList"]: buffers["buffersList"].append(buffers[key])
        for sampleName in metaData:
            samples = {}

            for sample in metaData[sampleName]["samples"]:
                idSample = sample["id"]
                
                if sampleName in buffers["buffersList"]:
                    buffer = False
                elif sampleName in buffers.keys():
                    buffer = metaData[buffers[sampleName]]["samples"]
                else:
                    buffer = metaData[buffers["default"]]["samples"] 
            
                samples[idSample] = Spectrum(sampleData=sample, name=sampleName, fileType=fileType, buffer=buffer)
            
            metaDataProc[sampleName] = samples
                
            
        return metaDataProc   

    @staticmethod
    def id_to_key(data):
        sampleDF = {}
        for sampleName in data:
            for sampleID in data[sampleName]:

                sampleDF[sampleID] = {"spectrum": data[sampleName][sampleID]}
                sampleDF[sampleID]["name"] = sampleName
        
        return sampleDF
            
