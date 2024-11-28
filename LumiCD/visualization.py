import matplotlib.pyplot as plt
import seaborn as sns

class Visualization:
    
    def __init__(self):
        ...
        
    @staticmethod
    def cdplot(wls,ramp_plot,label,color,linewidth):

        
        plt.plot(wls,ramp_plot,label=label, color=color, linewidth=linewidth,)

        plt.xlabel('Wavelength (nm)',fontsize=16),plt.ylabel('CD Absorbance (millidegrees)', fontsize=16)
        plt.legend()    
        