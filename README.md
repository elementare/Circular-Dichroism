
# LumiCD

**LumiCD** is a Python library for processing and analyzing Circular Dichroism (CD) samples. It provides tools for loading, processing, visualizing, and analyzing CD spectra, making it easier to handle experimental data.

## Features

- Load and process CD sample metadata and spectra.
- Apply custom adjustments and filters to spectra.
- Visualize CD spectra with customizable styles and colors.
- Normalize and preprocess experimental data.

## Installation

To install the library:

```bash
pip install LumiCD
```

Or clone the repository and install locally:

```bash
git clone https://github.com/yourusername/LumiCD.git
cd LumiCD
pip install .
```

## Example Usage

Here’s an example of how to use LumiCD to load metadata, process CD spectra, and plot the results:

```python
import LumiCD as cd
from LumiCD.utils import plotSamples, linearAdjustSamples

# Load metadata
metaData = cd.DataLoader.read_toml("metadata.toml")
metaData = metaData.json

# Extract and process spectra
data04_10 = cd.Spectrum.extract(metaData["04-10-2024"], fileType="filesTXT", buffers="NP")
samples04_10Proc = cd.Spectrum.id_to_key(data04_10)

# Plot samples
style = {"font.size": 12, "axes.titlesize": 14, "axes.labelsize": 12}
plotSamples(samples04_10Proc, samples=[1, 3, 5], show=False, style=style, 
            custom_colors=["#29f3e2", "#47ff4d", "#ffff47"])
```

## Input Format

The library expects data in the following format:
- A **TOML** metadata file with details about experiments and samples.
- Spectrum data files (`.txt`, `.csv`, etc.) with columns for wavelengths and ellipticity.

Example spectrum file format:
```
Wavelength (nm), Ellipticity (mdeg)
190, -1.25
191, -1.30
...
```

## Contributing

Contributions are welcome! To contribute:
1. Fork the repository.
2. Create a new branch: `git checkout -b feature-name`.
3. Commit your changes: `git commit -m "Add feature"`.
4. Push to the branch: `git push origin feature-name`.
5. Open a Pull Request.

## License

This project is licensed under the **MIT License**. See the [LICENSE](LICENSE) file for details.

## Contact

For questions or suggestions:
- **Author:** Carlos Daniel and Joaão Otávio
- **Email:** carlos23001@ilum.cnpem.br
- **Repository:** [LumiCD on GitHub](https://github.com/elementare/Circular-Dichroism)
